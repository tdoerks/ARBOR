#!/usr/bin/env python3
"""
patch_consensus.py

Build a new ARBOR reference from the A03 iVar consensus, filling the N-runs
(primer regions + sub-threshold-depth positions) with MP-12 reference bases.

WHY N-patching (industry standard for viral tiling-amplicon work, e.g. ARTIC):
    The consensus has N's wherever coverage fell below the min-depth threshold or
    where primers were trimmed. You have no real data there, so you borrow bases
    from a trusted reference to produce a complete, gap-free reference that a
    read aligner can index and that keeps downstream coordinates stable.

IMPORTANT quirk of ARBOR's per-segment iVar consensus files:
    Each ivar/A03_<seg>.fa file is NOT a single segment -- it is the WHOLE
    concatenated genome (M+L+S = 11979 bp, the header order of rvfv_reference.fa),
    with real bases only in that segment's slice and N everywhere else. So this
    script first SLICES each segment's real region out of its file using the
    wildtype reference (rvfv_reference.fa) to define the concat order and per-contig
    lengths, and only then patches. Without slicing you would emit three bloated
    11979-bp sequences and destroy the coordinate frame.

WHY the output is written in the MP-12 REFERENCE coordinate frame:
    MP-12 and the wildtype reference have identical per-segment lengths
    (S=1690, M=3885, L=6404), and the primer BED (rvfv_amplicons_v4.bed) is written
    in those per-contig coordinates. iVar consensus can add a base (the S file came
    out 11980 vs 11979), so writing the output in the MP-12 frame -- exactly one
    output base per non-gap MP-12 column -- GUARANTEES each segment is its exact
    wildtype length, keeping the BED and the pipeline's default --segments valid.
    The tradeoff: an A03 insertion relative to MP-12 is dropped, which is the right
    call for coordinate stability (a lone consensus insertion is low-confidence).
    Output keeps the segment IDs NC_014395_S / NC_014396_M / NC_014397_L.

WHAT it does, per segment (S, M, L):
    1. Slice the segment's real region out of the full-genome ivar file, using the
       wildtype reference for the concat order / offsets (last contig sliced to EOF
       to absorb iVar's extra base).
    2. MAFFT-align the MP-12 segment to that A03 consensus slice.
    3. Walk the alignment in MP-12 (reference) coordinates:
         - MP-12 has a gap ('-')      -> drop (consensus insertion; keeps ref frame)
         - consensus base is real     -> keep the A03 base (preserves real calls/SNVs)
         - consensus is N or gap      -> substitute the aligned MP-12 base
    4. Emit the patched segment (exact wildtype length) under its original NC_ id.

USAGE (on Beocat, from anywhere):
    module load MAFFT 2>/dev/null || module load mafft 2>/dev/null   # provide `mafft`
    python3 bin/patch_consensus.py

    Optional overrides (positional):
        python3 bin/patch_consensus.py <mp12_ref.fa> <consensus_dir> <out.fa>

OUTPUT:
    A03_patched_reference.fa  (3 segments, ready to pass as --reference)
    Prints per-segment: consensus length, output length (must match), N's patched.
"""

import subprocess
import os
import sys
import tempfile

# ---- Defaults (Beocat absolute paths) ---------------------------------------
DEFAULT_WT_REF = "/fastscratch/tylerdoe/ARBOR/assets/rvfv/rvfv_reference.fa"
DEFAULT_MP12_REF = "/fastscratch/tylerdoe/ARBOR/assets/rvfv/mp12_reference.fa"
DEFAULT_CONSENSUS_DIR = "/fastscratch/tylerdoe/ARBOR/tests/beocat/results/ivar"
DEFAULT_OUT = "/fastscratch/tylerdoe/ARBOR/assets/rvfv/A03_patched_reference.fa"

# output_id (== wildtype contig id / BED chrom / pipeline --segments),
# MP-12 header id, per-segment iVar consensus filename
SEGMENTS = [
    ("NC_014395_S", "DQ380154.1", "A03_NC_014395_S.fa"),
    ("NC_014396_M", "DQ380208.1", "A03_NC_014396_M.fa"),
    ("NC_014397_L", "DQ375404.1", "A03_NC_014397_L.fa"),
]


def read_fasta(path=None, text=None):
    """Return an ordered list of (header_first_token, sequence_upper)."""
    records = []
    header = None
    seq = []
    lines = open(path).readlines() if path else text.splitlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.upper())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def fasta_dict(path=None, text=None):
    """read_fasta as a {header: seq} dict."""
    return dict(read_fasta(path=path, text=text))


def concat_offsets(wt_records):
    """Cumulative start offset of each contig in wildtype file (concat) order."""
    offsets = {}
    pos = 0
    for cid, seq in wt_records:
        offsets[cid] = (pos, len(seq))
        pos += len(seq)
    return offsets, pos


def mafft_align(ref_id, ref_seq, cons_id, cons_seq):
    """Pairwise-align two sequences with MAFFT; return {id: aligned_seq}."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp.write(">%s\n%s\n>%s\n%s\n" % (ref_id, ref_seq, cons_id, cons_seq))
        tmp_path = tmp.name
    try:
        result = subprocess.run(
            ["mafft", "--quiet", "--auto", tmp_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            sys.exit("MAFFT failed for %s:\n%s" % (cons_id, result.stderr))
        return fasta_dict(text=result.stdout)
    finally:
        os.unlink(tmp_path)


def patch_segment(ref_aln, cons_aln):
    """Walk the alignment in REFERENCE (MP-12) coordinates.

    Output length == number of non-gap MP-12 columns == len(MP-12 segment), so the
    result is always the exact wildtype segment length and the primer BED stays
    valid regardless of insertions iVar may have introduced in the consensus.
        - MP-12 has a gap ('-')      -> drop (consensus insertion; keeps ref frame)
        - consensus base is real     -> keep the A03 base (preserves SNVs)
        - consensus is N or gap      -> fill from the MP-12 base
    """
    out = []
    for r, c in zip(ref_aln, cons_aln):
        if r == "-":
            continue                      # insertion in consensus -> drop, stay in ref frame
        if c in ("N", "-"):
            out.append(r)                 # no confident A03 base here -> fill from MP-12
        else:
            out.append(c)                 # keep the real A03 base
    return "".join(out)


def main():
    wt_ref = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_WT_REF
    mp12_ref = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_MP12_REF
    consensus_dir = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_CONSENSUS_DIR
    out_file = sys.argv[4] if len(sys.argv) > 4 else DEFAULT_OUT

    # Wildtype reference defines the consensus coordinate frame: concat order
    # (its header order) and each contig's length.
    wt_records = read_fasta(wt_ref)
    offsets, total_len = concat_offsets(wt_records)
    last_contig = wt_records[-1][0]   # slice this one to EOF (absorbs iVar's +1)
    print("Consensus frame (from %s):" % os.path.basename(wt_ref))
    for cid, seq in wt_records:
        start, length = offsets[cid]
        print("  %-13s [%d:%d]  len=%d" % (cid, start, start + length, length))
    print("  total = %d\n" % total_len)

    mp12_seqs = fasta_dict(mp12_ref)

    with open(out_file, "w") as out:
        for out_id, mp12_id, cons_file in SEGMENTS:
            if out_id not in offsets:
                sys.exit("%s not found in wildtype ref %s" % (out_id, wt_ref))
            if mp12_id not in mp12_seqs:
                sys.exit("MP-12 id %s not found in %s" % (mp12_id, mp12_ref))

            start, seg_len = offsets[out_id]
            mp12_seq = mp12_seqs[mp12_id]

            # The ivar file is the FULL concatenated genome; slice this segment out.
            # iVar may add a base (S came out 11980 vs 11979), so slice the LAST
            # contig to EOF and tolerate a small total-length drift -- the ref-frame
            # patching below forces the output back to the exact segment length.
            full = list(fasta_dict(os.path.join(consensus_dir, cons_file)).values())[0]
            drift = len(full) - total_len
            if abs(drift) > 5:
                sys.exit(
                    "%s: ivar file length %d differs from wildtype total %d by %d; "
                    "concat frame mismatch, aborting." % (cons_file, len(full), total_len, drift)
                )
            end = len(full) if out_id == last_contig else start + seg_len
            cons_seg = full[start:end]

            aln = mafft_align(mp12_id, mp12_seq, out_id, cons_seg)
            patched = patch_segment(aln[mp12_id], aln[out_id])

            # Ref-frame output must equal the MP-12 segment length (== wildtype len).
            if len(patched) != len(mp12_seq) or len(patched) != seg_len:
                sys.exit(
                    "%s: patched length %d != MP-12 %d / wildtype %d; BED coordinates "
                    "would break, aborting." % (out_id, len(patched), len(mp12_seq), seg_len)
                )

            n_before = cons_seg.count("N")
            n_after = patched.count("N")
            print(
                "%-13s segment_len=%d  output=%d  N: %d -> %d  (%d patched)"
                % (out_id, seg_len, len(patched),
                   n_before, n_after, n_before - n_after)
            )
            out.write(">%s\n%s\n" % (out_id, patched))

    print("\nDone. Patched reference written to:\n  %s" % out_file)
    print("Run with:  --reference %s  (BED and --segments unchanged)" % out_file)


if __name__ == "__main__":
    main()
