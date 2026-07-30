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

WHY the output stays in the CONSENSUS coordinate frame (this is the important bit):
    The A03 consensus was called against the wildtype reference (rvfv_reference.fa),
    so each segment slice already lives in that reference's coordinate system: same
    per-contig lengths, same headers (NC_014395_S / NC_014396_M / NC_014397_L). The
    primer BED (rvfv_amplicons_v4.bed) is written in those same per-contig
    coordinates, and iVar trim matches BED chrom names to the reference headers.
    Therefore the patched output:
        * keeps the segment IDs NC_014395_S / NC_014396_M / NC_014397_L
        * keeps each segment's length equal to the wildtype contig length
    so the existing BED and the pipeline's default --segments keep working with NO
    changes. MP-12 is only consulted to supply bases at N positions; it never
    changes the coordinate frame (MP-12 indels vs the consensus are absorbed by the
    pairwise alignment and do not shift output positions).

WHAT it does, per segment (S, M, L):
    1. Slice the segment's real region out of the full-genome ivar file, using the
       wildtype reference for the concat order / offsets.
    2. MAFFT-align the MP-12 segment to that A03 consensus slice.
    3. Walk the alignment in consensus coordinates:
         - consensus has a real base  -> keep the A03 base (preserves real calls/SNVs)
         - consensus has an N         -> substitute the aligned MP-12 base if present
         - consensus has a gap ('-')  -> drop (keeps output in consensus frame)
    4. Emit the patched segment under its original NC_ id.

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
    """Walk the alignment in consensus coordinates, filling N's from ref."""
    out = []
    for r, c in zip(ref_aln, cons_aln):
        if c == "-":
            continue                      # gap in consensus -> stay in consensus frame
        if c == "N" and r not in ("-", "N"):
            out.append(r)                 # fill N with the aligned MP-12 base
        else:
            out.append(c)                 # keep the real A03 base (or unfillable N)
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

            # The ivar file is the FULL concatenated genome; slice this segment out.
            full = list(fasta_dict(os.path.join(consensus_dir, cons_file)).values())[0]
            if len(full) != total_len:
                sys.exit(
                    "%s: ivar file length %d != wildtype total %d; concat frame "
                    "mismatch, aborting." % (cons_file, len(full), total_len)
                )
            cons_seg = full[start:start + seg_len]

            aln = mafft_align(mp12_id, mp12_seqs[mp12_id], out_id, cons_seg)
            patched = patch_segment(aln[mp12_id], aln[out_id])

            if len(patched) != seg_len:
                sys.exit(
                    "%s: patched length %d != expected segment length %d; BED "
                    "coordinates would break, aborting." % (out_id, len(patched), seg_len)
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
