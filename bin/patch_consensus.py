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

WHY the output stays in the CONSENSUS coordinate frame (this is the important bit):
    The A03 consensus was called against the wildtype reference (rvfv_reference.fa),
    so it already lives in that reference's coordinate system: same segment lengths,
    same headers (NC_014395_S / NC_014396_M / NC_014397_L). The primer BED
    (rvfv_amplicons_v4.bed) is written in those same coordinates, and iVar trim
    matches BED chrom names to the reference headers. Therefore the patched output:
        * keeps the segment IDs NC_014395_S / NC_014396_M / NC_014397_L
        * keeps each segment's length equal to the consensus length
    so the existing BED and the pipeline's default --segments keep working with NO
    changes. MP-12 is only consulted to supply bases at N positions; it never
    changes the coordinate frame (MP-12 indels vs the consensus are absorbed by the
    pairwise alignment and do not shift output positions).

WHAT it does, per segment (S, M, L):
    1. MAFFT-align the MP-12 segment to the A03 consensus segment.
    2. Walk the alignment in consensus coordinates:
         - consensus has a real base  -> keep the A03 base (preserves real calls/SNVs)
         - consensus has an N         -> substitute the aligned MP-12 base if present
         - consensus has a gap ('-')  -> drop (keeps output in consensus frame)
    3. Emit the patched segment under its original NC_ id.

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
DEFAULT_MP12_REF = "/fastscratch/tylerdoe/ARBOR/assets/rvfv/mp12_reference.fa"
DEFAULT_CONSENSUS_DIR = "/fastscratch/tylerdoe/ARBOR/tests/beocat/results/ivar"
DEFAULT_OUT = "/fastscratch/tylerdoe/ARBOR/assets/rvfv/A03_patched_reference.fa"

# output_id (== BED chrom / pipeline --segments), MP-12 header id, consensus filename
SEGMENTS = [
    ("NC_014395_S", "DQ380154.1", "A03_NC_014395_S.fa"),
    ("NC_014396_M", "DQ380208.1", "A03_NC_014396_M.fa"),
    ("NC_014397_L", "DQ375404.1", "A03_NC_014397_L.fa"),
]


def read_fasta(path=None, text=None):
    """Return {header_first_token: sequence_upper}."""
    seqs = {}
    header = None
    seq = []
    lines = open(path).readlines() if path else text.splitlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                seqs[header] = "".join(seq)
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.upper())
    if header is not None:
        seqs[header] = "".join(seq)
    return seqs


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
        return read_fasta(text=result.stdout)
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
    mp12_ref = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_MP12_REF
    consensus_dir = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_CONSENSUS_DIR
    out_file = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_OUT

    ref_seqs = read_fasta(mp12_ref)

    with open(out_file, "w") as out:
        for out_id, mp12_id, cons_file in SEGMENTS:
            if mp12_id not in ref_seqs:
                sys.exit("MP-12 id %s not found in %s" % (mp12_id, mp12_ref))
            ref_seq = ref_seqs[mp12_id]

            cons_seqs = read_fasta(os.path.join(consensus_dir, cons_file))
            cons_seq = list(cons_seqs.values())[0]

            aln = mafft_align(mp12_id, ref_seq, out_id, cons_seq)
            patched = patch_segment(aln[mp12_id], aln[out_id])

            n_before = cons_seq.count("N")
            n_after = patched.count("N")
            print(
                "%-13s consensus=%d  output=%d  N: %d -> %d  (%d patched)"
                % (out_id, len(cons_seq), len(patched),
                   n_before, n_after, n_before - n_after)
            )
            if len(patched) != len(cons_seq):
                print("  WARNING: output length != consensus length for %s; "
                      "primer BED coordinates may be affected." % out_id)

            out.write(">%s\n%s\n" % (out_id, patched))

    print("\nDone. Patched reference written to:\n  %s" % out_file)
    print("Run with:  --reference %s  (BED and --segments unchanged)" % out_file)


if __name__ == "__main__":
    main()
