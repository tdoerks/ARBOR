#!/usr/bin/env python3
"""
make_consensus_reference.py

Turn one sample's ARBOR consensus into a gap-free reference for a new ARBOR run
(e.g. an in-house stock consensus, so the stock's own fixed differences stop showing up as
"variants"). Generalised from the A03 remap script (patch_consensus.py).

Per segment:
  1. Slice the segment out of ARBOR's per-segment iVar consensus. Each results/ivar/<sample>_<seg>.fa
     holds the WHOLE concatenated genome with real bases only in that segment's slice; the bundled
     reference's record order/lengths define the slices (last contig to EOF: iVar may add a base).
  2. MAFFT-align the fill reference's segment (default: true MP-12, assets/rvfv/mp12_reference.fa)
     to the consensus slice.
  3. Walk the alignment in fill-reference coordinates:
       - fill reference has a gap       -> drop (consensus insertion; keeps the coordinate frame)
       - consensus base is real          -> keep it (the sample's own calls / fixed differences)
       - consensus is N or gap           -> use the fill reference's base (ends, primer sites)
  4. Write it under the bundled segment ID (NC_014395_S / NC_014396_M / NC_014397_L), at the
     exact bundled length, so the primer BED and the default --segments keep working.

Usage (needs `mafft` on PATH — on Beocat: module load MAFFT):
    python3 bin/make_consensus_reference.py <sample> [--results results] [--out <sample>_reference.fa]
        [--frame-ref assets/rvfv/rvfv_reference.fa] [--fill-ref assets/rvfv/mp12_reference.fa]
Then run ARBOR with --reference <out>.
"""
import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "assets", "rvfv")

# bundled segment ID  ->  MP-12 record ID in mp12_reference.fa (segments S, M, L)
FILL_IDS = {"NC_014395_S": "DQ380154.1", "NC_014396_M": "DQ380208.1", "NC_014397_L": "DQ375404.1"}


def read_fasta(path=None, text=None):
    recs, name = [], None
    for line in (open(path) if path else text.splitlines()):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs.append([name, []])
        elif name:
            recs[-1][1].append(line.upper())
    return [(n, "".join(s)) for n, s in recs]


def mafft_pair(a_id, a_seq, b_id, b_seq):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
        t.write(f">{a_id}\n{a_seq}\n>{b_id}\n{b_seq}\n")
        path = t.name
    try:
        r = subprocess.run(["mafft", "--quiet", "--auto", path], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"MAFFT failed for {b_id}:\n{r.stderr}")
        return dict(read_fasta(text=r.stdout))
    finally:
        os.unlink(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out")
    ap.add_argument("--frame-ref", default=os.path.join(ASSETS, "rvfv_reference.fa"))
    ap.add_argument("--fill-ref", default=os.path.join(ASSETS, "mp12_reference.fa"))
    args = ap.parse_args()
    out = args.out or f"{args.sample}_reference.fa"

    frame = read_fasta(args.frame_ref)
    offsets, pos = {}, 0
    for cid, seq in frame:
        offsets[cid] = (pos, len(seq))
        pos += len(seq)
    total = pos
    fill = dict(read_fasta(args.fill_ref))

    with open(out, "w") as fo:
        for seg, (start, seglen) in offsets.items():
            fill_id = FILL_IDS.get(seg)
            if fill_id not in fill:
                sys.exit(f"{seg}: fill-reference record {fill_id} not found in {args.fill_ref}")
            path = os.path.join(args.results, "ivar", f"{args.sample}_{seg}.fa")
            if not os.path.exists(path):
                sys.exit(f"{seg}: consensus not found: {path}")
            full = "".join(s for _, s in read_fasta(path))
            if abs(len(full) - total) > 5:
                sys.exit(f"{path}: length {len(full)} vs reference total {total} — wrong frame, aborting")
            end = len(full) if seg == frame[-1][0] else start + seglen
            cons = full[start:end]

            aln = mafft_pair(fill_id, fill[fill_id], seg, cons)
            out_seq, kept, filled, diffs, dropped = [], 0, 0, 0, 0
            for r, c in zip(aln[fill_id], aln[seg]):
                if r == "-":
                    dropped += c not in ("N", "-")
                    continue
                if c in ("N", "-"):
                    out_seq.append(r)
                    filled += 1
                else:
                    out_seq.append(c)
                    kept += 1
                    diffs += c != r
            seq = "".join(out_seq)
            if len(seq) != seglen or len(seq) != len(fill[fill_id]):
                sys.exit(f"{seg}: output {len(seq)} bp != reference {seglen} bp — aborting")
            if "N" in seq:
                sys.exit(f"{seg}: N remains in output — aborting")
            fo.write(f">{seg}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
            print(f"{seg}: {len(seq)} bp | {kept} from {args.sample} consensus, {filled} filled from "
                  f"{fill_id} | {diffs} position(s) differ from {fill_id}"
                  + (f" | {dropped} consensus insertion base(s) dropped" if dropped else ""))
    print(f"\n-> {out}  (use with: --reference {os.path.abspath(out)})")


if __name__ == "__main__":
    main()
