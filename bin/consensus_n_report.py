#!/usr/bin/env python3
"""
consensus_n_report.py

Where are the N's in a sample's ARBOR consensus? For each segment, list every run of N with its
position (1-based, segment coordinates), length, location (5' end / 3' end / internal) and any
primer it overlaps. Use it to decide whether a consensus can become a reference (end gaps are
normal for tiling amplicons and are usually filled from a reference) or has internal dropouts.

ARBOR's per-segment iVar consensus files (results/ivar/<sample>_<segment>.fa) each hold the WHOLE
concatenated genome, real bases only in that segment's slice — so the segment is sliced out using
the reference FASTA's record order and lengths first.

Usage:
    python3 bin/consensus_n_report.py <sample> [results_dir] [reference.fa] [primers.bed]
    e.g. python3 ../../bin/consensus_n_report.py MP12_merged results
Defaults: results_dir=results, reference/primers = the bundled assets/rvfv files.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sample = sys.argv[1] if len(sys.argv) > 1 else sys.exit(__doc__)
results = sys.argv[2] if len(sys.argv) > 2 else "results"
ref_fa = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "..", "assets", "rvfv", "rvfv_reference.fa")
bed = sys.argv[4] if len(sys.argv) > 4 else os.path.join(HERE, "..", "assets", "rvfv", "rvfv_amplicons_v4.bed")


def read_fasta(path):
    recs, name = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs.append([name, []])
        elif name:
            recs[-1][1].append(line.upper())
    return [(n, "".join(s)) for n, s in recs]


ref = read_fasta(ref_fa)
offsets, pos = {}, 0
for name, seq in ref:
    offsets[name] = (pos, len(seq))
    pos += len(seq)

primers = []
if os.path.exists(bed):
    for line in open(bed):
        f = line.split("\t")
        if len(f) >= 4:
            primers.append((f[0], int(f[1]), int(f[2]), f[3].strip()))   # BED: 0-based, end-exclusive

for seg, (off, seglen) in offsets.items():
    path = os.path.join(results, "ivar", f"{sample}_{seg}.fa")
    if not os.path.exists(path):
        print(f"{seg}: {path} not found")
        continue
    cons = "".join(s for _, s in read_fasta(path))
    last = seg == ref[-1][0]
    sl = cons[off:] if last else cons[off:off + seglen]     # last contig to EOF (iVar may add a base)
    runs, i = [], 0
    while i < len(sl):
        if sl[i] == "N":
            j = i
            while j < len(sl) and sl[j] == "N":
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    total = sum(j - i for i, j in runs)
    print(f"\n{seg}: {len(sl)} bp (reference {seglen}), {total} N in {len(runs)} run(s), "
          f"{100 * (1 - total / max(1, len(sl))):.1f}% called")
    for i, j in runs:
        where = "5' END" if i == 0 else "3' END" if j >= len(sl) else "internal"
        hit = [p[3] for p in primers if p[0] == seg and p[1] < j and p[2] > i]
        print(f"  {i + 1:>6}-{j:<6} {j - i:>5} N  {where:<8}  {('primer: ' + ', '.join(hit)) if hit else ''}")
