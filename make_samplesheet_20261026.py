import os
from collections import defaultdict

base    = "/fastscratch/tylerdoe/ARBOR/runs"
fq_dir  = base + "/fastq_20261023"
old     = base + "/RVFV_inhouse_MP12/samplesheet.csv"
out_ss  = base + "/samplesheet_combined_20261026.csv"
out_sh  = base + "/merge_fastq_20261026.sh"
mrg_dir = base + "/merged_fastq_20261026"

# Read old samplesheet grouped by sample name (skip header)
groups = defaultdict(list)   # sample -> [(r1, r2), ...]
order  = []
with open(old) as f:
    for i, line in enumerate(f):
        line = line.rstrip()
        if not line:
            continue
        parts = line.split(",")
        if i == 0 and parts[0].lower() in ("sample", "sample_name", "id"):
            continue
        s, r1, r2 = parts[0], parts[1], parts[2]
        if s not in groups:
            order.append(s)
        groups[s].append((r1, r2))

# Add new _10-23-26 samples
new_10_23 = []
for fname in sorted(os.listdir(fq_dir)):
    if not fname.endswith("_R1_001.fastq.gz"):
        continue
    s  = fname.split("_S")[0] + "_10-23-26"
    r1 = fq_dir + "/" + fname
    r2 = r1.replace("_R1_", "_R2_")
    if s not in groups:
        order.append(s)
        new_10_23.append(s)
    groups[s].append((r1, r2))

# Add _combined samples: original + 10-23-26 reads catted together (max coverage)
for s in new_10_23:
    base_name = s[:-len("_10-23-26")]
    if base_name in groups:
        combined = base_name + "_combined"
        if combined not in groups:
            order.append(combined)
        groups[combined] = groups[base_name] + groups[s]

# Build samplesheet + merge shell script
ss_rows  = ["sample,fastq_1,fastq_2"]
sh_lines = ["#!/bin/bash", "set -euo pipefail", "mkdir -p " + mrg_dir, ""]

for s in order:
    pairs = groups[s]
    if len(pairs) == 1:
        ss_rows.append(f"{s},{pairs[0][0]},{pairs[0][1]}")
    else:
        mr1 = f"{mrg_dir}/{s}_R1_001.fastq.gz"
        mr2 = f"{mrg_dir}/{s}_R2_001.fastq.gz"
        sh_lines.append(f"# {s}: merging {len(pairs)} inputs")
        sh_lines.append("cat " + " ".join(p[0] for p in pairs) + " > " + mr1)
        sh_lines.append("cat " + " ".join(p[1] for p in pairs) + " > " + mr2)
        sh_lines.append("")
        ss_rows.append(f"{s},{mr1},{mr2}")

with open(out_ss, "w") as f:
    f.write("\n".join(ss_rows) + "\n")
with open(out_sh, "w") as f:
    f.write("\n".join(sh_lines) + "\n")

dupes    = [s for s in order if len(groups[s]) > 1]
combined = [s for s in order if s.endswith("_combined")]
print(f"Wrote {len(ss_rows) - 1} samples (+ header) -> {out_ss}")
print(f"  {len(order) - len(new_10_23) - len(combined)} original  |  {len(new_10_23)} 10-23-26  |  {len(combined)} combined")
print(f"Merged {len(dupes)} multi-row samples: {', '.join(dupes)}")
print(f"Cat commands -> {out_sh}")
print(f"\nNext steps on Beocat:")
print(f"  bash {out_sh}")
print(f"  cp {out_ss} /fastscratch/tylerdoe/ARBOR/runs/RVFV_combined_20261026/samplesheet.csv")
print(f"  # then resubmit with -resume")
