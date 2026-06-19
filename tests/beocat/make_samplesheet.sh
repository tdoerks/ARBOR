#!/bin/bash
# Build an ARBOR samplesheet from a directory of Illumina BCLConvert fastqs.
# Handles names like:  R3D720_S5_L001_R1_001.fastq.gz
#   sample  = everything before _S<n>   (e.g. R3D720)
#
#   Usage:  bash make_samplesheet.sh /path/to/fastq_dir > samplesheet.csv
set -euo pipefail
dir=${1:?usage: make_samplesheet.sh <fastq_dir>}

echo "sample,fastq_1,fastq_2"
shopt -s nullglob
declare -A seen
for r1 in "$dir"/*_R1_001.fastq.gz; do
    base0=$(basename "$r1")
    # skip BCLConvert's Undetermined (unassigned barcodes) — not a real sample
    if [[ "$base0" == Undetermined_* ]]; then
        echo "INFO: skipping $base0 (unassigned reads)" >&2
        continue
    fi
    r2=${r1/_R1_001/_R2_001}
    if [[ ! -e "$r2" ]]; then
        echo "WARN: no R2 mate for $r1 — skipping" >&2
        continue
    fi
    base=$(basename "$r1")
    sample=${base%%_S[0-9]*}                 # strip _S<n>_L00x_R1_001.fastq.gz
    if [[ -n "${seen[$sample]:-}" ]]; then
        echo "WARN: '$sample' already seen (multiple lanes?). Only the first lane is used; merge lanes first if needed: $r1" >&2
        continue
    fi
    seen[$sample]=1
    printf '%s,%s,%s\n' "$sample" "$(readlink -f "$r1")" "$(readlink -f "$r2")"
done
