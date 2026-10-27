#!/bin/bash
# Combined run: 96 July 2026 samples + 15 new 10-23-26 samples (111 total)
# Reference: NCBI MP-12 (assets/rvfv/mp12_reference.fa)
#
#   1. mkdir -p /fastscratch/tylerdoe/ARBOR/runs/RVFV_combined_20261026
#   2. cp this file + samplesheet_combined_20261026.csv into that dir as samplesheet.csv
#   3. cd /fastscratch/tylerdoe/ARBOR/runs/RVFV_combined_20261026 && sbatch run_arbor_combined_20261026.sbatch
#
#SBATCH --job-name=arbor-combined-20261026
#SBATCH --partition=batch.q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=arbor_head_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tylerdoe@ksu.edu

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

export NXF_WORK="$PWD/work"
export NXF_SINGULARITY_CACHEDIR="/fastscratch/$USER/singularity_cache"
export NXF_HOME="$PWD/.nextflow_home"
export NXF_VER=26.04.3
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

module load Singularity 2>/dev/null || module load singularity 2>/dev/null || true
module load Java/17.0.15 2>/dev/null || module load Java/23.0.2 2>/dev/null \
    || module load Java/17 2>/dev/null || true

REPO_ROOT="$SLURM_SUBMIT_DIR"
RUN_DIR="/fastscratch/tylerdoe/ARBOR/runs/RVFV_combined_20261026"
NF="$REPO_ROOT/tests/beocat/nextflow"
if [ ! -x "$NF" ]; then
    echo "Installing standalone Nextflow ${NXF_VER}..."
    curl -s https://get.nextflow.io | bash
    mv nextflow "$NF"
fi

echo "== Java + Nextflow in use =="
java -version 2>&1 | head -1 || true
"$NF" -version 2>&1 | grep -i version || true
echo "============================"

"$NF" run "$REPO_ROOT" -profile beocat \
    --input   "$RUN_DIR/samplesheet.csv" \
    --outdir  "$RUN_DIR/results" \
    --reference "$REPO_ROOT/assets/rvfv/mp12_reference.fa" \
    -work-dir "$NXF_WORK" \
    -resume
