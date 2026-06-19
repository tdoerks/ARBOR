#!/bin/bash
# ARBOR head job for Beocat (KSU).
# The DRIVER: it stays alive and submits each pipeline step as its own SLURM job.
# Reference, primers, and segments default to the bundled RVFV set — you only need a samplesheet.
#
#   Submit from the dir holding samplesheet.csv (ARBOR/tests/beocat):
#     sbatch run_arbor.sbatch
#
#SBATCH --job-name=arbor
#SBATCH --partition=batch.q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=arbor_head_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tyler.tyler.d@gmail.com   # <- change to your preferred address (KSU email is most reliable on Beocat)

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# caches on scratch; reuse a shared per-user image cache so containers aren't re-pulled
export NXF_WORK="$PWD/work"
export NXF_SINGULARITY_CACHEDIR="/fastscratch/$USER/singularity_cache"
export NXF_HOME="$PWD/.nextflow_home"
export NXF_VER=26.04.3
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

# Environment. Beocat's bundled `Nextflow` module is a SELF-CONTAINED 24.04.2 that ignores
# NXF_VER and is too old for ARBOR's current nf-core modules (they need Nextflow >= 25.10).
# Load Singularity + Java + (for its Java dependency) the Nextflow module, then below we fall
# back to a standalone launcher if the available Nextflow is < 25.
module load Singularity 2>/dev/null || module load singularity 2>/dev/null || true
module load OpenJDK     2>/dev/null || module load Java 2>/dev/null || true
module load Nextflow    2>/dev/null || module load nextflow 2>/dev/null || true

# Decide which nextflow to use: prefer one >= 25; otherwise install a standalone launcher
# into scratch (the standalone launcher DOES honor NXF_VER and will run 26.04.3).
NF=nextflow
cur="$(nextflow -version 2>/dev/null | grep -oiE 'version [0-9]+' | grep -oE '[0-9]+' | head -1 || true)"
if [ -z "${cur:-}" ] || [ "${cur:-0}" -lt 25 ]; then
    echo "Available Nextflow is '${cur:-none}' (< 25) — installing standalone ${NXF_VER} into scratch..."
    [ -x "$PWD/nextflow" ] || curl -s https://get.nextflow.io | bash
    NF="$PWD/nextflow"
fi

echo "== Nextflow version in use =="
"$NF" -version 2>&1 | grep -i version || true
echo "============================="

# Run the LOCAL clone (two dirs up) so your `git pull` always takes effect — no remote-cache staleness.
REPO_ROOT="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"
"$NF" run "$REPO_ROOT" -profile beocat \
    --input   samplesheet.csv \
    --outdir  results \
    -work-dir "$NXF_WORK" \
    -resume

# Reference (RVFV MP-12 S/M/L), primer BED, and --segments are bundled defaults.
# For a different virus, add: --reference ref.fa --primer_bed primers.bed --segments id1,id2,...
# To compare against external strains in the trees, add: --context_fasta known_strains.fa
