#!/bin/bash
# ARBOR head job for Beocat (KSU).
# The DRIVER: it stays alive and submits each pipeline step as its own SLURM job.
# Reference, primers, and segments default to the bundled RVFV set — you only need a samplesheet.
#
#   Submit from the dir holding samplesheet.csv:
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

# keep caches on scratch; reuse a shared per-user image cache so containers aren't re-pulled
export NXF_WORK="$PWD/work"
export NXF_SINGULARITY_CACHEDIR="/fastscratch/$USER/singularity_cache"
export NXF_HOME="$PWD/.nextflow_home"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

# Pin a Nextflow >= 24.10 — current nf-core modules don't parse on older runtimes
# (Beocat's module default may be 24.04.x). Nextflow self-downloads this version.
export NXF_VER=26.04.3

module load Nextflow   2>/dev/null || module load nextflow   2>/dev/null || true
module load Singularity 2>/dev/null || module load singularity 2>/dev/null || true

# show what we're actually running (should report 26.04.3, NOT 24.04.2)
echo "== Nextflow version in use =="
nextflow -version 2>&1 | grep -i version || true
echo "============================="

nextflow run tdoerks/ARBOR -r main -profile beocat \
    --input   samplesheet.csv \
    --outdir  results \
    -work-dir "$NXF_WORK" \
    -resume

# Reference (RVFV MP-12 S/M/L), primer BED, and --segments are bundled defaults.
# For a different virus, add: --reference ref.fa --primer_bed primers.bed --segments id1,id2,...
# To compare against external strains in the trees, add: --context_fasta known_strains.fa
