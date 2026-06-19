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

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# keep caches on scratch
export NXF_WORK="$PWD/work"
export NXF_SINGULARITY_CACHEDIR="$PWD/singularity_cache"
export NXF_HOME="$PWD/.nextflow_home"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

module load Nextflow   2>/dev/null || module load nextflow   2>/dev/null || true
module load Singularity 2>/dev/null || module load singularity 2>/dev/null || true

nextflow run tdoerks/ARBOR -r main -profile beocat \
    --input   samplesheet.csv \
    --outdir  results \
    -work-dir "$NXF_WORK" \
    -resume

# Reference (RVFV MP-12 S/M/L), primer BED, and --segments are bundled defaults.
# For a different virus, add: --reference ref.fa --primer_bed primers.bed --segments id1,id2,...
# To compare against external strains in the trees, add: --context_fasta known_strains.fa
