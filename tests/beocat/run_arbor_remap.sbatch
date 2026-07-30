#!/bin/bash
# ARBOR head job for Beocat (KSU) -- REMAP run against the A03 patched reference.
# Same driver as run_arbor.sbatch, but:
#   * --reference points at assets/rvfv/A03_patched_reference.fa (A03 consensus with
#     primer/low-coverage N's filled from MP-12; built by bin/patch_consensus.py).
#   * outdir + work dir are dedicated to this run so the earlier A03-only `results/`
#     and `work/` are left untouched.
# The primer BED (rvfv_amplicons_v4.bed) and --segments stay at their defaults --
# the patched reference keeps the NC_014395_S / NC_014396_M / NC_014397_L ids and
# exact segment lengths, so those remain valid.
#
#   Submit from ARBOR/tests/beocat (the dir holding samplesheet.csv):
#     sbatch run_arbor_remap.sbatch
#
#SBATCH --job-name=arbor_remap
#SBATCH --partition=batch.q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=arbor_remap_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tylerdoe@ksu.edu   # KSU address — Beocat's mail relay delivers reliably here

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# Dedicated caches/work for the remap run so the A03-only results aren't clobbered.
export NXF_WORK="$PWD/work_remap_A03ref"
export NXF_SINGULARITY_CACHEDIR="/fastscratch/$USER/singularity_cache"
export NXF_HOME="$PWD/.nextflow_home"
export NXF_VER=26.04.3
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

# Environment.
#  - Beocat's Nextflow modules are only 23.10 / 24.04 — both too old for ARBOR's current nf-core
#    modules (need >= 25.10) AND they pin Java 11. So we DON'T load the Nextflow module.
#  - Instead: load Java 17+ (Nextflow 26 needs Java 17-26; Beocat default is Java 11) and use a
#    standalone Nextflow launcher (honors NXF_VER=26.04.3, uses the Java 17 we load here).
module load Singularity 2>/dev/null || module load singularity 2>/dev/null || true
module load Java/17.0.15 2>/dev/null || module load Java/23.0.2 2>/dev/null \
    || module load Java/17 2>/dev/null || true

NF="$PWD/nextflow"
if [ ! -x "$NF" ]; then
    echo "Installing standalone Nextflow ${NXF_VER} into scratch..."
    curl -s https://get.nextflow.io | bash
fi

echo "== Java + Nextflow in use (need Java 17+, Nextflow 26.04.3) =="
java -version 2>&1 | head -1 || true
"$NF" -version 2>&1 | grep -i version || true
echo "============================================================="

# Run the LOCAL clone (two dirs up) so your `git pull` always takes effect — no remote-cache staleness.
REPO_ROOT="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"
"$NF" run "$REPO_ROOT" -profile beocat \
    --input     samplesheet.csv \
    --reference "$REPO_ROOT/assets/rvfv/A03_patched_reference.fa" \
    --outdir    results_remap_A03ref \
    -work-dir   "$NXF_WORK" \
    -resume

# primer BED (rvfv_amplicons_v4.bed) and --segments (NC_014395_S,NC_014396_M,NC_014397_L)
# stay at their bundled defaults; the patched reference preserves those ids/lengths.
# To compare against external strains in the trees, add: --context_fasta known_strains.fa
