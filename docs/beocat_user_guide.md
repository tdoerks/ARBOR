# ARBOR — Beocat User Guide

How to run **ARBOR** (viral tiling-amplicon: primer-trimmed mapping → iSNVs → per-segment consensus →
phylogeny) on KSU [Beocat](https://www.k-state.edu/hpc/), from Illumina paired-end `fastq.gz`.

The RVFV MP-12 reference, primer BED, and `--segments` are **bundled defaults** — you only need a
samplesheet to run.

> Launch everything from a Beocat **login node** (e.g. `icr-helios`). A small head job drives Nextflow,
> which submits each step as its own SLURM job under the `beocat` profile.

---

## What it does

```
reads ─► FastQC ─► fastp ─► Bowtie2(--local) ─► sort+index ─┬─► QC (samtools stats, mosdepth)
                                                            ▼
                                          iVar trim (primer BED) ─► sort+index ─┬─► QC
                                                                                ▼
                          ┌──────────────── whole-reference BAM ───────────────┐
                          ▼                                                     ▼
              iVar variants + LoFreq  (iSNVs)        split BAM by segment (S/M/L)
                                                                  ▼
                                                       iVar consensus (per segment)
                                                                  ▼
                                         group by segment ─► MAFFT(+context) ─► IQ-TREE (per segment)
```

---

## 1. Upload your reads to Beocat

Your fastqs are local (BaseSpace/BCLConvert output, named like `R3D720_S5_L001_R1_001.fastq.gz`).
Get them onto scratch. From **Windows PowerShell** (or use WinSCP / MobaXterm):

```powershell
# adjust the local path and your Beocat username
scp -r "C:\Users\tdoerks\Downloads\BaseSpace\Nathali_run-507641134\*\*.fastq.gz" `
    tylerdoe@<beocat-login>:/fastscratch/tylerdoe/reads/
```

The per-sample BaseSpace subfolders are fine — copy all the `*.fastq.gz` into one `reads/` directory on
scratch. End state: `/fastscratch/$USER/reads/` full of `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`.

---

## 2. Get the pipeline + build the samplesheet

```bash
cd /fastscratch/$USER
git clone https://github.com/tdoerks/ARBOR.git
cd ARBOR/tests/beocat

# build the samplesheet from your reads folder (handles the _S#_L00#_R1_001 naming)
bash make_samplesheet.sh /fastscratch/$USER/reads > samplesheet.csv
column -t -s, samplesheet.csv          # eyeball: one row per sample, R1 + R2 filled
```

The generator derives the sample name as everything before `_S<n>` (so `R3D720_S5_L001_R1_001.fastq.gz`
→ `R3D720`). If a sample spans **multiple lanes** (L001 *and* L002), it warns and uses the first lane —
merge lanes beforehand if you need both.

---

## 3. Launch

```bash
sbatch run_arbor.sbatch
```

That's it — reference, primers, and segment names are bundled defaults. The head job runs:

```bash
nextflow run tdoerks/ARBOR -r main -profile beocat --input samplesheet.csv --outdir results -resume
```

---

## 4. Monitor

```bash
squeue -u $USER             # head job + the per-step jobs it spawns
tail -f arbor_head_*.log    # live Nextflow progress
```

---

## 5. Results (`results/`)

| Path | What |
|---|---|
| `results/multiqc/multiqc_report.html` | combined QC (FastQC, fastp, samtools stats, mosdepth) |
| `results/mosdepth/` | per-sample coverage — check for amplicon "valleys" (failed primer pairs) |
| `results/ivar/` | iVar variant tables (`.tsv`) — frequency-based iSNVs |
| `results/lofreq/` | LoFreq VCFs — significance-tested low-frequency variants |
| `results/ivarconsensus/` (or `ivar/`) | per-segment consensus FASTA per sample |
| `results/iqtree/*.treefile` | **one core-genome tree per segment** (S, M, L) |

---

## Defaults you can override

| Flag | Default | Notes |
|---|---|---|
| `--reference` | bundled RVFV MP-12 | swap for another virus |
| `--primer_bed` | bundled `rvfv_amplicons_v4.bed` | your tiling scheme |
| `--segments` | `NC_014395_S,NC_014396_M,NC_014397_L` | **must match reference FASTA headers** |
| `--context_fasta` | none | external NCBI strains to fold into the trees |
| `--ivar_variants_min_freq` / `_min_depth` / `_min_qual` | 0.01 / 250 / 20 | iVar iSNV thresholds |
| `--lofreq_sig` / `--lofreq_min_cov` | 0.005 / 100 | LoFreq thresholds |
| `--consensus_min_depth` | 10 | below this → N |
| `--iqtree_bootstrap` / `--iqtree_model` | 1000 / MFP | phylogeny |
| `--skip_lofreq` / `--skip_ivar_variants` / `--skip_phylogeny` | false | drop a stage |

---

## Troubleshooting

**`nextflow`/`singularity` not found** — Beocat module names differ; find them with `module spider
nextflow` / `module spider singularity` and update the two `module load` lines in `run_arbor.sbatch`.

**Network errors on first run** — container images download on first use. If compute nodes have no
internet, pre-pull on a login node (ask and we'll set up `NXF_SINGULARITY_CACHEDIR` staging).

**Empty/low coverage on a segment** — check `results/mosdepth/`; coverage valleys indicate a failed
primer pair in that amplicon. The trees still build per segment from whatever consensus is recovered.

**Re-running** — `-resume` is on; just `sbatch run_arbor.sbatch` again — completed steps are cached.
