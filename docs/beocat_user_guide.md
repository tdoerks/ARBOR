# ARBOR — Beocat User Guide

How to run **ARBOR** (viral tiling-amplicon: primer-trimmed mapping → iSNVs → per-segment consensus →
phylogeny) on KSU [Beocat](https://www.k-state.edu/hpc/), from Illumina paired-end `fastq.gz`.

The RVFV MP-12 reference, primer BED, and `--segments` are **bundled defaults** — you only need a
samplesheet to run.

> Launch everything from a Beocat **login node** (any `icr-*` host, e.g. `icr-helios`). A small head job
> drives Nextflow, which submits each step as its own SLURM job under the `beocat` profile.

> **First time on Beocat?** Skim [Good to know](#good-to-know-first-timers) at the bottom first — it
> covers the handful of things that surprise new users (where reads go, `Undetermined`, the `TDOERKS_`
> label, and copying results off scratch).

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
Put them on **scratch, not your home dir** — home has a small quota, and `/fastscratch` is the big, fast
space meant for job I/O. First make the folder (on Beocat):

```bash
mkdir -p /fastscratch/$USER/reads
```

Then copy from the machine where the data lives. From **Windows PowerShell** (or use WinSCP / MobaXterm),
replacing `<user>`, `<beocat-login>`, and the local path with your own:

```powershell
# If the *.fastq.gz sit directly in the BCLConvert output folder:
scp "C:\path\to\BCLConvert_output\*.fastq.gz" <user>@<beocat-login>:/fastscratch/<user>/reads/

# If they're in per-sample subfolders, flatten them into one reads/ dir instead:
scp -r "C:\path\to\BCLConvert_output\*\*.fastq.gz" <user>@<beocat-login>:/fastscratch/<user>/reads/
```

`<beocat-login>` is the same host you SSH into (e.g. `headnode.beocat.ksu.edu`). BCLConvert's
`Undetermined_*` files can come along — the samplesheet builder skips them automatically.
End state: `/fastscratch/$USER/reads/` full of `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`.

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

Optional — **email notifications**: edit `run_arbor.sbatch` and uncomment the two `##SBATCH --mail-*`
lines, setting `--mail-user` to your own KSU address. Left off by default. Then:

```bash
sbatch run_arbor.sbatch
```

That's it — reference, primers, and segment names are bundled defaults. The head job runs your **local
clone** (so a `git pull` always takes effect) under the `beocat` profile, roughly equivalent to:

```bash
nextflow run . -profile beocat --input samplesheet.csv --outdir results -resume
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
| `--context_dir` | none | dir of per-segment outgroup/reference FASTAs (`S.fasta`/`M.fasta`/`L.fasta`) folded into each segment's tree; see [`assets/context/README.md`](../assets/context/README.md) |
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

---

## Good to know (first-timers)

Things that surprise people on their first ARBOR run — none are errors:

- **Reads go on `/fastscratch/$USER/`, not your home dir.** Home has a small quota; scratch is the fast,
  roomy space for job I/O (work dir, container cache, and `results/` all live there too).
- **`/fastscratch` is auto-purged** after a period of inactivity. When the run finishes, **copy `results/`
  somewhere permanent** (home, lab storage) — scratch is for active runs, not keeping.
- **`Undetermined_*` fastqs are skipped automatically** — you'll see `INFO: skipping Undetermined_...`.
  That's expected (those are unassigned barcodes, not a sample).
- **The process list shows `TDOERKS_ARBOR:ARBOR:...`.** That `TDOERKS_` prefix is just the pipeline's
  *name* (from `manifest.name = 'tdoerks/arbor'` — the GitHub repo), like nf-core pipelines show
  `NFCORE_...`. It is **not** an account or a path: your job runs as **you** and writes under **your**
  `/fastscratch/$USER/`. Verify any time with `echo $USER`, `pwd`, and `squeue -u $USER`.
- **The first run is quiet for several minutes.** The head job installs Nextflow and pulls the
  Singularity containers before any step starts — the log will sit still, then fill in. Later runs reuse
  the cache and start fast.
- **`squeue` drip-feeds jobs.** The `beocat` profile throttles submission (≤100 queued, ≤10/min) so a
  big run doesn't flood the scheduler — you won't see all tasks at once. Normal.
- **Only per-user edit needed** is the optional email in `run_arbor.sbatch`; everything else keys off
  `$USER`, so the scripts work for anyone without renaming.

### Explore the results
The dashboard at [`dashboard/arbor_dashboard.html`](../dashboard/arbor_dashboard.html) opens in any
browser (no install) — drag in the run's `ivar/*.tsv`, `lofreq/*.vcf`, `mosdepth/*`, `iqtree/*.treefile`
(and optionally `pipeline_info/params.json`) to browse variants, coverage, trees, and the parameters used.
