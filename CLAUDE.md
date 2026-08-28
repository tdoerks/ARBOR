# Claude Code — Session Notes

Last active: (updated end of session 2026-08-28)

## Environment
- WSL Ubuntu-24.04 on Windows (KSU — tylerdoe)
- Claude Code runs inside Docker: `ghcr.io/tdoerks/claude-nextflow:latest`
- Launch: `cd ~/claude-workspace && docker run -it --rm -v "$PWD:/workspace" ghcr.io/tdoerks/claude-nextflow:latest`
- `/workspace` = `~/claude-workspace` on WSL host — only mounted dir, persists across restarts
- ARBOR repo: `/workspace/ARBOR`, COMPASS repo: `/workspace/COMPASS-pipeline`
- Beocat: `tylerdoe@icr-helios`, scratch at `/fastscratch/tylerdoe/`, bulk at `/bulk/tylerdoe/`
- NARMS bulk storage: `/bulk/tylerdoe/NARMS/` (uppercase)

## GitHub Auth (do once per container session)
```bash
gh auth login       # device flow → https://github.com/login/device
gh auth setup-git   # wires gh credentials to git
```

## Current Work

### ARBOR (github.com/tdoerks/ARBOR) — `main` branch
- **Pipeline fully working** — job 10854668 (MultiQC + dashboard rerun after SPAdes version fix)
  - SPAdes version string fix: `grep -oE '[0-9]+\.[0-9]+\.[0-9]+'` instead of raw `spades.py --version 2>&1` (commit 5aab167)
  - LOFREQ_CALL: 16GB × attempt, 4h × attempt, maxRetries=2 (commits af67bd5, 014f01e)
  - Assembly skips pooled samples: `!meta.id.endsWith('_pooled')` (commit 5fd4455)
  - Dashboard built ✔ — find with: `find /fastscratch/tylerdoe/ARBOR/tests/beocat/results_remap_A03ref -name "*.html" | grep -i dashboard`
- **In silico pooling**: D3/D7/D14 pooled samples for LoFreq variant frequency estimation
  - Submit dir: `/fastscratch/tylerdoe/ARBOR/tests/beocat`
  - Submit: `cd /fastscratch/tylerdoe/ARBOR/tests/beocat && sbatch run_arbor.sbatch`
- **Results rsync to bulk** — check if completed
  - `rsync -av /fastscratch/tylerdoe/ARBOR/tests/beocat/results_remap_A03ref/ /bulk/tylerdoe/ARBOR/results_remap_A03ref/`
- **Next**: Confirm dashboard HTML opens correctly; set up Globus share on /bulk/tylerdoe/ARBOR/

### COMPASS (github.com/tdoerks/COMPASS-pipeline) — `scratch` branch
- **Location on Beocat**: `/fastscratch/tylerdoe/COMPASS-1.1.0/`
- **Job 10765626** E. coli 6k run — status unknown (may have finished by now)
  - Watch: `tail -f compass*stdout*10765626*` from `/fastscratch/tylerdoe/COMPASS-1.1.0/`
  - Results: `/fastscratch/tylerdoe/COMPASS-1.1.0/results_mic_ecoli_6k/`
- **Viewer**: `bin/build_compass_viewer.py` — joins `compass_summary.tsv` + NCBI MIC/SIR metadata → self-contained HTML
  - Run: `git pull origin scratch && python3 bin/build_compass_viewer.py`
  - Copy to Windows: `scp tylerdoe@icr-helios:/fastscratch/tylerdoe/COMPASS-1.1.0/compass_mic_viewer.html /mnt/c/Users/tdoerks/Downloads/`
- **Shiga toxin (stx)** added (commit e723e6d): classifies isolates stx1/stx2/stx1+stx2/negative
- **Next**: Check job status → rebuild viewer if finished → review stx breakdown

### NARMS bulk storage cleanup — `/bulk/tylerdoe/NARMS/`
- **Goal**: Flatten FASTQs by year into `samples_clean/` folders, remove nested BaseSpace hash dirs
- **Structure**: `/bulk/tylerdoe/NARMS/2025/samples_clean/`, `/bulk/tylerdoe/NARMS/2026/`
- **In progress**: Moving 2025-collected samples (25KS prefix) from `2026/1-6-26_NARMS_WGS/` → `2025/samples_clean/`
  - FASTQs extracted from hash subdirs to `1-6-26_NARMS_WGS/` folder already
  - Move command: `mv /bulk/tylerdoe/NARMS/2026/1-6-26_NARMS_WGS/*.fastq.gz /bulk/tylerdoe/NARMS/2025/samples_clean/`
- **Duplicate issue found**: 26 samples appear in both top-level run folders AND `Runs/` subfolder
  - `Runs/` has 16 UNIQUE samples (March 3 requeue: `26KS0*_3_3_2026_*`) not elsewhere
  - Duplicates: top-level `2-23-26_NARMS_BWGS_2026-02-23/` and `3-5-26_NARMS_USDA_BWGS_EM_2026-03-05/` overlap with `Runs/`
  - Strategy: keep top-level organized versions; move unique Runs/ samples out before deleting Runs/
- **Next**: Confirm 1-6-26 FASTQs moved to samples_clean; continue with remaining 2026 run folders

### E-ink display (ordered 2026-08-24, arrives ~Aug 25-26)
- **Elecrow CrowPanel 2.13" ESP32-S3 e-ink** — ordered from Amazon ~$25
- Plan: Arduino sketch polls `https://raw.githubusercontent.com/tdoerks/open-source-lab-resources/main/checkin-board/status.json` every 10 min
- Displays: current status + note + QR code linking to web checkin-board
- Plugs in via USB-C (no built-in battery); 3D print case from Hale Library makerspace
- **Next**: Write Arduino sketch + flash via Arduino IDE

### BaseSpace download (pending)
- Run: `7-17-26_NARMS-NGS_BWGS_2026-07-17T21_35_42_9845731` (BCLConvert, 2GB)
- CLI: `C:\Users\tdoerks\Downloads\BaseSpace\bs.exe` — run from WSL terminal (not Docker)

## Key Patterns
- COMPASS log: `tail -f compass*stdout*<JOBID>*` (NOT `slurm-<JOBID>.out`)
- ARBOR log: `tail -f arbor_head_<JOBID>.log`
- Beocat jobs: `squeue -u tylerdoe`
- Lock file stuck: check `squeue`, cancel zombie job, then `rm -f .nextflow/cache/.../db/LOCK`
- SLURM preemption: just resubmit with `sbatch run_arbor.sbatch` (-resume handles it)
- Cancel duplicate job if submitted twice: `scancel <JOBID>`
- BaseSpace dirs have read-only permissions — `chmod u+w <dir>` before mv/rm

## How to Use This File
- At the start of a session, read this file to get context
- To save progress: "update CLAUDE.md with where we left off"
