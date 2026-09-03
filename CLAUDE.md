# Claude Code — Session Notes

Last active: (updated end of session 2026-09-03)

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
- **Structure**: `/bulk/tylerdoe/NARMS/2025/samples_clean/`, `/bulk/tylerdoe/NARMS/2026/` (flat)
- **Pass sample list**: 112 confirmed-pass 26KS samples (26KS01–26KS07) to move to `/bulk/tylerdoe/NARMS/2026/`
  - Script at `/workspace/copy_pass_samples.py` — finds each sample anywhere under NARMS and copies R1+R2 to flat 2026 folder
  - Get to Beocat: `wsl bash -c "cat ~/claude-workspace/copy_pass_samples.py | ssh tylerdoe@beocat.ksu.edu 'cat > /fastscratch/tylerdoe/copy_pass_samples.py'"`
  - Run on Beocat: `python3 /fastscratch/tylerdoe/copy_pass_samples.py`
- **BaseSpace permission issue**: dirs downloaded with `dr-x` — always `chmod -R u+w <folder>` before mv/rm
- **Progress so far**: 26KS01 + 26KS02 samples confirmed present in 2026 flat folder; 26KS03–07 still in unprocessed subfolders
- **Key finding**: Sxx slot numbers differ between pass list and disk filenames — that's normal, doesn't affect identity
- **One lib ID discrepancy**: pass list has 26KS02CL01-EC as `26DL039` but disk has `26DL038` — likely typo in pass list
- **Remaining run folders to process**: `2-23-26_NARMS_BWGS_2026-02-23`, `3_3_2026_NARMS_WGS_requeue_NEB`, `3-5-26_NARMS_USDA_BWGS_EM_2026-03-05`, `4-22-26_RVFV_req_NEB`, `ICA_Workflows_2026_04`, `NARMS-425401161`, `narms-435220242`, `Runs/` (has 16 unique March 3 requeue samples — don't delete without extracting those first)
- **Next**: Run copy_pass_samples.py on Beocat to bulk-copy all pass samples to flat 2026 folder

### Platinum-Calibration (github.com/tdoerks/Platinum-Calibration)
- **Live site**: https://tdoerks.github.io/Platinum-Calibration/ (serves from `main`)
- **Newest branch**: `unified-rebuild` (July 13) — adds DYMO label printer tab, 1 commit ahead of main, 5 behind
- **Google Drive integration**: `google-drive-integration` branch — Phase 1 auth started (May 3) but never merged
- **Goal**: Add Microsoft backend (instead of Google Drive) + pipette serial number registry with auto-fill
  - When serial # typed → XLOOKUP auto-fills manufacturer/model/max volume from registry
  - New pipettes: fill manually → add to registry
- **Data storage decision**: Local dedicated PC (SQLite) rather than cloud — data stays in building
- **Next**: Finish Google Drive integration branch OR rewrite for Microsoft Graph API

### Lab Server / Open-Source Lab Tools
- **Repo**: github.com/tdoerks/open-source-lab-resources
- **Tools**: bacterial-isolation-tracker, checkin-board, freezer-inventory, pipetting-tutorial, miseq-pooling, etc.
- **Goal**: Host tools on a dedicated lab PC, accessible to team from anywhere, client data secure
- **Architecture**:
  ```
  Internet → Cloudflare Access (login gate, free up to 50 users)
                  ↓
            Cloudflare Tunnel (outbound, works through guest WiFi/NAT)
                  ↓
            Dedicated lab PC (nginx serves HTML, SQLite stores data)
  ```
- **Tested 2026-09-02**: Cloudflare quick tunnel working on current PC
  - Temp URL was: `https://shapes-prizes-advertiser-atom.trycloudflare.com` (expired)
  - Test command (WSL): `python3 -m http.server 8080` + `cloudflared tunnel --url http://localhost:8080`
- **Guest WiFi OK**: Cloudflare Tunnel is outbound-only, works through client isolation
- **Data security plan**: Data on local PC (not cloud), Cloudflare Access gates login by email
- **Auth plan**: Personal Microsoft 365 or Google account (NOT KSU — KSU IT controls those)
- **Next**: Set up dedicated PC with Ubuntu + nginx + Cloudflare named tunnel + Access login wall

### E-ink display (Elecrow CrowPanel 2.13" ESP32-S3)
- **Sketch written**: `/workspace/eink_checkin/checkin_display.ino` — complete and ready to flash
  - WiFi → polls status.json every 10 min → deep sleep between refreshes
  - Renders: status lines (left) + QR code (right) + timestamp (bottom)
  - Libraries needed: GxEPD2, ArduinoJson, qrcode (Richard Moore)
  - Board setting: ESP32S3 Dev Module
  - Pins: EPD_CS=5, EPD_DC=17, EPD_RST=16, EPD_BUSY=4 (verify against Elecrow datasheet)
- **Flashing status: BLOCKED — USB driver issue**
  - Board plugged in via USB-C → shows as "billboard device" in Windows Device Manager (no COM port assigned)
  - ESP32-S3 uses native USB — needs ESP32-S3 USB CDC driver, user has no Windows admin rights
  - WSL only sees /dev/ttyS0-S7, no /dev/ttyUSB* or /dev/ttyACM*
  - esptool + python3-serial installed in WSL already: `sudo apt install esptool python3-serial -y`
- **Next steps to try**:
  1. Hold BOOT button on board, plug in USB, release BOOT after 2s → check Device Manager for new device
  2. Check if usbipd-win already installed in PowerShell: `usbipd list`
     - If yes: `usbipd bind --busid <id>` → `usbipd attach --wsl` → flash via WSL
  3. If neither works: try a different computer with admin rights to install ESP32 drivers + Arduino IDE

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
