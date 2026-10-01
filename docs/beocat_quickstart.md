# ARBOR on Beocat — command-by-command quickstart

Every command you need, in order, for running ARBOR under **your own** Beocat account. Nothing needs
renaming: wherever a path says `$USER`, Beocat fills in your username automatically. The only things you
type yourself are your **username** (step 3, when copying files from your computer) and a **run name**.

For background and options, see the full [Beocat User Guide](beocat_user_guide.md).

> **Paste one command at a time.** Copy one grey box, paste it into the Beocat terminal, press Enter,
> and wait for the prompt (`[you@icr-helios ...]$`) to come back before the next one. Pasting several
> lines at once can glue them together into one broken command.

---

## 0. Log in

From your computer (Windows PowerShell, Mac Terminal, or WSL Ubuntu), replacing `YOUR_EID` with your
K-State eID:

```bash
ssh YOUR_EID@beocat.ksu.edu
```

Check who you are — this is the `$USER` used in every path below:

```bash
echo $USER
```

---

## 1. Get the pipeline (first time only)

```bash
cd /fastscratch/$USER
```
```bash
git clone https://github.com/tdoerks/ARBOR.git
```

Already have it? Update it instead:

```bash
cd /fastscratch/$USER/ARBOR && git pull
```

---

## 2. Make a folder for this run

Pick a short run name with no spaces, e.g. `RVFV_2026-10`. Set it once per login session:

```bash
RUN=RVFV_2026-10
```
```bash
mkdir -p /fastscratch/$USER/ARBOR/runs/$RUN/reads
```

> Runs must live in `ARBOR/runs/<name>/` — the run script finds the pipeline two folders up.

---

## 3. Upload your reads (from your computer)

Download the run's FASTQ files from BaseSpace to your computer first. Then, **in a terminal on your
computer** (not Beocat), copy them up. Replace `YOUR_EID`, the run name, and the folder path:

**WSL Ubuntu / Mac:**
```bash
scp /path/to/fastq_folder/*.fastq.gz YOUR_EID@beocat.ksu.edu:/fastscratch/YOUR_EID/ARBOR/runs/RVFV_2026-10/reads/
```

**Windows PowerShell:**
```powershell
scp "C:\path\to\fastq_folder\*.fastq.gz" YOUR_EID@beocat.ksu.edu:/fastscratch/YOUR_EID/ARBOR/runs/RVFV_2026-10/reads/
```

If BaseSpace put each sample in its own subfolder, use `fastq_folder/*/*.fastq.gz` instead.
`Undetermined_*` files can come along — they are skipped automatically.

Back **on Beocat**, check they arrived (you should see `_R1_001` and `_R2_001` files):

```bash
ls /fastscratch/$USER/ARBOR/runs/$RUN/reads | head
```

---

## 4. Build the samplesheet

```bash
cd /fastscratch/$USER/ARBOR/runs/$RUN
```
```bash
bash ../../tests/beocat/make_samplesheet.sh $PWD/reads > samplesheet.csv
```
```bash
cat samplesheet.csv
```

Check the output:

- one row per sample, each with an R1 and an R2 path
- sample names come from the file names (`R3D720_S5_L001_R1_001.fastq.gz` → `R3D720`)
- **no spaces** in sample names and **no duplicates** — this prints nothing if names are unique:

```bash
cut -d, -f1 samplesheet.csv | sort | uniq -d
```

---

## 5. Get the run script

```bash
cp ../../tests/beocat/run_arbor.sbatch .
```

**Optional — email when the run finishes.** Open the script:

```bash
nano run_arbor.sbatch
```

Find the two lines starting `##SBATCH --mail-`, delete **one** `#` from each, and put your own address
after `--mail-user=`. Save with **Ctrl+O**, Enter, then exit with **Ctrl+X**.

---

## 6. Start the run

```bash
sbatch run_arbor.sbatch
```

Beocat replies `Submitted batch job 12345678` — note that number (your **job ID**).

Defaults: RVFV MP-12 reference, bundled primer scheme, S/M/L segments. Add options after the script
name if needed, e.g. a different reference:

```bash
sbatch run_arbor.sbatch --reference /path/to/reference.fa
```

---

## 7. Watch it

Your jobs (the head job named `arbor`, plus one job per step it launches):

```bash
squeue -u $USER
```

Live progress (replace the number with your job ID; **Ctrl+C** stops watching, not the run):

```bash
tail -f arbor_head_12345678.log
```

The first run in a new folder is quiet for a few minutes while it installs Nextflow and downloads
containers. It's finished when the log says `Pipeline completed successfully`.

---

## 8. Get the results

Everything is in `results/` inside your run folder. The two reports to open first:

```bash
ls results/dashboard/arbor_dashboard_loaded.html results/multiqc/multiqc_report.html
```

Copy them **to your computer** (run this on your computer, with your eID and run name):

```bash
scp YOUR_EID@beocat.ksu.edu:/fastscratch/YOUR_EID/ARBOR/runs/RVFV_2026-10/results/dashboard/arbor_dashboard_loaded.html .
```

Then double-click the HTML file to open it in a browser.

**Save the results somewhere permanent** — `/fastscratch` is purged after a period of inactivity. If you
have bulk storage (`/bulk/$USER`), on Beocat:

```bash
mkdir -p /bulk/$USER/ARBOR/$RUN && rsync -av results/ /bulk/$USER/ARBOR/$RUN/
```

---

## Fixing things

**Re-run after a failure or after `git pull`** — from the run folder, just submit again; finished steps
are reused:

```bash
sbatch run_arbor.sbatch
```

**Cancel a run** — cancel the head job, then check for leftover step jobs (cancelling the head does not
always stop them) and cancel those too:

```bash
scancel 12345678
```
```bash
squeue -u $USER
```

**"Unable to acquire lock"** — an earlier run in this folder is still active (check `squeue`) or was
killed uncleanly. Make sure nothing is running, then:

```bash
rm -f .nextflow/cache/*/db/LOCK
```

**A step failed** — the log names the step and its work folder; the error is in that folder's
`.command.err`:

```bash
tail -n 20 work/ab/cdef12…/.command.err
```

---

## Useful options

| Option | What it does |
|---|---|
| `--reference ref.fa` | Map to a different reference (default: RVFV MP-12) |
| `--skip_pooling false` | Turn on in silico D3/D7/D14 day pools (off by default) |
| `--pool_exclude '_rerun$'` | With pooling on, leave matching samples out of the pools |
| `--skip_lofreq` / `--skip_phylogeny` | Skip LoFreq variant calling / the trees |
| `--context_fasta strains.fa` | Add external strains to the trees |

Pass any of these after `run_arbor.sbatch`, e.g. `sbatch run_arbor.sbatch --skip_phylogeny`.
