#!/usr/bin/env python3
"""
make_samplesheet.py

Build an ARBOR samplesheet (sample,fastq_1,fastq_2) by discovering paired-end
fastqs across the RVFV read drops. Two directory layouts are handled:

  "folder" layout (reads_2026-07-02): one BaseSpace subdir per sample, named
      <SAMPLE>-ds.<hash>/  containing <SAMPLE>_S##_L001_R{1,2}_001.fastq.gz
      (A03 is here too, one extra level deep: A03-ds.<hash>/A03-ds.<hash>/...).
      Sample name = the folder prefix before "-ds.".

  "flat" layout (reads_2026-04-24): fastqs sit directly in the dir, named like
      S1-RVFV-MP-12-R3D701-repeat-4-22-2026_S4_L001_R1_001.fastq.gz
      Sample name = the RxDyyy token pulled from the filename.

Non-sample entries (BCL_Convert_*, Undetermined*, *.csv) are skipped. R1 files
are matched to their R2 partner; single-end / unpaired files are reported and
skipped (ARBOR here is run paired-end). Duplicate sample names abort the run.

USAGE (defaults target Beocat):
    python3 bin/make_samplesheet.py                    # writes tests/beocat/samplesheet.csv
    python3 bin/make_samplesheet.py <out.csv>          # custom output path
"""

import glob
import os
import re
import sys

# ---- Sources (directory, layout) --------------------------------------------
BASE = "/fastscratch/tylerdoe"
SOURCES = [
    (os.path.join(BASE, "reads_2026-07-02"), "folder"),
    (os.path.join(BASE, "reads_2026-04-24"), "flat"),
]
DEFAULT_OUT = "/fastscratch/tylerdoe/ARBOR/tests/beocat/samplesheet.csv"

EXCLUDE_PREFIXES = ("BCL_Convert", "Undetermined")
RXDY = re.compile(r"(R\d+D\d+[A-Z]?)")   # e.g. R3D701, R1D312F


def r2_for(r1):
    """Map an R1 fastq path to its expected R2 partner path."""
    return r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")


def is_excluded(name):
    return any(name.startswith(p) for p in EXCLUDE_PREFIXES)


def discover_folder(root):
    """One sample per <SAMPLE>-ds.<hash> subdir; recurse for the fastqs."""
    found = []
    for entry in sorted(os.listdir(root)):
        path = os.path.join(root, entry)
        if not os.path.isdir(path) or is_excluded(entry) or "-ds." not in entry:
            continue
        sample = entry.split("-ds.")[0]
        r1s = sorted(glob.glob(os.path.join(path, "**", "*_R1_001.fastq.gz"),
                               recursive=True))
        if not r1s:
            print("  WARN: no R1 fastq under %s -- skipped" % entry)
            continue
        if len(r1s) > 1:
            print("  WARN: %d R1 files under %s -- using %s"
                  % (len(r1s), entry, os.path.basename(r1s[0])))
        found.append((sample, r1s[0]))
    return found


def discover_flat(root):
    """Fastqs directly in the dir; sample name from the RxDyyy token."""
    found = []
    for r1 in sorted(glob.glob(os.path.join(root, "*_R1_001.fastq.gz"))):
        base = os.path.basename(r1)
        if is_excluded(base):
            continue
        m = RXDY.search(base)
        if not m:
            print("  WARN: no RxDyyy token in %s -- skipped" % base)
            continue
        found.append((m.group(1), r1))
    return found


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT

    rows = []          # (sample, r1, r2)
    seen = {}          # sample -> source (collision detection)
    for root, layout in SOURCES:
        if not os.path.isdir(root):
            sys.exit("Source directory not found: %s" % root)
        print("Scanning %s (%s layout)" % (root, layout))
        found = discover_folder(root) if layout == "folder" else discover_flat(root)

        n_ok = 0
        for sample, r1 in found:
            r2 = r2_for(r1)
            if not os.path.exists(r2):
                print("  WARN: %s has no R2 partner (%s) -- skipped"
                      % (sample, os.path.basename(r2)))
                continue
            if sample in seen:
                sys.exit("Duplicate sample name %r (in %s and %s); aborting."
                         % (sample, seen[sample], root))
            seen[sample] = root
            rows.append((sample, r1, r2))
            n_ok += 1
        print("  -> %d paired samples" % n_ok)

    rows.sort(key=lambda r: r[0])

    with open(out_path, "w") as out:
        out.write("sample,fastq_1,fastq_2\n")
        for sample, r1, r2 in rows:
            out.write("%s,%s,%s\n" % (sample, r1, r2))

    print("\nWrote %d samples to %s" % (len(rows), out_path))
    if "A03" in seen:
        print("A03 included (from %s)" % seen["A03"])
    else:
        print("NOTE: A03 was NOT found -- check its path.")


if __name__ == "__main__":
    main()
