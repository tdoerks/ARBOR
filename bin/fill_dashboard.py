#!/usr/bin/env python3
"""Pre-fill the ARBOR dashboard by embedding pipeline result files into the HTML.

Reads the ARBOR result files (iVar .tsv, LoFreq .vcf.gz, mosdepth .summary.txt /
.per-base.bed.gz, IQ-TREE .treefile, fastp .json, samtools .stats, metadata .csv),
base64-encodes each, and injects them into the dashboard template's
``__ARBOR_EMBEDDED_DATA__`` placeholder so the published HTML opens already populated.

Stdlib only -- runs in any python container and standalone on a results/ folder:

    python fill_dashboard.py --template dashboard/arbor_dashboard.html \\
        --metadata arbor_metadata.csv -o arbor_dashboard.html results/

The same script is invoked by the ARBOR_DASHBOARD Nextflow process.
"""
import argparse
import base64
import json
import os
import sys

PLACEHOLDER = "__ARBOR_EMBEDDED_DATA__"


def recognized(name):
    """True if a filename is a dashboard-relevant result file."""
    n = name.lower()
    return (
        n.endswith(".tsv")
        or n.endswith(".vcf")
        or n.endswith(".vcf.gz")
        or (n.endswith(".txt") and "summary" in n)
        or ("per-base" in n and n.endswith(".bed.gz"))
        or n.endswith(".treefile")
        or n.endswith(".nwk")
        or n.endswith(".json")
        or n.endswith(".stats")
        or n.endswith(".csv")
    )


def collect(paths):
    """Expand the given paths (files or directories) into a list of result files."""
    files = []
    for p in paths:
        if os.path.isdir(p):
            for root, _dirs, names in os.walk(p):
                for nm in names:
                    if recognized(nm):
                        files.append(os.path.join(root, nm))
        elif os.path.isfile(p):
            if recognized(os.path.basename(p)):
                files.append(p)
            else:
                print(f"[fill_dashboard] skipping unrecognized file: {p}", file=sys.stderr)
        else:
            print(f"[fill_dashboard] not found: {p}", file=sys.stderr)
    return files


def encode(path):
    """Return a {name, b64, gz} entry for one result file."""
    name = os.path.basename(path)
    with open(path, "rb") as fh:
        raw = fh.read()
    return {"name": name, "b64": base64.b64encode(raw).decode("ascii"), "gz": name.lower().endswith(".gz")}


def main():
    ap = argparse.ArgumentParser(description="Embed ARBOR results into the dashboard HTML.")
    ap.add_argument("--template", required=True, help="dashboard HTML template (with the placeholder)")
    ap.add_argument("--metadata", help="optional metadata CSV (sample,host,date,location ...)")
    ap.add_argument("-o", "--out", default="arbor_dashboard.html", help="output HTML path")
    ap.add_argument("files", nargs="*", help="result files or directories (default: scan cwd)")
    args = ap.parse_args()

    paths = list(args.files) or ["."]
    files = collect(paths)
    if args.metadata and os.path.isfile(args.metadata):
        files.append(args.metadata)

    # De-duplicate by basename (channels may stage a file more than once).
    seen, entries = set(), []
    for f in sorted(files):
        bn = os.path.basename(f)
        if bn in seen:
            continue
        seen.add(bn)
        entries.append(encode(f))

    if not entries:
        print("[fill_dashboard] WARNING: no recognized result files found to embed", file=sys.stderr)

    with open(args.template, "r", encoding="utf-8") as fh:
        html = fh.read()
    if PLACEHOLDER not in html:
        sys.exit(f"[fill_dashboard] template has no {PLACEHOLDER} placeholder: {args.template}")

    payload = json.dumps(entries, separators=(",", ":"))
    if "</script" in payload.lower():  # guard the inline <script> block; should never happen
        payload = payload.replace("</", "<\\/")
    html = html.replace(PLACEHOLDER, payload)

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"[fill_dashboard] embedded {len(entries)} files -> {args.out}")


if __name__ == "__main__":
    main()
