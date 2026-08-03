#!/usr/bin/env python3
"""
build_dashboard.py

Produce a PRE-LOADED ARBOR dashboard: take the standalone drop-in dashboard
(dashboard/arbor_dashboard.html) and a completed run's results dir, embed the
result files into a self-contained copy of the HTML, and auto-feed them to the
dashboard's own handleFiles() on page load. You open the output and it's already
populated -- no dragging -- with every tab the template supports.

This revives the old "ARBOR_DASHBOARD process" experience as a standalone script,
against the current template (so you keep the Parameters/Mapping/Consensus/SNP
tabs). It reuses the dashboard's existing JS parsers untouched: the embedded files
are reconstructed as File objects in the browser and passed to handleFiles(),
exactly as if you had dropped them.

Note: the template loads d3 + pako from a CDN, so rendering still needs internet
(pre-existing). We gunzip *.vcf.gz during embedding so the DATA doesn't need pako.

USAGE (defaults target the Beocat remap run):
    python3 bin/build_dashboard.py
    python3 bin/build_dashboard.py <results_dir> <template.html> <output.html>
"""

import base64
import gzip
import json
import os
import sys

DEFAULT_RESULTS = "/fastscratch/tylerdoe/ARBOR/tests/beocat/results_remap_A03ref"
DEFAULT_TEMPLATE = "/fastscratch/tylerdoe/ARBOR/dashboard/arbor_dashboard.html"
DEFAULT_OUTPUT = "/fastscratch/tylerdoe/ARBOR/tests/beocat/arbor_dashboard_A03_loaded.html"

MARKER = "__ARBOR_EMBED__"


def wanted(name):
    """True if this filename is one the dashboard's parsers consume.

    Curated (not every .log / .fasta) so a recursive walk of a full results tree
    picks up exactly the dashboard inputs and skips noise -- e.g. per-sample iVar
    *.fa full-genome consensus, the reference, FastQC HTML, BAMs. Works the same on
    a results dir (subdirs) and on Nextflow's flat-staged inputs.

    Includes *.per-base.bed.gz for the per-position coverage depth plot (Coverage
    tab). These are larger than the summary files (~1-5 MB each compressed) but
    are decompressed during embedding so the dashboard needs no pako for them.
    """
    n = name.lower()
    return (
        n.endswith(".tsv")                              # iVar variants
        or n.endswith(".vcf") or n.endswith(".vcf.gz") # LoFreq
        or n.endswith(".stats")                         # samtools stats (raw + primertrim)
        or n.endswith(".bowtie2.log")                   # bowtie2 alignment log
        or n.endswith(".mosdepth.summary.txt")          # mosdepth coverage summary (mean per segment)
        or n.endswith(".per-base.bed.gz")               # mosdepth per-position depth (coverage plot)
        or n.endswith(".per-base.bed")                  # mosdepth per-position depth (uncompressed)
        or n.endswith(".consensus.fasta")               # per-segment concatenated consensus
        or n.endswith(".treefile")                      # IQ-TREE
        or ("params" in n and n.endswith(".json"))
    )


def trim_consensus_ns(raw_bytes):
    """Strip flanking N-runs from every record in a FASTA byte string.

    iVar writes per-segment consensus against the full multi-contig reference, so
    each record is full-genome length (~11,979 bp for RVFV) with the real segment
    bases in the middle and N-padding on either side (e.g. NC_014395_S has ~10,289
    leading N's + 1,690 real bases). Trimming the flanks leaves each record at its
    natural segment length. Internal N's survive because .strip('Nn') only removes
    N's at the very start and end of the sequence.
    """
    out = []
    header = None
    seq_lines = []
    for line in raw_bytes.decode("utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if header is not None:
                seq = "".join(seq_lines).strip("Nn")
                out.append(header)
                out.append(seq)
            header = line
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        seq = "".join(seq_lines).strip("Nn")
        out.append(header)
        out.append(seq)
    return "\n".join(out).encode("utf-8")


def collect(results_dir):
    """Return [{name, b64}] for every dashboard input under results_dir.

    Walks recursively and classifies by filename, so it handles both a results
    tree (ivar/, lofreq/, ...) and a flat directory of Nextflow-staged files.
    *.gz inputs are decompressed and renamed (drop .gz) so the browser needs no
    pako to read the embedded data.
    """
    items = []
    seen = set()
    for root, _dirs, files in os.walk(results_dir):
        for fname in sorted(files):
            if not wanted(fname) or fname in seen:
                continue
            with open(os.path.join(root, fname), "rb") as fh:
                raw = fh.read()
            base = fname
            if base.endswith(".gz"):
                try:
                    raw = gzip.decompress(raw)
                    base = base[:-3]            # drop .gz so the parser takes the plain path
                except OSError:
                    pass                         # not actually gzipped; embed as-is
            if base in seen:
                continue
            # Strip flanking N-runs from consensus FASTA records. iVar writes per-segment
            # consensus against the full multi-contig reference, so each record is full-
            # genome length (~11,979 bp) with the real segment bases in the middle and
            # structural N-padding on either side. Trimming the flanking N's leaves each
            # record at its natural segment length with genuine coverage-gap N's intact
            # (they are flanked by real bases, so .strip('Nn') leaves them alone).
            if base.endswith(".consensus.fasta"):
                raw = trim_consensus_ns(raw)
            seen.add(base)
            items.append({"name": base, "b64": base64.b64encode(raw).decode("ascii")})
    items.sort(key=lambda it: it["name"])
    return items


def build(results_dir, template, output):
    if not os.path.isfile(template):
        sys.exit("Template not found: %s" % template)
    if not os.path.isdir(results_dir):
        sys.exit("Results dir not found: %s" % results_dir)

    items = collect(results_dir)
    if not items:
        sys.exit("No result files matched under %s -- nothing to embed." % results_dir)

    payload = json.dumps(items, separators=(",", ":"))
    # base64 content contains no '<', so it can't break out of the <script> tag.
    embed_block = (
        '<script id="%s" type="application/json">%s</script>\n'
        "<script>\n"
        "window.addEventListener('load',function(){\n"
        "  var el=document.getElementById('%s'); if(!el) return;\n"
        "  var items; try{items=JSON.parse(el.textContent);}catch(e){console.error('embed parse',e);return;}\n"
        "  var files=items.map(function(it){\n"
        "    var bin=atob(it.b64), arr=new Uint8Array(bin.length);\n"
        "    for(var i=0;i<bin.length;i++) arr[i]=bin.charCodeAt(i);\n"
        "    return new File([arr], it.name);\n"
        "  });\n"
        "  if(typeof handleFiles==='function'){ handleFiles(files); }\n"
        "  else{ console.error('handleFiles not defined -- template changed?'); }\n"
        "});\n"
        "</script>\n"
    ) % (MARKER, payload, MARKER)

    with open(template) as fh:
        html = fh.read()
    if "</body>" not in html:
        sys.exit("Template has no </body> -- cannot inject.")
    html = html.replace("</body>", embed_block + "</body>", 1)

    with open(output, "w") as fh:
        fh.write(html)

    size_mb = os.path.getsize(output) / 1e6
    print("Embedded %d files -> %s (%.1f MB)" % (len(items), output, size_mb))
    print("Open it in a browser (needs internet for d3/pako CDN); it self-loads.")


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_RESULTS
    template = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_TEMPLATE
    output = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_OUTPUT
    build(results_dir, template, output)


if __name__ == "__main__":
    main()
