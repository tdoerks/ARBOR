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
import glob
import gzip
import json
import os
import sys

DEFAULT_RESULTS = "/fastscratch/tylerdoe/ARBOR/tests/beocat/results_remap_A03ref"
DEFAULT_TEMPLATE = "/fastscratch/tylerdoe/ARBOR/dashboard/arbor_dashboard.html"
DEFAULT_OUTPUT = "/fastscratch/tylerdoe/ARBOR/tests/beocat/arbor_dashboard_A03_loaded.html"

# glob patterns (relative to results dir) the dashboard's parsers understand
PATTERNS = [
    "ivar/*.tsv",
    "lofreq/*.vcf.gz",
    "lofreq/*.vcf",
    "samtools/*.stats",
    "bowtie2/*.bowtie2.log",
    "mosdepth/*.summary.txt",
    "find/*.consensus.fasta",
    "iqtree/*.treefile",
    "pipeline_info/*params*.json",
]

MARKER = "__ARBOR_EMBED__"


def collect(results_dir):
    """Return [{name, b64}] for every matched result file (gz decompressed)."""
    items = []
    seen = set()
    for pat in PATTERNS:
        for path in sorted(glob.glob(os.path.join(results_dir, pat))):
            base = os.path.basename(path)
            if base in seen:
                continue
            with open(path, "rb") as fh:
                raw = fh.read()
            if base.endswith(".gz"):
                try:
                    raw = gzip.decompress(raw)
                    base = base[:-3]            # drop .gz so the parser takes the plain path
                except OSError:
                    pass                         # not actually gzipped; embed as-is
            if base in seen:
                continue
            seen.add(base)
            items.append({"name": base, "b64": base64.b64encode(raw).decode("ascii")})
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
