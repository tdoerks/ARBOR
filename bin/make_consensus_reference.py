#!/usr/bin/env python3
"""
make_consensus_reference.py

Turn one sample's ARBOR consensus into a gap-free reference for a new ARBOR run
(e.g. an in-house stock consensus, so the stock's own fixed differences stop showing up as
"variants"). Generalised from the A03 remap script (patch_consensus.py).

Per segment:
  1. Slice the segment out of ARBOR's per-segment iVar consensus. Each results/ivar/<sample>_<seg>.fa
     holds the WHOLE concatenated genome with real bases only in that segment's slice; the bundled
     reference's record order/lengths define the slices (last contig to EOF: iVar may add a base).
  2. MAFFT-align the fill reference's segment (default: true MP-12, assets/rvfv/mp12_reference.fa)
     to the consensus slice.
  3. Walk the alignment in fill-reference coordinates:
       - fill reference has a gap       -> drop (consensus insertion; keeps the coordinate frame)
       - consensus base is real          -> keep it (the sample's own calls / fixed differences)
       - consensus is N or gap           -> use the fill reference's base (ends, primer sites)
  4. Write it under the bundled segment ID (NC_014395_S / NC_014396_M / NC_014397_L), at the
     exact bundled length, so the primer BED and the default --segments keep working.

Usage (needs `mafft` on PATH — on Beocat: module load MAFFT):
    python3 bin/make_consensus_reference.py <sample> [--results results] [--out <sample>_reference.fa]
        [--frame-ref assets/rvfv/rvfv_reference.fa] [--fill-ref assets/rvfv/mp12_reference.fa]
Then run ARBOR with --reference <out>.

Also writes, next to <out>:
  <out base>.filled.bed   regions filled from the fill reference (no sample data) — load in IGV
  <out base>.diffs.tsv    positions where the sample consensus differs from the fill reference
  <out base>.map.html     downloadable map: per segment, where each base came from (consensus vs
                          filled), primer/amplicon tiling, differences, and a position lookup —
                          for checking whether a SNP position lies in real sequence or the stitch
"""
import argparse
import os
import subprocess
import sys
import html
import json
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "assets", "rvfv")

# bundled segment ID  ->  MP-12 record ID in mp12_reference.fa (segments S, M, L)
FILL_IDS = {"NC_014395_S": "DQ380154.1", "NC_014396_M": "DQ380208.1", "NC_014397_L": "DQ375404.1"}


def read_fasta(path=None, text=None):
    recs, name = [], None
    for line in (open(path) if path else text.splitlines()):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs.append([name, []])
        elif name:
            recs[-1][1].append(line.upper())
    return [(n, "".join(s)) for n, s in recs]


def mafft_pair(a_id, a_seq, b_id, b_seq):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
        t.write(f">{a_id}\n{a_seq}\n>{b_id}\n{b_seq}\n")
        path = t.name
    try:
        r = subprocess.run(["mafft", "--quiet", "--auto", path], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"MAFFT failed for {b_id}:\n{r.stderr}")
        return dict(read_fasta(text=r.stdout))
    finally:
        os.unlink(path)


FILL_LABEL = "MP12"


def to_runs(positions):
    runs = []
    for p in positions:
        if runs and p == runs[-1][1] + 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])
    return [tuple(r) for r in runs]


def read_primers(path):
    """BED -> {segment: [(start1, end1, name), ...]} in 1-based inclusive coordinates."""
    out = {}
    if path and os.path.exists(path):
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                out.setdefault(f[0], []).append((int(f[1]) + 1, int(f[2]), f[3]))
    return out


MAP_CSS = """
:root{--surface:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#8a8984;--line:#e4e3de;
--cons:#2a78d6;--fill:#eb6834;--amp:#b9b8b2;--primer:#52514e}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--surface:#1a1a19;--panel:#222221;--ink:#fff;
--ink2:#c3c2b7;--muted:#8f8e88;--line:#3a3a38;--cons:#3987e5;--fill:#d95926;--amp:#5a5955;--primer:#c3c2b7;color-scheme:dark}}
:root[data-theme="dark"]{--surface:#1a1a19;--panel:#222221;--ink:#fff;--ink2:#c3c2b7;--muted:#8f8e88;--line:#3a3a38;
--cons:#3987e5;--fill:#d95926;--amp:#5a5955;--primer:#c3c2b7;color-scheme:dark}
body{margin:0;padding:24px 16px;background:var(--surface);color:var(--ink);font:14px/1.5 system-ui,sans-serif}
main{max-width:1100px;margin:0 auto}h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 6px}
.sub{color:var(--ink2);margin:0 0 14px}.legend{display:flex;gap:18px;flex-wrap:wrap;color:var(--ink2);font-size:13px;margin:8px 0 4px}
.sw{display:inline-block;width:14px;height:10px;border-radius:2px;margin-right:6px;vertical-align:middle}
.card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px 14px;margin:10px 0}
svg{width:100%;height:auto;display:block}svg text{fill:var(--ink2);font:11px system-ui,sans-serif}
.stats{color:var(--ink2);font-size:13px}table{border-collapse:collapse;font-size:13px;width:100%}
th,td{text-align:left;padding:4px 8px;border-bottom:1px solid var(--line)}th{color:var(--ink2);font-weight:600}
.look{display:flex;gap:8px;flex-wrap:wrap;align-items:center}input,select{font:inherit;padding:4px 8px;border:1px solid var(--line);
border-radius:6px;background:var(--panel);color:var(--ink)}#lookres{font-weight:600}
.wrap{overflow-x:auto}
"""


def write_map(path, sample, segs, primers):
    """Self-contained HTML map: per segment, source of every base + amplicon tiling + differences."""
    W, X0, X1 = 1000, 70, 980
    parts = []
    for si in segs:
        L = si["length"]
        sx = lambda p: X0 + (p - 1) / max(1, L - 1) * (X1 - X0)
        svg = [f'<svg viewBox="0 0 {W} 168" role="img" aria-label="{si["seg"]} source map">']
        # track 1: base source — consensus (blue) with filled runs (orange) on top
        svg.append(f'<text x="0" y="34">source</text>')
        svg.append(f'<rect x="{X0}" y="22" width="{X1 - X0}" height="16" rx="3" fill="var(--cons)">'
                   f'<title>{html.escape(sample)} consensus (real sequenced bases)</title></rect>')
        for a, b in si["filled_runs"]:
            x, w = sx(a), max(3, sx(b) - sx(a))
            svg.append(f'<rect x="{x:.1f}" y="22" width="{w:.1f}" height="16" rx="2" fill="var(--fill)" '
                       f'stroke="var(--panel)" stroke-width="1"><title>{si["seg"]} {a}-{b}: filled from '
                       f'{si["fill_id"]} ({b - a + 1} bp, no sample data)</title></rect>')
        # track 2: amplicons (alternate two rows) + primers
        svg.append(f'<text x="0" y="74">amplicons</text>')
        amps = {}
        for s_, e_, name in primers.get(si["seg"], []):
            key = name.rsplit("_", 1)[0]
            amps.setdefault(key, []).append((s_, e_, name))
        for i, (key, pr) in enumerate(sorted(amps.items(), key=lambda kv: min(p[0] for p in kv[1]))):
            lo, hi = min(p[0] for p in pr), max(p[1] for p in pr)
            y = 58 + (i % 2) * 14
            svg.append(f'<line x1="{sx(lo):.1f}" x2="{sx(hi):.1f}" y1="{y + 5}" y2="{y + 5}" stroke="var(--amp)" '
                       f'stroke-width="2"><title>{key}: {lo}-{hi}</title></line>')
            for s_, e_, name in pr:
                svg.append(f'<rect x="{sx(s_):.1f}" y="{y + 1}" width="{max(2, sx(e_) - sx(s_)):.1f}" height="8" '
                           f'fill="var(--primer)"><title>{name}: {s_}-{e_}</title></rect>')
        # track 3: differences from the fill reference
        svg.append(f'<text x="0" y="110">vs {FILL_LABEL}</text>')
        label_ends = [-1e9, -1e9]          # right edge of the last label on each of two label rows
        for p_, r, c in si["diffs"]:
            x = sx(p_)
            label = f"{p_} {r}&#8594;{c}"
            width = 7 * (len(str(p_)) + 4)
            row = 0 if x >= label_ends[0] + 4 else 1 if x >= label_ends[1] + 4 else None
            svg.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="94" y2="130" stroke="var(--ink)" stroke-width="2"/>'
                       f'<rect x="{x - 6:.1f}" y="90" width="12" height="42" fill="transparent">'
                       f'<title>{si["seg"]} position {p_}: {FILL_LABEL} {r} &#8594; {html.escape(sample)} {c}</title></rect>')
            if row is not None:                # crowded beyond two rows: tooltip + table carry it
                svg.append(f'<text x="{x + 4:.1f}" y="{106 + 14 * row}" style="fill:var(--ink)">{label}</text>')
                label_ends[row] = x + 4 + width
        if not si["diffs"]:
            svg.append(f'<text x="{X0}" y="110">no differences</text>')
        # axis
        step = 500 if L <= 4000 else 1000
        for t in [1] + list(range(step, L, step)) + [L]:
            svg.append(f'<line x1="{sx(t):.1f}" x2="{sx(t):.1f}" y1="146" y2="150" stroke="var(--muted)"/>'
                       f'<text x="{sx(t):.1f}" y="164" text-anchor="middle">{t}</text>')
        svg.append(f'<line x1="{X0}" x2="{X1}" y1="146" y2="146" stroke="var(--line)"/></svg>')
        filled_txt = ", ".join(f"{a}–{b} ({b - a + 1} bp)" for a, b in si["filled_runs"]) or "none"
        parts.append(f'<h2>{si["seg"]} <span class="stats">· {L} bp · {si["kept"]} from consensus '
                     f'({100 * si["kept"] / L:.1f}%) · {si["filled"]} filled from {si["fill_id"]} · '
                     f'{len(si["diffs"])} difference(s)</span></h2><div class="card wrap">{"".join(svg)}</div>'
                     f'<p class="stats">Filled (stitch): {filled_txt}</p>')

    rows = "".join(f'<tr><td>{si["seg"]}</td><td>{a}</td><td>{b}</td><td>{b - a + 1}</td><td>filled from {si["fill_id"]}</td></tr>'
                   for si in segs for a, b in si["filled_runs"])
    drows = "".join(f'<tr><td>{si["seg"]}</td><td>{p_}</td><td>{r}</td><td>{c}</td></tr>'
                    for si in segs for p_, r, c in si["diffs"]) or '<tr><td colspan="4">none</td></tr>'
    lookup = {si["seg"]: {"length": si["length"], "filled": si["filled_runs"],
                          "diffs": {str(p_): [r, c] for p_, r, c in si["diffs"]},
                          "primers": [[s_, e_, n] for s_, e_, n in primers.get(si["seg"], [])]} for si in segs}
    opts = "".join(f'<option>{si["seg"]}</option>' for si in segs)
    doc = f"""<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(sample)} reference map</title><style>{MAP_CSS}</style><main>
<h1>{html.escape(sample)} in-house reference</h1>
<p class="sub">Where every base of this ARBOR reference came from. Blue = real {html.escape(sample)} consensus.
Orange = the "stitch": positions with no sample coverage, filled from {FILL_LABEL} so the reference has no gaps.
No SNP can be called in an orange region — reads do not reach it (it is outside the amplicon tiling).</p>
<div class="legend"><span><span class="sw" style="background:var(--cons)"></span>{html.escape(sample)} consensus</span>
<span><span class="sw" style="background:var(--fill)"></span>filled from {FILL_LABEL} (no data)</span>
<span><span class="sw" style="background:var(--primer)"></span>primer</span>
<span><span class="sw" style="background:var(--amp)"></span>amplicon</span>
<span>│ position differs from {FILL_LABEL}</span></div>
<div class="card"><div class="look"><b>Check a SNP position:</b> <select id="seg">{opts}</select>
<input id="pos" type="number" min="1" placeholder="position" style="width:120px">
<span id="lookres"></span></div></div>
{"".join(parts)}
<h2>Filled regions (stitch)</h2><div class="wrap"><table><tr><th>Segment</th><th>Start</th><th>End</th><th>Length</th><th>Source</th></tr>{rows}</table></div>
<h2>Differences from {FILL_LABEL}</h2><div class="wrap"><table><tr><th>Segment</th><th>Position</th><th>{FILL_LABEL}</th><th>{html.escape(sample)}</th></tr>{drows}</table></div>
</main><script>
const D={json.dumps(lookup)};
function look(){{const s=document.getElementById('seg').value,p=+document.getElementById('pos').value,o=document.getElementById('lookres'),d=D[s];
if(!p){{o.textContent='';return}} if(p<1||p>d.length){{o.textContent='outside '+s+' (1-'+d.length+')';return}}
const f=d.filled.some(r=>p>=r[0]&&p<=r[1]);const pr=d.primers.filter(x=>p>=x[0]&&p<=x[1]).map(x=>x[2]);
let t=f?'STITCH — filled from {FILL_LABEL}, no sample coverage (a SNP here is not real)':'Real consensus sequence';
if(d.diffs[p]) t+=' · this position differs from {FILL_LABEL} ('+d.diffs[p][0]+'→'+d.diffs[p][1]+')';
if(pr.length) t+=' · inside primer '+pr.join(', ')+' (primer-trimmed; lower confidence)';
o.textContent=t;}}
document.getElementById('pos').addEventListener('input',look);document.getElementById('seg').addEventListener('change',look);
</script>"""
    with open(path, "w") as f:
        f.write(doc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out")
    ap.add_argument("--frame-ref", default=os.path.join(ASSETS, "rvfv_reference.fa"))
    ap.add_argument("--fill-ref", default=os.path.join(ASSETS, "mp12_reference.fa"))
    ap.add_argument("--primers", default=os.path.join(ASSETS, "rvfv_amplicons_v4.bed"))
    args = ap.parse_args()
    out = args.out or f"{args.sample}_reference.fa"
    if os.path.dirname(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
    import shutil
    if not shutil.which("mafft"):
        sys.exit("mafft not found on PATH — on Beocat run:  module load MAFFT")

    frame = read_fasta(args.frame_ref)
    offsets, pos = {}, 0
    for cid, seq in frame:
        offsets[cid] = (pos, len(seq))
        pos += len(seq)
    total = pos
    fill = dict(read_fasta(args.fill_ref))
    segs_info = []      # for the map / BED / TSV

    with open(out, "w") as fo:
        for seg, (start, seglen) in offsets.items():
            fill_id = FILL_IDS.get(seg)
            if fill_id not in fill:
                sys.exit(f"{seg}: fill-reference record {fill_id} not found in {args.fill_ref}")
            path = os.path.join(args.results, "ivar", f"{args.sample}_{seg}.fa")
            if not os.path.exists(path):
                sys.exit(f"{seg}: consensus not found: {path}")
            full = "".join(s for _, s in read_fasta(path))
            if abs(len(full) - total) > 5:
                sys.exit(f"{path}: length {len(full)} vs reference total {total} — wrong frame, aborting")
            end = len(full) if seg == frame[-1][0] else start + seglen
            cons = full[start:end]

            aln = mafft_pair(fill_id, fill[fill_id], seg, cons)
            out_seq, kept, filled, diffs, dropped = [], 0, 0, 0, 0
            filled_pos, diff_list = [], []
            for r, c in zip(aln[fill_id], aln[seg]):
                if r == "-":
                    dropped += c not in ("N", "-")
                    continue
                if c in ("N", "-"):
                    out_seq.append(r)
                    filled += 1
                    filled_pos.append(len(out_seq))                  # 1-based
                else:
                    out_seq.append(c)
                    kept += 1
                    if c != r:
                        diffs += 1
                        diff_list.append((len(out_seq), r, c))       # pos, fill base, sample base
            seq = "".join(out_seq)
            if len(seq) != seglen or len(seq) != len(fill[fill_id]):
                sys.exit(f"{seg}: output {len(seq)} bp != reference {seglen} bp — aborting")
            if "N" in seq:
                sys.exit(f"{seg}: N remains in output — aborting")
            fo.write(f">{seg}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
            segs_info.append(dict(seg=seg, length=len(seq), fill_id=fill_id, kept=kept, filled=filled,
                                  filled_runs=to_runs(filled_pos), diffs=diff_list))
            print(f"{seg}: {len(seq)} bp | {kept} from {args.sample} consensus, {filled} filled from "
                  f"{fill_id} | {diffs} position(s) differ from {fill_id}"
                  + (f" | {dropped} consensus insertion base(s) dropped" if dropped else ""))
    base = out[:-3] if out.endswith(".fa") else out
    with open(f"{base}.filled.bed", "w") as fb:
        for si in segs_info:
            for a, b in si["filled_runs"]:
                fb.write(f"{si['seg']}\t{a - 1}\t{b}\tfilled_from_{si['fill_id']}\n")   # BED 0-based
    with open(f"{base}.diffs.tsv", "w") as fd:
        fd.write(f"segment\tposition\t{FILL_LABEL}_base\t{args.sample}_base\n")
        for si in segs_info:
            for p_, r, c in si["diffs"]:
                fd.write(f"{si['seg']}\t{p_}\t{r}\t{c}\n")
    write_map(f"{base}.map.html", args.sample, segs_info, read_primers(args.primers))
    print(f"\n-> {out}  (use with: --reference {os.path.abspath(out)})")
    print(f"-> {base}.map.html  (download and open in a browser)")
    print(f"-> {base}.filled.bed, {base}.diffs.tsv")


if __name__ == "__main__":
    main()
