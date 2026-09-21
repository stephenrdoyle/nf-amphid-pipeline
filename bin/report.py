#!/usr/bin/env python3
"""
Run report for nf-helminth-amplicon.

One self-contained HTML file per run. Self-contained is a requirement, not
a preference: this has to open on a laptop in a field site with no network,
so there are no CDN links, no web fonts and no JavaScript libraries. Charts
are inline SVG drawn from the data.

The report is arranged so that the things that would invalidate the results
appear before the results. A plate with a column effect, a run where the
standards did not scale, a set of wells whose negatives are unverified -
these decide whether the call table means anything, so they are not
relegated to an appendix.

Nothing here computes a new number. Every figure is read from the tables
the pipeline already wrote, so the report cannot disagree with them.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table  # noqa: E402

# A sequential scale that survives greyscale printing and the common forms
# of colour blindness. Ordered light to dark, so depth reads as intensity
# rather than as hue.
SCALE = ["#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6",
         "#4292c6", "#2171b5", "#08519c", "#08306b"]

# Okabe-Ito: eight hues chosen to stay distinguishable under the common
# forms of colour vision deficiency, which a default ggplot hue wheel does
# not. Grey is reserved for the pooled remainder, so it never competes with
# a named taxon.
TAXON_COLOURS = ["#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7",
                 "#56B4E9", "#F0E442", "#661100"]
OTHER_COLOUR = "#b9bdc2"

# Read fates. Deliberately not the taxon palette: this panel answers a
# different question, and the eye should not read a fate as a species.
# Ordered usable-first, so the good fraction sits at the bottom of the
# stack directly above the composition it becomes.
FATE_ORDER = ("target", "quant", "control", "vector",
              "ambiguous", "tag_unobservable", "unmatched", "no_primer")
FATE_COLOURS = {
    "target": "#3d6b8e",
    "quant": "#7fa8bf",
    "control": "#c98a5e",
    "vector": "#9b7fa0",
    "ambiguous": "#d9b365",
    "tag_unobservable": "#c6b89a",
    "unmatched": "#c9ced4",
    "no_primer": "#e8eaed",
}

CALL_COLOUR = {
    "detected": "#08519c",
    "not_detected": "#e8eaed",
    "assay_failed": "#b03030",
}

CSS = """
:root { --ink:#1a1c1e; --muted:#5f6368; --line:#d8dbdf; --bg:#ffffff;
        --panel:#f7f8fa; --warn:#8a5300; --bad:#b03030; --good:#0b6b3a; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
       font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",
       Roboto,Helvetica,Arial,sans-serif; }
.wrap { max-width:1100px; margin:0 auto; padding:32px 24px 80px; }
h1 { font-size:24px; margin:0 0 4px; letter-spacing:-0.01em; }
h2 { font-size:17px; margin:40px 0 10px; letter-spacing:-0.01em;
     padding-bottom:6px; border-bottom:1px solid var(--line); }
h3 { font-size:14px; margin:20px 0 6px; color:var(--muted);
     text-transform:uppercase; letter-spacing:0.06em; font-weight:600; }
p { margin:8px 0; max-width:70ch; }
.sub { color:var(--muted); margin:0 0 20px; }
.note { color:var(--muted); font-size:12.5px; max-width:72ch; margin:6px 0 0; }
.cards { display:flex; flex-wrap:wrap; gap:10px; margin:14px 0; }
.card { flex:1 1 150px; background:var(--panel); border:1px solid var(--line);
        border-radius:8px; padding:12px 14px; }
.card .n { font-size:22px; font-weight:650; letter-spacing:-0.02em; }
.card .k { font-size:11.5px; color:var(--muted); text-transform:uppercase;
           letter-spacing:0.05em; }
table { border-collapse:collapse; width:100%; font-size:13px; margin:10px 0; }
th,td { text-align:left; padding:6px 10px; border-bottom:1px solid var(--line);
        vertical-align:top; }
th { font-weight:600; color:var(--muted); font-size:11.5px;
     text-transform:uppercase; letter-spacing:0.05em; }
td.num,th.num { text-align:right; font-variant-numeric:tabular-nums; }
.scroll { max-height:420px; overflow:auto; border:1px solid var(--line);
          border-radius:8px; }
.scroll table { margin:0; }
.scroll thead th { position:sticky; top:0; background:var(--bg);
                   box-shadow:inset 0 -1px 0 var(--line); }
.flag { display:inline-block; padding:1px 7px; border-radius:99px;
        font-size:11.5px; font-weight:600; }
.flag.ok { background:#e6f4ec; color:var(--good); }
.flag.warn { background:#fdf1dd; color:var(--warn); }
.flag.bad { background:#fbe9e9; color:var(--bad); }
.legend { display:flex; align-items:center; gap:6px; font-size:12px;
          color:var(--muted); margin:6px 0 0; }
.sw { width:16px; height:10px; border:1px solid rgba(0,0,0,.12); }
.bar { height:9px; background:var(--panel); border-radius:5px;
       overflow:hidden; display:flex; }
.bar span { display:block; height:100%; }
figure { margin:12px 0; }
code { font:12.5px ui-monospace,SFMono-Regular,Menlo,monospace;
       background:var(--panel); padding:1px 5px; border-radius:4px; }
@media print { .scroll { max-height:none; overflow:visible; } }
"""


def esc(x) -> str:
    return html.escape(str(x))


def fmt(n) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def colour_for(v: float, lo: float, hi: float) -> str:
    if hi <= lo or v is None:
        return SCALE[0]
    t = (v - lo) / (hi - lo)
    return SCALE[max(0, min(len(SCALE) - 1, int(t * (len(SCALE) - 1) + 0.5)))]


# --------------------------------------------------------------------------
# charts
# --------------------------------------------------------------------------


def plate_svg(cells: dict, rows: list, cols: list, lo: float, hi: float,
              title: str, fmt_val=lambda v: f"{v:,.0f}") -> str:
    """A plate as a plate.

    Well position is the one QC axis that a table cannot show. A column
    running cold, a dead row, evaporation down one edge - these are obvious
    in the shape of the plate and invisible in a sorted list.
    """
    cw, ch, pad = 26, 22, 34
    w = pad + len(cols) * cw + 8
    h = pad + len(rows) * ch + 8
    out = [f'<svg viewBox="0 0 {w} {h}" width="100%" '
           f'style="max-width:{w}px" role="img" aria-label="{esc(title)}">']
    for j, c in enumerate(cols):
        out.append(f'<text x="{pad + j*cw + cw/2}" y="{pad-10}" '
                   f'font-size="9" fill="#5f6368" text-anchor="middle">{c}</text>')
    for i, r in enumerate(rows):
        out.append(f'<text x="{pad-8}" y="{pad + i*ch + ch/2 + 3}" '
                   f'font-size="9" fill="#5f6368" text-anchor="end">{r}</text>')
        for j, c in enumerate(cols):
            v = cells.get((r, c))
            fill = "#ffffff" if v is None else colour_for(v, lo, hi)
            stroke = "#e3e6ea" if v is None else "rgba(0,0,0,.10)"
            tip = f"{r}{c}" + ("" if v is None else f": {fmt_val(v)}")
            out.append(
                f'<rect x="{pad + j*cw}" y="{pad + i*ch}" width="{cw-3}" '
                f'height="{ch-3}" rx="3" fill="{fill}" stroke="{stroke}">'
                f'<title>{esc(tip)}</title></rect>'
            )
    out.append("</svg>")
    return "".join(out)


def _num(v: float) -> str:
    """Readable rather than scientific. A report is read by people, and
    2.051e+04 is a number nobody checks against a table."""
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1000:
        return f"{v:,.0f}"
    if a >= 1:
        return f"{v:,.3g}"
    return f"{v:.3g}"


def _axis_ticks(lo: float, hi: float, n: int = 5) -> list[float]:
    if hi <= lo:
        return [lo]
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def histogram_svg(values: list, title: str, log=False, bins=26,
                  unit: str = "") -> str:
    """A distribution, drawn so its counts can be read off it.

    Three things this refuses to do, each of which it used to.

    It will not bin a handful of distinct values into twenty-six bins. A
    set of twenty-one numbers that takes only two values is not a
    distribution with a spread; binned, it draws as two bars at opposite
    ends of a wide axis and reads as bimodal when the truth is "all the
    same except one". Few distinct values are drawn as themselves.

    It will not invent an axis. When every value is identical there is no
    range to plot, and stretching one to the next whole number shows a
    spread that is not in the data - so it says so in words instead.

    It will not mix scales. On a log axis a zero cannot be plotted, and
    silently leaving it linear drops it into an arbitrary bin. Such values
    are excluded and counted in the caption.

    And it draws a y-axis, so a bar height is a number the reader can
    check against the tables rather than a shape.
    """
    raw = [v for v in values if v is not None]
    if not raw:
        return ""
    dropped = 0
    if log:
        keep = [v for v in raw if v > 0]
        dropped = len(raw) - len(keep)
        raw = keep
        if not raw:
            return (f'<p class="note">No positive values to plot for '
                    f'{esc(title)}.</p>')
    vals = [math.log10(v) if log else v for v in raw]

    distinct = sorted(set(round(v, 12) for v in vals))
    w, h, pad_l, pad_b, pad_t = 560, 150, 46, 30, 12

    # nothing to distribute
    if len(distinct) == 1:
        shown = (10 ** distinct[0]) if log else distinct[0]
        return (
            f'<p class="note">All {len(raw)} values are '
            f'<b>{_num(shown)}</b>{esc(" " + unit if unit else "")} — '
            f'no distribution to plot.</p>'
        )

    # few enough distinct values that binning would misrepresent them
    categorical = len(distinct) <= 8
    if categorical:
        counts = [(d, sum(1 for v in vals if round(v, 12) == d))
                  for d in distinct]
        labels = [_num(10**d if log else d) for d, _ in counts]
        heights = [c for _, c in counts]
    else:
        lo, hi = min(vals), max(vals)
        counts = [0] * bins
        for v in vals:
            counts[min(bins - 1, int((v - lo) / (hi - lo) * bins))] += 1
        heights = counts
        labels = None

    top = max(heights) or 1
    plot_h = h - pad_t - pad_b
    n = len(heights)
    bw = (w - pad_l - 14) / n

    out = [f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" '
           f'role="img" aria-label="{esc(title)}">']

    # y axis: gridlines at whole counts, so a bar is a readable number
    steps = [0, top] if top <= 2 else [0, top // 2, top]
    for c in sorted(set(steps)):
        y = pad_t + plot_h - plot_h * c / top
        out.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{w-14}" '
                   f'y2="{y:.1f}" stroke="#eceef0"/>')
        out.append(f'<text x="{pad_l-6}" y="{y+3:.1f}" font-size="9" '
                   f'fill="#5f6368" text-anchor="end">{c}</text>')
    out.append(f'<text x="11" y="{pad_t + plot_h/2}" font-size="9" '
               f'fill="#5f6368" text-anchor="middle" '
               f'transform="rotate(-90 11 {pad_t + plot_h/2})">samples</text>')

    for i, c in enumerate(heights):
        bh = plot_h * c / top
        x = pad_l + i * bw
        if c:
            out.append(
                f'<rect x="{x:.1f}" y="{pad_t + plot_h - bh:.1f}" '
                f'width="{max(1.0, bw-1.5):.1f}" height="{bh:.1f}" '
                f'fill="#4292c6" rx="1"><title>{c}</title></rect>'
            )
            if categorical or c == top:
                out.append(
                    f'<text x="{x + bw/2:.1f}" '
                    f'y="{pad_t + plot_h - bh - 3:.1f}" font-size="9" '
                    f'fill="#3c4043" text-anchor="middle">{c}</text>'
                )

    out.append(f'<line x1="{pad_l}" y1="{pad_t+plot_h}" x2="{w-14}" '
               f'y2="{pad_t+plot_h}" stroke="#d8dbdf"/>')

    if categorical:
        for i, lab in enumerate(labels):
            out.append(f'<text x="{pad_l + i*bw + bw/2:.1f}" '
                       f'y="{pad_t+plot_h+13}" font-size="9" fill="#5f6368" '
                       f'text-anchor="middle">{esc(lab)}</text>')
    else:
        for t in _axis_ticks(min(vals), max(vals)):
            frac = ((t - min(vals)) / (max(vals) - min(vals)))
            x = pad_l + frac * (w - pad_l - 14)
            shown = 10 ** t if log else t
            out.append(f'<text x="{x:.1f}" y="{pad_t+plot_h+13}" '
                       f'font-size="9" fill="#5f6368" text-anchor="middle">'
                       f'{_num(shown)}</text>')
    if unit:
        out.append(f'<text x="{w-14}" y="{h-3}" font-size="9" '
                   f'fill="#5f6368" text-anchor="end">{esc(unit)}</text>')
    out.append("</svg>")
    svg = "".join(out)
    if dropped:
        svg += (f'<p class="note">{dropped} non-positive value(s) are not '
                f'shown: they cannot be placed on a logarithmic axis.</p>')
    return svg


def composition_svg(
    order: list, depths: dict, comp: dict, colours: dict,
    min_denom: int, denom: dict, fates: dict | None = None,
) -> str:
    """Depth, then where the reads went, then composition - one x axis.

    The three panels belong together, and the order is the argument.
    Depth says how much data a well produced. The fate panel says what
    fraction of it was usable, which a depth alone cannot: a well can
    carry five thousand reads and fifty targets. Composition then divides
    up that usable fraction - so the target slice of the middle panel is
    literally the denominator of the bottom one, sitting directly above
    it.

    A composition bar looks exactly the same whether it rests on eighty
    reads or eighty thousand, so shown alone it invites a confidence the
    data does not carry.

    Bars whose composition rests on fewer than `min_denom` target reads are
    drawn hatched. That is a stronger statement than low depth: a deep well
    can still have almost no target in it, and it is the target count, not
    the read count, that the proportions are divided by.
    """
    if not order:
        return ""
    n = len(order)
    bw = max(3.0, min(26.0, 900.0 / n))
    gap = 1.0 if bw > 5 else 0.4
    pad_l, pad_r, pad_t = 56, 12, 10
    top_h, mid_h, bot_h = 100, 16, 140
    fate_h = 90 if fates else 0
    lab_h = 58 if n <= 40 else 10
    w = pad_l + n * (bw + gap) + pad_r
    h = pad_t + top_h + mid_h + fate_h + (mid_h if fates else 0) + bot_h \
        + lab_h

    dmax = max([depths.get(s, 0) for s in order] + [1])
    step = 10 ** math.floor(math.log10(dmax)) / 2 or 1
    dtop = math.ceil(dmax / step) * step

    out = [f'<svg viewBox="0 0 {w:.0f} {h:.0f}" width="{w:.0f}" '
           f'height="{h:.0f}" role="img" '
           f'aria-label="read depth and composition per sample">',
           '<defs><pattern id="thin" width="4" height="4" '
           'patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
           '<rect width="4" height="4" fill="#ffffff"/>'
           '<line x1="0" y1="0" x2="0" y2="4" stroke="#c2c6cb" '
           'stroke-width="2"/></pattern></defs>']

    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        y = pad_t + top_h - frac * top_h
        out.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{w-pad_r}" '
                   f'y2="{y:.1f}" stroke="#eceef0"/>')
        out.append(f'<text x="{pad_l-6}" y="{y+3:.1f}" font-size="9" '
                   f'fill="#5f6368" text-anchor="end">'
                   f'{int(dtop*frac):,}</text>')
    out.append(f'<text x="10" y="{pad_t + top_h/2}" font-size="10" '
               f'fill="#3c4043" text-anchor="middle" '
               f'transform="rotate(-90 10 {pad_t + top_h/2})">Read depth</text>')

    for i, sid in enumerate(order):
        x = pad_l + i * (bw + gap)
        d = depths.get(sid, 0)
        bh = top_h * (d / dtop) if dtop else 0
        out.append(f'<rect x="{x:.1f}" y="{pad_t + top_h - bh:.1f}" '
                   f'width="{bw:.1f}" height="{bh:.1f}" fill="#8c9196">'
                   f'<title>{esc(sid)}: {d:,} sequences</title></rect>')

    # --- read fate panel
    if fates:
        yf = pad_t + top_h + mid_h
        for frac in (0, 0.5, 1.0):
            y = yf + fate_h - frac * fate_h
            out.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{w-pad_r}" '
                       f'y2="{y:.1f}" stroke="#eceef0"/>')
            out.append(f'<text x="{pad_l-6}" y="{y+3:.1f}" font-size="9" '
                       f'fill="#5f6368" text-anchor="end">'
                       f'{int(frac*100)}%</text>')
        out.append(f'<text x="10" y="{yf + fate_h/2}" font-size="10" '
                   f'fill="#3c4043" text-anchor="middle" '
                   f'transform="rotate(-90 10 {yf + fate_h/2})">'
                   f'Reads mapped</text>')
        for i, sid in enumerate(order):
            x = pad_l + i * (bw + gap)
            f = fates.get(sid) or {}
            tot = sum(f.values())
            if not tot:
                out.append(f'<rect x="{x:.1f}" y="{yf:.1f}" '
                           f'width="{bw:.1f}" height="{fate_h}" '
                           f'fill="#f4f5f6" stroke="#e6e8ea"/>')
                continue
            acc = 0.0
            for fate in FATE_ORDER:
                v = f.get(fate, 0)
                if not v:
                    continue
                fh_ = fate_h * v / tot
                out.append(
                    f'<rect x="{x:.1f}" y="{yf + fate_h - acc - fh_:.2f}" '
                    f'width="{bw:.1f}" height="{fh_:.2f}" '
                    f'fill="{FATE_COLOURS[fate]}">'
                    f'<title>{esc(sid)} — {esc(fate)}: {v:,} reads, '
                    f'{100*v/tot:.1f}% of {tot:,}</title></rect>'
                )
                acc += fh_

    y0 = pad_t + top_h + mid_h + (fate_h + mid_h if fates else 0)
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        y = y0 + bot_h - frac * bot_h
        out.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{w-pad_r}" '
                   f'y2="{y:.1f}" stroke="#eceef0"/>')
        out.append(f'<text x="{pad_l-6}" y="{y+3:.1f}" font-size="9" '
                   f'fill="#5f6368" text-anchor="end">{frac:.2f}</text>')
    out.append(f'<text x="10" y="{y0 + bot_h/2}" font-size="10" '
               f'fill="#3c4043" text-anchor="middle" '
               f'transform="rotate(-90 10 {y0 + bot_h/2})">'
               f'Relative abundance</text>')

    for i, sid in enumerate(order):
        x = pad_l + i * (bw + gap)
        shares = comp.get(sid, {})
        dn = denom.get(sid, 0)
        if not shares or not dn:
            out.append(f'<rect x="{x:.1f}" y="{y0:.1f}" width="{bw:.1f}" '
                       f'height="{bot_h}" fill="#f4f5f6" stroke="#e6e8ea">'
                       f'<title>{esc(sid)}: no target reads</title></rect>')
            continue
        acc = 0.0
        for taxon, share in sorted(shares.items(), key=lambda kv: -kv[1]):
            sh = bot_h * share
            yy = y0 + bot_h - acc - sh
            fill = colours.get(taxon, OTHER_COLOUR)
            out.append(
                f'<rect x="{x:.1f}" y="{yy:.1f}" width="{bw:.1f}" '
                f'height="{sh:.2f}" fill="{fill}">'
                f'<title>{esc(sid)} — {esc(taxon)}: {share*100:.1f}% '
                f'of {dn:,} target reads</title></rect>'
            )
            acc += sh
        if dn < min_denom:
            out.append(
                f'<rect x="{x:.1f}" y="{y0:.1f}" width="{bw:.1f}" '
                f'height="{bot_h}" fill="url(#thin)" opacity="0.72">'
                f'<title>{esc(sid)}: only {dn:,} target reads — these '
                f'proportions are not well determined</title></rect>'
            )

    if lab_h > 20:
        for i, sid in enumerate(order):
            x = pad_l + i * (bw + gap) + bw / 2
            out.append(f'<text x="{x:.1f}" y="{y0 + bot_h + 8}" '
                       f'font-size="9" fill="#5f6368" text-anchor="end" '
                       f'transform="rotate(-60 {x:.1f} {y0 + bot_h + 8})">'
                       f'{esc(sid)}</text>')
    out.append("</svg>")
    return "".join(out)


def stacked_bar(parts: list) -> str:
    total = sum(v for _, v, _ in parts) or 1
    segs = "".join(
        f'<span style="width:{100*v/total:.3f}%;background:{c}" '
        f'title="{esc(k)}: {fmt(v)}"></span>'
        for k, v, c in parts if v
    )
    return f'<div class="bar">{segs}</div>'


# --------------------------------------------------------------------------


def flag(level: str, text: str) -> str:
    return f'<span class="flag {level}">{esc(text)}</span>'


def build(args) -> str:
    calls = read_table(os.path.join(args.synthesis_dir, "calls.tsv"))
    with open(os.path.join(args.synthesis_dir, "run_qc.json")) as fh:
        qc = json.load(fh)
    cal_rows = read_table(os.path.join(args.synthesis_dir, "calibration.tsv"))
    fates: dict = {}
    fpath = os.path.join(args.synthesis_dir, "read_fates.tsv")
    if os.path.exists(fpath):
        for r in read_table(fpath):
            fates[r["sample_id"]] = {
                k: int(r.get(k) or 0) for k in FATE_ORDER
            }

    comparison = None
    if args.comparison and os.path.exists(args.comparison):
        with open(args.comparison) as fh:
            comparison = json.load(fh)

    run_id = qc.get("run_id") or args.run_id or "run"
    by_sample: dict = defaultdict(list)
    for c in calls:
        by_sample[c["sample_id"]].append(c)

    depths = {s: int(rows[0].get("depth") or 0)
              for s, rows in by_sample.items()}
    wells = {}
    for s, rows in by_sample.items():
        w = rows[0].get("well") or ""
        if w and w[0].isalpha():
            wells[s] = (w[0], int(w[1:]) if w[1:].isdigit() else 0)

    cs = qc.get("call_summary", {})
    ic = qc.get("internal_control_summary", {})
    d = qc.get("depth", {})
    parts = []

    # ---- header -------------------------------------------------------
    parts.append(f"<h1>{esc(run_id)}</h1>")
    parts.append(
        f'<p class="sub">{fmt(qc.get("n_samples", 0))} samples &middot; '
        f'{fmt(len(calls))} sample-by-target calls &middot; '
        f'generated by nf-helminth-amplicon</p>'
    )

    # ---- what would invalidate this -----------------------------------
    parts.append("<h2>Before the results</h2>")
    parts.append(
        "<p>Three things decide whether the call table below means anything. "
        "They are here rather than in an appendix because a run can fail "
        "any of them and still produce a full set of confident-looking "
        "numbers.</p>"
    )
    concerns = []

    failed = cs.get("assay_failed", 0)
    if failed:
        concerns.append((
            "bad" if failed > 0.1 * len(calls) else "warn",
            f"{fmt(failed)} calls are assay_failed",
            "Those wells cannot support a negative. They are not negatives.",
        ))
    unver = ic.get("unspiked", 0)
    if unver:
        concerns.append((
            "warn", f"{fmt(unver)} wells have no recovered standard",
            "Their negatives are unverified: nothing establishes that the "
            "well could have produced a positive.",
        ))
    sp = ic.get("spike_failed", 0)
    if sp:
        concerns.append((
            "bad", f"{fmt(sp)} wells failed their spike-in",
            "Standards were expected and did not appear.",
        ))
    slope = qc.get("calibration_slope", {})
    if slope.get("outside_tolerance"):
        concerns.append((
            "warn",
            f"{slope['outside_tolerance']} of {slope.get('n', 0)} wells "
            "have a calibration slope away from 1",
            "Copy numbers from those wells are compressed or stretched "
            "relative to the standards.",
        ))
    spread = d.get("fold_spread", 0)
    if spread and spread > 10:
        concerns.append((
            "warn", f"{spread}x spread between the deepest and shallowest well",
            "The same absolute detection limit is a much weaker claim in a "
            "shallow well than in a deep one.",
        ))
    geo = [g for g in qc.get("plate_geometry", [])
           if min(g["row_p"], g["col_p"], g["edge_p"]) < 0.01]
    if geo:
        concerns.append((
            "warn",
            f"{len(geo)} quantities vary with position on the plate",
            "Wells hold unrelated samples, so position should explain "
            "nothing. See the plate maps below.",
        ))
    fn = qc.get("read_funnel") or {}
    if fn.get("handoff_balances") is False:
        concerns.append((
            "bad", "reads were lost between preparation and assignment",
            "No reported filter accounts for them.",
        ))
    if fn.get("wells_losing_over_half_at_prep"):
        concerns.append((
            "warn",
            f"{len(fn['wells_losing_over_half_at_prep'])} well(s) lost over "
            "half their reads before assignment",
            "They appear merely shallow in the depth panel; the reaction "
            "is the problem, not the loading.",
        ))

    cf = qc.get("corr_factor", {})
    if cf.get("references_using_default"):
        concerns.append((
            "ok",
            f"{cf['references_using_default']} references use the default "
            "corr_factor of 1.0",
            "Copy numbers are relative to the standards' amplification "
            "efficiency rather than corrected to it - sound for comparing "
            "samples, not for comparing species.",
        ))

    if concerns:
        parts.append("<table><thead><tr><th></th><th>Finding</th>"
                     "<th>What it means</th></tr></thead><tbody>")
        for lvl, what, why in concerns:
            parts.append(
                f"<tr><td>{flag(lvl, lvl.upper())}</td><td>{esc(what)}</td>"
                f"<td>{esc(why)}</td></tr>"
            )
        parts.append("</tbody></table>")
    else:
        parts.append(f"<p>{flag('ok', 'CLEAR')} No run-level concern was "
                     "raised. That is not the same as a clean run: the "
                     "sensors below say what was actually checked.</p>")

    # ---- the funnel ----------------------------------------------------
    funnel = qc.get("read_funnel") or {}
    stages = funnel.get("stages") or []
    if stages:
        parts.append("<h2>Where the reads went</h2>")
        parts.append(
            "<p>Every read from the sequencer to a species call. The "
            "per-sample figures below begin at <i>reached assignment</i>, "
            "which is already net of everything preparation discarded — so "
            "this is the only place a well that produced plenty of data and "
            "lost it to primer dimer looks different from a well that was "
            "simply shallow. The two call for different fixes.</p>"
        )
        start = stages[0]["reads"] or 1
        parts.append("<table><thead><tr><th>Stage</th><th class='num'>Reads"
                     "</th><th class='num'>% of start</th><th></th>"
                     "<th>Lost here</th></tr></thead><tbody>")
        for st in stages:
            pct = st.get("percent_of_start")
            frac = st["reads"] / start
            bar = (f'<div class="bar"><span style="width:{100*frac:.2f}%;'
                   f'background:#3d6b8e"></span></div>')
            lost = ", ".join(f"{esc(k)} {fmt(v)}"
                             for k, v in st["lost_here"].items()) or "—"
            parts.append(
                f"<tr><td>{esc(st['stage'])}</td>"
                f"<td class='num'>{fmt(st['reads'])}</td>"
                f"<td class='num'>{'' if pct is None else f'{pct:.1f}%'}</td>"
                f"<td style='width:34%'>{bar}</td><td>{lost}</td></tr>"
            )
        parts.append("</tbody></table>")

        if funnel.get("handoff_balances") is False:
            parts.append(
                f"<p>{flag('bad', 'MISMATCH')} Reads left preparation that "
                "never arrived at assignment. Nothing else in the pipeline "
                "compares those two counts, so this is data lost between "
                "steps rather than by any reported filter.</p>"
            )
        bad = funnel.get("wells_losing_over_half_at_prep") or []
        if bad:
            parts.append(
                f"<p>{flag('warn', 'WELLS')} "
                f"{len(bad)} well(s) lost more than half their reads before "
                "assignment — the depth panel shows these as shallow, but "
                "the reaction rather than the loading is the problem:</p>"
            )
            parts.append("<table><thead><tr><th>Sample</th>"
                         "<th class='num'>Pairs in</th>"
                         "<th class='num'>Kept</th><th>Where they went</th>"
                         "</tr></thead><tbody>")
            for b in bad:
                where = ", ".join(f"{esc(k)} {fmt(v)}"
                                  for k, v in b["fates"].items() if v)
                parts.append(
                    f"<tr><td>{esc(b['sample_id'])}</td>"
                    f"<td class='num'>{fmt(b['pairs_in'])}</td>"
                    f"<td class='num'>{100*b['kept_fraction']:.0f}%</td>"
                    f"<td>{where}</td></tr>"
                )
            parts.append("</tbody></table>")
        if funnel.get("prep_loss_reasons"):
            parts.append(
                '<p class="note">'
                + " &middot; ".join(
                    f"<b>{esc(k)}</b> {esc(v)}"
                    for k, v in funnel["prep_loss_reasons"].items())
                + "</p>"
            )

    # ---- calls --------------------------------------------------------
    parts.append("<h2>Calls</h2>")
    parts.append('<div class="cards">')
    for k in ("detected", "not_detected", "assay_failed"):
        parts.append(f'<div class="card"><div class="n">{fmt(cs.get(k, 0))}'
                     f'</div><div class="k">{esc(k.replace("_", " "))}</div></div>')
    parts.append('</div>')
    parts.append(stacked_bar([
        (k, cs.get(k, 0), CALL_COLOUR[k])
        for k in ("detected", "not_detected", "assay_failed")
    ]))
    parts.append(
        '<p class="note">'
        + " &middot; ".join(f"<b>{esc(k)}</b> {esc(v)}"
                            for k, v in qc.get("call_meaning", {}).items())
        + "</p>"
    )

    # by taxon
    tax = defaultdict(Counter)
    for c in calls:
        tax[c.get("taxon") or c["seq_id"]][c["call"]] += 1
    parts.append("<h3>By taxon</h3><table><thead><tr><th>Taxon</th>"
                 "<th class='num'>Detected</th><th class='num'>Not detected</th>"
                 "<th class='num'>Assay failed</th><th>Share detected</th>"
                 "</tr></thead><tbody>")
    for t, cc in sorted(tax.items(), key=lambda kv: -kv[1]["detected"]):
        n = sum(cc.values()) or 1
        parts.append(
            f"<tr><td>{esc(t)}</td><td class='num'>{fmt(cc['detected'])}</td>"
            f"<td class='num'>{fmt(cc['not_detected'])}</td>"
            f"<td class='num'>{fmt(cc['assay_failed'])}</td>"
            f"<td>{stacked_bar([(k, cc[k], CALL_COLOUR[k]) for k in CALL_COLOUR])}</td></tr>"
        )
    parts.append("</tbody></table>")

    # ---- composition ---------------------------------------------------
    # Only speciation targets form a composition. A panel carrying a second
    # locus for variant calling would otherwise count a species twice, with
    # a denominator mixing loci at different copy numbers.
    spec = [c for c in calls
            if "speciation" in (c.get("purpose") or "speciation")]
    comp: dict = {}
    denom: dict = {}
    totals: Counter = Counter()
    for c in spec:
        totals[c.get("taxon") or c["seq_id"]] += int(c["reads"] or 0)
    ranked = [t for t, n in totals.most_common() if n > 0]
    named = ranked[:len(TAXON_COLOURS)]
    colours = {t: TAXON_COLOURS[i] for i, t in enumerate(named)}
    pooled = set(ranked[len(TAXON_COLOURS):])

    for s_id in by_sample:
        rows_s = [c for c in spec if c["sample_id"] == s_id]
        tot = sum(int(c["reads"] or 0) for c in rows_s)
        denom[s_id] = tot
        if not tot:
            continue
        acc: Counter = Counter()
        for c in rows_s:
            n = int(c["reads"] or 0)
            if not n:
                continue
            t = c.get("taxon") or c["seq_id"]
            acc["other" if t in pooled else t] += n
        comp[s_id] = {t: n / tot for t, n in acc.items()}
    if pooled:
        colours["other"] = OTHER_COLOUR

    # Plate order where it is known: a column effect is far easier to see
    # when neighbouring bars are neighbouring wells.
    order = sorted(by_sample, key=lambda s: (wells.get(s, ("Z", 999)), s))
    min_denom = int(args.min_composition_reads)
    if comp:
        parts.append("<h2>Composition</h2>")
        parts.append(
            "<p>Depth above, composition below, on one axis. The two belong "
            "together: a composition bar looks the same whether it rests on "
            "eighty reads or eighty thousand, and shown on its own it "
            "invites a confidence the data does not carry.</p>"
        )
        parts.append('<div class="scroll" style="max-height:none">'
                     + composition_svg(order, depths, comp, colours,
                                       min_denom, denom, fates)
                     + "</div>")
        if fates:
            used = [f for f in FATE_ORDER
                    if any(v.get(f, 0) for v in fates.values())]
            parts.append(
                '<div class="legend"><span>reads mapped:</span>'
                + "".join(
                    f'<span class="sw" style="background:{FATE_COLOURS[f]}">'
                    f'</span><span>{esc(f.replace("_", " "))}</span>'
                    for f in used)
                + "</div>"
            )
            tf = [v["target"] / sum(v.values()) for v in fates.values()
                  if sum(v.values())]
            if tf:
                tf.sort()
                lowest = [s_ for s_, v in sorted(fates.items())
                          if sum(v.values())
                          and v["target"] / sum(v.values()) < 0.5]
                parts.append(
                    f'<p class="note">The middle panel is where the reads '
                    f'went. Its target slice is the denominator of the '
                    f'composition below it, so a tall stack of anything '
                    f'else is a well whose composition rests on very '
                    f'little. Median target fraction '
                    f'<b>{100*tf[len(tf)//2]:.1f}%</b>'
                    + (f"; {len(lowest)} sample(s) below 50%."
                       if lowest else ".")
                    + ' Each fate means something different and calls for a '
                    'different fix: control is plasmid carryover, vector is '
                    'non-specific amplification, unmatched may be a gap in '
                    'the reference panel rather than noise.</p>'
                )

        parts.append(
            '<div class="legend"><span>composition:</span>'
            + "".join(
                f'<span class="sw" style="background:{c}"></span>'
                f'<span>{esc(t)}</span>'
                for t, c in colours.items())
            + '<span class="sw" style="background:'
              'repeating-linear-gradient(45deg,#fff,#fff 2px,#c2c6cb 2px,'
              '#c2c6cb 4px)"></span>'
            + f'<span>fewer than {min_denom:,} target reads</span></div>'
        )
        thin = [s for s in order if 0 < denom.get(s, 0) < min_denom]
        empty = [s for s in order if not denom.get(s, 0)]
        parts.append(
            f'<p class="note">Hatched bars rest on fewer than '
            f'{min_denom:,} target reads ({len(thin)} of {len(order)} '
            f'samples). That is a stronger warning than low depth: a deep '
            f'well can still hold almost no target, and it is the target '
            f'count, not the read count, that these proportions are '
            f'divided by. {len(empty)} samples have no target reads at all '
            f'and are drawn empty rather than as an arbitrary split.</p>'
        )
        if pooled:
            parts.append(
                f'<p class="note">{len(pooled)} less abundant taxa are '
                f'pooled as "other": {esc(", ".join(sorted(pooled)))}.</p>'
            )

    # ---- plate --------------------------------------------------------
    if wells:
        rows = sorted({r for r, _ in wells.values()})
        cols = sorted({c for _, c in wells.values()})
        parts.append("<h2>The plate</h2>")
        parts.append(
            "<p>Every well holds an unrelated sample, so these maps should "
            "look like noise. Structure in them is a statement about the "
            "plate, not about the samples.</p>"
        )
        cells = {wells[s]: v for s, v in depths.items() if s in wells}
        vals = [v for v in cells.values() if v]
        if vals:
            parts.append("<h3>Sequencing depth</h3>")
            parts.append(plate_svg(cells, rows, cols, min(vals), max(vals),
                                   "depth by well"))
            parts.append(
                f'<div class="legend"><span>{fmt(min(vals))}</span>'
                + "".join(f'<span class="sw" style="background:{c}"></span>'
                          for c in SCALE)
                + f'<span>{fmt(max(vals))} sequences</span></div>'
            )
        # one map per taxon actually detected somewhere
        detected_tax = [t for t, cc in tax.items() if cc["detected"]]
        for t in detected_tax[:4]:
            cells_t = {}
            for s, rws in by_sample.items():
                if s not in wells:
                    continue
                v = sum(int(r["reads"]) for r in rws
                        if (r.get("taxon") or r["seq_id"]) == t)
                cells_t[wells[s]] = v
            tv = [v for v in cells_t.values() if v]
            if not tv:
                continue
            parts.append(f"<h3>{esc(t)} &mdash; reads per well</h3>")
            parts.append(plate_svg(cells_t, rows, cols, 0, max(tv),
                                   f"{t} by well"))

        if qc.get("plate_geometry"):
            parts.append("<h3>Position effects</h3>")
            parts.append("<table><thead><tr><th>Quantity</th>"
                         "<th class='num'>Row p</th><th class='num'>Column p</th>"
                         "<th class='num'>Edge p</th></tr></thead><tbody>")
            for g in qc["plate_geometry"]:
                cellf = lambda p: (f"<td class='num'>{flag('warn', f'{p:.4f}')}</td>"
                                   if p < 0.01 else f"<td class='num'>{p:.4f}</td>")
                parts.append(
                    f"<tr><td>{esc(g['quantity'])}</td>"
                    + cellf(g["row_p"]) + cellf(g["col_p"]) + cellf(g["edge_p"])
                    + "</tr>"
                )
            parts.append("</tbody></table>")
            parts.append(f'<p class="note">{esc(qc.get("plate_geometry_note", ""))}</p>')

    # ---- quantification ------------------------------------------------
    parts.append("<h2>Quantification</h2>")
    slopes = [float(r["slope"]) for r in cal_rows if r.get("slope")]
    copies = [float(c["copies_per_unit_dna"]) for c in calls
              if c.get("copies_per_unit_dna")]
    parts.append('<div class="cards">')
    parts.append(f'<div class="card"><div class="n">'
                 f'{slope.get("median", "-")}</div>'
                 f'<div class="k">median slope</div></div>')
    parts.append(f'<div class="card"><div class="n">'
                 f'{slope.get("outside_tolerance", 0)}</div>'
                 f'<div class="k">slopes off unity</div></div>')
    parts.append(f'<div class="card"><div class="n">'
                 f'{fmt(qc.get("depth", {}).get("median", 0))}</div>'
                 f'<div class="k">median depth</div></div>')
    parts.append(f'<div class="card"><div class="n">{spread}x</div>'
                 f'<div class="k">depth spread</div></div>')
    parts.append("</div>")
    parts.append(
        "<p>A linear assay puts a slope of 1 through its two standards. "
        "A slope below 1 means the high standard is saturating, and every "
        "copy number from that well is compressed toward the middle.</p>"
    )
    if slopes:
        parts.append("<h3>Calibration slope, per well</h3>")
        parts.append(histogram_svg(slopes, "calibration slopes"))
    if copies:
        parts.append("<h3>Copies per unit DNA, where detected</h3>")
        parts.append(histogram_svg(copies, "copies", log=True,
                                   unit="copies per unit DNA"))
    if cf:
        parts.append(f'<p class="note">{esc(cf.get("note", ""))}</p>')

    # ---- detection limits ----------------------------------------------
    lods = qc.get("detection_limits", {})
    if lods:
        parts.append("<h2>Detection limits</h2>")
        parts.append(
            "<p>Measured from this run's own negative controls, which carry "
            "the same index hopping, carryover and reagent contamination as "
            "everything else. Two bounds are computed and the stricter "
            "applies per well.</p>"
        )
        nb = next(iter(lods.values()), {}).get("n_blanks", 0)
        if nb < 2:
            parts.append(f"<p>{flag('warn', 'THIN')} Only {nb} negative "
                         "control on this run. A limit from one blank is an "
                         "observation, not a distribution.</p>")
        parts.append("<table><thead><tr><th>Reference</th>"
                     "<th class='num'>Worst blank (reads)</th>"
                     "<th class='num'>Worst blank (fraction)</th>"
                     "<th class='num'>Floor</th></tr></thead><tbody>")
        for sid, v in sorted(lods.items(),
                             key=lambda kv: -kv[1]["max_reads_in_blank"])[:14]:
            parts.append(
                f"<tr><td>{esc(sid)}</td>"
                f"<td class='num'>{fmt(v['max_reads_in_blank'])}</td>"
                f"<td class='num'>{v['max_fraction_in_blank']:.2e}</td>"
                f"<td class='num'>{fmt(v['absolute_floor'])}</td></tr>"
            )
        parts.append("</tbody></table>")

    # ---- cross-run ------------------------------------------------------
    if comparison:
        parts.append("<h2>Agreement with the paired run</h2>")
        k = comparison.get("cohen_kappa")
        parts.append('<div class="cards">')
        parts.append(f'<div class="card"><div class="n">'
                     f'{comparison.get("agreement", 0):.3f}</div>'
                     f'<div class="k">raw agreement</div></div>')
        parts.append(f'<div class="card"><div class="n">'
                     f'{"-" if k is None else f"{k:.3f}"}</div>'
                     f'<div class="k">cohen kappa</div></div>')
        w2 = comparison.get("within_twofold")
        parts.append(f'<div class="card"><div class="n">'
                     f'{"-" if w2 is None else f"{100*w2:.0f}%"}</div>'
                     f'<div class="k">copies within 2x</div></div>')
        parts.append(f'<div class="card"><div class="n">'
                     f'{fmt(comparison.get("paired_observations", 0))}</div>'
                     f'<div class="k">paired observations</div></div>')
        parts.append("</div>")
        parts.append(f'<p class="note">{esc(comparison.get("kappa_note", ""))}</p>')
        only = comparison.get("detected_by_one_only", {})
        if only:
            dd = comparison.get("depth_of_disagreements") or {}
            parts.append(
                "<p>Found by one run only: "
                + ", ".join(f"<b>{esc(kk)}</b> {vv}" for kk, vv in only.items())
                + (f". Of {dd.get('n_disagreements_with_depth', 0)} such "
                   "disagreements with depth recorded, most were found by the "
                   "deeper run &mdash; which is a difference in sampling, not "
                   "between platforms." if dd else "")
                + "</p>"
            )

    # ---- per sample -----------------------------------------------------
    parts.append("<h2>Per sample</h2>")
    parts.append('<div class="scroll"><table><thead><tr>'
                 "<th>Sample</th><th>Type</th><th>Well</th>"
                 "<th class='num'>Depth</th><th>Internal control</th>"
                 "<th class='num'>Detected</th><th>Taxa detected</th>"
                 "</tr></thead><tbody>")
    cal_by = {r["sample_id"]: r for r in cal_rows}
    for s in sorted(by_sample):
        rws = by_sample[s]
        det = [r for r in rws if r["call"] == "detected"]
        state = rws[0].get("internal_control", "")
        lvl = {"spiked": "ok", "unspiked": "warn",
               "spike_failed": "bad"}.get(state, "warn")
        parts.append(
            f"<tr><td>{esc(s)}</td><td>{esc(rws[0].get('sample_type',''))}</td>"
            f"<td>{esc(rws[0].get('well',''))}</td>"
            f"<td class='num'>{fmt(rws[0].get('depth',0))}</td>"
            f"<td>{flag(lvl, state)}</td>"
            f"<td class='num'>{len(det)}</td>"
            f"<td>{esc(', '.join(sorted({r.get('taxon') or r['seq_id'] for r in det})) or '—')}</td></tr>"
        )
    parts.append("</tbody></table></div>")
    parts.append(
        '<p class="note">Copy numbers, detection limits and the reason '
        'behind every call are in <code>calls.tsv</code>; this table is a '
        'summary of it, not a replacement.</p>'
    )

    body = "\n".join(parts)
    return (
        "<!doctype html>\n<html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{esc(run_id)} — nf-helminth-amplicon</title>"
        f"<style>{CSS}</style></head><body><div class='wrap'>{body}"
        "</div></body></html>\n"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--synthesis-dir", required=True)
    ap.add_argument("--comparison", help="comparison.json, if there is a "
                                         "paired run")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--min-composition-reads", type=int, default=100,
                    help="target reads below which a composition bar is "
                         "hatched as poorly determined (default: 100)")
    ap.add_argument("--out", default="report.html")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    doc = build(args)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    if not args.quiet:
        print(f"{args.out}  ({len(doc)/1024:.0f} KB, self-contained)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
