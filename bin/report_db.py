#!/usr/bin/env python3
"""
Run report for database mode (--reference_database).

The same rules as the panel report (report.py), whose style and helpers
it shares: one self-contained HTML file, no network, no JavaScript
libraries, nothing computed that the tables do not already hold, and the
things that would undermine the results placed before the results.

Its shape follows VESPA's own figures (Owens et al. 2024, Figs 3-5): every
well's reads by organism group first, then each target group - helminths,
protists, microsporidia - on its own, shares recalculated within the group,
with a presence table and prevalence per taxon; then off-target.

What it puts first is different from the panel report, because what can
undermine a database run is different. There are no standards and no plasmids, so the
questions are: which release was this, were there blanks and what did
they show, how much of each well was not the target at all, and how
much of the target could be named only above species. Then the
composition, at the rank each call actually reached and labelled as
relative 18S abundance - never organisms, never copies.

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table  # noqa: E402
from report import (  # noqa: E402
    CSS, OTHER_COLOUR, TAXON_COLOURS, esc, flag, fmt, stacked_bar,
)

FUNNEL = (("clustered", "#3d6b8e", "classified"),
          ("unclustered", "#9fb4c4", "too few reads for a consensus"),
          ("too_short", "#c98a5e", "too short"),
          ("low_quality", "#d9a0a0", "low quality (Illumina)"),
          ("no_primer", "#c9c2b0", "no primer pair"))
# One colour per organism group, fixed so a group is the same colour in
# every run's report. A fixed-order categorical palette chosen to stay
# distinguishable under the common colour-vision deficiencies; the
# off-target remainder in greys so the eye goes to the targets.
GROUP_STYLE = (
    ("helminths", "#2a78d6", "helminths"),
    ("protists", "#1baf7a", "protists"),
    ("microsporidia", "#e87ba4", "microsporidia"),
    ("fungi", "#eda100", "fungi"),
    ("host", "#eb6834", "host"),
    ("plants", "#6250d6", "plants"),
    ("invertebrates", "#9085e9", "other invertebrates"),
    ("bacteria_archaea", "#888780", "bacteria & archaea"),
    ("organelle", "#b4b2a9", "organelles"),
    ("unassigned", "#d3d1c7", "unassigned"),
    ("chimera", "#e34948", "PCR chimeras"),
)
GROUP_COLOUR = {g: c for g, c, _l in GROUP_STYLE}
GROUP_LABEL = {g: lab for g, _c, lab in GROUP_STYLE}
TARGET = ("helminths", "protists", "microsporidia")
# Within a target group a run easily names more than the panel palette's
# eight taxa, and the ninth onwards went grey - one helminth in three on
# the fixture. Okabe-Ito, then Paul Tol's muted set, which stays
# distinguishable under the common colour-vision deficiencies.
WITHIN_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                  "#008300", "#6250d6", "#e34948", "#332288", "#88CCEE",
                  "#44AA99", "#999933", "#AA4499", "#882255"]
# Sequential: darkest = most specific.
RANK_COLOURS = {"species": "#0c447c", "genus": "#378add",
                "family": "#85b7eb", "above": "#d3d1c7"}
CHART_INK, CHART_GRID, CHART_AXIS = "#52514e", "#e1e0d9", "#c3c2b7"


def columns_svg(labels: list, series: list, height: int = 300,
                note: str = "", totals: list | None = None,
                totals_label: str = "reads") -> str:
    """Vertical stacked columns, one per sample, each scaled to 100%.

    `series` is [(name, colour, [value per label])]; each column is the
    values' shares of that column's total, stacked in series order from
    the bottom. Drawn as inline SVG rather than with a charting library:
    the report is one self-contained file that must open with no network.
    Gridlines every 10%, the axis in muted ink, a 2px gap between stacked
    segments, labels angled when there are many columns; every segment
    carries a hover title with its name, share and value. Columns with no
    data are drawn as an empty outline, so a blank well is visible as one.

    With `totals`, a small bar panel sits above the columns, aligned with
    them: the read count each column's percentages rest on. On a log
    scale, because counts run from a handful to tens of thousands and the
    shallow wells - where a 100% column means least - are the ones that
    must not vanish against the deepest.
    """
    n = len(labels)
    if not n:
        return ""
    th = 74 if totals is not None else 0       # height of the totals panel
    height += th
    left, top, bottom = 44, 8 + th, 64
    slot = max(14.0, min(56.0, 900.0 / n))
    width = left + slot * n + 8
    bar = min(24.0, slot * 0.72)
    ph = height - top - bottom
    out = [f"<svg viewBox='0 0 {width:.0f} {height}' width='100%' "
           f"style='max-width:{max(width, 360):.0f}px;display:block' "
           f"role='img' aria-label='{esc(note)}'>"]
    if totals is not None:
        import math
        tmax = max([t for t in totals if t] or [1])
        decades = max(1, math.ceil(math.log10(max(tmax, 10))))
        t0, tb = 6, 6 + th - 18                 # panel top and baseline
        tph = tb - t0
        def ty(v):
            return tb - tph * min(1.0, math.log10(max(v, 1)) / decades)
        for d in range(0, decades + 1):
            y = tb - tph * d / decades
            out.append(f"<line x1='{left}' x2='{width - 4:.0f}' y1='{y:.1f}' "
                       f"y2='{y:.1f}' stroke='{CHART_GRID if d else CHART_AXIS}' "
                       "stroke-width='1'/>")
            lab = {0: "1", 1: "10", 2: "100", 3: "1k", 4: "10k", 5: "100k",
                   6: "1M"}.get(d, f"1e{d}")
            if decades <= 4 or d % 2 == 0 or d == decades:
                out.append(f"<text x='{left - 6}' y='{y + 3.5:.1f}' "
                           f"text-anchor='end' font-size='10' "
                           f"fill='{CHART_INK}'>{lab}</text>")
        out.append(f"<text x='{left}' y='{t0 - 0.5:.1f}' font-size='10' "
                   f"fill='{CHART_INK}'>{esc(totals_label)} (log scale)</text>")
        for i, lab in enumerate(labels):
            v = totals[i] or 0
            x = left + slot * i + (slot - bar) / 2
            if v > 0:
                y = ty(v)
                out.append(f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar:.1f}' "
                           f"height='{max(tb - y, 0.5):.1f}' fill='#5f5e5a'>"
                           f"<title>{esc(lab)}: {fmt(v)} {esc(totals_label)}"
                           f"</title></rect>")
    for k in range(11):
        y = top + ph - ph * k / 10
        out.append(f"<line x1='{left}' x2='{width - 4:.0f}' y1='{y:.1f}' "
                   f"y2='{y:.1f}' stroke='{CHART_GRID if k else CHART_AXIS}' "
                   "stroke-width='1'/>")
        out.append(f"<text x='{left - 6}' y='{y + 3.5:.1f}' text-anchor='end' "
                   f"font-size='10.5' fill='{CHART_INK}'>{k * 10}%</text>")
    for i, lab in enumerate(labels):
        x = left + slot * i + (slot - bar) / 2
        tot = sum(vals[i] for _n, _c, vals in series)
        if tot <= 0:
            out.append(f"<rect x='{x:.1f}' y='{top}' width='{bar:.1f}' "
                       f"height='{ph}' fill='none' stroke='{CHART_GRID}'/>")
        else:
            y = top + ph
            for name, colour, vals in series:
                v = vals[i]
                if v <= 0:
                    continue
                h = ph * v / tot
                gap = 2 if h > 3 else 0
                out.append(
                    f"<rect x='{x:.1f}' y='{y - h:.1f}' width='{bar:.1f}' "
                    f"height='{max(h - gap, 0.5):.1f}' fill='{colour}'>"
                    f"<title>{esc(lab)} - {esc(name)}: {100 * v / tot:.1f}% "
                    f"({fmt(v)})</title></rect>")
                y -= h
        cx = x + bar / 2
        ly = top + ph + 10
        out.append(f"<text transform='translate({cx:.1f},{ly}) rotate(-55)' "
                   f"text-anchor='end' font-size='10.5' fill='{CHART_INK}'>"
                   f"{esc(lab)}</text>")
    out.append("</svg>")
    return "".join(out)
NEG = ("neg_extraction", "neg_ntc")


def pct(x) -> str:
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return "—"


def natural(sid: str):
    """Sort key: A2 before A10, plate wells in row then column order."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", sid)]


def legend(items) -> str:
    return ("<div class='legend' style='flex-wrap:wrap'>" + "".join(
        f"<span class='sw' style='background:{c}'></span>{esc(lab)}&nbsp;&nbsp;"
        for c, lab in items) + "</div>")


def build(args) -> str:
    d = args.synthesis_dir
    calls = read_table(os.path.join(d, "db_calls.tsv"))
    wells = sorted(read_table(os.path.join(d, "db_wells.tsv")),
                   key=lambda w: natural(w["sample_id"]))
    prevalence = read_table(os.path.join(d, "db_prevalence.tsv"))
    with open(os.path.join(d, "db_run_qc.json"), encoding="utf-8") as fh:
        qc = json.load(fh)
    run_id = args.run_id or qc.get("run_id", "")
    db = qc.get("database", {})
    neg = qc.get("negatives", {})
    P = []

    P.append(f"<h1>{esc(run_id)} — database mode</h1>")
    pr = db.get("primers", {})
    P.append(
        f"<p class='sub'>Classified against <b>{esc(db.get('name') or db.get('version'))}</b> "
        f"(release {esc(db.get('version'))}, built {esc(db.get('built', '')[:10])}) "
        f"with minimap2 {esc(qc.get('minimap2', ''))}. Primers "
        f"<code>{esc((pr.get('forward') or {}).get('sequence', ''))}</code> / "
        f"<code>{esc((pr.get('reverse') or {}).get('sequence', ''))}</code>.</p>")
    P.append(
        "<p><b>Every number here is relative 18S abundance</b> - a share of "
        "reads - and never organisms or copies. 18S copy number varies by "
        "orders of magnitude between taxa, and this protocol has no "
        "quantification standards to correct for it.</p>")

    det = [c for c in calls if c["call"] == "detected"]
    samples = [w for w in wells if w["sample_type"] not in NEG]
    P.append("<div class='cards'>" + "".join(
        f"<div class='card'><div class='n'>{v}</div><div class='k'>{k}</div></div>"
        for k, v in (("wells", len(wells)), ("blanks", neg.get("n", 0)),
                     ("helminth taxa", len({c["taxon"] for c in det
                                            if c["group"] == "helminths"})),
                     ("protist taxa", len({c["taxon"] for c in det
                                           if c["group"] == "protists"})),
                     ("calls above species", sum(
                         1 for c in det if c["group"] in TARGET
                         and c["rank"] != "species"))))
        + "</div>")
    P.append("<p class='note'>Sections: <a href='#funnel'>where the "
             "reads went</a> · <a href='#warn'>before reading the "
             "results</a> · <a href='#groups'>all reads by organism group</a> · "
             + " · ".join(f"<a href='#{g}'>{GROUP_LABEL[g]}</a>"
                          for g in TARGET)
             + " · <a href='#offtarget'>off-target</a> · <a href='#rank'>"
             "naming depth</a> · <a href='#calls'>every call</a></p>")

    # -- the read funnel, as in the panel report -------------------------
    # The per-well figures start at reads that reached a sample, already
    # net of splitting and demultiplexing; this is the one place a run
    # that lost its reads to barcodes, to primer dimer, or to bacteria
    # looks different from one that was simply shallow.
    funnel = qc.get("read_funnel") or {}
    stages = funnel.get("stages") or []
    if stages:
        P.append("<h2 id='funnel'>Where the reads went</h2>")
        P.append(
            "<p>Every read from the sequencer to a classified target "
            "organism, and what was lost at each step. The per-sample "
            "figures further down start at <i>demultiplexed to a sample</i>"
            ", so this is the only place a run that lost its reads to "
            "unreadable barcodes, to primer dimer or to off-target "
            "amplification looks different from one that was simply "
            "shallow.</p>")
        start = stages[0]["reads"] or 1
        P.append("<table><thead><tr><th>Stage</th><th class='num'>Reads"
                 "</th><th class='num'>% of start</th><th></th>"
                 "<th>Lost here</th></tr></thead><tbody>")
        for st in stages:
            p = st.get("percent_of_start")
            frac = min(1.0, st["reads"] / start)
            bar = (f'<div class="bar"><span style="width:{100 * frac:.2f}%;'
                   f'background:#3d6b8e"></span></div>')
            lost = "<br>".join(f"{esc(k)} {fmt(v)}"
                               for k, v in st["lost_here"].items()) or "—"
            P.append(
                f"<tr><td>{esc(st['stage'])}</td>"
                f"<td class='num'>{fmt(st['reads'])}</td>"
                f"<td class='num'>{'' if p is None else f'{p:.1f}%'}</td>"
                f"<td style='width:28%'>{bar}</td>"
                f"<td class='note'>{lost}</td></tr>")
        P.append("</tbody></table>")
        if funnel.get("handoff_balances") is False:
            h = funnel.get("handoff") or {}
            P.append(
                f"<p>{flag('bad', 'MISMATCH')} Demultiplexing wrote "
                f"{fmt(h.get('written_by_demultiplexing'))} reads to samples "
                f"but assignment read {fmt(h.get('read_by_assignment'))}. "
                "No filter accounts for the difference: data was lost "
                "between the steps.</p>")
        bad = funnel.get("wells_losing_over_half") or []
        if bad:
            P.append(
                f"<p>{flag('warn', 'wells')} {len(bad)} well(s) classified "
                "fewer than half the reads they received. They look shallow "
                "below, but the reaction rather than the loading is the "
                "likelier problem: "
                + esc(", ".join(
                    f"{b['sample_id']} ({100 * b['classified_fraction']:.0f}%)"
                    for b in sorted(bad, key=lambda b: natural(b["sample_id"]))))
                + ".</p>")

    # -- what could undermine the results -------------------------------
    P.append("<h2 id='warn'>Before reading the results</h2>")
    items = []
    if not neg.get("n"):
        items.append(flag("bad", "no blanks") + " No negative controls on this "
                     "run. With no plasmids or standards they are the only "
                     "measure of contamination and cross-talk, so every call "
                     "below rests on a fixed floor that measures nothing.")
    elif neg.get("taxa_in_blanks"):
        items.append(flag("warn", "blanks not empty") + " The blanks contain "
                     + esc(", ".join(neg["taxa_in_blanks"]))
                     + ". Calls of these taxa elsewhere must beat what the "
                     "blanks show; see <i>limit</i> in the call table.")
    else:
        items.append(flag("ok", "blanks clean") + " Nothing was classified "
                     "in the blanks.")
    if qc.get("low_depth_wells"):
        items.append(flag("warn", "low depth") + " Too few classified reads "
                     "to call anything in: "
                     + esc(", ".join(qc["low_depth_wells"])) + ".")
    if qc.get("empty_wells"):
        items.append(flag("warn", "empty") + " No reads demultiplexed to: "
                     + esc(", ".join(qc["empty_wells"])) + ".")
    ia = qc.get("internal_anchors", {})
    if ia.get("reads_in"):
        items.append(
            f"{flag('ok', 'adapters')} The splitter cut "
            f"{fmt(ia.get('reads_cut_by_splitter'))} of {fmt(ia['reads_in'])} "
            f"reads at internal adapters, and "
            f"{fmt(ia.get('quarantined_internal_adapter') or 0)} reads were "
            "set aside as more than one molecule. The release's look-alike "
            f"mask ({fmt(ia.get('mask_rows'))} sites) protects taxa in the "
            "database from both; a taxon missing from it is not protected, "
            "so a high rate here is worth a second look.")
    chq = qc.get("chimeras") or {}
    if chq.get("reads"):
        items.append(
            flag("warn" if chq.get("wells_over_5pct") else "ok", "chimeras")
            + f" {fmt(chq['reads'])} reads ({100 * chq['of_clustered']:.1f}% "
            "of classified) were PCR chimeras - joins of two more abundant "
            "sequences in the same well. They are counted but never named, "
            "and kept out of every composition."
            + (" Over 5% in: " + esc(", ".join(chq["wells_over_5pct"])) + "."
               if chq.get("wells_over_5pct") else ""))
    near = qc.get("near_3prime_calls") or []
    if near:
        items.append(flag("warn", "primer fit") + " Detected with a best hit "
                     "whose primer sites mismatch in the last five bases, so "
                     "their share may be held down by the PCR, not the "
                     "sample: " + esc("; ".join(near[:20]))
                     + (" …" if len(near) > 20 else "") + ".")
    P.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>")

    # -- all reads, by group (VESPA Figs 3a, 4a) --------------------------
    present = [g for g, _c, _l in GROUP_STYLE
               if any(int(w.get(f"{g}_reads") or 0) for w in wells)]
    P.append("<h2 id='groups'>All reads, by organism group</h2>")
    P.append(
        "<p>Each well's classified reads split by what they are, target or "
        "not - the first thing to know about a well before its helminths "
        "and protists are read off. The last column is the read funnel: "
        "how many of the reads that reached the sample were classified.</p>")
    P.append(legend([(GROUP_COLOUR[g], GROUP_LABEL[g]) for g in present]))
    P.append(columns_svg(
        [w["sample_id"] for w in wells],
        [(GROUP_LABEL[g], GROUP_COLOUR[g],
          [int(w.get(f"{g}_reads") or 0) for w in wells]) for g in present],
        note="Share of each well's classified reads by organism group",
        totals=[int(w.get("clustered") or 0) for w in wells],
        totals_label="classified reads"))
    rows = []
    for w in wells:
        funnel = [(lab, int(w.get(k) or 0), c) for k, c, lab in FUNNEL]
        rows.append(
            f"<tr><td>{esc(w['sample_id'])}<br><span class='note'>"
            f"{esc(w['sample_type'])}</span></td>"
            f"<td class='num'>{fmt(w['clustered'])}</td>"
            + "".join(f"<td class='num'>{pct(w.get(g + '_fraction'))}</td>"
                      for g in TARGET)
            + f"<td class='num'>{fmt(w['reads_in'])}</td>"
            f"<td style='min-width:80px'>{stacked_bar(funnel)}</td></tr>")
    P.append(
        "<div class='scroll'><table><thead><tr><th>sample</th>"
        "<th class='num'>classified reads</th>"
        + "".join(f"<th class='num'>{GROUP_LABEL[g]}</th>" for g in TARGET)
        + "<th class='num'>reads in</th><th>funnel</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>")
    P.append(legend([(c, f"funnel: {lab}") for _k, c, lab in FUNNEL]))

    # -- each target group on its own (VESPA Figs 4b-d, 5b-d) -------------
    by_sample = defaultdict(list)
    for c in calls:
        by_sample[c["sample_id"]].append(c)
    for g in TARGET:
        gcalls = [c for c in det if c["group"] == g and c["sample_type"] not in NEG]
        P.append(f"<h2 id='{g}'>{GROUP_LABEL[g].capitalize()}</h2>")
        if not any(c["group"] == g for c in calls):
            P.append(f"<p class='note'>No {GROUP_LABEL[g]} were classified "
                     "on this run.</p>")
            continue
        P.append(
            f"<p>Shares <b>within {GROUP_LABEL[g]}</b>: the other groups are "
            f"set aside and each well's {GROUP_LABEL[g]} reads are the "
            "denominator, so a well that is mostly host still shows its "
            f"{GROUP_LABEL[g]} composition. Few reads under a column means "
            "little confidence in it.</p>")
        totals = defaultdict(float)
        for c in gcalls:
            totals[(c["rank"], c["taxon"])] += float(c["group_share"] or 0)
        top = [k for k, _v in sorted(totals.items(), key=lambda kv: -kv[1])]
        colours = {k: WITHIN_COLOURS[i]
                   for i, k in enumerate(top[:len(WITHIN_COLOURS)])}
        # Legend by the readable label ("Chromadorea, unnamed family
        # (nearest Haemonchus)"), not PR2's placeholder name.
        label_of = {}
        for c in gcalls:
            label_of.setdefault((c["rank"], c["taxon"]), c["label"])
        def name(k):
            lab = label_of.get(k, k[1])
            return lab if lab != k[1] else f"{k[1]} ({k[0]})"
        P.append(legend([(colours[k], name(k)) for k in colours]
                        + ([(OTHER_COLOUR, "everything else")]
                           if len(top) > len(colours) else [])))
        # The chart: one column per sample, its detected taxa as shares of
        # the group. Taxa past the palette pool as "everything else".
        named = list(colours)
        series = []
        for k in named:
            series.append((name(k), colours[k], [
                sum(float(c["group_share"] or 0) for c in by_sample.get(w["sample_id"], [])
                    if c["group"] == g and c["call"] == "detected"
                    and (c["rank"], c["taxon"]) == k) for w in samples]))
        series.append(("everything else", OTHER_COLOUR, [
            sum(float(c["group_share"] or 0) for c in by_sample.get(w["sample_id"], [])
                if c["group"] == g and c["call"] == "detected"
                and (c["rank"], c["taxon"]) not in colours) for w in samples]))
        P.append(columns_svg([w["sample_id"] for w in samples], series,
                             note=f"Composition within {GROUP_LABEL[g]}, per sample",
                             totals=[int(w.get(f"{g}_reads") or 0) for w in samples],
                             totals_label=f"{GROUP_LABEL[g]} reads"))
        rr = []
        for w in samples:
            n = int(w.get(f"{g}_reads") or 0)
            rr.append(f"<tr><td>{esc(w['sample_id'])}</td>"
                      f"<td class='num'>{fmt(n)}</td>"
                      f"<td class='num'>{esc(w.get(g + '_richness', ''))}</td></tr>")
        P.append("<div class='scroll' style='max-height:240px'><table><thead>"
                 "<tr><th>sample</th>"
                 f"<th class='num'>{GROUP_LABEL[g]} reads</th>"
                 "<th class='num'>richness</th></tr></thead><tbody>"
                 + "".join(rr) + "</tbody></table></div>")
        P.append("<p class='note'>An empty column is a sample with no "
                 f"{GROUP_LABEL[g]} detected. The read count is what each "
                 "column rests on.</p>")

        # presence table: taxa x samples, with prevalence
        taxa = [r for r in prevalence if r["group"] == g]
        if taxa:
            sids = [w["sample_id"] for w in samples]
            look = {(c["sample_id"], c["rank"], c["taxon"]): c for c in calls
                    if c["group"] == g}
            head = "".join(f"<th class='num' style='writing-mode:vertical-rl;"
                           f"transform:rotate(180deg)'>{esc(s)}</th>"
                           for s in sids)
            body = []
            for t in taxa:
                cells = []
                for sid in sids:
                    c = look.get((sid, t["rank"], t["taxon"]))
                    if c is None:
                        cells.append("<td></td>")
                    elif c["call"] == "detected":
                        cells.append(
                            f"<td class='num' style='background:#e6f4ec' "
                            f"title='{esc(c['reads'])} reads'>"
                            f"{100 * float(c['group_share'] or 0):.0f}</td>")
                    else:
                        cells.append(f"<td class='num' title='{esc(c['call'])}'"
                                     " style='color:var(--muted)'>·</td>")
                body.append(
                    f"<tr><td>{esc(t['label'])}</td><td>{esc(t['rank'])}</td>"
                    f"<td class='num'>{esc(t['samples_detected'])}/"
                    f"{esc(t['samples_tested'])}</td>{''.join(cells)}</tr>")
            P.append(
                "<h3>Presence by sample</h3><div class='scroll'><table><thead>"
                "<tr><th>taxon</th><th>rank</th><th class='num'>prevalence</th>"
                + head + "</tr></thead><tbody>" + "".join(body)
                + "</tbody></table></div>")
            P.append(
                "<p class='note'>A number is the taxon's percentage within "
                f"the well's {GROUP_LABEL[g]}, where it was detected; a dot "
                "means it was seen but not called (below the blanks or the "
                "relative floor; hover for which). Prevalence counts samples "
                "deep enough to call anything, blanks excluded.</p>")

    # -- off-target --------------------------------------------------------
    off = [g for g, _c, _l in GROUP_STYLE if g not in TARGET]
    P.append("<h2 id='offtarget'>Off-target</h2>")
    P.append(
        "<p>What else the primers amplified, as shares of each well's "
        "classified reads, from the PR2 lineage each call reached; host is "
        f"<code>{esc(qc.get('settings', {}).get('host_clade', ''))}</code>. "
        "<i>Unassigned</i> stopped too high in the tree to place; "
        "<i>unclassified</i> is consensus that cleared no identity floor down "
        "to family - most often a taxon the database does not hold.</p>")
    rows = []
    for w in wells:
        rows.append(f"<tr><td>{esc(w['sample_id'])}</td>"
                    + "".join(f"<td class='num'>{pct(w.get(g + '_fraction'))}</td>"
                              for g in off)
                    + f"<td class='num'>{pct(w['unclassified_fraction'])}</td></tr>")
    P.append("<div class='scroll'><table><thead><tr><th>sample</th>"
             + "".join(f"<th class='num'>{GROUP_LABEL[g]}</th>" for g in off)
             + "<th class='num'>unclassified</th></tr></thead><tbody>"
             + "".join(rows) + "</tbody></table></div>")

    # -- rank reached -----------------------------------------------------
    P.append("<h2 id='rank'>How deep each well could be named</h2>")
    P.append(legend([(c, k) for k, c in RANK_COLOURS.items()]))
    P.append(columns_svg(
        [w["sample_id"] for w in samples],
        [(k, c, [int(w.get(f"rank_{k}") or 0) for w in samples])
         for k, c in RANK_COLOURS.items()],
        note="Share of each sample's classified reads by the rank reached",
        totals=[int(w.get("clustered") or 0) for w in samples],
        totals_label="classified reads"))
    fl = (qc.get("settings") or {}).get("floors") or {}
    P.append(
        "<p class='note'>A consensus is named at a rank only when the lowest "
        "common ancestor of its near-best database hits reaches it AND its "
        "best hit clears that rank's identity floor ("
        + esc(", ".join(f"{k} {v}" for k, v in fl.items()))
        + "). A genus call is not a failed species call: it says the "
        "database holds a close relative but not this sequence, or that V4 "
        "cannot tell the species apart.</p>")

    # -- the calls ---------------------------------------------------------
    P.append("<h2 id='calls'>Every call</h2>")
    lev = {"detected": "ok", "below_blanks": "warn", "below_floor": "warn",
           "low_depth": "bad", "chimera": "warn"}
    order = {g: i for i, (g, _c, _l) in enumerate(GROUP_STYLE)}
    cr = []
    for c in sorted(calls, key=lambda c: (c["sample_id"],
                                          order.get(c["group"], 99),
                                          -int(c["reads"]))):
        notes = "primer fit" if "near_3prime" in c["primer_fit"] else ""
        cr.append(
            f"<tr><td>{esc(c['sample_id'])}</td>"
            f"<td><span class='sw' style='display:inline-block;background:"
            f"{GROUP_COLOUR.get(c['group'], OTHER_COLOUR)}'></span> "
            f"{esc(GROUP_LABEL.get(c['group'], c['group']))}</td>"
            f"<td>{esc(c['rank'] or '—')}</td>"
            f"<td>{esc(c['label'] or '—')}</td>"
            f"<td class='num'>{fmt(c['reads'])}</td>"
            f"<td class='num'>{pct(c['rel_abundance'])}</td>"
            f"<td class='num'>{pct(c['group_share'])}</td>"
            f"<td class='num'>{100 * float(c['identity']):.1f}%</td>"
            f"<td>{flag(lev.get(c['call'], 'warn'), c['call'])}</td>"
            f"<td class='num' title='{esc(c['limit_basis'])}'>{fmt(c['limit_reads'])}</td>"
            f"<td>{esc(notes)}</td></tr>")
    P.append("<div class='scroll'><table><thead><tr><th>sample</th><th>group</th>"
             "<th>rank</th><th>taxon</th><th class='num'>reads</th>"
             "<th class='num'>of well</th><th class='num'>of group</th>"
             "<th class='num'>identity</th><th>call</th><th class='num'>limit</th>"
             "<th>notes</th></tr></thead><tbody>" + "".join(cr)
             + "</tbody></table></div>")
    P.append(
        "<p class='note'><i>Identity</i> is the best database hit's, read "
        "beside the rank: a species call at 99.6% and one at 100% are not the "
        "same claim. <i>Limit</i> is the fewest reads that beat the blanks "
        "(hover for how it was set). The cluster sequences, best hits and "
        "full lineages are in <code>assign_db/&lt;sample&gt;.clusters.tsv</code>; "
        "this table is <code>synthesis/db_calls.tsv</code>.</p>")

    body = "\n".join(P)
    return ("<!doctype html>\n<html lang='en'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{esc(run_id)} — amphid (database mode)</title>"
            f"<style>{CSS}</style></head><body><div class='wrap'>{body}"
            "</div></body></html>\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--synthesis-dir", required=True)
    ap.add_argument("--run-id", default="")
    ap.add_argument("--out", default="report.html")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    doc = build(args)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    if not args.quiet:
        print(f"{args.out}  ({len(doc) / 1024:.0f} KB, self-contained)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
