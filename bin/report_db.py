#!/usr/bin/env python3
"""
Run report for database mode (--reference_database).

The same rules as the panel report (report.py), whose style and helpers
it shares: one self-contained HTML file, no network, no JavaScript
libraries, nothing computed that the tables do not already hold, and the
things that would undermine the results placed before the results.

What it puts first is different, because what can undermine a database
run is different. There are no standards and no plasmids, so the
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
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table  # noqa: E402
from report import (  # noqa: E402
    CSS, OTHER_COLOUR, TAXON_COLOURS, esc, flag, fmt, stacked_bar,
)

FUNNEL = (("clustered", "#3d6b8e", "clustered - classified"),
          ("unclustered", "#9fb4c4", "unclustered - too few reads for a consensus"),
          ("too_short", "#c98a5e", "too short"),
          ("no_primer", "#c9c2b0", "no primer pair"))
OFF = (("bacteria_archaea", "#8c6d31", "bacteria & archaea"),
       ("organelle", "#b5cf6b", "organelle (plastid, mitochondrion)"),
       ("fungi", "#9e9ac8", "fungi"),
       ("host", "#d6616b", "host"))
RANK_COLOURS = {"species": "#0b6b3a", "genus": "#3d8f5f",
                "family": "#8a5300", "above": "#b03030"}


def pct(x) -> str:
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return "—"


def build(args) -> str:
    d = args.synthesis_dir
    calls = read_table(os.path.join(d, "db_calls.tsv"))
    wells = read_table(os.path.join(d, "db_wells.tsv"))
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
        "<p><b>Every number here is relative 18S abundance</b> - the share of "
        "a well's classified reads - and never organisms or copies. 18S copy "
        "number varies by orders of magnitude between taxa, and this "
        "protocol has no quantification standards to correct for it.</p>")

    det = [c for c in calls if c["call"] == "detected"]
    samples = [w for w in wells if w["sample_type"] not in
               ("neg_extraction", "neg_ntc")]
    P.append("<div class='cards'>" + "".join(
        f"<div class='card'><div class='n'>{v}</div><div class='k'>{k}</div></div>"
        for k, v in (("wells", len(wells)), ("blanks", neg.get("n", 0)),
                     ("detected calls", len(det)),
                     ("taxa at species", sum(1 for c in det
                                             if c["rank"] == "species")),
                     ("named above species", sum(1 for c in det
                                                 if c["rank"] != "species"))))
        + "</div>")

    # -- what could undermine the results -------------------------------
    P.append("<h2>Before reading the results</h2>")
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
                     "blanks show; see <i>limit</i> in the table.")
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
    near = qc.get("near_3prime_calls") or []
    if near:
        items.append(flag("warn", "primer fit") + " Detected with a best hit "
                     "whose primer sites mismatch in the last five bases, so "
                     "their share may be held down by the PCR, not the "
                     "sample: " + esc("; ".join(near[:20]))
                     + (" …" if len(near) > 20 else "") + ".")
    P.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>")

    # -- reads and off-target -------------------------------------------
    P.append("<h2>Where each well's reads went</h2>")
    P.append("<div class='legend'>" + "".join(
        f"<span class='sw' style='background:{c}'></span>{esc(lab)}&nbsp;&nbsp;"
        for _k, c, lab in FUNNEL) + "</div>")
    rows = []
    for w in wells:
        parts = [(lab, int(w.get(k) or 0), c) for k, c, lab in FUNNEL]
        rows.append(
            f"<tr><td>{esc(w['sample_id'])}</td><td>{esc(w['sample_type'])}</td>"
            f"<td class='num'>{fmt(w['reads_in'])}</td>"
            f"<td style='min-width:220px'>{stacked_bar(parts)}</td>"
            f"<td class='num'>{fmt(w['clustered'])}</td>"
            f"<td class='num'>{fmt(w['clusters'])}</td>"
            + "".join(f"<td class='num'>{pct(w.get(k + '_fraction'))}</td>"
                      for k, _c, _l in OFF)
            + f"<td class='num'>{pct(w['unclassified_fraction'])}</td></tr>")
    P.append(
        "<div class='scroll'><table><thead><tr><th>sample</th><th>type</th>"
        "<th class='num'>reads</th><th>fate</th><th class='num'>classified</th>"
        "<th class='num'>clusters</th>"
        + "".join(f"<th class='num'>{esc(lab.split(' (')[0])}</th>"
                  for _k, _c, lab in OFF)
        + "<th class='num'>unclassified</th></tr></thead><tbody>"
        + "".join(rows) + "</tbody></table></div>")
    P.append(
        "<p class='note'>Off-target columns are shares of each well's "
        "classified reads, from the PR2 lineage of each call's best hit; "
        f"host is the clade <code>{esc(qc.get('settings', {}).get('host_clade', ''))}</code>. "
        "<i>Unclassified</i> is consensus that cleared no identity floor down "
        "to family - most often a taxon the database does not hold.</p>")

    # -- rank reached -----------------------------------------------------
    P.append("<h2>How deep each well could be named</h2>")
    P.append("<div class='legend'>" + "".join(
        f"<span class='sw' style='background:{c}'></span>{k}&nbsp;&nbsp;"
        for k, c in RANK_COLOURS.items()) + "</div>")
    rr = []
    for w in samples:
        parts = [(k, int(w.get(f"rank_{k}") or 0), c)
                 for k, c in RANK_COLOURS.items()]
        rr.append(f"<tr><td>{esc(w['sample_id'])}</td>"
                  f"<td style='min-width:320px'>{stacked_bar(parts)}</td>"
                  + "".join(f"<td class='num'>{fmt(w.get('rank_' + k))}</td>"
                            for k in RANK_COLOURS) + "</tr>")
    P.append("<table><thead><tr><th>sample</th><th>reads by rank reached</th>"
             + "".join(f"<th class='num'>{k}</th>" for k in RANK_COLOURS)
             + "</tr></thead><tbody>" + "".join(rr) + "</tbody></table>")
    fl = (qc.get("settings") or {}).get("floors") or {}
    P.append(
        "<p class='note'>A consensus is named at a rank only when the lowest "
        "common ancestor of its near-best database hits reaches it AND its "
        "best hit clears that rank's identity floor ("
        + esc(", ".join(f"{k} {v}" for k, v in fl.items()))
        + "). A genus call is not a failed species call: it says the "
        "database holds a close relative but not this sequence, or that V4 "
        "cannot tell the species apart.</p>")

    # -- composition --------------------------------------------------------
    P.append("<h2>Composition: relative 18S abundance</h2>")
    totals = defaultdict(float)
    for c in det:
        if c["sample_type"] not in ("neg_extraction", "neg_ntc"):
            totals[(c["rank"], c["taxon"])] += float(c["rel_abundance"])
    top = [k for k, _v in sorted(totals.items(), key=lambda kv: -kv[1])]
    colours = {k: TAXON_COLOURS[i] for i, k in enumerate(top[:len(TAXON_COLOURS)])}
    P.append("<div class='legend' style='flex-wrap:wrap'>" + "".join(
        f"<span class='sw' style='background:{colours[k]}'></span>"
        f"{esc(k[1])} <span class='note'>({esc(k[0])})</span>&nbsp;&nbsp;"
        for k in colours) + f"<span class='sw' style='background:{OTHER_COLOUR}'>"
        "</span>everything else</div>")
    by_sample = defaultdict(list)
    for c in calls:
        by_sample[c["sample_id"]].append(c)
    cr = []
    for w in samples:
        mine = [c for c in by_sample.get(w["sample_id"], [])
                if c["call"] == "detected"]
        parts = [(f"{c['taxon']} ({c['rank']})", float(c["rel_abundance"]),
                  colours.get((c["rank"], c["taxon"]), OTHER_COLOUR))
                 for c in mine]
        cr.append(f"<tr><td>{esc(w['sample_id'])}</td>"
                  f"<td style='min-width:420px'>{stacked_bar(parts)}</td>"
                  f"<td class='num'>{len(mine)}</td></tr>")
    P.append("<table><thead><tr><th>sample</th><th>detected taxa</th>"
             "<th class='num'>n</th></tr></thead><tbody>" + "".join(cr)
             + "</tbody></table>")

    # -- the calls ---------------------------------------------------------
    P.append("<h2>Every call</h2>")
    lev = {"detected": "ok", "below_blanks": "warn", "below_floor": "warn",
           "low_depth": "bad"}
    cr = []
    for c in calls:
        tags = []
        if c["offtarget"]:
            tags.append(c["offtarget"].replace("_", " "))
        if "near_3prime" in c["primer_fit"]:
            tags.append("primer fit")
        cr.append(
            f"<tr><td>{esc(c['sample_id'])}</td><td>{esc(c['rank'] or '—')}</td>"
            f"<td>{esc(c['label'] or '—')}</td>"
            f"<td class='num'>{fmt(c['reads'])}</td>"
            f"<td class='num'>{pct(c['rel_abundance'])}</td>"
            f"<td class='num'>{100 * float(c['identity']):.1f}%</td>"
            f"<td>{flag(lev.get(c['call'], 'warn'), c['call'])}</td>"
            f"<td class='num' title='{esc(c['limit_basis'])}'>{fmt(c['limit_reads'])}</td>"
            f"<td>{esc(', '.join(tags))}</td></tr>")
    P.append("<div class='scroll'><table><thead><tr><th>sample</th><th>rank</th>"
             "<th>taxon</th><th class='num'>reads</th><th class='num'>share</th>"
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
            f"<title>{esc(run_id)} — nf-helminth-amplicon (database mode)</title>"
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
