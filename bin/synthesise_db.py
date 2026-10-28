#!/usr/bin/env python3
"""
Run-level synthesis for database mode (--reference_database).

Panel synthesis (synthesise.py) is built on what this protocol does not
have: quantification standards to turn shares into copies, and tagged
control plasmids to measure carryover. So database mode gets its own,
smaller synthesis, and what it can claim is narrower:

RELATIVE 18S ABUNDANCE, NOT ORGANISMS OR COPIES
    A read share is the share of 18S template, and 18S copy number varies
    widely between taxa - by orders of magnitude between a ciliate and a
    nematode. Every number here is labelled as what it is.

AT THE RANK EACH CALL REACHED
    Calls come from bin/assign_db.py at species, genus, family or above,
    and stay there. A species call and a genus call in the same genus
    are different statements and are kept as separate rows; nothing is
    rolled up or pushed down to make the table tidier.

JUDGED AGAINST THE BLANKS AND A RELATIVE FLOOR
    With no plasmids and no standards the negative controls are the only
    measure left of what an empty well shows. The detection limit is the
    panel's, per taxon from the run's own blanks: an absolute bound (the
    most reads the taxon reached in any blank - hopping-like) and a
    proportional one (its largest share of any blank deep enough to have
    a share - carryover-like), whichever is stricter, and never below a
    fixed floor. Then a relative-abundance floor, because a call that
    clears the blanks at 0.01% of a well is still one read in ten
    thousand. Without blanks the limit is the fixed floor, and every call
    says it is an assumption.

ORGANISM GROUPS
    Every call is put in one group, the way VESPA reports (Owens et al.
    2024, Figs 3-5): first every well as shares of all its classified
    reads by group - what the primers amplified, target or not - then
    each target group on its own, shares recalculated within the group,
    so a helminth composition is not squashed flat by a well that is 90%
    host. The groups, from the PR2 lineage down to the rank the call
    REACHED (a call that stopped at domain cannot be placed in a group,
    and is 'unassigned' rather than guessed), first match wins:

      bacteria_archaea  domain Bacteria or Archaea
      organelle         PR2's organelle domains (Eukaryota:plas, :mito...)
      host              --host-clade anywhere (Craniata: vertebrates)
      helminths         class Nematoda, Platyhelminthes or Acanthocephala
      microsporidia     order Microsporida
      fungi             subdivision Fungi
      invertebrates     any other Metazoa (arthropods, mites, annelids...)
      plants            supergroup Archaeplastida (land plants and algae)
      protists          any other eukaryote - Apicomplexa, Amoebozoa,
                        Metamonada, Blastocystis, ciliates...
      unassigned        not placed deep enough to say

    helminths, protists and microsporidia are the TARGET groups - what
    the protocol is for; the rest are off-target, reported as fractions
    of each well, as VESPA does.

PREDICTED AMPLIFICATION
    The release's primer_fit.tsv says how well the primers fit each
    record. A call whose best hit is classed near_3prime - a mismatch in
    the primer's last five bases - is flagged: its signal may be low
    because of the PCR, not the sample.

Standard library only.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    RANKS, load_release_manifest, read_table, release_version,
    restore_declared_ids, writer,
)

NEGATIVE_TYPES = {"neg_extraction", "neg_ntc"}

CALLS = ("detected", "below_blanks", "below_floor", "low_depth")
CALL_MEANING = {
    "detected": "above the run's blank-derived limit and the relative "
                "abundance floor",
    "below_blanks": "no more reads than the blanks show for this taxon",
    "below_floor": "above the blanks, but under the relative abundance floor",
    "low_depth": "the well has fewer clustered reads than --min-depth, so "
                 "nothing in it is called",
}

GROUPS = ("helminths", "protists", "microsporidia", "host", "fungi",
          "plants", "invertebrates", "bacteria_archaea", "organelle",
          "unassigned")
TARGET_GROUPS = ("helminths", "protists", "microsporidia")
OFFTARGET = tuple(g for g in GROUPS if g not in TARGET_GROUPS
                  and g != "unassigned")
HELMINTH_CLASSES = {"Nematoda", "Platyhelminthes", "Acanthocephala"}


def organism_group(lineage: list[str], host_clade: str) -> str:
    """The reporting group of a call, from its lineage down to the rank
    it reached (module docstring). 'unassigned' when that lineage stops
    too high to tell."""
    lin = lineage + [""] * (len(RANKS) - len(lineage))
    dom, sup, sub, cls, order = lin[0], lin[1], lin[3], lin[4], lin[5]
    if dom in ("Bacteria", "Archaea"):
        return "bacteria_archaea"
    if dom.startswith("Eukaryota:"):
        return "organelle"
    if dom != "Eukaryota":
        return "unassigned"
    if host_clade and host_clade in lin:
        return "host"
    if cls in HELMINTH_CLASSES:
        return "helminths"
    if order == "Microsporida":
        return "microsporidia"
    if sub == "Fungi":
        return "fungi"
    if sub == "Metazoa":
        # a Metazoa call that stopped above class could still be a
        # helminth: not placed
        return "invertebrates" if cls else "unassigned"
    if sup == "Archaeplastida":
        return "plants"
    if sup == "Obazoa":
        # Choanoflagellates, apusomonads and the like are protists; an
        # Obazoa call that stopped above subdivision could be an animal.
        return "protists" if sub else "unassigned"
    return "protists" if sup else "unassigned"


class Well:
    def __init__(self, sid: str):
        self.sid = sid
        self.sample_type = ""
        self.clusters: list[dict] = []
        self.summary: dict = {}

    @property
    def depth(self) -> int:
        return int(self.summary.get("reads", {}).get("clustered", 0))


def load_wells(assign_dir: str) -> dict[str, Well]:
    wells: dict[str, Well] = {}
    for path in sorted(glob.glob(os.path.join(assign_dir, "*.assign_db.json"))):
        with open(path, encoding="utf-8") as fh:
            js = json.load(fh)
        w = wells.setdefault(js["sample_id"], Well(js["sample_id"]))
        w.summary = js
    for path in sorted(glob.glob(os.path.join(assign_dir, "*.clusters.tsv"))):
        rows = read_table(path)
        sid = os.path.basename(path)[:-len(".clusters.tsv")]
        wells.setdefault(sid, Well(sid)).clusters = rows
    return wells


def taxon_rows(w: Well, host_clade: str, primer_fit: dict) -> list[dict]:
    """One row per (rank, taxon, label) in a well, clusters pooled."""
    agg: dict = {}
    for c in w.clusters:
        key = (c["rank"], c["taxon"], c["label"])
        a = agg.setdefault(key, {"reads": 0, "clusters": 0, "identity": 0.0,
                                 "lineage": c["lineage"], "best_hits": []})
        a["reads"] += int(c["reads"])
        a["clusters"] += 1
        a["identity"] = max(a["identity"], float(c["identity"] or 0))
        a["best_hits"].append(c["best_hit"])
    out = []
    for (rank, taxon, label), a in agg.items():
        fits = sorted({primer_fit.get(h, "") for h in a["best_hits"]} - {""})
        out.append({"rank": rank, "taxon": taxon, "label": label,
                    "reads": a["reads"], "clusters": a["clusters"],
                    "identity": a["identity"], "lineage": a["lineage"],
                    "group": organism_group(
                        [x for x in a["lineage"].split(";")] if a["lineage"]
                        else [], host_clade),
                    "primer_fit": ",".join(fits)})
    return out


def blank_limits(wells: dict, rows: dict, min_fraction_depth: int) -> dict:
    """Per (rank, taxon): the most reads, and the largest share, any blank
    showed. Shares only from blanks deep enough for a share to mean
    anything - an NTC of three reads with one stray is not a 33%
    carryover rate (the panel synthesis learnt that the hard way)."""
    blanks = [w for w in wells.values() if w.sample_type in NEGATIVE_TYPES]
    out: dict = defaultdict(lambda: {"reads": 0, "fraction": 0.0})
    for w in blanks:
        for r in rows[w.sid]:
            k = (r["rank"], r["taxon"])
            out[k]["reads"] = max(out[k]["reads"], r["reads"])
            if w.depth >= min_fraction_depth and w.depth:
                out[k]["fraction"] = max(out[k]["fraction"],
                                         r["reads"] / w.depth)
    return out, len(blanks)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--assign-dir", required=True,
                    help="assign_db outputs: *.clusters.tsv, *.assign_db.json")
    ap.add_argument("--preflight-dir", required=True)
    ap.add_argument("--database", required=True,
                    help="the release directory (for primer_fit.tsv)")
    ap.add_argument("--demux-summary")
    ap.add_argument("--split-summary")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--min-depth", type=int, default=100,
                    help="clustered reads a well needs before anything in "
                    "it is called")
    ap.add_argument("--absolute-floor", type=int, default=2,
                    help="reads a call needs whatever the blanks show")
    ap.add_argument("--min-rel-abundance", type=float, default=0.001,
                    help="share of a well's clustered reads a call needs")
    ap.add_argument("--host-clade", default="Craniata",
                    help="a name anywhere in the lineage that marks host "
                    "reads (PR2: Craniata = vertebrates)")
    ap.add_argument("--outdir", default="synthesis")
    args = ap.parse_args(argv)

    wells = load_wells(args.assign_dir)
    restore_declared_ids(wells, args.preflight_dir, args.run_id)
    for sid, w in wells.items():
        w.sid = sid
    declared = [r for r in read_table(os.path.join(args.preflight_dir,
                                                   "file_resolution.tsv"))
                if not args.run_id or r.get("run_id") == args.run_id]
    types = {r["sample_id"]: r.get("sample_type", "") for r in declared}
    for sid, w in wells.items():
        w.sample_type = types.get(sid, "")

    # Every declared sample must arrive, or the blanks this synthesis
    # leans on might be the ones that went missing. A well the
    # demultiplexer wrote no reads for has nothing to assign and is
    # reported as empty rather than lost.
    per_sample: dict = {}
    demux: dict = {}
    if args.demux_summary and os.path.exists(args.demux_summary):
        with open(args.demux_summary, encoding="utf-8") as fh:
            demux = json.load(fh)
        per_sample = demux.get("per_sample") or {}
    missing = [r for r in declared if r["sample_id"] not in wells]
    lost = [r["sample_id"] for r in missing
            if per_sample.get(r["sample_id"], 0) or
            per_sample.get(r.get("sanitised_id", ""), 0)]
    if lost:
        raise SystemExit(
            f"cohort incomplete: {', '.join(sorted(lost))} had reads but no "
            "assign_db output. Refusing to report on a partial run")
    empty = sorted(r["sample_id"] for r in missing)
    for sid in empty:
        w = wells[sid] = Well(sid)
        w.sample_type = types.get(sid, "")

    man = load_release_manifest(args.database)
    version = release_version(man)
    pf_path = os.path.join(args.database, "primer_fit.tsv")
    primer_fit = ({r["id"]: r.get("fit", "") for r in read_table(pf_path)}
                  if os.path.exists(pf_path) else {})
    with open(os.path.join(args.preflight_dir, "preflight_summary.json"),
              encoding="utf-8") as fh:
        pf_summary = json.load(fh)

    rows = {sid: taxon_rows(w, args.host_clade, primer_fit)
            for sid, w in wells.items()}
    limits, n_blanks = blank_limits(wells, rows, args.min_depth)

    os.makedirs(args.outdir, exist_ok=True)
    call_cols = ["run_id", "sample_id", "sample_type", "rank", "taxon",
                 "label", "reads", "rel_abundance", "call", "limit_reads",
                 "limit_basis", "identity", "clusters", "lineage",
                 "group", "group_reads", "group_share", "primer_fit",
                 "database_version"]
    calls = []
    for sid in sorted(wells):
        w = wells[sid]
        in_group = Counter()
        for r in rows[sid]:
            in_group[r["group"]] += r["reads"]
        for r in sorted(rows[sid], key=lambda r: -r["reads"]):
            rel = r["reads"] / w.depth if w.depth else 0.0
            lim = limits.get((r["rank"], r["taxon"]))
            if not n_blanks:
                limit, basis = args.absolute_floor, (
                    "fixed floor - no blanks on this run, so the limit is "
                    "an assumption, not a measurement")
            else:
                absolute = lim["reads"] if lim else 0
                prop = int(math.ceil((lim["fraction"] if lim else 0)
                                     * w.depth))
                limit = max(absolute + 1 if absolute else 0, prop + 1
                            if prop else 0, args.absolute_floor)
                basis = ("absolute: reads in a blank" if absolute + 1 == limit
                         and absolute else
                         "proportional: share of a blank" if prop + 1 == limit
                         and prop else "fixed floor - blanks were clean "
                         "for this taxon")
            if w.depth < args.min_depth:
                call = "low_depth"
            elif r["reads"] < limit:
                call = "below_blanks"
            elif rel < args.min_rel_abundance:
                call = "below_floor"
            else:
                call = "detected"
            calls.append({
                "run_id": args.run_id, "sample_id": sid,
                "sample_type": w.sample_type, "rank": r["rank"],
                "taxon": r["taxon"], "label": r["label"],
                "reads": r["reads"], "rel_abundance": f"{rel:.6f}",
                "call": call, "limit_reads": limit, "limit_basis": basis,
                "identity": f"{r['identity']:.4f}", "clusters": r["clusters"],
                "lineage": r["lineage"], "group": r["group"],
                "group_reads": in_group[r["group"]],
                # The share within its group, as VESPA's helminth-only and
                # protozoa-only panels show it.
                "group_share": f"{r['reads'] / in_group[r['group']]:.6f}"
                if in_group[r["group"]] else "",
                "primer_fit": r["primer_fit"], "database_version": version})
    with open(os.path.join(args.outdir, "db_calls.tsv"), "w",
              encoding="utf-8") as fh:
        wr = writer(fh, call_cols)
        wr.writeheader()
        wr.writerows(calls)

    # Per well: the read funnel and the off-target fractions.
    well_cols = ["run_id", "sample_id", "sample_type", "demuxed", "reads_in",
                 "no_primer", "too_short", "unclustered", "clustered",
                 "clusters", "haplotype_splits", "detected_taxa",
                 "unclassified_fraction"] + \
        [f"{g}_reads" for g in GROUPS] + \
        [f"{g}_fraction" for g in GROUPS] + \
        [f"{g}_richness" for g in TARGET_GROUPS] + \
        [f"rank_{r}" for r in ("species", "genus", "family", "above")] + \
        ["database_version"]
    wells_out = []
    for sid in sorted(wells):
        w = wells[sid]
        rd = w.summary.get("reads", {})
        off = Counter()
        unclass = 0
        ranks = Counter()
        for r in rows[sid]:
            off[r["group"]] += r["reads"]
            if r["label"].startswith("unclassified"):
                unclass += r["reads"]
            ranks[r["rank"] if r["rank"] in ("species", "genus", "family")
                  else "above"] += r["reads"]
        d = w.depth or 1
        wells_out.append({
            "run_id": args.run_id, "sample_id": sid,
            "sample_type": w.sample_type,
            "demuxed": per_sample.get(sid, ""),
            **{k: rd.get(k, 0) for k in ("reads_in", "no_primer", "too_short",
                                         "unclustered", "clustered")},
            "clusters": w.summary.get("clusters", 0),
            "haplotype_splits": w.summary.get("haplotype_splits", 0),
            "detected_taxa": sum(1 for c in calls if c["sample_id"] == sid
                                 and c["call"] == "detected"),
            "unclassified_fraction": f"{unclass / d:.4f}",
            **{f"{g}_reads": off[g] for g in GROUPS},
            **{f"{g}_fraction": f"{off[g] / d:.4f}" for g in GROUPS},
            # Detected taxa per target group: VESPA's richness.
            **{f"{g}_richness": sum(1 for c in calls if c["sample_id"] == sid
                                    and c["group"] == g
                                    and c["call"] == "detected")
               for g in TARGET_GROUPS},
            **{f"rank_{k}": ranks[k] for k in ("species", "genus", "family",
                                               "above")},
            "database_version": version})
    with open(os.path.join(args.outdir, "db_wells.tsv"), "w",
              encoding="utf-8") as fh:
        wr = writer(fh, well_cols)
        wr.writeheader()
        wr.writerows(wells_out)

    # Prevalence per taxon: samples it was detected in, of the samples
    # that could detect anything (not blanks, not low depth).
    tested = [w.sid for w in wells.values()
              if w.sample_type not in NEGATIVE_TYPES
              and w.depth >= args.min_depth]
    prev: dict = {}
    for c in calls:
        if c["sample_id"] not in tested:
            continue
        k = (c["group"], c["rank"], c["taxon"], c["label"])
        p = prev.setdefault(k, {"detected": 0, "max_share": 0.0})
        if c["call"] == "detected":
            p["detected"] += 1
            p["max_share"] = max(p["max_share"], float(c["group_share"] or 0))
    with open(os.path.join(args.outdir, "db_prevalence.tsv"), "w",
              encoding="utf-8") as fh:
        wr = writer(fh, ["group", "rank", "taxon", "label", "samples_detected",
                         "samples_tested", "prevalence",
                         "max_share_in_group", "database_version"])
        wr.writeheader()
        for (g, rank, taxon, label), p in sorted(
                prev.items(), key=lambda kv: (GROUPS.index(kv[0][0]),
                                              -kv[1]["detected"], kv[0][2])):
            if not p["detected"]:
                continue
            wr.writerow({"group": g, "rank": rank, "taxon": taxon,
                         "label": label, "samples_detected": p["detected"],
                         "samples_tested": len(tested),
                         "prevalence": f"{p['detected'] / len(tested):.4f}"
                         if tested else "",
                         "max_share_in_group": f"{p['max_share']:.4f}",
                         "database_version": version})

    split = {}
    if args.split_summary and os.path.exists(args.split_summary):
        with open(args.split_summary, encoding="utf-8") as fh:
            split = json.load(fh)
    blanks = [w for w in wells.values() if w.sample_type in NEGATIVE_TYPES]
    qc = {
        "run_id": args.run_id,
        "mode": "database",
        "unit": "relative 18S abundance - never organisms or copies; 18S "
                "copy number varies widely between taxa",
        "database": {"version": version, "name": man.get("name", ""),
                     "primers": man.get("primers", {}),
                     "built": man.get("built", "")},
        "minimap2": (pf_summary.get("reference_database") or {}).get(
            "minimap2", ""),
        "assay": pf_summary.get("assay", {}),
        "settings": {"min_depth": args.min_depth,
                     "absolute_floor": args.absolute_floor,
                     "min_rel_abundance": args.min_rel_abundance,
                     "host_clade": args.host_clade,
                     "floors": next((w.summary.get("settings", {})
                                     .get("floors") for w in wells.values()
                                     if w.summary), None)},
        "call_meaning": CALL_MEANING,
        "groups": {"order": list(GROUPS), "target": list(TARGET_GROUPS),
                   "helminth_classes": sorted(HELMINTH_CLASSES)},
        "samples_tested": len(tested),
        "calls": dict(Counter(c["call"] for c in calls)),
        "negatives": {
            "n": len(blanks),
            "note": ("no negative controls: every detection limit is the "
                     "fixed floor, an assumption" if not blanks else ""),
            "taxa_in_blanks": sorted(
                {f"{c['rank']} {c['taxon']}" for c in calls
                 if c["sample_type"] in NEGATIVE_TYPES and c["reads"] > 0}),
        },
        "empty_wells": empty,
        "low_depth_wells": sorted(w.sid for w in wells.values()
                                  if w.depth < args.min_depth),
        # The design note asks for the rate of internal anchor hits per
        # run: an organism missing from the database cannot be protected
        # by its look-alike mask, so its reads would be cut or
        # quarantined here and nowhere else would say so.
        "internal_anchors": {
            "reads_in": split.get("reads_in"),
            "reads_cut_by_splitter": split.get("reads_cut"),
            "quarantined_internal_adapter":
                (demux.get("by_tier") or {}).get("E"),
            "demuxed_reads": demux.get("total_reads"),
            "mask_rows": (pf_summary.get("reference_database") or {}).get(
                "lookalikes"),
        },
        "near_3prime_calls": sorted(
            f"{c['sample_id']}: {c['rank']} {c['taxon']}" for c in calls
            if "near_3prime" in c["primer_fit"] and c["call"] == "detected"),
    }
    with open(os.path.join(args.outdir, "db_run_qc.json"), "w",
              encoding="utf-8") as fh:
        json.dump(qc, fh, indent=2)
    print(f"database synthesis: {len(wells)} wells, "
          f"{qc['calls'].get('detected', 0)} detected calls, "
          f"{len(blanks)} blanks, release {version}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
