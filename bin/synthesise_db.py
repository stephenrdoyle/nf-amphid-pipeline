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
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    RANKS, load_rank_names, load_release_manifest, rank_meaning, read_table,
    release_version, restore_declared_ids, writer,
)

NEGATIVE_TYPES = {"neg_extraction", "neg_ntc"}
# The identity in a label ("97% to", "R1 74% to"), left out when pooling.
LABEL_PCT = re.compile(r"\d+% to ")

CALLS = ("detected", "below_blanks", "below_floor", "low_depth", "chimera")
CALL_MEANING = {
    "detected": "above the run's blank-derived limit and the relative "
                "abundance floor",
    "below_blanks": "no more reads than the blanks show for this taxon",
    "below_floor": "above the blanks, but under the relative abundance floor",
    "low_depth": "the well has fewer clustered reads than --min-depth, so "
                 "nothing in it is called",
    "chimera": "a PCR chimera of two more abundant sequences in the well; "
               "its reads are counted, never named",
}

GROUPS = ("helminths", "protists", "microsporidia", "host", "fungi",
          "plants", "invertebrates", "bacteria_archaea", "organelle",
          "unassigned", "chimera")
TARGET_GROUPS = ("helminths", "protists", "microsporidia")
OFFTARGET = tuple(g for g in GROUPS if g not in TARGET_GROUPS
                  and g not in ("unassigned", "chimera"))
HELMINTH_CLASSES = {"Nematoda", "Platyhelminthes", "Acanthocephala"}
# Groups the report shows taxon by taxon: the targets, and fungi, which
# are off-target but can dominate a sample (half of every primate sample
# in VESPA's public data) and say something about diet and collection.
DETAIL_GROUPS = TARGET_GROUPS + ("fungi",)


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
        if c.get("origin") == "chimera":
            key = ("", "", "chimera")      # one row per well, all pooled
        else:
            # Pooled on the label without its identity ("Hexamitidae (97%
            # to Enteromonas)" and "(96% to Enteromonas)" are one call);
            # the row keeps the label of its largest cluster.
            key = (c["rank"], c["taxon"], LABEL_PCT.sub("", c["label"]))
        a = agg.setdefault(key, {"reads": 0, "clusters": 0, "identity": 0.0,
                                 "lineage": c["lineage"], "best_hits": [],
                                 "margins": [], "label": c["label"],
                                 "top": -1})
        if key[2] != "chimera" and int(c["reads"]) > a["top"]:
            a["label"], a["top"] = c["label"], int(c["reads"])
        a["reads"] += int(c["reads"])
        a["clusters"] += 1
        a["identity"] = max(a["identity"], float(c["identity"] or 0))
        a["best_hits"].append(c["best_hit"])
        if c.get("db_margin"):
            a["margins"].append(c["db_margin"])
    out = []
    for (rank, taxon, _key), a in agg.items():
        label = "chimera" if _key == "chimera" else a["label"]
        fits = sorted({primer_fit.get(h, "") for h in a["best_hits"]} - {""})
        out.append({"rank": rank, "taxon": taxon, "label": label,
                    "reads": a["reads"], "clusters": a["clusters"],
                    "identity": a["identity"], "lineage": a["lineage"],
                    "group": "chimera" if label == "chimera" else organism_group(
                        [x for x in a["lineage"].split(";")] if a["lineage"]
                        else [], host_clade),
                    "primer_fit": ",".join(fits),
                    "margin": weakest_margin(a["margins"])})
    return out


def weakest_margin(margins: list[str]) -> str:
    """A species call pooled from several clusters is as fragile as its
    most fragile best hit: the fewest edits, or "large" when none has
    another species near."""
    nums = [int(m) for m in margins if m.isdigit()]
    return str(min(nums)) if nums else ("large" if margins else "")


def depth_bucket(rank: str, meaning: str) -> str:
    """Species, genus, family or "above family", for the naming-depth
    summaries, read through rank_names.tsv: PR2's family column holds an
    order in nematodes and a class in most fungi, so a call there is
    above family unless the table says the column is a family (or has no
    row for the clade, when it means roughly what it says)."""
    if rank in ("species", "genus"):
        return rank
    if rank == "family":
        real = [x for x in meaning.split("|") if x and x != "placeholder"]
        return "family" if not real or real[0] == "family" else "above family"
    return "above family"


def fragile_species(calls: list[dict]) -> list[dict]:
    """Detected species calls whose name rests on 1-2 edits, by species:
    the margin, the samples, and the group, most samples first."""
    by: dict = {}
    for c in calls:
        if c["call"] == "detected" and c["margin"] in ("1", "2"):
            f = by.setdefault(c["taxon"], {"taxon": c["taxon"],
                                           "label": c["label"],
                                           "group": c["group"],
                                           "margin": int(c["margin"]),
                                           "samples": []})
            f["margin"] = min(f["margin"], int(c["margin"]))
            f["samples"].append(c["sample_id"])
    return sorted(by.values(), key=lambda f: (-len(f["samples"]), f["taxon"]))


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


def read_funnel(wells_out: list[dict], demux: dict, split: dict,
                calls: list[dict] | None = None,
                family_floor: float | None = None) -> dict:
    """Every read from the sequencer to a classified target, and where
    the rest went.

    The per-well table starts at reads that reached a sample, which is
    already net of splitting and demultiplexing. A well that lost most
    of its reads to unreadable barcodes looks merely shallow there, and
    so does a well whose reads were mostly primer dimer or bacteria. The
    funnel is the one place those differ, and they call for different
    fixes: re-index or re-load, fix the reaction, or accept the sample.

    Stages, each with what was lost getting to it:
      reads off the sequencer       nanopore: the split's input
      amplicon units after splitting  a read can hold several molecules
      demultiplexed to a sample     by demux tier (B one end, C chimeric...)
      read pairs in                 Illumina: there is no demux here
      both primers found            no primer pair, too short, low quality
      in a consensus sequence       too few reads to form one
      classified                    PCR chimeras removed
        of which unclassified       under the family floor to every record
      target organisms              host, bacteria, plants... by group
        helminths, protists, microsporidia, each by the rank reached
    Rows marked "sub" are a breakdown of the stage above, not a further
    loss: their detail is what they hold rather than what was lost.

    It also checks the seam between demultiplexing and assignment: the
    reads the demultiplexer wrote must be the reads assignment read.
    Nothing else compares the two, so a sample file dropped between the
    steps would otherwise go unnoticed.
    """
    total = Counter()
    for w in wells_out:
        # chimera_reads is one of the group columns, so not listed again.
        for k in ("reads_in", "no_primer", "too_short", "low_quality",
                  "unclustered", "unmapped", "clustered",
                  *(f"{g}_reads" for g in GROUPS)):
            total[k] += int(w.get(k) or 0)
    stages = []
    if split.get("reads_in"):
        stages.append({"stage": "reads off the sequencer",
                       "reads": split["reads_in"], "lost_here": {}})
        dropped = split.get("pieces_dropped_below_min_length", 0)
        stages.append({
            "stage": "amplicon units after splitting at internal adapters",
            "reads": split.get("pieces_out", demux.get("total_reads", 0)),
            "lost_here": ({"pieces too short to hold an amplicon": dropped}
                          if dropped else {})})
    elif demux.get("total_reads"):
        stages.append({"stage": "reads in", "reads": demux["total_reads"],
                       "lost_here": {}})
    if demux.get("total_reads"):
        tiers = demux.get("by_tier") or {}
        meaning = demux.get("tier_meaning") or {}
        stages.append({
            "stage": "demultiplexed to a sample",
            "reads": demux.get("assigned_reads", tiers.get("A", 0)),
            # Name the tiers rather than lump them: B is a read with one
            # unreadable end, C a chimera of two wells, D no match - and
            # they point at different problems.
            "lost_here": {(f"{meaning[k].split(' - ')[0]} (tier {k})"
                           if k in meaning else f"tier {k}"): v
                          for k, v in sorted(tiers.items())
                          if k != "A" and v}})
    else:
        stages.append({"stage": "read pairs in", "reads": total["reads_in"],
                       "lost_here": {}})
    primed = (total["reads_in"] - total["no_primer"] - total["too_short"]
              - total["low_quality"])
    stages.append({
        "stage": "both primers found", "reads": primed,
        "lost_here": {k: v for k, v in (
            ("no primer pair (dimer, off-target or truncated)",
             total["no_primer"]),
            ("too short", total["too_short"]),
            ("low quality", total["low_quality"])) if v}})
    # Nanopore counts its unmapped reads inside unclustered; Illumina has
    # no unmapped count. Shown apart where known: a read with no hit in
    # the release is a different loss from a read in too small a group.
    unmapped = min(total["unmapped"], total["unclustered"])
    stages.append({
        "stage": "in a consensus sequence", "reads": total["clustered"],
        "lost_here": {k: v for k, v in (
            ("too few reads to form a consensus",
             total["unclustered"] - unmapped),
            ("no hit in the release", unmapped)) if v}})
    stages.append({
        "stage": "classified (PCR chimeras removed)",
        "reads": total["clustered"] - total["chimera_reads"],
        "lost_here": ({"PCR chimeras": total["chimera_reads"]}
                      if total["chimera_reads"] else {})})
    # From the calls: how much of the classified set is "unclassified,
    # NN% to X" (a sequence no rank's floor would name), by group, and
    # how deep each target group was named.
    unclass, depth = Counter(), defaultdict(Counter)
    for c in calls or []:
        if c["group"] == "chimera":
            continue
        n = int(c["reads"])
        if c["label"].startswith("unclassified"):
            unclass[c["group"]] += n
            depth[c["group"]]["unclassified"] += n
        else:
            depth[c["group"]][c.get("depth") or depth_bucket(
                c["rank"], c.get("rank_meaning", ""))] += n
    if calls is not None:
        stages.append({
            # Named by the floor it missed, from the run's own settings:
            # under the family floor a sequence is named to no family,
            # and is "unclassified, NN% to X". A tie stopped above
            # family is not counted here - the database holds it.
            "stage": "of which unclassified ("
                     + (f"under {100 * family_floor:.0f}% to any record"
                        if family_floor else "under the family floor")
                     + ", named to no family)",
            "reads": sum(unclass.values()), "sub": True,
            "lost_here": {g.replace("_", " / "): unclass[g]
                          for g in GROUPS if unclass[g]}})
    stages.append({
        "stage": "target organisms (" + ", ".join(TARGET_GROUPS) + ")",
        "reads": sum(total[f"{g}_reads"] for g in TARGET_GROUPS),
        "lost_here": {g.replace("_", " / "): total[f"{g}_reads"]
                      for g in GROUPS if g not in TARGET_GROUPS
                      and g != "chimera" and total[f"{g}_reads"]}})
    if calls is not None:
        for g in TARGET_GROUPS:
            stages.append({
                "stage": g, "reads": total[f"{g}_reads"], "sub": True,
                "lost_here": {f"{k}": depth[g][k] for k in
                              ("species", "genus", "family", "above family",
                               "unclassified") if depth[g][k]}})
    start = stages[0]["reads"] or 0
    for st in stages:
        st["percent_of_start"] = (round(100 * st["reads"] / start, 2)
                                  if start else None)
    out = {"stages": stages}
    if demux.get("assigned_reads") is not None:
        out["handoff_balances"] = demux["assigned_reads"] == total["reads_in"]
        out["handoff"] = {"written_by_demultiplexing": demux["assigned_reads"],
                          "read_by_assignment": total["reads_in"]}
    # Why reads reached only one barcode. On nanopore this is usually
    # the largest loss of the run, and the tier count alone does not say
    # whether it is the molecule (a read that stops before its second
    # adapter), the basecalls (an adapter or index too damaged to read)
    # or the plate design (one end cannot name a well when indexes are
    # combined). The demultiplexer records a reason for every end of
    # every read; these are whole-run counts, not a sample.
    if demux.get("i5_end_reasons") or demux.get("i7_end_reasons"):
        grid = demux.get("combination_grid") or {}
        out["barcode_ends"] = {
            "reads": demux.get("total_reads", 0),
            "both_ends": demux.get("assigned_reads"),
            "one_end": (demux.get("by_tier") or {}).get("B", 0),
            "i5": demux.get("i5_end_reasons") or {},
            "i7": demux.get("i7_end_reasons") or {},
            "i5_used": grid.get("i5_used"), "i7_used": grid.get("i7_used"),
            "wells_declared": grid.get("wells_declared"),
            "rescued_by_combination":
                (demux.get("combination_check") or {}).get("resolved", 0),
            "accept_single_end":
                (demux.get("settings") or {}).get("accept_single_end"),
            "basecall_model":
                (split.get("basecalling") or {}).get("basecall_model", ""),
        }
    # A well that lost most of its reads before classification looks
    # shallow in the per-well table; say which, and where they went.
    out["wells_losing_over_half"] = [
        {"sample_id": w["sample_id"], "reads_in": int(w["reads_in"] or 0),
         "classified_fraction": round(int(w["clustered"] or 0)
                                      / int(w["reads_in"]), 4),
         "fates": {k: int(w.get(k) or 0) for k in
                   ("no_primer", "too_short", "low_quality", "unclustered")}}
        for w in wells_out
        if int(w["reads_in"] or 0) and int(w["clustered"] or 0)
        < 0.5 * int(w["reads_in"])]
    return out


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
    # What the rank reached actually is, where PR2's column is not that
    # rank (rank_names.tsv, v0.8.0 on): "order" for a nematode's family
    # column. Blank where the table has no row or no table exists.
    rank_names = load_rank_names(args.database)
    for rs in rows.values():
        for r in rs:
            lin = r["lineage"].split(";") if r["lineage"] else []
            m = (rank_meaning(rank_names, lin, r["rank"])
                 if r["rank"] in ("subdivision", "class", "order", "family")
                 else "")
            r["rank_meaning"] = "" if m == r["rank"] else m
            r["depth"] = depth_bucket(r["rank"], m)
    limits, n_blanks = blank_limits(wells, rows, args.min_depth)

    os.makedirs(args.outdir, exist_ok=True)
    call_cols = ["run_id", "sample_id", "sample_type", "rank", "taxon",
                 "label", "reads", "rel_abundance", "call", "limit_reads",
                 "limit_basis", "identity", "clusters", "lineage",
                 "group", "group_reads", "group_share", "primer_fit",
                 "margin", "rank_meaning", "database_version"]
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
            if r["group"] == "chimera":
                call = "chimera"
            elif w.depth < args.min_depth:
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
                "primer_fit": r["primer_fit"],
                # Edits from the best record to the nearest other species
                # (species calls, v0.8.0 releases on). 1-2 is fragile.
                "margin": r["margin"],
                # What the rank column really is for this clade, where it
                # is not what it says (rank_names.tsv).
                "rank_meaning": r["rank_meaning"],
                "database_version": version})
    with open(os.path.join(args.outdir, "db_calls.tsv"), "w",
              encoding="utf-8") as fh:
        wr = writer(fh, call_cols)
        wr.writeheader()
        wr.writerows(calls)

    # Per well: the read funnel and the off-target fractions.
    well_cols = ["run_id", "sample_id", "sample_type", "platform",
                 "demuxed", "reads_in", "no_primer", "too_short",
                 "low_quality", "unclustered", "unmapped", "clustered",
                 "clusters", "haplotype_splits", "detected_taxa",
                 "unclassified_fraction"] + \
        [f"{g}_reads" for g in GROUPS] + \
        [f"{g}_fraction" for g in GROUPS] + \
        [f"{g}_richness" for g in DETAIL_GROUPS] + \
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
            ranks[{"above family": "above"}.get(r["depth"], r["depth"])] += \
                r["reads"]
        d = w.depth or 1
        wells_out.append({
            "run_id": args.run_id, "sample_id": sid,
            "sample_type": w.sample_type,
            "platform": w.summary.get("platform", "nanopore"
                                      if w.summary else ""),
            "demuxed": per_sample.get(sid, ""),
            **{k: rd.get(k, 0) for k in ("reads_in", "no_primer", "too_short",
                                         "low_quality", "unclustered",
                                         "unmapped", "clustered")},
            "clusters": w.summary.get("clusters", 0),
            "haplotype_splits": w.summary.get("haplotype_splits", 0),
            "detected_taxa": sum(1 for c in calls if c["sample_id"] == sid
                                 and c["call"] == "detected"),
            "unclassified_fraction": f"{unclass / d:.4f}",
            **{f"{g}_reads": off[g] for g in GROUPS},
            **{f"{g}_fraction": f"{off[g] / d:.4f}" for g in GROUPS},
            # Detected taxa per target group (VESPA's richness), and fungi.
            **{f"{g}_richness": sum(1 for c in calls if c["sample_id"] == sid
                                    and c["group"] == g
                                    and c["call"] == "detected")
               for g in DETAIL_GROUPS},
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
        # PCR chimeras found per well (assign_db.chimera_scan): counted,
        # never named, kept out of every composition.
        "chimeras": {"reads": sum(int(w2["chimera_reads"]) for w2 in wells_out),
                     "of_clustered": round(sum(int(w2["chimera_reads"])
                                               for w2 in wells_out)
                                           / max(1, sum(int(w2["clustered"] or 0)
                                                        for w2 in wells_out)), 4),
                     "wells_over_5pct": sorted(
                         w2["sample_id"] for w2 in wells_out
                         if float(w2["chimera_fraction"]) > 0.05)},
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
        "read_funnel": read_funnel(
            wells_out, demux, split, calls,
            (next((w.summary.get("settings", {}).get("floors") for w in
                   wells.values() if w.summary), None) or {}).get("family")),
        # Detected species calls whose name rests on 1-2 edits: the
        # nearest other species in the release is that close.
        "fragile_species_calls": fragile_species(calls),
        # Illumina pairs whose two windows point at different records
        # (assign_db_pairs): usually chimeras or artefacts, not gaps.
        "windows_disagree": {
            "reads": sum(int(c["reads"]) for c in calls
                         if c["label"].startswith("windows disagree")),
            "samples": len({c["sample_id"] for c in calls
                            if c["label"].startswith("windows disagree")})},
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
