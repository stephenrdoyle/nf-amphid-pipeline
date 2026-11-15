"""
vespa_v4: vespa-18s-refdb releases, 18S V4 - identity floors over
minimap2's hits.

A V4 record is a whole primer-to-primer insert, so a consensus (or a
pair's two windows) is compared with it end to end and an unaligned end
is a real difference. A sequence's call is the LCA of the records within
--lca-margin of its best alignment score, cut back to the deepest rank
whose identity floor the best record clears (DEFAULT_FLOORS for single
reads, the stricter ILLUMINA_FLOORS for pairs, which read only two
windows of the insert). Where minimap2 clipped a divergent end, the
identity is taken from a full alignment instead (refine_identity).

Pairs whose windows -x sr cannot align end to end are mapped again with
nanopore's penalties (SR_DIVERGENT), and a pair whose windows point at
different records is called at what the two agree on.

Read types: single_long_noisy (nanopore consensus), paired_short
(Illumina windows) and single_long_accurate (PacBio HiFi: the consensus
route with its PacBio steps, and the same floors). The floors were
calibrated on nanopore consensus sequences; a leave-out benchmark on
simulated HiFi reads (tests/bench_db_floors.py --error-model
pacbio_hifi, docs/CONTRACTS.md) measured them before HiFi reads were
taken, and again with the PacBio steps: the same calls.

Standard-library Python plus minimap2.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ampseqlib import (  # noqa: E402
    RANKS, apply_floors, bounded_edit, consensus_hits, hit_identity,
    lca_depth, map_windows, pair_hits, ranks_of, read_fasta, revcomp,
    write_fasta,
)
from protocols import Protocol  # noqa: E402

# Floors on the best hit's identity, per rank. A consensus is named at a
# rank only if its best hit is at least this close; otherwise it stops at
# the next rank up. From the leave-out benchmark at consensus error rates:
# species 0.995 kept 90% of correct species calls and cut wrong calls on
# missing species to 5-8%. Genus was 0.98 in the prototype, which let
# 27% of consensus from missing genera take a wrong genus on a 150-sample
# leave-out run of v0.4.0; 0.99 cut that to 20% at no cost to present
# taxa, the rest moving honestly to family (docs/REFERENCE_DATABASE.md
# section 11.3, tests/bench_db_floors.py).
DEFAULT_FLOORS = {"species": 0.995, "genus": 0.99, "family": 0.95}

# Floors for pairs, stricter than nanopore's. Over the two 2x150 windows
# (~260bp) the bases that separate sister taxa are mostly in the gap the
# pair never reads, so the nanopore floors - fractions calibrated on
# whole inserts - let a missing species take its sister's name. On a
# 150-sample leave-out run of v0.4.0 at 2x150 (tests/bench_db_floors.py
# --illumina 150), species 0.997 (an exact match over 260bp) cut wrong
# names on missing species from 16% to 9% at a cost of 2.6 points of
# species calls on present ones; genus 0.995 cut wrong genera on missing
# genera from 18% to 17%. docs/REFERENCE_DATABASE.md section 11.5.
ILLUMINA_FLOORS = {"species": 0.997, "genus": 0.995, "family": 0.95}


# Window mapping. SR_DIVERGENT keeps sr's seeds and chaining but scores
# as map-ont does (mismatch 4, gap 4+2k), so an alignment carries on
# through a divergent stretch instead of stopping at it (classify_pairs).
SR_DIVERGENT = ("-x", "sr", "-B4", "-O4,24", "-E2,1")
SR_CLIP = 0.10


def classify_hits(hits: list[dict], tax: dict, floors: dict,
                  margin: float) -> dict:
    """One consensus's call from its database hits."""
    ranks = ranks_of(tax)
    if not hits:
        return {"rank": "", "depth": -1, "identity": 0.0, "coverage": 0.0,
                "best_hit": "", "near_hits": 0, "lca_depth": -1}
    top = max(h["AS"] for h in hits)
    near = [h for h in hits if h["AS"] >= top * (1 - margin)]
    best = max(near, key=lambda h: (h["AS"], hit_identity(h), h["t"]))
    ident = hit_identity(best)
    names = sorted({h["t"] for h in near})
    # A record with no taxonomy at all - v0.4.0 has nine, conflicts its
    # build left unresolved - can vouch for nothing, and as a member of
    # the LCA it vetoed everything: human 18S tied with one stopped at
    # no rank. It is left out of the LCA, unless it is all there is.
    known = [t for t in names if tax.get(t, {}).get(ranks[0])]
    lca = lca_depth(known or names, tax)
    depth = apply_floors(lca, ident, floors, ranks)
    return {
        "depth": depth, "rank": ranks[depth] if depth >= 0 else "",
        "identity": ident,
        "coverage": (best["qe"] - best["qs"]) / max(1, best["qlen"]),
        "best_hit": best["t"], "near_hits": len(names),
        "lca_depth": lca,
    }


def aligned_identity(parts: list[tuple[str, str]], record: str,
                     max_frac: float = 0.35) -> float | None:
    """Identity of a consensus against one record, from a full alignment
    of each piece rather than minimap2's local block.

    hit_identity() charges every base minimap2 leaves unaligned at an end as
    a difference. For a sequence close to its record that is nothing -
    the alignment runs end to end - but -x sr stops extending at the
    first divergent stretch, so a sequence 6% from its nearest record
    was charged for a whole clipped end it would have matched most of:
    VESPA's Entamoeba-like candidates read 90% where BLAST gives 95% and
    an end-to-end alignment 94%, enough to fall under the family floor
    and read as unclassified. Here each piece is aligned in full
    (bounded_edit, unit costs):

      "whole"  against the whole record (a primer-to-primer insert);
      "start"  against the record's start, its tail free (R1's window);
      "end"    against the record's end, its head free (R2's window,
               given in the insert's orientation).

    Edits over the pieces' length; None when a piece is further than
    max_frac from the record, where no floor could be cleared anyway.
    """
    edits = length = 0
    for seq, where in parts:
        if not seq:
            continue
        cap = max(1, int(max_frac * len(seq)))
        if where == "whole":
            d = bounded_edit(seq, record, cap)
        elif where == "start":
            d = bounded_edit(seq, record, cap, free_b_suffix=True)
        else:
            d = bounded_edit(seq[::-1], record[::-1], cap, free_b_suffix=True)
        if d > cap:
            return None
        edits += d
        length += len(seq)
    return 1 - edits / length if length else None


def refine_identity(call: dict, parts: list[tuple[str, str]],
                    record: str | None, floors: dict,
                    ranks: list[str] = RANKS) -> dict:
    """Raise a call's identity to its full-alignment value, and re-apply
    the floors from its LCA. Never lowers it: where minimap2 aligned end
    to end the two measures agree, and the floors were calibrated on
    those calls, so they stay exactly as they were."""
    if not record or not call.get("best_hit"):
        return call
    new = aligned_identity(parts, record)
    if new is None or new <= call["identity"] + 1e-9:
        return call
    depth = (apply_floors(call["lca_depth"], new, floors, ranks)
             if call.get("lca_depth", -1) >= 0 else call["depth"])
    return {**call, "identity": new, "depth": depth,
            "rank": ranks[depth] if depth >= 0 else "",
            "identity_minimap2": call["identity"]}

class VespaV4(Protocol):
    name = "vespa_v4"
    title = "vespa-18s-refdb, 18S V4: identity floors"
    read_types = {
        "single_long_noisy": "reference-guided consensus per group, "
                             "classified by identity floors",
        "paired_short": "two fixed windows per pair, denoised, classified "
                        "by the stricter pair floors",
        "single_long_accurate": "the consensus route, as nanopore's, "
                                "classified by the same floors",
    }
    not_supported: dict = {}
    calling_rules = ("floors",)
    required_files = ("refdb.fasta", "taxonomy.tsv",
                      "resolution_nanopore.tsv")

    def classify_consensus(self, ctx, seqs: dict) -> dict:
        hits = consensus_hits(ctx.mm, ctx.index, ctx.tmp, seqs)
        calls = {}
        for k in seqs:
            calls[k] = classify_hits(hits.get(f"c{k}", []), ctx.tax,
                                     ctx.floors, ctx.args.lca_margin)
        # Identity from a full alignment to the best record, where
        # that is higher than minimap2's (aligned_identity).
        need = {c["best_hit"] for c in calls.values() if c["best_hit"]}
        recs = {i: s for i, s in read_fasta(ctx.db_fasta) if i in need}
        for k in seqs:
            calls[k] = refine_identity(
                calls[k], [(seqs[k], "whole")],
                recs.get(calls[k]["best_hit"]), ctx.floors, ctx.ranks)
        return calls

    def classify_pairs(self, ctx, centroids: list) -> tuple[dict, int]:
        mm, index, tmp, args = ctx.mm, ctx.index, ctx.tmp, ctx.args
        tax, floors, ranks = ctx.tax, ctx.floors, ctx.ranks
        calls: dict[int, dict] = {}
        hits, paths, extra = pair_hits(mm, index, tmp, centroids,
                                       args.max_secondary)
        # -x sr's penalties (mismatch 8, gap open 12) are made for reads
        # of a genome they belong to. Against a taxon ~15% from its
        # nearest record - VESPA's NHP gregarines, Hexamita,
        # Hypotrichomonas - the extension stops at the first divergent
        # stretch, the window aligns as two pieces, and the pair was
        # joined on one 60bp piece: identity ~0.5 over the pair, best hit
        # a record that merely shares the conserved start. A pair sr
        # leaves more than SR_CLIP of a window unaligned is mapped again,
        # BOTH its windows, with nanopore's penalties on the same seeds,
        # and that mapping replaces the first; the rest keep sr's, so
        # ordinary calls and the floors calibrated on them are unchanged.
        # Judged on the window's best-scoring hit, the one its call rests
        # on: sr will often also give some record a gappy end-to-end
        # alignment scoring below the clean partial ones. A window with
        # no hit at all counts as unaligned too.
        def poor(name, q):
            hs = hits[name].get(q)
            if not hs:
                return True
            h = max(hs.values(), key=lambda h: (h["AS"], h["t"]))
            return h["qs"] + h["qlen"] - h["qe"] > SR_CLIP * h["qlen"]
        redo = {
            f"p{k}" for k, (x, y, _n) in enumerate(centroids)
            if (poor("r1", f"p{k}") or poor("r2", f"p{k}") if y
                else poor("full", f"p{k}"))}
        if redo:
            seqs = {"r1": lambda x, y: x if y else "",
                    "r2": lambda x, y: y, "full": lambda x, y: "" if y else x}
            again_paths = {}
            for name, pick in seqs.items():
                again_paths[name] = paths[name] + ".divergent"
                write_fasta(again_paths[name], (
                    (f"p{k}", pick(x, y)) for k, (x, y, _n) in
                    enumerate(centroids) if f"p{k}" in redo and pick(x, y)))
            again = map_windows(mm, index, again_paths, extra, SR_DIVERGENT)
            for name in paths:
                for q in redo:
                    hits[name].pop(q, None)
                    if q in again[name]:
                        hits[name][q] = again[name][q]
        remapped = sum(centroids[int(q[1:])][2] for q in redo)
        for k, (x, y, _n) in enumerate(centroids):
            q = f"p{k}"
            if not y:
                calls[k] = classify_hits(
                    [h for h in hits["full"].get(q, {}).values()
                     if h["strand"] == "+"], tax, floors, args.lca_margin)
                continue
            joined = []
            a_hits, b_hits = hits["r1"].get(q, {}), hits["r2"].get(q, {})
            for t, h1 in a_hits.items():
                h2 = b_hits.get(t)
                # The shape a pair has on its insert: R1 forward from the
                # start, R2 reverse back from the end.
                if h2 is None or h1["strand"] != "+" or h2["strand"] != "-":
                    continue
                if h1["ts"] > 30 or h2["tlen"] - h2["te"] > 30:
                    continue
                c1 = h1["qs"] + h1["qlen"] - h1["qe"]
                c2 = h2["qs"] + h2["qlen"] - h2["qe"]
                joined.append({
                    "t": t, "AS": h1["AS"] + h2["AS"],
                    "match": h1["match"] + h2["match"],
                    "alen": h1["alen"] + h2["alen"],
                    "qlen": h1["qlen"] + h2["qlen"], "qs": c1 + c2,
                    "qe": h1["qlen"] + h2["qlen"]})
            if joined:
                calls[k] = classify_hits(joined, tax, floors,
                                         args.lca_margin)
                continue
            # No record carries both windows in a pair's shape: the two
            # halves point at different records (a chimera, or a lineage
            # the database has only parts of). Each window is classified
            # alone and the call is their LCA, at the lower identity - it
            # can say no more than both halves agree on.
            one = [classify_hits(list(hs.values()), tax, floors,
                                 args.lca_margin)
                   for hs in (a_hits, b_hits) if hs]
            if len(one) < 2:
                calls[k] = one[0] if one else classify_hits(
                    [], tax, floors, args.lca_margin)
                calls[k]["split_mates"] = True
                continue
            lin = [tax.get(c["best_hit"], {}) for c in one]
            depth = -1
            for di, r in enumerate(ranks):
                if di <= min(c["depth"] for c in one) and lin[0].get(r) \
                        and lin[0].get(r) == lin[1].get(r):
                    depth = di
                else:
                    break
            worse = min(one, key=lambda c: c["identity"])
            calls[k] = {**worse, "depth": depth,
                        "rank": ranks[depth] if depth >= 0 else "",
                        "lca_depth": depth, "split_mates": True,
                        # Each window's own call, for the label: two
                        # windows that disagree below what they share
                        # are usually a chimera or an artefact.
                        "windows": [(c["best_hit"], c["identity"],
                                     c["depth"]) for c in one]}
        # Identity from a full alignment of the windows to the best
        # record, where that is higher than minimap2's: -x sr clips a
        # divergent end and hit_identity() charges all of it
        # (aligned_identity). Joined pairs and whole inserts only; a
        # pair whose windows disagree has no one record to align to.
        need = {c["best_hit"] for c in calls.values()
                if c["best_hit"] and not c.get("split_mates")}
        need |= {h for c in calls.values() for h, _i, _d in
                 c.get("windows") or [] if h}
        recs = {i: s for i, s in read_fasta(
            os.path.join(args.database, "refdb.fasta")) if i in need}
        for k, (x, y, _n) in enumerate(centroids):
            c = calls.get(k)
            if c and c.get("windows") and len(c["windows"]) == 2:
                # Each window against its own best record, so the label's
                # "R1 NN%" is not minimap2's clipped figure: a gregarine
                # R1 read "37%" where it aligns 74% end to end.
                ws = []
                for (h, i, dd), part in zip(c["windows"], ((x, "start"),
                                                           (revcomp(y), "end"))):
                    v = aligned_identity([part], recs.get(h, ""), 0.5) \
                        if recs.get(h) else None
                    ws.append((h, max(i, v or 0.0), dd))
                calls[k] = {**c, "windows": ws,
                            "identity": min(w[1] for w in ws)}
            if not c or c.get("split_mates"):
                continue
            parts = [(x, "start"), (revcomp(y), "end")] if y else \
                [(x, "whole")]
            calls[k] = refine_identity(c, parts,
                                       recs.get(c["best_hit"]), floors,
                                       ranks)
        return calls, remapped


PROTOCOL = VespaV4()
