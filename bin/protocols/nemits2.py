"""
nemits2: nemabiome-its2-refdb releases, ITS2 - span scoring, the gap rule
and its tentative tier.

An ITS2 record stops where the ~47-base 28S stretch before NC2 starts,
and many records are partial, so a read always runs past its record.
Each candidate record minimap2 proposes is therefore scored exactly over
the record's own span, the query's overhang free (span_scores), and the
call follows the release's calling rule. The gap rule: a species needs
>= 99% identity AND the best record of any other species >= 4 edits
further off over the same span (span_call). The tiered rule, the
default: the gap rule's calls, plus a species named at >= 98% with the
same gap, flagged tentative, with the gap rule's own call beside it.
Every call says how close it came: the tied species and their identity,
the nearest other species and how many edits further (nearest, next).

On Illumina the release's windows count record bases, and R2 opens with
the 28S stretch, so R2's window takes that stretch too (tail: 130 + 83
at 2x150 is 130 + 130 read bases); --primer-offset shortens both
windows by the spacer bases some protocols put before each primer (0-3
for NC1/NC2). Mates that overlap (most ITS2 pairs at 2x250) are merged
into the insert and scored whole (merge_pair): in two overlapping
windows every difference in the overlap is counted twice, and the
AT-repeat stutter of an Oesophagostomum ITS2 (two bases) fell under the
99% a species call needs.

Read types: single_long_noisy (nanopore consensus, carrying the 28S
stretch), paired_short, and single_long_accurate (PacBio HiFi or CCS:
the consensus route with its PacBio steps, checked on S-BSST527 against
DADA2 ASVs of the same reads; docs/CONTRACTS.md).

Standard-library Python plus minimap2.
"""

from __future__ import annotations

import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ampseqlib import (  # noqa: E402
    anchored_distance, apply_floors, consensus_hits, infix_distance,
    lca_depth, pair_hits, parse_floors, ranks_of, read_fasta, revcomp,
)
from protocols import Protocol, scores_over_reference  # noqa: E402

# The gap rule, for releases scored over their records' span (ITS2):
# species needs >= 99% identity AND the best record of any other species
# >= 4 edits further off over the same span; genus >= 97%, family >= 95%.
# From nemabiome-its2-refdb's leave-one-out (3,164 queries, 1,210
# species): it named 52.9% of species present at 0.7% wrong, against
# 18.5% under the V4 floors; reads of missing species took a sister's
# name 6.6% of the time (V4 floors 2.4%, assignTaxonomy 42%). Per-genus
# or per-family limits did no better on held-out species.
GAP_FLOORS = {"species": 0.99, "genus": 0.97, "family": 0.95}
GAP_EDITS = 4
# ... or an exact match to the record (0 edits over its span) with the
# nearest other species >= 2 edits further off: the release's own
# classifier (its2_classify.py GAP_EXACT). In its leave-one-out an exact
# match with another species 1 edit off was right 60% of the time, at 2-3
# edits 91%; below 100% identity 2-3 edits stay too close.
GAP_EXACT = 2
# The tiered rule: the gap rule's species calls, plus tentative ones at
# >= 98% identity with the same gap (the release's its2_classify.py
# TENTATIVE_IDENTITY, and its recommended default). Half the present-
# species reads the gap rule stops at genus have no other species within
# 4 edits but are 1-3% from their own records - intraspecific variation.
# In the release's leave-one-out the 98% tier raised right species calls
# from 56% to 68% where the species is present, and wrong ones from 5.4%
# to 9.0% where it is hidden: a read 1-2% from its own species looks the
# same as one 1-2% from an unsampled sister, which is why the call is
# flagged and the gap rule's call is kept beside it. On real data it
# named the goat line of T. circumcincta and Oe. venulosum / Oe. radiatum
# reads the gap rule left at genus.
TENTATIVE_IDENTITY = 0.98
# Span scoring. A record must be aligned over at least MIN_SPAN bases to
# count as evidence at all (a 77-base fragment 32 edits off otherwise
# outranks full-length records by identity); among equally identical
# records, those within SPAN_SLACK bases of the longest span tie and the
# rest - partials matching over less of the read - drop out. A window's
# overhang past its record is free for at most WINDOW_TRIM bases beyond
# the 28S stretch. SPAN_CANDIDATES: records minimap2 proposes that are
# scored exactly, enough for every tie and the nearest other species.
MIN_SPAN = 100
SPAN_SLACK = 5
WINDOW_TRIM = 30
SPAN_CANDIDATES = 100


# span scoring and the gap rule (releases scored over their records' span)
# --------------------------------------------------------------------------


def hit_span_identity(h: dict) -> float:
    """minimap2's estimate of a hit's identity over the record's span:
    matches over the alignment plus the record's unaligned ends (the
    query's are free)."""
    return h["match"] / max(1, h["alen"] + h["ts"] + h["tlen"] - h["te"])


def span_candidates(hits: list[dict], n: int = SPAN_CANDIDATES) -> list[str]:
    """The records minimap2 proposes that are scored exactly: the n best
    by estimated identity ("ident", else hit_span_identity), then
    alignment score, and every record tied with the n-th. By identity,
    not score: a partial record scores lower on length alone, and a
    species with hundreds of records pushed a partial of the sister
    species 3 edits off - the one that decides the gap - to 247th by
    score (an H. placei sequence, 7% of a sheep set)."""
    best: dict[str, tuple] = {}
    for h in hits:
        k = (round(h["ident"] if "ident" in h else hit_span_identity(h), 6),
             h["AS"])
        best[h["t"]] = max(best.get(h["t"], k), k)
    order = sorted(best, key=lambda t: (-best[t][0], -best[t][1], t))
    if len(order) > n:
        cut = best[order[n - 1]]
        order = [t for t in order if best[t] >= cut]
    return order


def span_scores(parts: list[tuple[str, str]], cands, refs: dict,
                tail: int = 0) -> list[tuple]:
    """(edits, span, record) for each candidate record, scored over the
    record's span (manifest score_over), each piece in full:

      "whole"  the record placed whole inside the query, the query's
               bases either side free (a consensus, or an insert read
               through: ITS2 plus the 28S stretch);
      "start"  a window against the record's start, the record past the
               window free (R1);
      "end"    a window against the record's end, given in the insert's
               orientation; the window's last `tail` bases (R2's 28S
               stretch) run past the record and are free.

    A window may also overhang a partial record's outer end by up to
    WINDOW_TRIM bases for free; beyond that the overhang is charged.
    Edits and spans add up over the pieces. A record no piece fits is
    left out."""
    out = []
    for c in cands:
        r = refs.get(c)
        if not r:
            continue
        d = sp = 0
        for seq, where in parts:
            if not seq:
                continue
            if where == "whole":
                e, k = infix_distance(r, seq), len(r)
            elif where == "start":
                e, k = anchored_distance(
                    r, seq, min(len(r), max(1, len(seq) - WINDOW_TRIM)))
            else:
                e, k = anchored_distance(
                    r[::-1], seq[::-1],
                    min(len(r), max(1, len(seq) - tail - WINDOW_TRIM)))
            d += e
            sp += k
        if sp and d != float("inf"):
            out.append((int(d), sp, c))
    return out


def span_call(scores: list[tuple], tax, rule: str, floors: dict,
              gap_edits: int = GAP_EDITS, qlen: int = 0,
              min_species_span: int = MIN_SPAN,
              gap_exact: int = GAP_EXACT,
              tentative_identity: float = TENTATIVE_IDENTITY) -> dict:
    """One query's call from its span scores (span_scores), in
    classify_hits()'s shape plus what the gap rule adds.

    Records are ranked by identity over their span, not by raw edits: on
    edits a short partial, with fewer bases to get wrong, beat full-length
    records. The ties are every record at the best identity whose span is
    within SPAN_SLACK of the longest such span - a partial matching only
    over a shorter span does not tie with a full-length match (an H.
    placei sequence, 23% of one cattle set, tied with a partial H.
    contortus record and stopped at genus). The call is their LCA, cut
    back by the floors; under the gap rule a species also needs the
    nearest record of any OTHER species to be gap_edits further off over
    the same span - or, for an exact match, gap_exact (0: never). The
    tiered rule makes every call the gap rule makes, and also names a
    species the gap rule stops at genus when the identity is at least
    tentative_identity and the gap is as clear: that call is flagged
    tentative, and the gap rule's own call is kept beside it. A
    shared record (its members name two species) counts
    as another species; an unnamed one (hybrid, cf., environmental) does
    not - counting three hybrid records once cut a cattle set's species
    calls from 87% to 35%. With no other species among the candidates,
    the gap is the whole span. A tied record the release flags as
    contained exactly in a record of another species, or containing one
    (tax.contained), is 0 edits from that species: the two are identical
    wherever both have sequence, so a read that differs from the longer
    one only where the shorter has no sequence is no evidence for the
    shorter's species (the release's margin is 0 too). Only records in
    the database count.

    Adds: span and edits of the best record, gap (None unless the ties
    name exactly one species), next (record, edits further) for the
    nearest other species, ties (the tied records), tentative, and the
    gap rule's own call (gap_rule_depth, gap_rule_rank)."""
    ranks = ranks_of(tax)
    empty = {"rank": "", "depth": -1, "identity": 0.0, "coverage": 0.0,
             "best_hit": "", "near_hits": 0, "lca_depth": -1, "span": 0,
             "edits": None, "gap": None, "next": None, "ties": [],
             "tentative": False, "gap_rule_depth": -1, "gap_rule_rank": ""}
    scores = [x for x in scores if x[1] >= MIN_SPAN]
    if not scores:
        return empty
    top = max(1 - d / sp for d, sp, _c in scores)
    ties = [x for x in scores if 1 - x[0] / x[1] >= top - 1e-9]
    longest = max(x[1] for x in ties)
    ties = [x for x in ties if x[1] >= longest - SPAN_SLACK]
    d0, sp0, best = min(ties, key=lambda x: (x[0], -x[1], x[2]))
    species = {tax.get(c, {}).get("species", "") for _d, _s, c in ties}
    other = [(1 - d / sp, sp, c) for d, sp, c in scores
             if tax.get(c, {}).get("species", "") not in species
             and tax.get(c, {}).get("species_status", "shared") != "unnamed"]
    nxt = None
    if other:
        oi, _osp, oc = max(other)
        nxt = (oc, round((1 - oi) * sp0) - d0)
    gap = None
    if len(species) == 1 and "" not in species:
        gap = nxt[1] if nxt else sp0 - d0
        flagged = sorted({c for _d, _s, t in ties
                          for c in getattr(tax, "contained", {}).get(t, ())
                          if c in tax and tax[c].get("species")
                          and tax[c]["species"] not in species})
        if flagged and gap > 0:
            gap, nxt = 0, (flagged[0], 0)
    tied = sorted(c for _d, _s, c in ties)
    lca = lca_depth(tied, tax)
    ident = 1 - d0 / sp0
    depth = apply_floors(lca, ident, floors, ranks)
    k = ranks.index("species")
    clear = gap is not None and (
        gap >= gap_edits or (gap_exact > 0 and d0 == 0 and gap >= gap_exact))
    if depth >= k and (sp0 < min_species_span or
                       rule in ("gap", "tiered") and not clear):
        depth = k - 1
    gap_depth, tentative = depth, False
    if (rule == "tiered" and depth == k - 1 and lca >= k and clear
            and sp0 >= min_species_span and ident >= tentative_identity):
        depth, tentative = k, True
    return {"depth": depth, "rank": ranks[depth] if depth >= 0 else "",
            "identity": ident,
            "coverage": min(1.0, sp0 / qlen) if qlen else 1.0,
            "best_hit": best, "near_hits": len(tied), "lca_depth": lca,
            "span": sp0, "edits": d0, "gap": gap, "next": nxt,
            "ties": tied, "tentative": tentative, "gap_rule_depth": gap_depth,
            "gap_rule_rank": ranks[gap_depth] if gap_depth >= 0 else ""}


def record_label(rid: str, tax) -> str:
    """A record's species; for a shared record the species its members
    name ("Cooperia shared: punctata / spatulata"); else its deepest
    name and rank."""
    t = tax.get(rid, {})
    if t.get("species"):
        return t["species"]
    ranks = ranks_of(tax)
    rank, name = next(((r, t[r]) for r in reversed(ranks) if t.get(r)),
                      ("", ""))
    shared = getattr(tax, "shared", {}).get(rid)
    if shared and name:
        return f"{name} shared: " + " / ".join(
            s.split(" ", 1)[-1] for s in shared)
    return f"{name} ({rank})" if name else rid


def nearest_next(call: dict, tax) -> tuple[str, str]:
    """How close a span-scored call came: (nearest, next).

    nearest: the tied species and their identity ("Haemonchus contortus
    99.6%"); next: the nearest record of any other species and how many
    edits further it is ("Haemonchus placei +2"), or "none close" when
    minimap2 proposed none. A call stopped at genus then says whether it
    was "exactly X, but its sister is 2 edits away" (ITS2 cannot separate
    them) or "98% to X, nothing else close" (a variant the database does
    not hold, or a missing relative)."""
    if not call.get("span"):
        return "", ""
    shared = getattr(tax, "shared", {})
    names = sorted({n for t in call["ties"]
                    for n in (shared.get(t) or [record_label(t, tax)])})
    near = " / ".join(names[:3]) + (
        f" / +{len(names) - 3} more" if len(names) > 3 else "")
    nxt = call.get("next")
    return (f"{near} {100 * call['identity']:.1f}%",
            f"{record_label(nxt[0], tax)} +{nxt[1]}" if nxt else "none close")


def span_columns(call: dict, tax) -> dict:
    """The columns a span-scored call adds to clusters.tsv. Under the
    tiered rule, tentative says the species was named only by its lower
    identity floor, and gap_rule_rank / gap_rule_taxon are the gap rule's
    own call - the confident one - beside it."""
    near, nxt = nearest_next(call, tax)
    gap = call.get("gap")
    lin = tax.get(call.get("best_hit", ""), {})
    gd = call.get("gap_rule_depth", call.get("depth", -1))
    return {"span": call.get("span") or "",
            "edits": "" if call.get("edits") is None else call["edits"],
            "gap": "" if gap is None else gap,
            "nearest": near, "next": nxt,
            "grade": lin.get("grade", ""),
            "tentative": "yes" if call.get("tentative") else "",
            "gap_rule_rank": call.get("gap_rule_rank", call.get("rank", "")),
            "gap_rule_taxon": lin.get(ranks_of(tax)[gd], "") if gd >= 0 else ""}


SPAN_COLUMNS = ["span", "edits", "gap", "nearest", "next", "grade",
                "tentative", "gap_rule_rank", "gap_rule_taxon"]


MERGE_MIN_OVERLAP = 20
MERGE_MAX_MISMATCH = 0.05


def merge_pair(s1: str, q1: str, s2: str, q2: str):
    """(insert, quality) of a pair whose mates overlap, or None.

    The overlap is R1's end against R2's reverse complement's start, found
    from a 16-base seed at the start of R2's reverse complement (or the
    next 16, if the first holds an error) and accepted at
    MERGE_MIN_OVERLAP or more bases with at most MERGE_MAX_MISMATCH of
    them mismatched; the longest such overlap wins. In the overlap R1's
    bases are kept. As DADA2's mergePairs does in spirit, and as
    nemabiome-its2-refdb's its2_classify.py merges."""
    r = revcomp(s2)
    rq = q2[::-1]
    best = None
    for off in (0, 16):
        seed = r[off:off + 16]
        if len(seed) < 16:
            break
        i = s1.find(seed)
        while i >= 0:
            ov = len(s1) - i + off
            if MERGE_MIN_OVERLAP <= ov <= min(len(s1), len(r)) and \
                    (best is None or ov > best):
                a, b = s1[-ov:], r[:ov]
                if sum(x != y for x, y in zip(a, b)) <= MERGE_MAX_MISMATCH * ov:
                    best = ov
            i = s1.find(seed, i + 1)
        if best is not None:
            break
    if best is None:
        return None
    return s1 + r[best:], q1 + rq[best:]

class Nemits2(Protocol):
    name = "nemits2"
    title = ("nemabiome-its2-refdb, ITS2: span scoring, the gap rule and "
             "its tentative tier")
    read_types = {
        "single_long_noisy": "reference-guided consensus per group, each "
                             "candidate record scored over its span",
        "paired_short": "overlapping mates merged, others two windows "
                        "(R2 with its 28S stretch), scored over each "
                        "record's span",
        "single_long_accurate": "the consensus route, as nanopore's, each "
                                "candidate record scored over its span",
    }
    not_supported: dict = {}
    calling_rules = ("tiered", "gap", "floors")
    required_files = ("refdb.fasta", "taxonomy.tsv",
                      "resolution_nanopore.tsv")
    extra_columns = tuple(SPAN_COLUMNS)
    pair_fates = ("merged",)
    pair_note = (" merged counts the pairs (within clustered or "
                 "unclustered) whose mates overlapped and were merged "
                 "into the insert.")

    def release_problems(self, release, manifest, read_types=()):
        out = super().release_problems(release, manifest, read_types)
        if ("paired_short" in read_types and scores_over_reference(manifest)
                and not manifest.get("illumina_geometries")):
            out.append(("warn", (
                f"protocol {self.name}: the release lists no "
                "illumina_geometries, so the 28S stretch R2 reads before "
                "each record cannot be measured, and R2's window takes "
                "none of it")))
        return out

    def calling(self, args, manifest, tax):
        rule = super().calling(args, manifest, tax)[0]
        return rule, parse_floors(
            args.gap_floors if rule in ("gap", "tiered") else args.floors,
            ranks_of(tax))

    def calling_settings(self, rule, args):
        gapped = rule in ("gap", "tiered")
        return {"score_over": "reference span", "calling_rule": rule,
                "gap_edits": args.gap_edits if gapped else None,
                "gap_exact": args.gap_exact if gapped else None,
                "tentative_identity": args.tentative_identity
                if rule == "tiered" else None,
                "min_span": MIN_SPAN, "min_species_span": args.min_species_span,
                "span_slack": SPAN_SLACK}

    def column_values(self, call, tax):
        return span_columns(call, tax)

    def classify_consensus(self, ctx, seqs: dict) -> dict:
        # Scored over each record's span (span_call): the consensus
        # carries the 28S stretch its records lack.
        hits = consensus_hits(ctx.mm, ctx.index, ctx.tmp, seqs)
        cands = {k: span_candidates(hits.get(f"c{k}", [])) for k in seqs}
        need = {t for ts in cands.values() for t in ts}
        recs = {i: s for i, s in read_fasta(ctx.db_fasta) if i in need}
        calls = {}
        for k, cons in seqs.items():
            calls[k] = span_call(
                span_scores([(cons, "whole")], cands[k], recs),
                ctx.tax, ctx.rule, ctx.floors, ctx.args.gap_edits,
                len(cons), ctx.args.min_species_span, ctx.args.gap_exact,
                ctx.args.tentative_identity)
        return calls

    def tail(self, manifest, read_length=None):
        """Bases a primer-trimmed R2 reads before its record's span starts.

        0 for a release scored over whole inserts. For one scored over its
        records' span (ITS2), R2 opens with the 28S stretch the records stop
        short of, and the release's windows count record bases only (130 +
        83 at 2x150): a geometry's read length, less its reverse primer and
        its R2 window, is that stretch (150 - 20 - 83 = 47)."""
        if not scores_over_reference(manifest):
            return 0
        rev = len(((manifest.get("primers") or {}).get("reverse") or {})
                  .get("sequence", ""))
        geoms = [g for g in manifest.get("illumina_geometries") or []
                 if isinstance(g, dict) and len(g.get("windows") or []) == 2
                 and g.get("read_length")]
        if not geoms or not rev:
            return 0
        fits = [g for g in geoms
                if read_length and int(g["read_length"]) <= read_length]
        g = max(fits, key=lambda g: int(g["read_length"])) if fits else \
            min(geoms, key=lambda g: int(g["read_length"]))
        return max(0, int(g["read_length"]) - rev - int(g["windows"][1]))

    def window_settings(self, rw1, rw2, args):
        return {"read_windows": [rw1, rw2],
                "primer_offset": args.primer_offset}

    def merge_mates(self, s1, q1, s2, q2, args):
        return None if args.no_merge else merge_pair(s1, q1, s2, q2)

    def classify_pairs(self, ctx, centroids: list) -> tuple[dict, int]:
        args, tax, rule, floors, tail = (ctx.args, ctx.tax, ctx.rule,
                                         ctx.floors, ctx.tail)
        calls: dict[int, dict] = {}
        # R2's 28S stretch is unaligned in every pair, so these hits only
        # propose candidates; none is mapped again as on V4.
        hits, _paths, _extra = pair_hits(ctx.mm, ctx.index, ctx.tmp,
                                         centroids, args.max_secondary)
        # Every record either window proposed, scored exactly over
        # its span with both windows (span_scores): a pair is one
        # sequence with a gap in it, and a record that minimap2 saw
        # through one window only is still scored on both.
        cands = {}
        for k, (x, y, _n) in enumerate(centroids):
            q = f"p{k}"
            hs = [h for name in (("r1", "r2") if y else ("full",))
                  for h in hits[name].get(q, {}).values()
                  if h["strand"] == ("-" if name == "r2" else "+")]
            if y:
                # Ranked on both windows together: matches over the
                # aligned blocks (a window's record runs past it, so
                # the record's unaligned ends say nothing), then score.
                both: dict = defaultdict(lambda: [0, 0, 0])
                for h in hs:
                    b = both[h["t"]]
                    b[0] += h["match"]
                    b[1] += h["alen"]
                    b[2] += h["AS"]
                hs = [{"t": t, "ident": m / max(1, a), "AS": sc}
                      for t, (m, a, sc) in both.items()]
            cands[k] = span_candidates(hs)
        need = {t for ts in cands.values() for t in ts}
        recs = {i: s for i, s in read_fasta(
            os.path.join(args.database, "refdb.fasta")) if i in need}
        for k, (x, y, _n) in enumerate(centroids):
            parts = [(x, "start"), (revcomp(y), "end")] if y else \
                [(x, "whole")]
            scorable = len(x) + max(0, len(y) - tail) if y else \
                len(x) - tail
            calls[k] = span_call(
                span_scores(parts, cands[k], recs, tail), tax, rule,
                floors, args.gap_edits, scorable, args.min_species_span,
                args.gap_exact, args.tentative_identity)
        return calls, 0


PROTOCOL = Nemits2()
