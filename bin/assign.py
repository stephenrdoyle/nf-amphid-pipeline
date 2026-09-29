#!/usr/bin/env python3
"""
Competitive assignment for nf-helminth-amplicon.

Both platforms converge here. A read is scored against the whole panel at
once - targets, tagged control decoys, quantification spike-ins and vector
backbone in one index - and assigned to whichever it fits best, or to
nothing. Scoring against targets alone would give every read a best hit
whether or not it came from the panel, which is the failure this assay is
built to avoid: a plasmid carryover read has to be able to say so.

Three things decide the call.

  distance   Bounded edit distance to each candidate reference. Exact, not
             a similarity heuristic: the panel contains haplotypes three
             bases apart in three hundred, and an approximate score cannot
             be trusted to order those correctly.

  margin     How far the best fit beats the runner-up. A read that fits two
             references equally well has not identified either, and saying
             so is more useful than breaking the tie arbitrarily. Where the
             tied references share a taxon the species call still stands -
             only the haplotype is unresolved.

  tag        The 18bp insert that every control plasmid carries and no
             target does. This is a presence test, not a comparison, so it
             is independent of the alignment: a read carrying the tag is
             control-derived whatever it mapped to. Where the two disagree,
             that disagreement is reported rather than silently resolved -
             it is the sharpest evidence available that something in the
             panel or the library is not what it is believed to be.

Matching is staged for cost, not for accuracy: an exact lookup catches the
error-free majority, a k-mer prefilter narrows the field to a few plausible
references, and the exact distance is computed only against those. Every
read that reaches a call has been scored exactly.

Read accounting is preserved: each input sequence lands in exactly one
outcome, weighted by its dereplicated count, and the outcomes sum to the
input.

Standard library only.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    bounded_edit,
    find_anchor,
    load_assay,
    load_primers,
    open_reads,
    read_table,
    revcomp,
)

OUTCOMES = (
    "assigned",
    "taxon_only",
    "ambiguous",
    "unmatched",
    "tag_unobservable",
    "no_primer",
)

OUTCOME_MEANING = {
    "assigned": "one reference fits best by a clear margin",
    "taxon_only": "haplotype unresolved, but every candidate is one taxon",
    "ambiguous": "two or more references fit equally well across taxa",
    "unmatched": "no reference within the distance limit",
    "tag_unobservable": (
        "target and control fit equally - the tag lies between the mates, "
        "where neither read can see it"
    ),
    "no_primer": "no primer at the read start (nanopore input only)",
}

SPECIATION = "speciation"


# --------------------------------------------------------------------------
# the panel
# --------------------------------------------------------------------------


@dataclass
class Ref:
    seq_id: str
    sequence: str
    locus: str = ""
    type: str = "target"
    taxon: str = ""
    purpose: str = SPECIATION
    parent_id: str = ""


class Panel:
    """The reference panel, indexed three ways.

    The k-mer index is a filter, never a score. It decides which references
    are worth the cost of an exact comparison and nothing else, so a read
    whose k-mer profile is misleading is still measured properly against
    the candidates it reaches.
    """

    def __init__(self, refs: list[Ref], k: int = 15):
        self.refs = refs
        self.k = k
        self.exact: dict[str, int] = {}
        self.kmers: dict[str, list[int]] = defaultdict(list)
        for i, r in enumerate(refs):
            self.exact.setdefault(r.sequence, i)
            for j in range(0, len(r.sequence) - k + 1):
                self.kmers[r.sequence[j : j + k]].append(i)
        self.by_id = {r.seq_id: i for i, r in enumerate(refs)}

    def candidates(self, q: str, keep: int = 4) -> list[int]:
        k = self.k
        if len(q) < k:
            return list(range(len(self.refs)))
        hits: Counter = Counter()
        for j in range(0, len(q) - k + 1):
            for i in self.kmers.get(q[j : j + k], ()):
                hits[i] += 1
        if not hits:
            return []
        top = hits.most_common(1)[0][1]
        # Half the best count is a deliberately loose gate. The filter's job
        # is to avoid comparing against obviously unrelated references, not
        # to pick a winner - anything it discards is never measured, so it
        # must err toward keeping.
        return [i for i, n in hits.most_common(keep) if n >= top * 0.5]


def load_panel(fasta: str, meta_path: str | None) -> list[Ref]:
    seqs: dict[str, str] = {}
    key = None
    with open(fasta, encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith(">"):
                key = ln[1:].split()[0]
                seqs[key] = ""
            elif key:
                seqs[key] += ln.strip().upper()

    meta: dict[str, dict] = {}
    if meta_path:
        for r in read_table(meta_path):
            if r.get("seq_id"):
                meta[r["seq_id"]] = r

    refs = []
    for sid, seq in seqs.items():
        m = meta.get(sid, {})
        purpose = (m.get("purpose") or "").strip() or SPECIATION
        refs.append(
            Ref(
                seq_id=sid,
                sequence=seq,
                locus=m.get("locus", ""),
                type=m.get("type", "target"),
                taxon=(m.get("taxon") or sid),
                purpose=purpose,
                parent_id=m.get("parent_id", ""),
            )
        )
    # Controls report under their parent's taxon: a control read is evidence
    # about the plasmid, not about the species, and letting it inherit the
    # taxon keeps the two countable separately while still naming what it
    # is a control for.
    by_id = {r.seq_id: r for r in refs}
    for r in refs:
        if r.type == "control" and r.parent_id in by_id:
            r.taxon = by_id[r.parent_id].taxon
            r.purpose = by_id[r.parent_id].purpose
    return refs


# --------------------------------------------------------------------------
# calling
# --------------------------------------------------------------------------


@dataclass
class Call:
    outcome: str
    seq_id: str = ""
    taxon: str = ""
    distance: int = -1
    margin: int = -1
    tag: bool = False
    tied: tuple = ()
    ref_type: str = ""
    purpose: str = ""


def call_sequence(
    q: str, panel: Panel, tag: str, max_frac: float, min_margin: int,
    tag_max_edits: int, second: str | None = None,
) -> Call:
    """Assign one (dereplicated) sequence.

    `second` is the reverse mate of an unmerged pair. Both mates are scored
    against the same reference and their distances summed, so a pair only
    counts for a reference that explains both ends. One end can be wrong;
    both ends being wrong in the same direction is rare.
    """
    # Both mates are searched. A control differs from its target only by the
    # tag, so testing one mate and calling the other's silence evidence is
    # how a control read gets reported as its own target.
    has_tag = bool(tag) and (
        find_anchor(q, tag, tag_max_edits) is not None
        or (
            second is not None
            and find_anchor(revcomp(second), tag, tag_max_edits) is not None
        )
    )

    exact = panel.exact.get(q) if second is None else None
    if exact is not None:
        r = panel.refs[exact]
        # Even an exact hit needs its runner-up: a target and its control
        # differ by an insertion, so an exact hit on one is not by itself
        # evidence against the other.
        return _finish(panel, [(0, exact)], has_tag, min_margin, q, tag)

    max_d = max(2, int(len(q) * max_frac))
    scored: list[tuple[int, int]] = []
    for i in panel.candidates(q):
        ref = panel.refs[i]
        if second is None:
            d = bounded_edit(q, ref.sequence, max_d)
        else:
            d1 = bounded_edit(q, ref.sequence, max_d, free_b_suffix=True)
            if d1 > max_d:
                continue
            # The reverse mate covers the far end of the amplicon, so both
            # strings are reversed and the same prefix-anchored comparison
            # is used - the reference is then free to run past the mate at
            # the end it never reached.
            d2 = bounded_edit(
                revcomp(second)[::-1], ref.sequence[::-1], max_d,
                free_b_suffix=True,
            )
            d = d1 + d2
        if d <= max_d * (2 if second is not None else 1):
            scored.append((d, i))

    if not scored:
        return Call(outcome="unmatched", tag=has_tag)
    scored.sort()
    return _finish(panel, scored, has_tag, min_margin, q, tag)


def _finish(panel: Panel, scored, has_tag: bool, min_margin: int, q: str,
            tag: str) -> Call:
    best_d, best_i = scored[0]
    tied = [i for d, i in scored if d - best_d < min_margin]

    # The tag is deterministic where the alignment is not, but only in one
    # direction. Finding it proves the read is control-derived. NOT finding
    # it proves nothing unless the read actually covered the place the tag
    # would be - and a tie between a target and its control is itself the
    # proof that it did not, since an observed tag would have separated them
    # by its own length. Treating that silence as evidence of target is how
    # plasmid carryover gets reported as infection.
    mixed = len({panel.refs[i].type == "control" for i in tied}) > 1
    if len(tied) > 1 and tag and has_tag:
        narrowed = [i for i in tied if panel.refs[i].type == "control"]
        if narrowed:
            tied, best_i, mixed = narrowed, narrowed[0], False

    runner = next((d for d, i in scored if i not in tied), None)
    margin = -1 if runner is None else runner - best_d

    if len(tied) == 1:
        r = panel.refs[best_i]
        return Call("assigned", r.seq_id, r.taxon, best_d, margin, has_tag,
                    (), r.type, r.purpose)

    if mixed:
        return Call(
            "tag_unobservable", "", "", best_d, margin, has_tag,
            tuple(sorted(panel.refs[i].seq_id for i in tied)),
        )

    taxa = {panel.refs[i].taxon for i in tied}
    types = {panel.refs[i].type for i in tied}
    if len(taxa) == 1 and len(types) == 1:
        r = panel.refs[tied[0]]
        return Call(
            "taxon_only", "", r.taxon, best_d, margin, has_tag,
            tuple(sorted(panel.refs[i].seq_id for i in tied)),
            r.type, r.purpose,
        )
    return Call(
        "ambiguous", "", "", best_d, margin, has_tag,
        tuple(sorted(panel.refs[i].seq_id for i in tied)),
    )


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------


def read_derep_fasta(path: str) -> list[tuple[str, int]]:
    op = gzip.open if path.endswith(".gz") else open
    out, seq, size = [], [], 0
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if seq:
                    out.append(("".join(seq), size))
                seq = []
                size = int(ln.split("size=")[1]) if "size=" in ln else 1
            else:
                seq.append(ln.strip().upper())
    if seq:
        out.append(("".join(seq), size))
    return out


def read_pairs(path: str) -> list[tuple[str, str, int]]:
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        head = fh.readline()
        if not head:
            return out
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 4:
                out.append((parts[2], parts[3], int(parts[1])))
    return out


def derep_fastq(path: str, ps, cfg: dict) -> tuple[list[tuple[str, int]], int]:
    """Nanopore input: reads still carry their locus primers, because the
    demultiplexer trims the construct and the amplicon begins at the primer.
    They are trimmed here so both platforms present the same extent to the
    panel."""
    counts: Counter = Counter()
    no_primer = 0
    for _rid, seq, _q in open_reads(path):
        if ps is not None:
            f = ps.match(seq, "fwd")
            rc = revcomp(seq)
            r = ps.match(rc, "rev")
            if f is None or r is None:
                no_primer += 1
                continue
            # The reverse primer was located on the reverse complement, so
            # its end offset there is how much to remove from this end.
            seq = seq[f[1] : len(seq) - r[1]]
        counts[seq] += 1
    return list(counts.items()), no_primer


# --------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", required=True,
                    help="reference.insert.fasta from pre-flight")
    ap.add_argument("--reference", help="reference panel TSV, for metadata")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--merged", help="dereplicated merged FASTA (Illumina)")
    ap.add_argument("--pairs", help="unmerged pair table (Illumina)")
    ap.add_argument("--reads", help="per-sample FASTQ (nanopore)")
    ap.add_argument("--primers", help="required for --reads")
    ap.add_argument("--assay")
    ap.add_argument("--outdir", default="assign")
    ap.add_argument("--max-dist-frac", type=float, default=0.10,
                    help="distance limit as a fraction of read length")
    ap.add_argument("--min-margin", type=int, default=2,
                    help="edits by which the best fit must beat the next")
    ap.add_argument("--tag-max-edits", type=int, default=3)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    _T0 = time.time()

    if not (args.merged or args.reads):
        raise SystemExit("one of --merged (Illumina) or --reads (nanopore)")

    cfg = load_assay(args.assay)
    tag, tag_edits = "", args.tag_max_edits
    if args.assay and os.path.exists(args.assay):
        with open(args.assay, encoding="utf-8") as fh:
            ct = json.load(fh).get("control_tag", {}) or {}
        tag = ct.get("sequence", "") or ""
        # The assay owns the tolerance; the flag exists to override it, not
        # to shadow it. A tag matched at the wrong tolerance is the one
        # error that would make every control call quietly unreliable.
        if ct.get("max_edits") is not None and "--tag-max-edits" not in (
            argv if argv is not None else sys.argv[1:]
        ):
            tag_edits = int(ct["max_edits"])
    if not tag:
        raise SystemExit(
            "no control tag in the assay definition - target and control "
            "cannot be separated without it"
        )
    refs = load_panel(args.panel, args.reference)
    if not refs:
        raise SystemExit(f"panel is empty: {args.panel}")
    panel = Panel(refs)

    ps = None
    if args.reads:
        if not args.primers:
            raise SystemExit("--reads requires --primers")
        ps = load_primers(args.primers)

    items: list[tuple[str, str | None, int]] = []
    no_primer = 0
    if args.merged:
        items += [(s, None, n) for s, n in read_derep_fasta(args.merged)]
    if args.pairs and os.path.exists(args.pairs):
        items += [(a, b, n) for a, b, n in read_pairs(args.pairs)]
    if args.reads:
        uniq, no_primer = derep_fastq(args.reads, ps, cfg)
        items += [(s, None, n) for s, n in uniq]

    counts: Counter = Counter()
    outcomes: Counter = Counter()
    unique_by_ref: Counter = Counter()
    taxon_counts: Counter = Counter()
    # Composition and spike-ins are different measurements and must not be
    # added together. A control read is evidence about plasmid carryover
    # and a quant read is a dose of known concentration; neither is a
    # statement about what was in the stool.
    by_type: Counter = Counter()
    tag_vs_map: Counter = Counter()
    dists: list[int] = []
    tied_sets: Counter = Counter()

    for seq, mate, n in items:
        c = call_sequence(
            seq, panel, tag, args.max_dist_frac, args.min_margin,
            tag_edits, second=mate,
        )
        outcomes[c.outcome] += n
        if c.outcome == "assigned":
            counts[c.seq_id] += n
            unique_by_ref[c.seq_id] += 1
            by_type[c.ref_type] += n
            if c.ref_type == "target" and SPECIATION in c.purpose:
                taxon_counts[c.taxon] += n
            dists.append(c.distance)
            is_ctrl = panel.refs[panel.by_id[c.seq_id]].type == "control"
            tag_vs_map[
                ("tagged" if c.tag else "untagged")
                + " / "
                + ("control" if is_ctrl else "non-control")
            ] += n
        elif c.outcome == "taxon_only":
            by_type[c.ref_type] += n
            if c.ref_type == "target" and SPECIATION in c.purpose:
                taxon_counts[c.taxon] += n
            tied_sets[c.tied] += n
        elif c.outcome in ("ambiguous", "tag_unobservable"):
            tied_sets[c.tied] += n
    outcomes["no_primer"] += no_primer

    os.makedirs(args.outdir, exist_ok=True)
    sid = args.sample_id

    with open(os.path.join(args.outdir, f"{sid}.counts.tsv"), "w") as fh:
        fh.write("sample_id\tseq_id\tlocus\ttype\tpurpose\ttaxon\treads\t"
                 "unique_seqs\n")
        for r in refs:
            n = counts.get(r.seq_id, 0)
            fh.write(
                f"{sid}\t{r.seq_id}\t{r.locus}\t{r.type}\t{r.purpose}\t"
                f"{r.taxon}\t{n}\t{unique_by_ref.get(r.seq_id, 0)}\n"
            )

    # Composition is over speciation targets only. A panel carrying a second
    # locus for variant calling would otherwise count every species twice,
    # with a denominator mixing loci at different copy numbers.
    denom_ids = {
        r.seq_id for r in refs
        if r.type == "target" and SPECIATION in r.purpose
    }
    # The denominator is the reads in this table, not every assigned read -
    # so the column is a composition and sums to one, rather than a set of
    # fractions of something unstated.
    denom = sum(taxon_counts.values())
    with open(os.path.join(args.outdir, f"{sid}.taxa.tsv"), "w") as fh:
        fh.write("sample_id\ttaxon\treads\trel_abundance\n")
        for tx, n in sorted(taxon_counts.items(), key=lambda kv: -kv[1]):
            rel = f"{n / denom:.6f}" if denom else ""
            fh.write(f"{sid}\t{tx}\t{n}\t{rel}\n")

    total = sum(n for _s, _m, n in items) + no_primer
    # Wall clock per sample, and the read length that drives the cost.
    # Assignment aligns within a band proportional to read length, so its
    # cost grows as length squared - which makes "how long did this take"
    # uninterpretable without "how long were the reads". Recorded together
    # so a slow run can be explained rather than guessed at.
    lengths = sorted(len(i[0]) for i in items) if items else [0]
    summary = {
        "elapsed_s": round(time.time() - _T0, 2),
        "read_length_median": lengths[len(lengths) // 2],
        "read_length_p95": lengths[min(len(lengths) - 1, 19 * len(lengths) // 20)],
        "sample_id": sid,
        "sequences_in": total,
        "unique_in": len(items),
        "outcomes": {k: outcomes.get(k, 0) for k in OUTCOMES},
        "outcome_meaning": OUTCOME_MEANING,
        "composition_denominator": denom,
        "reads_by_reference_type": dict(by_type),
        "median_distance": sorted(dists)[len(dists) // 2] if dists else None,
        # A tagged read that mapped to a non-control, or an untagged read
        # that mapped to a control, is the panel contradicting itself.
        "tag_vs_mapping": dict(tag_vs_map),
        "tag_mapping_conflicts": sum(
            v for k, v in tag_vs_map.items()
            if k in ("tagged / non-control", "untagged / control")
        ),
        "top_ambiguous_sets": [
            {"references": list(k), "reads": v}
            for k, v in tied_sets.most_common(5)
        ],
        "accounting_balances": total == sum(outcomes.get(k, 0) for k in OUTCOMES),
    }
    with open(os.path.join(args.outdir, f"{sid}.assign.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    if not args.quiet:
        print(f"{sid}: {total:,} sequences ({len(items):,} unique)")
        for k in OUTCOMES:
            n = outcomes.get(k, 0)
            if n:
                print(f"    {k:12s} {n:8,}  {100*n/total:5.2f}%  "
                      f"{OUTCOME_MEANING[k]}")
        if summary["tag_mapping_conflicts"]:
            print(f"    tag/mapping conflicts: "
                  f"{summary['tag_mapping_conflicts']:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
