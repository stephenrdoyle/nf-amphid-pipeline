#!/usr/bin/env python3
"""
Illumina read preparation for amphid.

The sequencer has already demultiplexed, so unlike the nanopore path this
step never decides which sample a read came from. What it does decide is
what counts as the read: where the primer stops and the template starts,
whether the two mates describe the same molecule, and whether the pair is
good enough to keep.

Five things happen, in this order:

  1. orient       - R1 should open with a forward primer and R2 with a
                    reverse one. If they are the other way round the pair is
                    swapped rather than discarded.
  2. trim primers - the primer bases are synthesised, not observed. Whatever
                    the template carried under a primer, the read reports
                    the primer, so those bases are removed here and the
                    reference is trimmed to match (preflight writes
                    reference.insert.fasta for exactly this).
  3. read-through - when the insert is shorter than the read, the read runs
                    off the end of the insert, through the opposite primer
                    and into the opposite adapter. Both are removed, so the
                    read stops where the insert does.
  4. merge        - overlapping mates are collapsed into one sequence, with
                    disagreements resolved by quality.
  5. dereplicate  - identical sequences are collapsed to one record with a
                    count, which is what makes competitive assignment
                    affordable downstream.

Pairs whose amplicon is longer than the two reads together cannot merge.
They are not discarded: they are written to a separate table and assigned
as a pair, with both mates required to agree. That is the same concordance
requirement the nanopore demultiplexer applies to its two ends, for the same
reason - one end can be wrong, two ends being wrong the same way is rare.

Read accounting is preserved: every input pair ends in exactly one bucket,
and the buckets sum to the input. The input itself is checked too: R1 and
R2 must hold the same reads in the same order and end together, or the run
stops - pairs lost before the funnel are invisible to the funnel's sums.

Standard library only.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    PrimerSet,
    expected_errors,
    find_anchor_cost,
    load_assay,
    load_primers,
    open_reads,
    revcomp,
)

# How a pair can end up. Every input pair lands in exactly one of these.
FATES = (
    "merged",
    "paired",
    "no_primer",
    "too_short",
    "low_quality",
    # Never produced. Kept so every prep.json carries the same keys, which
    # synthesise.py reads by name; see FATE_MEANING for why.
    "length_mismatch",
)

FATE_MEANING = {
    "merged": "mates overlapped and were collapsed into one sequence",
    "paired": "no overlap - kept as a pair, both mates must agree later",
    "no_primer": "neither orientation opened with a primer pair",
    "too_short": "insert shorter than the minimum after trimming",
    "low_quality": "expected errors above the threshold",
    "length_mismatch": (
        "reserved, always 0 - mates whose overlap disagrees are kept as "
        "pairs, because a disagreeing seed hit on a long amplicon is far "
        "more often a chance match than evidence the pair is bad"
    ),
}


# --------------------------------------------------------------------------
# trimming and merging
# --------------------------------------------------------------------------


# Complement of every IUPAC code, so a degenerate primer can be turned
# round. ampseqlib.revcomp maps ACGTN only and would leave a Y in place
# where an R belongs - the right letter for the forward strand and the
# wrong one for the reverse, which is exactly the strand needed here.
_IUPAC_COMP = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def revcomp_iupac(seq: str) -> str:
    return seq.translate(_IUPAC_COMP)[::-1]


# The shortest read tail accepted as the START of the read-through junction
# when the read stops before the junction is complete. Ten bases with the
# error allowance below is a chance match about once in 10^5 reads, and
# what such a chance match costs is a few bases off the end of one mate -
# bases the other mate still covers wherever the pair overlaps.
MIN_PARTIAL_JUNCTION = 10


def trim_readthrough(
    seq: str, qual: str, adapter: str, opposite_primer: str = "",
    max_edits: int = 3, primer_edits: int = 2,
):
    """Remove everything past the end of the insert when the read runs off it.

    When the insert is shorter than the read, the read does not run
    straight from template into adapter. The molecule is

        fwd primer - insert - rc(rev primer) - rc(rev overhang) - ...

    so R1 runs through the reverse-complemented OPPOSITE primer before it
    reaches any adapter (and R2 likewise through rc(fwd primer)). Cutting
    at the adapter alone leaves that primer on the read: the merge seed
    then lands on primer rather than insert, the pair falls through to
    `paired`, and assignment charges each mate twenty-odd edits for
    synthesised bases - enough to leave the shortest amplicons unmatched.

    So the thing searched for is the junction itself: rc(opposite primer)
    followed by the first 20bp of the adapter, and the read is cut where
    the junction starts, i.e. where the insert ends. Which opposite primer
    is known - it is the one matched at the start of the other mate, the
    same molecule seen from the other end - so the junction is one
    sequence, not a search over the primer set. Searching for primer and
    adapter together is also more specific than either alone: ~40bp with
    the adapter's and the primer's edit allowance combined.

    A read can also stop part way into the junction - in the primer, or a
    few bases into the adapter - in which case the whole junction is not
    there to be found. The tail of the read is then tested against the
    START of the junction, held to the same error rate per base (not the
    same absolute allowance, which over ten bases would accept almost
    anything) and at least MIN_PARTIAL_JUNCTION bases long. A read that
    stops fewer than that many bases past the insert keeps them; that is
    a handful of bases, and the alternative - trimming any read whose last
    few bases happen to resemble a primer - would shorten reads that never
    ran through at all.

    Without `opposite_primer` this is the old adapter-only search, with
    the same partial-tail rule.
    """
    junction = revcomp_iupac(opposite_primer) + adapter[:20]
    # Myers' matcher packs the pattern into one 64-bit word.
    junction = junction[:64]
    m = len(junction)
    if not m or not seq:
        return seq, qual
    budget = max_edits + (primer_edits if opposite_primer else 0)

    hit = find_anchor_cost(seq, junction, budget)
    if hit is not None:
        start = hit[0]
        return seq[:start], qual[:start]

    # Partial junction at the read end. Only the last m bases can hold it,
    # and the read is padded with N - which the matcher accepts as any
    # base - so the part of the junction beyond the read end costs nothing.
    # The padding stops MIN_PARTIAL_JUNCTION short of the full junction, so
    # at least that many real bases must take part.
    tail_from = max(0, len(seq) - m)
    tail = seq[tail_from:]
    pad = m - MIN_PARTIAL_JUNCTION
    if pad <= 0 or len(tail) < MIN_PARTIAL_JUNCTION:
        return seq, qual
    hit = find_anchor_cost(tail + "N" * pad, junction, budget)
    if hit is None:
        return seq, qual
    start, end, cost = hit
    if end <= len(tail):
        # A complete junction inside the tail would have been found above;
        # this is a weaker duplicate of that search, not a partial.
        return seq, qual
    covered = len(tail) - start
    if covered < MIN_PARTIAL_JUNCTION or cost > (budget * covered) // m:
        return seq, qual
    cut = tail_from + start
    return seq[:cut], qual[:cut]


# The merge seed: rc(R2)'s first SEED_LEN bases located in R1 within
# SEED_EDITS. Overlaps shorter than the seed are tested at the same error
# rate per base (see _overlap_offsets).
SEED_LEN = 24
SEED_EDITS = 3


def _overlap_offsets(r1: str, rc2: str, min_overlap: int):
    """Candidate offsets of rc(R2) against R1, most likely first.

    An offset k >= 0 means rc(R2) starts at R1[k] - the ordinary overlap.
    k < 0 means rc(R2) starts -k bases BEFORE R1 does: a dovetail, which
    only happens when both reads ran a few bases past the insert into the
    junction (too few for trim_readthrough to recognise safely), so each
    mate's overhang beyond the other is read-through, not template.

    Three searches, cheapest-to-trust first:

      seed      rc2[:24] found anywhere in R1 - every overlap of 24bp or
                more, with indels tolerated in the seed.
      short     overlaps from 23bp down to min_overlap, tested directly at
                the end of R1. The seed has to fit inside R1, so it can
                never find these; without this --min-overlap quietly meant
                max(24, min_overlap). Worse, an overlap of 21-23bp could be
                "found" by the seed with its last bases scored as
                deletions, at an offset a few bases off - which then fails
                verification - so this is tried even when the seed hit.
                There are at most a few lengths to try, each one ungapped
                comparison, so it costs almost nothing.
      dovetail  R1's first 24bp found inside rc(R2).

    The allowance everywhere is the seed's own error rate, 3 in 24, scaled
    to the length compared. Each candidate is still verified in full by the
    caller, so a wrong candidate costs a comparison, not a bad merge.
    """
    probe = rc2[: min(SEED_LEN, len(rc2))]
    hit = find_anchor_cost(r1, probe, SEED_EDITS)
    if hit is not None:
        yield hit[0]
    for L in range(min(len(probe) - 1, len(r1), len(rc2)), min_overlap - 1, -1):
        mism = sum(1 for x, y in zip(r1[len(r1) - L:], rc2[:L]) if x != y)
        if mism <= (SEED_EDITS * L) // SEED_LEN:
            yield len(r1) - L
            break
    if hit is None:
        hit = find_anchor_cost(rc2, r1[: min(SEED_LEN, len(r1))], SEED_EDITS)
        if hit is not None and hit[0] > 0:
            yield -hit[0]


def merge_pair(
    r1: str, q1: str, r2: str, q2: str, min_overlap: int = 20,
    max_diff_frac: float = 0.25,
) -> tuple[str, str] | None:
    """Collapse overlapping mates, or return None if they do not overlap.

    The overlap is found by locating the start of the reverse-complemented
    R2 inside R1 with the shared approximate matcher, then verifying the
    whole overlap ungapped. Approximate for the seed, exact-offset for the
    verification: Illumina indels are rare enough that a gapped overlap is
    better left unmerged than patched over, and a pair that fails here is
    kept rather than thrown away. The candidate offsets, including overlaps
    too short for the seed and dovetailed mates, are in _overlap_offsets.

    Disagreements are resolved by quality, and the merged base takes the
    higher score only where the mates agree. Where they disagree the score
    is knocked down to the difference, so a contested base cannot inherit
    the confidence of the mate that happened to win.
    """
    rc2, rq2 = revcomp(r2), q2[::-1]
    if len(r1) < min_overlap or len(rc2) < min_overlap:
        return None

    for k in _overlap_offsets(r1, rc2, min_overlap):
        # Align the two on molecule coordinates: R1 from i1, rc(R2) from i2.
        i1, i2 = (k, 0) if k >= 0 else (0, -k)
        overlap = min(len(r1) - i1, len(rc2) - i2)
        if overlap < min_overlap:
            continue
        a, b = r1[i1 : i1 + overlap], rc2[i2 : i2 + overlap]
        diff = sum(1 for x, y in zip(a, b) if x != y)
        if diff > max_diff_frac * overlap:
            continue
        break
    else:
        return None

    cons, consq = [], []
    for i in range(overlap):
        x, y = a[i], b[i]
        qx, qy = ord(q1[i1 + i]) - 33, ord(rq2[i2 + i]) - 33
        if x == y:
            cons.append(x)
            consq.append(chr(33 + min(41, qx + qy)))
        elif qx >= qy:
            cons.append(x)
            consq.append(chr(33 + max(2, qx - qy)))
        else:
            cons.append(y)
            consq.append(chr(33 + max(2, qy - qx)))
    if k < 0:
        # Dovetail: the molecule is exactly the overlap. What either mate
        # has beyond it is junction it ran into, so both overhangs go.
        return "".join(cons), "".join(consq)
    return (
        r1[:k] + "".join(cons) + rc2[overlap:],
        q1[:k] + "".join(consq) + rq2[overlap:],
    )


# --------------------------------------------------------------------------
# per-sample prep
# --------------------------------------------------------------------------


@dataclass
class PrepResult:
    counts: Counter = field(default_factory=Counter)
    merged: Counter = field(default_factory=Counter)
    pairs: Counter = field(default_factory=Counter)
    primer_pairs: Counter = field(default_factory=Counter)
    overlap_len: list = field(default_factory=list)
    insert_len: list = field(default_factory=list)
    # Reads with no primer are not simply waste. Primer dimer, off-target
    # amplification and adapter-only clusters all land here, and the ratio
    # between them is a library-prep diagnostic. Dropping them silently
    # would turn "the PCR failed" into "the sample was negative".
    noprimer: Counter = field(default_factory=Counter)


def _mate_id(rid: str) -> str:
    """The read id with any old-style /1 or /2 mate suffix removed.

    Current Illumina headers carry the mate number after a space, which
    open_reads already drops; pre-Casava 1.8 files put it on the id itself,
    and comparing those verbatim would call every pair out of step.
    """
    return rid[:-2] if rid.endswith(("/1", "/2")) else rid


def paired_reads(r1_path: str, r2_path: str):
    """Yield ((id, seq, qual), (id, seq, qual)) mate by mate, refusing to
    pair reads that are not mates.

    zip() would stop at the shorter file without a word. A truncated R2 -
    an interrupted copy, a full disk - then drops every pair past the
    break, and because the dropped pairs never enter the funnel the read
    accounting still balances: the one check meant to catch lost reads
    is blind to exactly this loss. So both files must run out together.
    """
    end = object()
    it1, it2 = open_reads(r1_path), open_reads(r2_path)
    n = 0
    while True:
        a, b = next(it1, end), next(it2, end)
        if a is end and b is end:
            return
        if a is end or b is end:
            short, other = (r1_path, r2_path) if a is end else (r2_path, r1_path)
            raise SystemExit(
                f"{short} ran out after {n:,} reads while {other} still has "
                "reads - one mate file is truncated or the two are not a "
                "pair. Re-fetch or regenerate the files; nothing past the "
                "break could be paired."
            )
        if _mate_id(a[0]) != _mate_id(b[0]):
            raise SystemExit(
                f"R1/R2 out of step at {a[0]} vs {b[0]} - the files are not "
                "in the same order, which no downstream step can detect"
            )
        n += 1
        yield a, b


def prep_sample(r1_path: str, r2_path: str, ps: PrimerSet, cfg: dict, args) -> PrepResult:
    res = PrepResult()
    rt_r1 = revcomp(cfg["pcr1_overhang_rev"])
    rt_r2 = revcomp(cfg["pcr1_overhang_fwd"])
    # Primer sequence by name, for the read-through junction: the primer
    # matched at the start of one mate is the one whose reverse complement
    # the other mate runs into. First entry wins on a duplicated name, as
    # the primer set itself would have no way to tell them apart either.
    fwd_seq: dict[str, str] = {}
    rev_seq: dict[str, str] = {}
    for name, p in ps.fwd:
        fwd_seq.setdefault(name, p)
    for name, p in ps.rev:
        rev_seq.setdefault(name, p)

    for (rid1, s1, q1), (rid2, s2, q2) in paired_reads(r1_path, r2_path):
        # 1. orient
        f = ps.match(s1, "fwd")
        r = ps.match(s2, "rev")
        if f is None or r is None:
            f2, r2m = ps.match(s2, "fwd"), ps.match(s1, "rev")
            if f2 is not None and r2m is not None:
                s1, q1, s2, q2 = s2, q2, s1, q1
                f, r = f2, r2m
            else:
                res.counts["no_primer"] += 1
                res.noprimer[s1[:120]] += 1
                continue

        # 2. trim primers
        s1, q1 = s1[f[1] :], q1[f[1] :]
        s2, q2 = s2[r[1] :], q2[r[1] :]
        res.primer_pairs[f"{f[0]}/{r[0]}"] += 1

        # 3. read-through: R1 runs into rc(the reverse primer R2 opened
        # with) and then the reverse adapter, R2 the mirror image.
        s1, q1 = trim_readthrough(s1, q1, rt_r1, rev_seq.get(r[0], ""),
                                  primer_edits=ps.max_edits)
        s2, q2 = trim_readthrough(s2, q2, rt_r2, fwd_seq.get(f[0], ""),
                                  primer_edits=ps.max_edits)

        if len(s1) < args.min_length or len(s2) < args.min_length:
            res.counts["too_short"] += 1
            continue
        if (
            expected_errors(q1) > args.max_ee
            or expected_errors(q2) > args.max_ee
        ):
            res.counts["low_quality"] += 1
            continue

        # 4. merge
        m = merge_pair(s1, q1, s2, q2, args.min_overlap, args.max_diff)
        if m is not None:
            seq, _q = m
            if len(seq) < args.min_length:
                res.counts["too_short"] += 1
                continue
            res.counts["merged"] += 1
            res.merged[seq] += 1
            res.insert_len.append(len(seq))
            res.overlap_len.append(len(s1) + len(s2) - len(seq))
        else:
            res.counts["paired"] += 1
            res.pairs[(s1, s2)] += 1

    return res


# --------------------------------------------------------------------------


def write_derep(path: str, counts: Counter, prefix: str) -> None:
    with gzip.open(path, "wt") as fh:
        for i, (seq, n) in enumerate(
            sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])), 1
        ):
            fh.write(f">{prefix}{i};size={n}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j : j + 60] + "\n")


def write_pairs(path: str, counts: Counter) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("uid\tsize\tr1\tr2\n")
        for i, ((a, b), n) in enumerate(
            sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])), 1
        ):
            fh.write(f"P{i}\t{n}\t{a}\t{b}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--primers", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--outdir", default="prep")
    ap.add_argument("--min-length", type=int, default=100)
    ap.add_argument(
        "--min-overlap", type=int, default=20,
        help="shortest mate overlap that is merged; shorter stays a pair",
    )
    ap.add_argument("--max-diff", type=float, default=0.25)
    ap.add_argument(
        "--max-ee", type=float, default=2.0,
        help="maximum expected errors per mate after trimming",
    )
    ap.add_argument("--primer-max-edits", type=int, default=2)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    cfg = load_assay(args.assay)
    ps = load_primers(args.primers)
    ps.max_edits = args.primer_max_edits
    if not ps.fwd or not ps.rev:
        raise SystemExit("primers: need at least one fwd and one rev primer")

    os.makedirs(args.outdir, exist_ok=True)
    res = prep_sample(args.r1, args.r2, ps, cfg, args)

    sid = args.sample_id
    write_derep(os.path.join(args.outdir, f"{sid}.merged.fasta.gz"), res.merged, "M")
    write_pairs(os.path.join(args.outdir, f"{sid}.pairs.tsv.gz"), res.pairs)
    write_derep(
        os.path.join(args.outdir, f"{sid}.noprimer.fasta.gz"), res.noprimer, "N"
    )

    total = sum(res.counts.values())
    kept = res.counts["merged"] + res.counts["paired"]
    summary = {
        "sample_id": sid,
        "pairs_in": total,
        "pairs_kept": kept,
        "fates": {k: res.counts.get(k, 0) for k in FATES},
        "fate_meaning": FATE_MEANING,
        "unique_merged": len(res.merged),
        "unique_pairs": len(res.pairs),
        "unique_noprimer": len(res.noprimer),
        "primer_pairs": dict(res.primer_pairs.most_common()),
        "median_insert": (
            sorted(res.insert_len)[len(res.insert_len) // 2]
            if res.insert_len else None
        ),
        "median_overlap": (
            sorted(res.overlap_len)[len(res.overlap_len) // 2]
            if res.overlap_len else None
        ),
        # The invariant, asserted rather than assumed: a funnel that does not
        # balance has lost reads somewhere it does not report.
        "accounting_balances": total == sum(res.counts[k] for k in FATES),
    }
    with open(os.path.join(args.outdir, f"{sid}.prep.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    if not args.quiet:
        print(f"{sid}: {total:,} pairs in")
        for k in FATES:
            n = res.counts.get(k, 0)
            if n:
                pct = 100 * n / total if total else 0
                print(f"    {k:16s} {n:8,}  {pct:5.2f}%  {FATE_MEANING[k]}")
        print(
            f"    unique: {len(res.merged):,} merged, {len(res.pairs):,} paired"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
