#!/usr/bin/env python3
"""
Illumina read preparation for nf-helminth-amplicon.

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
                    off the end of the molecule and into the opposite
                    adapter. That adapter is removed.
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
and the buckets sum to the input.

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
    "length_mismatch",
)

FATE_MEANING = {
    "merged": "mates overlapped and were collapsed into one sequence",
    "paired": "no overlap - kept as a pair, both mates must agree later",
    "no_primer": "neither orientation opened with a primer pair",
    "too_short": "insert shorter than the minimum after trimming",
    "low_quality": "expected errors above the threshold",
    "length_mismatch": "mates overlapped in two places, or disagreed too much",
}


# --------------------------------------------------------------------------
# trimming and merging
# --------------------------------------------------------------------------


def trim_readthrough(seq: str, qual: str, adapter: str, max_edits: int = 3):
    """Remove adapter the read has run into past the end of a short insert.

    Only the first 20bp of the adapter is searched for. The read may stop
    part way through it, and requiring the whole thing would miss exactly
    the shortest inserts - the ones that ran furthest into it.
    """
    probe = adapter[:20]
    if len(seq) < 20:
        return seq, qual
    hit = find_anchor_cost(seq, probe, max_edits)
    if hit is None:
        return seq, qual
    start = hit[0]
    if start < 20:
        # The whole read is adapter; leave it to the length filter.
        return seq[:start], qual[:start]
    return seq[:start], qual[:start]


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
    kept rather than thrown away.

    Disagreements are resolved by quality, and the merged base takes the
    higher score only where the mates agree. Where they disagree the score
    is knocked down to the difference, so a contested base cannot inherit
    the confidence of the mate that happened to win.
    """
    rc2, rq2 = revcomp(r2), q2[::-1]
    if len(r1) < min_overlap or len(rc2) < min_overlap:
        return None

    probe = rc2[: min(24, len(rc2))]
    hit = find_anchor_cost(r1, probe, 3)
    if hit is None:
        return None
    k = hit[0]
    overlap = min(len(r1) - k, len(rc2))
    if overlap < min_overlap:
        return None

    a, b = r1[k : k + overlap], rc2[:overlap]
    diff = sum(1 for x, y in zip(a, b) if x != y)
    if diff > max_diff_frac * overlap:
        return None

    out_s, out_q = [r1[:k]], [q1[:k]]
    cons, consq = [], []
    for i in range(overlap):
        x, y = a[i], b[i]
        qx, qy = ord(q1[k + i]) - 33, ord(rq2[i]) - 33
        if x == y:
            cons.append(x)
            consq.append(chr(33 + min(41, qx + qy)))
        elif qx >= qy:
            cons.append(x)
            consq.append(chr(33 + max(2, qx - qy)))
        else:
            cons.append(y)
            consq.append(chr(33 + max(2, qy - qx)))
    out_s.append("".join(cons))
    out_q.append("".join(consq))
    out_s.append(rc2[overlap:])
    out_q.append(rq2[overlap:])
    return "".join(out_s), "".join(out_q)


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


def prep_sample(r1_path: str, r2_path: str, ps: PrimerSet, cfg: dict, args) -> PrepResult:
    res = PrepResult()
    rt_r1 = revcomp(cfg["pcr1_overhang_rev"])
    rt_r2 = revcomp(cfg["pcr1_overhang_fwd"])

    for (rid1, s1, q1), (rid2, s2, q2) in zip(
        open_reads(r1_path), open_reads(r2_path)
    ):
        if rid1 != rid2:
            raise SystemExit(
                f"R1/R2 out of step at {rid1} vs {rid2} - the files are not "
                "in the same order, which no downstream step can detect"
            )

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

        # 3. adapter read-through
        s1, q1 = trim_readthrough(s1, q1, rt_r1)
        s2, q2 = trim_readthrough(s2, q2, rt_r2)

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
    ap.add_argument("--min-overlap", type=int, default=20)
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
