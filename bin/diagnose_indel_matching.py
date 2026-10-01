#!/usr/bin/env python3
"""
What indel-aware barcode matching would change, measured on real reads.

The demultiplexer matches barcodes by HAMMING distance, deliberately:
the index sets are minimum edit distance 2, so one indel is not safely
correctable for the SET, and allowing indels blanket-wise would
misassign. That reasoning is right about the set and wrong about the
rule, in the same way the margin rule was: safety is a property of the
set only for a blanket tolerance. For a particular read, a window 1
indel from one declared index and 4 from every other is unambiguous,
and declining it throws away a read nothing else was competing for.

The cost of getting this wrong is not symmetric. A read declined is a
read lost; a read misassigned is a wrong answer in somebody's
prevalence estimate. So nothing here changes an assignment. It reports
what three rules WOULD do on the same reads:

  A  Hamming per copy, copies must agree          (what runs today)
  B  edit distance per copy with a margin, copies must agree
  C  edit distance summed over both copies, one decision

and it measures the one quantity that decides whether C is safe:

  RHO  how often the two tandem copies are corrupted IDENTICALLY.

C treats the copies as two observations and sums their evidence. That
is only worth more than one observation if the copies fail
independently. They are the same 8bp in the same context, so a
basecaller error systematic for that sequence corrupts both the same
way - and to the extent that happens, C is summing one observation
twice and its margin claims more confidence than it has. Simulation put
C ahead of B on volume at every rho, but its misassignment rose from
0.02% at rho=0 to 0.77% at rho=1, while B stayed at or below 0.05%
throughout. Which rule to adopt therefore depends on a number that can
only be read off the data.

Nothing is written but a report.

  diagnose_indel_matching.py --reads split/reads \\
      --samplesheet samples.csv --run-id RUN01 \\
      --assay assay.json --limit 200000
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ampseqlib import load_assay, open_reads_many, resolve_read_files  # noqa: E402
import demux_ont as dx  # noqa: E402


def seqlev(a: str, b: str) -> int:
    """Sequence-Levenshtein: Levenshtein with a free end gap on the right.

    The right measure for a barcode read IN CONTEXT. Plain Levenshtein
    assumes the barcode's boundary is known, but in a read it is not: a
    deletion inside the index pulls a base in from the adapter beyond
    it, and an insertion pushes one out. Scoring that as a plain edit
    overstates how far apart two indices really are - on these sets by
    exactly one step, which is the difference between "one indel
    correctable" and "none".

    Taking the minimum over the last row AND the last column keeps it
    symmetric, which a one-sided semi-global search is not.
    """
    prev = list(range(len(b) + 1))
    rows = [prev]
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
        rows.append(cur)
    return min(min(rows[-1]), min(r[-1] for r in rows))


def edit_match(window: str, by_seq, margin: int, cap: int) -> str | None:
    """Nearest declared index by edit distance, if it is clearly nearest.

    Same shape as the existing margin rule, different metric. The cap
    bounds how far a match may be at all; the margin is what makes it
    safe for THIS read regardless of how tightly the set is packed.
    """
    ds = sorted((seqlev(window, c), c) for c in by_seq)
    if not ds:
        return None
    best_d, best = ds[0]
    runner = ds[1][0] if len(ds) > 1 else 99
    if best_d > cap or runner - best_d < margin:
        return None
    return best


def pooled_match(windows, by_seq, margin: int, cap: int) -> str | None:
    """One decision from every copy, instead of one decision per copy."""
    if not windows:
        return None
    tot = sorted((sum(seqlev(w, c) for w in windows), c) for c in by_seq)
    best_d, best = tot[0]
    runner = tot[1][0] if len(tot) > 1 else 99
    if best_d > cap or runner - best_d < margin:
        return None
    return best


def windows_of(read: str, end, length: int, repeats: int) -> list[str]:
    """The barcode copy windows the demultiplexer actually looked at.

    Taken from the EndCall's own window_start rather than re-derived, so
    this measures the windows that failed rather than windows of this
    script's choosing. Re-deriving the position is how a diagnostic
    comes to describe a different problem from the one being diagnosed.
    """
    if end.window_start < 0:
        return []
    out = []
    for k in range(repeats):
        a = end.window_start + k * length
        if a >= 0 and a + length <= len(read):
            out.append(read[a:a + length])
    return out


def verdict(current: str | None, alt: str | None) -> str:
    if current is None and alt is None:
        return "both_declined"
    if current is None:
        return "RECOVERED"
    if alt is None:
        return "lost"
    return "same" if current == alt else "DISAGREE"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--limit", type=int, default=200000,
                    help="reads to examine. Edit distance over every "
                         "candidate for every copy is far more work than "
                         "the Hamming pass it is costing out, so this is "
                         "a sample by default rather than a full run")
    ap.add_argument("--margin", type=int, default=2)
    ap.add_argument("--cap", type=int, default=3,
                    help="furthest an edit-distance match may be, per copy")
    ap.add_argument("--pool-cap", type=int, default=6,
                    help="same, summed over the copies, for rule C")
    args = ap.parse_args(argv)

    cfg = load_assay(args.assay)
    samples, i5_set, i7_set, pair_to_sample = dx._build_sets(
        args.samplesheet, args.run_id, 1)
    repeats = int(cfg.get("barcode_repeats", 2) or 2)
    window = int(cfg.get("search_window", 200) or 200)

    files, _file_report = resolve_read_files(args.reads)
    if not files:
        raise SystemExit(f"no reads under {args.reads}")

    sets = {"i5": i5_set, "i7": i7_set}
    stats = {e: Counter() for e in ("i5", "i7")}
    rho = {e: Counter() for e in ("i5", "i7")}
    recovered_dist = {e: Counter() for e in ("i5", "i7")}
    seen = 0

    for _rid, seq, _q in open_reads_many(files):
        seen += 1
        if seen > args.limit:
            break
        # Both orientations, keeping whichever resolves more ends - the
        # same choice assign_read makes. Scanning only the forward
        # strand looks like it halves the work; what it actually does is
        # record every reverse-oriented read as having no barcode
        # window, which is the majority of the input on this platform
        # and would make the diagnostic describe orientation rather than
        # matching.
        best = None
        for s in (seq, dx.revcomp(seq)):
            five = dx.call_end(s, cfg["i5_inner_anchor"], i5_set, repeats,
                               window, barcode_before_anchor=True)
            tail_from = max(0, len(s) - window)
            three = dx.call_end(s[tail_from:], cfg["i7_inner_anchor"],
                                i7_set, repeats, window,
                                barcode_before_anchor=False)
            score = (1 if five.window_start >= 0 else 0) + \
                    (1 if three.window_start >= 0 else 0)
            if best is None or score > best[0]:
                best = (score, five, three, s, tail_from)
        _score, five, three, s, tail_from = best

        for end_name, end, text in (("i5", five, s),
                                ("i7", three, s[tail_from:])):
            bc = sets[end_name]
            if not bc.length:
                continue
            ws = windows_of(text, end, bc.length, repeats)
            if not ws:
                # No anchor, or the barcode ran off the read. Neither
                # rule can touch those, and counting them would
                # dilute the comparison with reads no matching
                # change can reach.
                stats[end_name]["no_window"] += 1
                continue

            # RHO. Among ends where the copies are not both a clean
            # exact match, how often are the two windows byte
            # identical - which is what systematic miscalling of the
            # same 8bp in the same context produces.
            if len(ws) >= 2:
                exact = all(w in bc.by_seq for w in ws)
                if not exact:
                    rho[end_name]["imperfect"] += 1
                    if len(set(ws)) == 1:
                        rho[end_name]["identical"] += 1

            b_calls = [edit_match(w, bc.by_seq, args.margin, args.cap)
                       for w in ws]
            b_seen = [c for c in b_calls if c]
            b = b_seen[0] if b_seen and len(set(b_seen)) == 1 else None
            c = pooled_match(ws, bc.by_seq, args.margin, args.pool_cap)

            stats[end_name]["B:" + verdict(end.barcode, b)] += 1
            stats[end_name]["C:" + verdict(end.barcode, c)] += 1
            if end.barcode is None and b is not None:
                recovered_dist[end_name][
                    min(seqlev(w, b) for w in ws)] += 1

    print(f"\n{min(seen, args.limit):,} reads examined from "
          f"{len(files)} file(s)\n")
    print("  A = Hamming per copy, copies agree   (what runs today)")
    print("  B = edit distance per copy + margin, copies agree")
    print("  C = edit distance summed over copies, one decision")
    print(f"      margin {args.margin}, cap {args.cap} per copy / "
          f"{args.pool_cap} pooled\n")

    for end in ("i5", "i7"):
        st = stats[end]
        usable = sum(v for k, v in st.items() if k.startswith("B:"))
        if not usable:
            continue
        print(f"--- {end}   {usable:,} ends with a readable barcode window "
              f"({st['no_window']:,} had none)")
        for rule in ("B", "C"):
            rec = st[f"{rule}:RECOVERED"]
            lost = st[f"{rule}:lost"]
            dis = st[f"{rule}:DISAGREE"]
            print(f"    {rule}:  recovered {rec:>9,} ({100*rec/usable:5.2f}%)"
                  f"   lost {lost:>7,}   DISAGREE {dis:>7,}"
                  f" ({100*dis/usable:.3f}%)")
        if recovered_dist[end]:
            d = recovered_dist[end]
            spread = ", ".join(f"d={k}: {d[k]:,}" for k in sorted(d))
            print(f"    B recovered, by edit distance - {spread}")
        r = rho[end]
        if r["imperfect"]:
            frac = r["identical"] / r["imperfect"]
            print(f"    RHO {frac:.3f}  ({r['identical']:,} of "
                  f"{r['imperfect']:,} imperfect ends had both copies "
                  f"byte-identical)")
            print("         low  -> the copies fail independently, so C's "
                  "pooling is real evidence")
            print("         high -> they fail together, so C sums one "
                  "observation twice and its margin overstates confidence")
        print()

    print("DISAGREE is the number that matters. A rule that recovers reads "
          "is only\nworth having if it does not also move reads that were "
          "already assigned.\nNothing was changed by this run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
