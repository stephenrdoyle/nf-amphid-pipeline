#!/usr/bin/env python3
"""
Why declared indices fail to match, and what would recover them.

On the reference MinION run, 10.7M ends (19.2% of i5, 9.9% of i7) had
their barcode window read and matched no declared index. The index is
physically there. Something about the matching rejected it, and the
tiers cannot say what.

Two mechanisms explain it, and they have opposite conclusions:

  MISPOSITIONED   the window is in the wrong place. Both tandem copies
                  are located by a fixed offset from the anchor with
                  +/- 3 slack, so a drift larger than that moves both
                  together and both fail. For the i5 end the offset is
                  derived as (alignment end - anchor length), which is
                  wrong by the net indel count INSIDE the anchor - so
                  an indel there displaces the window systematically.
                  Fixable in software.

  MISCALLED       the window is in the right place and the bases are
                  wrong. Both copies are the same 8bp sequence in the
                  same context, so a basecaller error that is
                  systematic for that sequence corrupts both
                  identically. Not fixable by matching.

The two are distinguishable: under MISPOSITIONED the correct index sits
nearby at a shifted offset, under MISCALLED it does not sit anywhere.
This measures which, and then costs each candidate remedy in ends
recovered against ambiguity introduced.

Nothing here changes any output. It reads reads and prints.

  diagnose_barcodes.py --reads split/reads --samplesheet samples.tsv \\
      --run-id RUN01 --assay assay.json --limit 200000
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ampseqlib import (  # noqa: E402
    find_anchor_cost,
    load_assay,
    open_reads_many,
    resolve_read_files,
)
import demux_ont as dx  # noqa: E402


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def set_geometry(name: str, seqs: list[str]) -> dict:
    """The distances that decide what tolerance is safe.

    k substitutions are correctable only if every pair is more than 2k
    apart in Hamming distance; k indels likewise in edit distance. A
    tolerance beyond that does not lose reads, it misassigns them.
    """
    if len(seqs) < 2:
        # One index imposes no constraint at all, so there is no safe
        # tolerance to report. Saying "49 substitutions" because the
        # minimum over an empty set of pairs is infinity would be worse
        # than saying nothing.
        return {"n": len(seqs), "min_hamming": None, "closest_hamming": ("", ""),
                "min_edit": None, "closest_edit": ("", ""),
                "safe_sub": None, "safe_indel": None}
    min_h, min_e = 99, 99
    pair_h = pair_e = ("", "")
    for i, a in enumerate(sorted(seqs)):
        for b in sorted(seqs)[i + 1:]:
            h = hamming(a, b)
            if h < min_h:
                min_h, pair_h = h, (a, b)
            hit = find_anchor_cost(a, b, len(b))
            e = hit[2] if hit else len(b)
            if e < min_e:
                min_e, pair_e = e, (a, b)
    return {"n": len(seqs), "min_hamming": min_h, "closest_hamming": pair_h,
            "min_edit": min_e, "closest_edit": pair_e,
            "safe_sub": max(0, (min_h - 1) // 2),
            "safe_indel": max(0, (min_e - 1) // 2)}


def landmark_gap(read: str, outer: str, anchor: str, before: bool,
                 max_edits: int) -> int | None:
    """Net indel count in the index region, measured not guessed.

    The index does not float freely: it sits between two pieces of
    FIXED, known sequence.

        i5 end:  [p5_tail 29][i5 x n][i5 anchor 14]
        i7 end:  [i7 anchor 15][i7 x n][p7_tail 24]

    So the distance between the inner edge of the outer tail and the
    outer edge of the anchor should be exactly the index span. It is
    not a parameter to be swept - it is a length that can be read off
    the molecule. Measuring 14 where 16 is expected means two net
    deletions in the index; 18 means two insertions.

    That does two things a wider search cannot. It places AND sizes the
    window, so a 7bp index from a deletion can be compared as 7bp
    rather than charged a mismatch at every base after the deletion.
    And being a measurement it adds no ambiguity, where a looser
    tolerance trades lost reads for misassigned ones.

    Returns the SIGNED difference from the expected span, or None when
    either landmark is missing - which is the common case at a
    truncated 3' end, where the outer tail is the first thing lost.
    """
    # Tolerance scales with the landmark's LENGTH. The outer tail is
    # 24-29bp, and at nanopore error rates a fixed budget of 2 edits
    # over 29bp finds it only in the cleanest reads - which biases the
    # measurement towards exactly the reads that did not fail. The
    # first run of this measured 1,207 of 33,978 failed ends and would
    # have been read as if it described all of them.
    a = find_anchor_cost(read, outer, max(max_edits, len(outer) // 5))
    b = find_anchor_cost(read, anchor, max(max_edits, len(anchor) // 5))
    if a is None or b is None:
        return None
    if before:
        # outer tail, then index, then anchor
        gap = b[0] - a[1]
    else:
        # anchor, then index, then outer tail
        gap = a[0] - b[1]
    return gap


def probe_window(read: str, expect: int, length: int, declared: list[str],
                 wide: int) -> dict:
    """Everything worth knowing about one failed barcode position.

    Returns the best Hamming match within the CURRENT slack, the best
    anywhere within +/- wide, and the best allowing indels. The offset
    of the wide best is the discriminator: 0 means the window was
    already right and the bases are wrong; non-zero means it was in the
    wrong place.
    """
    out = {"best": 99, "runner_up": 99, "best_offset": None,
           "wide_best": 99, "wide_offset": None, "indel_best": 99,
           "indel_runner_up": 99}
    for shift in range(-wide, wide + 1):
        a = expect + shift
        if a < 0 or a + length > len(read):
            continue
        w = read[a:a + length]
        ds = sorted(hamming(w, c) for c in declared)
        if not ds:
            continue
        if ds[0] < out["wide_best"]:
            out["wide_best"], out["wide_offset"] = ds[0], shift
        if abs(shift) <= 3 and ds[0] < out["best"]:
            out["best"] = ds[0]
            out["runner_up"] = ds[1] if len(ds) > 1 else 99
            out["best_offset"] = shift
    # Indels, at the EXPECTED position only. The pad is +/- 1, just
    # enough for the index to be one base shorter or longer. A wider
    # pad measures the wrong thing: it finds a correctly-placed copy
    # that the offset search missed and reports it as an indel, so a
    # mispositioned end is credited to indel tolerance. That happened
    # on the validation set - 400 shifted reads all came back "indel
    # distance 0" - which is why the window is this tight.
    pad = read[max(0, expect - 1):expect + length + 1]
    if pad:
        costs = []
        for c in declared:
            hit = find_anchor_cost(pad, c, 3)
            costs.append(hit[2] if hit else 99)
        costs.sort()
        out["indel_best"] = costs[0]
        # The runner-up matters more than the best. The index sets here
        # are only 2 edits apart, so (min_edit - 1) // 2 = 0 indels are
        # safely correctable by tolerance alone. A MARGIN rule is safe
        # whatever the set geometry: accept the nearest only when the
        # next nearest is clearly further, and decline otherwise. This
        # is the number that says how much of the indel-recoverable
        # pool can be taken without inventing misassignments.
        out["indel_runner_up"] = costs[1] if len(costs) > 1 else 99
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--barcode-repeats", type=int, default=2)
    ap.add_argument("--search-window", type=int, default=200)
    ap.add_argument("--anchor-max-edits", type=int, default=2)
    ap.add_argument("--wide", type=int, default=12,
                    help="how far either side of the expected offset to "
                         "look for the index the matcher missed")
    ap.add_argument("--limit", type=int, default=200000,
                    help="reads to examine. A subsample answers this; the "
                         "whole run is hours and would say the same")
    args = ap.parse_args(argv)

    cfg = load_assay(args.assay)
    _, i5_set, i7_set, pair_to_sample = dx._build_sets(
        args.samplesheet, args.run_id, 1)

    print("index set geometry")
    for name, bs in (("i5", i5_set), ("i7", i7_set)):
        g = set_geometry(name, list(bs.by_seq))
        if g["min_hamming"] is None:
            print(f"  {name}: n={g['n']} - a single index constrains "
                  "nothing, so no safe tolerance can be derived")
            continue
        print(f"  {name}: n={g['n']}  min Hamming {g['min_hamming']} "
              f"({g['closest_hamming'][0]} vs {g['closest_hamming'][1]})  "
              f"min edit {g['min_edit']}")
        print(f"      safely correctable: {g['safe_sub']} substitution(s), "
              f"{g['safe_indel']} indel(s)   [in use: "
              f"{bs.max_mismatch} substitution(s), 0 indels]")

    ends = {"i5": Counter(), "i7": Counter()}
    dist = {"i5": Counter(), "i7": Counter()}
    offs = {"i5": Counter(), "i7": Counter()}
    indel = {"i5": Counter(), "i7": Counter()}
    margin2 = {"i5": Counter(), "i7": Counter()}
    imargin = {"i5": Counter(), "i7": Counter()}
    # The gap, split by whether the end matched. The comparison is the
    # point: if failing ends have a different gap distribution from
    # matching ones, the index region really is being displaced.
    gap_ok = {"i5": Counter(), "i7": Counter()}
    gap_bad = {"i5": Counter(), "i7": Counter()}
    seen = 0

    read_files, _ = resolve_read_files(args.reads)
    if not read_files:
        raise SystemExit(f"no FASTQ found under {', '.join(args.reads)}")
    for read_id, seq, qual in open_reads_many(read_files):
        if seen >= args.limit:
            break
        seen += 1
        call = dx.assign_read(read_id, seq, cfg, i5_set, i7_set,
                              pair_to_sample, args.barcode_repeats,
                              args.search_window, False)
        s = seq if call.orientation == "+" else dx.revcomp(seq)

        for end, reason, anchor, bset, before in (
            ("i5", call.i5_reason, cfg["i5_inner_anchor"], i5_set, True),
            ("i7", call.i7_reason, cfg["i7_inner_anchor"], i7_set, False),
        ):
            ends[end][reason] += 1

            # Measured for EVERY end, matched or not, because the
            # failing distribution means nothing without the matching
            # one to compare it against.
            outer = cfg["p5_tail"] if before else dx.revcomp(cfg["p7_tail"])
            span = bset.length * args.barcode_repeats
            g = landmark_gap(s, outer, anchor, before, args.anchor_max_edits)
            if g is not None:
                bucket = gap_ok if reason.startswith("ok") else gap_bad
                bucket[end][max(-9, min(9, g - span))] += 1

            if reason != "no_copy_matched":
                continue

            # Re-locate the anchor exactly as the demultiplexer does, so
            # the window examined is the window that was rejected.
            if before:
                w = s[:args.search_window]
                base = 0
            else:
                base = max(0, len(s) - args.search_window)
                w = s[base:]
            exact = w.find(anchor)
            if exact >= 0:
                a_start, a_end = exact, exact + len(anchor)
            else:
                hit = find_anchor_cost(w, anchor, args.anchor_max_edits)
                if hit is None:
                    continue
                a_start, a_end = hit[0], hit[1]

            span = bset.length * args.barcode_repeats
            start = (a_start - span) if before else a_end
            declared = list(bset.by_seq)
            # Copy 0 only: the copies share an offset, so their fate is
            # the same and counting both would double every tally.
            p = probe_window(w, start, bset.length, declared, args.wide)
            dist[end][min(p["best"], 9)] += 1
            indel[end][min(p["indel_best"], 9)] += 1
            if p["wide_offset"] is not None and p["wide_best"] <= 1:
                offs[end][p["wide_offset"]] += 1
            if p["indel_best"] <= 2:
                imargin[end]["clear" if
                             p["indel_runner_up"] - p["indel_best"] >= 2
                             else "ambiguous"] += 1
            if p["best"] == 2:
                margin2[end]["clear" if p["runner_up"] - p["best"] >= 2
                             else "ambiguous"] += 1

    print(f"\n{seen:,} reads examined")
    for end in ("i5", "i7"):
        n = ends[end]["no_copy_matched"]
        tot = sum(ends[end].values())
        print(f"\n{end}: {n:,} of {tot:,} ends failed as no_copy_matched "
              f"({100*n/max(1,tot):.2f}%)")
        if not n:
            continue
        print("  Hamming distance to the nearest declared index, at the "
              "offset now used:")
        for d in sorted(dist[end]):
            c = dist[end][d]
            lab = f"{d}" if d < 9 else ">=9"
            print(f"    {lab:>4}  {c:>8,}  {100*c/n:5.1f}%")
        print("  with indels allowed (length-flexible window):")
        for d in sorted(indel[end]):
            c = indel[end][d]
            lab = f"{d}" if d < 9 else ">=9"
            print(f"    {lab:>4}  {c:>8,}  {100*c/n:5.1f}%")
        if offs[end]:
            print(f"  where a near-exact match (<=1) actually sits, as an "
                  f"offset from the expected position:")
            for o in sorted(offs[end]):
                c = offs[end][o]
                flag = "" if abs(o) <= 3 else "   <- OUTSIDE current slack"
                print(f"    {o:+3d}  {c:>8,}  {100*c/n:5.1f}%{flag}")
        if margin2[end]:
            clear = margin2[end]["clear"]
            amb = margin2[end]["ambiguous"]
            print(f"  of those at distance 2: {clear:,} have a clear "
                  f"runner-up gap (>=2), {amb:,} are ambiguous")

    print("\nnet indel count in the index region, from the distance "
          "between\nthe outer tail and the inner anchor (0 = nominal length)")
    for end in ("i5", "i7"):
        ok, bad = gap_ok[end], gap_bad[end]
        n_ok, n_bad = sum(ok.values()), sum(bad.values())
        if not (n_ok or n_bad):
            print(f"  {end}: neither landmark located - nothing to measure")
            continue
        print(f"  {end}:  matched ends n={n_ok:,}   failed ends n={n_bad:,}")
        print(f"    {'delta':>6}  {'matched':>16}  {'failed':>16}")
        for d in sorted(set(ok) | set(bad)):
            lab = f"{d:+d}" if abs(d) < 9 else ("<=-9" if d < 0 else ">=+9")
            a = ok.get(d, 0)
            b = bad.get(d, 0)
            print(f"    {lab:>6}  {a:>9,} {100*a/max(1,n_ok):5.1f}%  "
                  f"{b:>9,} {100*b/max(1,n_bad):5.1f}%")
        # The asymmetry is the answer to "are insertions and deletions
        # equally frequent". For this run, on this basecall model.
        for lab, c, n in (("matched", ok, n_ok), ("failed", bad, n_bad)):
            if not n:
                continue
            dele = sum(v for d, v in c.items() if d < 0)
            ins = sum(v for d, v in c.items() if d > 0)
            exact = c.get(0, 0)
            print(f"    {lab}: {100*exact/n:.1f}% nominal, "
                  f"{100*dele/n:.1f}% net deletion, "
                  f"{100*ins/n:.1f}% net insertion"
                  + (f"  (del:ins = {dele/ins:.1f}:1)" if ins else ""))

    print("\nwhat each remedy would recover, of the ends above")
    for end in ("i5", "i7"):
        n = ends[end]["no_copy_matched"]
        if not n:
            continue
        d = dist[end]
        i = indel[end]
        wide_out = sum(c for o, c in offs[end].items() if abs(o) > 3)
        print(f"  {end}:")
        print(f"    tolerance 2 substitutions      {d.get(2,0):>8,}  "
              f"{100*d.get(2,0)/n:5.1f}%  (needs min Hamming >= 5)")
        print(f"    distance 2 with a margin rule  "
              f"{margin2[end]['clear']:>8,}  "
              f"{100*margin2[end]['clear']/n:5.1f}%  (this window only - "
              "with the +/-3bp slack it cost more ends than it saved on "
              "simulated plates; see BarcodeSet.match)")
        print(f"    1-2 indels, blanket tolerance  {i.get(1,0)+i.get(2,0):>8,}  "
              f"{100*(i.get(1,0)+i.get(2,0))/n:5.1f}%  UNSAFE on this set "
              f"(min edit 2 allows 0)")
        print(f"    1-2 indels WITH a margin rule  "
              f"{imargin[end]['clear']:>8,}  "
              f"{100*imargin[end]['clear']/n:5.1f}%  safe: nearest beats "
              f"runner-up by >=2")
        print(f"      (of which declined as ambiguous: "
              f"{imargin[end]['ambiguous']:,})")
        print(f"    wider offset search            {wide_out:>8,}  "
              f"{100*wide_out/n:5.1f}%  (mispositioned, not miscalled)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
