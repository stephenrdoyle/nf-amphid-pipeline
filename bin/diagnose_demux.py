#!/usr/bin/env python3
"""Find out why nanopore reads are failing to demultiplex.

`demux_ont.py` reports the outcome - what fraction of reads reached each
tier - and that is the right thing to report when a run is healthy. It is
the wrong thing to have when a run is not: "52% of reads resolved only one
end" tells you the run is broken and nothing about what to change.

This tool answers the next question. It replays the same end-caller over a
subsample of reads and records, for each end separately, the step at which
the end was lost:

    anchor_not_found    the construct's anchor was not within the search
                        window at the required tolerance. Either the
                        window is too short, the tolerance too tight, or
                        the sequence is not what the assay config says
    barcode_off_read    anchor found, but the barcode window runs past the
                        end of the read - the molecule is truncated right
                        where the index sits
    no_copy_matched     the barcode window is there and neither copy is
                        within the index set's tolerance
    copies_disagreed    both copies read, and they disagree - a frame slip
    ok / ok_single_copy the end yielded a barcode

Then it sweeps the two parameters that are tuning rather than physics -
the anchor search window and the anchor edit budget - so the cost of
loosening either is a measurement rather than an argument. And it scans
the whole read for the anchor, which separates "the anchor is present but
outside the window" from "the anchor is not in these reads at all". Those
two look identical in the tier counts and have nothing in common as
problems.

Read-only: it writes a report and touches no pipeline output.

    bin/diagnose_demux.py --reads <minknow dir|fastq> \
        --samplesheet samplesheet.csv --run-id RUN01 \
        --assay assets/assay/ampseq_helminth.json

Self-contained: standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    find_anchor_cost,
    hamming,
    load_assay,
    open_reads_many,
    read_table,
    resolve_read_files,
    revcomp,
)
from demux_ont import BarcodeSet, EndCall, call_end  # noqa: E402
from split_concatemers import CLOSES, OPENS  # noqa: E402

# This tool reads instrumentation that lives in demux_ont.py, so the two
# files have to be the same vintage. They are separate files and get
# copied to compute separately, which means they can and do drift - and
# an AttributeError three minutes into a run is a poor way to find out.
if not hasattr(EndCall(), "reason"):
    raise SystemExit(
        f"{os.path.join(os.path.dirname(os.path.abspath(__file__)), 'demux_ont.py')}\n"
        "is older than this tool: its EndCall has no 'reason' field, which "
        "is where\nthe per-end failure cause is recorded. Copy the whole "
        "bin/ directory across\nrather than single files - diagnose_demux.py "
        "cannot work without it."
    )

REASONS = (
    "ok",
    "ok_single_copy",
    "anchor_not_found",
    "barcode_off_read",
    "no_copy_matched",
    "copies_disagreed",
    "no_anchor_configured",
)


def pct(n: int, d: int) -> str:
    return f"{100 * n / d:5.1f}%" if d else "    -"


def both_ends(seq, cfg, i5_set, i7_set, repeats, window, edits):
    """Call both ends in the orientation that resolves more of them.

    Mirrors assign_read's orientation choice: the construct is not
    palindromic, so the i5 end is only at the start on one strand, and
    picking the wrong strand loses both ends rather than one.
    """
    best = None
    for orient, s in (("+", seq), ("-", revcomp(seq))):
        five = call_end(s, cfg["i5_inner_anchor"], i5_set, repeats, window,
                        barcode_before_anchor=True, anchor_max_edits=edits)
        tail_from = max(0, len(s) - window)
        three = call_end(s[tail_from:], cfg["i7_inner_anchor"], i7_set,
                         repeats, window, barcode_before_anchor=False,
                         anchor_max_edits=edits)
        score = (1 if five.barcode else 0) + (1 if three.barcode else 0)
        if best is None or score > best[0]:
            best = (score, orient, five, three)
    return best[1], best[2], best[3]


def scan_all(text: str, pattern: str, max_edits: int, limit: int = 12):
    """Every non-overlapping occurrence of pattern, not just the first.

    The demultiplexer looks for one anchor per end because a read is
    supposed to be one molecule. Counting them instead is how you find out
    that it is not.
    """
    hits = []
    off = 0
    while len(hits) < limit and off < len(text):
        h = find_anchor_cost(text[off:], pattern, max_edits)
        if h is None:
            break
        s, e, c = h
        hits.append((off + s, off + e, c))
        off += max(e, s + 1)
    return hits


def census(seq: str, a5: str, a7: str, edits: int) -> tuple[int, int, bool]:
    """How many complete ends of each kind this read carries, either way round.

    A read with two i5 ends and two i7 ends is two amplicons ligated into
    one molecule. Scanning both strands separately also distinguishes
    head-to-tail joins from head-to-head ones, which need different
    handling: the first splits cleanly, the second has an inverted copy.
    """
    r5 = revcomp(a5)
    r7 = revcomp(a7)
    f5, b5 = len(scan_all(seq, a5, edits)), len(scan_all(seq, r5, edits))
    f7, b7 = len(scan_all(seq, a7, edits)), len(scan_all(seq, r7, edits))
    # Mixed strands for the same anchor means an inverted repeat.
    inverted = (f5 > 0 and b5 > 0) or (f7 > 0 and b7 > 0)
    return f5 + b5, f7 + b7, inverted


def why_no_match(window: str, declared: set[str], bl: int) -> tuple[int, int, int]:
    """Is this window a damaged index, or not an index at all?

    Two measurements on the same piece of sequence, which answer
    different questions.

    The substitution distance is what the demultiplexer itself uses. The
    edit distance additionally allows indels, which is the specific thing
    the current matcher refuses - so comparing the two says whether
    indels are what is breaking it.

    And the floor matters as much as either: 8bp of unrelated sequence
    sits about 4-5 substitutions from the nearest member of a set this
    size, purely by chance. A window that is 2 away is a damaged index; a
    window that is 5 away is not an index at all, and no tolerance will
    turn it into one.

    Returns (best substitutions, best edits, how many indices tie at the
    best edit distance). The tie count is the safety question: a window
    within one indel of two different indices cannot be assigned to
    either, so "recoverable" is not the same as "recoverable safely".
    """
    if not declared or not window:
        return 9, 9, 0
    best_h = min(hamming(window[:bl], d) if len(window) >= bl else 9
                 for d in declared)
    best_e, ties = 9, 0
    for d in declared:
        # Myers, not a straight Levenshtein: the index may sit anywhere
        # in the padded window, and the padding makes the two strings
        # different lengths. A length-sensitive comparison rejects that
        # outright - bounded_edit returns its ceiling as soon as the
        # lengths differ by more than the budget, which silently turns
        # every measurement into "no match" and hides the answer.
        hit = find_anchor_cost(window, d, 3)
        if hit is None:
            continue
        e = hit[2]
        if e < best_e:
            best_e, ties = e, 1
        elif e == best_e:
            ties += 1
    return best_h, best_e, ties


def junction_gaps(seq: str, a5: str, a7: str, edits: int,
                  head_len: dict, tail_len: dict):
    """How much sequence sits between one amplicon and the next.

    Two amplicons can end up in one read by two entirely different
    routes, and the fix differs completely:

      LIGATED     joined at the bench during library prep. The two
                  molecules' outer tails abut, so the gap is ~0.
      READ-THROUGH  two separate molecules traversed the pore in quick
                  succession and were recorded as one read. Each still
                  carries its own sequencing adapter, so tens of bases
                  that belong to neither amplicon sit between them.

    One is fixed with end-prep and adapter ratios, the other by read
    splitting in the basecaller. They look identical in every summary
    statistic computed so far, and the gap distinguishes them.

    Yields (gap, sequence) for each junction found.
    """
    hits = []
    for kind, pat in (("5+", a5), ("7+", a7),
                      ("5-", revcomp(a5)), ("7-", revcomp(a7))):
        off = 0
        while off < len(seq):
            h = find_anchor_cost(seq[off:], pat, edits)
            if h is None:
                break
            hits.append((off + h[0], off + h[1], kind))
            off += max(h[1], h[0] + 1)
    hits.sort()
    kept = []
    for h in hits:
        if kept and h[0] < kept[-1][1]:
            continue
        kept.append(h)

    out = []
    for i in range(1, len(kept)):
        _p, prev_end, prev_kind = kept[i - 1]
        pos, _e, kind = kept[i]
        if prev_kind in CLOSES and kind in OPENS:
            # The only unambiguous junction: one unit finishes and the
            # next begins. Its extent is from the end of the first unit's
            # tail to the start of the second unit's head.
            a = prev_end + tail_len[prev_kind]
            b = pos - head_len[kind]
            out.append((b - a, seq[max(0, a):max(0, b)]))
    return out


def percentiles(values: list[int], points=(0, 5, 25, 50, 75, 95, 100)) -> dict:
    if not values:
        return {}
    v = sorted(values)
    out = {}
    for p in points:
        i = min(len(v) - 1, max(0, round(p / 100 * (len(v) - 1))))
        out[f"p{p}"] = v[i]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--barcode-repeats", type=int, default=2)
    ap.add_argument("--max-mismatch", type=int, default=1)
    ap.add_argument("--search-window", type=int, default=120)
    ap.add_argument("--anchor-max-edits", type=int, default=2)
    ap.add_argument("--include-failed", action="store_true")
    ap.add_argument("--n", type=int, default=20000,
                    help="reads to sample; the whole point is that this is "
                         "quick, and 20k settles every rate below to well "
                         "under a percent")
    ap.add_argument("--sweep", action="store_true", default=True)
    ap.add_argument("--no-sweep", dest="sweep", action="store_false")
    ap.add_argument("--examples", type=int, default=3,
                    help="failing reads to print in full context")
    ap.add_argument("--census", type=int, default=5000,
                    help="reads to scan end to end for repeated adapter "
                         "structures - the concatemer count. Slower per "
                         "read than the rest, so sampled separately")
    ap.add_argument("--out", help="also write the report as JSON")
    args = ap.parse_args(argv)

    cfg = load_assay(args.assay)
    samples = [s for s in read_table(args.samplesheet)
               if s["run_id"] == args.run_id]
    if not samples:
        print(f"no samples for run '{args.run_id}' in the samplesheet",
              file=sys.stderr)
        return 1
    i5_set = BarcodeSet.build({s["index_i5"].upper() for s in samples},
                              args.max_mismatch)
    i7_set = BarcodeSet.build({s["index_i7"].upper() for s in samples},
                              args.max_mismatch)

    read_files, file_report = resolve_read_files(
        args.reads, prefer_pass=not args.include_failed)
    if not read_files:
        print(f"no FASTQ found under {', '.join(args.reads)}", file=sys.stderr)
        return 1

    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    span = i5_set.length * args.barcode_repeats
    head_len = {"5+": span + len(cfg["p5_tail"]),
                "7-": span + len(cfg["p7_tail"])}
    tail_len = {"7+": span + len(cfg["p7_tail"]),
                "5-": span + len(cfg["p5_tail"])}

    # The shortest molecule that could physically carry both ends. A read
    # below this cannot resolve both no matter how the search is tuned, so
    # it is worth knowing before tuning anything.
    floor = (len(cfg["p5_tail"]) + span + len(cfg["pcr1_overhang_fwd"])
             + len(cfg["pcr1_overhang_rev"]) + span + len(cfg["p7_tail"]))

    lengths: list[int] = []
    reasons = {"i5": Counter(), "i7": Counter()}
    anchor_edits = {"i5": Counter(), "i7": Counter()}
    anchor_pos = {"i5": [], "i7": []}
    # Whole-read scan at exactly the budget the real search uses. Using a
    # looser one here would be worse than useless: a 14bp pattern at 5
    # edits matches random sequence readily, so it would report chance
    # hits as "the anchor is there, you just missed it". Same stringency,
    # different extent, is the only comparison that isolates position.
    elsewhere = {"i5": 0, "i7": 0}
    absent = {"i5": 0, "i7": 0}
    examples: list[dict] = []
    leaders: list[str] = []
    unexpected = {"i5": Counter(), "i7": Counter()}
    unexpected_noisy = {"i5": 0, "i7": 0}
    noisy_windows = {"i5": Counter(), "i7": Counter()}
    both_ends_no_index = 0
    nomatch_sub = {"i5": Counter(), "i7": Counter()}
    nomatch_edit = {"i5": Counter(), "i7": Counter()}
    indel_recoverable = {"i5": 0, "i7": 0}
    indel_ambiguous = {"i5": 0, "i7": 0}
    n5: Counter = Counter()
    n7: Counter = Counter()
    inverted = 0
    neither = 0
    one_only = 0
    gaps: list[int] = []
    gap_seqs: Counter = Counter()
    multi_reads = 0
    units = 0
    splittable = 0
    census_n = 0
    n = 0

    for read_id, seq, _qual in open_reads_many(read_files):
        if n >= args.n:
            break
        n += 1
        lengths.append(len(seq))
        orient, five, three = both_ends(
            seq, cfg, i5_set, i7_set, args.barcode_repeats,
            args.search_window, args.anchor_max_edits)
        s = seq if orient == "+" else revcomp(seq)
        tail_from = max(0, len(s) - args.search_window)

        if census_n < args.census:
            census_n += 1
            c5, c7, inv = census(seq, a5, a7, args.anchor_max_edits)
            n5[min(c5, 4)] += 1
            n7[min(c7, 4)] += 1
            inverted += inv
            units += max(c5, c7)
            # No adapter structure of either kind, anywhere on either
            # strand. Such a read cannot carry an index, so it is not a
            # demultiplexing failure - there is nothing there to read.
            if c5 == 0 and c7 == 0:
                neither += 1
            elif c5 == 0 or c7 == 0:
                one_only += 1
            # A read is multi-unit if EITHER end structure repeats. Taking
            # the larger of the two per-end totals instead would undercount
            # the union, and could report fewer multi-unit reads than
            # inverted ones - which are a subset of them.
            if c5 >= 2 or c7 >= 2:
                multi_reads += 1
            if c5 == c7 and c5 >= 2:
                splittable += 1
            if max(c5, c7) >= 2:
                for g, gs in junction_gaps(seq, a5, a7, args.anchor_max_edits,
                                           head_len, tail_len):
                    gaps.append(g)
                    if 5 <= g <= 200:
                        gap_seqs[gs[:30]] += 1

        if (five.anchor_pos >= 0 and three.anchor_pos >= 0
                and not five.barcode and not three.barcode):
            both_ends_no_index += 1

        for tag, end, anchor, sub in (("i5", five, a5, s),
                                      ("i7", three, a7, s[tail_from:])):
            reasons[tag][end.reason] += 1
            # When the anchor is found and the barcode window is in place
            # but nothing in the index set matches, the question is whether
            # the sequence there is noise or a real index this samplesheet
            # does not list. Tandem copies that agree with each other
            # answer it: noise does not repeat itself exactly.
            if end.reason == "no_copy_matched" and end.window_start >= 0:
                bl = i5_set.length if tag == "i5" else i7_set.length
                copies = [sub[end.window_start + k * bl:
                              end.window_start + (k + 1) * bl]
                          for k in range(args.barcode_repeats)]
                if all(len(c) == bl for c in copies):
                    if len(set(copies)) == 1:
                        unexpected[tag][copies[0]] += 1
                    else:
                        unexpected_noisy[tag] += 1
                        # Keep the actual sequence. If a few constant
                        # sequences dominate here, the window is sitting
                        # on something fixed - adapter, or a PCR1
                        # overhang where an index should be - which is a
                        # different diagnosis from the window landing on
                        # insert at an internal junction.
                        noisy_windows[tag][copies[0]] += 1
                # The indel question, measured on the better of the two
                # copies. A couple of bases of context either side, so an
                # insertion has somewhere to go.
                decl = i5_set.by_seq if tag == "i5" else i7_set.by_seq
                best = (9, 9, 0)
                for k in range(args.barcode_repeats):
                    lo = max(0, end.window_start + k * bl - 2)
                    w = sub[lo:end.window_start + (k + 1) * bl + 2]
                    cand = why_no_match(w, set(decl), bl)
                    if cand[1] < best[1]:
                        best = cand
                sub_d, edit_d, ties = best
                nomatch_sub[tag][min(sub_d, 6)] += 1
                nomatch_edit[tag][min(edit_d, 6)] += 1
                if edit_d <= 1 and sub_d > i5_set.max_mismatch:
                    if ties == 1:
                        indel_recoverable[tag] += 1
                    else:
                        indel_ambiguous[tag] += 1
            if end.anchor_edits >= 0:
                anchor_edits[tag][end.anchor_edits] += 1
                # Distance from the read end this anchor belongs to: the i5
                # anchor is placed from the 5' end, the i7 anchor from the
                # 3'. Reporting both as absolute offsets would make the i7
                # figure track read length instead of construct geometry.
                if tag == "i5":
                    anchor_pos[tag].append(end.anchor_pos)
                else:
                    absolute = max(0, len(s) - args.search_window) + end.anchor_pos
                    anchor_pos[tag].append(len(s) - absolute)
            if (tag == "i5" and end.barcode
                    and end.window_start >= 0 and len(leaders) < 4000):
                # The i5 anchor is preceded by the index copies and then
                # P5. Anything before THAT is not part of the construct.
                lead_end = end.window_start - len(cfg["p5_tail"])
                if lead_end > 0:
                    leaders.append(sub[:lead_end])
                elif lead_end <= 0:
                    leaders.append("")
            if end.reason == "anchor_not_found":
                hit = find_anchor_cost(s, anchor, args.anchor_max_edits)
                if hit is None:
                    absent[tag] += 1
                else:
                    elsewhere[tag] += 1
                    if len(examples) < args.examples:
                        pos = hit[0]
                        examples.append({
                            "read_id": read_id, "end": tag,
                            "read_length": len(s),
                            "found_at": pos if tag == "i5" else pos - len(s),
                            "edits": hit[2],
                            "context": s[max(0, pos - 25): pos + len(anchor) + 25],
                        })

    if not n:
        print("no reads read", file=sys.stderr)
        return 1

    # ---------------------------------------------------------------- report
    print(f"\n{'='*66}")
    print(f"demultiplexing diagnosis: run {args.run_id}")
    print(f"{'='*66}")
    print(f"{n:,} reads sampled from {file_report['n_files']} file(s)")
    print(f"settings: window {args.search_window}  anchor edits "
          f"<={args.anchor_max_edits}  repeats {args.barcode_repeats}  "
          f"index tolerance i5<={i5_set.max_mismatch} i7<={i7_set.max_mismatch}")

    q = percentiles(lengths)
    short = sum(1 for L in lengths if L < floor)
    print(f"\nread lengths  p5 {q['p5']}  p25 {q['p25']}  median {q['p50']}  "
          f"p75 {q['p75']}  p95 {q['p95']}  max {q['p100']}")
    print(f"  shortest molecule that can carry both ends: {floor}bp "
          f"(adapters and indices alone, before any insert)")
    print(f"  reads below that: {short:,} ({pct(short, n)}) - these cannot "
          "resolve both ends however the search is tuned")

    for tag, anchor in (("i5", a5), ("i7", a7)):
        r = reasons[tag]
        ok = r["ok"] + r["ok_single_copy"]
        print(f"\n{tag} end   anchor {anchor} ({len(anchor)}bp)")
        print(f"  resolved                {ok:>8,}  {pct(ok, n)}")
        if r["ok_single_copy"]:
            print(f"    of which one copy only{r['ok_single_copy']:>8,}  "
                  f"{pct(r['ok_single_copy'], n)}")
        for why in ("anchor_not_found", "barcode_off_read",
                    "no_copy_matched", "copies_disagreed"):
            if r[why]:
                print(f"  {why:<22}  {r[why]:>8,}  {pct(r[why], n)}")
        if r["anchor_not_found"]:
            print(f"    anchor absent from the read entirely "
                  f"{absent[tag]:>6,}  {pct(absent[tag], n)}")
            print(f"    anchor present but outside the window "
                  f"{elsewhere[tag]:>5,}  {pct(elsewhere[tag], n)}")
        if anchor_edits[tag]:
            tot = sum(anchor_edits[tag].values())
            dist = "  ".join(
                f"{e}:{pct(anchor_edits[tag][e], tot).strip()}"
                for e in sorted(anchor_edits[tag]))
            print(f"  anchor edit distance when found:  {dist}")
            if anchor_pos[tag]:
                p = percentiles(anchor_pos[tag], (5, 50, 95, 100))
                side = "5'" if tag == "i5" else "3'"
                expect = len(cfg["p5_tail"]) + span if tag == "i5" else (
                    len(anchor) + span + len(cfg["p7_tail"]))
                print(f"  anchor starts this far from the {side} end:  "
                      f"median {p['p50']}  p95 {p['p95']}  max {p['p100']}"
                      f"   (construct predicts {expect})")

    declared = {"i5": set(i5_set.by_seq), "i7": set(i7_set.by_seq)}
    for tag in ("i5", "i7"):
        if not unexpected[tag] and not unexpected_noisy[tag]:
            continue
        clean = sum(unexpected[tag].values())
        noisy = unexpected_noisy[tag]
        print(f"\n{tag}: barcode window read, but no index in the samplesheet "
              "matched")
        print(f"  tandem copies agree with each other  {clean:>8,}  "
              f"{pct(clean, n)}   <- a real index, not in this samplesheet")
        # Not "noise". Two disagreeing copies mean the 16bp being read is
        # not a tandem pair at all, and a window landing somewhere other
        # than an index looks exactly like this - which is what happens
        # when the anchor found was an internal junction rather than the
        # end of the molecule. Basecalling quality is measurable from the
        # anchor edit distances above; if those are mostly exact, this is
        # not a quality problem and should not be labelled as one.
        print(f"  tandem copies disagree               {noisy:>8,}  "
              f"{pct(noisy, n)}   <- the window is not on an index")
        # An index set declared in the wrong orientation is the classic
        # way to get this whole section populated, and it is invisible
        # unless you look for it: every observed index would be the
        # reverse complement of a declared one. Checking is free.
        rc_hits = sum(c for s_, c in unexpected[tag].items()
                      if revcomp(s_) in declared[tag])
        if rc_hits:
            print(f"  ** {rc_hits:,} of them are the REVERSE COMPLEMENT of a "
                  "declared index.\n     The index set for this end is "
                  "declared in the wrong orientation.")
        if noisy_windows[tag]:
            top = noisy_windows[tag].most_common(5)
            share = sum(c for _s, c in top) / max(1, sum(noisy_windows[tag].values()))
            print(f"  the 5 commonest sequences at that position cover "
                  f"{100*share:.1f}% of them:")
            for s_, c in top:
                print(f"    {s_}  {c:>7,}  {pct(c, n)}")
            print("    A few sequences dominating means the window is on "
                  "something CONSTANT\n    - an adapter, or a PCR1 overhang "
                  "where an index should be. All\n    distinct means it is "
                  "on insert, at an internal junction.")
        for seq_, c in unexpected[tag].most_common(8):
            near = min((sum(a != b for a, b in zip(seq_, d)), d)
                       for d in declared[tag]) if declared[tag] else (0, "")
            note = (f"   rc is declared: {revcomp(seq_)}"
                    if revcomp(seq_) in declared[tag]
                    else f"   nearest declared {near[1]} "
                         f"({near[0]} mismatches)")
            print(f"    {seq_}  {c:>7,}  {pct(c, n)}{note}")

    print(f"\n  reads with an anchor at BOTH ends but an index at neither: "
          f"{both_ends_no_index:,} ({100*both_ends_no_index/n:.1f}%)")
    print("    That is the shape of a product amplified by PCR1 primers at "
          "both ends:\n    the overhangs, and so the anchors, are present "
          "but no index ever\n    got added.")

    if leaders:
        print(f"\n{'-'*66}")
        print("what sits in front of the construct")
        print(f"{'-'*66}")
        lq = percentiles([len(x) for x in leaders], (5, 50, 95, 100))
        print(f"  {len(leaders):,} reads with a resolved i5 end. Sequence "
              "before P5:")
        print(f"  length  median {lq['p50']}  p95 {lq['p95']}  max {lq['p100']}"
              "   (the construct predicts 0)")
        common = Counter(x[-20:] for x in leaders if len(x) >= 8)
        if common:
            print("  most common last 20bp before P5 - if one sequence "
                  "dominates,\n  it is an untrimmed adapter and belongs in "
                  "the assay config:")
            for s_, c in common.most_common(5):
                print(f"    {s_:<22} {c:>6,}  {pct(c, len(leaders))}")

    if any(nomatch_sub[t] for t in ("i5", "i7")):
        print(f"\n{'-'*66}")
        print("is it indels? distance from the failed window to the "
              "nearest index")
        print(f"{'-'*66}")
        print("  8bp of unrelated sequence sits ~4-5 substitutions from the")
        print("  nearest member of a set this size by chance. A window 1-2 "
              "away is\n  a damaged index; one 4+ away was never an index.")
        for tag in ("i5", "i7"):
            tot = sum(nomatch_sub[tag].values())
            if not tot:
                continue
            print(f"\n  {tag}  ({tot:,} windows that matched nothing)")
            sub_row = "  ".join(
                f"{d}{'+' if d == 6 else ''}:{100*nomatch_sub[tag][d]/tot:.0f}%"
                for d in sorted(nomatch_sub[tag]))
            edit_row = "  ".join(
                f"{d}{'+' if d == 6 else ''}:{100*nomatch_edit[tag][d]/tot:.0f}%"
                for d in sorted(nomatch_edit[tag]))
            print(f"    substitutions only (what the matcher uses):  {sub_row}")
            print(f"    edits, indels allowed:                       {edit_row}")
            rec, amb = indel_recoverable[tag], indel_ambiguous[tag]
            print(f"    within 1 indel of exactly one index:   {rec:>7,}  "
                  f"{pct(rec, n)}  <- recoverable, and safely")
            print(f"    within 1 indel of more than one index: {amb:>7,}  "
                  f"{pct(amb, n)}  <- recoverable, but not safely")
        print("\n  If the edit row is shifted well left of the substitution "
              "row, indels\n  are what the matcher is refusing. If the two "
              "are alike and both sit\n  at 4+, the window is not on an "
              "index and tolerance is not the issue.")

    if gaps:
        print(f"\n{'-'*66}")
        print("what sits between one amplicon and the next")
        print(f"{'-'*66}")
        q = percentiles(gaps, (5, 25, 50, 75, 95))
        near0 = sum(1 for g in gaps if -10 <= g <= 10)
        wide = sum(1 for g in gaps if g > 10)
        print(f"  {len(gaps):,} junctions measured")
        print(f"  gap between the two units' outer tails:")
        print(f"    p5 {q['p5']}  p25 {q['p25']}  median {q['p50']}  "
              f"p75 {q['p75']}  p95 {q['p95']}")
        print(f"    within 10bp of zero  {near0:>8,}  {pct(near0, len(gaps))}"
              "   <- ligated at the bench")
        print(f"    more than 10bp       {wide:>8,}  {pct(wide, len(gaps))}"
              "   <- something between them")
        if gap_seqs:
            print("  most common sequence in a non-zero gap (first 30bp):")
            for s_, c in gap_seqs.most_common(4):
                print(f"    {s_:<32} {c:>6,}")
            print("  If one sequence dominates, it is the sequencing adapter "
                  "and these\n  are separate molecules the read splitter "
                  "missed - a basecalling\n  setting, not a bench problem.")

    if census_n:
        print(f"\n{'-'*66}")
        print("adapter structures per read - is one read one molecule?")
        print(f"{'-'*66}")
        print(f"{'ends found':>12} {'i5':>10} {'i7':>10}")
        for k in range(5):
            label = "4+" if k == 4 else str(k)
            print(f"{label:>12} {pct(n5[k], census_n):>10} "
                  f"{pct(n7[k], census_n):>10}")
        multi = multi_reads
        print(f"\n  reads carrying more than one end structure: "
              f"{multi:,} of {census_n:,} ({pct(multi, census_n).strip()})")
        print(f"  of those, i5 and i7 counts agree:           "
              f"{splittable:,}")
        print(f"  reads with an inverted repeat:              "
              f"{inverted:,} ({pct(inverted, census_n).strip()})")
        both_n = census_n - neither - one_only
        print(f"\n  reads with NO adapter of either kind   {neither:>8,}  "
              f"{pct(neither, census_n)}  <- cannot carry an index")
        print(f"  reads with one kind of end only        {one_only:>8,}  "
              f"{pct(one_only, census_n)}")
        print(f"  reads with both kinds present          {both_n:>8,}  "
              f"{pct(both_n, census_n)}")
        print(f"\n  amplicon units seen:  {units:,} in {census_n:,} reads "
              f"({units / census_n:.2f} per read)")
        print("  A read is supposed to be one molecule. Every unit beyond "
              "the first is\n  an amplicon the demultiplexer cannot see, "
              "because it only looks at the\n  two ends of the read.")

    if examples:
        print("\nreads where the anchor exists but the window missed it")
        for e in examples:
            print(f"  {e['read_id'][:24]:<24} {e['end']} at {e['found_at']:>6} "
                  f"of {e['read_length']}bp, {e['edits']} edits")
            print(f"    {e['context']}")

    sweep: list[dict] = []
    if args.sweep:
        print(f"\n{'-'*66}")
        print("parameter sweep - what loosening the search would recover")
        print(f"{'-'*66}")
        print(f"{'window':>8} {'edits':>6} {'i5 ok':>8} {'i7 ok':>8} "
              f"{'both':>8}")
        sample = []
        for i, (read_id, seq, _q) in enumerate(open_reads_many(read_files)):
            if i >= min(args.n, 5000):
                break
            sample.append(seq)
        for window in (args.search_window, 200, 400):
            for edits in (args.anchor_max_edits, args.anchor_max_edits + 2,
                          args.anchor_max_edits + 4):
                c5 = c7 = cb = 0
                for seq in sample:
                    _o, f, t = both_ends(seq, cfg, i5_set, i7_set,
                                         args.barcode_repeats, window, edits)
                    c5 += bool(f.barcode)
                    c7 += bool(t.barcode)
                    cb += bool(f.barcode and t.barcode)
                row = {"window": window, "anchor_max_edits": edits,
                       "i5": c5 / len(sample), "i7": c7 / len(sample),
                       "both": cb / len(sample)}
                sweep.append(row)
                mark = " <- current" if (
                    window == args.search_window
                    and edits == args.anchor_max_edits) else ""
                print(f"{window:>8} {edits:>6} {pct(c5, len(sample)):>8} "
                      f"{pct(c7, len(sample)):>8} {pct(cb, len(sample)):>8}"
                      f"{mark}")
        print(f"\n  n={len(sample):,} reads. 'both' is the ceiling on tier A: "
              "a read\n  needs both ends before the index pair can be "
              "looked up at all.")

    report = {
        "run_id": args.run_id,
        "reads_sampled": n,
        "settings": {
            "search_window": args.search_window,
            "anchor_max_edits": args.anchor_max_edits,
            "barcode_repeats": args.barcode_repeats,
            "i5_tolerance": i5_set.max_mismatch,
            "i7_tolerance": i7_set.max_mismatch,
        },
        "read_lengths": q,
        "adapter_floor_bp": floor,
        "reads_anchor_both_ends_no_index": both_ends_no_index,
        "reads_below_floor": short,
        "per_end": {
            tag: {
                "reasons": dict(reasons[tag]),
                "anchor_absent_from_read": absent[tag],
                "anchor_outside_window": elsewhere[tag],
                "anchor_edit_distance": dict(anchor_edits[tag]),
                "undeclared_indices": dict(unexpected[tag].most_common(20)),
                "undeclared_noisy": unexpected_noisy[tag],
                "noisy_window_top": dict(noisy_windows[tag].most_common(10)),
            }
            for tag in ("i5", "i7")
        },
        "census": {
            "reads_scanned": census_n,
            "i5_ends_per_read": dict(n5),
            "i7_ends_per_read": dict(n7),
            "inverted_repeat_reads": inverted,
            "reads_without_any_adapter": neither,
            "reads_with_one_end_kind_only": one_only,
            "amplicon_units": units,
            "units_per_read": units / census_n if census_n else 0,
            "multi_end_reads": multi_reads,
            "multi_end_reads_splittable": splittable,
        },
        "sweep": sweep,
    }
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
        print(f"\nJSON report: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
