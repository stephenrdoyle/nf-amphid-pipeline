#!/usr/bin/env python3
"""
Cross-run and paired-platform comparison for nf-helminth-amplicon.

The same sample can appear twice: as a replicate on another plate, or
sequenced on both Illumina and nanopore from the same library. Either way
the second measurement is the only external check the assay has. Everything
else in this pipeline is the assay agreeing with itself.

Three things are measured, and the third is the one that matters.

  agreement   How often the two runs make the same call. Raw agreement
              flatters any assay whose answer is usually "not detected" -
              two coins that both say no ninety percent of the time agree
              eighty-two percent of the time - so Cohen's kappa is reported
              beside it, which asks how much of the agreement survives
              after chance is removed.

  copies      Where both runs detected, the ratio between their copy
              numbers. Reported as a spread, not a correlation: a
              correlation across four logs of burden will look excellent
              while individual samples are out by tenfold, because the
              range does the work rather than the agreement.

  asymmetry   Which run found it when only one did, and how much of that is
              explained by depth. A disagreement between a deep run and a
              shallow one is not evidence that the platforms differ.

Pairing is by sample_id. Replicates usually carry a suffix, so --strip
takes a regular expression whose match is removed to recover the subject.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table  # noqa: E402


def subject_of(sample_id: str, strip: str | None) -> str:
    if not strip:
        return sample_id
    return re.sub(strip, "", sample_id)


def load_calls(paths: list[str]) -> list[dict]:
    rows = []
    for p in paths:
        for r in read_table(p):
            r["_source"] = p
            rows.append(r)
    return rows


def cohen_kappa(pairs: list[tuple[str, str]]) -> float | None:
    """Agreement beyond what the marginal rates would produce by chance.

    Undefined when both raters always say the same thing - there is no
    chance-corrected agreement to measure when there is no disagreement
    possible. Returning None says so instead of returning 1.0, which would
    read as perfect agreement rather than as an unanswerable question.
    """
    if not pairs:
        return None
    labels = sorted({x for p in pairs for x in p})
    if len(labels) < 2:
        return None
    n = len(pairs)
    observed = sum(1 for a, b in pairs if a == b) / n
    ma = Counter(a for a, _ in pairs)
    mb = Counter(b for _, b in pairs)
    expected = sum((ma[k] / n) * (mb[k] / n) for k in labels)
    if abs(1 - expected) < 1e-12:
        return None
    return (observed - expected) / (1 - expected)


def quantiles(xs: list[float], ps=(0.1, 0.5, 0.9)) -> dict:
    if not xs:
        return {}
    xs = sorted(xs)
    return {f"p{int(p*100)}": xs[min(len(xs) - 1, int(p * (len(xs) - 1)))]
            for p in ps}


def compare(
    rows_a: list[dict], rows_b: list[dict], strip: str | None,
    label_a: str, label_b: str,
) -> dict:
    def index(rows):
        out = {}
        for r in rows:
            key = (subject_of(r["sample_id"], strip), r["seq_id"])
            # A subject measured twice within one side (two replicates on
            # the same run) is not a pairing - it is ambiguous, and taking
            # whichever came first would hide that.
            out.setdefault(key, []).append(r)
        return out

    ia, ib = index(rows_a), index(rows_b)
    shared = sorted(set(ia) & set(ib))
    ambiguous = [k for k in shared if len(ia[k]) > 1 or len(ib[k]) > 1]
    shared = [k for k in shared if k not in set(ambiguous)]

    all_pairs: list[tuple[str, str]] = []
    call_pairs: list[tuple[str, str]] = []
    n_failed = 0
    ratios: list[float] = []
    only: Counter = Counter()
    # (side that detected, depth A, depth B) for each one-sided detection
    disagreement_depth: list[tuple[str, int, int]] = []
    per_taxon: dict = defaultdict(lambda: {"n": 0, "agree": 0})

    for key in shared:
        a, b = ia[key][0], ib[key][0]
        ca, cb = a["call"], b["call"]
        all_pairs.append((ca, cb))
        # assay_failed is not a call - it is the absence of one. Two
        # failed wells "agreeing" says nothing about whether the runs see
        # the same biology, and a failed well against a negative is not a
        # disagreement about it either. Counting them made a run that
        # failed on both sides look concordant, and fed a third label into
        # kappa whose marginals were set by bench failures rather than by
        # detection. They stay in the call matrix, where they are visible,
        # and the number left out is reported.
        if "assay_failed" in (ca, cb):
            n_failed += 1
            continue
        call_pairs.append((ca, cb))
        t = a.get("taxon") or key[1]
        per_taxon[t]["n"] += 1
        per_taxon[t]["agree"] += int(ca == cb)

        if ca == "detected" and cb == "detected":
            va, vb = a.get("copies_per_unit_dna"), b.get("copies_per_unit_dna")
            if va and vb and float(va) > 0 and float(vb) > 0:
                ratios.append(math.log10(float(vb) / float(va)))
        elif "detected" in (ca, cb):
            finder = label_a if ca == "detected" else label_b
            only[finder] += 1
            try:
                disagreement_depth.append(
                    (finder, int(a["depth"]), int(b["depth"])))
            except (KeyError, ValueError):
                pass

    n = len(call_pairs)
    agree = sum(1 for x, y in call_pairs if x == y)
    kappa = cohen_kappa(call_pairs)

    # Is the side that found it simply the deeper side? If so the platforms
    # are not disagreeing about biology.
    #
    # Asked per disagreement: was the side that DETECTED the deeper one?
    # Counting how often A was deeper, regardless of who detected, answers
    # a different question - A could be deeper every time and find nothing
    # B missed. Equal depths are counted on their own rather than handed
    # to either side, because a tie is evidence for neither reading.
    deeper_found = None
    if disagreement_depth:
        by_deeper = by_shallower = ties = 0
        for finder, da, db in disagreement_depth:
            mine, theirs = (da, db) if finder == label_a else (db, da)
            if mine > theirs:
                by_deeper += 1
            elif mine < theirs:
                by_shallower += 1
            else:
                ties += 1
        deeper_found = {
            "n_disagreements_with_depth": len(disagreement_depth),
            "found_by_deeper_run": by_deeper,
            "found_by_shallower_run": by_shallower,
            "equal_depth": ties,
            # Only true when it is true: a strict majority of all the
            # one-sided detections with depth, ties included in the
            # denominator so that they cannot manufacture a majority.
            "mostly_found_by_deeper_run": (
                by_deeper * 2 > len(disagreement_depth)),
        }

    return {
        "label_a": label_a,
        "label_b": label_b,
        # Every pairing found, failed wells included - so "nothing
        # paired" still means the ids did not match, not that the wells
        # failed.
        "paired_observations": len(all_pairs),
        # The denominator of agreement, kappa and per-taxon agreement.
        "pairs_compared": n,
        "excluded_assay_failed": n_failed,
        "ambiguous_pairings": len(ambiguous),
        "agreement": round(agree / n, 4) if n else None,
        "cohen_kappa": None if kappa is None else round(kappa, 4),
        "kappa_note": (
            "Raw agreement is inflated by how often both runs say "
            "not_detected. Kappa removes that; a value near zero means the "
            "agreement is no better than the marginal rates predict. Both "
            "leave out pairs where either well was assay_failed, which is "
            "no call at all."
        ),
        "call_matrix": {f"{x} / {y}": c
                        for (x, y), c in Counter(all_pairs).most_common()},
        "detected_by_one_only": dict(only),
        "depth_of_disagreements": deeper_found,
        "log10_copy_ratio": quantiles(ratios),
        "copy_ratio_note": (
            "log10 of B/A where both detected. Zero is agreement; 0.3 is a "
            "twofold difference. Reported as a spread rather than a "
            "correlation, which a wide dynamic range would flatter."
        ),
        "n_both_detected": len(ratios),
        "within_twofold": (
            round(sum(1 for x in ratios if abs(x) < math.log10(2))
                  / len(ratios), 4) if ratios else None
        ),
        "per_taxon_agreement": {
            k: round(v["agree"] / v["n"], 4)
            for k, v in sorted(per_taxon.items()) if v["n"]
        },
    }


def side_labels(rows_a: list[dict], rows_b: list[dict], label_a: str,
                label_b: str) -> tuple[str, str]:
    """Names for the two sides: as given, else each file's own run_id.

    They must differ. Every per-side key - call_<label>, detected_by_one_only[<label>] -
    is built from them, so two sides sharing a name (a run compared with
    a re-run of itself, both RUN01) would silently overwrite each other.
    """
    def own(rows: list[dict], fallback: str) -> str:
        ids = {r.get("run_id", "") for r in rows} - {""}
        return ids.pop() if len(ids) == 1 else fallback

    a = label_a or own(rows_a, "A")
    b = label_b or own(rows_b, "B")
    if a == b:
        a, b = f"{a}_A", f"{b}_B"
    return a, b


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--a", nargs="+", required=True,
                    help="calls.tsv for the first run (or runs)")
    ap.add_argument("--b", nargs="+", required=True,
                    help="calls.tsv for the second")
    ap.add_argument("--label-a", default="",
                    help="name for the first side; defaults to the run_id "
                         "its calls.tsv records, else 'A'")
    ap.add_argument("--label-b", default="",
                    help="name for the second side; defaults as --label-a, "
                         "else 'B'")
    ap.add_argument("--strip", default="",
                    help="regex removed from sample_id to recover the "
                         "subject, e.g. '-r[0-9]+$' for replicates")
    ap.add_argument("--outdir", default="comparison")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    rows_a = load_calls(args.a)
    rows_b = load_calls(args.b)
    if not rows_a or not rows_b:
        raise SystemExit("one side has no calls to compare")
    args.label_a, args.label_b = side_labels(
        rows_a, rows_b, args.label_a, args.label_b)

    res = compare(rows_a, rows_b, args.strip or None,
                  args.label_a, args.label_b)

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "comparison.json"), "w") as fh:
        json.dump(res, fh, indent=2)

    # the paired rows themselves, so a disagreement can be looked at rather
    # than only counted
    def index(rows):
        out = {}
        for r in rows:
            out.setdefault(
                (subject_of(r["sample_id"], args.strip or None), r["seq_id"]),
                []).append(r)
        return out

    ia, ib = index(rows_a), index(rows_b)
    with open(os.path.join(args.outdir, "paired_calls.tsv"), "w",
              newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["subject", "seq_id", "taxon",
                    f"call_{args.label_a}", f"call_{args.label_b}",
                    f"reads_{args.label_a}", f"reads_{args.label_b}",
                    f"depth_{args.label_a}", f"depth_{args.label_b}",
                    f"copies_{args.label_a}", f"copies_{args.label_b}",
                    "log10_ratio", "agree"])
        for key in sorted(set(ia) & set(ib)):
            if len(ia[key]) > 1 or len(ib[key]) > 1:
                continue
            a, b = ia[key][0], ib[key][0]
            va, vb = a.get("copies_per_unit_dna"), b.get("copies_per_unit_dna")
            ratio = ""
            if va and vb and float(va) > 0 and float(vb) > 0:
                ratio = f"{math.log10(float(vb) / float(va)):.4f}"
            w.writerow([key[0], key[1], a.get("taxon", ""),
                        a["call"], b["call"], a.get("reads", ""),
                        b.get("reads", ""), a.get("depth", ""),
                        b.get("depth", ""), va or "", vb or "", ratio,
                        # Same rule as the summary: a failed well is not
                        # a call, so it neither agrees nor disagrees.
                        ("excluded" if "assay_failed" in (a["call"], b["call"])
                         else "yes" if a["call"] == b["call"] else "no")])

    if not args.quiet:
        print(f"{args.label_a} vs {args.label_b}: "
              f"{res['paired_observations']:,} paired observations")
        if res["ambiguous_pairings"]:
            print(f"    {res['ambiguous_pairings']} skipped as ambiguous "
                  "(a subject measured more than once on one side)")
        if not res["paired_observations"]:
            print("    NOTHING PAIRED: no sample_id and seq_id is present in "
                  "both runs, so there is no agreement to measure. Check "
                  "that the ids match, or set --strip (replicate_suffix) "
                  "to remove a replicate suffix")
            return 0
        if res["excluded_assay_failed"]:
            print(f"    {res['excluded_assay_failed']} left out of agreement "
                  "and kappa: assay_failed on at least one side")
        if res["agreement"] is None:
            print("    agreement undefined: every pairing had a failed well")
        else:
            print(f"    agreement {res['agreement']:.3f}", end="")
            print(f"   kappa {res['cohen_kappa']:.3f}"
                  if res["cohen_kappa"] is not None
                  else "   kappa undefined (no disagreement possible)")
        if res["n_both_detected"]:
            q = res["log10_copy_ratio"]
            print(f"    both detected: {res['n_both_detected']}, "
                  f"{res['within_twofold']:.1%} within twofold, "
                  f"median log10 ratio {q.get('p50', 0):+.3f}")
        if res["detected_by_one_only"]:
            print("    found by one side only: "
                  + ", ".join(f"{k} {v}" for k, v in
                              sorted(res["detected_by_one_only"].items())))
            d = res["depth_of_disagreements"]
            if d:
                print(f"       of {d['n_disagreements_with_depth']} with "
                      "depth recorded, the deeper run found it "
                      f"{d['found_by_deeper_run']} times, the shallower "
                      f"{d['found_by_shallower_run']}, at equal depth "
                      f"{d['equal_depth']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
