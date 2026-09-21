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

    call_pairs: list[tuple[str, str]] = []
    ratios: list[float] = []
    only: Counter = Counter()
    disagreement_depth: list[tuple[int, int]] = []
    per_taxon: dict = defaultdict(lambda: {"n": 0, "agree": 0})

    for key in shared:
        a, b = ia[key][0], ib[key][0]
        ca, cb = a["call"], b["call"]
        call_pairs.append((ca, cb))
        t = a.get("taxon") or key[1]
        per_taxon[t]["n"] += 1
        per_taxon[t]["agree"] += int(ca == cb)

        if ca == "detected" and cb == "detected":
            va, vb = a.get("copies_per_unit_dna"), b.get("copies_per_unit_dna")
            if va and vb and float(va) > 0 and float(vb) > 0:
                ratios.append(math.log10(float(vb) / float(va)))
        elif "detected" in (ca, cb) and "assay_failed" not in (ca, cb):
            only[label_a if ca == "detected" else label_b] += 1
            try:
                disagreement_depth.append((int(a["depth"]), int(b["depth"])))
            except (KeyError, ValueError):
                pass

    n = len(call_pairs)
    agree = sum(1 for x, y in call_pairs if x == y)
    kappa = cohen_kappa(call_pairs)

    # Is the side that found it simply the deeper side? If so the platforms
    # are not disagreeing about biology.
    deeper_found = None
    if disagreement_depth:
        wins = sum(1 for da, db in disagreement_depth if da > db)
        deeper_found = {
            "n_disagreements_with_depth": len(disagreement_depth),
            f"{label_a}_deeper": wins,
            f"{label_b}_deeper": len(disagreement_depth) - wins,
        }

    return {
        "label_a": label_a,
        "label_b": label_b,
        "paired_observations": n,
        "ambiguous_pairings": len(ambiguous),
        "agreement": round(agree / n, 4) if n else None,
        "cohen_kappa": None if kappa is None else round(kappa, 4),
        "kappa_note": (
            "Raw agreement is inflated by how often both runs say "
            "not_detected. Kappa removes that; a value near zero means the "
            "agreement is no better than the marginal rates predict."
        ),
        "call_matrix": {f"{x} / {y}": c
                        for (x, y), c in Counter(call_pairs).most_common()},
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--a", nargs="+", required=True,
                    help="calls.tsv for the first run (or runs)")
    ap.add_argument("--b", nargs="+", required=True,
                    help="calls.tsv for the second")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
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
                        "yes" if a["call"] == b["call"] else "no"])

    if not args.quiet:
        print(f"{args.label_a} vs {args.label_b}: "
              f"{res['paired_observations']:,} paired observations")
        if res["ambiguous_pairings"]:
            print(f"    {res['ambiguous_pairings']} skipped as ambiguous "
                  "(a subject measured more than once on one side)")
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
                      f"{max(d[f'{args.label_a}_deeper'], d[f'{args.label_b}_deeper'])} times")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
