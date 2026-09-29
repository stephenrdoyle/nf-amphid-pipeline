#!/usr/bin/env python3
"""
Choose an index set for a platform whose dominant error is indels.

The index set caps what demultiplexing can ever recover, and the cap is
set before a single read is sequenced. On the reference MinION run it
bound hard: both sets were minimum EDIT distance 2, so zero indels were
safely correctable, and of 10.7M ends whose index was read but matched
nothing, 96% were equidistant from two declared indices once indels
were allowed. No matching rule could take them.

The measure that matters is EDIT distance, not Hamming. k errors of
either kind are safely correctable only when every pair is more than 2k
apart, so:

    min edit 2  ->  0 indels        min Hamming 4  ->  1 substitution
    min edit 3  ->  1 indel         min Hamming 5  ->  2 substitutions
    min edit 5  ->  2 indels

A set chosen for Hamming distance can be 5 apart by substitution and 2
by edit, because a single deletion shifts every base after it. That is
what happened here: AAGGAGTA vs ACTGCATA is Hamming 4, edit 2 - nearly
a frameshift of each other, which is the worst case on this platform.

Two modes:

  evaluate   report the geometry of a set you already have
  select     find the largest subset of a pool meeting a distance floor,
             which often means better indexing from primers already on
             the shelf rather than new oligos

  design_indices.py --index-map assets/index_map.tsv --min-edit 3
  design_indices.py --sequences AAAACCCC,TTTTGGGG,... --evaluate
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def edit(a: str, b: str) -> int:
    """True Levenshtein distance. Symmetric, and it must be.

    The first version of this used find_anchor_cost, the pipeline's
    read-level matcher. That is a SEMI-GLOBAL search - it aligns the
    pattern anywhere within the text with free end gaps - so it is
    asymmetric and far too permissive between two equal-length
    barcodes: edit(GTAGAGGA, AAGAGGCA) came back 2 one way and 3 the
    other. Selection and reporting then disagreed, and a subset chosen
    at "min edit 3" contained a pair 2 apart.

    Distance between two declared indices is a property of the pair,
    not of a search, so it needs the symmetric measure.
    """
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def geometry(seqs: list[str]) -> dict:
    if len(seqs) < 2:
        return {"n": len(seqs), "min_edit": None, "min_hamming": None,
                "safe_indels": None, "safe_subs": None, "closest": ("", "")}
    worst_e, worst_h, pair = 99, 99, ("", "")
    for i, a in enumerate(seqs):
        for b in seqs[i + 1:]:
            e = edit(a, b)
            if e < worst_e:
                worst_e, pair = e, (a, b)
            worst_h = min(worst_h, hamming(a, b))
    return {"n": len(seqs), "min_edit": worst_e, "min_hamming": worst_h,
            "safe_indels": max(0, (worst_e - 1) // 2),
            "safe_subs": max(0, (worst_h - 1) // 2),
            "closest": pair}


def largest_subset(pool: list[tuple[str, str]], min_edit: int
                   ) -> list[tuple[str, str]]:
    """Greedy, seeded from every starting index.

    The maximum-clique problem this really is, is NP-hard; seeding the
    greedy pass from each member in turn gets close enough for sets of
    this size and never overstates, because whatever it returns has
    been checked pairwise.
    """
    best: list[tuple[str, str]] = []
    for start in range(len(pool)):
        chosen = [pool[start]]
        for cand in pool:
            if cand in chosen:
                continue
            if all(edit(cand[1], c[1]) >= min_edit for c in chosen):
                chosen.append(cand)
        if len(chosen) > len(best):
            best = chosen
    return best


def load_pool(path: str) -> dict[str, list[tuple[str, str]]]:
    rows = [l for l in open(path, encoding="utf-8") if not l.startswith("#")]
    out: dict[str, list[tuple[str, str]]] = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        out.setdefault(r["role"], []).append(
            (r["name"], r["sequence"].upper()))
    return out


def report(label: str, seqs: list[tuple[str, str]]) -> None:
    g = geometry([s for _, s in seqs])
    if g["min_edit"] is None:
        print(f"  {label}: {g['n']} index - a single index constrains nothing")
        return
    print(f"  {label}: n={g['n']}  min edit {g['min_edit']}  "
          f"min Hamming {g['min_hamming']}")
    print(f"      safely correctable: {g['safe_indels']} indel(s), "
          f"{g['safe_subs']} substitution(s)")
    print(f"      closest pair: {g['closest'][0]} vs {g['closest'][1]}")
    if g["safe_indels"] < 1:
        print("      NOT USABLE on nanopore: cannot correct a single "
              "deletion, which is the commonest error there")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--index-map", help="TSV with role, name, sequence")
    ap.add_argument("--sequences",
                    help="comma-separated sequences, instead of a map")
    ap.add_argument("--min-edit", type=int, default=3,
                    help="distance floor for subset selection. 3 permits 1 "
                         "indel, 5 permits 2")
    ap.add_argument("--evaluate", action="store_true",
                    help="only report the geometry; select nothing")
    ap.add_argument("--need", type=int, default=0,
                    help="wells required. Reports whether the selected "
                         "i5 x i7 grid reaches it")
    args = ap.parse_args(argv)

    if args.sequences:
        seqs = [(s, s) for s in args.sequences.upper().split(",") if s]
        print("supplied set")
        report("set", seqs)
        if not args.evaluate:
            sub = largest_subset(seqs, args.min_edit)
            print(f"\nlargest subset at min edit >= {args.min_edit}: "
                  f"{len(sub)} of {len(seqs)}")
            report(f"subset", sub)
            print("      " + " ".join(s for _, s in sub))
        return 0

    if not args.index_map:
        raise SystemExit("give --index-map or --sequences")

    pool = load_pool(args.index_map)
    print(f"pool: {args.index_map}")
    for role in sorted(pool):
        report(role, pool[role])

    if args.evaluate:
        return 0

    print(f"\nlargest subset at min edit >= {args.min_edit} "
          f"(permits {max(0, (args.min_edit - 1) // 2)} indel(s))")
    picked = {}
    for role in sorted(pool):
        sub = largest_subset(pool[role], args.min_edit)
        picked[role] = sub
        print(f"\n  {role}: {len(sub)} of {len(pool[role])}")
        report(f"{role} subset", sub)
        for name, seq in sub:
            print(f"        {name}\t{seq}")

    roles = sorted(picked)
    if len(roles) == 2:
        wells = len(picked[roles[0]]) * len(picked[roles[1]])
        print(f"\n  grid: {len(picked[roles[0]])} x "
              f"{len(picked[roles[1]])} = {wells} wells")
        if args.need:
            ok = "enough" if wells >= args.need else "NOT ENOUGH"
            print(f"  {args.need} wells required: {ok}")
            if wells < args.need:
                print("  Options: a lower --min-edit (less indel tolerance), "
                      "longer indices, or more oligos.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
