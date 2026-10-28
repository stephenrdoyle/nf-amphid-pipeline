#!/usr/bin/env python3
"""
Database-mode assignment for Illumina: one sample's read pairs against a
reference database release.

The nanopore route (assign_db.py) builds a consensus because a nanopore
read is too noisy to name a species: sister taxa differ by ~1% over V4
and a read errs at 2-3%. An Illumina read errs at ~0.1%, so the noise is
not in the bases but in the counting - a handful of reads carrying one
error each become a handful of extra "sequences". What this route needs
is the amplicon-sequence-variant step VESPA took with DADA2: collapse the
error copies onto the sequence they came from, then name each surviving
sequence. Then the same naming rules as nanopore - LCA over hits within
0.5% of the best score, identity floors per rank - and the same outputs,
so synthesis and the report do not know which platform they are reading.

Why the pairs are not merged. VESPA's own pairs did not overlap at 2x250
(they used DADA2's justConcatenate), and at 2x150 99.9% of V4 inserts
leave an unread gap between the mates (docs/REFERENCE_DATABASE.md, §4).
Dropping pairs that fail to merge would bias silently against long-V4
lineages - cestodes reach 970bp. So each pair stays two windows: the
first `--window-r1` bases of the insert from R1 and the last
`--window-r2` from R2, a fixed length each so that error copies line up
base for base. A pair whose insert is SHORTER than the read runs through
into the opposite primer; that is detected (prep_illumina's read-through
search), and the pair then holds the whole insert in R1 and is kept as
one sequence - short-V4 lineages are not lost either.

Per sample:

  1. Pair the mates, find the forward primer at the start of one and the
     reverse at the start of the other (either way round), trim both.
  2. Cut read-through; cut each mate to its window; drop pairs over
     --max-ee expected errors in the windows.
  3. Dereplicate; then denoise, UNOISE-style: in decreasing abundance, a
     sequence d mismatches from a more abundant one and at most
     1/2^(alpha*d + 1) of its abundance is that sequence's error copy
     (alpha = 2: an eighth at one mismatch, a 32nd at two). Sequences
     seen fewer than --min-size times found nothing of their own; they
     join the nearest surviving sequence within --max-diffs or stay
     unclustered.
  4. Map each window with minimap2 (-x sr), keep records hit by BOTH
     windows in the orientation a pair has (R1 forward near the record's
     start, R2 reverse near its end), and score a record by the two
     alignment scores together.
  5. LCA and floors as on nanopore, on the identity of both windows
     together (unaligned ends counted as differences).

Expected errors come from the quality string; on SRA-Lite data, where
every base reads '?' (Q30), that filter is blind and says so.

Standard-library Python plus minimap2.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
import tempfile
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    RANKS, load_assay, load_primers, load_release_manifest, load_taxonomy,
    read_table, release_version, revcomp, writer,
)
import assign_db as adb  # noqa: E402
from prep_illumina import paired_reads, trim_readthrough  # noqa: E402


# Floors for pairs, stricter than nanopore's. Over the two 2x150 windows
# (~260bp) the bases that separate sister taxa are mostly in the gap the
# pair never reads, so the nanopore floors - fractions calibrated on
# whole inserts - let a missing species take its sister's name. On a
# 150-sample leave-out run of v0.4.0 at 2x150 (tests/bench_db_floors.py
# --illumina 150), species 0.997 (an exact match over 260bp) cut wrong
# names on missing species from 16% to 9% at a cost of 2.6 points of
# species calls on present ones; genus 0.995 cut wrong genera on missing
# genera from 18% to 17%. docs/REFERENCE_DATABASE.md section 11.5.
ILLUMINA_FLOORS = {"species": 0.997, "genus": 0.995, "family": 0.95}


def expected_errors(qual: str) -> float:
    return sum(10 ** (-(ord(c) - 33) / 10) for c in qual)


def orient(r1, r2, ps):
    """((r1 seq, qual), (r2 seq, qual)) with primers trimmed, R1 the
    forward-primer mate; None when no primer pair is found."""
    for a, b in ((r1, r2), (r2, r1)):
        f = ps.match(a[1].upper(), "fwd")
        if f is None:
            continue
        r = ps.match(b[1].upper(), "rev")
        if r is None:
            continue
        return ((a[1][f[1]:].upper(), a[2][f[1]:]),
                (b[1][r[1]:].upper(), b[2][r[1]:]), f[0], r[0])
    return None


class Denoiser:
    """UNOISE-style denoising over fixed-length pairs (module docstring).

    Distances are Hamming within a length class - Illumina errors are
    substitutions, and the windows are cut to fixed lengths so they line
    up. Candidates are found exactly by the pigeonhole principle: a
    sequence within D mismatches of another shares at least one of D + 1
    position-anchored pieces with it, so only centroids sharing a piece
    are compared at all.
    """

    def __init__(self, max_diffs: int = 3, alpha: float = 2.0,
                 min_size: int = 2):
        self.D, self.alpha, self.min_size = max_diffs, alpha, min_size
        self.centroids: list[str] = []
        self.size: list[int] = []
        self.index: dict = defaultdict(list)

    def _pieces(self, s: str):
        n, k = len(s), self.D + 1
        cuts = [n * i // k for i in range(k + 1)]
        return [(n, i, s[cuts[i]:cuts[i + 1]]) for i in range(k)]

    def nearest(self, s: str):
        """(centroid, mismatches) of the closest centroid within D, the
        most abundant on a tie; None when there is none."""
        seen, best = set(), None
        for key in self._pieces(s):
            for c in self.index.get(key, ()):
                if c in seen:
                    continue
                seen.add(c)
                t = self.centroids[c]
                d = sum(1 for x, y in zip(s, t) if x != y)
                if d <= self.D and (best is None or (d, -self.size[c])
                                    < (best[1], -self.size[best[0]])):
                    best = (c, d)
        return best

    def run(self, counts: Counter) -> tuple[dict, int]:
        """Map every unique to a centroid index; unclustered reads are
        counted, not mapped. Returns (unique -> centroid, unclustered)."""
        member: dict[str, int] = {}
        unclustered = 0
        for s, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            near = self.nearest(s)
            if n >= self.min_size:
                if near is not None:
                    c, d = near
                    if n <= self.size[c] / 2 ** (self.alpha * d + 1):
                        member[s] = c
                        continue
                c = len(self.centroids)
                self.centroids.append(s)
                self.size.append(n)
                for key in self._pieces(s):
                    self.index[key].append(c)
                member[s] = c
            elif near is not None:
                member[s] = near[0]
            else:
                unclustered += n
        return member, unclustered


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--database", required=True)
    ap.add_argument("--index", help="minimap2 index of refdb.fasta "
                    "(default: build on the fly)")
    ap.add_argument("--primers", required=True)
    ap.add_argument("--assay", help="for the adapter, to cut read-through")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", default="assign_db")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--window-r1", type=int,
                    help="insert bases kept from R1 (default: the release's "
                    "illumina_windows)")
    ap.add_argument("--window-r2", type=int)
    ap.add_argument("--max-ee", type=float, default=2.0,
                    help="expected errors allowed over both windows")
    ap.add_argument("--min-size", type=int, default=2,
                    help="reads a sequence needs to stand on its own")
    ap.add_argument("--max-diffs", type=int, default=3,
                    help="mismatches within which error copies and rare "
                    "sequences are joined to a more abundant one")
    ap.add_argument("--alpha", type=float, default=2.0)
    ap.add_argument("--lca-margin", type=float, default=0.005)
    ap.add_argument("--floors", default=",".join(
        f"{k}={v}" for k, v in ILLUMINA_FLOORS.items()),
        help="identity floor per rank over both windows together")
    ap.add_argument("--max-secondary", type=int, default=5000,
                    help="hits reported per window; the first window of V4 "
                    "is conserved enough to tie with thousands of records")
    args = ap.parse_args(argv)
    t0 = time.time()

    sid = args.sample_id
    os.makedirs(args.outdir, exist_ok=True)
    man = load_release_manifest(args.database)
    version = release_version(man)
    win = man.get("illumina_windows") or [132, 128]
    w1 = args.window_r1 or int(win[0])
    w2 = args.window_r2 or int(win[1])
    tax = load_taxonomy(os.path.join(args.database, "taxonomy.tsv"))
    floors = adb.parse_floors(args.floors)
    ps = load_primers(args.primers)
    cfg = load_assay(args.assay)
    # The adapter a read runs into past the insert: from R1, the reverse
    # overhang's reverse complement; from R2, the forward one's. Without
    # the assay, prep_illumina's defaults.
    ad1 = revcomp(cfg.get("pcr1_overhang_rev", ""))
    ad2 = revcomp(cfg.get("pcr1_overhang_fwd", ""))
    fwd_seq = {n: s for n, s in ps.fwd}
    rev_seq = {n: s for n, s in ps.rev}

    fate: Counter = Counter()
    counts: Counter = Counter()
    qual_values: set = set()
    for a, b in paired_reads(args.r1, args.r2):
        fate["reads_in"] += 1
        if len(qual_values) < 4:
            qual_values.update(a[2], b[2])
        o = orient(a, b, ps)
        if o is None:
            fate["no_primer"] += 1
            continue
        (s1, q1), (s2, q2), fn, rn = o
        t1, q1t = trim_readthrough(s1, q1, ad1, rev_seq.get(rn, ""))
        if len(t1) < len(s1):
            # The insert ended inside R1: R1 holds all of it.
            if len(t1) < 50:
                fate["too_short"] += 1
                continue
            if expected_errors(q1t) > args.max_ee:
                fate["low_quality"] += 1
                continue
            counts[(t1, "")] += 1
            fate["whole_insert"] += 1
            continue
        if len(s1) < w1 or len(s2) < w2:
            fate["too_short"] += 1
            continue
        if expected_errors(q1[:w1]) + expected_errors(q2[:w2]) > args.max_ee:
            fate["low_quality"] += 1
            continue
        counts[(s1[:w1], s2[:w2])] += 1

    # Denoise within each shape: windowed pairs together, whole inserts by
    # length (a pair key is joined with '|' so the two windows never
    # compare across their boundary).
    by_shape: dict = defaultdict(Counter)
    for (x, y), n in counts.items():
        by_shape[(len(x), len(y))][x + "|" + y] += n
    centroids: list[tuple[str, str, int]] = []
    unclustered = 0
    for shape, cnt in sorted(by_shape.items()):
        dn = Denoiser(args.max_diffs, args.alpha, args.min_size)
        member, lost = dn.run(cnt)
        unclustered += lost
        size = Counter()
        for s, c in member.items():
            size[c] += cnt[s]
        for c, s in enumerate(dn.centroids):
            x, _, y = s.partition("|")
            centroids.append((x, y, size[c]))
    fate["unclustered"] = unclustered
    fate["clustered"] = sum(n for _x, _y, n in centroids)

    mm = adb.Minimap2(args.minimap2, args.threads)
    index = args.index or os.path.join(args.database, "refdb.fasta")
    tmp = tempfile.mkdtemp(prefix=f"assign_db_pairs_{sid}_", dir=args.outdir)
    calls: dict[int, dict] = {}
    try:
        f1, f2, ff = (os.path.join(tmp, n) for n in ("r1.fa", "r2.fa", "full.fa"))
        adb.write_fasta(f1, ((f"p{k}", x) for k, (x, y, _n) in enumerate(centroids) if y))
        adb.write_fasta(f2, ((f"p{k}", y) for k, (x, y, _n) in enumerate(centroids) if y))
        adb.write_fasta(ff, ((f"p{k}", x) for k, (x, y, _n) in enumerate(centroids) if not y))
        # -p 0.5, not the nanopore route's 0.9: each window is kept with
        # every record that scores within half of its best, so that a
        # record second-best for BOTH windows - which is what a taxon
        # missing from the database often looks like - is still there to
        # be joined.
        extra = ("--secondary=yes", "-N", str(args.max_secondary), "-p", "0.5")
        hits: dict[str, dict] = {"r1": defaultdict(dict), "r2": defaultdict(dict),
                                 "full": defaultdict(dict)}
        for name, path in (("r1", f1), ("r2", f2), ("full", ff)):
            if os.path.getsize(path):
                for h in mm.map(index, path, preset=("-x", "sr"), extra=extra):
                    prev = hits[name][h["q"]].get(h["t"])
                    if prev is None or h["AS"] > prev["AS"]:
                        hits[name][h["q"]][h["t"]] = h
        for k, (x, y, _n) in enumerate(centroids):
            q = f"p{k}"
            if not y:
                calls[k] = adb.classify_hits(
                    [h for h in hits["full"].get(q, {}).values()
                     if h["strand"] == "+"], tax, floors, args.lca_margin)
                continue
            joined = []
            a_hits, b_hits = hits["r1"].get(q, {}), hits["r2"].get(q, {})
            for t, h1 in a_hits.items():
                h2 = b_hits.get(t)
                # The shape a pair has on its insert: R1 forward from the
                # start, R2 reverse back from the end.
                if h2 is None or h1["strand"] != "+" or h2["strand"] != "-":
                    continue
                if h1["ts"] > 30 or h2["tlen"] - h2["te"] > 30:
                    continue
                c1 = h1["qs"] + h1["qlen"] - h1["qe"]
                c2 = h2["qs"] + h2["qlen"] - h2["qe"]
                joined.append({
                    "t": t, "AS": h1["AS"] + h2["AS"],
                    "match": h1["match"] + h2["match"],
                    "alen": h1["alen"] + h2["alen"],
                    "qlen": h1["qlen"] + h2["qlen"], "qs": c1 + c2,
                    "qe": h1["qlen"] + h2["qlen"]})
            if joined:
                calls[k] = adb.classify_hits(joined, tax, floors,
                                             args.lca_margin)
                continue
            # No record carries both windows in a pair's shape: the two
            # halves point at different records (a chimera, or a lineage
            # the database has only parts of). Each window is classified
            # alone and the call is their LCA, at the lower identity - it
            # can say no more than both halves agree on.
            one = [adb.classify_hits(list(hs.values()), tax, floors,
                                     args.lca_margin)
                   for hs in (a_hits, b_hits) if hs]
            if len(one) < 2:
                calls[k] = one[0] if one else adb.classify_hits(
                    [], tax, floors, args.lca_margin)
                calls[k]["split_mates"] = True
                continue
            lin = [tax.get(c["best_hit"], {}) for c in one]
            depth = -1
            for di, r in enumerate(RANKS):
                if di <= min(c["depth"] for c in one) and lin[0].get(r) \
                        and lin[0].get(r) == lin[1].get(r):
                    depth = di
                else:
                    break
            worse = min(one, key=lambda c: c["identity"])
            calls[k] = {**worse, "depth": depth,
                        "rank": RANKS[depth] if depth >= 0 else "",
                        "lca_depth": depth, "split_mates": True}
        mm_version = mm.version()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # Outputs in assign_db.py's shape, so synthesis reads either platform.
    resolution = {}
    rpath = os.path.join(args.database,
                         f"resolution_illumina_{w1 + len(next(iter(fwd_seq.values()), ''))}"
                         f"x{w2 + len(next(iter(rev_seq.values()), ''))}.tsv")
    if os.path.exists(rpath):
        resolution = {r["id"]: r.get("deepest_rank", "") for r in read_table(rpath)}
    total = fate["clustered"]
    order = sorted(range(len(centroids)), key=lambda k: (-centroids[k][2], k))
    rows = []
    label_depth = RANKS.index("family")
    for n, k in enumerate(order, 1):
        x, y, cnt = centroids[k]
        call = calls[k]
        taxon, label = adb.label_for(call, tax, label_depth)
        lin = tax.get(call["best_hit"], {})
        d = call["depth"]
        rows.append({
            "sample_id": sid, "cluster": f"{sid}_p{n}", "reads": cnt,
            "rel_abundance": f"{cnt / total:.6f}" if total else "",
            "rank": call["rank"], "taxon": taxon, "label": label,
            "identity": f"{call['identity']:.4f}",
            "coverage": f"{call['coverage']:.3f}",
            "lineage": ";".join(lin.get(r, "") for r in RANKS[:d + 1]),
            "best_hit": call["best_hit"],
            "best_hit_lineage": ";".join(lin.get(r, "") for r in RANKS),
            "near_hits": call["near_hits"],
            "lca_rank": RANKS[call["lca_depth"]] if call["lca_depth"] >= 0 else "",
            "db_resolution": resolution.get(call["best_hit"], ""),
            "origin": ("whole_insert" if not y else
                       "pair_split" if call.get("split_mates") else "pair"),
            "merged": 0,
            "length": len(x) + len(y),
            "database_version": version,
            # The insert as far as it was read, gap marked by Ns.
            "consensus": x + ("N" * 10 + revcomp(y) if y else ""),
        })
    cols = ["sample_id", "cluster", "reads", "rel_abundance", "rank", "taxon",
            "label", "identity", "coverage", "lineage", "best_hit",
            "best_hit_lineage", "near_hits", "lca_rank", "db_resolution",
            "origin", "merged", "length", "database_version", "consensus"]
    with open(os.path.join(args.outdir, f"{sid}.clusters.tsv"), "w",
              encoding="utf-8") as fh:
        w = writer(fh, cols)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.outdir, f"{sid}.consensus.fasta"), "w",
              encoding="utf-8") as fh:
        for r in rows:
            fh.write(f">{r['cluster']} reads={r['reads']} rank={r['rank'] or 'none'} "
                     f"taxon={r['taxon'] or 'none'} identity={r['identity']}\n"
                     f"{r['consensus']}\n")
    agg: dict = {}
    for r in rows:
        a = agg.setdefault((r["rank"], r["taxon"], r["label"]),
                           {"reads": 0, "clusters": 0, "identity": 0.0,
                            "lineage": r["lineage"]})
        a["reads"] += r["reads"]
        a["clusters"] += 1
        a["identity"] = max(a["identity"], float(r["identity"]))
    with open(os.path.join(args.outdir, f"{sid}.taxa.tsv"), "w",
              encoding="utf-8") as fh:
        w = writer(fh, ["sample_id", "rank", "taxon", "label", "reads",
                        "rel_abundance", "clusters", "best_identity",
                        "lineage", "database_version"])
        w.writeheader()
        for (rank, taxon, label), a in sorted(agg.items(),
                                              key=lambda kv: -kv[1]["reads"]):
            w.writerow({"sample_id": sid, "rank": rank, "taxon": taxon,
                        "label": label, "reads": a["reads"],
                        "rel_abundance": f"{a['reads'] / total:.6f}" if total else "",
                        "clusters": a["clusters"],
                        "best_identity": f"{a['identity']:.4f}",
                        "lineage": a["lineage"], "database_version": version})
    summary = {
        "sample_id": sid, "mode": "database", "platform": "illumina",
        "database": {"version": version, "name": man.get("name", ""),
                     "path": os.path.abspath(args.database)},
        "minimap2": mm_version,
        "elapsed_s": round(time.time() - t0, 2),
        "reads": {k: fate.get(k, 0) for k in
                  ("reads_in", "no_primer", "too_short", "low_quality",
                   "clustered", "unclustered", "whole_insert")},
        "reads_note": (
            "Read pairs. reads_in = no_primer + too_short + low_quality + "
            "clustered + unclustered; whole_insert counts the pairs (within "
            "clustered or unclustered) whose insert was shorter than R1. "
            "rel_abundance is over clustered pairs: relative 18S abundance."),
        # SRA-Lite keeps two quality values (Q30 for a pass, Q3 for a
        # failed read): the expected-error filter then only removes the
        # failed reads, and cannot tell good bases from bad.
        "quality_values": len(qual_values),
        "quality_filter_blind": len(qual_values) <= 2,
        "clusters": len(rows),
        "haplotype_splits": 0,
        "settings": {"window_r1": w1, "window_r2": w2, "max_ee": args.max_ee,
                     "min_size": args.min_size, "max_diffs": args.max_diffs,
                     "alpha": args.alpha, "floors": floors,
                     "lca_margin": args.lca_margin,
                     "resolution_table": os.path.basename(rpath)
                     if resolution else ""},
        "rank_reached": dict(Counter(r["rank"] or "none" for r in rows)),
    }
    with open(os.path.join(args.outdir, f"{sid}.assign_db.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"{sid}: {fate['reads_in']} pairs, {len(rows)} sequences, "
          f"{fate['clustered']} clustered, {unclustered} unclustered, "
          f"{summary['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
