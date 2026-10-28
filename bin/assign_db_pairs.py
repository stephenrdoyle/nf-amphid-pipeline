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
     alignment scores together. A pair whose best hit on either window
     -x sr leaves more than 10% unaligned (a lineage too divergent for
     sr's penalties) is mapped again with nanopore's penalties.
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


# Window mapping. SR_DIVERGENT keeps sr's seeds and chaining but scores
# as map-ont does (mismatch 4, gap 4+2k), so an alignment carries on
# through a divergent stretch instead of stopping at it (main()).
SR = ("-x", "sr")
SR_DIVERGENT = ("-x", "sr", "-B4", "-O4,24", "-E2,1")
SR_CLIP = 0.10


def map_windows(mm, index, paths: dict, extra, preset) -> dict:
    """{window: {query: {record: best hit}}} for each FASTA in paths."""
    hits: dict[str, dict] = {name: defaultdict(dict) for name in paths}
    for name, path in paths.items():
        if os.path.getsize(path):
            for h in mm.map(index, path, preset=preset, extra=extra):
                prev = hits[name][h["q"]].get(h["t"])
                if prev is None or h["AS"] > prev["AS"]:
                    hits[name][h["q"]][h["t"]] = h
    return hits


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


_LOW_BITS: dict[int, int] = {}


def as_int(s: str) -> int:
    return int.from_bytes(s.encode("ascii"), "big")


def hamming_int(x: int, y: int, n: int) -> int:
    """Mismatches between two equal-length strings given as as_int().

    Python compares 450-character strings a character at a time; the
    denoiser does that hundreds of thousands of times per sample. XOR of
    the bytes as one integer, folded so each differing byte leaves one
    bit, counts the same thing in C - ~30x faster, same answer.
    """
    m = _LOW_BITS.get(n)
    if m is None:
        m = _LOW_BITS[n] = int.from_bytes(b"\x01" * n, "big")
    z = x ^ y
    z |= z >> 4
    z |= z >> 2
    z |= z >> 1
    return bin(z & m).count("1")


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
        self.cint: list[int] = []
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
        si, n = as_int(s), len(s)
        for key in self._pieces(s):
            for c in self.index.get(key, ()):
                if c in seen:
                    continue
                seen.add(c)
                d = hamming_int(si, self.cint[c], n)
                if d <= self.D and (best is None or (d, -self.size[c])
                                    < (best[1], -self.size[best[0]])):
                    best = (c, d)
        return best

    def _add(self, s: str, n: int) -> int:
        c = len(self.centroids)
        self.centroids.append(s)
        self.cint.append(as_int(s))
        self.size.append(n)
        for key in self._pieces(s):
            self.index[key].append(c)
        return c

    def rescue(self, left: Counter, member: dict) -> int:
        """Second chance for sequences no centroid took.

        A rare taxon with error-prone reads - AT-rich, long homopolymers,
        read in a 2x300 run's noisy R2 tails - can have no read seen twice:
        every one is unique, none may stand on its own, and the whole taxon
        went unclustered. Found on VESPA's human samples, where Endolimax
        nana carried 1.5% of a sample's pairs and vanished. So among the
        leftovers, the sequence with the most others within D mismatches
        founds a cluster of them if there are at least min_size + 1, and the
        cluster's sequence is their per-position majority - a consensus,
        as on nanopore, so no single read's errors become the call. One
        greedy pass, densest first. Returns reads still unclustered.
        """
        seqs = list(left)
        # D + 2 pieces, so a neighbour within D shares at least TWO of them
        # exactly; requiring two before comparing at all is what keeps
        # this affordable on 40,000 leftovers, where one shared piece -
        # a conserved stretch - made nearly every pair a candidate.
        def pieces2(x):
            n, k = len(x), self.D + 2
            cuts = [n * i // k for i in range(k + 1)]
            return [(n, i, x[cuts[i]:cuts[i + 1]]) for i in range(k)]
        idx: dict = defaultdict(list)
        for i, x in enumerate(seqs):
            for key in pieces2(x):
                idx[key].append(i)
        ints = [as_int(x) for x in seqs]
        near: dict[int, set] = {}
        for i, x in enumerate(seqs):
            shared = Counter(j for key in pieces2(x) for j in idx[key])
            L = len(x)
            near[i] = {j for j, k in shared.items() if k >= 2 and
                       hamming_int(ints[i], ints[j], L) <= self.D}
        # One greedy pass, densest neighbourhood first (weights taken once,
        # not re-scored after every cluster, which was quadratic).
        weight = {i: sum(left[seqs[j]] for j in near[i]) for i in near}
        alive = set(range(len(seqs)))
        for best in sorted(alive, key=lambda i: (-weight[i], i)):
            if best not in alive or weight[best] < self.min_size + 1:
                continue
            group = near[best] & alive
            reads = sum(left[seqs[j]] for j in group)
            if reads < self.min_size + 1:
                continue
            cols = zip(*(seqs[j] for j in group))
            cons = "".join(Counter(col).most_common(1)[0][0] for col in cols)
            c = self._add(cons, reads)
            for j in group:
                member[seqs[j]] = c
            alive -= group
        return sum(left[seqs[j]] for j in alive)

    def run(self, counts: Counter) -> tuple[dict, int]:
        """Map every unique to a centroid index; unclustered reads are
        counted, not mapped. Returns (unique -> centroid, unclustered)."""
        member: dict[str, int] = {}
        unclustered = 0
        left: Counter = Counter()
        for s, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            near = self.nearest(s)
            if n >= self.min_size:
                if near is not None:
                    c, d = near
                    if n <= self.size[c] / 2 ** (self.alpha * d + 1):
                        member[s] = c
                        continue
                member[s] = self._add(s, n)
            elif near is not None:
                member[s] = near[0]
            else:
                left[s] += n
        left = self.absorb_far(left, member)
        unclustered = self.rescue(left, member)
        return member, unclustered

    def absorb_far(self, left: Counter, member: dict) -> Counter:
        """Leftovers within 2D of a surviving sequence join it.

        Most of what no centroid took within D is not a rare taxon but a
        noisy copy of an abundant one - a 2x300 run's R2 tail can carry a
        dozen errors. On VESPA's largest sample that was 41,000 reads, and
        searching them against each other for rescue() took over ten
        minutes. Against the few hundred surviving sequences it takes
        seconds, and what is left for rescue() is what it is for: reads
        far from everything found so far."""
        if not self.centroids or not left:
            return left
        D2 = 2 * self.D
        def pieces(x):
            n, k = len(x), D2 + 1
            cuts = [n * i // k for i in range(k + 1)]
            return [(n, i, x[cuts[i]:cuts[i + 1]]) for i in range(k)]
        idx: dict = defaultdict(list)
        for c, x in enumerate(self.centroids):
            for key in pieces(x):
                idx[key].append(c)
        rest: Counter = Counter()
        for s, n in left.items():
            best = None
            si, L = as_int(s), len(s)
            for c in {c for key in pieces(s) for c in idx.get(key, ())}:
                d = hamming_int(si, self.cint[c], L)
                if d <= D2 and (best is None or d < best[0]):
                    best = (d, c)
            if best is None:
                rest[s] = n
            else:
                member[s] = best[1]
        return rest


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
    ap.add_argument("--max-diffs", type=int, default=None,
                    help="mismatches within which error copies and rare "
                    "sequences are joined to a more abundant one (default: "
                    "1.5%% of the sequence, at least 3)")
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
        # The radius grows with the sequence: three mismatches is 1% of a
        # 2x150 pair's 260bp, but 0.7% of 450bp from 2x300 reads, whose
        # R2 tails err far more - and there a fixed three left 28% of the
        # VESPA pairs unclustered against 9% at 2x150.
        D = args.max_diffs if args.max_diffs is not None else \
            max(3, round(0.015 * (shape[0] + shape[1])))
        dn = Denoiser(D, args.alpha, args.min_size)
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
        paths = {"r1": f1, "r2": f2, "full": ff}
        hits = map_windows(mm, index, paths, extra, SR)
        # -x sr's penalties (mismatch 8, gap open 12) are made for reads
        # of a genome they belong to. Against a taxon ~15% from its
        # nearest record - VESPA's NHP gregarines, Hexamita,
        # Hypotrichomonas - the extension stops at the first divergent
        # stretch, the window aligns as two pieces, and the pair was
        # joined on one 60bp piece: identity ~0.5 over the pair, best hit
        # a record that merely shares the conserved start. A pair sr
        # leaves more than SR_CLIP of a window unaligned is mapped again,
        # BOTH its windows, with nanopore's penalties on the same seeds,
        # and that mapping replaces the first; the rest keep sr's, so
        # ordinary calls and the floors calibrated on them are unchanged.
        # Judged on the window's best-scoring hit, the one its call rests
        # on: sr will often also give some record a gappy end-to-end
        # alignment scoring below the clean partial ones. A window with
        # no hit at all counts as unaligned too.
        def poor(name, q):
            hs = hits[name].get(q)
            if not hs:
                return True
            h = max(hs.values(), key=lambda h: (h["AS"], h["t"]))
            return h["qs"] + h["qlen"] - h["qe"] > SR_CLIP * h["qlen"]
        redo = {f"p{k}" for k, (x, y, _n) in enumerate(centroids)
                if (poor("r1", f"p{k}") or poor("r2", f"p{k}") if y
                    else poor("full", f"p{k}"))}
        if redo:
            seqs = {"r1": lambda x, y: x if y else "",
                    "r2": lambda x, y: y, "full": lambda x, y: "" if y else x}
            again_paths = {}
            for name, pick in seqs.items():
                again_paths[name] = paths[name] + ".divergent"
                adb.write_fasta(again_paths[name], (
                    (f"p{k}", pick(x, y)) for k, (x, y, _n) in
                    enumerate(centroids) if f"p{k}" in redo and pick(x, y)))
            again = map_windows(mm, index, again_paths, extra, SR_DIVERGENT)
            for name in paths:
                for q in redo:
                    hits[name].pop(q, None)
                    if q in again[name]:
                        hits[name][q] = again[name][q]
        fate["remapped"] = sum(centroids[int(q[1:])][2] for q in redo)
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
        # Pairs whose windows -x sr could not align end to end, mapped
        # again with softer penalties (SR_DIVERGENT).
        "remapped_divergent": fate.get("remapped", 0),
        # SRA-Lite keeps two quality values (Q30 for a pass, Q3 for a
        # failed read): the expected-error filter then only removes the
        # failed reads, and cannot tell good bases from bad.
        "quality_values": len(qual_values),
        "quality_filter_blind": len(qual_values) <= 2,
        "clusters": len(rows),
        "haplotype_splits": 0,
        "settings": {"window_r1": w1, "window_r2": w2, "max_ee": args.max_ee,
                     "min_size": args.min_size,
                     "max_diffs": args.max_diffs if args.max_diffs is not None
                     else "1.5% of length, at least 3",
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
