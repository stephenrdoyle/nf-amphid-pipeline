#!/usr/bin/env python3
"""
Database-mode assignment for Illumina: one sample's read pairs against a
reference database release.

This is the paired route of the protocol layer (docs/CONTRACTS.md): it
takes the pairs a platform adapter delivered for one sample (read type
paired_short), orients, windows and denoises them, and asks the
release's protocol (bin/protocols/, --protocol) to classify the result.
Pairing, windowing and denoising are shared; what a protocol changes is
said below, and lives in its module.

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

The nemits2 protocol (nemabiome ITS2, bin/protocols/nemits2.py) differs
in three places. Its windows count record bases, and R2
opens with the 28S stretch the records stop short of, so R2's window
takes that stretch too (illumina_tail: 130 + 83 at 2x150 is 130 + 130
read bases). Step 4 only proposes candidate records; each is scored
exactly over its own span, the stretch free (nemits2.span_scores).
And step 5 is the release's calling rule, the gap rule for ITS2
(nemits2.span_call). --primer-offset shortens both windows by the
spacer bases some protocols put before each primer (0-3 for NC1/NC2),
so a pair with a 3-base spacer still fills them. And mates that overlap
(most ITS2 pairs at 2x250) are merged into the insert and scored whole
(merge_pair): in two overlapping windows every difference in the overlap
is counted twice, and the AT-repeat stutter of an Oesophagostomum ITS2
(two bases) fell under the 99% a species call needs.

Expected errors come from the quality string; on SRA-Lite data, where
every base reads '?' (Q30), that filter is blind and says so.

Standard-library Python plus minimap2.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
import time
from collections import Counter, defaultdict
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    Minimap2, illumina_windows, load_assay, load_primers, load_rank_names,
    load_release_manifest, load_taxonomy, ranks_of, read_table,
    release_version, revcomp, write_assignment,
)
import assign_db as adb  # noqa: E402
from prep_illumina import paired_reads, trim_readthrough  # noqa: E402
from protocols.vespa_v4 import ILLUMINA_FLOORS  # noqa: E402


def expected_errors(qual: str) -> float:
    return sum(10 ** (-(ord(c) - 33) / 10) for c in qual)


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
    ap.add_argument("--read-type", choices=["paired_short"],
                    default="paired_short",
                    help="what the platform delivered (contract A)")
    ap.add_argument("--outdir", default="assign_db")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--window-r1", type=int,
                    help="insert bases kept from R1 (default: the windows of "
                    "the release's longest illumina_geometries entry within "
                    "--read-length, else its illumina_windows)")
    ap.add_argument("--window-r2", type=int)
    ap.add_argument("--read-length", type=int,
                    help="sequenced read length: picks the release's windows "
                    "and resolution table (2x301 reads take the 300 geometry)")
    ap.add_argument("--max-ee", type=float, default=2.0,
                    help="expected errors allowed over both windows")
    ap.add_argument("--min-size", type=int, default=2,
                    help="reads a sequence needs to stand on its own")
    ap.add_argument("--max-diffs", type=int, default=None,
                    help="mismatches within which error copies and rare "
                    "sequences are joined to a more abundant one (default: "
                    "1.5%% of the sequence, at least 3)")
    ap.add_argument("--alpha", type=float, default=2.0)
    ap.add_argument("--primer-max-edits", type=int, default=3,
                    help="edits allowed when finding each primer (see "
                    "assign_db.py)")
    ap.add_argument("--no-chimera-check", action="store_true",
                    help="do not test sequences for PCR chimeras")
    ap.add_argument("--lca-margin", type=float, default=0.005)
    ap.add_argument("--floors", default=",".join(
        f"{k}={v}" for k, v in ILLUMINA_FLOORS.items()),
        help="identity floor per rank over both windows together")
    adb.add_calling_args(ap)
    ap.add_argument("--primer-offset", type=int, default=0,
                    help="spacer bases a protocol may put before each primer "
                    "(nemabiome NC1/NC2: 3). Both windows are this much "
                    "shorter, so every pair fills them whatever its spacer")
    ap.add_argument("--no-merge", action="store_true",
                    help="span-scored releases: keep overlapping mates as two "
                    "windows instead of merging them into the insert")
    ap.add_argument("--no-indel-fold", action="store_true",
                    help="span-scored releases: keep sequences an indel or two "
                    "from a much more abundant one as their own, as the pair "
                    "denoiser alone does")
    ap.add_argument("--max-secondary", type=int, default=5000,
                    help="hits reported per window; the first window of V4 "
                    "is conserved enough to tie with thousands of records")
    args = ap.parse_args(argv)
    t0 = time.time()

    sid = args.sample_id
    os.makedirs(args.outdir, exist_ok=True)
    man = load_release_manifest(args.database)
    version = release_version(man)
    proto = adb.protocol_of(args, man)
    # The windows of the longest geometry the release resolved that the
    # reads can hold (250/200 at 2x300), unless set here.
    w1, w2, rtable = illumina_windows(man, args.read_length,
                                      args.window_r1, args.window_r2)
    # Read bases kept from each mate: the windows, R2's past-the-record
    # stretch (0 for V4), less any spacer allowance.
    rw1, rw2, tail = proto.read_windows(man, args.read_length, w1, w2,
                                        args.primer_offset)
    tax = load_taxonomy(os.path.join(args.database, "taxonomy.tsv"), man)
    ranks = ranks_of(tax)
    rank_names = load_rank_names(args.database)
    rule, floors = proto.calling(args, man, tax)
    ps = load_primers(args.primers)
    ps.max_edits = args.primer_max_edits
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
        o = proto.orient(a, b, ps, fate)
        if o is None:
            fate["no_primer"] += 1
            continue
        (s1, q1), (s2, q2), fn, opposite = o
        t1, q1t = trim_readthrough(s1, q1, ad1, opposite)
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
        # Mates the protocol scores as one merged insert (nemits2).
        j = proto.merge_mates(s1, q1, s2, q2, args)
        if j is not None:
            if expected_errors(j[1]) > args.max_ee:
                fate["low_quality"] += 1
                continue
            counts[(j[0], "")] += 1
            fate["merged"] += 1
            continue
        if len(s1) < rw1 or len(s2) < rw2:
            fate["too_short"] += 1
            continue
        if expected_errors(q1[:rw1]) + expected_errors(q2[:rw2]) > args.max_ee:
            fate["low_quality"] += 1
            continue
        counts[(s1[:rw1], s2[:rw2])] += 1

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
    # The protocol's own second pass (nemits2: variants an indel apart,
    # which the denoiser above, comparing within one length, never sees).
    centroids, folded = proto.fold_variants(centroids, args)
    if folded:
        fate["indel_folded"] = folded
    fate["unclustered"] = unclustered
    fate["clustered"] = sum(n for _x, _y, n in centroids)

    mm = Minimap2(args.minimap2, args.threads)
    index = args.index or os.path.join(args.database, "refdb.fasta")
    tmp = tempfile.mkdtemp(prefix=f"assign_db_pairs_{sid}_", dir=args.outdir)
    try:
        # The protocol's step: map, classify, refine.
        ctx = SimpleNamespace(mm=mm, index=index, tmp=tmp, tax=tax,
                              ranks=ranks, rule=rule, floors=floors,
                              args=args, tail=tail,
                              db_fasta=os.path.join(args.database,
                                                    "refdb.fasta"))
        calls, fate["remapped"] = proto.classify_pairs(ctx, centroids)
        mm_version = mm.version()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # Outputs in assign_db.py's shape, so synthesis reads either platform.
    resolution = {}
    # The resolution table of the windows used. A release that lists no
    # geometries names its tables by READ length: --read-length, or else
    # window + primer, which is right only when the windows fill the reads.
    if rtable is None and not man.get("illumina_geometries"):
        if args.read_length:
            rl1 = rl2 = args.read_length
        else:
            rl1 = w1 + len(next(iter(fwd_seq.values()), ""))
            rl2 = w2 + len(next(iter(rev_seq.values()), ""))
        rtable = f"resolution_illumina_{rl1}x{rl2}.tsv"
    rpath = os.path.join(args.database, rtable or "")
    if rtable and os.path.exists(rpath):
        # id -> (deepest rank, margin); see assign_db.margin_of.
        resolution = {r["id"]: (r.get("deepest_rank", ""), r.get("margin"))
                      for r in read_table(rpath)}
    total = fate["clustered"]
    order = sorted(range(len(centroids)), key=lambda k: (-centroids[k][2], k))
    cid = {k: f"{sid}_p{n}" for n, k in enumerate(order, 1)}
    # PCR chimeras (assign_db.chimera_scan), within each shape so the
    # windows line up base for base. A join in the unread gap - R1 from
    # one parent, R2 from another - is the commonest kind on Illumina,
    # and the concatenated windows show it as a breakpoint at the '|'.
    chim: dict = {}
    if not args.no_chimera_check:
        shapes: dict = defaultdict(list)
        for k in order:
            x, y, _n = centroids[k]
            shapes[(len(x), len(y))].append(k)
        for ks in shapes.values():
            found = adb.chimera_scan(
                [centroids[k][0] + "|" + centroids[k][1] for k in ks],
                [centroids[k][2] for k in ks], fixed=True)
            for q, (a, b, bp, *_r) in found.items():
                chim[ks[q]] = (ks[a], ks[b], bp)
    # Only joins of two genera are excluded (assign_db.parents_differ).
    kept_joins = {k: v for k, v in chim.items()
                  if not adb.parents_differ(calls[v[0]], calls[v[1]], tax)}
    chim = {k: v for k, v in chim.items() if k not in kept_joins}
    rows = []
    label_depth = ranks.index("family")
    for n, k in enumerate(order, 1):
        x, y, cnt = centroids[k]
        call = calls[k]
        taxon, label = adb.label_for(call, tax, label_depth, rank_names)
        ws = call.get("windows") or []
        if len(ws) == 2 and call["depth"] < min(w[2] for w in ws):
            # The windows point at different records, and the call is
            # only what they share - often just "Eukaryota". That is not
            # a tie the database cannot settle. BLASTing the largest of
            # these on VESPA's data found yeast and Blastocystis fragments
            # joined in the gap, and one apparent nematode that was no
            # nematode. But on the primate samples most are gregarines
            # (R2 97% to Leidyana) whose R1 aligns 74-81% to the same
            # record - barely above the 64% (p99 74%) any unrelated 18S
            # record gives, so a divergent region of one organism cannot
            # be told from a join. The label says what is known.
            label = "windows disagree: " + ", ".join(
                f"{r} {100 * i:.0f}% to "
                f"{adb.nearest_name(tax.get(h, {})) or h}"
                for r, (h, i, _d) in zip(("R1", "R2"), ws)) + (
                " (a chimera, an artefact, or a stretch the database has "
                "nothing close to)")
        lin = tax.get(call["best_hit"], {})
        d = call["depth"]
        chimera_of = ""
        if k in kept_joins:
            a, b, bp = kept_joins[k]
            chimera_of = f"{cid[a]}+{cid[b]}@{bp} within genus, kept"
        if k in chim:
            a, b, bp = chim[k]
            # In the gap when every base telling the parents apart takes
            # the left parent's in R1's window and the right one's in R2's.
            q = x + "|" + y
            pa = centroids[a][0] + "|" + centroids[a][1]
            pb = centroids[b][0] + "|" + centroids[b][1]
            inf = [i for i in range(len(q)) if pa[i] != pb[i]]
            gap = y and inf and all(
                (q[i] == pa[i]) == (i < len(x)) for i in inf
                if q[i] in (pa[i], pb[i]))
            where = "gap" if gap else str(bp)
            chimera_of = f"{cid[a]}+{cid[b]}@{where}"
            taxon, label, d = "", f"chimera of {cid[a]} and {cid[b]}", -1
            call = {**call, "rank": ""}
        rows.append({
            "sample_id": sid, "cluster": f"{sid}_p{n}", "reads": cnt,
            "rel_abundance": f"{cnt / total:.6f}" if total else "",
            "rank": call["rank"], "taxon": taxon, "label": label,
            "identity": f"{call['identity']:.4f}",
            "coverage": f"{call['coverage']:.3f}",
            "lineage": ";".join(lin.get(r, "") for r in ranks[:d + 1]),
            "best_hit": call["best_hit"],
            "best_hit_lineage": ";".join(lin.get(r, "") for r in ranks),
            "near_hits": call["near_hits"],
            "lca_rank": ranks[call["lca_depth"]] if call["lca_depth"] >= 0 else "",
            "db_resolution": resolution.get(call["best_hit"], ("",))[0],
            "db_margin": adb.margin_of(call["rank"], call["best_hit"],
                                       resolution),
            "origin": ("chimera" if k in chim else "whole_insert" if not y
                       else "pair_split" if call.get("split_mates") else "pair"),
            "merged": 0,
            "chimera_of": chimera_of,
            "length": len(x) + len(y),
            "database_version": version,
            # The insert as far as it was read, gap marked by Ns.
            "consensus": x + ("N" * 10 + revcomp(y) if y else ""),
            **proto.column_values(call, tax),
        })
    summary = {
        "sample_id": sid, "mode": "database", "protocol": proto.name,
        "platform": "illumina",
        "database": {"version": version, "name": man.get("name", ""),
                     "path": os.path.abspath(args.database)},
        "minimap2": mm_version,
        "elapsed_s": round(time.time() - t0, 2),
        "reads": {k: fate.get(k, 0) for k in
                  ("reads_in", "no_primer", "too_short", "low_quality",
                   "clustered", "unclustered", "whole_insert")
                  + proto.pair_fates},
        "reads_note": (
            "Read pairs. reads_in = no_primer + too_short + low_quality + "
            "clustered + unclustered; whole_insert counts the pairs (within "
            "clustered or unclustered) whose insert was shorter than R1. "
            "rel_abundance is over clustered pairs: relative 18S abundance."
            + proto.pair_note),
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
        "chimeras": {"sequences": len(chim),
                     "reads": sum(centroids[k][2] for k in chim),
                     "within_genus_kept": {
                         "sequences": len(kept_joins),
                         "reads": sum(centroids[k][2] for k in kept_joins)}},
        "settings": {"window_r1": w1, "window_r2": w2, "max_ee": args.max_ee,
                     "min_size": args.min_size,
                     "max_diffs": args.max_diffs if args.max_diffs is not None
                     else "1.5% of length, at least 3",
                     "alpha": args.alpha, "floors": floors,
                     "lca_margin": args.lca_margin,
                     **proto.window_settings(rw1, rw2, args),
                     **proto.calling_settings(rule, args),
                     "resolution_table": os.path.basename(rpath)
                     if resolution else ""},
        "rank_reached": dict(Counter(r["rank"] or "none" for r in rows)),
    }
    write_assignment(args.outdir, sid, rows, proto.extra_columns, version,
                     summary)
    print(f"{sid}: {fate['reads_in']} pairs, {len(rows)} sequences, "
          f"{fate['clustered']} clustered, {unclustered} unclustered, "
          f"{summary['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
