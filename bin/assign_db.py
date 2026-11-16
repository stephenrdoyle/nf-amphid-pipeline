#!/usr/bin/env python3
"""
Database-mode assignment of single reads: one sample's reads against a
taxonomically annotated reference database release.

This is the single-read route of the protocol layer (docs/CONTRACTS.md):
it takes the reads a platform adapter delivered for one sample, builds
consensus sequences, and asks the release's protocol (bin/protocols/,
--protocol) to classify them. Steps 1-6 below are shared by every
protocol; step 7 is the protocol's - what follows describes vespa_v4,
and nemits2 scores each record over its own span under the gap rule
instead (bin/protocols/nemits2.py). Read types: single_long_noisy, and
whatever else a protocol declares.

Panel mode asks "which of these few targets is this read?", and a read
can only ever be called something on the panel. Database mode asks "what
is in this sample?" against tens of thousands of 18S V4 records, and the
answer has to be honest about two limits that a panel never meets: the
database does not contain everything that will be amplified, and V4 does
not separate every pair of sister species. So the unit of classification
here is not the read, and the answer is not always a species.

Why not classify reads. Benchmarked with minimap2 against simulated reads
from genera and species held OUT of the index, an LCA over equally best
hits named 79-85% of reads from missing species or genera at species
level - a confident wrong name. No identity floor can fix that at read
level: at 3% error a read from a species that IS present sits at ~97%
identity, and a read from its missing sister at ~96%, because sister taxa
differ by about 1% over V4. The two distributions overlap. A consensus of
many reads has ~0.2% error, and there a species floor of ~0.995 kept 90%
of correct species calls while wrong calls on missing taxa fell to 5-8%.

So, per sample:

  1. Trim each read at the primers, which orients it and gives it the
     extent of a database insert.
  2. Map every read to the database (minimap2) and group reads by their
     primary hit.
  3. Build a consensus for each group, guided by its reference: pile the
     reads up on the reference through minimap2's cs tags, take the
     majority base at each position, re-align to that draft and repeat
     until it stops changing.
  4. Split hidden haplotypes by linkage. A group whose reads come from two
     sequences gives a consensus of the majority, and the minority is lost
     inside it. Their differences show as alleles that are too frequent to
     be read errors AND that occur together in the same reads: an error is
     independent of every other error, a haplotype's differences travel
     together. The largest linked set is the second haplotype's
     signature, its reads form a new group, and both are resolved again.
  5. Merge consensus sequences within one edit of each other: groups that
     were one taxon split across near-identical references converge on
     the same sequence. For PacBio reads (single_long_accurate) also
     absorb a consensus that differs from a much larger one only in the
     lengths of a few homopolymer runs (merge_homopolymer): CCS errors
     are homopolymer indels, reads that share one map to the same
     database record, and that group's consensus copies the error. Then
     give a smaller consensus the run lengths of the nearest much larger
     one where they differ only by such indels and substitutions
     (correct_homopolymer): a real haplotype whose group shared a CCS
     indel keeps its substitutions and loses the indel.
  6. Count reads per consensus by mapping every read to the final set, so
     reads from groups too small for a consensus of their own are counted
     too - or left unclustered, if they fit their own database hit better
     than any consensus (a taxon too rare here to build one).
  7. Classify each consensus: LCA over the database hits within 0.5% of
     the best alignment score, then cut back to the deepest rank whose
     identity floor the best hit clears. The identity is reported beside
     the rank, so a near-identical relative is visible, not hidden, and a
     consensus that clears no floor down to family is reported as
     "unclassified, NN% to X" rather than given X's name.

The idea in step 4 - that diversity within a group is evidence of a
second sequence when it is LINKED - is Stephen Doyle's. On simulated
samples it raised hidden-haplotype recovery from 0% to 37% and exact
consensus for present taxa from 87% to 98%, for about three false splits
per 30 samples; merging cut spurious consensus sequences from 33 to 5-10.

Output is relative 18S abundance, never organisms or copies: 18S copy
number varies widely between taxa.

Standard-library Python plus minimap2, the pipeline's one external
dependency, used in database mode only.
"""

from __future__ import annotations

import argparse
import math
import os
import random
import re
import shutil
import sys
import tempfile
import time
from collections import Counter, defaultdict
from itertools import combinations
from operator import add
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    INDEX_ARGS, RANKS, Minimap2, bounded_edit,
    iter_trimmed, load_primers, load_rank_names, load_release_manifest,
    load_taxonomy, rank_word, ranks_of, read_fasta, read_table,
    release_version, write_assignment, write_fasta,
)
import protocols  # noqa: E402
from protocols.nemits2 import (  # noqa: E402
    GAP_EDITS, GAP_EXACT, GAP_FLOORS, MIN_SPAN, TENTATIVE_IDENTITY)
from protocols.vespa_v4 import DEFAULT_FLOORS  # noqa: E402

CS = re.compile(r"(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)")




# --------------------------------------------------------------------------
# reference-guided consensus
# --------------------------------------------------------------------------


class Pileup:
    """The reads of one group laid over a draft sequence.

    base[p]   Counter of what the reads have at draft position p: a base,
              or '-' for a deletion;
    ins[p]    Counter of insertions placed before position p;
    cover[p]  how many reads span p;
    alleles   per read, its differences from the draft - substitutions
              and deletions only, keyed by draft position - which is what
              the linkage test compares between reads;
    span      per read, the draft interval it covers.
    """

    def __init__(self, ref: str, alns: list[dict]):
        n = len(ref)
        self.ref = ref
        self.base = [Counter() for _ in range(n)]
        self.ins = [Counter() for _ in range(n + 1)]
        self.cover = [0] * (n + 1)
        self.alleles: dict[str, dict[int, str]] = {}
        self.span: dict[str, tuple[int, int]] = {}
        for a in alns:
            p, al = a["ts"], {}
            start = p
            for op in CS.findall(a["cs"]):
                kind = op[0]
                if kind == ":":
                    for _ in range(int(op[1:])):
                        self.base[p][ref[p]] += 1
                        self.cover[p] += 1
                        p += 1
                elif kind == "*":
                    b = op[2].upper()
                    self.base[p][b] += 1
                    self.cover[p] += 1
                    al[p] = b
                    p += 1
                elif kind == "-":
                    for _ in range(len(op) - 1):
                        self.base[p]["-"] += 1
                        self.cover[p] += 1
                        al[p] = "-"
                        p += 1
                else:
                    self.ins[p][op[1:].upper()] += 1
            self.alleles[a["q"]] = al
            self.span[a["q"]] = (start, p)

    def consensus(self, n_reads: int) -> str:
        """Majority base at every position at least half the reads span,
        with an insertion wherever most reads spanning it carry one."""
        out = []
        need = max(1, n_reads // 2)
        for p in range(len(self.ref)):
            if self.cover[p] < need:
                continue
            if self.ins[p] and sum(self.ins[p].values()) > self.cover[p] / 2:
                out.append(self.ins[p].most_common(1)[0][0])
            b = self.base[p].most_common(1)[0][0]
            if b != "-":
                out.append(b)
        end = len(self.ref)
        if self.ins[end] and sum(self.ins[end].values()) > \
                self.cover[end - 1] / 2:
            out.append(self.ins[end].most_common(1)[0][0])
        return "".join(out)

    def error_rates(self) -> tuple[float, float]:
        """(substitution, deletion) rates of the reads against the draft.

        These set what counts as too frequent to be an error in the
        linkage test, and they are measured, not assumed: the run's
        chemistry, basecaller and quality are all in them. A second
        haplotype inflates them slightly, which only makes the test more
        conservative."""
        sub = dele = tot = 0
        for p, c in enumerate(self.base):
            tot += self.cover[p]
            for b, k in c.items():
                if b == "-":
                    dele += k
                elif b != self.ref[p]:
                    sub += k
        if not tot:
            return 0.02, 0.01
        return (min(0.2, max(0.002, sub / tot)),
                min(0.2, max(0.001, dele / tot)))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), summed in log space."""
    if k <= 0:
        return 1.0
    if p <= 0:
        return 0.0
    lp, lq = math.log(p), math.log1p(-p)
    tot = 0.0
    for i in range(k, n + 1):
        tot += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1)
                        - math.lgamma(n - i + 1) + i * lp + (n - i) * lq)
    return tot


def in_homopolymer(ref: str, p: int) -> bool:
    return (p > 0 and ref[p - 1] == ref[p]) or \
        (p + 1 < len(ref) and ref[p + 1] == ref[p])


class Settings:
    """The consensus and splitting parameters, in one place, so the
    summary can record exactly what a sample was resolved with."""

    def __init__(self, **kw):
        self.min_reads = kw.get("min_reads", 5)
        self.max_reads = kw.get("max_reads", 300)
        self.hap_min_frac = kw.get("hap_min_frac", 0.10)
        self.hap_p = kw.get("hap_p", 1e-6)
        self.min_cooccur = kw.get("min_cooccur", 3)
        self.cooccur_ratio = kw.get("cooccur_ratio", 3.0)
        self.max_split_depth = kw.get("max_split_depth", 2)
        self.max_iter = kw.get("max_iter", 4)
        self.merge_edits = kw.get("merge_edits", 1)
        self.split = kw.get("split", True)

    def as_dict(self) -> dict:
        return dict(vars(self))


def signature(pu: Pileup, st: Settings) -> list[tuple[int, str]] | None:
    """The linked set of non-consensus alleles marking a second
    haplotype, or None.

    A candidate allele must be at least hap_min_frac of the reads at its
    position, and improbable as read error under the measured rate
    (binomial tail below hap_p). Deletions in homopolymers are passed
    over: that is where nanopore indels concentrate and where their
    placement smears between neighbouring positions.

    Two candidates are linked when at least min_cooccur reads spanning
    both carry both, and that is at least cooccur_ratio times what
    independence would give. Read errors are independent of each other;
    a haplotype's differences are not. The largest connected set of
    linked alleles is the signature - two or more alleles, because one
    frequent allele on its own cannot be told from a systematic error.
    """
    sub_rate, del_rate = pu.error_rates()
    ref = pu.ref
    floor = max(st.min_reads * 2, 10)
    cand = []
    for p in range(len(ref)):
        c = pu.cover[p]
        if c < floor:
            continue
        for b, k in pu.base[p].items():
            if b == ref[p] or k < max(st.min_cooccur, st.hap_min_frac * c):
                continue
            if b == "-" and in_homopolymer(ref, p):
                continue
            if binom_sf(k, c, sub_rate / 3 if b != "-" else del_rate) \
                    > st.hap_p:
                continue
            cand.append((p, b))
    if len(cand) < 2:
        return None
    carriers = {x: {q for q, al in pu.alleles.items() if al.get(x[0]) == x[1]}
                for x in cand}
    adj: dict = defaultdict(set)
    for x, y in combinations(cand, 2):
        lo, hi = min(x[0], y[0]), max(x[0], y[0])
        cov = {q for q, (s, e) in pu.span.items() if s <= lo and e > hi}
        if len(cov) < 10:
            continue
        a, b = carriers[x] & cov, carriers[y] & cov
        both = len(a & b)
        expect = len(a) * len(b) / len(cov)
        if both >= st.min_cooccur and both >= st.cooccur_ratio * max(expect, 0.5):
            adj[x].add(y)
            adj[y].add(x)
    best: list = []
    seen: set = set()
    for v in sorted(adj):
        if v in seen:
            continue
        comp, stack = set(), [v]
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            stack.extend(adj[u] - comp)
        seen |= comp
        if len(comp) > len(best):
            best = sorted(comp)
    return best if len(best) >= 2 else None


class ConsensusBuilder:
    """Reference-guided consensus for the groups of one sample."""

    def __init__(self, mm: Minimap2, reads: dict[str, str], tmp: str,
                 st: Settings):
        self.mm, self.reads, self.tmp, self.st = mm, reads, tmp, st
        self.splits = 0

    def build(self, qs: list[str], start: str) -> tuple[str, Pileup | None]:
        """Iterate pileup consensus from `start` until it is stable;
        return it and the reads' pileup against it."""
        draft = start
        qfa = os.path.join(self.tmp, "group.fa")
        tfa = os.path.join(self.tmp, "draft.fa")
        write_fasta(qfa, ((q, self.reads[q]) for q in qs))
        pu = None
        for _ in range(self.st.max_iter + 1):
            write_fasta(tfa, [("draft", draft)])
            alns = [a for a in self.mm.map(tfa, qfa, cs=True,
                                           extra=("--secondary=no",))
                    if a["tp"] == "P"]
            if not alns:
                return draft, None
            pu = Pileup(draft, alns)
            new = pu.consensus(len(alns))
            if new == draft or not new:
                return draft, pu
            draft = new
        return draft, pu

    def resolve(self, qs: list[str], start: str, depth: int = 0,
                origin: str = "") -> list[dict]:
        cons, pu = self.build(qs, start)
        if pu is None:
            return []
        if self.st.split and depth < self.st.max_split_depth:
            sig = signature(pu, self.st)
            if sig:
                need = (len(sig) + 1) // 2
                hap = [q for q in qs
                       if sum(1 for p, b in sig
                              if pu.alleles.get(q, {}).get(p) == b) >= need]
                hs = set(hap)
                rest = [q for q in qs if q not in hs]
                if (len(hap) >= max(self.st.min_reads,
                                    self.st.hap_min_frac * len(qs))
                        and len(rest) >= self.st.min_reads):
                    self.splits += 1
                    tag = f"{origin}split{depth + 1}"
                    return (self.resolve(hap, cons, depth + 1, tag + "a.")
                            + self.resolve(rest, cons, depth + 1, tag + "b."))
        return [{"consensus": cons, "reads": list(qs),
                 "origin": origin.rstrip(".") or "group"}]


def merge_close(clusters: list[dict], max_edits: int) -> list[dict]:
    """Merge consensus sequences within max_edits of each other, largest
    first, so the survivor is the better-supported sequence."""
    merged: list[dict] = []
    for c in sorted(clusters, key=lambda x: (-len(x["reads"]), x["consensus"])):
        for m in merged:
            if abs(len(m["consensus"]) - len(c["consensus"])) <= max_edits and \
                    bounded_edit(c["consensus"], m["consensus"], max_edits) \
                    <= max_edits:
                m["reads"].extend(c["reads"])
                m["merged"] += 1
                if "weight" in m:
                    m["weight"] += c["weight"]
                break
        else:
            merged.append({**c, "reads": list(c["reads"]), "merged": 0})
    return merged


# Absorbing homopolymer-length variants (merge_homopolymer), per read type:
# (largest total run-length difference, least abundance ratio). Only
# single_long_accurate, measured on S-BSST527 (RSII CCS, docs/CONTRACTS.md):
# there, reads grouped by their primary hit share CCS homopolymer indels, so
# a small group's consensus reproduces them. Nanopore consensus is untouched
# - its homopolymer errors are spread over every read of a group, not
# gathered into groups of their own, and its floors were set without this.
HOMOPOLYMER_MERGE = {"single_long_accurate": (4, 4.0)}
# ... and how far (edits) a smaller consensus may be from the larger one
# whose run lengths correct_homopolymer gives it; 0: no correction.
HOMOPOLYMER_REACH = {"single_long_accurate": 8}
# Dropping consensus sequences no read supports (unsupported), per read
# type: (fraction of the sample's best-supported exact-read rate, binomial
# p below which a consensus falls short of it). Only where reads are
# accurate enough for a real sequence to be read exactly by many of them.
EXACT_SUPPORT = {"single_long_accurate": (1 / 3, 1e-3)}


def runs(seq: str) -> tuple[str, list[int]]:
    """A sequence as its run bases and run lengths: AAGCCC -> AGC, 2 1 3."""
    bases, lens = [], []
    for b in seq:
        if bases and bases[-1] == b:
            lens[-1] += 1
        else:
            bases.append(b)
            lens.append(1)
    return "".join(bases), lens


def homopolymer_edits(a: str, b: str) -> int | None:
    """The indels between a and b when they differ only in the lengths of
    their homopolymer runs (a run never vanishes, no base changes); else
    None. A run of one counts: A -> AA is a homopolymer indel."""
    (ba, la), (bb, lb) = runs(a), runs(b)
    if ba != bb:
        return None
    return sum(abs(x - y) for x, y in zip(la, lb))


def merge_homopolymer(clusters: list[dict], max_edits: int,
                      min_ratio: float) -> list[dict]:
    """Absorb a consensus that differs from a much larger one only by
    homopolymer-length indels: at most max_edits of them, and the larger
    carrying at least min_ratio times its reads ("weight", the reads its
    group stood for, not just those drawn to build it). Largest first,
    against the survivors, so a variant joins the best-supported sequence
    it is a copy of. A substitution keeps two sequences apart whatever their
    abundance, and so does a ratio below min_ratio: a real haplotype
    differing only in a run length survives when it is common enough."""
    out: list[dict] = []
    for c in sorted(clusters, key=lambda x: (-x["weight"], x["consensus"])):
        for m in out:
            if m["weight"] < min_ratio * c["weight"] or \
                    abs(len(m["consensus"]) - len(c["consensus"])) > max_edits:
                continue
            e = homopolymer_edits(c["consensus"], m["consensus"])
            if e is not None and e <= max_edits:
                m["reads"].extend(c["reads"])
                m["merged"] += c["merged"] + 1
                m["weight"] += c["weight"]
                break
        else:
            out.append(c)
    return out


def align_ops(a: str, b: str, band: int) -> list[tuple[str, int, int]]:
    """A global unit-cost alignment of a to b within a diagonal band, as its
    edits: ('S', i, j) a[i] for b[j]; ('I', i, j) a[i] not in b, before
    b[j]; ('D', i, j) b[j] not in a, before a[i]. Empty when a == b."""
    n, m = len(a), len(b)
    band = max(band, abs(n - m)) + 1
    INF = n + m + 1
    D = [dict() for _ in range(n + 1)]
    for j in range(0, min(m, band) + 1):
        D[0][j] = j
    for i in range(1, n + 1):
        Di, Dp = D[i], D[i - 1]
        for j in range(max(0, i - band), min(m, i + band) + 1):
            best = Dp.get(j, INF) + 1
            if j:
                best = min(best, Dp.get(j - 1, INF) + (a[i - 1] != b[j - 1]),
                           Di.get(j - 1, INF) + 1)
            Di[j] = best
    ops = []
    i, j = n, m
    while i or j:
        here = D[i].get(j, INF)
        if i and j and here == D[i - 1].get(j - 1, INF) + (a[i - 1] != b[j - 1]):
            if a[i - 1] != b[j - 1]:
                ops.append(("S", i - 1, j - 1))
            i, j = i - 1, j - 1
        elif i and here == D[i - 1].get(j, INF) + 1:
            ops.append(("I", i - 1, j))
            i -= 1
        else:
            ops.append(("D", i, j - 1))
            j -= 1
    return ops[::-1]


def correct_runs(a: str, b: str, band: int) -> tuple[str, int]:
    """a with the homopolymer-length indels between it and b undone - each
    run they share takes b's length - and its substitutions kept; with how
    many bases were changed. Only when every indel between them is such a
    run length: the base gained or lost sits in a run of two or more in
    the sequence that has it, beside the same base in the other. Any other
    indel leaves a as it is. On S-BSST527, a looser rule (beside the same
    base in the other sequence only, which is how an alignment places an
    indel next to a substitution) moved as many reads' consensus away from
    DADA2's ASVs as towards them."""
    ins, dels = set(), {}
    for op, i, j in align_ops(a, b, band):
        if op == "S":
            continue
        if op == "I" and in_homopolymer(a, i) and (
                (j < len(b) and b[j] == a[i]) or (j > 0 and b[j - 1] == a[i])):
            ins.add(i)
        elif op == "D" and in_homopolymer(b, j) and (
                (i < len(a) and a[i] == b[j]) or (i > 0 and a[i - 1] == b[j])):
            dels.setdefault(i, []).append(b[j])
        else:
            return a, 0
    out = []
    for i, x in enumerate(a):
        out.extend(dels.get(i, ()))
        if i not in ins:
            out.append(x)
    out.extend(dels.get(len(a), ()))
    return "".join(out), len(ins) + sum(map(len, dels.values()))


def correct_homopolymer(clusters: list[dict], max_edits: int,
                        min_ratio: float, reach: int,
                        stats: Counter | None = None) -> list[dict]:
    """Undo, in a smaller consensus, the homopolymer-length indels that set
    it apart from the nearest much larger one (at least min_ratio times
    its reads, within reach edits), keeping its substitutions: a real
    haplotype whose group shared a CCS indel. At most max_edits bases are
    changed. A consensus corrected into one already present joins it.
    A "known" consensus - its group's database record, over the record's
    span - is left as it is: a run length a record carries is not a read
    error. (On S-BSST527 that kept a consensus identical to a DADA2 ASV
    whose TTT a 6x larger neighbour had as TT.) The test is the record
    found anywhere in the consensus, so where the record is partial, a
    run outside its span is protected too: fewer corrections, never a
    wrong one. stats["corrected"] counts every consensus corrected,
    including those that then join one already present."""
    out: list[dict] = []
    for c in sorted(clusters, key=lambda x: (-x["weight"], x["consensus"])):
        q = c["consensus"]
        near = None
        for m in out:
            if m["weight"] < min_ratio * c["weight"]:
                continue
            d = bounded_edit(q, m["consensus"], reach)
            if d <= reach and (near is None or d < near[0]):
                near = (d, m)
        if near and not c.get("known"):
            fixed, k = correct_runs(q, near[1]["consensus"], reach)
            if 0 < k <= max_edits:
                c = {**c, "consensus": fixed,
                     "corrected": c.get("corrected", 0) + k}
                if stats is not None:
                    stats["corrected"] += 1
        same = next((m for m in out if m["consensus"] == c["consensus"]),
                    None)
        if same:
            same["reads"].extend(c["reads"])
            same["merged"] += c["merged"] + 1
            same["weight"] += c["weight"]
        else:
            out.append(c)
    return out


def exact_reads(reads_fa: str, seqs: list[str]) -> Counter:
    """Per consensus (by index), the reads that are it exactly: the read,
    trimmed at its primers as the consensus was, begins with it (a read
    runs on past a consensus cut to its record's span). A read that begins
    with two counts to the longer."""
    by_len: dict[int, dict[str, int]] = defaultdict(dict)
    for k, q in enumerate(seqs):
        by_len[len(q)].setdefault(q, k)
    lens = sorted(by_len, reverse=True)
    out: Counter = Counter()
    for _q, seq in read_fasta(reads_fa):
        for n in lens:
            k = by_len[n].get(seq[:n])
            if k is not None:
                out[k] += 1
                break
    return out


def unsupported(counts: Counter, exact: Counter, frac: float,
                p: float) -> set[int]:
    """The consensus sequences too few of whose reads are them exactly.

    A sequence in the sample is read exactly by some share of its reads -
    the run's accuracy, which the sample's largest consensus measures. A
    consensus built from reads that share errors (noisy reads that one
    database record fits best, their errors falling where its alleles
    differ) is read exactly by almost none: on S-BSST527 (RSII CCS) 3% of
    their reads against 40%. One is dropped when its exact reads are
    improbably few (binomial, below p) at frac of the largest's rate. A
    small consensus cannot fall that short, so it stays."""
    if not counts:
        return set()
    top = max(counts, key=lambda k: (counts[k], -k))
    rate = frac * exact.get(top, 0) / counts[top]
    if rate <= 0:
        return set()
    return {k for k, n in counts.items()
            if n and 1 - binom_sf(exact.get(k, 0) + 1, n, rate) < p}


# --------------------------------------------------------------------------
# chimeras
# --------------------------------------------------------------------------


def match_flags(q: str, p: str, band: int = 25) -> list[bool]:
    """For each base of q, whether a global alignment to p matches it.

    Banded edit-distance alignment (unit costs) with a traceback; an
    inserted query base, or a substitution, is a mismatch at that base, and
    a deletion from the query is charged to the base after it. Equal-length
    sequences with no band needed are compared base for base, which is
    what Illumina's fixed windows are.
    """
    n, m = len(q), len(p)
    if n == m and band <= 0:
        return [a == b for a, b in zip(q, p)]
    band = max(band, abs(n - m) + 2)
    INF = n + m + 1
    # dp rows over j in [i - band, i + band]
    prev = {j: j for j in range(0, min(m, band) + 1)}
    back = [dict() for _ in range(n + 1)]
    for j in prev:
        back[0][j] = "L" if j else ""
    for i in range(1, n + 1):
        cur = {}
        lo, hi = max(0, i - band), min(m, i + band)
        for j in range(lo, hi + 1):
            best, how = INF, ""
            if j == 0:
                best, how = i, "U"
            else:
                d = prev.get(j - 1, INF) + (q[i - 1] != p[j - 1])
                if d < best:
                    best, how = d, "D" if q[i - 1] == p[j - 1] else "S"
                l = cur.get(j - 1, INF) + 1
                if l < best:
                    best, how = l, "L"
            u = prev.get(j, INF) + 1
            if u < best:
                best, how = u, "U"
            cur[j] = best
            back[i][j] = how
        prev = cur
    flags = [False] * n
    i, j = n, m
    pending_del = False
    while i > 0 or j > 0:
        how = back[i].get(j, "U" if i else "L")
        if how == "D":
            flags[i - 1] = not pending_del
            pending_del = False
            i, j = i - 1, j - 1
        elif how == "S":
            flags[i - 1] = False
            pending_del = False
            i, j = i - 1, j - 1
        elif how == "U":
            flags[i - 1] = False
            pending_del = False
            i -= 1
        else:
            pending_del = True
            j -= 1
    return flags


def chimera_scan(seqs: list[str], abund: list[int], fixed: bool = False,
                 skew: float = 2.0, max_model_diffs: int = 2,
                 min_gain: int = 3, max_parents: int = 6) -> dict:
    """Which of a sample's sequences are chimeras of two others.

    De novo, in the manner of UCHIME: a PCR chimera is made after its
    parents, so it is rarer than they are, and it is explained better by
    a join of two of them than by either alone. For each sequence, in
    decreasing abundance, the candidate parents are those at least `skew`
    times as abundant, the closest few by shared 8-mers. Each is aligned
    to it base for base (match_flags), and the best single parent is
    compared with the best two-parent model - left of a breakpoint from
    one, right of it from the other, either way round. A chimera when

      * the two-parent model leaves at most max_model_diffs differences;
      * it explains at least min_gain more bases than the best single
        parent - otherwise the sequence is simply a variant of one;
      * and each side of the breakpoint holds a base its own parent has
        and the other does not - both parents are really used.

    Returns {index: (parent_left, parent_right, breakpoint, diffs_single,
    diffs_model)} for each chimera, parents indexing `seqs`. A parent is
    never itself a chimera found here.
    """
    order = sorted(range(len(seqs)), key=lambda k: (-abund[k], k))
    kmers = [set(s[i:i + 8] for i in range(len(s) - 7)) for s in seqs]
    found: dict = {}
    for rank_q, q in enumerate(order):
        cand = [p for p in order[:rank_q]
                if p not in found and abund[p] >= skew * abund[q]]
        if len(cand) < 2:
            continue
        cand.sort(key=lambda p: -len(kmers[q] & kmers[p]))
        cand = cand[:max_parents]
        qs = seqs[q]
        L = len(qs)
        vecs = {p: (match_flags(qs, seqs[p], 0) if fixed and
                    len(seqs[p]) == L else match_flags(qs, seqs[p]))
                for p in cand}
        single = max(sum(v) for v in vecs.values())
        if L - single <= max_model_diffs:
            continue            # explained by one parent: not a chimera
        # Prefix sums (left parent) and suffix sums (right parent) once
        # per candidate, breakpoints 1..L-1; each pair's best breakpoint
        # is its first maximum, the first pair to reach the best wins.
        pre, suf = {}, {}
        for p in cand:
            v = vecs[p]
            run = [0]
            for x in v:
                run.append(run[-1] + x)
            pre[p] = run[1:L]
            suf[p] = [run[L] - c for c in run[1:L]]
        best = None
        for a in cand:
            pa = pre[a]
            for b in cand:
                if a == b:
                    continue
                tot = list(map(add, pa, suf[b]))
                sc = max(tot)
                if best is None or sc > best[0]:
                    best = (sc, a, b, tot.index(sc) + 1)
        sc, a, b, i = best
        model_diffs = L - sc
        if model_diffs > max_model_diffs or sc - single < min_gain:
            continue
        va, vb = vecs[a], vecs[b]
        left_ok = any(va[k] and not vb[k] for k in range(i))
        right_ok = any(vb[k] and not va[k] for k in range(i, L))
        if left_ok and right_ok:
            found[q] = (a, b, i, L - single, model_diffs)
    return found


def parents_differ(call_a: dict, call_b: dict, tax: dict,
                   by_best_hit: bool = False) -> bool:
    """Do a chimera's two parents belong to different genera (or higher)?

    Compared rank by rank down to genus, over the ranks BOTH parents were
    named to. A join of two genera is a PCR artefact that would otherwise
    read as a taxon, and is excluded. A join within one genus is kept: the
    test cannot tell it from real variation (18S is multi-copy, and a
    sample can carry genuine recombinant-looking haplotypes), V4 barely
    separates species, and the database's species labels are its least
    reliable part - Blastocystis ST2 x "Blastocystis sp." would otherwise
    be excluded purely because one parent's species is missing. A parent
    named only above genus counts as possibly the same genus, so doubt
    about the parents never removes reads - with one exception: when BOTH
    parents are "unclassified" (stopped above family because no record is
    close), the join is excluded. Nothing places either parent, so there
    is no genus to share, and on a real flow cell these were mostly joins
    between unplaced gut bacteria (decided 2026-10-29).

    by_best_hit (a protocol's chimera_parents_by_best_hit; nemits2): a
    parent named only above genus is judged by its best record's genus
    instead. On ITS2 such a parent is most often another chimera or an
    error copy of a known worm, "Rhabditida" only because its call fell
    under the floors, and the rule above kept every Haemonchus x
    Trichostrongylus join the check found in the zoo antelopes
    (PRJNA1137257) because of it.
    """
    ranks = ranks_of(tax)
    if by_best_hit:
        ga = tax.get(call_a.get("best_hit", ""), {}).get("genus", "")
        gb = tax.get(call_b.get("best_hit", ""), {}).get("genus", "")
        if ga and gb:
            return ga != gb
    lab_depth = ranks.index("family")
    if all(label_for(c, tax, lab_depth)[1].startswith("unclassified")
           for c in (call_a, call_b)):
        return True
    la = tax.get(call_a.get("best_hit", ""), {})
    lb = tax.get(call_b.get("best_hit", ""), {})
    da, db = call_a.get("depth", -1), call_b.get("depth", -1)
    for d in range(0, min(da, db, ranks.index("genus")) + 1):
        r = ranks[d]
        if la.get(r) and lb.get(r) and la[r] != lb[r]:
            return True
    return False


# --------------------------------------------------------------------------
# labels: what a call says, whichever protocol made it
# --------------------------------------------------------------------------


def deepest_named(lineage: dict, ranks: list[str] = RANKS) -> str:
    for r in reversed(ranks):
        if lineage.get(r):
            return lineage[r]
    return ""


# PR2 fills ranks that have no name with its parent's name and an _X per
# level (Chromadorea_X, Craniata_XX), and leaves a species unnamed as
# Genus_sp.
_PLACEHOLDER = re.compile(r"_X+$")


def is_placeholder(name: str) -> bool:
    return bool(_PLACEHOLDER.search(name or ""))


def nearest_name(lin: dict, ranks: list[str] = RANKS) -> str:
    """A record named by its deepest real name: genus if it has one, else
    up the tree past any placeholder ("90% to Prevotellaceae", not "90%
    to Prevotellaceae_X")."""
    above = tuple(reversed(ranks[:ranks.index("genus")]))
    return next((lin.get(r) for r in ("genus",) + above
                 if lin.get(r) and not is_placeholder(lin.get(r))),
                "") or deepest_named(lin, ranks)


def margin_of(call_rank: str, best_hit: str, resolution: dict) -> str:
    """How far a species call's best record sits from the nearest record
    of another species, in edits over this platform's view (the
    release's `margin` column, from v0.8.0).

    Only a species call has one. A species row the release leaves blank
    has no other species near, which is "large", not unknown. Releases
    before v0.8.0 have no column, and every margin is blank. About 30%
    of species-level records rest on 1-2 edits, so a species name on one
    of those is one base from being another species.
    """
    if call_rank != "species":
        return ""
    rank, margin = resolution.get(best_hit, ("", None))
    if margin is None:          # pre-v0.8.0 release: no margin column
        return ""
    return margin or ("large" if rank == "species" else "")



def label_for(call: dict, tax: dict, label_depth: int,
              rank_names: list | None = None) -> tuple[str, str]:
    """(taxon at the rank reached, label for the report).

    Down to label_depth (family) a call is the taxon's name. Above it
    the consensus cleared no floor that would let it be named: it is
    "unclassified, NN% to X", where X is the best hit's genus (or its
    deepest named rank) - so a gap in the database reads as a gap, and
    how far off the nearest record is stays on the page.

    Where the name at the rank reached is one of PR2's placeholders, the
    label says what it means instead of repeating it. "Chromadorea_X" in
    a report reads as a taxon; it is "no family name recorded, somewhere
    in Chromadorea" - on a real VESPA run it was mostly Haemonchus at
    100%, stopped at family by a tie with a second genus, and the label
    said nothing a reader could use. Now: "Chromadorea, unnamed family
    (nearest Haemonchus)". The taxon column keeps PR2's name, so tables
    still join on it.

    A species the tiered rule names but the gap rule does not is
    "<species> (tentative)", so it is never pooled with confident calls
    of the same species; a species or genus call whose nearest records
    all carry one of the release's lineages says it too ("Teladorsagia
    circumcincta (goat line)", "Necator (Necator sp. type II)").

    And a call stopped above family by a TIE between records - the floors
    let it through, the LCA did not - is "X, unresolved below <rank>
    (NN% to Y)", not "unclassified": the database holds the sequence and
    cannot choose between records, which is not the same as a gap.

    With the release's rank_names.tsv (v0.8.0 on) the rank in a label is
    the real one, not PR2's column: a nematode's "Chromadorea_X" in the
    family column is an unnamed ORDER, and an ascomycete stopped at the
    order column is "unresolved below subphylum".""" 
    ranks = ranks_of(tax)
    lin = tax.get(call["best_hit"], {})
    d = call["depth"]
    taxon = lin.get(ranks[d], "") if d >= 0 else ""
    lineage = [lin.get(r, "") for r in ranks]

    def word(k: int) -> str:
        return rank_word(rank_names or [], lineage, ranks[k])
    near = nearest_name(lin, ranks) or call["best_hit"] or "nothing"
    # Above genus, the nearest genus and how close: a family call reads
    # "Hexamitidae (97% to Enteromonas)". Without it the reader cannot
    # tell a sequence one base from Enteromonas, stopped by the 99.5%
    # Illumina genus floor, from one 5% away. Naming the genus outright
    # was weighed and refused: at a 95% genus floor, 41% of sequences
    # whose genus is missing from the release get a wrong one.
    pct = f"{100 * call.get('identity', 0.0):.0f}%"
    above_genus = d < ranks.index("genus")
    if d >= label_depth:
        if not is_placeholder(taxon):
            if above_genus and near and near != taxon and not is_placeholder(near):
                return taxon, f"{taxon} ({pct} to {near})"
            # Named by the tiered rule's lower floor only, or nearest records
            # of one line (the release's lineage: "goat line", "Necator sp.
            # type II"): said in the label, so each is its own row in every
            # table, never pooled with the species' or genus' other calls.
            notes = (["tentative"] if call.get("tentative") else []) + (
                [call["release_lineage"]] if call.get("release_lineage")
                and not above_genus else [])
            return taxon, f"{taxon} ({'; '.join(notes)})" if notes else taxon
        named = next((lin.get(ranks[k]) for k in range(d - 1, -1, -1)
                      if lin.get(ranks[k]) and not is_placeholder(lin.get(ranks[k]))),
                     "")
        nearest = "" if not near or is_placeholder(near) or near == named \
            else f" ({pct} to {near})"
        return taxon, f"{named or taxon}, unnamed {word(d)}{nearest}"
    if not call["best_hit"]:
        return taxon, "unclassified, no database hit"
    if call.get("lca_depth", -2) == d and d >= 0:
        # Not stopped by an identity floor but by records disagreeing
        # (or the two halves of a pair pointing at different ones):
        # the database has the sequence, it cannot say which it is.
        # "Unclassified" would claim a gap that is not there.
        named = next((lin.get(ranks[k]) for k in range(d, -1, -1)
                      if lin.get(ranks[k]) and not is_placeholder(lin.get(ranks[k]))),
                     taxon)
        return taxon, (f"{named}, unresolved below {word(d)} "
                       f"({100 * call['identity']:.0f}% to {near})")
    return taxon, f"unclassified, {100 * call['identity']:.0f}% to {near}"


# --------------------------------------------------------------------------
# the sample
# --------------------------------------------------------------------------


def add_calling_args(ap) -> None:
    """The protocol and calling-rule options both routes share."""
    ap.add_argument("--protocol", default="auto",
                    help="the database protocol (bin/protocols/): "
                    f"{', '.join(protocols.PROTOCOLS)}. auto: the release's "
                    "manifest 'protocol', else nemits2 for a release scored "
                    "over its records' span or called by the gap or tiered "
                    "rule, else vespa_v4")
    ap.add_argument("--calling-rule",
                    choices=("auto", "floors", "gap", "tiered"),
                    default="auto",
                    help="floors: an identity floor per rank (--floors); gap: "
                    "species only when the next species is --gap-edits "
                    "further off (--gap-floors); tiered: the gap rule's "
                    "calls plus a species at --tentative-identity with the "
                    "same gap, flagged tentative. auto: the release's "
                    "calling_rule, else tiered for a release scored over "
                    "its records' span (ITS2), else floors (V4)")
    ap.add_argument("--gap-floors", default=",".join(
        f"{k}={v}" for k, v in GAP_FLOORS.items()),
        help="identity floor per rank under the gap rule")
    ap.add_argument("--gap-edits", type=int, default=GAP_EDITS,
                    help="edits the nearest other species must be further "
                    "off for a species call under the gap rule")
    ap.add_argument("--gap-exact", type=int, default=GAP_EXACT,
                    help="under the gap rule, the edits the nearest other "
                    "species must be further off for a species call on an "
                    "exact match (0: exact matches need --gap-edits too)")
    ap.add_argument("--tentative-identity", type=float,
                    default=TENTATIVE_IDENTITY,
                    help="under the tiered rule, the identity at which a "
                    "species the gap rule stops at genus is named, flagged "
                    "tentative (the gap must be as clear)")
    ap.add_argument("--no-repeats", action="store_true",
                    help="span-scored releases: score simple repeats base by "
                    "base, not by copy number (repeat-aware scoring is the "
                    "release's default)")
    ap.add_argument("--min-species-span", type=int, default=MIN_SPAN,
                    help="record bases a span-scored species call must rest "
                    f"on (default {MIN_SPAN}, the least any record must "
                    "cover to count)")


def protocol_of(args, man: dict):
    """The release's protocol under these options (protocols.resolve);
    SystemExit naming the problem when there is none, or when it cannot
    take this read type."""
    try:
        proto = protocols.resolve(man, args.protocol, args.calling_rule)
    except ValueError as exc:
        raise SystemExit(f"--protocol: {exc}") from None
    why = proto.unsupported(args.read_type)
    if why:
        raise SystemExit(why)
    return proto


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", required=True, nargs="+",
                    help="the sample's FASTQ (several are read as one)")
    ap.add_argument("--read-type", choices=[t for t in protocols.READ_TYPES
                                            if t.startswith("single")],
                    default="single_long_noisy",
                    help="what the platform delivered (contract A)")
    ap.add_argument("--database", required=True,
                    help="reference database release directory")
    ap.add_argument("--index", help="minimap2 index of the release's "
                    "refdb.fasta, built with " + " ".join(INDEX_ARGS) +
                    " (default: build from refdb.fasta on the fly)")
    ap.add_argument("--primers", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", default="assign_db")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-cluster-reads", type=int, default=5,
                    help="reads a group needs for a consensus of its own")
    ap.add_argument("--max-consensus-reads", type=int, default=300,
                    help="reads drawn from a larger group to build its "
                    "consensus; every read is still counted")
    ap.add_argument("--haplotype-min-frac", type=float, default=0.10)
    ap.add_argument("--no-split", action="store_true",
                    help="do not split groups into haplotypes")
    ap.add_argument("--merge-edits", type=int, default=1)
    ap.add_argument("--homopolymer-merge-edits", type=int,
                    help="absorb a consensus differing from a larger one "
                    "only by this many homopolymer-length indels or fewer "
                    "(0: never). Default by read type: "
                    + ", ".join(f"{k} {v[0]}" for k, v in
                                HOMOPOLYMER_MERGE.items()) + ", else 0")
    ap.add_argument("--homopolymer-correct-reach", type=int,
                    help="give a smaller consensus the run lengths of a "
                    "larger one within this many edits, keeping its other "
                    "differences (0: never). Independent of "
                    "--homopolymer-merge-edits, but it takes the merge's "
                    "ratio, and its cap on bases changed (the read type's "
                    "default when the merge is off). Default by read type: "
                    + ", ".join(f"{k} {v}" for k, v in
                                HOMOPOLYMER_REACH.items()) + ", else 0")
    ap.add_argument("--exact-support", type=float,
                    help="drop a consensus whose reads are it exactly far "
                    "less often than this fraction of the sample's largest "
                    "consensus' rate (0: never). Default by read type: "
                    + ", ".join(f"{k} {v[0]:.3g}" for k, v in
                                EXACT_SUPPORT.items()) + ", else 0")
    ap.add_argument("--homopolymer-merge-ratio", type=float,
                    help="...when the larger has at least this many times "
                    "its reads. Default by read type: "
                    + ", ".join(f"{k} {v[1]}" for k, v in
                                HOMOPOLYMER_MERGE.items()))
    ap.add_argument("--lca-margin", type=float, default=0.005,
                    help="hits within this fraction of the best alignment "
                    "score take part in the LCA")
    ap.add_argument("--floors", default=",".join(
        f"{k}={v}" for k, v in DEFAULT_FLOORS.items()),
        help="identity floor per rank, e.g. species=0.995,genus=0.99")
    add_calling_args(ap)
    ap.add_argument("--recruit-margin", type=int, default=10,
                    help="a read is counted to a consensus unless its own "
                    "database hit scores more than this much better "
                    "(minimap2 alignment score)")
    ap.add_argument("--primer-max-edits", type=int, default=3,
                    help="edits allowed when finding each primer (panel mode "
                    "uses 2). 3, because off-target bacterial 16S differs "
                    "from the forward primer at its 3' end, a quarter of it "
                    "by three or more: at 2 a fifth of a real run's reads "
                    "were dropped as no_primer instead of counted as "
                    "bacteria")
    ap.add_argument("--no-chimera-check", action="store_true",
                    help="do not test consensus sequences for PCR chimeras")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)
    t0 = time.time()

    sid = args.sample_id
    os.makedirs(args.outdir, exist_ok=True)
    man = load_release_manifest(args.database)
    version = release_version(man)
    tax = load_taxonomy(os.path.join(args.database, "taxonomy.tsv"), man)
    ranks = ranks_of(tax)
    rank_names = load_rank_names(args.database)
    # id -> (deepest rank, margin). margin is None for a release without
    # the column (before v0.8.0), so "no margin recorded" and "no other
    # species near" stay different things.
    resolution = {r["id"]: (r.get("deepest_rank", ""), r.get("margin"))
                  for r in read_table(os.path.join(args.database,
                                                   "resolution_nanopore.tsv"))}
    db_fasta = os.path.join(args.database, "refdb.fasta")
    index = args.index or db_fasta
    proto = protocol_of(args, man)
    rule, floors = proto.calling(args, man, tax)
    st = Settings(min_reads=args.min_cluster_reads,
                  max_reads=args.max_consensus_reads,
                  hap_min_frac=args.haplotype_min_frac,
                  merge_edits=args.merge_edits, split=not args.no_split)
    hp_edits, hp_ratio = HOMOPOLYMER_MERGE.get(args.read_type, (0, 0.0))
    if args.homopolymer_merge_edits is not None:
        hp_edits = args.homopolymer_merge_edits
    if args.homopolymer_merge_ratio is not None:
        hp_ratio = args.homopolymer_merge_ratio
    support = EXACT_SUPPORT.get(args.read_type)
    if args.exact_support is not None:
        support = (args.exact_support, 1e-3) if args.exact_support > 0 \
            else None
    hp_reach = HOMOPOLYMER_REACH.get(args.read_type, 0) \
        if args.homopolymer_correct_reach is None \
        else args.homopolymer_correct_reach
    # The correction's cap on bases changed and its ratio are the merge's;
    # with the merge turned off they stay the read type's defaults, so
    # --homopolymer-merge-edits 0 does not silently turn the correction
    # off too.
    hp_cap, hp_cap_ratio = HOMOPOLYMER_MERGE.get(args.read_type, (4, 4.0))
    if hp_edits > 0:
        hp_cap = hp_edits
    if args.homopolymer_merge_ratio is not None or hp_edits > 0:
        hp_cap_ratio = hp_ratio
    mm = Minimap2(args.minimap2, args.threads)

    primers = load_primers(args.primers)
    primers.max_edits = args.primer_max_edits
    fate: Counter = Counter()
    tmp = tempfile.mkdtemp(prefix=f"assign_db_{sid}_", dir=args.outdir)
    clusters: list[dict] = []
    hp_merges = hp_corrected = 0
    dropped = {"sequences": 0, "reads": 0}
    # Per read, only the alignment score of its primary database hit: the
    # one field used later (recruitment, below). A deep well holds millions
    # of these, so a full hit record each was the bulk of its memory.
    db_score: dict[str, int] = {}
    try:
        # Trimmed reads go straight to disk; minimap2 reads them from there,
        # and only those chosen to build a consensus come back into memory.
        rfa = os.path.join(tmp, "reads.fa")
        n_reads = 0
        with open(rfa, "w", encoding="utf-8") as fh:
            for q, seq in iter_trimmed(args.reads, primers, fate):
                fh.write(f">{q}\n{seq}\n")
                n_reads += 1
        groups: dict[str, list[str]] = defaultdict(list)
        if n_reads:
            for a in mm.map(index, rfa, preset=INDEX_ARGS,
                            extra=("--secondary=no",)):
                if a["tp"] == "P" and a["q"] not in db_score:
                    db_score[a["q"]] = a["AS"]
                    groups[a["t"]].append(a["q"])
        fate["unmapped"] = n_reads - len(db_score)

        # Consensus per group large enough to have one. The reference
        # sequence is the starting draft; the database is read once, for
        # just the records that head a group.
        heads = {t for t, qs in groups.items() if len(qs) >= st.min_reads}
        refseq = {i: s for i, s in read_fasta(db_fasta) if i in heads}
        rng = random.Random(f"{args.seed}:{sid}")
        chosen = []
        # reads per read drawn: a group's consensus is built from at most
        # max_reads of them, but stands for them all.
        scale: dict[str, float] = {}
        for t in sorted(heads, key=lambda t: (-len(groups[t]), t)):
            qs = sorted(groups[t])
            if len(qs) > st.max_reads:
                qs = sorted(rng.sample(qs, st.max_reads))
            chosen.append((t, qs))
            scale[t] = len(groups[t]) / len(qs)
        # Only the reads that build a consensus come back into memory;
        # counting maps them all again from reads.fa.
        need = {q for _t, qs in chosen for q in qs}
        reads = {q: seq for q, seq in read_fasta(rfa) if q in need}
        del groups, need
        cb = ConsensusBuilder(mm, reads, tmp, st)
        built = []
        for t, qs in chosen:
            for c in cb.resolve(qs, refseq[t]):
                c["group"] = t
                c["weight"] = len(c["reads"]) * scale[t]
                c["known"] = refseq[t] in c["consensus"]
                built.append(c)
        clusters = merge_close(built, st.merge_edits)
        n_built = len(clusters)
        if hp_edits > 0:
            clusters = merge_homopolymer(clusters, hp_edits, hp_ratio)
        if hp_reach > 0:
            hp_stats: Counter = Counter()
            clusters = correct_homopolymer(clusters, hp_cap, hp_cap_ratio,
                                           hp_reach, hp_stats)
            hp_corrected = hp_stats["corrected"]
        # Every consensus absorbed: merged as a copy, or corrected into
        # one already present.
        hp_merges = n_built - len(clusters)

        # Count every read against the final set.
        # Each read is decided on its first primary hit as it streams
        # past; only the ids already decided are kept.
        def count(clusters: list[dict]) -> Counter:
            counts: Counter = Counter()
            seen: set[str] = set()
            if clusters:
                cfa = os.path.join(tmp, "consensus.fa")
                write_fasta(cfa, ((f"c{k}", c["consensus"])
                                  for k, c in enumerate(clusters)))
                for a in mm.map(cfa, rfa, extra=("--secondary=no",)):
                    if a["tp"] != "P" or a["q"] in seen:
                        continue
                    seen.add(a["q"])
                    own = db_score.get(a["q"])
                    if own is not None and \
                            own - a["AS"] > args.recruit_margin:
                        continue
                    counts[int(a["t"][1:])] += 1
            return counts
        counts = count(clusters)
        # A consensus no read supports is dropped, and its reads counted
        # again against the rest.
        if support:
            exact = exact_reads(rfa, [c["consensus"] for c in clusters])
            drop = unsupported(counts, exact, *support)
            if drop:
                dropped = {"sequences": len(drop),
                           "reads": sum(counts[k] for k in drop)}
                clusters = [c for k, c in enumerate(clusters)
                            if k not in drop]
                counts = count(clusters)
        n_clustered = sum(counts.values())
        fate["clustered"] = n_clustered
        fate["unclustered"] = n_reads - n_clustered

        # Classify the consensus sequences that kept any reads: the
        # protocol's step.
        keep = [k for k in range(len(clusters)) if counts[k] > 0]
        calls: dict[int, dict] = {}
        if keep:
            ctx = SimpleNamespace(mm=mm, index=index, db_fasta=db_fasta,
                                  tmp=tmp, tax=tax, ranks=ranks, rule=rule,
                                  floors=floors, args=args)
            calls = proto.classify_consensus(
                ctx, {k: clusters[k]["consensus"] for k in keep})
        mm_version = mm.version()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- outputs ---------------------------------------------------------
    label_depth = ranks.index("family")
    total = sum(counts.values())
    order = sorted(keep, key=lambda k: (-counts[k], clusters[k]["consensus"]))
    cid = {k: f"{sid}_c{n}" for n, k in enumerate(order, 1)}
    # PCR chimeras among the final consensus sequences (chimera_scan): a
    # chimera is kept and its reads counted, but it is not named - its
    # label says which two clusters it joins - and synthesis keeps it out
    # of every composition.
    chim = {} if args.no_chimera_check else chimera_scan(
        [clusters[k]["consensus"] for k in order], [counts[k] for k in order])
    chim = {order[q]: (order[a], order[b], bp) for q, (a, b, bp, *_r)
            in chim.items()}
    # Only joins of two genera are excluded (parents_differ); a join
    # within a genus is kept, named as usual, and marked.
    kept_joins = {k: v for k, v in chim.items()
                  if not parents_differ(calls[v[0]], calls[v[1]], tax)}
    chim = {k: v for k, v in chim.items() if k not in kept_joins}
    rows = []
    for n, k in enumerate(order, 1):
        c, call = clusters[k], calls[k]
        taxon, label = label_for(call, tax, label_depth, rank_names)
        lin = tax.get(call["best_hit"], {})
        d = call["depth"]
        chimera_of = ""
        if k in kept_joins:
            a, b, bp = kept_joins[k]
            chimera_of = f"{cid[a]}+{cid[b]}@{bp} within genus, kept"
        if k in chim:
            a, b, bp = chim[k]
            chimera_of = f"{cid[a]}+{cid[b]}@{bp}"
            taxon, label, d = "", f"chimera of {cid[a]} and {cid[b]}", -1
            call = {**call, "rank": ""}
        rows.append({
            "sample_id": sid, "cluster": f"{sid}_c{n}",
            "reads": counts[k],
            "rel_abundance": f"{counts[k] / total:.6f}" if total else "",
            "rank": call["rank"], "taxon": taxon, "label": label,
            "identity": f"{call['identity']:.4f}",
            "coverage": f"{call['coverage']:.3f}",
            "lineage": ";".join(lin.get(r, "") for r in ranks[:d + 1]),
            "best_hit": call["best_hit"],
            "best_hit_lineage": ";".join(lin.get(r, "") for r in ranks),
            "near_hits": call["near_hits"],
            "lca_rank": ranks[call["lca_depth"]] if call["lca_depth"] >= 0 else "",
            "db_resolution": resolution.get(call["best_hit"], ("",))[0],
            "db_margin": margin_of(call["rank"], call["best_hit"],
                                   resolution),
            "origin": "chimera" if k in chim else c["origin"],
            "merged": c["merged"],
            "chimera_of": chimera_of,
            "length": len(c["consensus"]),
            "database_version": version,
            "consensus": c["consensus"],
            **proto.column_values(call, tax),
        })
    summary = {
        "sample_id": sid,
        "mode": "database",
        # The protocol that made these calls (bin/protocols/), so a result
        # says how it was reached without the run's logs.
        "protocol": proto.name,
        "database": {"version": version, "name": man.get("name", ""),
                     "path": os.path.abspath(args.database)},
        "minimap2": mm_version,
        "elapsed_s": round(time.time() - t0, 2),
        "reads": {k: fate.get(k, 0) for k in
                  ("reads_in", "no_primer", "too_short", "unmapped",
                   "clustered", "unclustered")},
        "reads_note": (
            "reads_in = no_primer + too_short + clustered + unclustered. "
            "unmapped reads (no database hit) are counted within "
            "clustered or unclustered, depending on whether a consensus "
            "took them. rel_abundance in the cluster and taxon tables is "
            f"over clustered reads: relative {proto.marker} abundance, never "
            "organisms or copies."),
        "clusters": len(rows),
        "groups": len(heads),
        "haplotype_splits": cb.splits,
        "chimeras": {"sequences": len(chim),
                     "reads": sum(counts[k] for k in chim),
                     "within_genus_kept": {
                         "sequences": len(kept_joins),
                         "reads": sum(counts[k] for k in kept_joins)}},
        "minimap2_calls": mm.calls,
        # Only where it ran, so other read types' summaries are unchanged.
        **({"homopolymer_merges": hp_merges}
           if hp_edits > 0 or hp_reach > 0 else {}),
        **({"homopolymer_corrected": hp_corrected} if hp_reach > 0 else {}),
        **({"unsupported_dropped": dropped} if support else {}),
        "settings": {**st.as_dict(), "floors": floors,
                     "lca_margin": args.lca_margin,
                     "recruit_margin": args.recruit_margin,
                     "index_args": INDEX_ARGS,
                     **({"homopolymer_merge_edits": hp_edits,
                         "homopolymer_merge_ratio": hp_ratio}
                        if hp_edits > 0 else {}),
                     **({"homopolymer_correct_reach": hp_reach,
                         **({"homopolymer_correct_cap": hp_cap,
                             "homopolymer_correct_ratio": hp_cap_ratio}
                            if hp_edits <= 0 else {})}
                        if hp_reach > 0 else {}),
                     **({"exact_support": {"fraction": support[0],
                                           "p": support[1]}}
                        if support else {}),
                     **proto.calling_settings(rule, args)},
        "rank_reached": dict(Counter(r["rank"] or "none" for r in rows)),
    }
    write_assignment(args.outdir, sid, rows, proto.extra_columns, version,
                     summary)
    print(f"{sid}: {fate['reads_in']} reads, {len(rows)} clusters, "
          f"{fate['clustered']} clustered, {fate['unclustered']} unclustered, "
          f"{cb.splits} haplotype split(s), {summary['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
