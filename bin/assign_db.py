#!/usr/bin/env python3
"""
Database-mode assignment for nanopore: one sample's reads against a
taxonomically annotated reference database release.

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
     the same sequence.
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
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    RANKS, anchored_distance, bounded_edit, calling_rule, infix_distance,
    load_primers, load_rank_names, load_release_manifest, load_taxonomy,
    open_reads, rank_word, read_fasta, read_table, release_version,
    scores_over_reference, writer,
)

# minimap2's index parameters for the database. k=19, w=10 gave the same
# accuracy as the map-ont default k=15 on the leave-out benchmark and ran
# 1.6x faster; the index builds in 1-2 s, so it is built once per run.
INDEX_ARGS = ["-x", "map-ont", "-k", "19", "-w", "10"]

# Floors on the best hit's identity, per rank. A consensus is named at a
# rank only if its best hit is at least this close; otherwise it stops at
# the next rank up. From the leave-out benchmark at consensus error rates:
# species 0.995 kept 90% of correct species calls and cut wrong calls on
# missing species to 5-8%. Genus was 0.98 in the prototype, which let
# 27% of consensus from missing genera take a wrong genus on a 150-sample
# leave-out run of v0.4.0; 0.99 cut that to 20% at no cost to present
# taxa, the rest moving honestly to family (docs/REFERENCE_DATABASE.md
# section 11.3, tests/bench_db_floors.py).
DEFAULT_FLOORS = {"species": 0.995, "genus": 0.99, "family": 0.95}

# The gap rule, for releases scored over their records' span (ITS2):
# species needs >= 99% identity AND the best record of any other species
# >= 4 edits further off over the same span; genus >= 97%, family >= 95%.
# From nemabiome-its2-refdb's leave-one-out (3,164 queries, 1,210
# species): it named 52.9% of species present at 0.7% wrong, against
# 18.5% under the V4 floors; reads of missing species took a sister's
# name 6.6% of the time (V4 floors 2.4%, assignTaxonomy 42%). Per-genus
# or per-family limits did no better on held-out species.
GAP_FLOORS = {"species": 0.99, "genus": 0.97, "family": 0.95}
GAP_EDITS = 4
# Span scoring. A record must be aligned over at least MIN_SPAN bases to
# count as evidence at all (a 77-base fragment 32 edits off otherwise
# outranks full-length records by identity); among equally identical
# records, those within SPAN_SLACK bases of the longest span tie and the
# rest - partials matching over less of the read - drop out. A window's
# overhang past its record is free for at most WINDOW_TRIM bases beyond
# the 28S stretch. SPAN_CANDIDATES: records minimap2 proposes that are
# scored exactly, enough for every tie and the nearest other species.
MIN_SPAN = 100
SPAN_SLACK = 5
WINDOW_TRIM = 30
SPAN_CANDIDATES = 100

CS = re.compile(r"(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)")


# --------------------------------------------------------------------------
# minimap2
# --------------------------------------------------------------------------


class Minimap2:
    """Runs minimap2 and parses its PAF. One instance per sample, so the
    executable and thread count are decided once."""

    def __init__(self, exe: str = "minimap2", threads: int = 1):
        path = shutil.which(exe) or (exe if os.path.isfile(exe) else None)
        if not path:
            raise SystemExit(
                f"minimap2 not found ('{exe}'). Database mode needs it; see "
                "docs/USAGE.md, field setup")
        self.exe, self.threads = path, max(1, threads)
        self.calls = 0

    def version(self) -> str:
        out = subprocess.run([self.exe, "--version"], capture_output=True,
                             text=True)
        return out.stdout.strip()

    def map(self, target: str, query: str, extra=(), preset=("-x", "map-ont"),
            cs: bool = False):
        """Yield minimap2's hits one at a time, as it writes them.

        Streamed rather than collected: a deep well has millions of reads,
        and holding minimap2's whole output as text, then as a list of
        dicts, cost several GB before a single hit was used (STOP23 on
        vespa_run1 was killed at 16 GB). Callers keep only what they need.
        stderr goes to a file so a chatty minimap2 cannot fill a pipe and
        stall while stdout is being read."""
        args = [self.exe, *preset, "-c", "-t", str(self.threads),
                *(["--cs"] if cs else []), *extra, target, query]
        self.calls += 1
        with tempfile.TemporaryFile("w+") as err:
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=err,
                                    text=True, bufsize=1 << 20)
            stopped = False
            try:
                for line in proc.stdout:
                    f = line.rstrip("\n").split("\t")
                    tags = {t[:2]: t[5:] for t in f[12:]}
                    yield {
                        "q": f[0], "qlen": int(f[1]), "qs": int(f[2]),
                        "qe": int(f[3]), "strand": f[4], "t": f[5],
                        "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]),
                        "match": int(f[9]), "alen": int(f[10]),
                        "AS": int(tags.get("AS", 0)), "cs": tags.get("cs", ""),
                        "tp": tags.get("tp", "P"),
                    }
            finally:
                # A caller that stops early must not leave minimap2 running.
                if proc.poll() is None:
                    stopped = True
                    proc.stdout.close()
                    proc.kill()
                proc.wait()
            if proc.returncode != 0 and not stopped:
                err.seek(0)
                raise SystemExit(f"minimap2 failed ({' '.join(args)}):\n"
                                 f"{err.read()[-2000:]}")


def write_fasta(path: str, recs) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n{seq}\n")


def ranks_of(tax) -> list[str]:
    """The release's rank columns (ampseqlib.Taxonomy), PR2's by default."""
    return getattr(tax, "ranks", RANKS)


def identity(hit: dict) -> float:
    """Matching bases over the alignment block PLUS any of the query left
    unaligned at either end.

    minimap2's own ratio (column 10 over column 11) ignores the ends it
    clipped, so a consensus aligned over only part of its length would
    look as close as one aligned end to end. A consensus and a record are
    both primer-to-primer inserts, so an unaligned end is a real
    difference and is charged as one."""
    clip = hit["qs"] + (hit["qlen"] - hit["qe"])
    return hit["match"] / max(1, hit["alen"] + clip)


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
                break
        else:
            merged.append({**c, "reads": list(c["reads"]), "merged": 0})
    return merged


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
        best = None
        for a in cand:
            va = vecs[a]
            pre = [0]
            for x in va:
                pre.append(pre[-1] + x)
            for b in cand:
                if a == b:
                    continue
                vb = vecs[b]
                suf = [0] * (L + 1)
                for i in range(L - 1, -1, -1):
                    suf[i] = suf[i + 1] + vb[i]
                for i in range(1, L):
                    sc = pre[i] + suf[i]
                    if best is None or sc > best[0]:
                        best = (sc, a, b, i)
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


def parents_differ(call_a: dict, call_b: dict, tax: dict) -> bool:
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
    """
    ranks = ranks_of(tax)
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
# classification
# --------------------------------------------------------------------------


def lca_depth(names, tax: dict) -> int:
    """Index into RANKS of the deepest rank at which every record agrees;
    -1 when they disagree even at domain. A blank rank (a PR2
    placeholder) breaks agreement: no record can vouch for a name it
    does not have."""
    depth = -1
    for d, r in enumerate(ranks_of(tax)):
        vals = {tax.get(t, {}).get(r, "") for t in names}
        if len(vals) == 1 and "" not in vals:
            depth = d
        else:
            break
    return depth


def apply_floors(depth: int, ident: float, floors: dict,
                 ranks: list[str] = RANKS) -> int:
    """Cut a rank back until the identity clears its floor. Ranks
    without a floor are kept as the LCA gives them."""
    for r in ("species", "genus", "family", "order", "class"):
        if r in floors and depth >= ranks.index(r) and ident < floors[r]:
            depth = ranks.index(r) - 1
    return depth


def classify_hits(hits: list[dict], tax: dict, floors: dict,
                  margin: float) -> dict:
    """One consensus's call from its database hits."""
    ranks = ranks_of(tax)
    if not hits:
        return {"rank": "", "depth": -1, "identity": 0.0, "coverage": 0.0,
                "best_hit": "", "near_hits": 0, "lca_depth": -1}
    top = max(h["AS"] for h in hits)
    near = [h for h in hits if h["AS"] >= top * (1 - margin)]
    best = max(near, key=lambda h: (h["AS"], identity(h), h["t"]))
    ident = identity(best)
    names = sorted({h["t"] for h in near})
    # A record with no taxonomy at all - v0.4.0 has nine, conflicts its
    # build left unresolved - can vouch for nothing, and as a member of
    # the LCA it vetoed everything: human 18S tied with one stopped at
    # no rank. It is left out of the LCA, unless it is all there is.
    known = [t for t in names if tax.get(t, {}).get(ranks[0])]
    lca = lca_depth(known or names, tax)
    depth = apply_floors(lca, ident, floors, ranks)
    return {
        "depth": depth, "rank": ranks[depth] if depth >= 0 else "",
        "identity": ident,
        "coverage": (best["qe"] - best["qs"]) / max(1, best["qlen"]),
        "best_hit": best["t"], "near_hits": len(names),
        "lca_depth": lca,
    }


def aligned_identity(parts: list[tuple[str, str]], record: str,
                     max_frac: float = 0.35) -> float | None:
    """Identity of a consensus against one record, from a full alignment
    of each piece rather than minimap2's local block.

    identity() charges every base minimap2 leaves unaligned at an end as
    a difference. For a sequence close to its record that is nothing -
    the alignment runs end to end - but -x sr stops extending at the
    first divergent stretch, so a sequence 6% from its nearest record
    was charged for a whole clipped end it would have matched most of:
    VESPA's Entamoeba-like candidates read 90% where BLAST gives 95% and
    an end-to-end alignment 94%, enough to fall under the family floor
    and read as unclassified. Here each piece is aligned in full
    (bounded_edit, unit costs):

      "whole"  against the whole record (a primer-to-primer insert);
      "start"  against the record's start, its tail free (R1's window);
      "end"    against the record's end, its head free (R2's window,
               given in the insert's orientation).

    Edits over the pieces' length; None when a piece is further than
    max_frac from the record, where no floor could be cleared anyway.
    """
    edits = length = 0
    for seq, where in parts:
        if not seq:
            continue
        cap = max(1, int(max_frac * len(seq)))
        if where == "whole":
            d = bounded_edit(seq, record, cap)
        elif where == "start":
            d = bounded_edit(seq, record, cap, free_b_suffix=True)
        else:
            d = bounded_edit(seq[::-1], record[::-1], cap, free_b_suffix=True)
        if d > cap:
            return None
        edits += d
        length += len(seq)
    return 1 - edits / length if length else None


def refine_identity(call: dict, parts: list[tuple[str, str]],
                    record: str | None, floors: dict,
                    ranks: list[str] = RANKS) -> dict:
    """Raise a call's identity to its full-alignment value, and re-apply
    the floors from its LCA. Never lowers it: where minimap2 aligned end
    to end the two measures agree, and the floors were calibrated on
    those calls, so they stay exactly as they were."""
    if not record or not call.get("best_hit"):
        return call
    new = aligned_identity(parts, record)
    if new is None or new <= call["identity"] + 1e-9:
        return call
    depth = (apply_floors(call["lca_depth"], new, floors, ranks)
             if call.get("lca_depth", -1) >= 0 else call["depth"])
    return {**call, "identity": new, "depth": depth,
            "rank": ranks[depth] if depth >= 0 else "",
            "identity_minimap2": call["identity"]}


# --------------------------------------------------------------------------
# span scoring and the gap rule (releases scored over their records' span)
# --------------------------------------------------------------------------


def hit_span_identity(h: dict) -> float:
    """minimap2's estimate of a hit's identity over the record's span:
    matches over the alignment plus the record's unaligned ends (the
    query's are free)."""
    return h["match"] / max(1, h["alen"] + h["ts"] + h["tlen"] - h["te"])


def span_candidates(hits: list[dict], n: int = SPAN_CANDIDATES) -> list[str]:
    """The records minimap2 proposes that are scored exactly: the n best
    by estimated identity ("ident", else hit_span_identity), then
    alignment score, and every record tied with the n-th. By identity,
    not score: a partial record scores lower on length alone, and a
    species with hundreds of records pushed a partial of the sister
    species 3 edits off - the one that decides the gap - to 247th by
    score (an H. placei sequence, 7% of a sheep set)."""
    best: dict[str, tuple] = {}
    for h in hits:
        k = (round(h["ident"] if "ident" in h else hit_span_identity(h), 6),
             h["AS"])
        best[h["t"]] = max(best.get(h["t"], k), k)
    order = sorted(best, key=lambda t: (-best[t][0], -best[t][1], t))
    if len(order) > n:
        cut = best[order[n - 1]]
        order = [t for t in order if best[t] >= cut]
    return order


def span_scores(parts: list[tuple[str, str]], cands, refs: dict,
                tail: int = 0) -> list[tuple]:
    """(edits, span, record) for each candidate record, scored over the
    record's span (manifest score_over), each piece in full:

      "whole"  the record placed whole inside the query, the query's
               bases either side free (a consensus, or an insert read
               through: ITS2 plus the 28S stretch);
      "start"  a window against the record's start, the record past the
               window free (R1);
      "end"    a window against the record's end, given in the insert's
               orientation; the window's last `tail` bases (R2's 28S
               stretch) run past the record and are free.

    A window may also overhang a partial record's outer end by up to
    WINDOW_TRIM bases for free; beyond that the overhang is charged.
    Edits and spans add up over the pieces. A record no piece fits is
    left out."""
    out = []
    for c in cands:
        r = refs.get(c)
        if not r:
            continue
        d = sp = 0
        for seq, where in parts:
            if not seq:
                continue
            if where == "whole":
                e, k = infix_distance(r, seq), len(r)
            elif where == "start":
                e, k = anchored_distance(
                    r, seq, min(len(r), max(1, len(seq) - WINDOW_TRIM)))
            else:
                e, k = anchored_distance(
                    r[::-1], seq[::-1],
                    min(len(r), max(1, len(seq) - tail - WINDOW_TRIM)))
            d += e
            sp += k
        if sp and d != float("inf"):
            out.append((int(d), sp, c))
    return out


def span_call(scores: list[tuple], tax, rule: str, floors: dict,
              gap_edits: int = GAP_EDITS, qlen: int = 0,
              min_species_span: int = MIN_SPAN) -> dict:
    """One query's call from its span scores (span_scores), in
    classify_hits()'s shape plus what the gap rule adds.

    Records are ranked by identity over their span, not by raw edits: on
    edits a short partial, with fewer bases to get wrong, beat full-length
    records. The ties are every record at the best identity whose span is
    within SPAN_SLACK of the longest such span - a partial matching only
    over a shorter span does not tie with a full-length match (an H.
    placei sequence, 23% of one cattle set, tied with a partial H.
    contortus record and stopped at genus). The call is their LCA, cut
    back by the floors; under the gap rule a species also needs the
    nearest record of any OTHER species to be gap_edits further off over
    the same span. A shared record (its members name two species) counts
    as another species; an unnamed one (hybrid, cf., environmental) does
    not - counting three hybrid records once cut a cattle set's species
    calls from 87% to 35%. With no other species among the candidates,
    the gap is the whole span.

    Adds: span and edits of the best record, gap (None unless the ties
    name exactly one species), next (record, edits further) for the
    nearest other species, and ties (the tied records)."""
    ranks = ranks_of(tax)
    empty = {"rank": "", "depth": -1, "identity": 0.0, "coverage": 0.0,
             "best_hit": "", "near_hits": 0, "lca_depth": -1, "span": 0,
             "edits": None, "gap": None, "next": None, "ties": []}
    scores = [x for x in scores if x[1] >= MIN_SPAN]
    if not scores:
        return empty
    top = max(1 - d / sp for d, sp, _c in scores)
    ties = [x for x in scores if 1 - x[0] / x[1] >= top - 1e-9]
    longest = max(x[1] for x in ties)
    ties = [x for x in ties if x[1] >= longest - SPAN_SLACK]
    d0, sp0, best = min(ties, key=lambda x: (x[0], -x[1], x[2]))
    species = {tax.get(c, {}).get("species", "") for _d, _s, c in ties}
    other = [(1 - d / sp, sp, c) for d, sp, c in scores
             if tax.get(c, {}).get("species", "") not in species
             and tax.get(c, {}).get("species_status", "shared") != "unnamed"]
    nxt = None
    if other:
        oi, _osp, oc = max(other)
        nxt = (oc, round((1 - oi) * sp0) - d0)
    gap = None
    if len(species) == 1 and "" not in species:
        gap = nxt[1] if nxt else sp0 - d0
    tied = sorted(c for _d, _s, c in ties)
    lca = lca_depth(tied, tax)
    ident = 1 - d0 / sp0
    depth = apply_floors(lca, ident, floors, ranks)
    k = ranks.index("species")
    if depth >= k and (sp0 < min_species_span or rule == "gap" and (
            gap is None or gap < gap_edits)):
        depth = k - 1
    return {"depth": depth, "rank": ranks[depth] if depth >= 0 else "",
            "identity": ident,
            "coverage": min(1.0, sp0 / qlen) if qlen else 1.0,
            "best_hit": best, "near_hits": len(tied), "lca_depth": lca,
            "span": sp0, "edits": d0, "gap": gap, "next": nxt,
            "ties": tied}


def record_label(rid: str, tax) -> str:
    """A record's species; for a shared record the species its members
    name ("Cooperia shared: punctata / spatulata"); else its deepest
    name and rank."""
    t = tax.get(rid, {})
    if t.get("species"):
        return t["species"]
    ranks = ranks_of(tax)
    rank, name = next(((r, t[r]) for r in reversed(ranks) if t.get(r)),
                      ("", ""))
    shared = getattr(tax, "shared", {}).get(rid)
    if shared and name:
        return f"{name} shared: " + " / ".join(
            s.split(" ", 1)[-1] for s in shared)
    return f"{name} ({rank})" if name else rid


def nearest_next(call: dict, tax) -> tuple[str, str]:
    """How close a span-scored call came: (nearest, next).

    nearest: the tied species and their identity ("Haemonchus contortus
    99.6%"); next: the nearest record of any other species and how many
    edits further it is ("Haemonchus placei +2"), or "none close" when
    minimap2 proposed none. A call stopped at genus then says whether it
    was "exactly X, but its sister is 2 edits away" (ITS2 cannot separate
    them) or "98% to X, nothing else close" (a variant the database does
    not hold, or a missing relative)."""
    if not call.get("span"):
        return "", ""
    shared = getattr(tax, "shared", {})
    names = sorted({n for t in call["ties"]
                    for n in (shared.get(t) or [record_label(t, tax)])})
    near = " / ".join(names[:3]) + (
        f" / +{len(names) - 3} more" if len(names) > 3 else "")
    nxt = call.get("next")
    return (f"{near} {100 * call['identity']:.1f}%",
            f"{record_label(nxt[0], tax)} +{nxt[1]}" if nxt else "none close")


def span_columns(call: dict, tax) -> dict:
    """The columns a span-scored call adds to clusters.tsv."""
    near, nxt = nearest_next(call, tax)
    gap = call.get("gap")
    return {"span": call.get("span") or "",
            "edits": "" if call.get("edits") is None else call["edits"],
            "gap": "" if gap is None else gap,
            "nearest": near, "next": nxt,
            "grade": tax.get(call.get("best_hit", ""), {}).get("grade", "")}


SPAN_COLUMNS = ["span", "edits", "gap", "nearest", "next", "grade"]


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
            return taxon, taxon
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


def parse_floors(text: str, ranks: list[str] = RANKS) -> dict:
    out = {}
    for part in (text or "").split(","):
        if not part.strip():
            continue
        k, _, v = part.partition("=")
        k = k.strip()
        if k not in ranks:
            raise SystemExit(f"--floors: '{k}' is not a rank ({', '.join(ranks)})")
        out[k] = float(v)
    return out


def add_calling_args(ap) -> None:
    """The calling-rule options both routes share."""
    ap.add_argument("--calling-rule", choices=("auto", "floors", "gap"),
                    default="auto",
                    help="floors: an identity floor per rank (--floors); gap: "
                    "species only when the next species is --gap-edits "
                    "further off (--gap-floors). auto: the release's "
                    "calling_rule, else gap for a release scored over its "
                    "records' span (ITS2), else floors (V4)")
    ap.add_argument("--gap-floors", default=",".join(
        f"{k}={v}" for k, v in GAP_FLOORS.items()),
        help="identity floor per rank under the gap rule")
    ap.add_argument("--gap-edits", type=int, default=GAP_EDITS,
                    help="edits the nearest other species must be further "
                    "off for a species call under the gap rule")
    ap.add_argument("--min-species-span", type=int, default=MIN_SPAN,
                    help="record bases a span-scored species call must rest "
                    f"on (default {MIN_SPAN}, the least any record must "
                    "cover to count)")


def calling_setup(args, man: dict, tax) -> tuple:
    """(exact, rule, floors): whether calls come from span scoring, the
    calling rule, and its floors over the release's ranks."""
    rule = calling_rule(man, args.calling_rule)
    exact = scores_over_reference(man) or rule == "gap"
    floors = parse_floors(args.gap_floors if rule == "gap" else args.floors,
                          ranks_of(tax))
    return exact, rule, floors


def calling_settings(exact: bool, rule: str, args) -> dict:
    """What the summary records of how calls were made. Nothing for a V4
    release under its floors, so those summaries read as they did."""
    if not exact:
        return {}
    return {"score_over": "reference span", "calling_rule": rule,
            "gap_edits": args.gap_edits if rule == "gap" else None,
            "min_span": MIN_SPAN, "min_species_span": args.min_species_span,
            "span_slack": SPAN_SLACK}


def trim_reads(path: str, primers) -> tuple[dict[str, str], Counter]:
    """Oriented, primer-trimmed inserts by read id, and what was lost.

    The release's records are primer-trimmed inserts, so reads must be
    too: identity is measured over the insert, and a primer left on one
    side and not the other would be charged as a difference."""
    fate: Counter = Counter()
    return dict(iter_trimmed(path, primers, fate)), fate


def iter_trimmed(path: str, primers, fate: Counter):
    """trim_reads() one read at a time, so main() can write them to disk
    without holding a deep well's reads in memory; fate is filled in as
    they pass. Only the ids are kept, to make duplicates unique."""
    out: set[str] = set()
    from ampseqlib import revcomp
    for rid, seq, _q in open_reads(path):
        fate["reads_in"] += 1
        seq = seq.upper()
        hit = None
        for s in (seq, revcomp(seq)):
            f = primers.match(s, "fwd")
            r = primers.match(revcomp(s), "rev")
            if f is not None and r is not None:
                hit = s[f[1]:len(s) - r[1]]
                break
        if hit is None:
            fate["no_primer"] += 1
            continue
        if len(hit) < 50:
            fate["too_short"] += 1
            continue
        # Duplicate ids happen when a split read keeps its parent's name;
        # minimap2 would merge them, so each gets its own.
        key = rid if rid not in out else f"{rid}_{fate['reads_in']}"
        out.add(key)
        yield key, hit


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", required=True, help="the sample's FASTQ")
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
    exact, rule, floors = calling_setup(args, man, tax)
    st = Settings(min_reads=args.min_cluster_reads,
                  max_reads=args.max_consensus_reads,
                  hap_min_frac=args.haplotype_min_frac,
                  merge_edits=args.merge_edits, split=not args.no_split)
    mm = Minimap2(args.minimap2, args.threads)

    primers = load_primers(args.primers)
    primers.max_edits = args.primer_max_edits
    fate: Counter = Counter()
    tmp = tempfile.mkdtemp(prefix=f"assign_db_{sid}_", dir=args.outdir)
    clusters: list[dict] = []
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
        for t in sorted(heads, key=lambda t: (-len(groups[t]), t)):
            qs = sorted(groups[t])
            if len(qs) > st.max_reads:
                qs = sorted(rng.sample(qs, st.max_reads))
            chosen.append((t, qs))
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
                built.append(c)
        clusters = merge_close(built, st.merge_edits)

        # Count every read against the final set.
        # Each read is decided on its first primary hit as it streams
        # past; only the ids already decided are kept.
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
                if own is not None and own - a["AS"] > args.recruit_margin:
                    continue
                counts[int(a["t"][1:])] += 1
        n_clustered = sum(counts.values())
        fate["clustered"] = n_clustered
        fate["unclustered"] = n_reads - n_clustered

        # Classify the consensus sequences that kept any reads.
        keep = [k for k in range(len(clusters)) if counts[k] > 0]
        calls: dict[int, dict] = {}
        if keep:
            kfa = os.path.join(tmp, "classify.fa")
            write_fasta(kfa, ((f"c{k}", clusters[k]["consensus"]) for k in keep))
            hits: dict[str, list] = defaultdict(list)
            for a in mm.map(index, kfa, preset=INDEX_ARGS,
                            extra=("--secondary=yes", "-N", "1000",
                                   "-p", "0.9")):
                hits[a["q"]].append(a)
            if exact:
                # Scored over each record's span (span_call): the
                # consensus carries the 28S stretch its records lack.
                cands = {k: span_candidates(hits.get(f"c{k}", []))
                         for k in keep}
                need = {t for ts in cands.values() for t in ts}
                recs = {i: s for i, s in read_fasta(db_fasta) if i in need}
                for k in keep:
                    cons = clusters[k]["consensus"]
                    calls[k] = span_call(
                        span_scores([(cons, "whole")], cands[k], recs),
                        tax, rule, floors, args.gap_edits, len(cons),
                        args.min_species_span)
            else:
                for k in keep:
                    calls[k] = classify_hits(hits.get(f"c{k}", []), tax,
                                             floors, args.lca_margin)
                # Identity from a full alignment to the best record, where
                # that is higher than minimap2's (aligned_identity).
                need = {c["best_hit"] for c in calls.values() if c["best_hit"]}
                recs = {i: s for i, s in read_fasta(db_fasta) if i in need}
                for k in keep:
                    calls[k] = refine_identity(
                        calls[k], [(clusters[k]["consensus"], "whole")],
                        recs.get(calls[k]["best_hit"]), floors, ranks)
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
            **(span_columns(call, tax) if exact else {}),
        })
    cols = ["sample_id", "cluster", "reads", "rel_abundance", "rank", "taxon",
            "label", "identity", "coverage", "lineage", "best_hit",
            "best_hit_lineage", "near_hits", "lca_rank", "db_resolution",
            "db_margin", "origin", "merged", "chimera_of", "length",
            "database_version", "consensus"] + (SPAN_COLUMNS if exact else [])
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

    # Per taxon at the rank each call reached. Two clusters of one
    # species (haplotypes) are one row; a species call and a genus call
    # in the same genus are NOT merged - they are different statements.
    agg: dict = {}
    for r in rows:
        key = (r["rank"], r["taxon"], r["label"])
        a = agg.setdefault(key, {"reads": 0, "clusters": 0, "identity": 0.0,
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
        "sample_id": sid,
        "mode": "database",
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
            "over clustered reads: relative 18S abundance, never organisms "
            "or copies."),
        "clusters": len(rows),
        "groups": len(heads),
        "haplotype_splits": cb.splits,
        "chimeras": {"sequences": len(chim),
                     "reads": sum(counts[k] for k in chim),
                     "within_genus_kept": {
                         "sequences": len(kept_joins),
                         "reads": sum(counts[k] for k in kept_joins)}},
        "minimap2_calls": mm.calls,
        "settings": {**st.as_dict(), "floors": floors,
                     "lca_margin": args.lca_margin,
                     "recruit_margin": args.recruit_margin,
                     "index_args": INDEX_ARGS,
                     **calling_settings(exact, rule, args)},
        "rank_reached": dict(Counter(r["rank"] or "none" for r in rows)),
    }
    with open(os.path.join(args.outdir, f"{sid}.assign_db.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"{sid}: {fate['reads_in']} reads, {len(rows)} clusters, "
          f"{fate['clustered']} clustered, {fate['unclustered']} unclustered, "
          f"{cb.splits} haplotype split(s), {summary['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
