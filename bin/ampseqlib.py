#!/usr/bin/env python3
"""
Shared primitives for amphid.

Everything in here is used by at least two of the executables, and the
reason it is shared rather than copied is that the alternative has already
bitten this pipeline once: the approximate matcher had a one-character bug
that silently halved assignment, and a second copy of it would have had to
be found and fixed separately.

Standard library only, so the field profile has nothing to install.
"""

from __future__ import annotations

import csv
import glob
import gzip
import json
import math
import os
import random
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import compress
from operator import sub

# --------------------------------------------------------------------------
# construct defaults, overridden by the assay config
# --------------------------------------------------------------------------

DEFAULTS = {
    "p5_tail": "AATGATACGGCGACCACCGAGATCTACAC",
    "p7_tail": "CAAGCAGAAGACGGCATACGAGAT",
    "i5_inner_anchor": "TCGTCGGCAGCGTC",
    "i7_inner_anchor": "CCGAGCCCACGAGAC",
    # The anchors are only the DISCRIMINATING part of each overhang - the
    # part that differs between the two ends. Both overhangs also carry the
    # 19bp Nextera mosaic end, which is shared and so useless for telling the
    # ends apart, but is still adapter.
    "pcr1_overhang_fwd": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "pcr1_overhang_rev": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
}


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


# IUPAC ambiguity, as the set of concrete bases each code stands for. Real
# primers are routinely degenerate, and a matcher that treats Y as a literal
# character fails to find the primer in every sequence that carries either
# base it was designed to cover.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


# --------------------------------------------------------------------------
# approximate matching
# --------------------------------------------------------------------------


def find_anchor_cost(
    text: str, pattern: str, max_edits: int
) -> tuple[int, int, int] | None:
    """As find_anchor, but also returns the edit cost.

    Primer trimming needs the cost: several primers compete for the same
    read start, and the one that fits best is the one that was used.
    """
    m = len(pattern)
    if not m or m > 64 or not text:
        return None

    # Bits per concrete base, expanding any degenerate code in the pattern,
    # then one mask per IUPAC code the text might use - so degeneracy on
    # either side matches, and an N in either matches anything.
    base: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i, c in enumerate(pattern):
        for b in IUPAC.get(c, ""):
            base[b] |= 1 << i
    peq = {
        code: (base["A"] if "A" in bases else 0)
        | (base["C"] if "C" in bases else 0)
        | (base["G"] if "G" in bases else 0)
        | (base["T"] if "T" in bases else 0)
        for code, bases in IUPAC.items()
    }

    full = (1 << m) - 1
    top = 1 << (m - 1)
    pv, mv, score = full, 0, m
    best_end, best_cost = -1, max_edits + 1

    for j, c in enumerate(text, 1):
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) ^ pv) | eq) & full
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = (mh | ~(xv | ph)) & full
        mv = ph & xv
        if score < best_cost:
            best_cost, best_end = score, j

    if best_end < 0 or best_cost > max_edits:
        return None
    return max(0, best_end - m), best_end, best_cost


def find_anchor_leftmost(
    text: str, pattern: str, max_edits: int
) -> tuple[int, int, int] | None:
    """The FIRST occurrence of `pattern` within tolerance, not the best.

    find_anchor_cost answers "where does this pattern fit best in this
    window", which is the right question when the window is known to
    hold one copy - a primer at a read start, an anchor at a read end.
    It is the wrong question for a caller that wants every copy, left to
    right. Such a caller searches the remainder after each hit, and if
    the remainder holds a one-edit copy followed by an exact one, the
    best hit is the exact one: the search jumps straight past the
    earlier copy and never comes back for it. In the splitter that loses
    exactly the anchor that marks a junction, and the read goes on to be
    demultiplexed on one molecule's i5 and another's i7.

    So: the first end position whose cost is within tolerance. Not that
    position alone, though - an exact copy ending at j also matches
    with k deletions ending at j - k, so the first position within
    tolerance is usually a worse-scoring early view of the same copy.
    The scan carries on for max_edits further columns and takes the
    lowest cost in that window (the earliest, on a tie), which is the
    end of the copy it found rather than a ragged prefix of it. Beyond
    that window any further hit is a different copy, and is left for
    the caller's next search.

    Returns (start, end, cost) like find_anchor_cost, with start taken
    the same way - end minus pattern length - so the two agree on any
    window holding a single copy.
    """
    m = len(pattern)
    if not m or m > 64 or not text:
        return None

    base: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i, c in enumerate(pattern):
        for b in IUPAC.get(c, ""):
            base[b] |= 1 << i
    peq = {
        code: (base["A"] if "A" in bases else 0)
        | (base["C"] if "C" in bases else 0)
        | (base["G"] if "G" in bases else 0)
        | (base["T"] if "T" in bases else 0)
        for code, bases in IUPAC.items()
    }

    full = (1 << m) - 1
    top = 1 << (m - 1)
    pv, mv, score = full, 0, m
    best_end, best_cost, stop_at = -1, max_edits + 1, -1

    for j, c in enumerate(text, 1):
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) ^ pv) | eq) & full
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = (mh | ~(xv | ph)) & full
        mv = ph & xv
        if score < best_cost:
            best_cost, best_end = score, j
            if stop_at < 0:
                stop_at = j + max_edits
        if 0 <= stop_at <= j:
            break

    if best_end < 0 or best_cost > max_edits:
        return None
    return max(0, best_end - m), best_end, best_cost


def find_anchor(text: str, pattern: str, max_edits: int) -> tuple[int, int] | None:
    """Locate `pattern` in `text` allowing indels; return (start, end) of the
    best alignment, or None.

    Edit distance rather than Hamming: a fixed-width comparison fails
    outright on an indel, and indels are nanopore's dominant error - roughly
    one read in eight is lost on a 14bp anchor that way.

    Myers' bit-parallel algorithm rather than the obvious dynamic program.
    The DP is O(window x pattern) per call and, at four calls per read across
    a hundred thousand reads, runs to hundreds of millions of Python-level
    operations. Myers is O(window) integer operations for any pattern up to
    64bp, which every anchor and primer here is.

    Returning the END of the alignment matters: an indel inside the pattern
    shifts everything after it, so a downstream offset must be measured from
    where the match actually finished, not from its nominal start plus its
    length.
    """
    hit = find_anchor_cost(text, pattern, max_edits)
    return None if hit is None else (hit[0], hit[1])


# --------------------------------------------------------------------------
# seed-and-verify: the same answers as the matchers above, sooner
# --------------------------------------------------------------------------
#
# The splitter runs find_anchor_leftmost four times across every read, and
# the look-alike mask runs find_anchor_cost once per known look-alike at
# every hit. Both are pure-Python loops over every base, and on a real
# MinION run together they were two-thirds of the pipeline's runtime -
# about 35 hours on a laptop for one flow cell.
#
# Almost all of those bases are nowhere near a match, and that can be
# known without aligning them. Cut the pattern into max_edits + 1 pieces:
# an alignment with at most max_edits edits cannot touch every piece, so
# at least one piece appears in the text EXACTLY (the pigeonhole
# argument). str.find locates pieces at C speed, and every end position
# whose cost could be within tolerance lies within max_edits of where
# that piece's alignment would put the pattern's end. Only those short
# stretches are handed to the bit-parallel alignment.
#
# The answers are not approximations of the old ones; they are the old
# ones. An alignment of an m-base pattern ending at j with cost c <= k
# starts no earlier than j - m - k, so running Myers from there (or from
# the caller's own left edge, if that is later) gives every within-
# tolerance cost exactly, and a cost above tolerance stays above it.
# Only ends whose cost is within tolerance ever decide a result, so
# visiting nothing else changes nothing.
#
# The argument needs exact matching to mean equality, which IUPAC codes
# break: an N in a read matches anything. A pattern or text that holds
# any code other than A, C, G and T goes the old way. Basecallers emit
# neither, so that path is for completeness, not for speed.

_AMBIGUOUS_CODES = tuple(c for c in IUPAC if c not in "ACGT")


def is_plain_acgt(text: str) -> bool:
    """True when no character of `text` is an IUPAC code other than A, C,
    G or T - so an exact substring match is the only kind possible.

    Characters outside IUPAC altogether (lower case, gaps) match nothing
    in the aligners either, so they do not disqualify a text.
    """
    return not any(c in text for c in _AMBIGUOUS_CODES)


def seed_pieces(pattern: str, max_edits: int) -> list[tuple[int, str]] | None:
    """`pattern` cut into max_edits + 1 contiguous pieces, as (offset,
    piece), or None where the pigeonhole argument does not hold.
    """
    m, n = len(pattern), max_edits + 1
    if max_edits < 0 or m < n or set(pattern) - set("ACGT"):
        return None
    cuts = [m * i // n for i in range(n + 1)]
    return [(cuts[i], pattern[cuts[i]:cuts[i + 1]]) for i in range(n)]


def _pattern_masks(pattern: str) -> dict[str, int]:
    # Built exactly as find_anchor_cost builds its own.
    base: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i, c in enumerate(pattern):
        for b in IUPAC.get(c, ""):
            base[b] |= 1 << i
    return {
        code: (base["A"] if "A" in bases else 0)
        | (base["C"] if "C" in bases else 0)
        | (base["G"] if "G" in bases else 0)
        | (base["T"] if "T" in bases else 0)
        for code, bases in IUPAC.items()
    }


class SeededAnchorScan:
    """One pattern, searched seed-and-verify.

    scan() gives what calling find_anchor_leftmost on each remainder in
    turn gives - the splitter's anchor scan - and best() gives what
    find_anchor_cost gives - the look-alike mask's test. Neither aligns
    more of the text than the seeds point at.
    """

    def __init__(self, pattern: str, max_edits: int):
        self.pattern, self.k, self.m = pattern, max_edits, len(pattern)
        self.seeds = seed_pieces(pattern, max_edits) if self.m <= 64 else None
        self.peq = _pattern_masks(pattern)

    def candidate_ends(self, text: str) -> list[list[int]]:
        """Sorted, disjoint [lo, hi] ranges of 1-based end positions
        holding every end whose cost could be within tolerance."""
        m, k, n = self.m, self.k, len(text)
        spans = []
        for a, piece in self.seeds:
            shift = m - a
            i = text.find(piece)
            while i >= 0:
                spans.append((i + shift - k, i + shift + k))
                i = text.find(piece, i + 1)
        spans.sort()
        merged: list[list[int]] = []
        for lo, hi in spans:
            lo, hi = max(lo, 1), min(hi, n)
            if lo > hi:
                continue
            if merged and lo <= merged[-1][1] + 1:
                if hi > merged[-1][1]:
                    merged[-1][1] = hi
            else:
                merged.append([lo, hi])
        return merged

    def _costs(self, text: str, w: int, stop: int) -> list[int]:
        """Myers from a fresh start at text[w], as find_anchor_leftmost
        runs it; the cost at each end j in w+1..stop, at index j - w - 1.
        """
        peq, m = self.peq, self.m
        full = (1 << m) - 1
        top = 1 << (m - 1)
        pv, mv, score = full, 0, m
        out = []
        for c in text[w:stop]:
            eq = peq.get(c, 0)
            xv = eq | mv
            xh = ((((eq & pv) + pv) ^ pv) | eq) & full
            ph = (mv | ~(xh | pv)) & full
            mh = pv & xh
            if ph & top:
                score += 1
            elif mh & top:
                score -= 1
            ph = (ph << 1) & full
            mh = (mh << 1) & full
            pv = (mh | ~(xv | ph)) & full
            mv = ph & xv
            out.append(score)
        return out

    def best(self, text: str,
             plain: bool | None = None) -> tuple[int, int, int] | None:
        """What find_anchor_cost(text, pattern, k) returns: the lowest-cost
        end within tolerance, the earliest on a tie, as (start, end, cost).
        """
        m, k = self.m, self.k
        if plain is None:
            plain = is_plain_acgt(text)
        if self.seeds is None or not plain:
            return find_anchor_cost(text, self.pattern, k)
        best_end, best_cost = -1, k + 1
        for lo, hi in self.candidate_ends(text):
            w = max(0, lo - m - k)
            costs = self._costs(text, w, hi)
            for j in range(lo, hi + 1):
                if costs[j - w - 1] < best_cost:
                    best_end, best_cost = j, costs[j - w - 1]
        if best_end < 0:
            return None
        return max(0, best_end - m), best_end, best_cost

    def scan(self, text: str, limit: int = 24,
             plain: bool | None = None) -> list[tuple[int, int]]:
        """(start, end) of every non-overlapping occurrence, left to right,
        exactly as the loop

            off = 0
            while len(hits) < limit and off < len(text):
                s, e, _ = find_anchor_leftmost(text[off:], pattern, k)
                hits.append((off + s, off + e)); off += max(e, s + 1)

        would give. `plain` is is_plain_acgt(text), if already known.
        """
        m, k, n = self.m, self.k, len(text)
        if plain is None:
            plain = is_plain_acgt(text)
        if self.seeds is None or not plain:
            hits, off = [], 0
            while len(hits) < limit and off < n:
                h = find_anchor_leftmost(text[off:], self.pattern, k)
                if h is None:
                    break
                hits.append((off + h[0], off + h[1]))
                off += max(h[1], h[0] + 1)
            return hits

        ranges = self.candidate_ends(text)
        hits: list[tuple[int, int]] = []
        off, r = 0, 0
        while len(hits) < limit and off < n and r < len(ranges):
            lo, hi = ranges[r]
            if hi <= off:
                r += 1
                continue
            # find_anchor_leftmost(text[off:]) sees no base before off, so
            # neither may the alignment here.
            lo = max(lo, off + 1)
            w = max(off, lo - m - k)
            costs = self._costs(text, w, min(n, hi + k))
            first = -1
            for j in range(lo, hi + 1):
                if costs[j - w - 1] <= k:
                    first = j
                    break
            if first < 0:
                r += 1
                continue
            # find_anchor_leftmost's settling window: the lowest cost in
            # the next k columns, earliest on a tie. Columns outside every
            # candidate range cost more than k and never win.
            best_end, best_cost = first, costs[first - w - 1]
            for j in range(first + 1, min(first + k, n) + 1):
                c = costs[j - w - 1]
                if c < best_cost:
                    best_end, best_cost = j, c
            start = max(off, best_end - m)
            hits.append((start, best_end))
            off = max(best_end, start + 1)
            # The same range may hold the next copy; look again from off.
        return hits


# --------------------------------------------------------------------------
# primers
# --------------------------------------------------------------------------


@dataclass
class PrimerSet:
    """The primers, kept as flat lists because a read does not know its locus
    yet - that is what assignment is for. Each read start is tried against
    every primer of the right direction and the best fit wins."""

    fwd: list[tuple[str, str]] = field(default_factory=list)
    rev: list[tuple[str, str]] = field(default_factory=list)
    max_edits: int = 2

    def match(self, seq: str, direction: str) -> tuple[str, int] | None:
        """Best primer at the start of `seq`; returns (name, end offset).

        The search window is the primer length plus a little slack, not the
        whole read: a primer sequence occurring by chance in the middle of an
        amplicon is not evidence of anything, and allowing the match to float
        would let it trim away real template.
        """
        best = None
        for name, p in self.fwd if direction == "fwd" else self.rev:
            window = seq[: len(p) + 6]
            hit = find_anchor_cost(window, p, self.max_edits)
            if hit is None:
                continue
            start, end, cost = hit
            # A primer that does not start at (or very near) the read start
            # is a coincidence, not the primer.
            if start > 4:
                continue
            if best is None or cost < best[2]:
                best = (name, end, cost)
        return None if best is None else (best[0], best[1])


def load_primers(path: str) -> PrimerSet:
    ps = PrimerSet()
    for r in read_table(path):
        d = r.get("direction", "").strip().lower()
        seq = r.get("sequence", "").strip().upper()
        if not seq:
            continue
        if d.startswith("f"):
            ps.fwd.append((r.get("name", "?"), seq))
        elif d.startswith("r"):
            ps.rev.append((r.get("name", "?"), seq))
    return ps


# --------------------------------------------------------------------------
# bounded edit distance
# --------------------------------------------------------------------------


def _strip_affixes(a: str, b: str) -> tuple[str, str]:
    """Remove the shared prefix and suffix.

    Safe for Levenshtein - an optimal alignment can always be arranged to
    match identical flanks - and it is what makes competitive scoring
    affordable on Illumina data, where a read differs from its reference in
    a handful of places and is identical either side of them.
    """
    i, n = 0, min(len(a), len(b))
    while i < n and a[i] == b[i]:
        i += 1
    a, b = a[i:], b[i:]
    j, n = 0, min(len(a), len(b))
    while j < n and a[len(a) - 1 - j] == b[len(b) - 1 - j]:
        j += 1
    return (a[: len(a) - j], b[: len(b) - j]) if j else (a, b)


def _bounded_edit_dp(a: str, b: str, max_d: int, free_b_suffix: bool = False) -> int:
    """Levenshtein distance between a and b, or max_d + 1 if it exceeds it.

    The original banded dynamic program, kept as the reference that
    bounded_edit is proven against (tests/test_assign.py) and as the path
    for a negative max_d, which nothing in the pipeline asks for. Its
    return values ARE the contract; bounded_edit only computes them faster.

    Banded: only cells within max_d of the diagonal are computed, because a
    path that strays further has already spent more than the budget. That
    turns an O(n*m) comparison into O(n*max_d), which is what makes it
    usable against a whole panel.

    With free_b_suffix, b may run past the end of a at no cost - the mode
    needed when a read covers only the start of a reference, as the mates of
    an unmerged pair do. Without it a 250bp mate would be charged for every
    base of the amplicon it never reached.
    """
    if not free_b_suffix:
        a, b = _strip_affixes(a, b)
    n, m = len(a), len(b)
    if n == 0:
        return m if (m <= max_d or free_b_suffix) else max_d + 1
    if m == 0:
        return n if n <= max_d else max_d + 1
    if not free_b_suffix and abs(n - m) > max_d:
        return max_d + 1

    inf = max_d + 1
    prev = [j if j <= max_d else inf for j in range(m + 1)]
    for i in range(1, n + 1):
        lo, hi = max(1, i - max_d), min(m, i + max_d)
        cur = [inf] * (m + 1)
        if i <= max_d:
            cur[0] = i
        best = inf
        ai = a[i - 1]
        for j in range(lo, hi + 1):
            c = prev[j - 1] + (ai != b[j - 1])
            d = prev[j] + 1
            if d < c:
                c = d
            d = cur[j - 1] + 1
            if d < c:
                c = d
            cur[j] = c
            if c < best:
                best = c
        if best > max_d:
            return max_d + 1
        prev = cur

    out = min(prev[1:]) if free_b_suffix else prev[m]
    return out if out <= max_d else max_d + 1



def bounded_edit(a: str, b: str, max_d: int, free_b_suffix: bool = False) -> int:
    """Levenshtein distance between a and b, or max_d + 1 if it exceeds it.

    With free_b_suffix, b may run past the end of a at no cost - the mode
    needed when a read covers only the start of a reference, as the mates of
    an unmerged pair do. Without it a 250bp mate would be charged for every
    base of the amplicon it never reached.

    The result is exactly that of _bounded_edit_dp, the banded dynamic
    program this replaced, edge cases included; only the arithmetic differs.
    That program spent 90% of assignment in its inner loop - about 5 ms a
    comparison on a 400bp nanopore amplicon, ~12 h of a laptop field run
    on a full flow cell - because it touched every cell of a 400 x (2*max_d
    + 1) band one Python operation at a time. Here a whole column of the
    alignment matrix is one integer (Myers 1999, in Hyyro's form for global
    distance), so a comparison is ~15 big-integer operations per base of b
    instead of ~10 Python operations per cell.

    Why the two agree, so that nobody has to take the tests on trust:

    - Without free_b_suffix the band changes nothing that is returned. A
      path leaving the band has spent more than max_d by then, the band
      fills the cells outside it with max_d + 1 (a lower bound on their true
      value), and every number that grows from them is already over the
      limit. So the DP returns min(distance, max_d + 1) - which is what the
      full, unbanded matrix gives when capped.
    - With free_b_suffix it returns min over j >= 1 of D[len(a)][j], capped,
      by the same argument. j starts at 1, not 0: the DP never looks at
      D[len(a)][0], and nor does this. (It only matters for an empty a, and
      that case returns before either computation.)
    - Characters are compared with plain equality in both: an N matches an
      N and nothing else, as do IUPAC codes and lower case. The match masks
      are built from the characters of a as they are, so no alphabet is
      assumed.
    - The empty-string and length-difference returns are the DP's own
      lines, run before either computation, including its one uncapped
      value: an empty a in free_b_suffix mode returns len(b) whatever the
      budget.
    """
    if max_d < 0:
        return _bounded_edit_dp(a, b, max_d, free_b_suffix)
    if not free_b_suffix:
        a, b = _strip_affixes(a, b)
    n, m = len(a), len(b)
    if n == 0:
        return m if (m <= max_d or free_b_suffix) else max_d + 1
    if m == 0:
        return n if n <= max_d else max_d + 1
    if not free_b_suffix and abs(n - m) > max_d:
        return max_d + 1

    # a down the rows, b across the columns, because the free end belongs
    # to b: the bottom row D[n][j] is then read off one column at a time,
    # as the score, and its minimum is the free-suffix answer.
    #
    # Bit i of pv / mv says that D[i+1][j] - D[i][j] is +1 / -1 in the
    # current column (neither: 0). Column 0 is D[i][0] = i, all +1.
    # peq[c] marks the rows whose base equals c.
    peq: dict[str, int] = {}
    bit = 1
    for c in a:
        peq[c] = peq.get(c, 0) | bit
        bit <<= 1
    mask = bit - 1
    top = bit >> 1
    # Every base of b gets an entry, 0 where it occurs nowhere in a, so
    # the loop below is a plain list walk with no lookups in it.
    for c in set(b).difference(peq):
        peq[c] = 0
    eqs = list(map(peq.__getitem__, b))
    pv, mv, score = mask, 0, n
    inf = max_d + 1

    # One step per column j of the matrix. ph / mh mark the +1 / -1
    # horizontal steps D[i][j] - D[i][j-1]; the top row's is the change in
    # the score. ph is left unmasked (it may go negative): every operation
    # here carries only upwards, so bits above the top row never reach the
    # ones that are read, and the one mask on pv keeps the state bounded.
    # The "| 1" is the global-distance boundary: D[0][j] = j, so the step
    # along row 0 is always +1. The two modes differ only in what they
    # read from the bottom row, and are written out twice to keep that
    # test out of the loop.
    if free_b_suffix:
        best = inf
        for eq in eqs:
            xv = eq | mv
            xh = (((eq & pv) + pv) ^ pv) | eq
            ph = mv | ~(xh | pv)
            mh = pv & xh
            if ph & top:
                score += 1
            elif mh & top:
                score -= 1
            if score < best:
                best = score
            ph = (ph << 1) | 1
            pv = ((mh << 1) | ~(xv | ph)) & mask
            mv = ph & xv
        return best

    # D[n][m] is at least D[n][j] - (m - j): the bottom row can fall by at
    # most one per column. Once even that cannot bring it back within
    # max_d, stop. The DP would have stopped too or finished over the
    # limit; either way the answer is max_d + 1.
    lim = max_d + m
    for j, eq in enumerate(eqs, 1):
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh
        if ph & top:
            score += 1
            if score + j > lim:
                return inf
        elif mh & top:
            score -= 1
        ph = (ph << 1) | 1
        pv = ((mh << 1) | ~(xv | ph)) & mask
        mv = ph & xv
    return score if score <= max_d else inf

# --------------------------------------------------------------------------
# io
# --------------------------------------------------------------------------


def writer(fh, fieldnames=None, delimiter: str = "\t"):
    """A csv writer that ends lines the way Unix does.

    csv defaults to CRLF regardless of platform. Every table here is read
    by awk, cut, pandas and shell loops, and a trailing carriage return
    turns the last column of every row into a value with an invisible
    character on the end - which surfaces as a file that cannot be opened
    and a path that looks exactly right in the error message.
    """
    if fieldnames is None:
        return csv.writer(fh, delimiter=delimiter, lineterminator="\n")
    return csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter,
                          lineterminator="\n")


# --------------------------------------------------------------------------
# anchor look-alikes inside the panel's own amplicons
# --------------------------------------------------------------------------

# Sequence kept either side of a look-alike to recognise it again in a
# read. 24 + a 14-15bp anchor stays inside the 64bp the bit-parallel
# matcher handles, and is long enough that no adapter junction resembles
# it: a junction's flanks are index and Nextera overhang, not insert.
LOOKALIKE_FLANK = 24

# Found at one edit looser than the splitter and demultiplexer search at,
# so a look-alike that a sequencing error carries inside their tolerance
# is still in the mask. Masking a site that never triggers costs nothing:
# it is insert either way.
LOOKALIKE_MAX_EDITS = 3


def find_anchor_lookalikes(entries, anchors: dict, max_edits: int,
                           flank: int = LOOKALIKE_FLANK) -> list[dict]:
    """Every place an adapter anchor is nearly spelled out inside an amplicon.

    The splitter cuts a read wherever it finds an anchor, and the chimera
    check quarantines a read with an anchor in its interior - both on the
    assumption that an amplicon's insert contains no adapter. A panel
    entry whose insert happens to contain a near-copy of an anchor breaks
    that assumption for every one of its reads: each is cut in two, or set
    aside as a chimera, and the target is quietly under-counted. Found on
    a real panel, where a Trichuris incognita amplicon carried the i5
    anchor at two edits and lost most of its reads.

    entries: (seq_id, sequence) pairs; anchors: name -> sequence. Both
    strands are searched, since reads arrive in either orientation.
    """
    rows = []
    for seq_id, seq in entries:
        seq = (seq or "").upper()
        for name, anchor in anchors.items():
            if not anchor:
                continue
            for strand, s in (("+", seq), ("-", revcomp(seq))):
                off = 0
                while off < len(s):
                    h = find_anchor_leftmost(s[off:], anchor, max_edits)
                    if h is None:
                        break
                    st, en, cost = h[0] + off, h[1] + off, h[2]
                    rows.append({
                        "seq_id": seq_id, "anchor": name, "strand": strand,
                        "position": st, "edits": cost,
                        "context": s[max(0, st - flank):en + flank],
                    })
                    off = max(en, st + 1)
    return rows


class LookalikeMask:
    """Recognise a panel look-alike in a read, so it is not taken for an
    adapter.

    Built from pre-flight's anchor_lookalikes.tsv. An anchor hit in a read
    is a look-alike when the sequence around it matches a known
    look-alike's surroundings - insert on both sides - rather than the
    index and overhang that flank a real junction. Only those hits are
    ignored; every other anchor is treated exactly as before, so a real
    junction, including one with no index, is still cut and still
    quarantined.
    """

    def __init__(self, contexts=()):
        self.contexts: list[str] = []
        # The membership test is a set, not the list. With a panel's dozen
        # contexts the list was fine; a database release lists ~46,000
        # rows, and testing each against a growing list took six seconds
        # in every splitting shard before a single read was read. The
        # order kept - and so every search below - is unchanged.
        seen: set[str] = set()
        for c in contexts:
            c = (c or "").upper()
            if c and c not in seen:
                seen.add(c)
                self.contexts.append(c)
                rc = revcomp(c)
                if rc not in seen:
                    seen.add(rc)
                    self.contexts.append(rc)
        # One seeded search per context, at the tolerance covers() uses.
        self.searches = [SeededAnchorScan(c, max(4, len(c) // 8))
                         for c in self.contexts]
        self.widest = max(map(len, self.contexts), default=0)
        self._index = (_ContextIndex(self.contexts)
                       if len(self.contexts) > LOOKALIKE_INDEX_MIN else None)

    @classmethod
    def load(cls, path: str | None) -> "LookalikeMask":
        if not path or not os.path.exists(path):
            return cls()
        return cls(r.get("context", "") for r in read_table(path))

    def __bool__(self) -> bool:
        return bool(self.contexts)

    def covers(self, read: str, start: int, end: int) -> bool:
        """Does a known look-alike's context sit over read[start:end]?

        The context is sought in a stretch of read around the hit, at a
        nanopore-sized error budget (one edit in eight); a match counts
        when it overlaps the hit.

        Most hits are real adapters, so most contexts are absent, and each
        used to cost a full alignment of the stretch around the hit. The
        alignment is now seed-and-verify (SeededAnchorScan.best): a stretch
        holding no exact piece of a context cannot match it and is passed
        over, and one that does is aligned only near those pieces. The
        answer is find_anchor_cost's, unchanged.
        """
        widest = self.widest
        plain = is_plain_acgt(read[max(0, start - widest):end + widest])
        if self._index is not None:
            ws = max(0, start - widest)
            todo = self._index.candidates(read[ws:end + widest],
                                          start - ws, end - ws, widest,
                                          plain)
            if todo is not None:
                return self._covers_by(todo, read, start, end, plain)
        return self._covers_by(range(len(self.contexts)), read, start, end,
                               plain)

    def _covers_by(self, which, read: str, start: int, end: int,
                   plain: bool) -> bool:
        """covers() over the contexts numbered in `which`. With every
        context, this is the original loop; the index only ever narrows
        `which` to the contexts that could possibly match."""
        for i in which:
            ctx, search = self.contexts[i], self.searches[i]
            lo = max(0, start - len(ctx))
            hi = min(len(read), end + len(ctx))
            h = search.best(read[lo:hi], plain)
            if h is not None and lo + h[0] < end and lo + h[1] > start:
                return True
        return False


# Above this many contexts (both strands) covers() consults an index
# rather than trying every context. A panel has a handful and the plain
# loop is quicker than building the index; a database release has tens of
# thousands, where the loop took 0.3 s per anchor hit - weeks for a
# flow cell.
LOOKALIKE_INDEX_MIN = 64


class _ContextIndex:
    """Which look-alike contexts could possibly match a stretch of read.

    An exact filter, not a heuristic: every context it leaves out provably
    has no alignment within its tolerance anywhere in the stretch, so
    covers() gives the same answer with it as without it. It only saves
    the work of finding that out one context at a time.

    The argument is the pigeonhole principle, taken one step further than
    SeededAnchorScan takes it. Cut a context into k + 2 disjoint pieces,
    where k is its tolerance. An alignment with at most k edits can touch
    at most k of them, so at least TWO pieces occur exactly in the read -
    and in the right places relative to each other: between two exact
    pieces the alignment can have gained or lost at most k bases, so their
    offsets in the read differ from their offsets in the context by at
    most k. Each piece occurrence gives a diagonal (read offset minus
    context offset); a context is a candidate only when two of its piece
    occurrences fall on diagonals within k of each other.

    Why two pieces and not one. One exact piece - what SeededAnchorScan
    uses - is the standard filter, and over the 31,000 contexts of a
    database release it leaves ~1,200 candidates per anchor hit, because
    the pieces are only 7-8bp and the stretch around a hit is ~140bp: a
    random 7-mer is present about one time in 120. Requiring two on a
    consistent diagonal leaves a median of ~35, none of which match on
    real adapters. Three pieces (k + 3 shorter ones) was tried and is
    slower: more, shorter pieces cost more to look up than they save.

    How it is fast in pure Python. Every piece occurrence becomes one
    integer, context * STRIDE + OFFSET + diagonal, produced by map() and
    list.extend() rather than a Python loop. Sorting those integers puts
    each context's occurrences side by side in diagonal order, and two
    neighbours differing by at most kmax are a candidate pair - found with
    map() and compress(), again without a Python loop. STRIDE keeps
    different contexts' integers far enough apart that they can never
    pair. A release's integers stay below 2**30, so the sort compares
    them on the fast single-digit path. ~1 ms per anchor hit on an M1, against
    ~300 ms for the plain loop over the same contexts.

    Contexts the argument does not cover - too short to
    cut, or too ambiguous to spell out (see _cut) - are always candidates.
    """

    STRIDE = 512

    def __init__(self, contexts: list[str]):
        post: dict[str, list[int]] = {}
        self.always: list[int] = []
        self.kmax = 0
        self.offset = 0
        cut = []
        for i, ctx in enumerate(contexts):
            k = max(4, len(ctx) // 8)
            # Any length: the index only narrows the field, and a context
            # over 64bp is still verified by find_anchor_cost, as before.
            pieces = self._cut(ctx, k + 2)
            if not pieces:
                self.always.append(i)
                continue
            self.kmax = max(self.kmax, k)
            cut.append((i, pieces))
            self.offset = max(self.offset, max(a for a, _p in pieces))
        for i, pieces in cut:
            base = i * self.STRIDE + self.offset
            for a, spellings in pieces:
                for p in spellings:
                    post.setdefault(p, []).append(base - a)
        self.post = {p: tuple(v) for p, v in post.items()}
        self.lengths = sorted({len(p) for p in self.post})
        # The widest stretch whose integers stay clear of the next
        # context's: diagonals run from -offset to the stretch length.
        self.max_text = self.STRIDE - self.offset - self.kmax - 2

    @staticmethod
    def _cut(ctx: str, n: int):
        """`ctx` in n contiguous pieces, as (offset, every plain-ACGT
        spelling of the piece), or None.

        A piece holding an ambiguity code is still an exact match where
        the read has any base the code allows - that is how the aligner
        scores it - so it is indexed under each such spelling. PR2 keeps
        its ambiguity codes, and ~1,200 of a release's contexts carry
        one; left out of the index they were tried on every call and
        made it a hundred times slower. A piece with more than 1,024
        spellings (five Ns) is not worth indexing, and its context is
        always tried instead.
        """
        m = len(ctx)
        if m < n:
            return None
        cuts = [m * i // n for i in range(n + 1)]
        out = []
        for i in range(n):
            piece = ctx[cuts[i]:cuts[i + 1]]
            spellings = [""]
            for c in piece:
                bases = IUPAC.get(c, "")
                if not bases:
                    return None
                spellings = [s + b for s in spellings for b in bases]
                if len(spellings) > 1024:
                    return None
            out.append((cuts[i], spellings))
        return out

    def candidates(self, text: str, start: int, end: int,
                   widest: int, plain: bool = True) -> list[int] | None:
        """Context numbers that could align within tolerance in `text`
        OVERLAPPING text[start:end], in index order; None when `text` is
        too long for the integer packing, so the caller tries every
        context instead.

        Why the overlap matters. covers() only answers yes for an
        alignment that overlaps the anchor hit, and the stretch around a
        hit is mostly something else - on a VESPA read, the 18S insert
        begins ~37bp past the i5 anchor, and hundreds of near-identical
        release contexts genuinely align there. Every one of them was a
        candidate, was aligned, and was then rejected for not overlapping:
        14 ms per hit. An alignment that does overlap has a bounded
        position, and so do its exact pieces. If it starts at x, with at
        most k edits its pieces' diagonals lie within k of x, and it ends
        within k of x + m; overlapping [start, end) needs x < end and
        x + m + k > start. So every piece diagonal d satisfies
            start - m - 2k  <  d  <  end + k,
        and a pair outside that range cannot belong to an overlapping
        alignment. The widest context and largest tolerance make the
        bound hold for all of them at once. ~2 ms per hit again.

        A read with an ambiguity code near the hit (`plain` False) is
        still indexed: the aligner lets a read's N match any base, so a
        piece of read holding one is looked up under every spelling it
        allows. Without this, one N sent the call to every context in
        turn - three seconds for a single hit.
        """
        if len(text) > self.max_text:
            return None
        post, codes = self.post, []
        extend = codes.extend
        for L in self.lengths:
            for i in range(len(text) - L + 1):
                piece = text[i:i + L]
                if plain:
                    hit = post.get(piece)
                    if hit:
                        extend(map(i.__add__, hit))
                    continue
                spellings = [""]
                for c in piece:
                    spellings = [sp + b for sp in spellings
                                 for b in IUPAC.get(c, "")]
                    if len(spellings) > 256:
                        return None
                seen = set()
                for sp in spellings:
                    for code in post.get(sp, ()):
                        if code not in seen:
                            seen.add(code)
                            codes.append(code + i)
        out = set(self.always)
        if len(codes) > 1:
            codes.sort()
            near = map((self.kmax + 1).__gt__, map(sub, codes[1:], codes))
            stride, k = self.STRIDE, self.kmax
            # Diagonal bounds, shifted into the packed integers.
            lo = start - widest - 2 * k + self.offset
            hi = end + k + self.offset
            out.update(c // stride for c in compress(codes, near)
                       if lo < c % stride < hi)
        return sorted(out)


# ASCII only, deliberately. str.isalnum() also accepts letters such as
# 'é', and demultiplexing once named its files that way while pre-flight
# recorded '_' in their place - so the sanitised id in file_resolution.tsv
# and the filename that reached assignment disagreed, and the restore
# below matched nothing.
ID_SAFE = re.compile(r"[^A-Za-z0-9._-]")


def sanitise_id(value: str) -> str:
    """Path-safe form of an identifier. The original is always what gets
    reported; this is only what appears in directory and file names.

    The ONE implementation. Pre-flight records it in file_resolution.tsv
    and demultiplexing names its per-sample files with it, and the two
    must agree character for character for restore_declared_ids to find
    the declared id again."""
    return ID_SAFE.sub("_", value.strip())


def restore_declared_ids(keyed: dict, preflight_dir: str | None,
                         run_id: str = "") -> list[tuple[str, str]]:
    """Re-key a per-sample dict from filename-safe ids back to declared ids.

    Sample ids reach the later stages as FILENAMES. Demultiplexing writes
    one FASTQ per sample, the workflow takes the id from that file's
    name, and assignment stamps it into counts.tsv - so anything
    sanitise_id() had to replace arrives changed. Every downstream join
    is then against the samplesheet's ORIGINAL id and quietly matches
    nothing.

    Found on a plate whose sample ids contained '=':

        ampure_F1=400_F2=100_rep1    declared
        ampure_F1_400_F2_100_rep1    what arrived

    The earlier plates used hyphens, which survive sanitisation, so this
    sat unexposed until a naming convention changed. The symptom was a
    cohort check reporting "88 declared, 88 reached" and then listing all
    88 as lost - incoherent on its face, which is what made it visible.
    What it would have done unnoticed is worse: sample_type, indices,
    plate position and metadata all join on the declared id, and all of
    them silently no-op on a miss, so negative controls stop being
    negative controls and detection limits are computed from nothing.

    file_resolution.tsv carries both spellings side by side for exactly
    this purpose, so the mapping is read rather than recomputed - a
    second implementation of sanitise_id here is one that can disagree
    with the first.

    Mutates `keyed` in place and returns the renames made, so the caller
    can report them. Ambiguity is declined rather than guessed: a
    sanitised form claimed by two declared ids is left alone (pre-flight
    already refuses those), and so is one whose declared id is in use.
    """
    if not keyed or not preflight_dir:
        return []
    path = os.path.join(preflight_dir, "file_resolution.tsv")
    if not os.path.exists(path):
        return []

    owners: dict[str, list[str]] = {}
    for r in read_table(path):
        if run_id and r.get("run_id", "") != run_id:
            continue
        safe, declared = r.get("sanitised_id", ""), r.get("sample_id", "")
        if safe and declared and safe != declared:
            owners.setdefault(safe, []).append(declared)

    renamed: list[tuple[str, str]] = []
    for safe, declared_ids in sorted(owners.items()):
        if len(set(declared_ids)) != 1:
            continue                      # ambiguous: two ids, one filename
        declared = declared_ids[0]
        if safe not in keyed or declared in keyed:
            continue
        value = keyed.pop(safe)
        if hasattr(value, "sample_id"):
            value.sample_id = declared
        keyed[declared] = value
        renamed.append((safe, declared))
    return renamed


def delimiter_for(path: str) -> str:
    """Comma for .csv, tab for everything else. The ONE copy of this rule.

    It lived in three files - here, preflight.py and simulate.py - as
    three identical lines that had to agree. They did, which is the only
    reason it was not already a bug: a reader and a writer that disagree
    about a delimiter do not raise, they produce one wide column whose
    header matches nothing, and the complaint that comes back names
    every required field as missing. That sends the reader to check
    their headers, where nothing is wrong.

    The rule is the file EXTENSION rather than sniffing the content,
    which makes it per-file and therefore mixable: a TSV samplesheet
    beside a CSV one both work, and neither has to be converted. That
    is what lets the samplesheet move to TSV without invalidating
    every sheet already written.
    """
    return "," if path.lower().endswith(".csv") else "\t"


def read_table(path: str) -> list[dict[str, str]]:
    delim = delimiter_for(path)
    # utf-8-sig, not utf-8: Excel's "CSV UTF-8" writes a byte-order mark,
    # and read as plain utf-8 it becomes part of the first header - so
    # '\ufeffrun_id' is a column nobody asks for and every lookup on
    # run_id misses. utf-8-sig drops a BOM if there is one and reads plain
    # UTF-8 identically.
    with open(path, newline="", encoding="utf-8-sig") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    return [
        {k.strip(): (v or "").strip() for k, v in row.items()}
        for row in csv.DictReader(lines, delimiter=delim)
    ]


def load_assay(path: str | None) -> dict:
    cfg = dict(DEFAULTS)
    if path and os.path.exists(path):
        # utf-8-sig for the same reason as read_table: an editor that saves
        # a BOM would otherwise make json.load refuse the file outright.
        with open(path, encoding="utf-8-sig") as fh:
            data = json.load(fh)
        lib = {
            k: v for k, v in data.get("library", {}).items() if not k.startswith("_")
        }
        cfg.update({k: v for k, v in lib.items() if k in DEFAULTS})
    return cfg


FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def resolve_read_files(paths, prefer_pass: bool = True) -> tuple[list, dict]:
    """Work out which FASTQ files a nanopore run actually consists of.

    A MinION does not produce a FASTQ. It produces a directory tree of
    them, written in batches while the run proceeds, split into
    fastq_pass and fastq_fail, and split again into per-barcode folders
    if MinKNOW was asked to demultiplex. Every one of those splits is a
    decision about which reads enter the analysis, so each is made
    explicitly here and reported rather than absorbed by a glob.

    Two of them matter enough to be worth stating.

    fastq_fail is excluded. Those reads failed the basecaller's own
    quality filter, and including them silently would import a different
    error distribution into the middle of a demultiplexer whose barcode
    tolerances were set for passing reads. They are counted and reported
    so the exclusion is visible.

    MinKNOW's own barcode folders are a trap for this assay. The barcodes
    here are in-line and tandem-duplicated inside the amplicon; MinKNOW
    knows nothing about them and will have been matching its own kit
    barcodes, which are not present. Its `unclassified` folder is
    therefore where the real data is, not where the failures are - so
    everything is taken and the situation is flagged.

    Returns (files, report). Files are sorted, so a run processes in the
    same order every time.
    """
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]

    found: list[str] = []
    report = {
        "inputs": [str(p) for p in paths],
        "n_files": 0,
        "from_fastq_pass": 0,
        "excluded_fastq_fail": 0,
        "ignored_outside_fastq_pass": [],
        "minknow_barcode_dirs": [],
        "notes": [],
    }

    def is_fastq(name: str) -> bool:
        return name.lower().endswith(FASTQ_SUFFIXES)

    for p in paths:
        p = str(p)
        if os.path.isfile(p):
            found.append(p)
            continue
        if not os.path.isdir(p):
            hits = sorted(glob.glob(p))
            if not hits:
                report["notes"].append(f"nothing matched '{p}'")
            found.extend(h for h in hits if is_fastq(h))
            continue

        passed, failed, other = [], [], []
        for root, dirs, files in os.walk(p):
            dirs.sort()
            rel = os.path.relpath(root, p).lower()
            parts = set(rel.split(os.sep))
            # Note the barcode folders BEFORE deciding whether this level
            # holds files. In a MinKNOW tree fastq_pass contains only
            # directories, so a check that runs after "no files here" is a
            # check that never runs.
            for d in dirs:
                if d.lower().startswith("barcode") or d.lower() == "unclassified":
                    report["minknow_barcode_dirs"].append(
                        os.path.relpath(os.path.join(root, d), p))

            here = [os.path.join(root, f) for f in sorted(files)
                    if is_fastq(f)]
            if not here:
                continue
            if "fastq_fail" in parts:
                failed.extend(here)
            elif "fastq_pass" in parts:
                passed.extend(here)
            else:
                other.extend(here)

        if passed:
            found.extend(passed)
            report["from_fastq_pass"] += len(passed)
            # Once a run has a fastq_pass, fastq_pass IS the run. Anything
            # else in the tree is something made from it afterwards - most
            # often an all_pass.fastq.gz someone built with cat to have one
            # file to hand round. Taking that too reads every read twice,
            # and deduplicating by path cannot catch it because it is a
            # different file holding the same reads. So they are left out
            # and named, so a file that genuinely belongs can be pointed
            # at directly.
            if other:
                report["ignored_outside_fastq_pass"].extend(
                    os.path.relpath(f, p) for f in other)
                shown = ", ".join(os.path.relpath(f, p) for f in other[:3])
                more = f" and {len(other) - 3} more" if len(other) > 3 else ""
                report["notes"].append(
                    f"{len(other)} FASTQ outside fastq_pass were ignored "
                    f"({shown}{more}) - with a fastq_pass present those are "
                    "usually copies of it, and reading them would count "
                    "every read twice. Pass a file directly to use it"
                )
        else:
            found.extend(other)
        if failed:
            report["excluded_fastq_fail"] += len(failed)
            if prefer_pass:
                report["notes"].append(
                    f"{len(failed)} file(s) under fastq_fail were excluded - "
                    "those reads failed the basecaller's quality filter and "
                    "carry a different error profile from the one the "
                    "barcode tolerances were set for"
                )
            else:
                found.extend(failed)

    # A run tree can present the same file by more than one route.
    seen, uniq = set(), []
    for f in found:
        key = os.path.realpath(f)
        if key not in seen:
            seen.add(key)
            uniq.append(f)
    uniq.sort()
    report["n_files"] = len(uniq)
    if report["minknow_barcode_dirs"]:
        report["notes"].append(
            "MinKNOW barcode folders are present, so its demultiplexing was "
            "on. It does not know the in-line tandem barcodes this assay "
            "uses, so its assignment is not meaningful here and its "
            "'unclassified' reads are ordinary data. All folders were read; "
            "check that MinKNOW barcoding was not also trimming"
        )
    return uniq, report


def open_reads_many(paths):
    """Stream reads from any number of FASTQ files as one sequence."""
    for p in paths:
        yield from open_reads(p)


def open_reads(path: str):
    """Yield (read_id, sequence, quality) from a FASTQ, gzipped or not."""
    # Lower-cased, as is_fastq does when deciding the file is a FASTQ at
    # all. Checked case-sensitively here, a READS.FASTQ.GZ was accepted
    # by the one and then opened as text by the other.
    op = gzip.open if str(path).lower().endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            # A blank line where a header should be - a trailing empty
            # line, most often, left by an editor or a careless cat - is
            # not a record. Taken as one it raised IndexError on the id.
            if not head.strip():
                continue
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            ident = head[1:].split()
            yield (ident[0] if ident else ""), seq, qual


def read_run_metadata(paths) -> dict:
    """What produced these reads, taken from the first read's header.

    Dorado writes SAM-style tags after the read id, and one of them says
    which basecalling model was used:

        RG:Z:<run uuid>_dna_r10.4.1_e8.2_400bps_hac@v5.2.0

    That belongs in the outputs. A basecalling model is not a neutral
    detail - it changes the error profile, and therefore the barcode
    recovery, the anchor hit rate and ultimately how many reads reach
    each sample. Two runs basecalled with different models can differ in
    composition for reasons that have nothing to do with biology, and
    nothing downstream would reveal it. Recording the model is what makes
    that check possible later.

    Best effort by design: the tags are absent from reads written by
    anything that rewrote the headers, and absent entirely from older
    files. An empty dict is a fine answer and never an error.
    """
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        try:
            op = gzip.open if str(path).lower().endswith(".gz") else open
            with op(path, "rt") as fh:
                head = fh.readline().rstrip("\n")
        except OSError:
            continue
        if not head.startswith("@"):
            continue
        fields = head[1:].split("\t")[1:]
        tags = {}
        for f in fields:
            parts = f.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
        if not tags:
            continue
        out = {}
        rg = tags.get("RG", "")
        if rg:
            # RG is "<36-char run uuid>_<model>". Split on the uuid width
            # rather than on the first underscore, because model names
            # contain underscores of their own.
            if len(rg) > 37 and rg[36] == "_":
                out["acquisition_id"], out["basecall_model"] = rg[:36], rg[37:]
            else:
                out["basecall_model"] = rg
        for tag, name in (("PU", "flow_cell"), ("LB", "library"),
                          ("DT", "run_started")):
            if tags.get(tag):
                out[name] = tags[tag]
        if out:
            # Deliberately NOT recording which file this came from. The
            # same reads sharded differently would then report different
            # provenance, and a run must not describe itself differently
            # according to how its input happened to be batched.
            return out
    return {}


def phred(qual: str) -> list[int]:
    return [ord(c) - 33 for c in qual]


def expected_errors(qual: str) -> float:
    """Sum of per-base error probabilities.

    A better filter than mean quality, which is an average of logarithms and
    so lets a read through on the strength of its good bases while a handful
    of terrible ones decide the outcome.
    """
    total = 0.0
    for c in qual:
        total += 10.0 ** (-(ord(c) - 33) / 10.0)
    return total


# --------------------------------------------------------------------------
# reference database releases (--reference_database)
# --------------------------------------------------------------------------
#
# A release is one directory built by the vespa-18s-refdb repository, and
# its layout is the contract between the two (docs/REFERENCE_DATABASE.md
# section 4.1). Pre-flight validates it and assignment reads it; both read
# it through here, so the contract is written down once in code.

# The nine PR2 ranks, shallowest first. taxonomy.tsv carries one column
# for each, blank where PR2 has only a placeholder.
RANKS = ["domain", "supergroup", "division", "subdivision", "class",
         "order", "family", "genus", "species"]
# NCBI's eight, as nemabiome-its2-refdb releases carry them: real ranks,
# no placeholders and no rank_names.tsv.
NCBI_RANKS = ["superkingdom", "kingdom", "phylum", "class", "order",
              "family", "genus", "species"]
# The rank sets a taxonomy.tsv header is recognised by, PR2's first so a
# vespa release reads exactly as it always has.
KNOWN_RANKS = (RANKS, NCBI_RANKS)

# What the pipeline reads. A release missing any of these cannot be used.
REFDB_REQUIRED = ["manifest.json", "refdb.fasta", "taxonomy.tsv",
                  "anchor_lookalikes.tsv", "resolution_nanopore.tsv"]
# The rest of the contract: not read by the nanopore route, but their
# absence means the release was not built the agreed way.
REFDB_EXPECTED = ["members.tsv", "conflicts.tsv", "primer_fit.tsv"]

LOOKALIKE_COLUMNS = ["seq_id", "anchor", "strand", "position", "edits",
                     "context"]


def load_rank_names(release: str) -> list[dict]:
    """The release's rank_names.tsv (v0.8.0 on): what each positional
    taxonomy column holds, by clade. [] when the release has none.

    The nine columns of taxonomy.tsv are PR2's positional levels, not
    Linnaean ranks. In nematodes `family` holds an order (Rhabditida) or
    a placeholder; in Ascomycota `order` holds the subphylum and `family`
    the class. Read as named, "Chromadorea, unnamed family" is really an
    unnamed order, and a call "at family" in the fungi is at class. Rows
    are matched in order on one lineage column (`match_column` equal to
    `match_value`); the values may list several ranks, commonest first
    ("class|order"), or "placeholder" where the column carries nothing.
    """
    path = os.path.join(release, "rank_names.tsv")
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def rank_meaning(rank_names: list[dict], lineage: list[str],
                 column: str) -> str:
    """What `column` holds for a lineage, from rank_names.tsv: "order",
    "class|order", "placeholder"... or "" where no row applies (the
    column then means roughly what it says)."""
    for row in rank_names:
        k = RANKS.index(row["match_column"]) if row.get(
            "match_column") in RANKS else -1
        if 0 <= k < len(lineage) and lineage[k] == row.get("match_value"):
            return row.get(column, "") or ""
    return ""


def rank_word(rank_names: list[dict], lineage: list[str], column: str) -> str:
    """The rank to SAY for a positional column: its commonest real rank
    from rank_names.tsv, or the column's own name. "family" for a
    nematode is "order"; for Apicomplexa, or with no table, "family"."""
    m = rank_meaning(rank_names, lineage, column)
    real = [x for x in m.split("|") if x and x != "placeholder"]
    return real[0] if real else column


def load_release_manifest(release: str) -> dict:
    """manifest.json of a release, or ValueError saying what is wrong."""
    path = os.path.join(release, "manifest.json")
    if not os.path.isfile(path):
        raise ValueError(f"no manifest.json in {release}")
    try:
        with open(path, encoding="utf-8-sig") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path} is not valid JSON - {exc}") from None
    if not isinstance(data, dict):
        raise ValueError(f"{path} is not a JSON object")
    return data


def release_version(manifest: dict) -> str:
    """The name every result carries: the build version, or the release
    name when an older manifest has no version field."""
    return str(manifest.get("version") or manifest.get("name") or "unknown")


def illumina_windows(manifest: dict, read_length: int | None = None,
                     window_r1: int | None = None,
                     window_r2: int | None = None) -> tuple:
    """(window_r1, window_r2, resolution table or None) for Illumina pairs
    of `read_length`.

    A release lists the windows it resolved at each read length
    (illumina_geometries: 132/128 at 2x150, 250/200 at 2x300). The
    windows are those of the longest geometry the reads can hold, so a
    2x301 run reads past the 2x150 windows - Onchocerca and Loa differ by
    one base at insert 347, between 132 and 128 but inside 250/200. A
    read length under every geometry, or none given, takes the manifest's
    illumina_windows, as do releases that list no geometries. Explicit
    windows override either mate. The table is the geometry whose windows
    are the ones used; None when no geometry has them, since another
    geometry's table would describe reads that were not taken."""
    geoms = [g for g in manifest.get("illumina_geometries") or []
             if isinstance(g, dict) and len(g.get("windows") or []) == 2]
    base = list(manifest.get("illumina_windows") or [132, 128])
    fits = [g for g in geoms
            if read_length and int(g.get("read_length", 0)) <= read_length]
    if fits:
        base = max(fits, key=lambda g: int(g["read_length"]))["windows"]
    w1 = int(window_r1 or base[0])
    w2 = int(window_r2 or base[1])
    table = next((g.get("file") for g in geoms
                  if [int(x) for x in g["windows"]] == [w1, w2]), None)
    return w1, w2, table


def read_fasta(path: str):
    """(id, sequence) pairs; the id is the header up to the first space."""
    name, parts = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def open_text(path: str):
    """A text handle on a plain or gzipped file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def ranks_in(header, manifest: dict | None = None) -> list[str] | None:
    """The rank columns of a taxonomy.tsv header, shallowest first: the
    manifest's `ranks` when it gives them, else the first known set
    (PR2's nine, NCBI's eight) whose every column is present. None when
    neither fits."""
    given = (manifest or {}).get("ranks")
    if isinstance(given, list) and given and all(r in header for r in given):
        return [str(r) for r in given]
    return next((list(rs) for rs in KNOWN_RANKS
                 if all(r in header for r in rs)), None)


# Per-record columns kept beside the ranks when a release has them
# (nemabiome-its2-refdb): species_status says whether a record with no
# species names one (named), two or more (shared) or none (unnamed);
# grade and coverage how much of the marker it holds.
RECORD_FIELDS = ("species_status", "grade", "coverage")


class Taxonomy(dict):
    """id -> {rank: name}, which is what every caller indexes, plus the
    release's rank columns as `ranks` and, for shared records, the
    species their members name (`shared`)."""

    def __init__(self, ranks: list[str]):
        super().__init__()
        self.ranks = ranks
        self.shared: dict[str, list[str]] = {}


def load_taxonomy(path: str, manifest: dict | None = None) -> Taxonomy:
    """id -> {rank: name} over the release's ranks, blanks kept as '';
    RECORD_FIELDS too where the release has them. The ranks are read from
    the header (ranks_in), so PR2 and NCBI releases both load; a
    members.tsv beside the table gives what each shared record's members
    name."""
    with open_text(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        ranks = ranks_in(rdr.fieldnames or [], manifest) or RANKS
        extra = [k for k in RECORD_FIELDS if k in (rdr.fieldnames or [])]
        out = Taxonomy(ranks)
        for r in rdr:
            rec = {k: (r.get(k) or "").strip() for k in ranks}
            for k in extra:
                rec[k] = (r.get(k) or "").strip()
            out[r["id"]] = rec
    members = os.path.join(os.path.dirname(path), "members.tsv")
    shared = {i for i, r in out.items() if r.get("species_status") == "shared"}
    if shared and os.path.exists(members):
        with open_text(members) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("id") in shared and r.get("role") == "kept":
                    sp = (r.get("label") or "").split(";")[-1].strip()
                    names = out.shared.setdefault(r["id"], [])
                    if sp and sp not in names:
                        names.append(sp)
        for names in out.shared.values():
            names.sort()
    return out


# --------------------------------------------------------------------------
# plate position: row, column and edge effects
# --------------------------------------------------------------------------
#
# Shared by the panel synthesis and the database synthesis. Wells hold
# unrelated samples, so where a well sits should explain nothing about
# it; when it does, the cause is the plate - a multichannel channel, a
# thermal gradient, evaporation at the edges, a bad column of tips.


def kruskal_h(groups: list[list[float]]) -> float:
    """Rank-based test statistic for 'do these groups differ'.

    Ranks rather than values because read counts are heavily skewed and a
    handful of heavy infections would otherwise drive any comparison of
    means. No distributional assumption is made, and none is needed - the
    p-value comes from permutation below.
    """
    vals = [v for g in groups for v in g]
    n = len(vals)
    if n < 3 or len(groups) < 2:
        return 0.0
    order = sorted(range(n), key=lambda i: vals[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    idx = 0
    total = 0.0
    for g in groups:
        rs = sum(ranks[idx: idx + len(g)])
        if g:
            total += rs * rs / len(g)
        idx += len(g)
    return 12.0 / (n * (n + 1)) * total - 3 * (n + 1)


def _midranks(vals: list[float]) -> list[float]:
    """Each value's rank, ties averaged - kruskal_h's ranks, in input order."""
    n = len(vals)
    order = sorted(range(n), key=lambda i: vals[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def permutation_p(
    values: list[float], labels: list, rng: random.Random, n_perm: int
) -> tuple[float, float]:
    """How unusual is the grouping, if well position means nothing?

    Wells are meant to be independent, so the null is simply that the
    labels could have been attached to any well. Shuffling them is a more
    honest null than any parametric one, and it costs nothing here.

    Shuffling labels never changes a value's rank, so the ranks are taken
    once and each permutation only re-sums them by label - the statistic
    kruskal_h would give, to the bit (half-integer sums are exact), at a
    fraction of the cost, which is what lets a 384-well plate be tested
    taxon by taxon.
    """
    by = defaultdict(list)
    for v, lab in zip(values, labels):
        by[lab].append(v)
    groups = [g for g in by.values() if g]
    if len(groups) < 2:
        return 0.0, 1.0
    observed = kruskal_h(groups)

    n = len(values)
    ranks = _midranks(values)
    scale = 12.0 / (n * (n + 1)) if n >= 3 else 0.0
    shuffled = list(labels)
    worse = 0
    for _ in range(n_perm):
        rng.shuffle(shuffled)
        sums: dict = {}
        sizes: dict = {}
        for r, lab in zip(ranks, shuffled):
            if lab in sums:
                sums[lab] += r
                sizes[lab] += 1
            else:
                sums[lab] = r
                sizes[lab] = 1
        if n < 3:
            h = 0.0
        else:
            total = 0.0
            for lab, rs in sums.items():
                total += rs * rs / sizes[lab]
            h = scale * total - 3 * (n + 1)
        if h >= observed:
            worse += 1
    # +1 on both sides: a permutation test can never justify p = 0, and
    # reporting one would overstate what the resampling can support.
    return observed, (worse + 1) / (n_perm + 1)


def plate_effects(
    rows: list, cols: list, series: list[tuple[str, list[float]]],
    rng: random.Random, n_perm: int, min_wells: int = 12,
) -> list[dict]:
    """Row, column and edge p-values for each (name, per-well values).

    `rows` and `cols` give each well's position, in the order of the
    values. A quantity with fewer than three distinct values is skipped:
    there is nothing for position to explain. Fewer than `min_wells`
    placed wells and nothing is tested at all - a permutation test on a
    handful of wells cannot tell a plate effect from chance.
    """
    if len(rows) < min_wells:
        return []
    rs, cs = sorted(set(rows)), sorted(set(cols))
    edge = ["edge" if (r in (rs[0], rs[-1]) or c in (cs[0], cs[-1]))
            else "interior" for r, c in zip(rows, cols)]
    findings = []
    for name, vals in series:
        if len(set(vals)) < 3:
            continue
        _h, row_p = permutation_p(vals, list(rows), rng, n_perm)
        _h, col_p = permutation_p(vals, list(cols), rng, n_perm)
        _h, edge_p = permutation_p(vals, edge, rng, n_perm)
        findings.append({
            "quantity": name,
            "n_wells": len(rows),
            "row_p": round(row_p, 4),
            "col_p": round(col_p, 4),
            "edge_p": round(edge_p, 4),
        })
    return findings


# --------------------------------------------------------------------------
# index misassignment on combinatorial plates
# --------------------------------------------------------------------------
#
# On a plate indexed by combination - each i5 a row, each i7 a column -
# a read whose i7 is misread lands in another well of the SAME ROW (it
# keeps its i5), and a read whose i5 is misread lands in the SAME
# COLUMN. A well sharing neither index needs both misread, and gets
# next to nothing. So the reads of taxon T that well w receives by
# misassignment are about
#
#     rate_row x (T reads in the other wells of w's row)
#   + rate_col x (T reads in the other wells of w's column)
#
# where each rate is per PAIR of wells: the share of one well's reads
# that land in one particular other well along that axis. Per pair,
# because a misread lands on whichever declared index is nearest, so the
# rate belongs to the index chemistry rather than to how many wells the
# plate has - and it is comparable from run to run.
#
# The rates are measured on each run from the run's own data, pooled
# over every taxon: flow cell, basecaller and library all move them. For
# a (well, taxon) pair, true presence can only ADD to what misassignment
# put there, so every pair's ratio reads / (reads next door) is at least
# the rate, and the pairs where the taxon is truly absent pile up AT it.
# The rate is that pile: the lowest well-supported mode of the ratios.
# Nothing here needs an unused index combination, so it works on a full
# plate; where there are unused combinations their reads give a direct
# measurement, reported beside it as a cross-check.

# Where the run cannot pin a rate, calls are still checked against this
# one, per pair. Twice what the first combinatorial nanopore plates
# showed (~5e-4), on the side of calling misassignment.
MISASSIGN_DEFAULT_RATE = 1e-3
# A pair informs the fit only when its sources hold this many reads of
# the taxon: under ~10 expected reads a leak is lost to Poisson noise
# and to clustering, which never forms a consensus from fewer than 5.
MISASSIGN_MIN_SOURCE = 10000
# Once a rate is found, a pair is kept only where it predicts this many
# leaked reads: at a quarter of that, the 5-read minimum loses nothing.
MISASSIGN_MIN_EXPECTED = 10.0
# ...and when one axis carries at least this multiple of the other's
# source reads, so the pair measures that axis.
MISASSIGN_DOMINANCE = 3.0
# Fewer informative pairs than this, or a mode carried by fewer pairs,
# and the axis is not fitted.
MISASSIGN_MIN_PAIRS = 20
MISASSIGN_MIN_SUPPORT = 8
# The placebo must be beaten this convincingly: with 100 shuffles, by
# all but one. Steady contamination from the whole plate - not
# misassignment - passed one test in sixteen at this level in
# simulation; every real leak simulated passed.
MISASSIGN_PLACEBO_P = 0.02
MISASSIGN_PLACEBOS = 100
# A fitted rate whose 95% interval spans more than this factor is not a
# measurement.
MISASSIGN_MAX_SPREAD = 10.0
# Nor is a rate above this, per pair: 0.5% of a well's reads into EACH
# other well of its row would be a tenth of the well gone on a 384-well
# plate, three times the worst seen (1.5e-3 on the first combinatorial
# nanopore plates; Illumina hopping at its worst, ~2% of reads spread
# over a row, is ~1e-3). A cluster that high is shared biology - wells
# laid out by group, a column of one host, a reagent taxon next to a
# blank full of it - lined up with an axis, not misreads.
MISASSIGN_MAX_RATE = 0.005
# Half-width of the window, in log10 units, used to find the mode.
_MODE_HALF = 0.25
_MODE_STEP = 0.05
# A mode must hold at least this share of the densest window's pairs.
_MODE_REL = 0.25
AXES = ("row", "col")


class IndexSharing:
    """Which wells share an index with which, and how many reads of each
    taxon sit next door along each axis.

    `geom` maps a well to its (row, column) - from the index map, the i5
    giving the row and the i7 the column - and `counts` maps a well to
    {taxon key: reads}. Wells without a position are simply absent.
    """

    def __init__(self, geom: dict, counts: dict):
        self.geom = geom
        self.counts = counts
        self.row_tot: dict = defaultdict(int)
        self.col_tot: dict = defaultdict(int)
        self.tot: dict = defaultdict(int)
        for sid, (r, c) in geom.items():
            for k, v in (counts.get(sid) or {}).items():
                self.row_tot[(r, k)] += v
                self.col_tot[(c, k)] += v
                self.tot[k] += v

    def sources(self, sid: str, key) -> tuple[int, int, int]:
        """Reads of `key` in the OTHER wells of sid's row, of its column,
        and in wells sharing neither index."""
        r, c = self.geom[sid]
        y = (self.counts.get(sid) or {}).get(key, 0)
        xr = self.row_tot[(r, key)] - y
        xc = self.col_tot[(c, key)] - y
        return xr, xc, self.tot[key] - xr - xc - y

    def expected(self, sid: str, key, rate_row: float,
                 rate_col: float) -> float:
        if sid not in self.geom:
            return 0.0
        xr, xc, _xn = self.sources(sid, key)
        return rate_row * xr + rate_col * xc

    def pairs(self, min_source: int, keys=None) -> list[tuple]:
        """(well, reads, row sources, column sources, taxon) for every
        placed well and every taxon with min_source reads on the run -
        zeros included, which are as much evidence as reads."""
        keys = [k for k, t in self.tot.items() if t >= min_source
                and (keys is None or k in keys)]
        out = []
        for sid in self.geom:
            cs = self.counts.get(sid) or {}
            for k in keys:
                xr, xc, _xn = self.sources(sid, k)
                out.append((sid, cs.get(k, 0), xr, xc, k))
        return out


def _quantile(vals: list[float], q: float) -> float:
    v = sorted(vals)
    x = q * (len(v) - 1)
    i = int(x)
    return v[i] if i + 1 >= len(v) else v[i] + (v[i + 1] - v[i]) * (x - i)


def lowest_mode(logs: list[float], min_support: int = MISASSIGN_MIN_SUPPORT
                ) -> tuple[float, int] | None:
    """The lowest well-supported cluster in a set of log10 ratios.

    Density is counted in a sliding window +-_MODE_HALF wide; the first
    window from the bottom that holds at least min_support pairs and a
    quarter of the densest window's is climbed to its local peak, and the
    median of the pairs within twice the half-width of that peak is the
    answer. Returns (log10 rate, pairs in the window), or None.

    The lowest, not the densest: pairs where the taxon is really present
    can outnumber those where it is absent (a gut bacterium in every
    stool), and they sit above the rate, never below it.
    """
    import bisect
    if len(logs) < min_support:
        return None
    s = sorted(logs)
    lo, hi = s[0], s[-1]
    grid = [lo + i * _MODE_STEP
            for i in range(int((hi - lo) / _MODE_STEP) + 2)]
    dens = [bisect.bisect_right(s, g + _MODE_HALF)
            - bisect.bisect_left(s, g - _MODE_HALF) for g in grid]
    need = max(min_support, _MODE_REL * max(dens))
    for i, d in enumerate(dens):
        if d < need:
            continue
        j = i
        while j + 1 < len(dens) and dens[j + 1] >= dens[j]:
            j += 1
        c = grid[j]
        win = [x for x in s if abs(x - c) <= 2 * _MODE_HALF]
        return win[len(win) // 2], len(win)
    return None


def _axis_logs(pairs: list[tuple], axis: str, other_rate: float,
               min_source: int, dominance: float
               ) -> tuple[list[tuple[str, float | None]], list[int]]:
    """(well, log10 ratio or None for nothing seen) for the pairs that
    measure `axis`, with what the other axis is expected to have put
    there taken off first; and their source read counts."""
    out, xs = [], []
    for sid, y, xr, xc, *_k in pairs:
        x, xo = (xr, xc) if axis == "row" else (xc, xr)
        if x < min_source or x < dominance * xo:
            continue
        v = y - other_rate * xo
        out.append((sid, math.log10(v / x) if v > 0 else None))
        xs.append(x)
    return out, xs


def _mostly_empty(logs: list, xs: list, rate: float) -> bool:
    """Half or more of the pairs where `rate` predicts 20+ leaked reads
    show none at all: the rate is not what put reads in the others."""
    sure = [v for (_s, v), x in zip(logs, xs) if rate * x >= 20]
    return len(sure) >= MISASSIGN_MIN_SUPPORT and \
        sum(1 for v in sure if v is None) >= 0.5 * len(sure)


def _low_tail_share(pairs: list[tuple], axis: str, rate: float,
                    min_source: float, dominance: float) -> tuple[float, int]:
    """Among the pairs whose ratio is at most half a decade above `rate`
    (nothing seen included), the share within half a decade of it; and
    how many pairs that is."""
    logs, _xs = _axis_logs(pairs, axis, 0.0, min_source, dominance)
    lr = math.log10(rate)
    low = [v for _s, v in logs if v is None or v <= lr + 0.5]
    if not low:
        return 0.0, 0
    return sum(1 for v in low if v is not None and v >= lr - 0.5) / len(low), \
        len(low)


def _placebo_layout(geom: dict, rng: random.Random) -> dict:
    """The plate's wells dealt out again over its positions at random."""
    sids = sorted(geom)
    pos = [geom[s] for s in sids]
    rng.shuffle(pos)
    return dict(zip(sids, pos))


def fit_misassignment(
    geom: dict, counts: dict, rng: random.Random, *,
    min_source: int = MISASSIGN_MIN_SOURCE,
    dominance: float = MISASSIGN_DOMINANCE,
    min_pairs: int = MISASSIGN_MIN_PAIRS,
    min_support: int = MISASSIGN_MIN_SUPPORT,
    max_spread: float = MISASSIGN_MAX_SPREAD,
    max_rate: float = MISASSIGN_MAX_RATE,
    default_rate: float = MISASSIGN_DEFAULT_RATE,
    n_boot: int = 200, n_placebo: int = MISASSIGN_PLACEBOS,
    placebo_p: float = MISASSIGN_PLACEBO_P, taxa=None,
    groups: dict | None = None,
) -> dict:
    """The run's per-pair misassignment rate along rows and columns.

    For each axis the result says how the rate was arrived at:

      fitted       measured on this run: the lowest mode of the ratios,
                   with a 95% interval from resampling destination wells
      upper_bound  wells that should have shown a leak mostly show
                   nothing: the rate is below what this run can see
      default      the run cannot measure it (too few wells next to an
                   abundant source, no clear mode, an interval too wide,
                   a low cluster that does not follow index sharing, or
                   one too high to be misreading at all)

    and gives two rates. `for_calls` checks each call against what
    misassignment could have put there: the fitted rate, nothing under a
    bound (the run showed none where it would show), else the
    conservative default. `for_blanks` takes misassignment off
    the blanks before they set a limit - and is the fitted rate only,
    never an assumed one: subtracting an assumption from a blank would
    lower every limit on a guess.

    A pair is used once its sources would leak 20+ reads at the axis's
    rate (and never under min_source reads): a leak of a few reads is
    either lost - a consensus needs 5 - or kept only when it happened to
    come out high, and either way it misstates the rate. The fit starts
    from the deepest pairs for the same reason.

    With `groups` ({taxon: group}), each group with enough pairs is also
    fitted on its own. Presence can only lift a group's ratios, so where
    one group's whole 95% interval sits below the pooled rate, the pool
    is being lifted by taxa really present next door - diet plants and
    gut flora in most stool wells - and the lowest such group's rate is
    used instead (`from_group`). Where the groups agree, as they do on
    simulated plates, the pooled fit stands unchanged.

    The validity check is a placebo. The wells are dealt out again over
    the plate's positions at random: leakage that really follows index
    sharing makes a tighter low cluster with the true neighbours than
    with shuffled ones. A cluster no tighter than chance - contamination
    from the whole plate at a steady level, say - is not misassignment,
    and is not fitted. It cannot be told apart where the plate itself
    is balanced so that every row (or column) holds the same heavy
    sources: then a well's neighbours and the whole plate predict the
    same reads, and the data cannot say which put them there.
    """
    out: dict = {"applicable": False, "wells_placed": len(geom),
                 "settings": {"min_source_reads": min_source,
                              "dominance": dominance,
                              "min_pairs": min_pairs,
                              "min_support": min_support,
                              "max_spread": max_spread,
                              "max_rate": max_rate,
                              "default_rate": default_rate,
                              "bootstrap": n_boot, "placebos": n_placebo,
                              "placebo_p": placebo_p}}
    rows = {r for r, _c in geom.values()}
    cols = {c for _r, c in geom.values()}
    if len(rows) < 2 or len(cols) < 2:
        out["note"] = ("no combinatorial layout: misassignment along rows "
                       "and columns needs at least two of each"
                       if geom else
                       "no plate positions (the index map gave none), so "
                       "no well is known to share an index with another")
        return out
    out["applicable"] = True
    pairs = IndexSharing(geom, counts).pairs(min_source, taxa)

    def other(ax):
        return "col" if ax == "row" else "row"

    def threshold(r):
        return max(min_source, MISASSIGN_MIN_EXPECTED / r) if r > 0 \
            else min_source

    # Start from the deepest pairs, where clustering's 5-read minimum
    # cannot bias anything, and let the rate found there set how shallow
    # a pair may be. Started shallow instead, the pairs that survive the
    # minimum only by coming out high pull the rate up, a higher rate
    # admits shallower pairs, and the two feed each other.
    start = {}
    for ax in AXES:
        xs = sorted((x for _s, _y, xr, xc, *_k in pairs
                     for x, xo in [((xr, xc) if ax == "row" else (xc, xr))]
                     if x >= dominance * xo), reverse=True)
        k = 4 * min_pairs
        start[ax] = max(min_source, xs[k - 1] if len(xs) >= k else 0)
    thr = dict(start)

    # Both axes together: each is measured with what the other is
    # expected to have added taken off - but only an axis with enough
    # pairs to be measured at all is subtracted from the other.
    rate = {ax: 0.0 for ax in AXES}
    modes: dict = {}
    logs: dict = {}
    used: dict = {}
    for _ in range(20):
        new = {}
        for ax in AXES:
            o = other(ax)
            sub = rate[o] if len(logs.get(o, ([], []))[0]) >= min_pairs else 0.0
            used[ax] = thr[ax]
            logs[ax] = _axis_logs(pairs, ax, sub, thr[ax], dominance)
            m = lowest_mode([v for _s, v in logs[ax][0] if v is not None],
                            min_support)
            modes[ax] = m
            new[ax] = 10 ** m[0] if m else 0.0
        done = all(abs(new[a] - rate[a]) <= 1e-3 * max(new[a], 1e-12)
                   for a in AXES)
        rate = new
        thr = {ax: threshold(rate[ax]) if rate[ax] else start[ax]
               for ax in AXES}
        if done:
            break

    for ax in AXES:
        lg, xs = logs[ax]
        zeros = sum(1 for _s, v in lg if v is None)
        res = {"pairs": len(lg), "pairs_with_nothing": zeros,
               "destination_wells": len({s for s, _v in lg}),
               "min_source_reads": round(used[ax]),
               "rate": None, "ci95": None, "basis": "default", "reason": ""}
        m = modes.get(ax)
        if len(lg) < min_pairs:
            res["reason"] = (
                f"too few wells next to an abundant source: {len(lg)} "
                f"well-taxon pairs with {res['min_source_reads']:,}+ "
                f"source reads along this axis, {min_pairs} needed")
        elif (m is None and zeros >= 0.5 * len(lg)) or (
                m is not None and _mostly_empty(lg, xs, rate[ax])):
            # Wells that should show a leak mostly show nothing: the rate
            # is below what the run can resolve, and any low cluster is
            # something else. Bounded by the median pair: a consensus
            # needs ~5 reads, so a leak under that is unseen.
            bound = 5.0 / _quantile(xs, 0.5)
            res.update(basis="upper_bound", rate=bound,
                       reason=f"{zeros} of {len(lg)} informative pairs "
                       "show nothing where a leak would be visible: the "
                       f"rate is below about {bound:.1e}")
        elif m is None:
            res["reason"] = ("no clear low cluster of ratios to measure a "
                             "rate from")
        elif rate[ax] > max_rate:
            res["reason"] = (
                f"the lowest cluster of ratios is at {rate[ax]:.1e} per "
                f"pair, above the {max_rate:.0e} misreading could reach: "
                "the taxa shared along this axis are more likely the "
                "samples' own layout (wells grouped by host or site) "
                "than misassignment")
        else:
            res["rate"] = rate[ax]
            res["mode_support"] = m[1]
            byw: dict = defaultdict(list)
            for s, v in lg:
                byw[s].append(v)
            ws = sorted(byw)
            boots = []
            for _ in range(n_boot):
                pick = [v for w in (rng.choice(ws) for _w in ws)
                        for v in byw[w] if v is not None]
                bm = lowest_mode(pick, min_support)
                boots.append(10 ** bm[0] if bm else math.inf)
            boots.sort()
            lo = boots[int(0.025 * n_boot)]
            hi = boots[min(n_boot - 1, int(math.ceil(0.975 * n_boot)) - 1)]
            res["ci95"] = [lo, hi]
            # The placebo: of the pairs at or below the rate's window, how
            # many sit IN it? Where the rate is real, a well's reads of a
            # taxon it lacks track what ITS neighbours hold, so those
            # pairs pile up there; deal the wells out over the plate at
            # random and the same reads, divided by the wrong neighbours,
            # scatter. Reads that arrive some other way - a steady level
            # from the whole plate, say - track the true neighbours no
            # better than shuffled ones. Presence sits above the window
            # either way, so it is left out rather than diluting this.
            # Nothing is taken off for the other axis, which a shuffled
            # plate has no rate for.
            tight, n_low = _low_tail_share(pairs, ax, rate[ax], used[ax],
                                           dominance)
            beaten = 0
            for _ in range(n_placebo):
                pp = IndexSharing(_placebo_layout(geom, rng), counts
                                  ).pairs(min_source, taxa)
                if _low_tail_share(pp, ax, rate[ax], used[ax],
                                   dominance)[0] >= tight:
                    beaten += 1
            if n_low < min_support:
                beaten = n_placebo
            res["placebo_p"] = round((beaten + 1) / (n_placebo + 1), 4)
            res["low_tail_share"] = round(tight, 4)
            if not math.isfinite(hi) or hi > max_spread * lo:
                res.update(rate=None, reason=(
                    "the run's data do not pin the rate: its 95% interval "
                    f"spans {lo:.1e} to {hi:.1e}"))
            elif res["placebo_p"] > placebo_p:
                res.update(rate=None, reason=(
                    "the low cluster of ratios is no tighter with each "
                    "well's true neighbours than with the wells shuffled "
                    f"over the plate (placebo p {res['placebo_p']}): what "
                    "it measures does not follow index sharing"))
            else:
                res["basis"] = "fitted"
        out[ax] = res

    # A bound means the wells that should show a leak show nothing, so
    # there is nothing to check calls against: zero, as for the blanks.
    if groups:
        _lowest_group(out, geom, counts, rng, groups,
                      dict(min_source=min_source, dominance=dominance,
                           min_pairs=min_pairs, min_support=min_support,
                           max_spread=max_spread, max_rate=max_rate,
                           default_rate=default_rate, n_boot=n_boot,
                           n_placebo=n_placebo, placebo_p=placebo_p))

    out["for_calls"] = {ax: (0.0 if out[ax]["basis"] == "upper_bound"
                             else out[ax]["rate"] if out[ax]["rate"] is not None
                             else default_rate) for ax in AXES}
    out["for_calls_upper"] = {
        ax: (out[ax]["ci95"][1] if out[ax]["basis"] == "fitted"
             else out["for_calls"][ax]) for ax in AXES}
    out["for_blanks"] = {ax: (out[ax]["rate"] if out[ax]["basis"] == "fitted"
                              else 0.0) for ax in AXES}
    return out


def _lowest_group(out: dict, geom: dict, counts: dict, rng: random.Random,
                  groups: dict, kw: dict) -> None:
    """Fit each group alone; where a group's 95% interval lies wholly
    below the pooled fitted rate, take the lowest such group's fit for
    that axis. Each axis records every group's result."""
    by: dict = defaultdict(set)
    for k, g in groups.items():
        by[g].add(k)
    fits = {}
    for g in sorted(by):
        # Cheap screen before the full fit: enough source-rich taxa?
        if len(IndexSharing(geom, counts).pairs(kw["min_source"], by[g])) \
                < kw["min_pairs"]:
            continue
        fits[g] = fit_misassignment(geom, counts, rng, taxa=by[g], **kw)
    for ax in AXES:
        pooled = out[ax]
        pooled["by_group"] = {
            g: {"basis": f[ax]["basis"], "rate": f[ax]["rate"],
                "ci95": f[ax]["ci95"], "pairs": f[ax]["pairs"]}
            for g, f in fits.items()}
        if pooled["basis"] != "fitted":
            continue
        lower = [(f[ax]["rate"], g) for g, f in fits.items()
                 if f[ax]["basis"] == "fitted"
                 and f[ax]["ci95"][1] < pooled["rate"]]
        if lower:
            _r, g = min(lower)
            pooled["pooled"] = {k: pooled[k] for k in
                                ("rate", "ci95", "pairs", "placebo_p")}
            keep = pooled["by_group"]
            pooled.update({k: v for k, v in fits[g][ax].items()})
            pooled["by_group"] = keep
            pooled["from_group"] = g


def _wls_through_origin(cells: list, fixed: dict) -> dict:
    """Rates per axis for unused cells, reads ~ row x X_row + col x X_col,
    by least squares weighted for Poisson counts (1 / (reads + 1)) and
    held non-negative. An axis in `fixed` is taken off first and not
    fitted; an axis whose sources are all zero is not fitted (None)."""
    w = [1.0 / (n + 1) for n, _x, _z in cells]
    free = [ax for ax in AXES if ax not in fixed
            and any((x if ax == "row" else z) for _n, x, z in cells)]
    y = [n - fixed.get("row", 0.0) * x - fixed.get("col", 0.0) * z
         for n, x, z in cells]
    col = {"row": [x for _n, x, _z in cells], "col": [z for _n, _x, z in cells]}

    def sse(rates):
        return sum(wi * (yi - sum(rates[a] * col[a][k] for a in rates)) ** 2
                   for k, (wi, yi) in enumerate(zip(w, y)))

    def one(ax):
        sxx = sum(wi * v * v for wi, v in zip(w, col[ax]))
        sxy = sum(wi * v * yi for wi, v, yi in zip(w, col[ax], y))
        return max(0.0, sxy / sxx) if sxx else 0.0

    out: dict = {ax: None for ax in AXES if ax not in fixed}
    if len(free) == 1:
        out[free[0]] = one(free[0])
    elif len(free) == 2:
        a, b = col["row"], col["col"]
        saa = sum(wi * v * v for wi, v in zip(w, a))
        sbb = sum(wi * v * v for wi, v in zip(w, b))
        sab = sum(wi * u * v for wi, u, v in zip(w, a, b))
        say = sum(wi * u * yi for wi, u, yi in zip(w, a, y))
        sby = sum(wi * v * yi for wi, v, yi in zip(w, b, y))
        det = saa * sbb - sab * sab
        r = c = -1.0
        if det > 1e-9 * saa * sbb:
            r, c = (say * sbb - sby * sab) / det, (sby * saa - say * sab) / det
        if r >= 0 and c >= 0:
            out.update(row=r, col=c)
        else:
            # the constraint binds: the better of each axis alone
            cand = [{"row": one("row"), "col": 0.0},
                    {"row": 0.0, "col": one("col")}]
            out.update(min(cand, key=sse))
    return out


def sentinel_rates(matrix: list[dict], fitted: dict | None = None,
                   n_boot: int = 200, seed: int = 1) -> dict | None:
    """The direct measurement, where the plate has unused combinations.

    A read on an index pair no well claims can only be misassignment, so
    its count against the reads of its row (same i5) and its column (same
    i7) measures the per-pair rates outright. `matrix` is the
    demultiplexer's index_matrix.tsv. Each axis is measured with the
    other taken off at its rate in `fitted` (measured from the calls)
    where there is one; otherwise the two are fitted together over the
    unused cells by non-negative least squares. Unused cells in
    a part-filled column share that column's sources while their rows'
    differ, so the two separate - a median of row ratios alone ran 20-75%
    high on such a plate. The 95% interval is from resampling the cells.
    Counted over every demultiplexed read, classified or not, so it can
    differ somewhat from the fitted rates, which are over classified
    reads. None without unused combinations.
    """
    rows: dict = defaultdict(int)
    cols: dict = defaultdict(int)
    raw = []
    for r in matrix:
        n = int(r.get("reads") or 0)
        if r.get("assigned") == "yes":
            rows[r["index_i5"]] += n
            cols[r["index_i7"]] += n
        else:
            raw.append((r["index_i5"], r["index_i7"], n))
    if not raw:
        return None
    cells = [(n, rows[a], cols[b]) for a, b, n in raw]
    known = {ax: v for ax, v in (fitted or {}).items() if v and v > 0}

    def measure(cs):
        # each axis with the OTHER held at its fitted rate where there is
        # one, else both together
        joint = None
        res = {}
        for ax in AXES:
            other = "col" if ax == "row" else "row"
            if other in known:
                res[ax] = _wls_through_origin(
                    cs, {other: known[other]}).get(ax)
            else:
                if joint is None:
                    joint = _wls_through_origin(cs, {})
                res[ax] = joint.get(ax)
        return res

    est = measure(cells)
    rng = random.Random(seed)
    boots: dict = {ax: [] for ax in est}
    for _ in range(n_boot):
        b = measure([rng.choice(cells) for _c in cells])
        for ax in est:
            if b.get(ax) is not None:
                boots[ax].append(b[ax])
    out: dict = {"unused_combinations": len(cells),
                 "reads_on_unused": sum(n for n, _x, _z in cells)}
    for ax in AXES:
        idx = 1 if ax == "row" else 2
        n_cells = sum(1 for c in cells if c[idx])
        if est.get(ax) is None or not n_cells:
            out[ax] = None
            continue
        bs = sorted(boots[ax])
        ci = ([bs[int(0.025 * len(bs))],
               bs[min(len(bs) - 1, int(math.ceil(0.975 * len(bs))) - 1)]]
              if bs else None)
        out[ax] = {"rate": est[ax], "ci95": ci, "cells": n_cells,
                   "other_axis": ("taken off at the fitted rate"
                                  if ("col" if ax == "row" else "row")
                                  in known else "fitted together")}
    return out


# --------------------------------------------------------------------------
# span distances (releases scored over their records' span)
# --------------------------------------------------------------------------


def _myers_columns(ref: str, read: str):
    """Myers' bit-parallel edit distance of `ref` against `read`, the
    read's start free (no boundary bit). Yields, after each read base,
    (score, pv, mv): the cost of all of ref ending there, and the column's
    vertical deltas (bit i: D[i+1] - D[i] is +1 in pv, -1 in mv)."""
    m = len(ref)
    peq: dict[str, int] = {}
    for i, c in enumerate(ref):
        peq[c] = peq.get(c, 0) | (1 << i)
    full, top = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score = full, 0, m
    for c in read:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & full)
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = mh | (~(xv | ph) & full)
        mv = ph & xv
        yield score, pv, mv


def infix_distance(ref: str, read: str) -> int:
    """Edits to place all of `ref` somewhere inside `read`, the read's
    bases either side free: a whole insert scored over its record's span
    (as nemabiome-its2-refdb's bin/build.py)."""
    if not ref:
        return 0
    best = len(ref)
    for score, _pv, _mv in _myers_columns(ref, read):
        best = min(best, score)
    return best


def anchored_distance(ref: str, read: str, min_span: int = 1) -> tuple:
    """(edits, span) of a read window against the START of a record.

    The record's first base is aligned (a partial record that starts late
    leaves the read's bases before it free). Then either the window is
    read to its last base and the record past it is free, or - a window
    longer than the record, R1 reading into the 28S stretch at 2x250 -
    the whole record is aligned and the window past it is free. span is
    the record bases aligned. Of the possible end points, the one with
    the best identity over at least min_span bases. For a window against
    a record's END, pass both reversed: R2's 28S stretch, past the
    record's end, is then the free overhang. (inf, 0) when nothing fits.
    """
    m = len(ref)
    if not m or not read:
        return float("inf"), 0
    pv = mv = 0
    whole = m
    for s, pv, mv in _myers_columns(ref, read):
        whole = min(whole, s)
    best, d = None, 0
    for i in range(1, m + 1):
        d += (pv >> (i - 1) & 1) - (mv >> (i - 1) & 1)
        # The record ending inside the window: all of it, read past free.
        if i == m:
            d = min(d, whole)
        # The best identity; on a tie, the longer span.
        if i >= min_span and (best is None or
                              d * best[1] <= best[0] * i):
            best = (d, i)
    return best if best else (float("inf"), 0)


# --------------------------------------------------------------------------
# database mode: helpers every protocol may use
# --------------------------------------------------------------------------
#
# Database mode has three layers (docs/CONTRACTS.md): a platform adapter
# gets each sample's reads, a protocol (bin/protocols/) turns them into
# calls, and synthesis and the report are shared. What follows is what
# the protocols and the read-type routes share - minimap2, primer
# trimming, the LCA and identity floors - so that no protocol needs to
# import another.

# minimap2's index parameters for the database. k=19, w=10 gave the same
# accuracy as the map-ont default k=15 on the leave-out benchmark and ran
# 1.6x faster; the index builds in 1-2 s, so it is built once per run.
INDEX_ARGS = ["-x", "map-ont", "-k", "19", "-w", "10"]


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


def hit_identity(hit: dict) -> float:
    """Matching bases over the alignment block PLUS any of the query left
    unaligned at either end.

    minimap2's own ratio (column 10 over column 11) ignores the ends it
    clipped, so a consensus aligned over only part of its length would
    look as close as one aligned end to end. A consensus and a record are
    both primer-to-primer inserts, so an unaligned end is a real
    difference and is charged as one."""
    clip = hit["qs"] + (hit["qlen"] - hit["qe"])
    return hit["match"] / max(1, hit["alen"] + clip)


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


def trim_reads(path, primers) -> tuple[dict[str, str], Counter]:
    """Oriented, primer-trimmed inserts by read id, and what was lost.

    The release's records are primer-trimmed inserts, so reads must be
    too: identity is measured over the insert, and a primer left on one
    side and not the other would be charged as a difference."""
    fate: Counter = Counter()
    return dict(iter_trimmed(path, primers, fate)), fate


def iter_trimmed(path, primers, fate: Counter):
    """trim_reads() one read at a time, so main() can write them to disk
    without holding a deep well's reads in memory; fate is filled in as
    they pass. Only the ids are kept, to make duplicates unique. `path`
    may be a list: a sample delivered as several FASTQs (a predemultiplexed
    run) is read as one."""
    out: set[str] = set()
    paths = path if isinstance(path, (list, tuple)) else [path]
    for rid, seq, _q in open_reads_many(paths):
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


def orient_pair(r1, r2, ps):
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


def consensus_hits(mm, index: str, tmp: str, seqs: dict) -> dict:
    """{"c<key>": [hits]} for sequences {key: sequence} against the
    database, every secondary hit within 90% of the best kept: what a
    consensus is classified from."""
    kfa = os.path.join(tmp, "classify.fa")
    write_fasta(kfa, ((f"c{k}", s) for k, s in seqs.items()))
    hits: dict[str, list] = defaultdict(list)
    for a in mm.map(index, kfa, preset=INDEX_ARGS,
                    extra=("--secondary=yes", "-N", "1000",
                           "-p", "0.9")):
        hits[a["q"]].append(a)
    return hits


# Window mapping, with minimap2's short-read preset.
SR = ("-x", "sr")


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


def pair_hits(mm, index: str, tmp: str, centroids: list,
              max_secondary: int) -> tuple[dict, dict, tuple]:
    """(hits, paths, extra): each denoised pair's windows mapped to the
    database (map_windows), R1 and R2 windows as "r1" and "r2", whole
    inserts as "full", queries named p<k>; and the FASTAs and options
    used, for a protocol that maps some of them again."""
    f1, f2, ff = (os.path.join(tmp, n) for n in ("r1.fa", "r2.fa", "full.fa"))
    write_fasta(f1, ((f"p{k}", x) for k, (x, y, _n) in enumerate(centroids) if y))
    write_fasta(f2, ((f"p{k}", y) for k, (x, y, _n) in enumerate(centroids) if y))
    write_fasta(ff, ((f"p{k}", x) for k, (x, y, _n) in enumerate(centroids) if not y))
    # -p 0.5, not the nanopore route's 0.9: each window is kept with
    # every record that scores within half of its best, so that a
    # record second-best for BOTH windows - which is what a taxon
    # missing from the database often looks like - is still there to
    # be joined.
    extra = ("--secondary=yes", "-N", str(max_secondary), "-p", "0.5")
    paths = {"r1": f1, "r2": f2, "full": ff}
    return map_windows(mm, index, paths, extra, SR), paths, extra


# --------------------------------------------------------------------------
# contract B: what a protocol hands to synthesis (docs/CONTRACTS.md)
# --------------------------------------------------------------------------

# <sample>.clusters.tsv: these columns, in this order, from every
# protocol; then the protocol's own declared extra columns
# (Protocol.extra_columns), if any.
CLUSTER_COLUMNS = [
    "sample_id", "cluster", "reads", "rel_abundance", "rank", "taxon",
    "label", "identity", "coverage", "lineage", "best_hit",
    "best_hit_lineage", "near_hits", "lca_rank", "db_resolution",
    "db_margin", "origin", "merged", "chimera_of", "length",
    "database_version", "consensus"]
# <sample>.taxa.tsv: exactly these.
TAXA_COLUMNS = ["sample_id", "rank", "taxon", "label", "reads",
                "rel_abundance", "clusters", "best_identity", "lineage",
                "database_version"]
# <sample>.assign_db.json: at least these keys; "reads" at least
# ASSIGN_READ_FATES, which balance (docs/CONTRACTS.md).
ASSIGN_SUMMARY_KEYS = ["sample_id", "mode", "database", "minimap2",
                       "elapsed_s", "reads", "reads_note", "clusters",
                       "haplotype_splits", "chimeras", "settings",
                       "rank_reached"]
ASSIGN_READ_FATES = ["reads_in", "no_primer", "too_short", "clustered",
                     "unclustered"]
# synthesis/db_calls.tsv, the shared layer's own output: these columns,
# then any of DB_CALL_OPTIONAL that a protocol's clusters carry.
DB_CALL_COLUMNS = [
    "run_id", "sample_id", "sample_type", "rank", "taxon", "label", "reads",
    "rel_abundance", "call", "limit_reads", "limit_basis",
    "misassign_expected", "misassign_excess", "identity", "clusters",
    "lineage", "group", "group_reads", "group_share", "primer_fit",
    "margin", "rank_meaning", "database_version"]
DB_CALL_OPTIONAL = ["nearest", "next", "grade"]


def write_assignment(outdir: str, sid: str, rows: list[dict],
                     extra_columns, version: str, summary: dict) -> None:
    """Write one sample's contract-B files: clusters.tsv (rows, most
    reads first, with the protocol's extra columns after the required
    ones), consensus.fasta, taxa.tsv and assign_db.json."""
    cols = CLUSTER_COLUMNS + list(extra_columns)
    with open(os.path.join(outdir, f"{sid}.clusters.tsv"), "w",
              encoding="utf-8") as fh:
        w = writer(fh, cols)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(outdir, f"{sid}.consensus.fasta"), "w",
              encoding="utf-8") as fh:
        for r in rows:
            fh.write(f">{r['cluster']} reads={r['reads']} rank={r['rank'] or 'none'} "
                     f"taxon={r['taxon'] or 'none'} identity={r['identity']}\n"
                     f"{r['consensus']}\n")

    # Per taxon at the rank each call reached. Two clusters of one
    # species (haplotypes) are one row; a species call and a genus call
    # in the same genus are NOT merged - they are different statements.
    total = sum(r["reads"] for r in rows)
    agg: dict = {}
    for r in rows:
        key = (r["rank"], r["taxon"], r["label"])
        a = agg.setdefault(key, {"reads": 0, "clusters": 0, "identity": 0.0,
                                 "lineage": r["lineage"]})
        a["reads"] += r["reads"]
        a["clusters"] += 1
        a["identity"] = max(a["identity"], float(r["identity"]))
    with open(os.path.join(outdir, f"{sid}.taxa.tsv"), "w",
              encoding="utf-8") as fh:
        w = writer(fh, TAXA_COLUMNS)
        w.writeheader()
        for (rank, taxon, label), a in sorted(agg.items(),
                                              key=lambda kv: -kv[1]["reads"]):
            w.writerow({"sample_id": sid, "rank": rank, "taxon": taxon,
                        "label": label, "reads": a["reads"],
                        "rel_abundance": f"{a['reads'] / total:.6f}" if total else "",
                        "clusters": a["clusters"],
                        "best_identity": f"{a['identity']:.4f}",
                        "lineage": a["lineage"], "database_version": version})
    with open(os.path.join(outdir, f"{sid}.assign_db.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
