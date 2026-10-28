#!/usr/bin/env python3
"""
Shared primitives for nf-helminth-amplicon.

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
import os
import re
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

    Contexts the argument does not cover - longer than 64bp, too short to
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
            pieces = self._cut(ctx, k + 2) if len(ctx) <= 64 else None
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
        made it a hundred times slower. A piece with more than 64
        spellings (four Ns) is not worth indexing, and its context is
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
                if len(spellings) > 64:
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

# What the pipeline reads. A release missing any of these cannot be used.
REFDB_REQUIRED = ["manifest.json", "refdb.fasta", "taxonomy.tsv",
                  "anchor_lookalikes.tsv", "resolution_nanopore.tsv"]
# The rest of the contract: not read by the nanopore route, but their
# absence means the release was not built the agreed way.
REFDB_EXPECTED = ["members.tsv", "conflicts.tsv", "primer_fit.tsv"]

LOOKALIKE_COLUMNS = ["seq_id", "anchor", "strand", "position", "edits",
                     "context"]


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


def load_taxonomy(path: str) -> dict[str, dict[str, str]]:
    """id -> {rank: name} over RANKS, blanks kept as ''."""
    out = {}
    with open_text(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["id"]] = {k: (r.get(k) or "").strip() for k in RANKS}
    return out
