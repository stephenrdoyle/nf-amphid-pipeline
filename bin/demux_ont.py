#!/usr/bin/env python3
"""
Nanopore demultiplexing for nf-helminth-amplicon.

The nanopore runs sequence Illumina libraries, so this is in-line barcode
demultiplexing on reads carrying full Nextera adapter structure - not
kit-based ONT barcoding. Four things shape the implementation, each of them
a decision that could reasonably have gone the other way:

1. ANCHOR, NEVER SCAN. The i5 and i7 index sets are only one substitution
   apart at their closest (S502 CTCTCTAT vs N707 CTCTCTAC), so an unanchored
   search cannot tell which end of a molecule it is reading. The barcode
   window is located by the adapter sequence abutting it, and only that
   window is tested.

2. SUBSTITUTIONS ONLY, INDELS DISABLED. The index sets have minimum Hamming
   distance 4 but minimum *edit* distance 2 (i7 N721 vs N723). One
   substitution is safely correctable; one indel is not. Since indels are
   nanopore's dominant error, allowing them would misassign.

3. THE TANDEM DUPLICATION IS INDEL INSURANCE, NOT EXTRA TOLERANCE. An indel
   makes one copy unmatchable and the other carries the read; a frame slip
   produces two copies that disagree, which is what detects it. With
   barcode_repeats=1 there is no second observation, and the rotational
   collisions in these sets become live misassignment routes.

4. TIERS, NOT A SINGLE PASS. Both ends agreeing is the confident case;
   single-ended reads are assigned but flagged, so the cost of accepting
   them can be measured rather than assumed.

Self-contained: standard library only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    DEFAULTS,
    sanitise_id,
    open_reads_many,
    resolve_read_files,
    find_anchor,
    find_anchor_cost,
    hamming,
    load_assay,
    open_reads,
    read_run_metadata,
    read_table,
    revcomp,
)

# --------------------------------------------------------------------------
# construct defaults, overridden by the assay config
# --------------------------------------------------------------------------

# Assignment tiers. A is the confident case; D is the residue.
TIERS = ("A", "B", "C", "D", "E")

TIER_MEANING = {
    "A": "both ends agree and match a samplesheet pair",
    "B": "one end only - the other truncated or unreadable",
    "C": "ends disagree - chimeric",
    "D": "no confident match",
    "E": "internal adapter - the read is more than one molecule",
}


# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# barcode matching
# --------------------------------------------------------------------------


@dataclass
class BarcodeSet:
    """The declared barcodes for one end, plus the tolerance the set itself
    can support. Tolerance is derived from the set rather than supplied:
    a number someone typed cannot know the set's minimum distance."""

    by_seq: dict[str, str] = field(default_factory=dict)  # sequence -> name
    length: int = 0
    max_mismatch: int = 1
    # How much further the runner-up must be for a match beyond
    # max_mismatch to be accepted, and how far the best may be. 0
    # disables the margin rule entirely and restores tolerance-only
    # behaviour.
    margin: int = 2
    margin_max_dist: int = 2

    @classmethod
    def build(cls, seqs: set[str], requested: int = 1,
              margin: int = 2, margin_max_dist: int = 2) -> "BarcodeSet":
        seqs = {s for s in seqs if s}
        if not seqs:
            return cls()
        lengths = {len(s) for s in seqs}
        length = max(lengths) if lengths else 0
        by_seq = {s: s for s in seqs}

        safe = requested
        ordered = sorted(seqs)
        for i, a in enumerate(ordered):
            for b in ordered[i + 1 :]:
                if len(a) != len(b):
                    continue
                # k mismatches are only safely correctable if every pair is
                # more than 2k apart.
                safe = min(safe, (hamming(a, b) - 1) // 2)
        return cls(by_seq=by_seq, length=length, max_mismatch=max(0, safe),
                   margin=margin, margin_max_dist=margin_max_dist)

    def match(self, window: str) -> str | None:
        """Exact first, then the nearest - by tolerance, or by margin.

        Two acceptance rules, and the second is what makes this worth
        more than max_mismatch alone.

        TOLERANCE is a property of the SET: k substitutions are safely
        correctable only when every pair is more than 2k apart, so
        max_mismatch is derived from the set's minimum distance rather
        than chosen. On an index set whose closest pair is 4 apart that
        caps at 1, whatever the data looks like.

        MARGIN is a property of the READ. A window 2 away from one
        declared index and 5 from every other is unambiguous - not
        because the set is widely spaced, but because this particular
        window is not near anything else. Accepting it introduces no
        misassignment that the set's geometry would otherwise have
        prevented, because the alternative was not close either.

        Measured on the reference run: the margin rule recovers 5.0% of
        failed i5 ends and 18.2% of i7, all with a runner-up at least 2
        further away. The ambiguous remainder is declined, as before -
        guessing would be worse than declining, and that has not
        changed.
        """
        if window in self.by_seq:
            return window
        if len(window) != self.length:
            return None
        if self.max_mismatch < 1 and self.margin < 1:
            return None

        ds = sorted((hamming(window, c), c) for c in self.by_seq)
        if not ds:
            return None
        best_d, best = ds[0]
        runner_up = ds[1][0] if len(ds) > 1 else 99

        # A tie is never a match, under either rule.
        if best_d == runner_up:
            return None
        if best_d <= self.max_mismatch:
            return best
        if (self.margin >= 1 and best_d <= self.margin_max_dist
                and runner_up - best_d >= self.margin):
            return best
        return None

    def candidates(self, window: str, max_dist: int = 4
                   ) -> tuple[tuple[str, ...], int]:
        """Every declared barcode tied for nearest, and that distance.

        What `match` throws away. `match` resolves one end in isolation,
        so a window equidistant from two declared barcodes is
        unresolvable and is correctly declined. But the read has a
        SECOND end, and if the plate is SPARSE the two together can
        resolve what neither can alone: a window tied between two i7
        candidates still resolves when only one of them forms a declared
        well with the i5 this read carries. That is MiniBar's trick.

        How much it can deliver is a property of the PLATE, not of this
        function, and it is worth knowing before reading any number it
        produces. On the reference titration plate the grid is 8 i5 x 12
        i7 = 96 combinations with 94 declared - 97.9% saturated. A tie
        there almost always has both candidates forming real wells, so
        the trick has nearly no room to work. MiniBar's own example
        plate was three samples drawn from a much larger index space.

        Returns the TIED set, not the k nearest: a candidate further
        away than another is already excluded by distance and needs no
        combination to rule it out. A set of ONE therefore means nothing
        was tied, and the caller must not treat it as a combination
        win - see combination_check(). Empty when nothing is within
        max_dist, which is the genuinely unreadable case.
        """
        if not self.by_seq or len(window) != self.length:
            return (), -1
        ds = sorted((hamming(window, c), c) for c in self.by_seq)
        if not ds or ds[0][0] > max_dist:
            return (), -1
        best_d = ds[0][0]
        return tuple(c for d, c in ds if d == best_d), best_d


@dataclass
class EndCall:
    barcode: str | None = None
    copies_seen: int = 0
    copies_agreed: bool = False
    anchor_pos: int = -1
    anchor_end: int = -1
    window_start: int = -1
    # Why this end did or did not yield a barcode. The tier a read lands in
    # says how much evidence survived; this says what destroyed it, which is
    # the part you need to fix a run where most reads fail. Costs a string
    # assignment per end and is read by bin/diagnose_demux.py.
    reason: str = "ok"
    anchor_edits: int = -1
    # The declared barcodes tied for nearest, when this end was declined
    # and the window was at least readable. Populated only on failure -
    # a resolved end has nothing to disambiguate - and used by the
    # combination check, which asks whether the OTHER end narrows these
    # to one declared well. Measurement only: nothing here changes a
    # barcode, a tier or a sample.
    candidates: tuple[str, ...] = ()
    candidate_dist: int = -1


def call_end(
    read: str,
    anchor: str,
    bc: BarcodeSet,
    repeats: int,
    search_window: int,
    barcode_before_anchor: bool,
    slack: int = 3,
    anchor_max_edits: int = 2,
) -> EndCall:
    """Locate the anchor, read the tandem barcode copies beside it, and
    require them to agree.

    Which side matters, and the two ends differ. The construct is

        P5 - i5 x n - [i5 anchor]...amplicon...[i7 anchor] - i7 x n - rc(P7)

    so the i5 copies sit immediately BEFORE their anchor (the anchor is the
    5' end of the PCR1 forward overhang) while the i7 copies sit immediately
    AFTER theirs (the anchor is the 3' end of the reverse-complemented PCR1
    reverse overhang). Assuming one side for both silently loses every i7.

    The copies are matched independently rather than concatenated. Joining
    them pools the error budget: two mismatches across 16bp, concentrated in
    one half, can push a read onto the wrong barcode. Two independent
    observations that agree is stronger evidence than one longer noisy one.
    """
    out = EndCall()
    if not anchor or not bc.length:
        out.reason = "no_anchor_configured"
        return out

    # The anchor is searched only within a window from the read end it
    # belongs to. The construct puts it at a known offset, and scanning the
    # whole read would invite a chance match elsewhere.
    window = read[:search_window]
    exact = window.find(anchor)
    if exact >= 0:
        a_start, a_end = exact, exact + len(anchor)
        out.anchor_edits = 0
    else:
        hit = find_anchor_cost(window, anchor, anchor_max_edits)
        if hit is None:
            out.reason = "anchor_not_found"
            return out
        a_start, a_end, out.anchor_edits = hit
    out.anchor_pos = a_start
    out.anchor_end = a_end

    span = bc.length * repeats
    # Measured from the alignment boundary, so an indel inside the anchor
    # does not displace the barcode window.
    start = a_start - span if barcode_before_anchor else a_end
    if start + slack < 0 or start + span - slack > len(read):
        # the anchor was found but the barcode runs off the end of the read
        out.reason = "barcode_off_read"
        return out
    out.window_start = start

    # Each copy is sought at its expected offset PLUS OR MINUS a little
    # slack, rather than at a fixed position. An indel between the anchor and
    # the barcode shifts every downstream base, so a fixed window turns one
    # indel into two failed copies - which is exactly the case the tandem
    # duplication exists to survive. Slack is small and anchored, so this is
    # still a handful of tested windows rather than a scan.
    calls = []
    # The best window per copy, kept whether or not it matched, so a
    # declined end can still report what it was nearest to. Taken at the
    # EXPECTED offset only: the shift search exists to find a match, and
    # reusing it here would let the candidate set drift to a different
    # position than the one that failed.
    windows = []
    for k in range(repeats):
        expect = start + k * bc.length
        if 0 <= expect and expect + bc.length <= len(read):
            windows.append(read[expect : expect + bc.length])
        # No distance ceiling here: bc.match has already decided what is
        # acceptable, by tolerance OR by margin. This used to start at
        # max_mismatch + 1, which silently discarded every hit the margin
        # rule accepted beyond max_mismatch - so the rule was measured,
        # reported in the summary as active, and never applied.
        best, best_d = None, None
        for shift in range(-slack, slack + 1):
            a = expect + shift
            if a < 0 or a + bc.length > len(read):
                continue
            hit = bc.match(read[a : a + bc.length])
            if hit is None:
                continue
            d = hamming(read[a : a + bc.length], hit)
            if best_d is None or d < best_d:
                best, best_d = hit, d
        calls.append(best)

    def record_candidates() -> None:
        """Nearest declared barcodes over the copies, for the combination
        check. Pooled across copies by taking the closest: the copies are
        the same 8bp in the same context, so the better-read one is the
        better evidence about which barcodes are in contention."""
        pooled: tuple[str, ...] = ()
        pooled_d = -1
        for w in windows:
            cand, d = bc.candidates(w)
            if not cand:
                continue
            if pooled_d < 0 or d < pooled_d:
                pooled, pooled_d = cand, d
            elif d == pooled_d:
                pooled = tuple(sorted(set(pooled) | set(cand)))
        out.candidates, out.candidate_dist = pooled, pooled_d

    seen = [c for c in calls if c]
    out.copies_seen = len(seen)
    if not seen:
        out.reason = "no_copy_matched"
        record_candidates()
        return out

    if len(set(seen)) == 1:
        out.barcode = seen[0]
        out.copies_agreed = len(seen) == repeats
        out.reason = "ok" if out.copies_agreed else "ok_single_copy"
    else:
        # copies disagree: a frame slip from an upstream indel produces
        # exactly this, and guessing which copy is right would defeat the
        # point of duplicating them
        out.barcode = None
        out.copies_agreed = False
        out.reason = "copies_disagreed"
        record_candidates()
    return out


# --------------------------------------------------------------------------
# per-read assignment
# --------------------------------------------------------------------------


@dataclass
class ReadCall:
    read_id: str
    tier: str
    sample_id: str | None
    i5: str | None
    i7: str | None
    orientation: str
    insert_start: int = -1
    insert_end: int = -1
    # Why each end did or did not yield a barcode, for the orientation
    # that was kept. The tier says how much evidence survived; these say
    # what destroyed it - and that is the difference between "your PCR2
    # is inefficient" and "the barcode matcher is too strict". A run
    # that comes back half single-ended cannot be diagnosed from tiers
    # alone, and re-reading 30M reads with a separate tool to find out
    # is not a reasonable answer.
    i5_reason: str = "ok"
    i7_reason: str = "ok"
    # Verdict of the combination check, which costs nothing and decides
    # nothing. See combination_check().
    combo: str = "not_applicable"
    combo_sample: str | None = None
    combo_dist: int = -1


COMBO_MEANING = {
    "not_applicable": "both ends resolved, or neither - no tie to narrow",
    "unreadable": "the tied end had no declared index within reach, so "
                  "there was no candidate set to narrow",
    "no_declared": "none of the tied candidates forms a declared well with "
                   "the resolved end, so declining this read cost nothing",
    "resolved": "the end was genuinely TIED between two or more declared "
                "indices, and exactly one of them forms a declared well "
                "with the end that did resolve. This is the safe recovery: "
                "the combination supplied the evidence, not a looser "
                "distance, so no misassignment risk is added",
    "distance_only": "there was only ONE candidate, so nothing was tied and "
                     "the combination contributed nothing. The read was "
                     "declined because that candidate sits further away "
                     "than the index set safely supports, and accepting it "
                     "would be a tolerance increase wearing a combination "
                     "check as a disguise. Counted separately, and by "
                     "distance, because it is not the same offer",
    "still_ambiguous": "more than one tied candidate forms a declared "
                       "well, so declining it was right",
}


END_REASON_MEANING = {
    "ok": "barcode read, and the tandem copies agreed",
    "ok_single_copy": "barcode read from one copy; the other was absent "
                      "or unreadable, so the replicate could not confirm it",
    "anchor_not_found": "no adapter anchor at this end - the molecule most "
                        "likely never carried this index, so this is PCR2 "
                        "efficiency rather than a matching problem",
    "barcode_off_read": "the anchor was found but the barcode runs past the "
                        "end of the read - the read is truncated here",
    "no_copy_matched": "the barcode window was read but matched no declared "
                       "index within tolerance; the index is present, so "
                       "this is a matching or index-set problem",
    "copies_disagreed": "the tandem copies read as different indices, which "
                        "is the frame slip the duplication exists to catch",
    "no_anchor_configured": "the assay defines no anchor for this end",
}


def has_internal_adapter(
    s: str, start: int, end: int, anchors: tuple[str, ...], max_edits: int
) -> bool:
    """Is there another molecule's end between this read's two ends?

    A read is assumed to be one molecule, and the two indices found at its
    ends are assumed to belong to the same one. Ligation does not respect
    that assumption. Two amplicons joined head to tail present the first
    one's i5 at the 5' end and the SECOND one's i7 at the 3' - a pair of
    real, declared indices that very probably belongs to some well on the
    plate, just not to either molecule. It demultiplexes cleanly into the
    wrong sample and there is nothing downstream that could notice.

    The junction is the evidence. A clean amplicon's interior is insert;
    it contains no adapter. So finding an anchor between the two outer
    anchors means the read spans more than one molecule and its ends
    cannot be trusted as a pair.

    Both strands are searched, because the join may be head to head.
    """
    if start < 0 or end < 0 or end <= start:
        return False
    interior = s[start:end]
    return any(find_anchor(interior, p, max_edits) is not None for p in anchors)


def assign_read(
    read_id: str,
    seq: str,
    cfg: dict,
    i5_set: BarcodeSet,
    i7_set: BarcodeSet,
    pair_to_sample: dict[tuple[str, str], str],
    repeats: int,
    search_window: int,
    accept_single_end: bool,
    junction_anchors: tuple[str, ...] = (),
    anchor_max_edits: int = 2,
) -> ReadCall:
    """Try both orientations; keep the one that resolves more ends.

    Nanopore reads arrive in either orientation, and the construct is not
    palindromic, so the i5 end is only at the start on one strand.
    """
    best: ReadCall | None = None
    for orient, s in (("+", seq), ("-", revcomp(seq))):
        five = call_end(
            s, cfg["i5_inner_anchor"], i5_set, repeats, search_window,
            barcode_before_anchor=True,
        )
        # The i7 anchor sits near the 3' end, so search the tail - and its
        # barcode follows the anchor rather than preceding it.
        tail_from = max(0, len(s) - search_window)
        three = call_end(
            s[tail_from:], cfg["i7_inner_anchor"], i7_set, repeats,
            search_window, barcode_before_anchor=False,
        )

        resolved = (1 if five.barcode else 0) + (1 if three.barcode else 0)
        call = _classify(
            read_id, five, three, orient, pair_to_sample, accept_single_end
        )
        # Trim past the whole overhang, not just the anchor: the mosaic end
        # is shared between the two overhangs, so it cannot identify an end,
        # but it is adapter all the same and does not belong in the insert.
        five_tail = max(0, len(cfg["pcr1_overhang_fwd"]) - len(cfg["i5_inner_anchor"]))
        three_head = max(0, len(cfg["pcr1_overhang_rev"]) - len(cfg["i7_inner_anchor"]))
        if five.anchor_end >= 0:
            call.insert_start = five.anchor_end + five_tail
        if three.anchor_pos >= 0:
            call.insert_end = tail_from + three.anchor_pos - three_head

        # Only a read that resolved both ends is at risk of a fabricated
        # pair, and only those pay for the interior scan. A read that
        # already failed has nothing to protect.
        if junction_anchors and call.tier == "A":
            if has_internal_adapter(s, call.insert_start, call.insert_end,
                                    junction_anchors, anchor_max_edits):
                call.tier = "E"
                call.sample_id = None

        call.i5_reason = five.reason
        call.i7_reason = three.reason
        # Measurement only. The verdict is recorded and counted; the
        # tier, the sample and the output files are untouched. Whether
        # to act on it is a decision that needs this number first.
        call.combo, call.combo_sample, _ = combination_check(
            five, three, pair_to_sample)
        # The distance matters only for distance_only: "unique nearest
        # at 2" and "unique nearest at 4" are different risks, and one
        # number covering both would hide that.
        tied_end = three if five.barcode else five
        call.combo_dist = (tied_end.candidate_dist
                           if call.combo == "distance_only" else -1)
        if best is None or resolved > _resolved(best):
            best = call
            best._resolved = resolved  # type: ignore[attr-defined]
    return best  # type: ignore[return-value]


def _resolved(call: ReadCall) -> int:
    return getattr(call, "_resolved", 0)


def combination_check(
    five: EndCall,
    three: EndCall,
    pair_to_sample: dict[tuple[str, str], str],
) -> tuple[str, str | None, str | None]:
    """Does the resolved end narrow the other end's tie to one well?

    The question MiniBar answers and this demultiplexer currently does
    not. Each end is resolved in isolation, so a window equidistant from
    two declared barcodes is declined - correctly, for one end alone.
    But the plate is sparse. 16 i5 x 24 i7 is 384 possible combinations
    and 94 were used, so roughly three in four combinations are not a
    well at all, and a tie between two candidates is more likely than
    not to contain only one that pairs with the other end.

    This MEASURES the prize and changes nothing. It returns a verdict,
    and the caller counts verdicts. Deciding whether to act on it needs
    a number first, and the number has never been taken:

      not_applicable  neither end tied, or both failed - nothing to do
      unreadable      the tied end had no candidate within reach
      no_declared     none of the tied candidates forms a declared well
                      with the resolved end, so the read was never
                      assignable and declining it lost nothing
      resolved        exactly one does - this read is recoverable, and
                      recoverable with no new misassignment risk,
                      because the alternatives are not wells
      still_ambiguous more than one does - correctly declined

    On `no_declared` being good news: it means the decline was free. A
    large no_declared count would say the current rule is already
    taking everything there is to take.
    """
    i5, i7 = five.barcode, three.barcode
    # Exactly one end resolved, and the other must have candidates to
    # narrow. Both-resolved is already tier A or C; both-failed has no
    # fixed end to pivot on and the declared set cannot narrow a tie
    # against another tie without inventing evidence.
    if bool(i5) == bool(i7):
        return "not_applicable", None, None
    tied = three if i5 else five
    if not tied.candidates:
        return "unreadable", None, None

    # THE SPLIT THAT MATTERS, and the first version of this got it
    # wrong by not making it.
    #
    # A candidate set of ONE is not a tie. `match` declines a unique
    # nearest candidate whenever its distance exceeds what the set's
    # geometry safely supports and the margin rule does not rescue it.
    # The combination check then "confirms" the single candidate forms a
    # declared well and contributes NOTHING - the only thing that
    # recovered the read was looking further out. That is a tolerance
    # increase, and on this index set a tolerance increase is exactly
    # what is unsafe: minimum Hamming 4 caps safe substitutions at 1,
    # and i7 N721 vs N723 are edit distance 2 apart.
    #
    # Reporting both as "resolved" put the safe mechanism and the unsafe
    # one in one number. On the reference plate that number was 25% of
    # tier B and almost none of it was the combination.
    if len(tied.candidates) < 2:
        c = tied.candidates[0]
        pair = (i5, c) if i5 else (c, i7)
        if pair not in pair_to_sample:
            return "no_declared", None, None
        return "distance_only", pair_to_sample[pair], c

    fits = [c for c in tied.candidates
            if (i5, c) in pair_to_sample or (c, i7) in pair_to_sample]
    if not fits:
        return "no_declared", None, None
    if len(fits) > 1:
        return "still_ambiguous", None, None
    c = fits[0]
    pair = (i5, c) if i5 else (c, i7)
    return "resolved", pair_to_sample[pair], c


def _classify(
    read_id: str,
    five: EndCall,
    three: EndCall,
    orient: str,
    pair_to_sample: dict[tuple[str, str], str],
    accept_single_end: bool,
) -> ReadCall:
    i5, i7 = five.barcode, three.barcode

    if i5 and i7:
        sample = pair_to_sample.get((i5, i7))
        if sample:
            return ReadCall(read_id, "A", sample, i5, i7, orient)
        # both ends read cleanly but the pair is not one that was used:
        # that is a hopped or chimeric molecule, and the combination itself
        # is the evidence
        return ReadCall(read_id, "C", None, i5, i7, orient)

    if (i5 or i7) and accept_single_end:
        # one end only. Unique resolution is required: if several samples
        # share that index, a single end cannot pick between them.
        matches = [
            smp
            for (a, b), smp in pair_to_sample.items()
            if (i5 and a == i5 and not i7) or (i7 and b == i7 and not i5)
        ]
        if len(set(matches)) == 1:
            return ReadCall(read_id, "B", matches[0], i5, i7, orient)
        return ReadCall(read_id, "D", None, i5, i7, orient)

    if i5 or i7:
        return ReadCall(read_id, "B", None, i5, i7, orient)

    # copies disagreed at one end and nothing resolved at the other
    if five.copies_seen > 1 or three.copies_seen > 1:
        return ReadCall(read_id, "C", None, None, None, orient)
    return ReadCall(read_id, "D", None, None, None, orient)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def write_index_matrix(grid, samples, path: str) -> dict:
    """The full observed i7 x i5 grid, assigned combinations and not.

    Every read carries a pair of indices, and only some pairs belong to a
    sample. Reads landing on a combination no sample claims cannot be
    anything but hopping or misreading - there is no library with that
    pair - so those cells measure the run's cross-talk directly rather
    than by inference. The assigned cells are here too, because the rate
    only means something against the depth it is a fraction of.

    `grid` is a {(i5, i7): reads} tally accumulated as reads stream past,
    not a list of calls: the grid is bounded by the index sets, a list of
    calls is bounded only by the size of the run.
    """
    assigned = {(s["index_i5"], s["index_i7"]) for s in samples}
    i5s = sorted({s["index_i5"] for s in samples})
    i7s = sorted({s["index_i7"] for s in samples})

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["index_i5", "index_i7", "reads", "assigned"])
        for a in i5s:
            for b in i7s:
                w.writerow([a, b, grid.get((a, b), 0),
                            "yes" if (a, b) in assigned else "no"])

    sentinel = [(a, b) for a in i5s for b in i7s if (a, b) not in assigned]
    on_sentinels = sum(grid.get(k, 0) for k in sentinel)
    on_assigned = sum(v for k, v in grid.items() if k in assigned)
    return {
        "n_i5": len(i5s), "n_i7": len(i7s),
        "combinations_possible": len(i5s) * len(i7s),
        "combinations_assigned": len(assigned),
        "sentinel_combinations": len(sentinel),
        "reads_on_sentinels": on_sentinels,
        "reads_on_assigned": on_assigned,
        # The hopping rate this run can actually see. Hopping onto an
        # assigned pair is invisible, so the true rate is higher by the
        # reciprocal of the observable fraction.
        "sentinel_rate": (on_sentinels / (on_sentinels + on_assigned)
                          if (on_sentinels + on_assigned) else 0.0),
        "observable_fraction": (len(sentinel) / (len(i5s) * len(i7s))
                                if i5s and i7s else 0.0),
        "note": (
            "Reads on a combination no sample claims can only be index "
            "cross-talk. Hopping onto an assigned pair is indistinguishable "
            "from a real read, so the observed rate is a floor: divide by "
            "the observable fraction for an estimate of the total."
        ),
    }


class SampleWriter:
    """Per-sample FASTQ output that holds neither the run nor 384 handles.

    Reads arrive interleaved, one sample at a time being impossible, so
    there are only two obvious ways to write them out and both scale
    badly. Holding every read until the end makes memory grow with the
    run - about 1.4 KB per read, so a 20M-read flow cell wants ~28 GB and
    gets killed after an hour of work. Keeping one open gzip stream per
    sample instead costs a compressor state and a file descriptor per
    sample, on a plate that may have 384 of them.

    So each sample gets a small buffer which is appended to its file and
    dropped. gzip members concatenate, and every gzip reader handles a
    multi-member file transparently, so a file written in a thousand
    flushes reads back as one stream. Memory is then bounded by
    (samples x flush_records) rather than by the size of the run, and at
    most one file is open at a time.
    """

    def __init__(self, outdir: str, flush_records: int = 2000,
                 stats_only: bool = False):
        # stats_only counts reads per sample without writing any FASTQ.
        #
        # It saves DISK and avoids a second copy of the per-sample reads;
        # it does NOT save much time. Measured on 51k reads: 15.2s with
        # writing, 16.3s without - the write is about 2% of this stage,
        # and the per-read barcode scanning is the rest. It was added on
        # the assumption that writing dominated, which the measurement
        # disproved, so the claim is recorded here rather than repeated.
        #
        # The counters are unaffected, so the summary is identical to a
        # normal run apart from reads_written.
        self.stats_only = stats_only
        self.dir = outdir
        if not stats_only:
            os.makedirs(self.dir, exist_ok=True)
        self.flush_records = flush_records
        self._buf: dict[str, list[str]] = defaultdict(list)
        self._started: set[str] = set()
        self.counts: Counter = Counter()
        self.written = 0

    def path_for(self, sample_id: str) -> str:
        return os.path.join(self.dir, f"{sanitise_id(sample_id)}.fastq.gz")

    def add(self, sample_id: str, read_id: str, seq: str, qual: str) -> None:
        self.counts[sample_id] += 1
        if self.stats_only:
            return
        self._buf[sample_id].append(f"@{read_id}\n{seq}\n+\n{qual}\n")
        if len(self._buf[sample_id]) >= self.flush_records:
            self._flush(sample_id)

    def _flush(self, sample_id: str) -> None:
        recs = self._buf.get(sample_id)
        if not recs:
            return
        # 'wt' the first time so a re-run truncates rather than appends to
        # the output of the run before it; 'at' after that.
        mode = "at" if sample_id in self._started else "wt"
        with gzip.open(self.path_for(sample_id), mode) as fh:
            fh.write("".join(recs))
        self._started.add(sample_id)
        self.written += len(recs)
        recs.clear()

    def close(self) -> None:
        for sample_id in list(self._buf):
            self._flush(sample_id)


def code_fingerprint(read_files, **settings) -> str:
    """Identify what produced a part: code, settings, and inputs.

    Resuming reuses any file whose part already exists, which is only
    safe while the thing being resumed is the same thing. On this
    pipeline that assumption failed in practice: a first run finished 65
    of 74 files, the cut rule was then fixed, and the rerun reused all 65
    - producing a run that was 88% one algorithm and 12% another, with a
    summary that averaged the two and no sign anything was wrong.

    The inputs are fingerprinted as well as the code, because this stage
    consumes the SPLITTER's output. If splitting changes, every part here
    must be recomputed even though nothing in this file did.
    """
    h = hashlib.sha256()
    here = os.path.dirname(os.path.abspath(__file__))
    for name in ("demux_ont.py", "ampseqlib.py"):
        try:
            with open(os.path.join(here, name), "rb") as fh:
                h.update(fh.read())
        except OSError:
            h.update(name.encode())
    for k in sorted(settings):
        h.update(f"{k}={settings[k]}".encode())
    for f in sorted(read_files):
        try:
            st = os.stat(f)
            h.update(f"{os.path.basename(f)}:{st.st_size}".encode())
        except OSError:
            h.update(str(f).encode())
    return h.hexdigest()[:16]


def _build_sets(samplesheet: str, run_id: str, max_mismatch: int,
                margin: int = 2, margin_max_dist: int = 2):
    """Rebuild the index sets from files. Workers get paths, not objects:
    what the workers agree on is then the input, not a pickled snapshot of
    one process's interpretation of it."""
    samples = [s for s in read_table(samplesheet) if s["run_id"] == run_id]
    i5_set = BarcodeSet.build({s["index_i5"].upper() for s in samples},
                              max_mismatch, margin, margin_max_dist)
    i7_set = BarcodeSet.build({s["index_i7"].upper() for s in samples},
                              max_mismatch, margin, margin_max_dist)
    pair_to_sample = {
        (s["index_i5"].upper(), s["index_i7"].upper()): s["sample_id"]
        for s in samples
    }
    return samples, i5_set, i7_set, pair_to_sample


def demux_one_file(job: dict) -> dict:
    """Demultiplex one input file into one part directory.

    Every output this produces is per-part: its own per-sample FASTQ, its
    own slice of the per-read table, its own tallies. Nothing is shared
    between workers and nothing is appended to in common, so the parts can
    be produced in any order, on any number of cores, and merged after.

    The merge is cheap because of how the pieces are shaped. gzip members
    concatenate, so per-sample reads are joined by copying bytes rather
    than by decompressing and recompressing tens of gigabytes; the tables
    are appended; the tallies are summed.
    """
    cfg = load_assay(job["assay"])
    _samples, i5_set, i7_set, pair_to_sample = _build_sets(
        job["samplesheet"], job["run_id"], job["max_mismatch"],
        job.get("barcode_margin", 2), job.get("barcode_margin_max_dist", 2))

    junction_anchors: tuple[str, ...] = ()
    if job["chimera_check"]:
        a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
        junction_anchors = (a5, a7, revcomp(a5), revcomp(a7))

    part_reads = job["part_reads"]
    if os.path.isdir(part_reads):
        for stale in os.listdir(part_reads):
            if stale.endswith(".fastq.gz"):
                os.remove(os.path.join(part_reads, stale))
    writer = SampleWriter(part_reads, flush_records=job["flush_records"],
                          stats_only=job.get("stats_only", False))

    tiers: Counter = Counter()
    orients: Counter = Counter()
    grid: Counter = Counter()
    i5_why: Counter = Counter()
    i7_why: Counter = Counter()
    combo: Counter = Counter()
    combo_dist: Counter = Counter()
    # Which samples the recoverable reads would land in. A recovery that
    # all went to one well is a different finding from one spread over
    # the plate - the first looks like one bad index, the second like a
    # rule that is simply too strict.
    combo_by_sample: Counter = Counter()
    total = assigned = 0

    tmp_tsv = job["assign_tsv"] + ".partial"
    with open(tmp_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for read_id, seq, qual in open_reads(job["src"]):
            call = assign_read(
                read_id, seq, cfg, i5_set, i7_set, pair_to_sample,
                job["barcode_repeats"], job["search_window"],
                job["accept_single_end"],
                junction_anchors=junction_anchors,
                anchor_max_edits=job["anchor_max_edits"],
            )
            total += 1
            tiers[call.tier] += 1
            orients[call.orientation] += 1
            i5_why[call.i5_reason] += 1
            i7_why[call.i7_reason] += 1
            combo[call.combo] += 1
            if call.combo == "distance_only" and call.combo_dist >= 0:
                combo_dist[call.combo_dist] += 1
            if call.combo == "resolved" and call.combo_sample:
                combo_by_sample[call.combo_sample] += 1
            if call.i5 and call.i7:
                grid[(call.i5, call.i7)] += 1
            w.writerow([call.read_id, call.tier, call.sample_id or "",
                        call.i5 or "", call.i7 or "", call.orientation])
            if call.sample_id:
                assigned += 1
                if writer.stats_only:
                    # Counted through the writer, not beside it: two
                    # places incrementing the same tally is how they
                    # come to disagree.
                    writer.add(call.sample_id, read_id, "", "")
                    continue
                s = seq if call.orientation == "+" else revcomp(seq)
                q = qual if call.orientation == "+" else qual[::-1]
                if 0 <= call.insert_start < call.insert_end <= len(s):
                    s = s[call.insert_start:call.insert_end]
                    q = q[call.insert_start:call.insert_end]
                writer.add(call.sample_id, read_id, s, q)
    writer.close()
    os.replace(tmp_tsv, job["assign_tsv"])

    stats = {
        "source": job["src"],
        "total_reads": total,
        "assigned_reads": assigned,
        "by_tier": dict(tiers),
        "orientation": dict(orients),
        "i5_end_reasons": dict(i5_why),
        "i7_end_reasons": dict(i7_why),
        "combination_check": dict(combo),
        "combination_distance_only_by_distance": {
            str(d): combo_dist[d] for d in sorted(combo_dist)},
        "combination_recoverable_by_sample": dict(combo_by_sample),
        # Tuple keys are not JSON; the pair is joined with a tab, which
        # cannot occur in a barcode.
        "grid": {f"{a}\t{b}": c for (a, b), c in grid.items()},
        "per_sample": dict(writer.counts),
        "reads_written": writer.written,
        "code": job["code"],
    }
    with open(job["stats"], "w", encoding="utf-8") as fh:
        json.dump(stats, fh)
    return stats


def merge_sample_reads(part_dirs: list[str], dest_dir: str) -> dict[str, int]:
    """Join each sample's per-part FASTQ into one file, without recompressing.

    A gzip file may be several members back to back, and every reader
    handles that transparently, so concatenating the compressed bytes is
    a valid gzip file of the concatenated content. At this scale that
    matters: the alternative is decompressing and recompressing the whole
    run to achieve nothing.
    """
    os.makedirs(dest_dir, exist_ok=True)
    for stale in os.listdir(dest_dir):
        if stale.endswith(".fastq.gz"):
            os.remove(os.path.join(dest_dir, stale))
    sizes: dict[str, int] = {}
    for part in part_dirs:
        if not os.path.isdir(part):
            continue
        for name in sorted(os.listdir(part)):
            if not name.endswith(".fastq.gz"):
                continue
            src = os.path.join(part, name)
            with open(src, "rb") as fin, \
                    open(os.path.join(dest_dir, name), "ab") as fout:
                while True:
                    chunk = fin.read(1 << 20)
                    if not chunk:
                        break
                    fout.write(chunk)
            sizes[name] = sizes.get(name, 0) + os.path.getsize(src)
    return sizes


def _progress(done: int, total: int, results: list, started: float) -> None:
    reads = sum(r["total_reads"] for r in results)
    dt = max(1e-6, time.time() - started)
    left = (total - done) * (dt / max(1, done))
    sys.stderr.write(
        f"\r  {done}/{total} files  {reads:,} reads  "
        f"{reads / dt:,.0f} reads/s  ~{left / 60:.0f} min left    ")
    sys.stderr.flush()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Demultiplex nanopore reads carrying in-line Illumina indices."
    )
    ap.add_argument(
        "--reads", required=True, nargs="+",
        help=(
            "pooled reads: a FASTQ, a glob, or a MinKNOW run directory. A "
            "directory is walked, fastq_pass is used and fastq_fail is "
            "excluded and reported"
        ),
    )
    ap.add_argument(
        "--include-failed", action="store_true",
        help="also read fastq_fail - those reads failed basecalling QC and "
             "carry a different error profile",
    )
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--barcode-repeats", type=int, default=2)
    ap.add_argument(
        "--max-mismatch",
        type=int,
        default=1,
        help="requested substitutions per barcode copy; capped by what the "
        "index set can actually support",
    )
    ap.add_argument(
        "--search-window",
        type=int,
        default=120,
        help="how far from each read end to look for the anchor",
    )
    ap.add_argument("--accept-single-end", action="store_true",
                    help="fold tier B into the output; measure the cost first")
    ap.add_argument(
        "--no-chimera-check", dest="chimera_check", action="store_false",
        help="stop quarantining reads that carry an adapter between their "
             "two ends. Those reads span more than one molecule, so the "
             "index at one end need not belong to the molecule at the "
             "other - they demultiplex into a real sample that never "
             "contained them. Disabling this trades correctness for yield",
    )
    ap.add_argument("--anchor-max-edits", type=int, default=2,
                    help="edit budget for the anchor search, including the "
                         "interior scan for junctions")
    ap.add_argument(
        "--flush-records", type=int, default=2000,
        help="reads buffered per sample before appending to its file. "
             "Peak memory is roughly samples x this x read length, so it "
             "trades memory against the number of gzip members written. "
             "Lower it only to test the flush boundary",
    )
    ap.add_argument("--threads", type=int, default=0,
                    help="worker processes; 0 picks min(8, cores). Work\n                         is divided by input file and resumes per file")
    ap.add_argument("--barcode-margin", type=int, default=2,
                    help="accept a barcode beyond max_mismatch when the "
                         "runner-up is at least this much further away. "
                         "The set's spacing caps the tolerance; the "
                         "margin is a property of the individual read, so "
                         "it adds no misassignment the geometry would "
                         "have prevented. 0 disables it")
    ap.add_argument("--barcode-margin-max-dist", type=int, default=2,
                    help="the furthest a barcode may be and still be "
                         "accepted under the margin rule")
    ap.add_argument("--stats-only", action="store_true",
                    help="compute the tiers, the per-end reasons and the "
                         "index matrix WITHOUT writing per-sample FASTQ. "
                         "For re-analysing a run already demultiplexed: "
                         "saves ~8 GB and leaves the existing reads "
                         "untouched. Saves little time - the write is "
                         "about 2%% of this stage, the scanning is the "
                         "rest")
    ap.add_argument("--outdir", default="demux")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    cfg = load_assay(args.assay)

    samples = [
        s for s in read_table(args.samplesheet) if s["run_id"] == args.run_id
    ]
    if not samples:
        print(f"no samples for run '{args.run_id}' in the samplesheet", file=sys.stderr)
        return 1

    i5_set = BarcodeSet.build({s["index_i5"].upper() for s in samples}, args.max_mismatch)
    i7_set = BarcodeSet.build({s["index_i7"].upper() for s in samples}, args.max_mismatch)
    pair_to_sample = {
        (s["index_i5"].upper(), s["index_i7"].upper()): s["sample_id"] for s in samples
    }

    print(
        f"{len(samples)} samples | i5 set n={len(i5_set.by_seq)} "
        f"tolerance {i5_set.max_mismatch} | i7 set n={len(i7_set.by_seq)} "
        f"tolerance {i7_set.max_mismatch} | repeats {args.barcode_repeats}"
    )
    if i5_set.max_mismatch < args.max_mismatch or i7_set.max_mismatch < args.max_mismatch:
        print(
            f"  tolerance capped below the requested {args.max_mismatch}: the "
            "index sets in this run are too close to support it"
        )

    # A MinION run is a tree of batched FASTQ, not a file. Resolving it
    # here, and printing what was resolved, keeps the choice of which
    # reads entered the analysis visible rather than buried in a glob.
    read_files, file_report = resolve_read_files(
        args.reads, prefer_pass=not args.include_failed
    )
    if not read_files:
        raise SystemExit(
            f"no FASTQ found under {', '.join(args.reads)} - expected "
            ".fastq/.fq, optionally gzipped"
        )
    print(f"  {file_report['n_files']} FASTQ file(s)"
          + (f", {file_report['from_fastq_pass']} from fastq_pass"
             if file_report["from_fastq_pass"] else ""))
    for note in file_report["notes"]:
        print(f"  note: {note}")

    # Nothing below accumulates per-read state. Reads stream past once;
    # the per-read table is written as they go, the per-sample FASTQ are
    # buffered in small chunks, and everything else is a tally whose size
    # is set by the number of samples and indices, not the size of the
    # run. Peak memory is therefore flat in the read count - which is the
    # only property that makes a full flow cell safe to run.
    # Both anchors and both their reverse complements: a head-to-tail join
    # leaves the junction in the same orientation, a head-to-head join
    # leaves it inverted, and both fabricate a pair.
    # One shard per input file, exactly as the splitter divides its work.
    # The parts are independent, so they can run on any number of cores,
    # and a part already finished is skipped - which is what makes a stage
    # this long survive a wall-clock limit.
    parts_dir = os.path.join(args.outdir, "parts")
    os.makedirs(parts_dir, exist_ok=True)

    fingerprint = code_fingerprint(
        read_files,
        barcode_repeats=args.barcode_repeats, max_mismatch=args.max_mismatch,
        search_window=args.search_window,
        accept_single_end=args.accept_single_end,
        chimera_check=args.chimera_check,
        anchor_max_edits=args.anchor_max_edits,
        run_id=args.run_id, samplesheet=args.samplesheet,
        # In the fingerprint, and it has to be. A stats-only part records
        # itself as finished having written no FASTQ. Without this, a
        # later real run resumes onto those parts, believes the work is
        # done, finds nothing to merge, and produces an empty output
        # while reporting success - the precise failure this project has
        # already had once, when two algorithm versions shared a cache.
        stats_only=args.stats_only,
        # Changes which reads are assigned, so parts computed with a
        # different margin must not be mixed with these.
        barcode_margin=args.barcode_margin,
        barcode_margin_max_dist=args.barcode_margin_max_dist,
        assay=args.assay or "")

    jobs, results, part_dirs, tsvs, stale = [], [], [], [], 0
    for i, src in enumerate(read_files):
        stem = f"{i:05d}"
        stats_path = os.path.join(parts_dir, stem + ".json")
        tsv_path = os.path.join(parts_dir, stem + ".assignments.tsv")
        part_reads = os.path.join(parts_dir, stem + ".reads")
        part_dirs.append(part_reads)
        tsvs.append(tsv_path)
        if os.path.exists(stats_path) and os.path.exists(tsv_path):
            with open(stats_path, encoding="utf-8") as fh:
                prev = json.load(fh)
            if prev.get("code") == fingerprint:
                results.append(prev)
                continue
            stale += 1
        jobs.append({
            "src": src, "stats": stats_path, "assign_tsv": tsv_path,
            "part_reads": part_reads, "assay": args.assay,
            "samplesheet": args.samplesheet, "run_id": args.run_id,
            "max_mismatch": args.max_mismatch,
            "barcode_repeats": args.barcode_repeats,
            "search_window": args.search_window,
            "accept_single_end": args.accept_single_end,
            "chimera_check": args.chimera_check,
            "anchor_max_edits": args.anchor_max_edits,
            "flush_records": max(1, args.flush_records),
            "barcode_margin": args.barcode_margin,
            "barcode_margin_max_dist": args.barcode_margin_max_dist,
            "stats_only": args.stats_only,
            "code": fingerprint,
        })

    threads = args.threads or min(8, os.cpu_count() or 1)
    print(f"  {len(read_files)} file(s), {len(results)} already done, "
          f"{len(jobs)} to do on {threads} thread(s)")
    if stale:
        print(f"  {stale} previously-finished file(s) were made by different "
              "code, settings or inputs\n  and will be redone - a resumed run "
              "must not blend two versions")

    started = time.time()
    if jobs:
        if threads > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                for k, stats in enumerate(pool.map(demux_one_file, jobs), 1):
                    results.append(stats)
                    _progress(k, len(jobs), results, started)
        else:
            for k, job in enumerate(jobs, 1):
                results.append(demux_one_file(job))
                _progress(k, len(jobs), results, started)
        sys.stderr.write("\n")

    # -- merge --------------------------------------------------------------
    # Parts are merged in input order, so the output is the same whatever
    # order the workers finished in and whatever thread count was used.
    reads_dir = os.path.join(args.outdir, "reads")
    # Nothing to merge in stats-only mode: no part wrote a FASTQ.
    if not args.stats_only:
        merge_sample_reads(part_dirs, reads_dir)

    assign_path = os.path.join(args.outdir, "assignments.tsv")
    with open(assign_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["read_id", "tier", "sample_id", "index_i5", "index_i7",
                    "orientation"])
        for tsv in tsvs:
            if not os.path.exists(tsv):
                continue
            with open(tsv, encoding="utf-8") as part:
                for line in part:
                    fh.write(line)

    by_source = {r["source"]: r for r in results}
    ordered = [by_source[s] for s in read_files if s in by_source]
    tiers: Counter = Counter()
    orients: Counter = Counter()
    grid: Counter = Counter()
    per_sample: Counter = Counter()
    i5_why: Counter = Counter()
    i7_why: Counter = Counter()
    combo: Counter = Counter()
    combo_dist: Counter = Counter()
    combo_by_sample: Counter = Counter()
    total = assigned = written = 0
    for r in ordered:
        total += r["total_reads"]
        assigned += r["assigned_reads"]
        written += r["reads_written"]
        tiers.update(r["by_tier"])
        orients.update(r["orientation"])
        per_sample.update(r["per_sample"])
        # Absent from parts written by an older version of this tool, so
        # a resumed run that blends vintages reports what it has rather
        # than dying on a missing key.
        i5_why.update(r.get("i5_end_reasons", {}))
        i7_why.update(r.get("i7_end_reasons", {}))
        combo.update(r.get("combination_check", {}))
        combo_dist.update({int(k): v for k, v in r.get(
            "combination_distance_only_by_distance", {}).items()})
        combo_by_sample.update(
            r.get("combination_recoverable_by_sample", {}))
        for key, c in r["grid"].items():
            a, b = key.split("\t")
            grid[(a, b)] += c

    # -- summary ------------------------------------------------------------
    # Every read lands in exactly one of the known tiers, and the tiers sum
    # to the input. Asserted rather than reported: if reads have vanished,
    # the results may be wrong for reasons unrelated to biology. Summing
    # only over TIERS rather than over everything seen also catches a tier
    # that is not one of the four, which a total-vs-total check cannot.
    accounted = sum(tiers.get(t, 0) for t in TIERS)

    summary = {
        "run_id": args.run_id,
        "total_reads": total,
        "assigned_reads": assigned,
        "by_tier": {t: tiers.get(t, 0) for t in TIERS},
        "tier_meaning": TIER_MEANING,
        # Per END, not per read: a read has two ends and they fail for
        # different reasons. Counts sum to total_reads in each map.
        #
        # This is the breakdown that says whether a disappointing run is
        # chemistry or software. anchor_not_found dominating means the
        # molecule never carried that end - PCR2 did not add the index,
        # and no amount of tolerance will recover it. no_copy_matched or
        # copies_disagreed dominating means the index IS there and the
        # matcher is rejecting it, which is a parameter away from being
        # fixed. Those two have the same symptom in the tiers and
        # opposite remedies.
        "i5_end_reasons": dict(i5_why),
        "i7_end_reasons": dict(i7_why),
        # Per READ, and a measurement only - nothing here has changed a
        # tier, a sample or a file. `resolved` is the size of the prize:
        # reads whose tied end is narrowed to exactly one declared well
        # by the end that did resolve. See combination_check().
        "combination_check": dict(combo),
        "combination_meaning": COMBO_MEANING,
        # How saturated the plate is, because that is what decides
        # whether the combination can resolve anything at all. A full
        # grid leaves no room: every tie has both candidates forming
        # real wells. Reported beside the verdicts so the two are never
        # read apart.
        "combination_grid": {
            "i5_used": len({a for a, _ in pair_to_sample}),
            "i7_used": len({b for _, b in pair_to_sample}),
            "combinations_possible": (
                len({a for a, _ in pair_to_sample})
                * len({b for _, b in pair_to_sample})),
            "wells_declared": len(pair_to_sample),
            "note": (
                "The combination check can only resolve a tie when some "
                "of the tied candidates do NOT form declared wells. The "
                "closer wells_declared is to combinations_possible, the "
                "less it can ever do."
            ),
        },
        "combination_distance_only_by_distance": {
            str(d): combo_dist[d] for d in sorted(combo_dist)},
        "combination_recoverable_by_sample": dict(
            sorted(combo_by_sample.items(), key=lambda kv: -kv[1])[:20]),
        "end_reason_meaning": END_REASON_MEANING,
        "orientation": dict(orients),
        "index_matrix": write_index_matrix(
            grid, samples, os.path.join(args.outdir, "index_matrix.tsv")),
        "input_files": file_report,
        # Empty when the reads have already been through the splitter,
        # which rewrites headers - split_summary.json carries it then.
        "basecalling": read_run_metadata(read_files),
        "settings": {
            "barcode_repeats": args.barcode_repeats,
            "max_mismatch_requested": args.max_mismatch,
            "max_mismatch_i5": i5_set.max_mismatch,
            "max_mismatch_i7": i7_set.max_mismatch,
            "barcode_margin": args.barcode_margin,
            "barcode_margin_max_dist": args.barcode_margin_max_dist,
            "allow_indels": False,
            "require_copy_agreement": args.barcode_repeats >= 2,
            "accept_single_end": args.accept_single_end,
            "chimera_check": args.chimera_check,
            "anchor_max_edits": args.anchor_max_edits,
        },
        # How many parts this invocation actually computed. The driver
        # uses it to decide whether per-sample assignment downstream is
        # still valid, since its inputs have just changed.
        "parts_recomputed": len(jobs),
        "read_accounting_balanced": accounted == total,
        # A second, independent check: the reads written to per-sample
        # files must equal the reads called to a sample. The tier tally and
        # the FASTQ output are produced by different code paths, so this
        # catches a dropped buffer that the tier check cannot see.
        "reads_written": written,
        # Not asserted in stats-only mode, because nothing was written.
        # Reported as None rather than True: a check that did not run is
        # not a check that passed.
        "reads_written_balanced": None if args.stats_only else written == assigned,
        "stats_only": args.stats_only,
        "per_sample": dict(sorted(per_sample.items())),
    }
    with open(os.path.join(args.outdir, "demux_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n{total:,} reads")
    for t in TIERS:
        n = tiers.get(t, 0)
        print(f"  tier {t}  {n:>8,}  {100*n/max(1,total):5.2f}%  {TIER_MEANING[t]}")
    print(f"  assigned {assigned:,} ({100*assigned/max(1,total):.2f}%) "
          f"to {len(per_sample)} samples")

    # Printed, not merely recorded. The tiers say how many reads failed;
    # this says what failed, and it is the difference between a bench
    # problem and a parameter problem. A number that only exists in a
    # JSON file nobody opens is half a diagnostic.
    for end, why in (("i5", i5_why), ("i7", i7_why)):
        if not why:
            continue
        print(f"\n  {end} ends")
        for reason, n in sorted(why.items(), key=lambda kv: -kv[1]):
            print(f"    {n:>9,}  {100*n/max(1,total):5.2f}%  {reason}"
                  f"{'' if reason.startswith('ok') else ' - ' + END_REASON_MEANING.get(reason, '')}")
    lost = i5_why.get("anchor_not_found", 0) + i7_why.get("anchor_not_found", 0)
    matched_but_rejected = sum(
        w.get(k, 0) for w in (i5_why, i7_why)
        for k in ("no_copy_matched", "copies_disagreed"))
    if lost or matched_but_rejected:
        print(f"\n  ends with no anchor at all      {lost:>9,}  "
              "(index never added - PCR2 efficiency)")
        print(f"  ends read but not matched       {matched_but_rejected:>9,}  "
              "(index present - matching or index set)")

    # The combination check. Printed because it is the one number that
    # says whether the current per-end rule is leaving anything on the
    # table, and a number nobody sees does not get acted on.
    if combo:
        decided = sum(v for k, v in combo.items() if k != "not_applicable")
        if decided:
            print(f"\n  single-ended reads whose tie the other end could "
                  f"narrow ({decided:,} with a tie to narrow)")
            for k in ("resolved", "distance_only", "still_ambiguous",
                      "no_declared", "unreadable"):
                n = combo.get(k, 0)
                if n:
                    print(f"    {k:16s} {n:>9,}  {100*n/decided:5.1f}%  "
                          f"{COMBO_MEANING[k]}")
            r = combo.get("resolved", 0)
            do = combo.get("distance_only", 0)
            print(f"\n    SAFE  {r:,} reads ({100*r/max(1,total):.2f}% of "
                  f"input) recoverable from the PAIR - genuinely tied, "
                  f"narrowed to one declared well.")
            if do:
                print(f"    NOT   {do:,} reads ({100*do/max(1,total):.2f}% "
                      f"of input) had a single candidate too far out. "
                      f"Taking those is a tolerance increase, which this "
                      f"index set does not support - not a combination win.")
                if combo_dist:
                    lo = ", ".join(f"d={d}: {combo_dist[d]:,}"
                                   for d in sorted(combo_dist))
                    print(f"            by distance - {lo}")
            print("    NOTHING WAS CHANGED - this is a measurement.")

    if accounted != total:
        print(
            f"\nread accounting does not balance: {accounted:,} of {total:,}",
            file=sys.stderr,
        )
        return 1
    # Not applicable in stats-only mode: nothing was written, by
    # request. Left unguarded this printed a hard error and returned 1
    # on a run that had worked perfectly, which any scheduler reports as
    # a failed job.
    if written != assigned and not args.stats_only:
        print(
            f"\n{written:,} reads were written to per-sample files but "
            f"{assigned:,} were called to a sample - the FASTQ output is "
            "incomplete and must not be analysed",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
