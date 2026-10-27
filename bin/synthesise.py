#!/usr/bin/env python3
"""
Run-level synthesis for nf-helminth-amplicon.

Counts become calls here, and this is the step where a pipeline most easily
starts asserting more than it measured. Three decisions are made, and each
one is reported with the evidence behind it rather than as a bare verdict.

QUANTIFICATION
    Two spike-ins, two orders of magnitude apart, go into every sample after
    extraction at a fixed dose per unit of faecal DNA. A read count is
    compositional - it says what share of the library a reference took, not
    how much was there - so the spike-ins are what convert a share back into
    an amount. The unit is copies per unit DNA, not per gram of stool,
    because a fixed mass of DNA is what was actually spiked.

    Two points also make the calibration checkable. If the assay is linear
    the fitted slope through them is 1. It is reported every time, because
    a slope well below 1 means the high standard is saturating and every
    copy number on that plate is compressed toward the middle.

DETECTION
    The limit is measured, not assumed. Negative controls on the same plate
    carry the same index hopping, the same carryover and the same reagent
    contamination as everything else, so what appears in them is what a
    well can show while containing nothing. Two models are computed - an
    absolute floor, because hopped reads arrive at a rate set by the run
    rather than by the receiving well, and a proportional one, because
    carryover scales with depth - and the stricter applies.

INTERPRETATION
    A negative is only meaningful if the sample could have produced a
    positive. The spike-ins establish that: if they are absent, the well
    failed and `assay_failed` is the honest answer, not `not_detected`.
    Where no spike was used at all the negative is reported as unverified
    rather than quietly promoted to a true negative.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table, restore_declared_ids  # noqa: E402

CALLS = ("detected", "not_detected", "assay_failed")

CALL_MEANING = {
    "detected": "above the plate's measured detection limit",
    "not_detected": "below it, in a well whose spike-ins worked",
    "assay_failed": "the well cannot support either answer",
}

INTERNAL = ("spiked", "unspiked", "spike_failed", "slope_implausible")

INTERNAL_MEANING = {
    "spiked": "both standards recovered, within the expected ratio",
    "unspiked": "no standard present - a negative here is unverified",
    "spike_failed": "standards present but wrong, so the well is unreliable",
    "slope_implausible": (
        "both standards recovered, so the well amplified and its calls "
        "stand, but they did not scale with their nominal concentrations "
        "closely enough to convert reads to copies - no copy number is "
        "given"
    ),
}

# The band of standard-curve slopes from which a copy number is still
# worth emitting. Separate from --slope-tolerance, which only decides when
# a slope earns a note: a slope of 0.8 is worth a warning but still gives
# copy numbers in the right order of magnitude, whereas outside this band
# the inversion stops being a measurement at all. The arithmetic is why:
# copies are 10 ** (log-distance / slope), so a slope of 0.05 multiplies
# every log-distance from the standards by twenty, and a negative slope
# makes copies FALL as reads rise. A slope of 0.5 already means a hundred-
# fold input difference produced a tenfold read difference; 2 means the
# reverse. Beyond either, the two standards are not describing a dose-
# response line and nothing read off it should be reported as an amount.
#
# Detection is not affected. Both standards were recovered, which is
# exactly the evidence that the well could have produced a positive, so
# its detected and not_detected calls stand; only quantification is
# withheld.
SLOPE_PLAUSIBLE = (0.5, 2.0)

NEGATIVE_TYPES = {"neg_ntc", "neg_extraction", "blank", "negative"}

# Wells whose content IS the tagged control plasmids. A mock is assembled
# from them, so in a mock the plasmids are called and quantified like
# targets; in every other well they are contamination.
PLASMID_WELLS = {"mock"}

# Amplification efficiency relative to the standards, when the panel does
# not state one. Unity is the right default - it is what "no correction"
# already meant arithmetically - but the output records that it was assumed
# rather than measured.
DEFAULT_CORR_FACTOR = 1.0

# Only speciation references form a composition; a variants-only locus is
# counted but never enters the denominator.
SPECIATION = "speciation"


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


@dataclass
class Sample:
    sample_id: str
    sample_type: str = "unknown"
    # Who the library came from. Carried so that replicates can be
    # grouped without a metadata file - several libraries of one subject
    # is the commonest grouping there is, and it is already declared in
    # the samplesheet.
    subject_id: str = ""
    index_i5: str = ""
    index_i7: str = ""
    plate_row: str = ""
    plate_col: int = 0
    counts: dict = field(default_factory=dict)
    total_target: int = 0
    depth: int = 0
    outcomes: dict = field(default_factory=dict)
    by_ref_type: dict = field(default_factory=dict)
    # Read-weighted edit-distance histograms per outcome, as assignment
    # measured them. Carried up because a per-sample JSON that nobody
    # opens is not an output: on a 94-well plate the choice was 94 files
    # by hand or nothing, and it was nothing.
    dist: dict = field(default_factory=dict)


def load_counts(assign_dir: str) -> dict[str, Sample]:
    out: dict[str, Sample] = {}
    for path in sorted(glob.glob(os.path.join(assign_dir, "*.counts.tsv"))):
        for r in read_table(path):
            sid = r["sample_id"]
            s = out.setdefault(sid, Sample(sample_id=sid))
            s.counts[r["seq_id"]] = int(r["reads"])
    for path in sorted(glob.glob(os.path.join(assign_dir, "*.assign.json"))):
        with open(path, encoding="utf-8") as fh:
            j = json.load(fh)
        s = out.get(j["sample_id"])
        if s:
            s.depth = j.get("sequences_in", 0)
            s.outcomes = j.get("outcomes", {}) or {}
            s.by_ref_type = j.get("reads_by_reference_type", {}) or {}
            s.dist = j.get("distance_by_outcome", {}) or {}
    return out


# The fates a read can meet once it reaches assignment. Together they
# partition the depth exactly: the first four come from what the read was
# assigned TO, the rest from why it was not assigned at all.
FATES = ("target", "quant", "control", "vector",
         "ambiguous", "tag_unobservable", "unmatched", "no_primer")

FATE_MEANING = {
    "target": "a panel target - the only reads a composition is built from",
    "quant": "a spike-in standard, which is how depth becomes an amount",
    "control": "a tagged control plasmid, so plasmid carryover",
    "vector": "plasmid backbone, so non-specific amplification",
    "ambiguous": "fits references from more than one taxon equally well",
    "tag_unobservable": "target and control tied, the tag between the mates",
    "unmatched": "nothing in the panel within the distance limit",
    "no_primer": "no primer at the read start",
}


# What prep does with a read pair before assignment ever sees it. The
# first two survive; the rest are losses, and each names a different
# problem at the bench.
PREP_KEPT = ("merged", "paired")
PREP_LOST = ("no_primer", "too_short", "low_quality", "length_mismatch")

PREP_MEANING = {
    "no_primer": "no primer at either read start - dimer, off-target "
                 "product, or a read too truncated to place",
    "too_short": "insert shorter than the minimum once trimmed",
    "low_quality": "expected errors above the threshold",
    "length_mismatch": "mates overlapped inconsistently",
}


def load_prep(prep_dir: str | None) -> dict[str, dict]:
    """Per-sample prep outcomes, if the directory was given."""
    out: dict[str, dict] = {}
    if not prep_dir:
        return out
    for path in sorted(glob.glob(os.path.join(prep_dir, "*.prep.json"))):
        with open(path, encoding="utf-8") as fh:
            j = json.load(fh)
        out[j["sample_id"]] = j
    return out


def read_funnel(
    samples: dict, prep: dict, demux: dict | None, fate_totals: Counter,
    split: dict | None = None,
) -> dict:
    """Every read from the sequencer to a species call, and where they go.

    The per-sample figures start at reads that reached assignment, which
    is already net of everything prep discarded. That hides a specific
    failure: a well that produced plenty of data and lost most of it to
    primer dimer looks merely shallow, and shallow calls for reloading the
    pool while dimer calls for fixing the reaction. The funnel is where
    the difference shows.

    It also checks across a tool boundary that nothing vanished silently -
    reads leaving prep must equal reads entering assignment. Nothing else
    in the pipeline compares those two numbers, so a file lost between the
    steps would otherwise go unnoticed.
    """
    stages, notes = [], []

    raw = sum(p.get("pairs_in", 0) for p in prep.values())
    kept = sum(p.get("fates", {}).get(k, 0)
               for p in prep.values() for k in PREP_KEPT)
    lost_by_reason = {
        k: sum(p.get("fates", {}).get(k, 0) for p in prep.values())
        for k in PREP_LOST
    }
    lost_by_reason = {k: v for k, v in lost_by_reason.items() if v}

    if demux:
        total = demux.get("total_reads", 0)
        tiers = demux.get("by_tier") or {}
        meaning = demux.get("tier_meaning") or {}
        assigned = demux.get("assigned_reads", tiers.get("A", 0))
        if total:
            # On the nanopore path a read may have been cut into several
            # amplicon units before demultiplexing, so the count entering
            # demux is not the count that came off the sequencer. Start
            # the funnel at the sequencer, or the split silently changes
            # the denominator every later percentage is taken against.
            if split and split.get("reads_in"):
                stages.append({"stage": "reads in the run",
                               "reads": split["reads_in"], "lost_here": {}})
                dropped = split.get("pieces_dropped_below_min_length", 0)
                stages.append({
                    "stage": "amplicon units after splitting",
                    "reads": split.get("pieces_out", total),
                    "lost_here": ({"pieces too short to carry an insert":
                                   dropped} if dropped else {}),
                })
            else:
                stages.append({"stage": "reads in the run", "reads": total,
                               "lost_here": {}})
            stages.append({
                "stage": "demultiplexed to a sample", "reads": assigned,
                # Name the tiers rather than lumping them: B is a read
                # whose second barcode was unreadable, C is a chimera, D
                # is neither - and they point at different problems.
                "lost_here": {
                    f"tier {k}"
                    + (f" ({meaning[k].split(' - ')[0]})" if k in meaning
                       else ""): v
                    for k, v in sorted(tiers.items()) if k != "A" and v
                },
            })
    elif raw:
        stages.append({"stage": "read pairs from the sequencer",
                       "reads": raw, "lost_here": {}})
        stages.append({"stage": "passed read preparation", "reads": kept,
                       "lost_here": lost_by_reason})

    entering = sum(s.depth for s in samples.values())
    stages.append({"stage": "reached assignment", "reads": entering,
                   "lost_here": {}})

    mapped = sum(fate_totals.get(k, 0)
                 for k in ("target", "quant", "control", "vector"))
    stages.append({
        "stage": "mapped to the panel", "reads": mapped,
        "lost_here": {k: fate_totals.get(k, 0)
                      for k in ("ambiguous", "tag_unobservable",
                                "unmatched", "no_primer")
                      if fate_totals.get(k, 0)},
    })
    stages.append({
        "stage": "mapped to a target", "reads": fate_totals.get("target", 0),
        "lost_here": {k: fate_totals.get(k, 0)
                      for k in ("quant", "control", "vector")
                      if fate_totals.get(k, 0)},
    })

    for st in stages:
        st["percent_of_start"] = (
            round(100 * st["reads"] / stages[0]["reads"], 2)
            if stages and stages[0]["reads"] else None
        )

    # The cross-boundary check, on whichever path produced the reads.
    # Both platforms have the same seam: a step hands reads to assignment
    # and nothing else in the pipeline compares the two counts, so a file
    # dropped between them would leave every stage balancing internally
    # while data quietly went missing.
    boundary_ok, left, whence = None, None, None
    if demux and demux.get("assigned_reads") is not None:
        left, whence = demux["assigned_reads"], "demultiplexing"
    elif kept:
        left, whence = kept, "preparation"
    if left is not None:
        boundary_ok = left == entering
        if not boundary_ok:
            notes.append(
                f"{left:,} reads left {whence} but {entering:,} entered "
                "assignment. Nothing else compares these two numbers, so "
                "the difference is reads lost between the steps rather "
                "than by any reported filter"
            )

    # a run median hides the one well that failed
    outliers = []
    for sid, p in sorted(prep.items()):
        n = p.get("pairs_in", 0)
        if n < 50:
            continue
        k = sum(p.get("fates", {}).get(x, 0) for x in PREP_KEPT)
        if k / n < 0.5:
            outliers.append({"sample_id": sid, "pairs_in": n,
                             "kept_fraction": round(k / n, 4),
                             "fates": p.get("fates", {})})
    if outliers:
        notes.append(
            f"{len(outliers)} well(s) lost more than half their reads "
            "before assignment. A low bar in the depth panel reads as a "
            "shallow well, but these are wells where the reaction, not "
            "the loading, is the problem"
        )

    return {
        "stages": stages,
        "prep_loss_reasons": PREP_MEANING,
        "handoff_balances": boundary_ok,
        "handoff_from": whence,
        "wells_losing_over_half_at_prep": outliers[:8],
        "notes": notes,
    }


def read_fates(s: Sample) -> dict:
    """Where a sample's reads went, as a partition of its depth.

    A depth and a composition together still leave the important question
    unanswered: what fraction of the run was usable at all. A well can
    carry five thousand reads and fifty targets, and the difference is
    diagnostic - dimer, off-target amplification, plasmid carryover and
    reference gaps each mean something different and each call for a
    different fix.
    """
    out = {f: 0 for f in FATES}
    for k in ("target", "quant", "control", "vector"):
        out[k] = int(s.by_ref_type.get(k, 0))
    for k in ("ambiguous", "tag_unobservable", "unmatched", "no_primer"):
        out[k] = int(s.outcomes.get(k, 0))
    return out


def well_of(row: str, col: int) -> str:
    return f"{row}{col}" if row and col else ""


# --------------------------------------------------------------------------
# calibration
# --------------------------------------------------------------------------


@dataclass
class Calibration:
    sample_id: str
    slope: float = float("nan")
    intercept: float = float("nan")
    reads_hi: int = 0
    reads_lo: int = 0
    conc_hi: float = 0.0
    conc_lo: float = 0.0
    state: str = "unspiked"
    note: str = ""

    def copies(self, reads: float) -> float | None:
        """Invert the fit. Returns None when there is no usable calibration,
        which is different from returning zero."""
        if reads <= 0 or self.state != "spiked":
            # slope_implausible lands here too: its standards were found,
            # but the line through them cannot be inverted honestly.
            return None
        if not (self.slope and math.isfinite(self.slope)):
            return None
        return 10 ** ((math.log10(reads) - self.intercept) / self.slope)


def run_was_spiked(samples: dict, quants: list, min_reads: int) -> bool:
    """Do the reads LOOK like a spiked run - a quarter of wells or more
    recovering a standard?

    Advisory only. This rule once decided whether a run was spiked, and it
    misreads a run whose spikes mostly failed and a run with widespread
    carryover. Whether the standards went in is now declared with
    --spiked; this only checks the declaration against the reads (see
    spike_declaration).
    """
    if not quants:
        return False
    hits = 0
    for s in samples.values():
        if any(s.counts.get(q["seq_id"], 0) >= min_reads for q in quants):
            hits += 1
    return hits >= max(1, 0.25 * len(samples))


def spike_declaration(spiked: bool, samples: dict, quants: list,
                      min_reads: int) -> dict:
    """Record whether the run was spiked, and whether the reads agree.

    Whether the standards went in is a fact of the protocol, so it is
    declared with --spiked and never inferred. Inferring it from the
    reads - spiked if a quarter of wells recovered a standard - was the
    only route once, and it fails both ways: a run whose spikes mostly
    failed reads as never spiked, so real failures pass as unverified
    negatives, and carryover across a quarter of the plate reads as a
    spike, so every clean well fails.

    The same rule survives here only as a cross-check. It decides
    nothing, but a declaration the reads disagree with is worth seeing:
    few standards on a run declared spiked is a run-wide spike failure,
    and many on a run declared unspiked is carryover - or a run that was
    spiked and is missing --spiked.
    """
    with_std = sum(
        1 for s in samples.values()
        if any(s.counts.get(q["seq_id"], 0) >= min_reads for q in quants))
    looks_spiked = run_was_spiked(samples, quants, min_reads)
    decision = {
        "spiked": spiked,
        "source": "--spiked" if spiked else "default: --spiked not given",
        "wells_with_standard": with_std,
        "wells": len(samples),
    }
    if spiked and not looks_spiked:
        decision["note"] = (
            f"declared spiked, but only {with_std} of {len(samples)} wells "
            "recovered a standard. If the spike really was added, most "
            "wells lost it - and those wells are assay_failed."
        )
    elif not spiked and looks_spiked:
        decision["note"] = (
            f"not declared spiked, but {with_std} of {len(samples)} wells "
            "carry standard reads. If the standards were added, rerun with "
            "--spiked to get copy numbers and verified negatives; if not, "
            "this is carryover."
        )
    return decision


def undeclared_standards(s: Sample, quants: list,
                         min_reads: int) -> "Calibration":
    """A well on an unspiked run (no --spiked).

    Never calibrated, whatever standard reads it holds. Fitting a curve to
    standards that were not added would turn carryover into copy numbers,
    so the reads are counted, named as contamination, and left alone.
    """
    stray = sum(s.counts.get(q["seq_id"], 0) for q in quants)
    note = "run not spiked (no --spiked) - every negative is unverified"
    if stray >= min_reads:
        note += (f"; {stray} standard reads present anyway, treated as "
                 "carryover rather than a calibration")
    return Calibration(s.sample_id, state="unspiked", note=note)


def calibrate(
    s: Sample, quants: list, min_reads: int, slope_tol: float,
    run_spiked: bool = True,
) -> Calibration:
    """Fit log10(reads) = intercept + slope * log10(copies) through the
    standards.

    With two standards the line is determined exactly, so the fit cannot
    fail - which is why the slope has to be inspected rather than trusted.
    A slope of 1 means a tenfold change in input produced a tenfold change
    in reads. Below that, the dynamic range is compressed and the high end
    is saturating; above it, the low standard is probably being lost to the
    noise floor rather than measured.
    """
    pts = []
    for q in quants:
        n = s.counts.get(q["seq_id"], 0)
        conc = float(q.get("default_conc") or 0)
        if conc > 0:
            pts.append((q["seq_id"], n, conc))
    if not pts:
        return Calibration(s.sample_id, state="unspiked",
                           note="no quant references in the panel")

    pts.sort(key=lambda p: -p[2])
    hi = pts[0]
    lo = pts[-1] if len(pts) > 1 else None
    cal = Calibration(
        s.sample_id, reads_hi=hi[1], conc_hi=hi[2],
        reads_lo=(lo[1] if lo else 0), conc_lo=(lo[2] if lo else 0.0),
    )

    present = [p for p in pts if p[1] >= min_reads]
    if not present:
        if run_spiked:
            cal.state = "spike_failed"
            cal.note = (
                "no standard recovered on a run declared spiked (--spiked) "
                "- this well failed rather than going unspiked, so nothing "
                "in it can be called negative"
            )
        else:
            cal.state = "unspiked"
            cal.note = (
                "no standard recovered anywhere on this run - the protocol "
                "appears not to have included a spike, so every negative "
                "here is unverified"
            )
        return cal
    if len(pts) > 1 and len(present) < len(pts):
        # The high standard missing is a failed well. The low one missing is
        # ordinary at shallow depth and only costs the bottom of the range.
        if hi[1] < min_reads:
            cal.state = "spike_failed"
            cal.note = (
                f"the high standard returned {hi[1]} reads - the well did "
                "not amplify, so nothing in it can be called negative"
            )
            return cal
        cal.state = "spiked"
        cal.slope = 1.0
        cal.intercept = math.log10(hi[1]) - math.log10(hi[2])
        cal.note = (
            "only the high standard was recovered; slope assumed 1 and the "
            "low end of the range is unverified"
        )
        return cal

    if len(present) == 1:
        cal.state = "spiked"
        cal.slope = 1.0
        cal.intercept = math.log10(present[0][1]) - math.log10(present[0][2])
        cal.note = "single standard; slope assumed 1"
        return cal

    x1, y1 = math.log10(lo[2]), math.log10(max(lo[1], 1))
    x2, y2 = math.log10(hi[2]), math.log10(max(hi[1], 1))
    if x2 == x1:
        cal.state = "spike_failed"
        cal.note = "standards have the same nominal concentration"
        return cal
    cal.slope = (y2 - y1) / (x2 - x1)
    cal.intercept = y1 - cal.slope * x1
    cal.state = "spiked"
    lo_ok, hi_ok = SLOPE_PLAUSIBLE
    if not (lo_ok <= cal.slope <= hi_ok):
        # The slope and intercept are kept and reported - they are the
        # evidence - but the state changes so that copies() declines to
        # invert them. Not spike_failed: the standards DID amplify, so the
        # well can support a negative, and failing its calls would discard
        # a detection result to punish a quantification problem.
        cal.state = "slope_implausible"
        cal.note = (
            f"slope {cal.slope:.2f}, outside the plausible "
            f"{lo_ok:g}-{hi_ok:g} - the standards did not describe a "
            "dose-response line, so no copy number is given for this "
            "well; detection calls still stand, since both standards "
            "were recovered"
        )
        return cal
    if abs(cal.slope - 1.0) > slope_tol:
        cal.note = (
            f"slope {cal.slope:.2f} rather than 1 - the standards did not "
            "scale with their nominal concentrations, so copy numbers from "
            "this well are compressed"
            if cal.slope < 1 else
            f"slope {cal.slope:.2f} rather than 1 - the low standard is "
            "probably at or below the noise floor here"
        )
    return cal


# --------------------------------------------------------------------------
# detection limit
# --------------------------------------------------------------------------


def detection_limits(
    samples: dict[str, Sample], seq_ids: list[str], floor: int,
    min_fraction_depth: int,
) -> dict[str, dict]:
    """What a well can show while containing nothing.

    Derived per reference from the run's own negative controls, which is
    the only place the answer is measurable. Two bounds, because the two
    things that put reads in an empty well behave differently: index
    hopping delivers a roughly fixed number of reads regardless of how deep
    the receiving well is, while carryover and cross-contamination scale
    with it.

    The proportional bound is taken only from blanks sequenced to at least
    min_fraction_depth. A blank is near-empty by design, and a fraction
    over a handful of reads is noise rather than a carryover rate: an NTC
    of 3 reads with 1 hopped read would otherwise read as 33%, and hold a
    20,000-read sample to 6,667 reads before anything could be detected.
    Shallow blanks still set the absolute bound - their stray reads are
    exactly the fixed-count events that bound exists for.
    """
    blanks = [s for s in samples.values() if s.sample_type in NEGATIVE_TYPES]
    deep = [s for s in blanks if s.depth >= min_fraction_depth and s.depth]
    out: dict[str, dict] = {}
    for sid in seq_ids:
        obs = [s.counts.get(sid, 0) for s in blanks]
        fracs = [s.counts.get(sid, 0) / s.depth for s in deep]
        out[sid] = {
            "n_blanks": len(blanks),
            "n_blanks_for_fraction": len(deep),
            "min_fraction_depth": min_fraction_depth,
            "max_reads_in_blank": max(obs) if obs else 0,
            "max_fraction_in_blank": max(fracs) if fracs else 0.0,
            "absolute_floor": floor,
        }
    return out


def limit_for(sid: str, s: Sample, lods: dict) -> tuple[int, str]:
    d = lods.get(sid)
    if not d:
        return 1, "no blanks on this run - the floor is a single read, "\
                  "which is an assumption rather than a measurement"
    # detection_limits writes an entry for every reference whether or not
    # the run had blanks, so the branch above never sees a run without
    # them. Without this one, zero blanks fell through to "blanks were
    # clean" - a measurement claimed where none was made. The fixed floor
    # still applies; it is just not evidence of anything.
    if not d.get("n_blanks"):
        return d["absolute_floor"], (
            f"the fixed floor of {d['absolute_floor']} - there were no "
            "blanks on this run, so the limit is an assumption rather "
            "than a measurement")
    absolute = d["max_reads_in_blank"]
    proportional = int(math.ceil(d["max_fraction_in_blank"] * s.depth))
    limit = max(absolute, proportional, d["absolute_floor"])
    which = (
        "absolute (hopping-like, flat across wells)"
        if absolute >= proportional and absolute >= d["absolute_floor"]
        else "proportional to depth (carryover-like)"
        if proportional >= d["absolute_floor"]
        else "the fixed floor - blanks were clean"
    )
    return limit, which


# --------------------------------------------------------------------------
# plate geometry
# --------------------------------------------------------------------------


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


def permutation_p(
    values: list[float], labels: list, rng: random.Random, n_perm: int
) -> tuple[float, float]:
    """How unusual is the grouping, if well position means nothing?

    Wells are meant to be independent, so the null is simply that the
    labels could have been attached to any well. Shuffling them is a more
    honest null than any parametric one, and it costs nothing here.
    """
    by = defaultdict(list)
    for v, lab in zip(values, labels):
        by[lab].append(v)
    groups = [g for g in by.values() if g]
    if len(groups) < 2:
        return 0.0, 1.0
    observed = kruskal_h(groups)

    shuffled = list(labels)
    worse = 0
    for _ in range(n_perm):
        rng.shuffle(shuffled)
        by2 = defaultdict(list)
        for v, lab in zip(values, shuffled):
            by2[lab].append(v)
        if kruskal_h([g for g in by2.values() if g]) >= observed:
            worse += 1
    # +1 on both sides: a permutation test can never justify p = 0, and
    # reporting one would overstate what the resampling can support.
    return observed, (worse + 1) / (n_perm + 1)


def plate_qc(
    samples: dict[str, Sample], taxa: list[str], rng: random.Random,
    n_perm: int,
) -> list[dict]:
    """Row, column and edge effects.

    Every well is an independent sample of a different person's stool.
    Position should therefore explain nothing, and if it does, the cause is
    the plate: a multichannel channel, a thermal gradient, evaporation at
    the edges, a bad column of tips. This does not say which - only that
    the signal is not coming from biology.
    """
    placed = [
        s for s in samples.values()
        if s.plate_row and s.plate_col and s.sample_type not in NEGATIVE_TYPES
    ]
    if len(placed) < 12:
        return []
    rows = sorted({s.plate_row for s in placed})
    cols = sorted({s.plate_col for s in placed})
    findings = []
    for what in ["depth"] + taxa:
        if what == "depth":
            vals = [float(s.depth) for s in placed]
        else:
            # Share of the library, not raw reads: otherwise a row effect in
            # depth would show up again as a row effect in every taxon.
            vals = [
                (s.counts.get(what, 0) / s.depth) if s.depth else 0.0
                for s in placed
            ]
        if len(set(vals)) < 3:
            continue
        row_h, row_p = permutation_p(
            vals, [s.plate_row for s in placed], rng, n_perm)
        col_h, col_p = permutation_p(
            vals, [s.plate_col for s in placed], rng, n_perm)
        edge_labels = [
            "edge" if (s.plate_row in (rows[0], rows[-1])
                       or s.plate_col in (cols[0], cols[-1])) else "interior"
            for s in placed
        ]
        edge_h, edge_p = permutation_p(vals, edge_labels, rng, n_perm)
        findings.append({
            "quantity": what,
            "n_wells": len(placed),
            "row_p": round(row_p, 4),
            "col_p": round(col_p, 4),
            "edge_p": round(edge_p, 4),
        })
    return findings


# --------------------------------------------------------------------------
# grouping by metadata
# --------------------------------------------------------------------------


# Matches pre-flight. A library with no metadata row is bucketed here
# rather than dropped, so the group counts always add up to the plate.
NO_VALUE = "(no value)"

# Fields that are identifiers rather than categories. Grouping by one
# gives a row per library, which is the per-sample table again with a
# worse layout. A cap rather than a name list, because the field names
# are the user's and cannot be anticipated.
MAX_GROUPS = 50

# Declared once, because every row must carry every column. A row built
# from whatever that group happened to have gives a TSV whose header
# depends on which group sorted first.
COLUMNS = (
    ["run_id", "field", "value", "n_libraries", "n_controls", "n_analysed",
     "n_subjects", "n_empty", "reads_assigned", "median_depth", "min_depth",
     "max_depth"]
    + [f"reads_{k}" for k in FATES]
    + ["target_fraction", "plasmid_share", "detected", "not_detected",
       "assay_failed", "detection_rate"]
    + [f"wells_{k}" for k in INTERNAL]
    + ["deepest_negative", "note"]
)


def load_sample_metadata(path: str | None, run_id: str
                         ) -> tuple[dict[str, dict], list[str]]:
    """The resolved per-library table pre-flight wrote.

    Flat by design: the keying question - subject-level or
    library-level, and what a subject with three libraries means - was
    settled in pre-flight, and nothing here has to know how it was
    answered. This reads sample_id -> {field: value} and the field
    order, and that is the whole interface.
    """
    if not path or not os.path.exists(path):
        return {}, []
    rows = read_table(path)
    if not rows:
        return {}, []
    skip = {"run_id", "sample_id", "sample_type"}
    fields = [c for c in rows[0] if c not in skip]
    out: dict[str, dict] = {}
    for r in rows:
        # sample_id is unique within a run, not across runs, so a
        # multi-run table must be filtered or two plates' rows collide
        # on one id.
        if run_id and r.get("run_id") and r["run_id"] != run_id:
            continue
        out[r["sample_id"]] = {f: (r.get(f) or NO_VALUE) for f in fields}
    return out, fields


def _blank(v) -> str:
    """None as an empty cell, not as the string 'None'.

    A well with no assigned reads has no distance, and 'None' in a TSV
    column of integers parses as a string and poisons whatever reads it.
    """
    return "" if v is None else str(v)


def hist_quantile(hist: dict, q: float) -> int | None:
    """Quantile of an integer-keyed histogram whose values are counts.

    Keys arrive as strings, because JSON has no integer keys. Sorting
    them as strings would put 10 before 2 and make every quantile
    wrong - which is the kind of bug that produces a plausible number
    rather than an error, so the cast is not optional.
    """
    pairs = sorted((int(k), v) for k, v in hist.items())
    total = sum(v for _, v in pairs)
    if not total:
        return None
    target, seen = q * total, 0
    for d, v in pairs:
        seen += v
        if seen >= target:
            return d
    return pairs[-1][0]


# How much higher than the trough the far mode must peak before a dip is
# called a trough at all. Not a biological threshold - it is a noise
# criterion for reading the SHAPE, and it is needed because a single
# stray read one bin out otherwise makes a trough out of a clean tail.
# 2x is deliberately weak: the signal it exists to find was 18x on the
# validation set, and a weak criterion that is documented beats a strong
# one that quietly discards real second modes.
TROUGH_RISE = 2.0


def find_trough(hist: dict, rise: float = TROUGH_RISE) -> int | None:
    """Where the error mode ends and something that is not error begins.

    The per-well summary this feeds used to be a 95th percentile, and
    that could not work. The population worth finding is a fraction of a
    percent - measured at 0.70% on a plate spiked with 1.5% chimeras -
    so it sits above p95 AND above p99, and the column read identically
    in spiked and clean wells. A quantile cannot see a tail smaller than
    its own complement; no choice of quantile fixes that.

    What does work is to find the boundary in the RUN's own pooled
    histogram and then measure each well against it. The boundary is
    measured, not chosen: a clean run decays monotonically from the
    error mode and has no trough, so the answer is None and the column
    is empty - which is itself the finding. A run carrying chimeras
    dips and rises again, and the dip is the boundary.

    Returns the last distance belonging to the error mode, so "beyond
    the trough" means strictly greater. None when the histogram is
    monotone, too short to have a shape, or the rise is too weak to
    distinguish from a noisy tail.
    """
    if not hist:
        return None
    pairs = sorted((int(k), int(v)) for k, v in hist.items())
    # Gaps are counts of zero and are part of the shape - a histogram
    # that jumps 4 -> 0 -> 37 has a trough at 5, and reading only the
    # keys present would miss it entirely.
    lo, hi = pairs[0][0], pairs[-1][0]
    if hi - lo < 3:
        return None
    counts = dict(pairs)
    series = [counts.get(d, 0) for d in range(lo, hi + 1)]

    mode = max(range(len(series)), key=lambda i: series[i])
    i = mode
    while i + 1 < len(series) and series[i + 1] <= series[i]:
        i += 1
    if i >= len(series) - 1:
        return None                      # monotone decay: a clean run
    trough_i, trough_v = i, series[i]
    far = series[trough_i + 1:]
    if max(far) < max(trough_v, 1) * rise:
        return None                      # a noisy tail, not a second mode
    return lo + trough_i


def pool_distances(samples: dict) -> dict:
    """One distance distribution for the run, per outcome.

    Pooled rather than averaged, because the thing worth seeing is the
    SHAPE and a mean of medians destroys it. A clean run puts assigned
    reads in a single mode at the platform's error rate. A second mode
    further out is reads that are not errors - a chimera between two
    panel members sits at about half their divergence, is assigned to
    one of them, and passes the margin rule, so it appears in no
    outcome count anywhere. This histogram is the only place it shows.

    Controls are included. The pooled distribution is a property of the
    chemistry rather than a per-sample metric, an NTC's handful of reads
    cannot move it, and if a negative IS full of chimeras that is worth
    seeing rather than excluding.
    """
    pooled: dict[str, Counter] = defaultdict(Counter)
    for s in samples.values():
        for outcome, block in (s.dist or {}).items():
            for d, n in (block.get("histogram") or {}).items():
                pooled[outcome][int(d)] += n
    out = {}
    for outcome, h in sorted(pooled.items()):
        as_str = {str(d): h[d] for d in sorted(h)}
        trough = find_trough(as_str)
        beyond = sum(v for d, v in h.items() if trough is not None
                     and d > trough)
        total = sum(h.values())
        out[outcome] = {
            "reads": total,
            "median": hist_quantile(as_str, 0.5),
            "p95": hist_quantile(as_str, 0.95),
            "p99": hist_quantile(as_str, 0.99),
            "max": max(h),
            # None means the histogram is monotone - one mode, no
            # boundary to draw, which is what a clean run looks like.
            "trough": trough,
            "reads_beyond_trough": beyond if trough is not None else None,
            "fraction_beyond_trough": (f"{beyond / total:.5f}"
                                       if trough is not None and total
                                       else None),
            "histogram": as_str,
        }
    return out


def median(vals: list) -> float:
    if not vals:
        return 0.0
    v = sorted(vals)
    n = len(v)
    return float(v[n // 2]) if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def group_summary(
    samples: dict, meta: dict[str, dict], fields: list[str],
    call_rows: list[dict], cals: dict, plasmid_share: dict, run_id: str,
) -> list[dict]:
    """Consolidate the run by each grouping field in turn.

    What this is for is the comparison a run is often designed around -
    uncleaned against AMPure-cleaned PCR1, one extraction batch against
    another, replicates of one subject. The question is always the same
    shape: did the arms differ, and in what.

    It consolidates and does NOT pool. The arms of a designed comparison
    are the experiment; summing their reads into one number answers a
    question nobody asked and destroys the one that was. So each group
    gets its own row and they sit side by side.

    Two things that decide whether a row can be read at all:

    CONTROLS ARE EXCLUDED from every metric, and counted in their own
    column. A group's median depth is meaningless if an NTC is in the
    denominator - and worse, it moves when the number of NTCs changes,
    which makes two runs incomparable for a reason invisible in the
    table. What the negatives are for is a different question, so it
    gets its own column: the deepest negative in the group, which is
    what says whether that arm stayed clean.

    EVERY LIBRARY APPEARS. A library whose subject has no metadata row
    is grouped under '(no value)', never dropped, so n_libraries summed
    over a field's rows equals the plate. Without that a reader cannot
    tell a small group from a broken join.

    No test statistic. With two arms and a handful of libraries each, a
    p-value here would be read as an answer when it is barely evidence;
    the counts are the honest output, and a comparison worth testing is
    worth designing for.
    """
    # Every row, with no filter on reference_type - and the absence of
    # that filter is the point.
    #
    # calls.tsv is built from the references called in each well - the
    # targets, plus the control plasmids in a mock well, whose plasmids
    # are its content. So the decision about what counts as a result has
    # already been made upstream, once, where it belongs. Re-applying it here looked like
    # caution and was a bug: on a plate whose TEMPLATES are the control
    # plasmids, filtering to reference_type == "target" discards every
    # real call and reports a detection rate built from the handful of
    # genuine targets that happen to be present. Measured on such a run,
    # that read 4-11% where the plasmids themselves were never scored.
    detected_by: dict[str, Counter] = defaultdict(Counter)
    for r in call_rows:
        detected_by[r["sample_id"]][r["call"]] += 1

    rows: list[dict] = []
    # subject_id first, and unconditionally: replicates are a grouping
    # that needs no metadata file, and a run of replicates should be
    # summarisable with nothing supplied at all.
    all_fields = ["subject_id"] + [f for f in fields if f != "subject_id"]

    def value_of(sid: str, field: str) -> str:
        # subject_id is read off the sample itself, not the metadata
        # table, so grouping by subject works with no --metadata at all.
        # A library whose samplesheet declared no subject is its own
        # subject, which is what pre-flight already assumes.
        if field == "subject_id":
            return samples[sid].subject_id or sid
        return meta.get(sid, {}).get(field, NO_VALUE)

    for field in all_fields:
        members: dict[str, list[str]] = defaultdict(list)
        for sid in samples:
            members[value_of(sid, field)].append(sid)

        if field == "subject_id" and all(len(v) < 2 for v in members.values()):
            # One library per subject: grouping by it reproduces the
            # per-sample table, so it is not a grouping at all.
            continue
        if len(members) > MAX_GROUPS:
            # Still a row, with the full column set. A skipped field
            # that leaves no trace looks the same as a field that was
            # never in the file.
            rows.append({
                **{k: "" for k in COLUMNS},
                "run_id": run_id, "field": field, "value": "(not grouped)",
                "n_libraries": len(samples),
                "note": f"{len(members)} distinct values, above the limit of "
                        f"{MAX_GROUPS} - this looks like an identifier "
                        f"rather than a category",
            })
            continue

        for value in sorted(members):
            sids = members[value]
            ctrl = [s for s in sids if samples[s].sample_type in NEGATIVE_TYPES]
            keep = [s for s in sids
                    if samples[s].sample_type not in NEGATIVE_TYPES]

            fates = Counter()
            for s in keep:
                for k, v in read_fates(samples[s]).items():
                    fates[k] += v
            depths = [samples[s].depth for s in keep]
            total = sum(depths)
            calls = Counter()
            for s in keep:
                calls.update(detected_by.get(s, Counter()))
            states = Counter(cals[s].state for s in keep if s in cals)
            shares = [plasmid_share[s] for s in keep
                      if plasmid_share.get(s) is not None
                      and samples[s].sample_type not in PLASMID_WELLS]

            rows.append({
                "run_id": run_id,
                "field": field,
                "value": value,
                "n_libraries": len(sids),
                "n_controls": len(ctrl),
                "n_analysed": len(keep),
                "n_subjects": len({samples[s].subject_id or s for s in keep}),
                "n_empty": sum(1 for d in depths if d == 0),
                "reads_assigned": total,
                "median_depth": f"{median(depths):g}",
                "min_depth": min(depths) if depths else 0,
                "max_depth": max(depths) if depths else 0,
                **{f"reads_{k}": fates[k] for k in FATES},
                "target_fraction": (f"{fates['target'] / total:.4f}"
                                    if total else ""),
                "plasmid_share": (f"{median(shares):.4f}" if shares else ""),
                "detected": calls["detected"],
                "not_detected": calls["not_detected"],
                "assay_failed": calls["assay_failed"],
                "detection_rate": (
                    f"{calls['detected'] / (calls['detected'] + calls['not_detected']):.4f}"
                    if (calls["detected"] + calls["not_detected"]) else ""),
                **{f"wells_{k}": states[k] for k in INTERNAL},
                # The one control number that matters for a bench
                # comparison: did this arm's negatives stay clean.
                "deepest_negative": (max((samples[s].depth for s in ctrl),
                                         default="") if ctrl else ""),
                "note": "",
            })
    return rows


# --------------------------------------------------------------------------


def check_cohort_complete(samples: dict, preflight_dir: str, run_id: str,
                          demux_summary: str, allow_missing: bool) -> list[str]:
    """Every sample the run declared must reach synthesis.

    The one failure mode that produces a confident, wrong report rather
    than an error: some samples never get assigned, and synthesis
    happily summarises the ones that did. It happened - a Nextflow
    fan-out collapsed eight wells to one, and nothing downstream
    objected, because nothing downstream knew how many wells there were
    meant to be.

    Synthesis has always known. file_resolution.tsv is the run's
    declared cohort and load_counts() is what turned up; the two were
    simply never compared. They are now, and a shortfall is fatal.

    The one legitimate shortfall is a nanopore well that received no
    reads at all - an empty well demultiplexes to nothing and so has no
    FASTQ to assign. That is distinguishable, because the
    demultiplexing summary counts reads per sample: a sample it never
    saw is empty, a sample it wrote reads for and that then vanished is
    a bug. Only the second is fatal.
    """
    geom = os.path.join(preflight_dir or "", "file_resolution.tsv")
    if not os.path.exists(geom):
        return []

    declared = [r for r in read_table(geom)
                if not run_id or r.get("run_id", "") == run_id]
    if not declared:
        return []

    missing = [r for r in declared if r.get("sample_id", "") not in samples]
    if not missing:
        return []

    # Which of the missing ones the demultiplexer actually wrote reads
    # for. Keys may be the sample id or the filename-safe version of
    # it, so both are accepted.
    seen: dict[str, int] = {}
    if demux_summary and os.path.exists(demux_summary):
        try:
            with open(demux_summary, encoding="utf-8") as fh:
                seen = json.load(fh).get("per_sample", {}) or {}
        except (OSError, ValueError):
            seen = {}

    def had_reads(row: dict) -> bool:
        for key in (row.get("sample_id", ""), row.get("sanitised_id", "")):
            if key and seen.get(key, 0) > 0:
                return True
        return False

    # With no demultiplexing summary there is nothing to excuse a gap:
    # on Illumina pre-flight already proved every file exists, so a
    # missing sample can only be a lost task.
    lost = [r["sample_id"] for r in missing if not seen or had_reads(r)]
    empty = [r["sample_id"] for r in missing if seen and not had_reads(r)]

    if lost and not allow_missing:
        raise SystemExit(
            f"cohort incomplete: {len(declared)} samples declared for run "
            f"'{run_id or 'all'}', {len(samples)} reached synthesis.\n"
            f"  lost: {', '.join(sorted(lost)[:12])}"
            f"{' ...' if len(lost) > 12 else ''}\n"
            "These samples have reads but no *.counts.tsv, so assignment "
            "never ran for them.\n"
            "On Nextflow this is usually a fan-out capped by a single-item "
            "queue channel - see tests/test_nextflow.py::"
            "test_fanout_is_not_capped.\n"
            "Refusing to report on a partial cohort. Pass --allow-missing "
            "to override, knowing the run QC below will be computed from "
            "an incomplete plate."
        )
    return sorted(empty) + (sorted(lost) if allow_missing else [])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assign-dir", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--prep-dir",
                    help="directory of *.prep.json, so the funnel can start "
                         "at the reads that came off the sequencer")
    ap.add_argument("--demux-summary",
                    help="demux_summary.json, for the nanopore path")
    ap.add_argument("--split-summary",
                    help="split_summary.json, so the funnel starts at reads "
                         "off the sequencer rather than at amplicon units")
    ap.add_argument(
        "--sample-metadata",
        help="preflight/sample_metadata.tsv - the resolved one-row-per-"
             "library join. Given it, synthesis writes group_summary.tsv: "
             "the run consolidated by each grouping field in turn. Absent, "
             "grouping by subject_id still happens, since that is declared "
             "in the samplesheet")
    ap.add_argument("--preflight-dir",
                    help="for plate geometry derived at pre-flight")
    ap.add_argument("--outdir", default="synthesis")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--min-quant-reads", type=int, default=10)
    ap.add_argument("--min-depth", type=int, default=100)
    # Detection is strictly above the limit (n > limit), and the floor is
    # one of the terms the limit is the maximum of. So a well needs MORE
    # than this many reads - with the default of 2, three - however clean
    # the blanks were. The help says that rather than "below which", which
    # described a rule one read looser than the one applied.
    ap.add_argument("--absolute-floor", type=int, default=2,
                    help="the smallest detection limit any reference can "
                         "have: a detection needs MORE reads than this, "
                         "however clean the blanks were")
    # Default matches params.slope_tolerance, which the workflow passes;
    # this default only matters when the script is run by hand.
    ap.add_argument("--slope-tolerance", type=float, default=0.15,
                    help="how far a well's standard-curve slope may sit "
                         "from 1 before it is noted. Copy numbers are only "
                         "withheld outside the wider plausible band "
                         f"{SLOPE_PLAUSIBLE[0]:g}-{SLOPE_PLAUSIBLE[1]:g}")
    ap.add_argument("--spiked", action="store_true",
                    help="the quantification standards were added to every "
                         "sample on this run. Calibrates each well against "
                         "them, gives copy numbers, and treats a well "
                         "missing its standards as FAILED (assay_failed). "
                         "Without it the run is unspiked: no copy numbers, "
                         "every negative unverified, and any standard reads "
                         "are reported as carryover")
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--plasmid-warn", type=float, default=0.05,
                    help="warn when this fraction of a well's amplicon reads "
                         "come from control plasmid")
    ap.add_argument("--allow-missing", action="store_true",
                    help="report even though samples are missing from the "
                         "cohort. Off by default: a partial plate gives a "
                         "confident wrong answer rather than an error.")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    refs = read_table(args.reference)
    by_id = {r["seq_id"]: r for r in refs}
    quants = [r for r in refs if r.get("type") == "quant"]
    targets = [r for r in refs if r.get("type") == "target"]
    controls = [r for r in refs if r.get("type") == "control"]
    # What is called depends on the WELL, not the run. A mock is built
    # from tagged control plasmids, so in a mock well the plasmids are the
    # content: called, counted and quantified like targets. In any other
    # well the same tagged reads are contamination - never called as a
    # finding, measured instead as plasmid_share and plasmid_copies. A
    # run-wide switch (--controls-as-targets, removed) could only get one
    # of the two right on a plate that holds both.
    def called_in(s_) -> list:
        return targets + controls if s_.sample_type in PLASMID_WELLS \
            else targets

    samples = load_counts(args.assign_dir)
    if not samples:
        raise SystemExit(f"no *.counts.tsv found in {args.assign_dir}")

    # Before anything is compared: put the declared ids back.
    #
    # Assignment keys its output on the filename-safe form of the sample
    # id, because that is what reached it - demultiplexing writes one
    # FASTQ per sample and the workflow reads the id off the file. Every
    # join below is against the samplesheet's ORIGINAL id, so this has to
    # happen before the cohort check and before anything is attached to a
    # sample, or all of it silently matches nothing.
    renamed = restore_declared_ids(samples, args.preflight_dir, args.run_id)
    if renamed and not args.quiet:
        print(f"restored {len(renamed)} declared sample id(s) from their "
              f"filename-safe form, e.g. {renamed[0][0]} -> {renamed[0][1]}")

    # Before anything is computed: is this the whole plate?
    empty_wells = check_cohort_complete(
        samples, args.preflight_dir, args.run_id,
        args.demux_summary, args.allow_missing)

    for r in read_table(args.samplesheet):
        s = samples.get(r.get("sample_id", ""))
        if s:
            s.sample_type = r.get("sample_type", "unknown")
            s.subject_id = (r.get("subject_id") or "").strip()
            s.index_i5 = r.get("index_i5", "")
            s.index_i7 = r.get("index_i7", "")

    # Plate geometry, as pre-flight derived it from the index primer names.
    # Without it the position sensors are simply skipped - they are not
    # guessed at from sample order, which would invent a layout.
    geom = os.path.join(args.preflight_dir or "", "file_resolution.tsv")
    if os.path.exists(geom):
        for r in read_table(geom):
            s = samples.get(r.get("sample_id", ""))
            if s:
                s.plate_row = r.get("plate_row", "")
                # file_resolution is pre-flight's resolved view, so it
                # carries the subject even when the samplesheet column
                # was absent and pre-flight defaulted it to the library.
                s.subject_id = (r.get("subject_id") or s.subject_id or "")
                try:
                    s.plate_col = int(r.get("plate_col") or 0)
                except ValueError:
                    s.plate_col = 0

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    # -- calibration ------------------------------------------------------
    spiked_run = args.spiked
    spike_decision = spike_declaration(
        spiked_run, samples, quants, args.min_quant_reads)
    if spiked_run:
        cals = {sid: calibrate(s, quants, args.min_quant_reads,
                               args.slope_tolerance, True)
                for sid, s in samples.items()}
    else:
        cals = {sid: undeclared_standards(s, quants, args.min_quant_reads)
                for sid, s in samples.items()}

    # -- detection limits -------------------------------------------------
    seq_ids = [r["seq_id"] for r in targets + controls]
    lods = detection_limits(samples, seq_ids, args.absolute_floor,
                            args.min_depth)

    # -- composition, always both denominators ----------------------------
    # Two questions, two denominators, and a switch must not silently
    # decide which one a column means:
    #   within a group - which target, or which plasmid, and in what ratio
    #   plasmid share  - how much of this well is plasmid at all
    # The second is the contamination measure and stays meaningful even
    # when the plasmids are the experiment.
    target_ids = [t["seq_id"] for t in targets]
    control_ids = [c["seq_id"] for c in controls]
    group_total, plasmid_share = {}, {}
    for sid, s_ in samples.items():
        tsum = sum(s_.counts.get(i, 0) for i in target_ids)
        csum = sum(s_.counts.get(i, 0) for i in control_ids)
        group_total[sid] = {"target": tsum, "control": csum}
        plasmid_share[sid] = (csum / (tsum + csum)) if (tsum + csum) else None

    # Contaminating plasmid in copies, on the same scale as the targets -
    # only where there is a scale (a calibrated well) and only where the
    # plasmid IS contamination (not a mock). Summed per reference rather
    # than inverted once from the pooled reads, because the curve is a
    # power law and the inversion does not distribute over a sum.
    plasmid_copies: dict[str, float | None] = {}
    for sid, s_ in samples.items():
        if s_.sample_type in PLASMID_WELLS or cals[sid].state != "spiked":
            plasmid_copies[sid] = None
            continue
        tot = 0.0
        for c in controls:
            n = s_.counts.get(c["seq_id"], 0)
            raw = cals[sid].copies(n) if n > 0 else None
            if raw is not None:
                cf = float((c.get("corr_factor") or "").strip()
                           or DEFAULT_CORR_FACTOR)
                tot += raw / cf
        plasmid_copies[sid] = tot

    # -- calls ------------------------------------------------------------
    rows = []
    summary_calls: Counter = Counter()
    for sid in sorted(samples):
        s = samples[sid]
        cal = cals[sid]
        for t in called_in(s):
            tid = t["seq_id"]
            kind = t.get("type") or "target"
            n = s.counts.get(tid, 0)
            limit, basis = limit_for(tid, s, lods)

            if s.depth < args.min_depth:
                call = "assay_failed"
                why = (f"only {s.depth} sequences reached assignment, below "
                       f"the minimum of {args.min_depth}")
            elif cal.state == "spike_failed":
                call = "assay_failed"
                why = cal.note
            elif n > limit:
                call = "detected"
                why = f"{n} reads against a limit of {limit}, set by {basis}"
            else:
                call = "not_detected"
                why = f"{n} reads, at or below the limit of {limit}"
                if cal.state == "unspiked":
                    why += " - and unverified, since no standard was recovered"

            raw = cal.copies(n) if call == "detected" else None
            if call == "detected" and cal.state == "slope_implausible":
                # Said on the row, because an empty copies column on a
                # detected row otherwise looks like a bug rather than a
                # decision.
                why += (f" - no copy number: calibration slope "
                        f"{cal.slope:.2f} is outside the plausible "
                        f"{SLOPE_PLAUSIBLE[0]:g}-{SLOPE_PLAUSIBLE[1]:g}")
            # corr_factor defaults to 1.0: unity is the sensible assumption
            # when efficiency has not been measured, and it keeps every
            # reference on the same arithmetic footing. The source is
            # carried alongside the value, because "measured and found to
            # be 1" and "assumed to be 1" produce the same number from
            # different amounts of evidence.
            supplied = (t.get("corr_factor") or "").strip()
            cf = float(supplied) if supplied else DEFAULT_CORR_FACTOR
            copies = (raw / cf) if raw is not None else None

            rows.append({
                "run_id": args.run_id,
                "sample_id": sid,
                "sample_type": s.sample_type,
                "well": well_of(s.plate_row, s.plate_col),
                "seq_id": tid,
                "reference_type": kind,
                "locus": t.get("locus", ""),
                "purpose": t.get("purpose") or SPECIATION,
                "taxon": t.get("taxon") or tid,
                "reads": n,
                "depth": s.depth,
                # Proportion within this reference's own group, so a
                # plasmid proportion is never diluted by target reads and
                # the reverse.
                "fraction_of_group": (
                    "" if not group_total[sid].get(kind)
                    else f"{n / group_total[sid][kind]:.4f}"),
                # And the contamination number, identical on every row of
                # a sample: what share of its amplicons are plasmid.
                "plasmid_share": ("" if plasmid_share[sid] is None
                                  else f"{plasmid_share[sid]:.4f}"),
                # The same contamination in copies, on a spiked run. Empty
                # in a mock, where the plasmids are the content and each
                # has its own copies on its own row.
                "plasmid_copies": ("" if plasmid_copies[sid] is None
                                   else f"{plasmid_copies[sid]:.3f}"),
                "detection_limit": limit,
                "call": call,
                "internal_control": cal.state,
                "copies_per_unit_dna": ("" if copies is None
                                        else f"{copies:.3f}"),
                "corr_factor": f"{cf:g}",
                "corr_factor_source": ("panel" if supplied else "default"),
                "calibration_slope": ("" if not math.isfinite(cal.slope)
                                      else f"{cal.slope:.3f}"),
                "basis": why,
            })
            summary_calls[call] += 1

    with open(os.path.join(args.outdir, "calls.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # -- consolidation by metadata group ----------------------------------
    # Written whenever there is a grouping to write: a metadata file, or
    # simply a subject with more than one library. Absent otherwise,
    # rather than an empty file, so its presence means something.
    meta, meta_fields = load_sample_metadata(args.sample_metadata,
                                             args.run_id)
    group_rows = group_summary(samples, meta, meta_fields, rows, cals,
                               plasmid_share, args.run_id)
    if group_rows:
        with open(os.path.join(args.outdir, "group_summary.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                               lineterminator="\n", restval="")
            w.writeheader()
            w.writerows(group_rows)

    with open(os.path.join(args.outdir, "calibration.tsv"), "w",
              newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "state", "slope", "intercept", "reads_hi",
                    "conc_hi", "reads_lo", "conc_lo", "note"])
        for sid in sorted(cals):
            c = cals[sid]
            w.writerow([
                sid, c.state,
                "" if not math.isfinite(c.slope) else f"{c.slope:.4f}",
                "" if not math.isfinite(c.intercept) else f"{c.intercept:.4f}",
                c.reads_hi, c.conc_hi, c.reads_lo, c.conc_lo, c.note,
            ])

    # -- read fates -------------------------------------------------------
    # Pooled first, because the per-well column is measured against the
    # RUN's boundary rather than its own. A trough found inside one well
    # would move from well to well and the column would not be
    # comparable across the plate - which is the only thing it is for.
    pooled_dist = pool_distances(samples)
    trough = (pooled_dist.get("assigned") or {}).get("trough")

    fate_rows, fate_totals, unbalanced = [], Counter(), []
    for sid in sorted(samples):
        s = samples[sid]
        f = read_fates(s)
        tot = sum(f.values())
        if s.depth and tot != s.depth:
            unbalanced.append(sid)
        for k, v in f.items():
            fate_totals[k] += v
        # Distance summary for the reads that were ASSIGNED. Here rather
        # than in a diagnostic file because this is the one column set
        # that can betray a chimera: the read is assigned, to one parent,
        # with a passing margin, so every other count in this table looks
        # clean. Only how far it sat from its reference does not.
        #
        # The share beyond the run's trough, NOT a quantile. A p95 column
        # shipped here first and was useless: the population it exists to
        # find was 0.70% of reads, so it sat above p95 and above p99 both,
        # and the column read 2 in every well whether spiked or clean. A
        # quantile cannot see a tail smaller than its own complement.
        assigned_d = (s.dist or {}).get("assigned") or {}
        hist = assigned_d.get("histogram") or {}
        n_assigned = sum(int(v) for v in hist.values())
        over = (sum(int(v) for d, v in hist.items() if int(d) > trough)
                if trough is not None else None)
        fate_rows.append({
            "run_id": args.run_id, "sample_id": sid,
            "sample_type": s.sample_type, "well": well_of(s.plate_row,
                                                          s.plate_col),
            "depth": s.depth, **f,
            "target_fraction": (f"{f['target']/s.depth:.6f}"
                                if s.depth else ""),
            "assigned_dist_median": _blank(hist_quantile(hist, 0.5)),
            "assigned_dist_max": _blank(max((int(k) for k in hist),
                                            default=None)),
            # Empty for the whole plate when the run is monotone, which
            # is the clean case and is meant to look like nothing.
            "assigned_dist_over_trough": _blank(over),
            "assigned_dist_over_trough_frac": (
                f"{over / n_assigned:.5f}"
                if over is not None and n_assigned else ""),
        })
    with open(os.path.join(args.outdir, "read_fates.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fate_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(fate_rows)

    prep = load_prep(args.prep_dir)
    # Prep summaries are keyed by the filename-safe id the workflow ran
    # with, exactly as the counts were; bring them back to the declared
    # id the same way so the outlier list names samples as declared.
    restore_declared_ids(prep, args.preflight_dir, args.run_id)
    demux = None
    if args.demux_summary and os.path.exists(args.demux_summary):
        with open(args.demux_summary, encoding="utf-8") as fh:
            demux = json.load(fh)
    split = None
    if args.split_summary and os.path.exists(args.split_summary):
        with open(args.split_summary, encoding="utf-8") as fh:
            split = json.load(fh)
    funnel = read_funnel(samples, prep, demux, fate_totals, split)

    # -- QC ---------------------------------------------------------------
    # Mock wells are left out of every contamination figure: plasmid is
    # what they are made of, and counting it would flag every mock and
    # inflate the run's share with the experiment itself.
    contam = [sid for sid in samples
              if samples[sid].sample_type not in PLASMID_WELLS]
    flagged = sorted(
        ({"sample_id": sid,
          "well": well_of(samples[sid].plate_row, samples[sid].plate_col),
          "sample_type": samples[sid].sample_type,
          "plasmid_share": round(plasmid_share[sid], 4),
          "plasmid_copies": (None if plasmid_copies[sid] is None
                             else round(plasmid_copies[sid], 3))}
         for sid in contam
         if plasmid_share[sid] is not None
         and plasmid_share[sid] > args.plasmid_warn),
        key=lambda d: -d["plasmid_share"])
    all_t = sum(group_total[sid]["target"] for sid in contam)
    all_c = sum(group_total[sid]["control"] for sid in contam)
    plasmid_qc = {
        "controls_called_in": sorted(PLASMID_WELLS),
        "warn_above": args.plasmid_warn,
        "run_plasmid_share": (round(all_c / (all_t + all_c), 4)
                              if (all_t + all_c) else None),
        "n_wells_above_threshold": len(flagged),
        "wells_above_threshold": flagged[:20],
        "note": (
            "Share of a well's panel amplicons coming from tagged control "
            "plasmid, in every well that is not a mock - there the plasmid "
            "is contamination, recognised by its tag, and never called as "
            "a finding. On a spiked run it is also given in copies "
            "(plasmid_copies in calls.tsv). Mock wells are left out: they "
            "are built from the plasmids, which are called in them like "
            "targets."
        ),
    }
    taxa_for_qc = [t["seq_id"] for t in targets][:6]
    geometry = plate_qc(samples, taxa_for_qc, rng, args.permutations)
    depths = sorted(s.depth for s in samples.values())
    slopes = [c.slope for c in cals.values() if math.isfinite(c.slope)]

    qc = {
        "run_id": args.run_id,
        "n_samples": len(samples),
        # Declared minus reached. Empty nanopore wells are legitimate and
        # named here rather than left as a silent difference between the
        # samplesheet and the report; anything else has already been
        # fatal by this point unless --allow-missing was given.
        "samples_absent": empty_wells,
        "cohort_complete": not empty_wells,
        "call_summary": dict(summary_calls),
        "call_meaning": CALL_MEANING,
        "assignment_distance": pooled_dist,
        "assignment_distance_note": (
            "Edit distance from each read to the reference it was "
            "assigned to, pooled over the run and weighted by reads. "
            "Read the shape, not the median. One mode at the platform's "
            "error rate is a clean run. A second mode further out is "
            "reads that are not errors: a chimera between two panel "
            "members sits at about half their divergence, is assigned to "
            "one of them, and passes the margin rule - so it is counted "
            "as a confident call and appears in no outcome tally "
            "anywhere. This is the only place it is visible. "
            "`trough` is where the first mode ends, FOUND in this run's "
            "own histogram rather than chosen: a monotone histogram has "
            "none and reports null, which is the clean case. Per-well "
            "shares beyond it are in read_fates.tsv. Nothing is "
            "filtered on this."
        ),
        "groups": {
            "fields": sorted({r["field"] for r in group_rows}),
            "n_rows": len(group_rows),
            "metadata_supplied": bool(meta),
            "note": (
                "Consolidated per group in group_summary.tsv. Negative "
                "controls are excluded from every metric and counted "
                "separately, and a library with no metadata row is "
                "grouped under '(no value)' rather than dropped, so "
                "n_libraries sums to the plate within each field."
            ),
        },
        "run_was_spiked": spiked_run,
        "spike_declaration": spike_decision,
        "internal_control_summary": dict(
            Counter(c.state for c in cals.values())),
        "internal_control_meaning": INTERNAL_MEANING,
        "depth": {
            "min": depths[0] if depths else 0,
            "median": depths[len(depths) // 2] if depths else 0,
            "max": depths[-1] if depths else 0,
            # A run is not described by its mean depth. The spread is what
            # decides whether the shallowest wells could have detected
            # anything at all.
            "fold_spread": (round(depths[-1] / max(1, depths[0]), 1)
                            if depths else 0),
            "below_minimum": sum(1 for d in depths if d < args.min_depth),
        },
        "corr_factor": {
            "default_applied": DEFAULT_CORR_FACTOR,
            "references_using_default": sum(
                1 for t in targets if not (t.get("corr_factor") or "").strip()
            ),
            "references_with_measured_value": sum(
                1 for t in targets if (t.get("corr_factor") or "").strip()
            ),
            "note": (
                "Where no efficiency was measured, 1.0 is assumed. Copy "
                "numbers for those references are therefore relative to the "
                "standards' efficiency, not corrected to it - a reference "
                "that amplifies twice as well reads twice as abundant."
            ),
        },
        "calibration_slope": {
            "n": len(slopes),
            "median": round(sorted(slopes)[len(slopes) // 2], 3) if slopes
            else None,
            "outside_tolerance": sum(
                1 for x in slopes if abs(x - 1) > args.slope_tolerance),
            "tolerance": args.slope_tolerance,
            # The stricter question: how many wells had copy numbers
            # withheld altogether, not merely flagged.
            "plausible_band": list(SLOPE_PLAUSIBLE),
            "outside_plausible_band": sum(
                1 for c in cals.values() if c.state == "slope_implausible"),
        },
        "detection_limits": lods,
        "read_funnel": funnel,
        "read_fates": {
            "totals": dict(fate_totals),
            "meaning": FATE_MEANING,
            # The fates have to account for every read that reached
            # assignment. A sample where they do not means a stage lost
            # reads without reporting it.
            "partition_balances": not unbalanced,
            "unbalanced_samples": unbalanced[:5],
            "target_fraction_median": (
                sorted(f["target"] / f["depth"] for f in fate_rows
                       if f["depth"])[
                    len([f for f in fate_rows if f["depth"]]) // 2]
                if any(f["depth"] for f in fate_rows) else None
            ),
        },
        # What produced the reads, carried up from whichever stage still
        # had the original headers. Recorded because a basecalling model
        # changes the error profile, and therefore barcode recovery and
        # what reaches each sample: two runs basecalled differently can
        # differ in composition for reasons that are not biological, and
        # nothing else in these outputs would show it.
        "basecalling": (
            (split or {}).get("basecalling")
            or (demux or {}).get("basecalling")
            or {}
        ),
        # Plasmid contamination, reported whatever the flag says. A well
        # whose amplicons are largely plasmid is suspect even on a run
        # where the plasmids ARE the subject: it means something carried
        # across, not that the titration worked.
        "plasmid_contamination": plasmid_qc,
        "plate_geometry": geometry,
        "plate_geometry_note": (
            "Wells hold unrelated samples, so position should explain "
            "nothing. A small p-value is evidence about the plate - a "
            "channel, a gradient, evaporation at the edges - not about the "
            "samples. These are unadjusted for multiple testing; with one "
            "test per taxon per axis, expect roughly one in twenty to fall "
            "below 0.05 by chance alone."
        ),
    }
    with open(os.path.join(args.outdir, "run_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2)

    if not args.quiet:
        print(f"{args.run_id or 'run'}: {len(samples)} samples, "
              f"{len(targets)} targets")
        for k in CALLS:
            n = summary_calls.get(k, 0)
            if n:
                print(f"    {k:14s} {n:6,}   {CALL_MEANING[k]}")
        ic = qc["internal_control_summary"]
        print("    internal control: "
              + ", ".join(f"{k} {v}" for k, v in sorted(ic.items())))
        for st in funnel["stages"]:
            pct = ("" if st["percent_of_start"] is None
                   else f"  {st['percent_of_start']:5.1f}%")
            print(f"    {st['stage']:<30s} {st['reads']:>10,}{pct}")
            for k, v in st["lost_here"].items():
                print(f"        lost to {k:<22s} {v:>10,}")
        for n in funnel["notes"]:
            print(f"    ! {n}")
        d = qc["depth"]
        print(f"    depth {d['min']:,}-{d['max']:,} "
              f"({d['fold_spread']}x spread), "
              f"{d['below_minimum']} below minimum")
        flagged = [
            g for g in geometry
            if min(g["row_p"], g["col_p"], g["edge_p"]) < 0.01
        ]
        if flagged:
            print(f"    plate structure in {len(flagged)} quantities:")
            for g in flagged:
                print(f"       {g['quantity']:16s} row p={g['row_p']:.4f}  "
                      f"col p={g['col_p']:.4f}  edge p={g['edge_p']:.4f}")
        else:
            print("    no plate-position structure detected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
