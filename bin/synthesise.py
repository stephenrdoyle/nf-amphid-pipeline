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
from ampseqlib import read_table  # noqa: E402

CALLS = ("detected", "not_detected", "assay_failed")

CALL_MEANING = {
    "detected": "above the plate's measured detection limit",
    "not_detected": "below it, in a well whose spike-ins worked",
    "assay_failed": "the well cannot support either answer",
}

INTERNAL = ("spiked", "unspiked", "spike_failed")

INTERNAL_MEANING = {
    "spiked": "both standards recovered, within the expected ratio",
    "unspiked": "no standard present - a negative here is unverified",
    "spike_failed": "standards present but wrong, so the well is unreliable",
}

NEGATIVE_TYPES = {"neg_ntc", "neg_extraction", "blank", "negative"}

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
    index_i5: str = ""
    index_i7: str = ""
    plate_row: str = ""
    plate_col: int = 0
    counts: dict = field(default_factory=dict)
    total_target: int = 0
    depth: int = 0
    outcomes: dict = field(default_factory=dict)
    by_ref_type: dict = field(default_factory=dict)


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
            return None
        if not (self.slope and math.isfinite(self.slope)):
            return None
        return 10 ** ((math.log10(reads) - self.intercept) / self.slope)


def run_was_spiked(samples: dict, quants: list, min_reads: int) -> bool:
    """Was this run spiked at all?

    A single well cannot answer this. No standard reads means either the
    spike was omitted or it failed, and those call for opposite
    conclusions: the first makes a negative unverified, the second makes
    the well unusable. The run can answer it - if most wells recovered
    their standards, the few that did not have failed, whereas if none did,
    the protocol simply did not include a spike.
    """
    if not quants:
        return False
    hits = 0
    for s in samples.values():
        if any(s.counts.get(q["seq_id"], 0) >= min_reads for q in quants):
            hits += 1
    return hits >= max(1, 0.25 * len(samples))


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
                "no standard recovered, on a run where most wells recovered "
                "theirs - this well failed rather than going unspiked, so "
                "nothing in it can be called negative"
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
    samples: dict[str, Sample], seq_ids: list[str], floor: int
) -> dict[str, dict]:
    """What a well can show while containing nothing.

    Derived per reference from the run's own negative controls, which is
    the only place the answer is measurable. Two bounds, because the two
    things that put reads in an empty well behave differently: index
    hopping delivers a roughly fixed number of reads regardless of how deep
    the receiving well is, while carryover and cross-contamination scale
    with it.
    """
    blanks = [s for s in samples.values() if s.sample_type in NEGATIVE_TYPES]
    out: dict[str, dict] = {}
    for sid in seq_ids:
        obs = [s.counts.get(sid, 0) for s in blanks]
        fracs = [
            (s.counts.get(sid, 0) / s.depth) if s.depth else 0.0
            for s in blanks
        ]
        out[sid] = {
            "n_blanks": len(blanks),
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
    ap.add_argument("--preflight-dir",
                    help="for plate geometry derived at pre-flight")
    ap.add_argument("--outdir", default="synthesis")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--min-quant-reads", type=int, default=10)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--absolute-floor", type=int, default=2,
                    help="reads below which nothing is ever called detected")
    ap.add_argument("--slope-tolerance", type=float, default=0.25)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    refs = read_table(args.reference)
    by_id = {r["seq_id"]: r for r in refs}
    quants = [r for r in refs if r.get("type") == "quant"]
    targets = [r for r in refs if r.get("type") == "target"]
    controls = [r for r in refs if r.get("type") == "control"]

    samples = load_counts(args.assign_dir)
    if not samples:
        raise SystemExit(f"no *.counts.tsv found in {args.assign_dir}")

    for r in read_table(args.samplesheet):
        s = samples.get(r.get("sample_id", ""))
        if s:
            s.sample_type = r.get("sample_type", "unknown")
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
                try:
                    s.plate_col = int(r.get("plate_col") or 0)
                except ValueError:
                    s.plate_col = 0

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    # -- calibration ------------------------------------------------------
    spiked_run = run_was_spiked(samples, quants, args.min_quant_reads)
    cals = {sid: calibrate(s, quants, args.min_quant_reads,
                           args.slope_tolerance, spiked_run)
            for sid, s in samples.items()}

    # -- detection limits -------------------------------------------------
    seq_ids = [r["seq_id"] for r in targets + controls]
    lods = detection_limits(samples, seq_ids, args.absolute_floor)

    # -- calls ------------------------------------------------------------
    rows = []
    summary_calls: Counter = Counter()
    for sid in sorted(samples):
        s = samples[sid]
        cal = cals[sid]
        for t in targets:
            tid = t["seq_id"]
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
                "locus": t.get("locus", ""),
                "purpose": t.get("purpose") or SPECIATION,
                "taxon": t.get("taxon") or tid,
                "reads": n,
                "depth": s.depth,
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
    fate_rows, fate_totals, unbalanced = [], Counter(), []
    for sid in sorted(samples):
        s = samples[sid]
        f = read_fates(s)
        tot = sum(f.values())
        if s.depth and tot != s.depth:
            unbalanced.append(sid)
        for k, v in f.items():
            fate_totals[k] += v
        fate_rows.append({
            "run_id": args.run_id, "sample_id": sid,
            "sample_type": s.sample_type, "well": well_of(s.plate_row,
                                                          s.plate_col),
            "depth": s.depth, **f,
            "target_fraction": (f"{f['target']/s.depth:.6f}"
                                if s.depth else ""),
        })
    with open(os.path.join(args.outdir, "read_fates.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fate_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(fate_rows)

    prep = load_prep(args.prep_dir)
    demux = None
    if args.demux_summary and os.path.exists(args.demux_summary):
        with open(args.demux_summary, encoding="utf-8") as fh:
            demux = json.load(fh)
    funnel = read_funnel(samples, prep, demux, fate_totals)

    # -- QC ---------------------------------------------------------------
    taxa_for_qc = [t["seq_id"] for t in targets][:6]
    geometry = plate_qc(samples, taxa_for_qc, rng, args.permutations)
    depths = sorted(s.depth for s in samples.values())
    slopes = [c.slope for c in cals.values() if math.isfinite(c.slope)]

    qc = {
        "run_id": args.run_id,
        "n_samples": len(samples),
        "call_summary": dict(summary_calls),
        "call_meaning": CALL_MEANING,
        "run_was_spiked": spiked_run,
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
