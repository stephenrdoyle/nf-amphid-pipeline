#!/usr/bin/env python3
"""
Amplification efficiency from mock communities.

A read count is not an amount. Two templates present in equal numbers do
not yield equal reads: they prime differently, they differ in length and
base composition, and they compete for the same polymerase. `corr_factor`
is the number that converts one into the other, and this is the only place
it can be measured - from mixtures whose composition was decided before
sequencing.

THE QUESTION THAT COMES FIRST

Whether a single number per species is even the right shape of answer.

If bias is a fixed property of each template, then for any pair

    log(observed A/B) = log(input A/B) + log(eA/eB)

which is a line of slope 1, and its intercept is the correction. If
instead the templates compete, a dominant species is capped and the
observed ratio compresses toward the middle: the slope falls below 1, and
no constant is correct anywhere except the one mixture it was fitted at.
Real samples are mostly one dominant species plus a minority, so a
constant fitted at 50/50 fails exactly where it is needed.

So the slope is tested before any correction is emitted, and a constant is
emitted only when the data support one. Where they do not, this reports
the fitted relationship and declines to collapse it to a number - because
a wrong correction applied confidently is worse than none, and it would
inflate precisely the minority calls that matter.

Confidence intervals come from resampling wells rather than from a
distributional assumption, matching how the plate-effect test works
elsewhere in this pipeline. A correction without an interval cannot be
told from noise.

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
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table  # noqa: E402

MODELS = ("constant", "composition_dependent", "insufficient_data")

MODEL_MEANING = {
    "constant": "a single corr_factor is supported and is given",
    "composition_dependent": (
        "efficiency changes with the mixture, so no single corr_factor is "
        "valid - the fitted slope is reported instead"
    ),
    "insufficient_data": (
        "too few mixtures or replicates to decide; corr_factor left unset "
        "and the default of 1.0 will apply"
    ),
}


# --------------------------------------------------------------------------
# small statistics, resampled rather than assumed
# --------------------------------------------------------------------------


def geomean(xs: list[float]) -> float:
    xs = [x for x in xs if x > 0]
    if not xs:
        return float("nan")
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def boot_ci(
    values: list[float], stat, rng: random.Random, n: int = 2000,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Percentile bootstrap over whole observations.

    Resampling the wells, not the reads: the uncertainty that matters is
    how much the answer moves between pipettings and PCRs, not the
    multinomial noise within one library.
    """
    if len(values) < 2:
        return (float("nan"), float("nan"))
    out = []
    k = len(values)
    for _ in range(n):
        out.append(stat([values[rng.randrange(k)] for _ in range(k)]))
    out = sorted(x for x in out if x == x)
    if not out:
        return (float("nan"), float("nan"))
    lo = out[max(0, int(alpha / 2 * len(out)))]
    hi = out[min(len(out) - 1, int((1 - alpha / 2) * len(out)))]
    return lo, hi


def ols(points: list[tuple[float, float]]) -> tuple[float, float]:
    n = len(points)
    if n < 2:
        return (float("nan"), float("nan"))
    mx = sum(p[0] for p in points) / n
    my = sum(p[1] for p in points) / n
    sxx = sum((p[0] - mx) ** 2 for p in points)
    if sxx <= 0:
        return (float("nan"), float("nan"))
    sxy = sum((p[0] - mx) * (p[1] - my) for p in points)
    slope = sxy / sxx
    return slope, my - slope * mx


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def load_mocks(path: str) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    for r in read_table(path):
        sid = r.get("sample_id", "")
        if not sid:
            continue
        try:
            amount = float(r["amount"])
        except (KeyError, ValueError):
            continue
        if amount <= 0:
            continue
        out[sid].append({
            "seq_id": r.get("seq_id", ""),
            "amount": amount,
            "units": (r.get("amount_units") or "copies").lower(),
            "mix": r.get("mix") or sid,
            "replicate": r.get("replicate", ""),
            "pcr_batch": r.get("pcr_batch", ""),
        })
    return dict(out)


def load_observed(assign_dir: str) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = defaultdict(dict)
    for path in sorted(glob.glob(os.path.join(assign_dir, "*.counts.tsv"))):
        for r in read_table(path):
            out[r["sample_id"]][r["seq_id"]] = int(r["reads"])
    return dict(out)


# --------------------------------------------------------------------------
# efficiency
# --------------------------------------------------------------------------


def per_well_efficiency(
    rows: list[dict], obs: dict[str, int], quants: set[str], min_reads: int,
) -> dict[str, float]:
    """Efficiency of each species in one well, relative to the standards.

    Anchoring to the spike-ins rather than to another target is what makes
    the number transferable: field samples are measured on the same
    standards, so an efficiency expressed against them can be applied
    there. An efficiency expressed against another target only holds while
    that target is present.
    """
    anchor = [r for r in rows if r["seq_id"] in quants]
    anchor_obs = sum(obs.get(r["seq_id"], 0) for r in anchor)
    anchor_in = sum(r["amount"] for r in anchor)
    if not anchor or anchor_obs < min_reads or anchor_in <= 0:
        return {}
    scale = (anchor_obs / anchor_in)
    out = {}
    for r in rows:
        if r["seq_id"] in quants:
            continue
        n = obs.get(r["seq_id"], 0)
        if n < min_reads:
            continue
        out[r["seq_id"]] = (n / r["amount"]) / scale
    return out


def pair_slopes(
    mocks: dict, observed: dict, quants: set[str], min_reads: int,
    rng: random.Random, tol: float,
) -> dict:
    """For every species pair seen at three or more input ratios, fit
    log(observed ratio) against log(input ratio) and ask whether the slope
    is 1."""
    series: dict[tuple, list] = defaultdict(list)
    for sid, rows in mocks.items():
        obs = observed.get(sid)
        if not obs:
            continue
        targets = [r for r in rows if r["seq_id"] not in quants]
        for i in range(len(targets)):
            for j in range(i + 1, len(targets)):
                a, b = targets[i], targets[j]
                na, nb = obs.get(a["seq_id"], 0), obs.get(b["seq_id"], 0)
                if na < min_reads or nb < min_reads:
                    continue
                key = tuple(sorted([a["seq_id"], b["seq_id"]]))
                first = key[0]
                ca, cb = (a, b) if a["seq_id"] == first else (b, a)
                oa = obs.get(ca["seq_id"], 0)
                ob = obs.get(cb["seq_id"], 0)
                series[key].append((
                    math.log10(ca["amount"] / cb["amount"]),
                    math.log10(oa / ob),
                    sid,
                ))

    out = {}
    for key, pts in sorted(series.items()):
        xs = {round(p[0], 4) for p in pts}
        if len(pts) < 3 or len(xs) < 3:
            out[f"{key[0]}|{key[1]}"] = {
                "n_points": len(pts), "n_ratios": len(xs),
                "verdict": "insufficient_data",
            }
            continue
        slope, intercept = ols([(p[0], p[1]) for p in pts])
        idx = list(range(len(pts)))
        lo, hi = boot_ci(
            idx,
            lambda s: ols([(pts[i][0], pts[i][1]) for i in s])[0],
            rng, n=1000,
        )
        excludes_one = (lo == lo and hi == hi) and not (lo <= 1.0 <= hi)
        far = abs(slope - 1.0) > tol
        out[f"{key[0]}|{key[1]}"] = {
            "n_points": len(pts),
            "n_ratios": len(xs),
            "slope": round(slope, 4),
            "slope_ci": [None if lo != lo else round(lo, 4),
                         None if hi != hi else round(hi, 4)],
            "log10_efficiency_ratio": round(intercept, 4),
            "verdict": ("composition_dependent"
                        if (excludes_one and far) else "constant"),
        }
    return out


# --------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mocks", required=True)
    ap.add_argument("--assign-dir", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--outdir", default="calibration")
    ap.add_argument("--min-reads", type=int, default=50,
                    help="reads below which a species is not used in a well")
    ap.add_argument("--slope-tolerance", type=float, default=0.15,
                    help="how far the slope may sit from 1 and still be "
                         "treated as a constant model")
    ap.add_argument("--update-reference",
                    help="write a NEW panel file with corr_factor filled in "
                         "(never edits the input in place)")
    ap.add_argument("--seed", type=int, default=20260920)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    refs = read_table(args.reference)
    quants = {r["seq_id"] for r in refs if r.get("type") == "quant"}
    mocks = load_mocks(args.mocks)
    observed = load_observed(args.assign_dir)
    if not mocks:
        raise SystemExit(f"no mock definitions in {args.mocks}")
    usable = [s for s in mocks if s in observed]
    if not usable:
        raise SystemExit(
            "no mock sample has counts in the assignment directory - check "
            "the sample ids match the samplesheet"
        )

    # per-well efficiencies, kept with their batch so the interval can
    # reflect between-PCR variation rather than only pipetting
    per_species: dict[str, list[float]] = defaultdict(list)
    batches: dict[str, set] = defaultdict(set)
    recovery = []
    for sid in sorted(usable):
        rows = mocks[sid]
        eff = per_well_efficiency(rows, observed[sid], quants, args.min_reads)
        tot_in = sum(r["amount"] for r in rows if r["seq_id"] not in quants)
        tot_obs = sum(observed[sid].get(r["seq_id"], 0)
                      for r in rows if r["seq_id"] not in quants)
        for r in rows:
            n = observed[sid].get(r["seq_id"], 0)
            recovery.append({
                "sample_id": sid, "mix": r["mix"],
                "replicate": r["replicate"], "pcr_batch": r["pcr_batch"],
                "seq_id": r["seq_id"],
                "input_amount": r["amount"],
                "input_fraction": (f"{r['amount']/tot_in:.6f}"
                                   if tot_in and r["seq_id"] not in quants
                                   else ""),
                "observed_reads": n,
                "observed_fraction": (f"{n/tot_obs:.6f}"
                                      if tot_obs and r["seq_id"] not in quants
                                      else ""),
                "efficiency": (f"{eff[r['seq_id']]:.4f}"
                               if r["seq_id"] in eff else ""),
            })
        for k, v in eff.items():
            per_species[k].append(v)
            batches[k].add(rows[0]["pcr_batch"])

    slopes = pair_slopes(mocks, observed, quants, args.min_reads, rng,
                         args.slope_tolerance)

    # A species inherits the verdict of the pairs it took part in. If any
    # pair it belongs to shows composition dependence, a constant for that
    # species is not defensible.
    involved = defaultdict(set)
    for k, v in slopes.items():
        a, b = k.split("|")
        involved[a].add(v.get("verdict"))
        involved[b].add(v.get("verdict"))

    factors = []
    for seq_id in sorted(per_species):
        vals = per_species[seq_id]
        gm = geomean(vals)
        lo, hi = boot_ci(vals, geomean, rng)
        verdicts = involved.get(seq_id, set())
        if "composition_dependent" in verdicts:
            model = "composition_dependent"
        elif len(vals) < 2 or "constant" not in verdicts:
            model = "insufficient_data"
        else:
            model = "constant"
        factors.append({
            "seq_id": seq_id,
            "corr_factor": (f"{gm:.4f}" if model == "constant" else ""),
            "ci_low": ("" if lo != lo else f"{lo:.4f}"),
            "ci_high": ("" if hi != hi else f"{hi:.4f}"),
            "n_wells": len(vals),
            "n_pcr_batches": len(batches.get(seq_id, set())),
            "model": model,
            "note": MODEL_MEANING[model],
        })

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "corr_factors.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(factors[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(factors)
    with open(os.path.join(args.outdir, "mock_recovery.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recovery[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(recovery)

    n_const = sum(1 for f in factors if f["model"] == "constant")
    summary = {
        "n_mock_wells": len(usable),
        "n_pcr_batches": len({r["pcr_batch"] for rows in mocks.values()
                              for r in rows}),
        "species_calibrated": n_const,
        "species_composition_dependent": sum(
            1 for f in factors if f["model"] == "composition_dependent"),
        "species_insufficient_data": sum(
            1 for f in factors if f["model"] == "insufficient_data"),
        "pair_slopes": slopes,
        "model_meaning": MODEL_MEANING,
        "measures": (
            "Amplification bias only. Templates were mixed before PCR, so "
            "extraction and lysis bias are not included; applied to eggs "
            "in stool this corrects one step of three."
        ),
    }
    with open(os.path.join(args.outdir, "calibration.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    if args.update_reference:
        got = {f["seq_id"]: f["corr_factor"] for f in factors
               if f["corr_factor"]}
        for r in refs:
            if r["seq_id"] in got:
                r["corr_factor"] = got[r["seq_id"]]
        with open(args.update_reference, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(refs[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(refs)

    if not args.quiet:
        print(f"{len(usable)} mock wells, "
              f"{summary['n_pcr_batches']} PCR batch(es)")
        for f in factors:
            cf = f["corr_factor"] or "—"
            ci = (f" [{f['ci_low']}–{f['ci_high']}]"
                  if f["ci_low"] else "")
            print(f"    {f['seq_id']:22s} {cf:>8s}{ci:<20s} {f['model']}")
        cd = [k for k, v in slopes.items()
              if v.get("verdict") == "composition_dependent"]
        if cd:
            print("\n    composition-dependent pairs "
                  "(a single corr_factor is not valid):")
            for k in cd:
                v = slopes[k]
                print(f"       {k}  slope {v['slope']:.3f} "
                      f"CI {v['slope_ci'][0]}–{v['slope_ci'][1]}")
            print("    A constant fitted here would overstate every "
                  "minority species.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
