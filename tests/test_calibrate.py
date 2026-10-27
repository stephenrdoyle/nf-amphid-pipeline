#!/usr/bin/env python3
"""
Tests for mock-community calibration.

The failure this suite guards against is a fitter that always returns a
number. A correction factor is easy to compute and hard to justify: fit a
line through any set of mocks and an answer appears, whether or not a
single number is the right shape of answer at all. So the calibrator is
tested twice on the same design - once against data where a constant is
correct, and once against data where it is not - and it has to get both
right. Recovering the right efficiency proves it can measure; refusing to
emit one under competition proves it knows when not to.

Counts are constructed here rather than sequenced. The read-level path is
tested in the other suites; what matters here is the inference, so the
tests build count tables from a known composition and a known efficiency
model and ask whether the calibrator recovers what went in.

Run: python3 tests/test_calibrate.py
"""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import random
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


K = _load("calibrate_mocks", "bin/calibrate_mocks.py")
pf = _load("preflight", "bin/preflight.py")
sim = _load("simulate", "bin/simulate.py")

ASSAY = os.path.join(ROOT, "assets", "assay", "ampseq_helminth.json")
PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


# --------------------------------------------------------------------------
# the small statistics
# --------------------------------------------------------------------------


def test_stats() -> None:
    print("\nstatistics")
    check("geometric mean of a symmetric ratio pair is one",
          abs(K.geomean([0.5, 2.0]) - 1.0) < 1e-12)
    check("geometric mean ignores non-positive values",
          abs(K.geomean([0.0, 4.0, 9.0]) - 6.0) < 1e-9,
          str(K.geomean([0.0, 4.0, 9.0])))
    check("an empty sample has no mean", K.geomean([]) != K.geomean([]))

    slope, icept = K.ols([(0.0, 1.0), (1.0, 3.0), (2.0, 5.0)])
    check("least squares recovers an exact line",
          abs(slope - 2.0) < 1e-9 and abs(icept - 1.0) < 1e-9,
          f"{slope}, {icept}")
    check("a vertical set has no slope",
          K.ols([(1.0, 1.0), (1.0, 2.0)])[0] != K.ols([(1.0, 1.0), (1.0, 2.0)])[0])

    rng = random.Random(3)
    lo, hi = K.boot_ci([1.0] * 20, K.geomean, rng, n=200)
    check("a sample with no spread gives a zero-width interval",
          abs(lo - 1.0) < 1e-9 and abs(hi - 1.0) < 1e-9, f"{lo}, {hi}")
    lo, hi = K.boot_ci([0.5, 1.0, 2.0, 1.2, 0.8], K.geomean, rng, n=500)
    check("a spread sample gives an interval containing its own mean",
          lo < K.geomean([0.5, 1.0, 2.0, 1.2, 0.8]) < hi, f"{lo}, {hi}")
    check("one observation supports no interval",
          K.boot_ci([1.0], K.geomean, rng)[0] != K.boot_ci([1.0], K.geomean, rng)[0])


# --------------------------------------------------------------------------
# a mock run, built to order
# --------------------------------------------------------------------------

RATIOS = [(99, 1), (90, 10), (75, 25), (50, 50), (25, 75), (10, 90), (1, 99)]
TRUTH = {"SppA": 1.80, "SppB": 0.55}


def build_run(tmp: str, competition: float, depth: int = 40000,
              ratios=RATIOS, batches=("B1", "B2", "B3"), seed: int = 11):
    """Mocks at known ratios, counts drawn through the real bias model."""
    rng = random.Random(seed)

    ref = os.path.join(tmp, "reference.tsv")
    with open(ref, "w") as fh:
        fh.write("seq_id\tlocus\ttype\tpurpose\ttaxon\tparent_id\t"
                 "default_conc\tconc_units\tcorr_factor\tsequence\n")
        for s, cf in TRUTH.items():
            fh.write(f"{s}\tITS2\ttarget\tspeciation\t{s}\t\t\t\t{cf}\tACGT\n")
        fh.write("QUANT_hi\tITS2\tquant\t\t\t\t10000\tcopies_per_rxn\t1.0\tACGT\n")
        fh.write("QUANT_lo\tITS2\tquant\t\t\t\t100\tcopies_per_rxn\t1.0\tACGT\n")

    class R:
        def __init__(self, cf):
            self.corr_factor = cf

    by_id = {s: R(cf) for s, cf in TRUTH.items()}
    by_id["QUANT_hi"] = R(1.0)
    by_id["QUANT_lo"] = R(1.0)

    mock_rows, adir = [], os.path.join(tmp, "assign")
    os.makedirs(adir, exist_ok=True)
    for bi, batch in enumerate(batches):
        for a, b in ratios:
            sid = f"MOCK_{a:02d}_{b:02d}_{batch}"
            counts = {"SppA": a * 2000.0, "SppB": b * 2000.0,
                      "QUANT_hi": 10000.0, "QUANT_lo": 100.0}
            for seq, amt in counts.items():
                mock_rows.append({
                    "sample_id": sid, "seq_id": seq, "amount": amt,
                    "amount_units": "copies", "mix": f"Pair_{a}_{b}",
                    "replicate": bi + 1, "pcr_batch": batch, "note": "",
                })
            comp = sim.Composition(counts=dict(counts))
            comp.amplify(by_id, competition)
            total = sum(comp.weights.values())
            obs = {}
            for seq, w in comp.weights.items():
                # multinomial-ish draw with a little well-to-well noise
                mu = depth * w / total
                obs[seq] = max(0, int(rng.gauss(mu, math.sqrt(max(mu, 1)))))
            with open(os.path.join(adir, f"{sid}.counts.tsv"), "w") as fh:
                fh.write("sample_id\tseq_id\tlocus\ttype\tpurpose\ttaxon\t"
                         "reads\tunique_seqs\n")
                for seq, n in obs.items():
                    fh.write(f"{sid}\t{seq}\tITS2\t-\t-\t{seq}\t{n}\t1\n")

    mocks = os.path.join(tmp, "mocks.tsv")
    with open(mocks, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mock_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(mock_rows)
    return ref, mocks, adir


def run_cal(tmp: str, ref, mocks, adir, extra=None):
    out = os.path.join(tmp, "cal")
    with contextlib.redirect_stdout(io.StringIO()):
        K.main(["--mocks", mocks, "--assign-dir", adir, "--reference", ref,
                "--outdir", out, "--quiet"] + (extra or []))
    factors = {r["seq_id"]: r for r in csv.DictReader(
        open(os.path.join(out, "corr_factors.tsv")), delimiter="\t")}
    summary = json.load(open(os.path.join(out, "calibration.json")))
    return factors, summary, out


def test_constant_case() -> None:
    print("\nwhen a single factor is the right answer")
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, adir = build_run(tmp, competition=0.0)
        factors, summary, out = run_cal(tmp, ref, mocks, adir)

        for s, true in TRUTH.items():
            f = factors[s]
            check(f"{s}: model is constant", f["model"] == "constant",
                  f["model"])
            got = float(f["corr_factor"]) if f["corr_factor"] else float("nan")
            check(f"{s}: recovers the true efficiency of {true}",
                  abs(math.log(got / true)) < math.log(1.15),
                  f"got {got:.3f}")
            check(f"{s}: the interval covers the truth",
                  float(f["ci_low"]) <= true <= float(f["ci_high"]),
                  f"[{f['ci_low']}, {f['ci_high']}]")

        pair = summary["pair_slopes"]["SppA|SppB"]
        check("the slope is one, as a multiplicative model implies",
              abs(pair["slope"] - 1.0) < 0.1, f"{pair['slope']}")
        check("and the pair is judged constant",
              pair["verdict"] == "constant", pair["verdict"])
        check("every species is calibrated",
              summary["species_calibrated"] == len(TRUTH),
              str(summary["species_calibrated"]))
        check("what the number does not cover is stated",
              "extraction" in summary["measures"])

        # recovery table, for plotting and for checking by hand
        rec = list(csv.DictReader(
            open(os.path.join(out, "mock_recovery.tsv")), delimiter="\t"))
        check("a per-well recovery table is written",
              len(rec) == 21 * 4, str(len(rec)))
        check("it carries both input and observed fractions",
              all(r["input_fraction"] for r in rec if r["seq_id"] in TRUTH))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_composition_dependent_case() -> None:
    """The one that matters: a constant must be refused, not fitted."""
    print("\nwhen a single factor is the wrong answer")
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, adir = build_run(tmp, competition=0.5)
        factors, summary, out = run_cal(tmp, ref, mocks, adir)

        pair = summary["pair_slopes"]["SppA|SppB"]
        check("the slope is well below one", pair["slope"] < 0.8,
              f"{pair['slope']}")
        lo, hi = pair["slope_ci"]
        check("and its interval excludes one", hi is not None and hi < 1.0,
              f"[{lo}, {hi}]")
        # the slope should recover the strength of the competition injected
        check("the slope measures how strong the competition was",
              abs(pair["slope"] - 0.5) < 0.12, f"{pair['slope']} vs 0.5")

        for s in TRUTH:
            f = factors[s]
            check(f"{s}: reported as composition dependent",
                  f["model"] == "composition_dependent", f["model"])
            check(f"{s}: NO corr_factor is emitted", f["corr_factor"] == "",
                  f"emitted {f['corr_factor']}")
        check("none are counted as calibrated",
              summary["species_calibrated"] == 0,
              str(summary["species_calibrated"]))

        # And nothing constant survives into the updated panel. The input
        # panel carries a factor for both species; this run showed no
        # constant is valid for either, so carrying the old one over -
        # which is what happened - wrote a contradicted number into a
        # file called "calibrated", and synthesis applied it.
        newref = os.path.join(tmp, "updated.tsv")
        _f, summary2, _o = run_cal(tmp, ref, mocks, adir,
                                   ["--update-reference", newref])
        rows = {r["seq_id"]: r for r in csv.DictReader(
            open(newref), delimiter="\t")}
        check("an updated panel clears a composition-dependent factor",
              all(rows[s]["corr_factor"] == "" for s in TRUTH),
              str({s: rows[s]["corr_factor"] for s in TRUTH}))
        cleared = summary2.get("reference_update", {}).get(
            "cleared_composition_dependent", [])
        check("and records what it cleared, and from what",
              sorted(c["seq_id"] for c in cleared) == sorted(TRUTH)
              and all(c["was"] == str(TRUTH[c["seq_id"]]) for c in cleared),
              str(cleared))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_thin_designs() -> None:
    print("\ndesigns that cannot answer the question")
    tmp = tempfile.mkdtemp()
    try:
        # one ratio only: a slope cannot be fitted
        ref, mocks, adir = build_run(tmp, 0.0, ratios=[(50, 50)],
                                     batches=("B1", "B2"))
        factors, summary, _ = run_cal(tmp, ref, mocks, adir)
        pair = summary["pair_slopes"]["SppA|SppB"]
        check("a single ratio yields no slope",
              pair["verdict"] == "insufficient_data", str(pair))
        check("and no corr_factor is emitted from it",
              all(factors[s]["corr_factor"] == "" for s in TRUTH),
              str({s: factors[s]["corr_factor"] for s in TRUTH}))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # and the design checks that run before any sequencing
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, _ = build_run(tmp, 0.0, ratios=[(50, 50)],
                                  batches=("B1",))
        # strip the anchors to trigger the warning
        rows = [r for r in csv.DictReader(open(mocks), delimiter="\t")
                if not r["seq_id"].startswith("QUANT")]
        thin = os.path.join(tmp, "thin.tsv")
        with open(thin, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            pf.main(["--reference", ref, "--mocks", thin, "--assay", ASSAY,
                     "--outdir", os.path.join(tmp, "pf")])
        log = out.getvalue()
        for want, label in (
            ("At least three are needed", "one ratio is flagged"),
            ("no quant standard", "a missing anchor is flagged"),
            ("one PCR batch", "a single PCR batch is flagged"),
            ("no replicate", "an unreplicated mixture is flagged"),
        ):
            check(f"pre-flight: {label}", want in log,
                  "not raised before sequencing")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_update_reference() -> None:
    print("\nwriting the factors back")
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, adir = build_run(tmp, competition=0.0)
        before = open(ref).read()
        newref = os.path.join(tmp, "updated.tsv")
        factors, _s, _o = run_cal(tmp, ref, mocks, adir,
                                  ["--update-reference", newref])
        check("the input panel is not edited in place",
              open(ref).read() == before)
        rows = {r["seq_id"]: r for r in csv.DictReader(
            open(newref), delimiter="\t")}
        check("the new panel carries the fitted factors",
              all(abs(float(rows[s]["corr_factor"])
                      - float(factors[s]["corr_factor"])) < 1e-9
                  for s in TRUTH),
              str({s: rows[s]["corr_factor"] for s in TRUTH}))
        check("and keeps every original column",
              set(rows["SppA"]) == set(csv.DictReader(
                  open(ref), delimiter="\t").fieldnames))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_refusals() -> None:
    """Inputs that used to produce a crash or a wrong number."""
    print("\ninputs that cannot be calibrated")

    # Spike-ins in copies, targets in ng. The efficiency would carry a
    # unit conversion of thousands-fold and look like any other number.
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, adir = build_run(tmp, competition=0.0)
        rows = list(csv.DictReader(open(mocks), delimiter="\t"))
        for r in rows:
            if not r["seq_id"].startswith("QUANT"):
                r["amount_units"] = "ng"
        with open(mocks, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        try:
            run_cal(tmp, ref, mocks, adir)
            got = "ran"
        except SystemExit as e:
            got = str(e)
        check("mixed amount_units in a mock are refused",
              "amount_units" in got and "ng" in got and "copies" in got,
              got[:200])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # Nothing reaches --min-reads: was an IndexError on factors[0].
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, adir = build_run(tmp, competition=0.0, depth=200)
        try:
            run_cal(tmp, ref, mocks, adir, ["--min-reads", "100000"])
            got = "ran"
        except SystemExit as e:
            got = str(e)
        except Exception as e:                      # noqa: BLE001
            got = f"raised {type(e).__name__}: {e}"
        check("no measurable efficiency is a clear refusal, not a crash",
              got.startswith("no efficiency could be measured"), got[:200])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # Too little design to decide: the panel's existing factor is kept -
    # this run is no evidence against it - and listed as unverified.
    tmp = tempfile.mkdtemp()
    try:
        ref, mocks, adir = build_run(tmp, 0.0, ratios=[(50, 50)],
                                     batches=("B1", "B2"))
        newref = os.path.join(tmp, "updated.tsv")
        _f, summary, _o = run_cal(tmp, ref, mocks, adir,
                                  ["--update-reference", newref])
        rows = {r["seq_id"]: r for r in csv.DictReader(
            open(newref), delimiter="\t")}
        check("an insufficient-data species keeps its earlier factor",
              all(rows[s]["corr_factor"] == str(TRUTH[s]) for s in TRUTH),
              str({s: rows[s]["corr_factor"] for s in TRUTH}))
        kept = summary.get("reference_update", {}).get("kept_unverified", [])
        check("and is listed as kept unverified",
              sorted(k["seq_id"] for k in kept) == sorted(TRUTH), str(kept))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_stats()
    test_constant_case()
    test_composition_dependent_case()
    test_thin_designs()
    test_update_reference()
    test_refusals()

    print(f"\n{'='*62}")
    print(f"{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("\nfailures:")
        for name, detail in FAIL:
            print(f"  {name}: {detail}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
