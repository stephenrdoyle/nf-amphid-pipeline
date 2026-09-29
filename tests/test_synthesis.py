#!/usr/bin/env python3
"""
Tests for run-level synthesis.

Two kinds of claim are made here and they need different evidence.

The quantification claims are checked against the simulator's composition
truth - the input copies, recorded before amplification and before
sequencing. Nothing downstream can see that file, so a recovery figure
measured against it is a measurement rather than a consistency check. The
amplification correction is tested by running the same data twice, once
with corr_factor and once without: a correction that cannot be shown to
change anything has not been shown to work.

The plate-geometry claim is statistical, so the test is too. A permutation
test is only worth reporting if its p-values mean what they say, so the
null is simulated directly - values with no relationship to well position -
and the false-positive rate is checked against nominal. Power is measured
as well, because a sensor that never fires is not evidence of a clean
plate.

Run: python3 tests/test_synthesis.py
"""

from __future__ import annotations

import contextlib
import csv
import glob
import importlib.util
import io
import json
import math
import os
import random
import shutil
import sys
import tempfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S = _load("synthesise", "bin/synthesise.py")
A = _load("assign", "bin/assign.py")
P = _load("prep_illumina", "bin/prep_illumina.py")
pf = _load("preflight", "bin/preflight.py")
sim = _load("simulate", "bin/simulate.py")
sys.path.insert(0, HERE)
import make_fixtures as mf  # noqa: E402

INDEX_MAP = os.path.join(ROOT, "assets", "index_map.tsv")
ASSAY = os.path.join(ROOT, "assets", "assay", "ampseq_helminth.json")

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else float("nan")


# --------------------------------------------------------------------------
# calibration
# --------------------------------------------------------------------------


QUANTS = [
    {"seq_id": "QUANT_hi", "default_conc": "10000"},
    {"seq_id": "QUANT_lo", "default_conc": "100"},
]


def sample_with(counts, depth=5000, stype="unknown"):
    s = S.Sample(sample_id="X", sample_type=stype, depth=depth)
    s.counts = dict(counts)
    return s


def test_calibration() -> None:
    print("\ntwo-point calibration")

    # a linear well: 100x more input gave 100x more reads
    s = sample_with({"QUANT_hi": 10000, "QUANT_lo": 100})
    c = S.calibrate(s, QUANTS, 10, 0.25)
    check("a linear well fits a slope of one",
          c.state == "spiked" and abs(c.slope - 1.0) < 1e-9, f"{c.slope}")
    check("and inverts back to the input it was given",
          abs(c.copies(10000) - 10000) < 1, str(c.copies(10000)))
    check("a target read count converts to copies",
          abs(c.copies(500) - 500) < 1, str(c.copies(500)))

    # saturation: the high standard did not keep up
    s = sample_with({"QUANT_hi": 2000, "QUANT_lo": 100})
    c = S.calibrate(s, QUANTS, 10, 0.25)
    check("a compressed range gives a slope below one",
          c.state == "spiked" and c.slope < 0.8, f"{c.slope:.3f}")
    check("and says so, rather than reporting the number quietly",
          "compressed" in c.note, c.note)

    # the well failed: no high standard
    # An empty well means different things depending on the run it sits in,
    # and only the run can say which.
    s = sample_with({"QUANT_hi": 0, "QUANT_lo": 0})
    c = S.calibrate(s, QUANTS, 10, 0.25, run_spiked=True)
    check("no standards, on a spiked run, is a failed well",
          c.state == "spike_failed", f"{c.state}")
    c = S.calibrate(s, QUANTS, 10, 0.25, run_spiked=False)
    check("the same well on an unspiked run is merely unverified",
          c.state == "unspiked", f"{c.state}")
    check("an uncalibrated well yields no copy number",
          c.copies(500) is None)

    s = sample_with({"QUANT_hi": 1, "QUANT_lo": 0})
    c = S.calibrate(s, QUANTS, 10, 0.25, run_spiked=True)
    check("a trace of the high standard is still a failed well",
          c.state == "spike_failed", f"{c.state}")

    # and the run-level inference itself
    spiked = {f"s{i}": sample_with({"QUANT_hi": 900, "QUANT_lo": 9})
              for i in range(10)}
    spiked["dead"] = sample_with({"QUANT_hi": 0, "QUANT_lo": 0})
    check("a run where most wells recovered standards counts as spiked",
          S.run_was_spiked(spiked, QUANTS, 10))
    none = {f"s{i}": sample_with({"QUANT_hi": 0, "QUANT_lo": 0})
            for i in range(10)}
    check("a run where none did counts as unspiked",
          not S.run_was_spiked(none, QUANTS, 10))

    # only the low standard recovered - shallow, but usable
    s = sample_with({"QUANT_hi": 9000, "QUANT_lo": 2})
    c = S.calibrate(s, QUANTS, 10, 0.25)
    check("losing the low standard costs the bottom of the range, not the well",
          c.state == "spiked" and "unverified" in c.note, c.note)


def test_calls() -> None:
    print("\nthree-state calling")
    tmp = tempfile.mkdtemp()
    try:
        ref = os.path.join(tmp, "reference.tsv")
        with open(ref, "w") as fh:
            fh.write("seq_id\tlocus\ttype\tpurpose\ttaxon\tparent_id\t"
                     "default_conc\tconc_units\tcorr_factor\tsequence\n")
            fh.write("T1\tITS2\ttarget\tspeciation\tSpecies_a\t\t\t\t\tACGT\n")
            fh.write("QUANT_hi\tITS2\tquant\t\t\t\t10000\tcopies_per_rxn\t1.0\tACGT\n")
            fh.write("QUANT_lo\tITS2\tquant\t\t\t\t100\tcopies_per_rxn\t1.0\tACGT\n")
        ss = os.path.join(tmp, "samplesheet.csv")
        with open(ss, "w") as fh:
            fh.write("run_id,platform,sample_id,index_i7,index_i5,sample_type\n")
            for sid, st in (("BLANK-01", "neg_ntc"), ("GOOD", "unknown"),
                            ("FAILED", "unknown"), ("NOSPIKE", "unknown")):
                fh.write(f"RUN01,illumina,{sid},A,B,{st}\n")

        adir = os.path.join(tmp, "assign")
        os.makedirs(adir)
        rows = {
            "BLANK-01": {"T1": 0, "QUANT_hi": 900, "QUANT_lo": 9},
            "GOOD":     {"T1": 400, "QUANT_hi": 1000, "QUANT_lo": 10},
            # standards absent: the PCR did not work in this well
            "FAILED":   {"T1": 0, "QUANT_hi": 0, "QUANT_lo": 0},
            "NOSPIKE":  {"T1": 0, "QUANT_hi": 0, "QUANT_lo": 0},
        }
        for sid, counts in rows.items():
            with open(os.path.join(adir, f"{sid}.counts.tsv"), "w") as fh:
                fh.write("sample_id\tseq_id\tlocus\ttype\tpurpose\ttaxon\t"
                         "reads\tunique_seqs\n")
                for k, v in counts.items():
                    fh.write(f"{sid}\t{k}\tITS2\t-\t-\t-\t{v}\t1\n")
            with open(os.path.join(adir, f"{sid}.assign.json"), "w") as fh:
                json.dump({"sample_id": sid,
                           "sequences_in": sum(counts.values())}, fh)

        out = os.path.join(tmp, "synth")
        with contextlib.redirect_stdout(io.StringIO()):
            S.main(["--assign-dir", adir, "--samplesheet", ss,
                    "--reference", ref, "--outdir", out,
                    "--run-id", "RUN01", "--permutations", "20", "--quiet"])
        calls = {r["sample_id"]: r for r in csv.DictReader(
            open(os.path.join(out, "calls.tsv")), delimiter="\t")}

        check("a well with signal above the limit is detected",
              calls["GOOD"]["call"] == "detected", str(calls["GOOD"]))
        check("a failed well is assay_failed, never not_detected",
              calls["FAILED"]["call"] == "assay_failed",
              f"{calls['FAILED']['call']} - "
              f"{calls['FAILED']['internal_control']}")
        check("the reason for the failure travels with the call",
              "did not amplify" in calls["FAILED"]["basis"]
              or "below the minimum" in calls["FAILED"]["basis"],
              calls["FAILED"]["basis"])
        check("a detected well reports copies",
              calls["GOOD"]["copies_per_unit_dna"] != "")
        check("a failed well reports none",
              calls["FAILED"]["copies_per_unit_dna"] == "")
        check("the blank is not called detected",
              calls["BLANK-01"]["call"] != "detected")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_detection_limit() -> None:
    print("\ndetection limit from the run's own blanks")
    samples = {
        "B1": sample_with({"T1": 5}, depth=1000, stype="neg_ntc"),
        "B2": sample_with({"T1": 2}, depth=1000, stype="neg_extraction"),
        "S1": sample_with({"T1": 4}, depth=1000),
        "S2": sample_with({"T1": 4}, depth=20000),
    }
    for sid, s in samples.items():
        s.sample_id = sid
    lods = S.detection_limits(samples, ["T1"], 2)
    check("the limit is taken from the worst blank, not the average",
          lods["T1"]["max_reads_in_blank"] == 5, str(lods["T1"]))

    lim1, basis1 = S.limit_for("T1", samples["S1"], lods)
    lim2, basis2 = S.limit_for("T1", samples["S2"], lods)
    check("a shallow well is held to the absolute floor",
          lim1 == 5, f"{lim1} ({basis1})")
    check("a deep well is held to a proportionally higher one",
          lim2 > lim1, f"{lim2} vs {lim1}")
    check("and the output says which bound applied",
          "proportional" in basis2, basis2)

    clean = {"B1": sample_with({"T1": 0}, 1000, "neg_ntc")}
    clean["B1"].sample_id = "B1"
    lods2 = S.detection_limits(clean, ["T1"], 2)
    lim, basis = S.limit_for("T1", sample_with({"T1": 1}, 1000), lods2)
    check("clean blanks still leave a floor - one read is never a detection",
          lim == 2, f"{lim} ({basis})")

    none = S.detection_limits({}, ["T1"], 2)
    lim, basis = S.limit_for("T2", sample_with({"T2": 1}), none)
    check("with no blanks at all the limit is declared an assumption",
          "assumption" in basis, basis)


# --------------------------------------------------------------------------
# the statistics
# --------------------------------------------------------------------------


def test_permutation_calibration() -> None:
    """A p-value that is not calibrated is decoration."""
    print("\nplate-effect test: calibration and power")
    rng = random.Random(1)
    ps = []
    for trial in range(240):
        vals = [rng.lognormvariate(0, 1.5) for _ in range(88)]
        labs = [chr(65 + i % 8) for i in range(88)]
        rng.shuffle(labs)
        ps.append(S.permutation_p(vals, labs, random.Random(trial), 200)[1])
    r05 = sum(1 for p in ps if p < 0.05) / len(ps)
    r10 = sum(1 for p in ps if p < 0.10) / len(ps)
    check("under the null, 5% of runs fall below p=0.05",
          0.02 <= r05 <= 0.09, f"{r05:.3f}")
    check("and 10% below p=0.10", 0.05 <= r10 <= 0.16, f"{r10:.3f}")
    check("a permutation test never claims p = 0", min(ps) > 0, str(min(ps)))

    hits = 0
    trials = 60
    for trial in range(trials):
        r = random.Random(500 + trial)
        vals, labs = [], []
        for i in range(88):
            row = chr(65 + i % 8)
            vals.append(r.lognormvariate(0, 0.8) * (2.0 if row in "AB" else 1.0))
            labs.append(row)
        hits += S.permutation_p(vals, labs, random.Random(trial), 200)[1] < 0.05
    power = hits / trials
    check("a twofold effect in two of eight rows is usually caught",
          power > 0.4, f"{power:.2f}")


# --------------------------------------------------------------------------
# end to end against composition truth
# --------------------------------------------------------------------------


def build_run(tmp: str, n_samples=30, reads=1500, corr=True):
    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    if corr:
        # amplification efficiencies that differ, so the correction has
        # something to correct
        factors = {"Tri": 0.55, "Asc": 1.80}
        for r in ref:
            if r["type"] == "target":
                r["corr_factor"] = str(
                    factors.get(r["seq_id"].split("_")[0], 1.0))
            elif r["type"] == "quant":
                r["corr_factor"] = "1.0"
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, n_samples, 2, "RUN01", "illumina")
    mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
    mf.write_primers(os.path.join(tmp, "primers.tsv"))
    mf.write_samplesheet(samples, os.path.join(tmp, "samplesheet.csv"))
    mf.write_runs(os.path.join(tmp, "runs.tsv"), "RUN01", "illumina",
                  os.path.join(tmp, "sim", "RUN01"))

    with contextlib.redirect_stdout(io.StringIO()):
        sim.main(["--reference", os.path.join(tmp, "reference.tsv"),
                  "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                  "--runs", os.path.join(tmp, "runs.tsv"), "--assay", ASSAY,
                  "--outdir", os.path.join(tmp, "sim"),
                  "--reads-per-sample", str(reads)])
        pf.main(["--reference", os.path.join(tmp, "reference.tsv"),
                 "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                 "--runs", os.path.join(tmp, "runs.tsv"),
                 "--primers", os.path.join(tmp, "primers.tsv"),
                 "--index-map", INDEX_MAP, "--assay", ASSAY,
                 "--outdir", os.path.join(tmp, "preflight")])
        for r1 in sorted(glob.glob(os.path.join(tmp, "sim", "RUN01",
                                                "*_R1_001.fastq.gz"))):
            sid = os.path.basename(r1).split("_S")[0]
            if sid == "Undetermined":
                continue
            P.main(["--r1", r1, "--r2", r1.replace("_R1_", "_R2_"),
                    "--sample-id", sid,
                    "--primers", os.path.join(tmp, "primers.tsv"),
                    "--assay", ASSAY, "--outdir", os.path.join(tmp, "prep"),
                    "--quiet"])
            A.main(["--panel", os.path.join(tmp, "preflight",
                                            "reference.insert.fasta"),
                    "--reference", os.path.join(tmp, "reference.tsv"),
                    "--sample-id", sid,
                    "--merged", os.path.join(tmp, "prep",
                                             f"{sid}.merged.fasta.gz"),
                    "--pairs", os.path.join(tmp, "prep",
                                            f"{sid}.pairs.tsv.gz"),
                    "--assay", ASSAY,
                    "--outdir", os.path.join(tmp, "assign"), "--quiet"])
    truth = defaultdict(dict)
    with open(os.path.join(tmp, "sim", "RUN01.composition.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            truth[r["sample_id"]][r["seq_id"]] = float(r["input_copies"])
    return truth


def synth(tmp: str, reference: str, out: str):
    with contextlib.redirect_stdout(io.StringIO()):
        S.main(["--assign-dir", os.path.join(tmp, "assign"),
                "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                "--reference", reference,
                "--preflight-dir", os.path.join(tmp, "preflight"),
                "--outdir", out, "--run-id", "RUN01",
                "--permutations", "200", "--quiet"])
    return list(csv.DictReader(open(os.path.join(out, "calls.tsv")),
                               delimiter="\t"))


def synth_with(tmp: str, reference: str, out: str, extra):
    """As synth(), plus extra flags. Kept separate so the existing
    callers keep their exact argument list."""
    with contextlib.redirect_stdout(io.StringIO()):
        S.main(["--assign-dir", os.path.join(tmp, "assign"),
                "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                "--reference", reference,
                "--preflight-dir", os.path.join(tmp, "preflight"),
                "--outdir", out, "--run-id", "RUN01",
                "--permutations", "200", "--quiet"] + list(extra))
    with open(os.path.join(out, "calls.tsv"), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(out, "run_qc.json"), encoding="utf-8") as fh:
        qc = json.load(fh)
    return rows, qc


def test_end_to_end() -> None:
    print("\nquantification against the composition truth")
    tmp = tempfile.mkdtemp()
    try:
        truth = build_run(tmp)
        ref = os.path.join(tmp, "reference.tsv")
        calls = synth(tmp, ref, os.path.join(tmp, "synth"))

        # accounting: one row per sample per target, no more, no fewer
        keys = [(c["sample_id"], c["seq_id"]) for c in calls]
        check("every sample x target appears exactly once",
              len(keys) == len(set(keys)), f"{len(keys)} rows, "
              f"{len(set(keys))} unique")
        check("every call is one of the three states",
              all(c["call"] in S.CALLS for c in calls))

        errs, by_ref = [], defaultdict(list)
        for c in calls:
            if c["call"] != "detected" or not c["copies_per_unit_dna"]:
                continue
            t = truth.get(c["sample_id"], {}).get(c["seq_id"], 0.0)
            if t <= 0:
                continue
            lr = math.log10(float(c["copies_per_unit_dna"]) / t)
            errs.append(lr)
            by_ref[c["seq_id"]].append(lr)
        check("enough detections to measure recovery", len(errs) > 20,
              str(len(errs)))
        within2 = sum(1 for x in errs if abs(x) < math.log10(2)) / len(errs)
        check("recovered copies land within twofold of the input",
              within2 > 0.85, f"{within2:.1%}")
        check("and are unbiased overall", abs(med(errs)) < 0.10,
              f"median log10 ratio {med(errs):+.3f}")

        # specificity: nothing in a blank or a mock
        bad = [c for c in calls
               if c["sample_type"] in ("neg_ntc", "neg_extraction", "mock")
               and c["call"] == "detected"]
        check("no target is detected in a blank or a mock", not bad,
              "; ".join(f"{c['sample_id']}/{c['seq_id']}" for c in bad[:3]))

        # sensitivity, against what was actually put in
        tp = fn = 0
        for c in calls:
            t = truth.get(c["sample_id"], {}).get(c["seq_id"], 0.0)
            if t <= 0 or c["call"] == "assay_failed":
                continue
            if c["call"] == "detected":
                tp += 1
            else:
                fn += 1
        sens = tp / max(1, tp + fn)
        check("most of what was put in is recovered", sens > 0.75,
              f"{sens:.3f} ({tp} of {tp + fn})")

        # the control experiment for corr_factor
        rows = list(csv.DictReader(open(ref), delimiter="\t"))
        for r in rows:
            r["corr_factor"] = ""
        nocf = os.path.join(tmp, "reference_nocf.tsv")
        with open(nocf, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        un = synth(tmp, nocf, os.path.join(tmp, "synth_nocf"))
        un_by_ref = defaultdict(list)
        for c in un:
            if c["call"] != "detected" or not c["copies_per_unit_dna"]:
                continue
            t = truth.get(c["sample_id"], {}).get(c["seq_id"], 0.0)
            if t > 0:
                un_by_ref[c["seq_id"]].append(
                    math.log10(float(c["copies_per_unit_dna"]) / t))

        check("the output records the factor it used and where it came from",
              {c["corr_factor_source"] for c in calls} == {"panel"}
              and {c["corr_factor_source"] for c in un} == {"default"},
              f"{sorted({c['corr_factor_source'] for c in calls})} vs "
              f"{sorted({c['corr_factor_source'] for c in un})}")
        check("an unstated corr_factor defaults to 1, not to nothing",
              {c["corr_factor"] for c in un} == {"1"},
              str(sorted({c["corr_factor"] for c in un})))

        # With every factor defaulted to 1 the arithmetic is unchanged, so
        # the comparison below is measured vs assumed-unity - which is
        # exactly the choice a panel without measured efficiencies makes.
        skewed = [k for k in by_ref
                  if k in un_by_ref and abs(med(un_by_ref[k])) > 0.15]
        check("a reference that amplifies differently is biased without "
              "the correction", len(skewed) >= 2,
              str({k: round(med(un_by_ref[k]), 3) for k in un_by_ref}))
        fixed = all(abs(med(by_ref[k])) < abs(med(un_by_ref[k])) / 2
                    for k in skewed)
        check("and the correction removes that bias", fixed,
              str({k: (round(med(un_by_ref[k]), 3), round(med(by_ref[k]), 3))
                   for k in skewed}))

        qc = json.load(open(os.path.join(tmp, "synth", "run_qc.json")))
        check("run QC records the depth spread",
              qc["depth"]["fold_spread"] >= 1, str(qc["depth"]))
        check("run QC reports the calibration slopes",
              qc["calibration_slope"]["n"] > 0,
              str(qc["calibration_slope"]))
        check("plate geometry was tested, the layout being known",
              len(qc["plate_geometry"]) > 0)
        check("plate-effect p-values carry their multiple-testing caveat",
              "multiple testing" in qc["plate_geometry_note"])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_read_funnel() -> None:
    """The funnel exists to separate two failures that look alike.

    A well that produced plenty of reads and lost them to primer dimer,
    and a well that was simply loaded light, both show a short bar in the
    depth panel. Only the funnel distinguishes them, and only it compares
    reads leaving preparation with reads entering assignment - a boundary
    no other part of the pipeline checks.
    """
    print("\nthe read funnel")
    from collections import Counter as _C

    def sample(depth):
        s = S.Sample(sample_id="x")
        s.depth = depth
        return s

    healthy = {
        f"S{i}": {"sample_id": f"S{i}", "pairs_in": 5000,
                  "fates": {"merged": 4700, "paired": 100,
                            "no_primer": 200}}
        for i in range(10)
    }
    samples = {k: sample(4800) for k in healthy}
    fates = _C({"target": 40000, "quant": 6000, "unmatched": 2000})
    f = S.read_funnel(samples, healthy, None, fates)

    names = [st["stage"] for st in f["stages"]]
    check("the funnel starts at the sequencer, not at assignment",
          names[0].startswith("read pairs"), str(names))
    check("and ends at reads on a target", names[-1].endswith("a target"),
          str(names))
    check("it is monotonic - no stage gains reads",
          all(a["reads"] >= b["reads"]
              for a, b in zip(f["stages"], f["stages"][1:])),
          str([st["reads"] for st in f["stages"]]))
    prep_stage = next(st for st in f["stages"] if "preparation" in st["stage"])
    check("prep losses are named, not lumped together",
          prep_stage["lost_here"].get("no_primer") == 2000,
          str(prep_stage["lost_here"]))
    check("percentages are of the starting count",
          abs(f["stages"][0]["percent_of_start"] - 100.0) < 1e-9)
    check("a clean run balances across the prep/assignment boundary",
          f["handoff_balances"] is True)
    check("and flags no wells", not f["wells_losing_over_half_at_prep"])

    # one well drowned in dimer, hidden inside a healthy run
    bad = dict(healthy)
    bad["BAD"] = {"sample_id": "BAD", "pairs_in": 5000,
                  "fates": {"merged": 900, "paired": 100,
                            "no_primer": 4000}}
    s2 = {k: sample(sum(v["fates"].get(x, 0) for x in ("merged", "paired")))
          for k, v in bad.items()}
    f2 = S.read_funnel(s2, bad, None, fates)
    flagged = [o["sample_id"] for o in f2["wells_losing_over_half_at_prep"]]
    check("a well that lost most reads to dimer is singled out",
          flagged == ["BAD"], str(flagged))
    check("with the fraction it kept",
          f2["wells_losing_over_half_at_prep"][0]["kept_fraction"] == 0.2,
          str(f2["wells_losing_over_half_at_prep"][0]["kept_fraction"]))
    check("and the note says why a depth panel would not show it",
          any("shallow" in n for n in f2["notes"]), str(f2["notes"]))

    # reads disappearing between two tools
    f3 = S.read_funnel({k: sample(100) for k in healthy}, healthy, None,
                       fates)
    check("reads lost between prep and assignment are caught",
          f3["handoff_balances"] is False)
    check("and the loss is quantified",
          any("left preparation" in n for n in f3["notes"]),
          str(f3["notes"]))

    # --- the nanopore path, which was silently broken ---------------
    # The demux summary keys its tiers under "by_tier". An earlier
    # version looked for "tiers", found nothing, and skipped the
    # demultiplexing stage without complaint - the funnel simply began
    # at assignment and nobody could tell it had lost a step.
    demux = {
        "run_id": "R", "total_reads": 8109, "assigned_reads": 7732,
        "by_tier": {"A": 7732, "B": 343, "C": 32, "D": 2},
        "tier_meaning": {
            "A": "both ends agree - ok", "B": "one end only - truncated",
            "C": "ends disagree - chimeric", "D": "no confident match",
        },
    }
    ont = {f"S{i}": sample(7732 // 12) for i in range(12)}
    # make the totals line up exactly, as a real run would
    ont["S0"] = sample(7732 - (7732 // 12) * 11)
    fo = S.read_funnel(ont, {}, demux, _C({"target": 2978, "quant": 3440,
                                           "control": 880, "unmatched": 17,
                                           "no_primer": 417}))
    names = [st["stage"] for st in fo["stages"]]
    check("the nanopore funnel starts at the whole run",
          names[0] == "reads in the run", str(names))
    check("and includes a demultiplexing stage",
          any("demultiplexed" in n for n in names), str(names))
    dx = next(st for st in fo["stages"] if "demultiplexed" in st["stage"])
    check("tier A is what survives demultiplexing", dx["reads"] == 7732,
          str(dx["reads"]))
    check("the other tiers are named individually, not pooled",
          len(dx["lost_here"]) == 3, str(dx["lost_here"]))
    check("and each carries a plain-language reason, not just a letter",
          all("(" in k and len(k) > len("tier X") + 3
              for k in dx["lost_here"]),
          str(list(dx["lost_here"])))
    check("the nanopore handoff is checked too",
          fo["handoff_balances"] is True and
          fo["handoff_from"] == "demultiplexing",
          f"{fo['handoff_from']}: {fo['handoff_balances']}")

    short = {f"S{i}": sample(10) for i in range(12)}
    fo2 = S.read_funnel(short, {}, demux, fates)
    check("reads lost between demultiplexing and assignment are caught",
          fo2["handoff_balances"] is False)
    check("and the note names demultiplexing as the source",
          any("left demultiplexing" in n for n in fo2["notes"]),
          str(fo2["notes"]))

    # with no prep directory the funnel still runs, just starting later
    f4 = S.read_funnel(samples, {}, None, fates)
    check("without prep data the funnel starts at assignment",
          f4["stages"][0]["stage"] == "reached assignment",
          f4["stages"][0]["stage"])
    check("and claims no boundary check it could not make",
          f4["handoff_balances"] is None)


def test_control_plasmid_reporting() -> None:
    """Plasmids are contamination or subject, and the flag decides only calls.

    A tagged control plasmid means two different things depending on the
    run. In a faecal sample it is contamination, and the number that
    matters is what share of the well it took. In a titration the
    plasmids ARE the experiment and each needs a call of its own.

    What must NOT depend on the flag is arithmetic. Both proportions are
    emitted either way - within a reference's own group, and plasmid as a
    share of all panel amplicons - because a column whose denominator
    changes with a command-line switch is a column that will be misread.
    """
    print("\ncontrol plasmid reporting")
    tmp = tempfile.mkdtemp()
    try:
        build_run(tmp)
        ref = os.path.join(tmp, "reference.tsv")
        off, qoff = synth_with(tmp, ref, os.path.join(tmp, "s_off"), [])
        on, qon = synth_with(tmp, ref, os.path.join(tmp, "s_on"),
                             ["--controls-as-targets"])

        check("by default only targets are called",
              {r["reference_type"] for r in off} == {"target"},
              str({r["reference_type"] for r in off}))
        check("with the flag, controls are called too",
              {r["reference_type"] for r in on} == {"target", "control"},
              str({r["reference_type"] for r in on}))
        check("and every target row is byte-identical either way",
              [r for r in on if r["reference_type"] == "target"] == off,
              "the flag altered target rows, which it must not")

        poff = qoff["plasmid_contamination"]
        pon = qon["plasmid_contamination"]
        check("the contamination measure does not depend on the flag",
              poff["run_plasmid_share"] == pon["run_plasmid_share"],
              f"{poff['run_plasmid_share']} vs {pon['run_plasmid_share']}")
        check("it is reported even when plasmids are the subject",
              pon["run_plasmid_share"] is not None)
        check("and the flag is recorded in the output, not merely applied",
              poff["controls_called_as_targets"] is False
              and pon["controls_called_as_targets"] is True)

        # A denominator error shows up here and almost nowhere else.
        for tag, rows in (("off", off), ("on", on)):
            by: dict = {}
            for r in rows:
                if r["fraction_of_group"]:
                    by.setdefault((r["sample_id"], r["reference_type"]),
                                  []).append(float(r["fraction_of_group"]))
            bad = [k for k, v in by.items() if abs(sum(v) - 1.0) > 1e-3]
            check(f"[{tag}] within-group proportions sum to 1 per sample",
                  not bad, f"{len(bad)} groups off, e.g. {bad[:2]}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_cohort_completeness() -> None:
    """Synthesis must refuse to report on part of a plate.

    This is the failure that motivated the check: a Nextflow fan-out
    collapsed eight samples to one, every stage exited zero, and
    synthesis produced a clean report on 1/8 of the run. Nothing in the
    pipeline objected, because nothing compared the cohort that was
    declared with the cohort that arrived - even though synthesis reads
    both files.
    """
    print("\nsynthesis refuses a partial cohort")
    # Which samples arrive is the subject here, not the
    # statistics, so the permutation test is turned down.
    CHEAP = ["--permutations", "5"]
    tmp = tempfile.mkdtemp()
    try:
        build_run(tmp, n_samples=8, reads=700)
        ref = os.path.join(tmp, "reference.tsv")
        assign = os.path.join(tmp, "assign")
        declared = sorted(glob.glob(os.path.join(assign, "*.counts.tsv")))
        check("the fixture built a full cohort", len(declared) == 8,
              str(len(declared)))

        # Whole plate: no complaint, and the QC says so.
        _, qc = synth_with(tmp, ref, os.path.join(tmp, "full"), CHEAP)
        check("a complete cohort is reported normally", qc["cohort_complete"])
        check("and records nothing absent", qc["samples_absent"] == [],
              str(qc["samples_absent"]))

        # Now lose seven of them, exactly as the fan-out bug did.
        kept = declared[0]
        for path in declared[1:]:
            os.remove(path)
            js = path.replace(".counts.tsv", ".assign.json")
            if os.path.exists(js):
                os.remove(js)

        try:
            synth_with(tmp, ref, os.path.join(tmp, "partial"), CHEAP)
            check("a partial cohort is fatal", False, "no error raised")
        except SystemExit as exc:
            msg = str(exc)
            check("a partial cohort is fatal", True)
            check("and the error counts both sides",
                  "8 samples declared" in msg and "1 reached" in msg, msg[:90])
            check("and names the samples that were lost",
                  "MOCK-A" in msg or "iCCH" in msg, msg[:90])

        # A nanopore well that got no reads is NOT a lost sample. The
        # demultiplexing summary is what tells the two apart, so the
        # check must consult it rather than refusing every gap.
        kept_id = os.path.basename(kept)[:-len(".counts.tsv")]
        demux = os.path.join(tmp, "demux_summary.json")
        with open(demux, "w", encoding="utf-8") as fh:
            json.dump({"per_sample": {kept_id: 1500}}, fh)
        _, qc = synth_with(tmp, ref, os.path.join(tmp, "empty"),
                           CHEAP + ["--demux-summary", demux])
        check("an empty nanopore well is not treated as lost",
              len(qc["samples_absent"]) == 7, str(qc["samples_absent"]))
        check("but the run is still not called complete",
              qc["cohort_complete"] is False)

        # ...whereas a sample the demultiplexer DID write reads for, and
        # which then failed to arrive, is a lost task and still fatal.
        gone = os.path.basename(declared[1])[:-len(".counts.tsv")]
        with open(demux, "w", encoding="utf-8") as fh:
            json.dump({"per_sample": {kept_id: 1500, gone: 9000}}, fh)
        try:
            synth_with(tmp, ref, os.path.join(tmp, "lost"),
                       CHEAP + ["--demux-summary", demux])
            check("reads written but not assigned is fatal", False,
                  "no error raised")
        except SystemExit as exc:
            check("reads written but not assigned is fatal", True)
            check("and only that sample is named as lost",
                  gone in str(exc), str(exc)[:90])

        # The override exists, and says what it cost.
        _, qc = synth_with(tmp, ref, os.path.join(tmp, "forced"),
                           CHEAP + ["--allow-missing"])
        check("--allow-missing reports anyway", qc["n_samples"] == 1)
        check("and the report admits the cohort was incomplete",
              qc["cohort_complete"] is False
              and len(qc["samples_absent"]) == 7,
              str(qc["samples_absent"]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_calibration()
    test_read_funnel()
    test_control_plasmid_reporting()
    test_detection_limit()
    test_calls()
    test_permutation_calibration()
    test_end_to_end()
    test_cohort_completeness()

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
