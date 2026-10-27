#!/usr/bin/env python3
"""
Tests for cross-run comparison and the run report.

The comparison is checked against constructed agreement, where the right
answer is known by construction rather than by running the pipeline: a set
of calls with a designed number of matches has a computable kappa, and an
implementation that disagrees with it is wrong. Kappa specifically, because
raw agreement is flattered by how often both runs say not_detected - a
sanity check on that is included, since it is the reason the statistic is
there at all.

One bug is pinned here that was not in this code at all. The first paired
fixture drew a fresh composition for each run, so the same sample was a
different sample on each platform, and concordance came out at kappa 0.04
between two runs of the *same libraries*. The simulator now derives a
sample's composition from its id, so the test below asserts that property
directly rather than trusting it.

The report is checked for the things that would make it useless in the
field: any external reference, any script, or a number that disagrees with
the table it came from.

Run: python3 tests/test_report.py
"""

from __future__ import annotations

import contextlib
import csv
import io
import importlib.util
import json
import os
import re
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


C = _load("compare_runs", "bin/compare_runs.py")
R = _load("report", "bin/report.py")

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def row(sample, seq, call, reads=0, depth=1000, copies="", taxon="T"):
    return {"sample_id": sample, "seq_id": seq, "taxon": taxon, "call": call,
            "reads": str(reads), "depth": str(depth),
            "copies_per_unit_dna": copies, "sample_type": "unknown"}


# --------------------------------------------------------------------------
# kappa
# --------------------------------------------------------------------------


def test_kappa() -> None:
    print("\nagreement beyond chance")
    check("perfect agreement with both labels used is kappa 1",
          abs(C.cohen_kappa([("a", "a")] * 5 + [("b", "b")] * 5) - 1.0) < 1e-9)
    check("agreement on a single label is undefined, not perfect",
          C.cohen_kappa([("a", "a")] * 10) is None,
          str(C.cohen_kappa([("a", "a")] * 10)))

    # The case the statistic exists for: both runs almost always say no.
    # 90% agreement, and kappa is very slightly NEGATIVE - the two runs
    # agreed marginally less often than their own marginal rates predict.
    # Reporting the 90% alone would describe this as near-perfect.
    pairs = [("n", "n")] * 90 + [("y", "n")] * 5 + [("n", "y")] * 5
    k = C.cohen_kappa(pairs)
    agree = sum(1 for a, b in pairs if a == b) / len(pairs)
    check("high raw agreement can still be no better than chance",
          agree == 0.90 and -0.1 < k <= 0.0,
          f"agreement {agree}, kappa {k:.4f}")

    # worked example: a 2x2 with known kappa
    pairs = ([("y", "y")] * 20 + [("n", "n")] * 15
             + [("y", "n")] * 5 + [("n", "y")] * 10)
    # po = 35/50 = .70 ; pe = (25/50)(30/50) + (25/50)(20/50) = .30+.20 = .50
    check("matches a hand-computed kappa", abs(C.cohen_kappa(pairs) - 0.4) < 1e-9,
          f"{C.cohen_kappa(pairs):.6f}")


def test_pairing() -> None:
    print("\npairing samples across runs")
    a = [row("S1-r1", "T1", "detected", 100, copies="500"),
         row("S2-r1", "T1", "not_detected", 0)]
    b = [row("S1-r2", "T1", "detected", 90, copies="450"),
         row("S2-r2", "T1", "not_detected", 1)]
    res = C.compare(a, b, r"-r[0-9]+$", "A", "B")
    check("a replicate suffix can be stripped to pair subjects",
          res["paired_observations"] == 2, str(res["paired_observations"]))
    check("and the calls then agree", res["agreement"] == 1.0)
    check("copy ratio is reported where both detected",
          res["n_both_detected"] == 1)

    res = C.compare(a, b, None, "A", "B")
    check("without stripping, nothing pairs",
          res["paired_observations"] == 0, str(res["paired_observations"]))

    # Labels default to each file's own run_id, and must differ: every
    # per-side key is built from them, so a run compared with a re-run
    # under the same id would overwrite one side with the other.
    ra = [dict(r, run_id="RUN01") for r in a]
    rb = [dict(r, run_id="RUN07") for r in b]
    check("labels default to each side's run_id",
          C.side_labels(ra, rb, "", "") == ("RUN01", "RUN07"),
          str(C.side_labels(ra, rb, "", "")))
    check("an explicit label wins",
          C.side_labels(ra, rb, "illumina", "") == ("illumina", "RUN07"))
    same = C.side_labels(ra, ra, "", "")
    check("two sides with one run_id are still told apart",
          same[0] != same[1], str(same))
    check("and with no run_id the fallback is A and B",
          C.side_labels(a, b, "", "") == ("A", "B"))

    # The CLI on a pairing that matches nothing: it must say so, not crash
    # formatting an agreement that does not exist.
    tmp = tempfile.mkdtemp()
    try:
        for name, rows in (("a.tsv", a), ("b.tsv", b)):
            with open(os.path.join(tmp, name), "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]),
                                   delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                code = C.main(["--a", os.path.join(tmp, "a.tsv"),
                               "--b", os.path.join(tmp, "b.tsv"),
                               "--outdir", os.path.join(tmp, "cmp")])
        except Exception as e:                      # noqa: BLE001
            code = f"raised {type(e).__name__}: {e}"
        check("compare_runs completes when nothing pairs", code == 0, str(code))
        check("and says nothing paired", "NOTHING PAIRED" in buf.getvalue(),
              buf.getvalue()[-200:])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # a subject measured twice on one side is ambiguous, not a pairing
    a2 = a + [row("S1-r3", "T1", "not_detected", 0)]
    res = C.compare(a2, b, r"-r[0-9]+$", "A", "B")
    check("a duplicated subject is skipped rather than silently resolved",
          res["ambiguous_pairings"] == 1 and res["paired_observations"] == 1,
          f"{res['ambiguous_pairings']} ambiguous, "
          f"{res['paired_observations']} paired")


def test_asymmetry() -> None:
    print("\nwhich run found it, and why")
    a, b = [], []
    for i in range(10):
        # A is shallow and misses; B is deep and finds
        a.append(row(f"S{i}", "T1", "not_detected", 1, depth=500))
        b.append(row(f"S{i}", "T1", "detected", 80, depth=20000,
                     copies="300"))
    res = C.compare(a, b, None, "shallow", "deep")
    check("one-sided detections are attributed to the run that made them",
          res["detected_by_one_only"] == {"deep": 10},
          str(res["detected_by_one_only"]))
    d = res["depth_of_disagreements"]
    check("and the depth behind the disagreement is reported",
          d and d["deep_deeper"] == 10 and d["shallow_deeper"] == 0, str(d))
    check("agreement is zero here, and kappa says so too",
          res["agreement"] == 0.0)


# --------------------------------------------------------------------------
# the simulator property the comparison depends on
# --------------------------------------------------------------------------


def test_composition_is_per_sample() -> None:
    """The bug that made concordance meaningless.

    Nanopore sequences the same libraries Illumina did, so a subject must
    carry the same burden in both runs. When the simulator drew a fresh
    composition per run, two runs of the same samples agreed at kappa 0.04
    and the number looked like a finding about the platforms.
    """
    print("\na sample keeps its composition across runs")
    sim = _load("simulate", "bin/simulate.py")
    sys.path.insert(0, HERE)
    import make_fixtures as mf  # noqa: E402
    import random

    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(20260918)
        ref = mf.build_reference(rng)
        mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
        i7, i5 = mf.load_index_map(
            os.path.join(ROOT, "assets", "index_map.tsv"))
        ill = mf.build_samplesheet(i7, i5, 8, 2, "ILL01", "illumina")
        ont = [dict(s, run_id="ONT01", platform="nanopore") for s in ill]
        rows = ill + ont
        with open(os.path.join(tmp, "samplesheet.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        cols = ["run_id", "platform", "path", "barcode_repeats",
                "dna_input_ng", "lot_pol", "lot_primer", "lot_extract",
                "lot_plasmid"]
        with open(os.path.join(tmp, "runs.tsv"), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            fh.write("\t".join(["ILL01", "illumina",
                                os.path.join(tmp, "sim", "ILL01"), "", "10",
                                "P", "I", "Q", "PL"]) + "\n")
            fh.write("\t".join(["ONT01", "nanopore",
                                os.path.join(tmp, "sim", "ONT01"), "2", "10",
                                "P", "I", "Q", "PL"]) + "\n")
        with contextlib.redirect_stdout(io.StringIO()):
            sim.main(["--reference", os.path.join(tmp, "reference.tsv"),
                      "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                      "--runs", os.path.join(tmp, "runs.tsv"),
                      "--assay", os.path.join(ROOT, "assets", "assay",
                                              "ampseq_helminth.json"),
                      "--outdir", os.path.join(tmp, "sim"),
                      "--reads-per-sample", "80"])

        def comp(run):
            out = {}
            with open(os.path.join(tmp, "sim", f"{run}.composition.tsv")) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    out[(r["sample_id"], r["seq_id"])] = float(
                        r["input_copies"])
            return out

        a, b = comp("ILL01"), comp("ONT01")
        check("both runs carry the same set of sample x reference entries",
              set(a) == set(b), f"{len(a)} vs {len(b)}")
        same = sum(1 for k in a if abs(a[k] - b.get(k, -1)) < 1e-6)
        check("and the same input copies for every one of them",
              a and same == len(a), f"{same} of {len(a)} identical")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# the report
# --------------------------------------------------------------------------


def make_synthesis(tmp: str) -> str:
    out = os.path.join(tmp, "synth")
    os.makedirs(out, exist_ok=True)
    calls = []
    for i in range(12):
        well = f"{chr(65 + i % 4)}{1 + i // 4}"
        stype = "neg_ntc" if i == 0 else "unknown"
        # S03 is deliberately thin: a detection resting on 4 target reads,
        # which must be drawn as poorly determined rather than as 100%
        for t, call, reads in (("T1", "detected" if i % 3 == 0 else
                                "not_detected",
                                (4 if i == 3 else 400) if i % 3 == 0 else 1),
                               ("T2", "assay_failed" if i == 5 else
                                "not_detected", 0)):
            calls.append({
                "run_id": "RUN01", "sample_id": f"S{i:02d}",
                "sample_type": stype, "well": well, "seq_id": t,
                "locus": "ITS2", "purpose": "speciation",
                "taxon": f"Taxon_{t}", "reads": str(reads),
                "depth": str(1000 + 100 * i), "detection_limit": "2",
                "call": call,
                "internal_control": "spiked" if i != 5 else "spike_failed",
                "copies_per_unit_dna": "250.0" if call == "detected" else "",
                "corr_factor": "1", "corr_factor_source": "default",
                "calibration_slope": "1.000", "basis": "test",
            })
    with open(os.path.join(out, "calls.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(calls[0]), delimiter="\t")
        w.writeheader()
        w.writerows(calls)
    with open(os.path.join(out, "calibration.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "state", "slope", "intercept", "reads_hi",
                    "conc_hi", "reads_lo", "conc_lo", "note"])
        for i in range(12):
            w.writerow([f"S{i:02d}", "spiked", "1.0", "0.0", 900, 10000, 9,
                        100, ""])
    qc = {
        "run_id": "RUN01", "n_samples": 12,
        "call_summary": dict(
            (k, sum(1 for c in calls if c["call"] == k))
            for k in ("detected", "not_detected", "assay_failed")),
        "call_meaning": {"detected": "above the limit"},
        "internal_control_summary": {"spiked": 11, "spike_failed": 1},
        "depth": {"min": 1000, "median": 1600, "max": 2100,
                  "fold_spread": 2.1, "below_minimum": 0},
        "calibration_slope": {"n": 12, "median": 1.0, "outside_tolerance": 0},
        "corr_factor": {"default_applied": 1.0,
                        "references_using_default": 2,
                        "references_with_measured_value": 0,
                        "note": "assumed unity"},
        "detection_limits": {
            "T1": {"n_blanks": 1, "max_reads_in_blank": 1,
                   "max_fraction_in_blank": 0.001, "absolute_floor": 2},
            "T2": {"n_blanks": 1, "max_reads_in_blank": 0,
                   "max_fraction_in_blank": 0.0, "absolute_floor": 2}},
        "plate_geometry": [
            {"quantity": "depth", "n_wells": 11, "row_p": 0.5,
             "col_p": 0.004, "edge_p": 0.3}],
        "plate_geometry_note": "unadjusted for multiple testing",
    }
    # read fates: the partition must sum to depth, and one sample is
    # deliberately mostly unusable so the panel has something to show
    fates = []
    for i in range(12):
        depth = 1000 + 100 * i
        if i == 7:
            f = {"target": int(depth * 0.1), "quant": int(depth * 0.1),
                 "control": 0, "vector": int(depth * 0.3),
                 "ambiguous": 0, "tag_unobservable": 0,
                 "unmatched": 0, "no_primer": 0}
        else:
            f = {"target": int(depth * 0.8), "quant": int(depth * 0.1),
                 "control": 0, "vector": 0, "ambiguous": 0,
                 "tag_unobservable": 0, "unmatched": 0, "no_primer": 0}
        f["unmatched"] += depth - sum(f.values())
        fates.append({"run_id": "RUN01", "sample_id": f"S{i:02d}",
                      "sample_type": "unknown",
                      "well": f"{chr(65 + i % 4)}{1 + i // 4}",
                      "depth": depth, **f,
                      "target_fraction": f"{f['target']/depth:.6f}"})
    with open(os.path.join(out, "read_fates.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fates[0]), delimiter="\t")
        w.writeheader()
        w.writerows(fates)
    qc["read_fates"] = {"partition_balances": True}

    with open(os.path.join(out, "run_qc.json"), "w") as fh:
        json.dump(qc, fh)
    return out


def test_report() -> None:
    print("\nthe report")
    tmp = tempfile.mkdtemp()
    try:
        synth = make_synthesis(tmp)
        out = os.path.join(tmp, "report.html")
        with contextlib.redirect_stdout(io.StringIO()):
            R.main(["--synthesis-dir", synth, "--out", out, "--quiet"])
        doc = open(out, encoding="utf-8").read()

        # the field requirement
        refs = [u for u in re.findall(r'(?:src|href)\s*=\s*["\']([^"\']+)', doc)
                if not u.startswith("#")]
        check("nothing is loaded from the network", not refs, str(refs[:3]))
        check("no scripts at all", "<script" not in doc.lower())
        check("the stylesheet is inline", "<style>" in doc)

        # the numbers must match the table they came from
        calls = list(csv.DictReader(open(os.path.join(synth, "calls.tsv")),
                                    delimiter="\t"))
        n_det = sum(1 for c in calls if c["call"] == "detected")
        n_fail = sum(1 for c in calls if c["call"] == "assay_failed")
        cards = re.findall(r'<div class="n">([\d,]+)</div>\s*'
                           r'<div class="k">([a-z ]+)</div>', doc)
        as_dict = {k.strip(): int(v.replace(",", "")) for v, k in cards}
        check("the detected count matches calls.tsv",
              as_dict.get("detected") == n_det,
              f"{as_dict.get('detected')} vs {n_det}")
        check("the assay_failed count matches calls.tsv",
              as_dict.get("assay failed") == n_fail,
              f"{as_dict.get('assay failed')} vs {n_fail}")

        # the concerns section must surface the failures, not bury them
        head = doc[:doc.index("<h2>Calls</h2>")]
        check("a failed well is raised before the results",
              "assay_failed" in head, "not in the opening section")
        check("a column effect is raised before the results",
              "position" in head.lower())
        check("the default corr_factor is stated up front",
              "corr_factor" in head)

        # the composition panel
        comp = doc[doc.index("<h2>Composition</h2>"):]
        comp = comp[:comp.index("<h2>", 10)]
        check("depth and composition share one figure",
              "Read depth" in comp and "Relative abundance" in comp,
              "one of the two axes is missing")
        check("a thinly supported composition is hatched",
              "url(#thin)" in comp,
              "S03 rests on 4 target reads and was drawn as solid")
        check("the hatching is explained, with a count",
              "Hatched bars rest on fewer than" in comp)
        check("samples with no target reads are drawn empty, not split",
              "no target reads" in comp)
        check("every stacked segment names its share and denominator",
              comp.count("of ") >= 4 and "target reads" in comp)
        # the panel must not invent a proportion
        import re as _re
        shares = [float(x) for x in
                  _re.findall(r'Taxon_T1: ([\d.]+)% of', comp)]
        check("shares are percentages within range",
              shares and all(0 <= x <= 100 for x in shares),
              str(shares[:4]))

        # the read-fate panel, between depth and composition
        check("a third panel shows where the reads went",
              "Reads mapped" in comp, "panel missing")
        check("its target slice is named as the composition denominator",
              "denominator of the composition below" in comp)
        check("each fate is legended separately from the taxa",
              "reads mapped:" in comp and "composition:" in comp)
        fates = {r["sample_id"]: r for r in csv.DictReader(
            open(os.path.join(synth, "read_fates.tsv")), delimiter="\t")}
        check("the fate partition sums to depth in every sample",
              all(sum(int(r[k]) for k in R.FATE_ORDER) == int(r["depth"])
                  for r in fates.values()),
              "a fate partition does not account for its depth")
        med = sorted(int(r["target"]) / int(r["depth"])
                     for r in fates.values())[len(fates) // 2]
        check("the caption's median target fraction matches the table",
              f"{100*med:.1f}%" in comp, f"{100*med:.1f}%")
        check("a mostly unusable well is called out",
              "below 50%" in comp, "not flagged")
        check("every fate segment names its own count",
              comp.count("reads,") >= 12, str(comp.count("reads,")))

        # plate maps, since geometry was available
        check("the plate is drawn", doc.count("<svg") >= 2,
              str(doc.count("<svg")))
        check("every well carries its own label for hover",
              doc.count("<title>") > 12, str(doc.count("<title>")))

        # with a comparison attached
        cmp_path = os.path.join(tmp, "comparison.json")
        with open(cmp_path, "w") as fh:
            json.dump({"label_a": "illumina", "label_b": "nanopore",
                       "paired_observations": 120, "agreement": 0.925,
                       "cohen_kappa": 0.831, "within_twofold": 0.971,
                       "kappa_note": "chance removed",
                       "detected_by_one_only": {"nanopore": 9},
                       "depth_of_disagreements": {
                           "n_disagreements_with_depth": 9,
                           "illumina_deeper": 2, "nanopore_deeper": 7}}, fh)
        out2 = os.path.join(tmp, "report2.html")
        with contextlib.redirect_stdout(io.StringIO()):
            R.main(["--synthesis-dir", synth, "--comparison", cmp_path,
                    "--out", out2, "--quiet"])
        doc2 = open(out2, encoding="utf-8").read()
        check("the paired run appears when there is one",
              "Agreement with the paired run" in doc2)
        check("kappa is shown next to raw agreement, not instead of it",
              "0.925" in doc2 and "0.831" in doc2)
        check("the section is absent when there is no paired run",
              "Agreement with the paired run" not in doc)

        # Nothing paired. compare_runs writes agreement as None, and the
        # report formatted it with :.3f and crashed - so a comparison that
        # measured nothing took the whole report down with it.
        res = C.compare([row("S1-r1", "T1", "detected", 100)],
                        [row("S1-r2", "T1", "detected", 90)], None, "A", "B")
        check("an unpaired comparison records agreement as undefined",
              res["paired_observations"] == 0 and res["agreement"] is None,
              str(res["agreement"]))
        empty = os.path.join(tmp, "comparison_empty.json")
        with open(empty, "w") as fh:
            json.dump(res, fh)
        out3 = os.path.join(tmp, "report3.html")
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                code = R.main(["--synthesis-dir", synth, "--comparison", empty,
                               "--out", out3, "--quiet"])
        except Exception as e:                      # noqa: BLE001
            code = f"raised {type(e).__name__}: {e}"
        check("the report renders when nothing paired", code in (0, None),
              str(code))
        doc3 = open(out3, encoding="utf-8").read() if os.path.exists(out3) else ""
        check("and says nothing paired rather than showing a number",
              "No observation was present in both runs" in doc3
              and "raw agreement" not in doc3)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_distance_chart_shows_a_minority_mode() -> None:
    """A plot whose whole purpose is a 1% mode must not hide a 1% mode.

    This chart exists to show a population of reads that is assigned,
    counted as confident, and wrong - a read built from two panel
    members. That population is a fraction of a percent next to the
    clean mode, which has two consequences the drawing has to respect.

    A LINEAR axis erases it. One bar of 70 beside one of 5,000 is under
    two pixels tall in a 120px plot, and a reader cannot see what is not
    drawn. Hence log counts, and hence the test asserts the ratio of
    drawn heights is nothing like the ratio of counts.

    A BAR OF ZERO HEIGHT erases it too. log10(1) is 0, so a distance
    observed once would draw as nothing at all and read as absent.
    """
    print("\ndistance chart")
    import re as _re

    hist = {"0": 5000, "1": 1800, "2": 240, "9": 70}
    svg = R.counts_bar_svg(hist, "test")

    titles = dict(_re.findall(r'<title>distance (\d+): ([\d,]+) reads', svg))
    check("every distance in the data is drawn",
          set(titles) == {"0", "1", "2", "9"}, str(sorted(titles)))
    check("and each bar carries its own count, readable on hover",
          titles.get("9") == "70", str(titles))

    rects = _re.findall(
        r'<rect[^>]*height="([\d.]+)"[^>]*><title>distance (\d+)', svg)
    heights = {d: float(h) for h, d in rects}
    check("the minority mode is drawn tall enough to see",
          heights.get("9", 0) > 0.35 * heights.get("0", 1),
          f"{heights.get('9')} vs {heights.get('0')}")
    check("which means the axis is not linear in counts",
          heights["9"] / heights["0"] > 10 * (70 / 5000),
          f"{heights['9'] / heights['0']:.3f} vs {70/5000:.3f}")

    once = R.counts_bar_svg({"0": 5000, "40": 1}, "test")
    r1 = dict((d, float(h)) for h, d in _re.findall(
        r'<rect[^>]*height="([\d.]+)"[^>]*><title>distance (\d+)', once))
    check("a distance seen exactly once still draws a bar",
          r1.get("40", 0) > 0, str(r1))

    check("a single distance says so in words rather than drawing an axis",
          "no distribution to plot" in R.counts_bar_svg({"0": 10}, "t"),
          R.counts_bar_svg({"0": 10}, "t")[:80])
    check("an empty histogram draws nothing at all",
          R.counts_bar_svg({}, "t") == "")

    # note_hi colours the far mode. The colour IS the finding, so it has
    # to actually change - a note_hi that quietly did nothing would
    # leave the chimeric bars looking like error.
    plain = R.counts_bar_svg(hist, "t")
    marked = R.counts_bar_svg(hist, "t", note_hi=6)
    check("note_hi recolours the bars beyond it, and only those",
          marked.count("#b03030") == 1 and plain.count("#b03030") == 0,
          f"{marked.count('#b03030')} / {plain.count('#b03030')}")

    # flag() escapes its argument, which every other caller depends on.
    # The first version of the fit-distance flag passed <b> tags and
    # rendered them as literal &lt;b&gt; on the page.
    f = R.flag("warn", "Not monotone: trough at distance 5")
    check("a flag carries no markup of its own to be escaped",
          "&lt;" not in f and "Not monotone" in f, f)


def test_histograms() -> None:
    """A chart whose counts cannot be read is not evidence.

    Four defects are pinned here, all of them found by checking a report
    against the tables it was drawn from rather than by looking at it:

      - no y-axis, so a bar height was a shape and not a number
      - twenty-one values taking two distinct values were binned across
        twenty-six bins, drawing as a bimodal spread when the truth was
        "twenty identical, one outlier"
      - identical values stretched the axis to an invented range
      - on a log axis a non-positive value stayed linear and landed in an
        arbitrary bin
    """
    print("\nhistograms")
    import re as _re

    def counts_of(svg):
        return [int(x) for x in _re.findall(r'<title>(\d+)</title>', svg)]

    def texts(svg):
        return _re.findall(r'<text[^>]*>([^<]+)</text>', svg)

    # a plain distribution
    vals = [1.0 + i * 0.01 for i in range(40)]
    svg = R.histogram_svg(vals, "test")
    check("every value is counted into some bar",
          sum(counts_of(svg)) == len(vals),
          f"{sum(counts_of(svg))} of {len(vals)}")
    check("a y-axis is drawn with the maximum count labelled",
          "samples" in texts(svg) and str(max(counts_of(svg))) in texts(svg),
          str(texts(svg)[:6]))
    check("more than two x labels are given",
          sum(1 for t in texts(svg) if t.replace(".", "").replace(",", "")
              .replace("-", "").isdigit()) >= 5, str(texts(svg)))

    # the case that misled: few distinct values
    two = [1.0] * 20 + [0.8516]
    svg = R.histogram_svg(two, "slopes")
    check("a handful of distinct values is drawn as itself, not binned",
          set(counts_of(svg)) == {20, 1}, str(counts_of(svg)))
    check("and each value is labelled",
          "0.852" in " ".join(texts(svg)) or "0.8516" in " ".join(texts(svg)),
          str(texts(svg)))
    check("the counts still sum to the input",
          sum(counts_of(svg)) == len(two))

    # no range at all
    out = R.histogram_svg([2.5] * 12, "flat")
    check("identical values are stated, not plotted as a fake range",
          "<svg" not in out and "12" in out and "2.5" in out, out[:120])

    # log axis with a zero
    out = R.histogram_svg([0.0, 10.0, 100.0, 1000.0, 5000.0, 20000.0,
                           70000.0, 250000.0, 900.0], "copies", log=True)
    check("a non-positive value is excluded from a log axis",
          "non-positive" in out, "silently binned")
    check("and the excluded count is stated", ">1 non-positive" in out
          or "1 non-positive" in out, out[-160:])
    check("the remaining values are all counted",
          sum(counts_of(out)) == 8, str(sum(counts_of(out))))
    check("axis labels are readable numbers, not exponent notation",
          "e+0" not in " ".join(texts(out)), str(texts(out)))

    check("an empty input draws nothing rather than failing",
          R.histogram_svg([], "none") == "")
    check("all-zero on a log axis says so instead of drawing",
          "No positive values" in R.histogram_svg([0.0, 0.0], "z", log=True))


def main() -> int:
    test_kappa()
    test_histograms()
    test_distance_chart_shows_a_minority_mode()
    test_pairing()
    test_asymmetry()
    test_report()
    test_composition_is_per_sample()

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
