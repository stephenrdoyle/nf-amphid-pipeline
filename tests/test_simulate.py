#!/usr/bin/env python3
"""
Tests for the read simulator.

A simulator that quietly produces the wrong thing is worse than none, since
every downstream test would then be validated against it. So these check
that what comes out matches what was asked for: injected rates land within
sampling error, the truth table accounts for every read, and the biological
constraints that make the assay work are actually honoured - a mock contains
no real target sequence, a blank contains only spike.

Run: python3 tests/test_simulate.py
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


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


def near(observed: float, expected: float, tol: float) -> bool:
    return abs(observed - expected) <= tol


# --------------------------------------------------------------------------
# unit
# --------------------------------------------------------------------------


def test_construct() -> None:
    print("\nconstruct assembly")
    a = sim.DEFAULTS
    mol = sim.build_molecule("AAAACCCCGGGG", "TAAGGCGA", "CTCTCTAT", 2, a)

    check("starts with the P5 tail", mol.startswith(a["p5_tail"]))
    check("ends with the reverse-complemented P7 tail",
          mol.endswith(sim.revcomp(a["p7_tail"])))
    check("i5 appears tandemly after P5",
          a["p5_tail"] + "CTCTCTAT" + "CTCTCTAT" in mol)
    check("i7 appears tandemly before P7",
          "TAAGGCGA" + "TAAGGCGA" + sim.revcomp(a["p7_tail"]) in mol)
    check("PCR1 forward overhang precedes the amplicon",
          a["pcr1_overhang_fwd"] + "AAAACCCCGGGG" in mol)
    check("amplicon is followed by the reverse-complemented rev overhang",
          "AAAACCCCGGGG" + sim.revcomp(a["pcr1_overhang_rev"]) in mol)

    single = sim.build_molecule("AAAA", "TAAGGCGA", "CTCTCTAT", 1, a)
    check("repeats=1 gives one copy of each index",
          single.count("CTCTCTAT") == 1 and single.count("TAAGGCGA") == 1)
    check("repeats=1 is shorter by exactly two indices",
          len(mol) - len(sim.build_molecule("AAAACCCCGGGG", "TAAGGCGA", "CTCTCTAT", 1, a)) == 16)


def test_error_model() -> None:
    print("\nerror model")
    import random

    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(5000))

    for name, expected in (("illumina", 0.001), ("nanopore_hac", 0.02)):
        model = sim.ERROR_MODELS[name]
        tot_sub = tot_indel = 0
        for _ in range(20):
            _, s, i = sim.apply_errors(rng, seq, model)
            tot_sub += s
            tot_indel += i
        observed = (tot_sub + tot_indel) / (len(seq) * 20)
        check(
            f"{name} error rate ~{expected}",
            near(observed, expected, expected * 0.35),
            f"{observed:.5f}",
        )
        if name == "nanopore_hac":
            frac = tot_indel / max(1, tot_sub + tot_indel)
            check("nanopore errors are roughly half indels",
                  near(frac, 0.5, 0.15), f"{frac:.2f}")

    # Homopolymers must attract indels, since that is what defeats a
    # Hamming-based matcher and justifies using edit distance.
    # Note "ACGT"*20 ends in T, so the run below is 9, not 8.
    hp = "ACGT" * 20 + "TTTTTTTT" + "GCA"
    runs = sim.homopolymer_lengths(hp)
    check("homopolymer run length computed", max(runs) == 9, str(max(runs)))
    check("run length resets after the homopolymer", runs[-1] == 1, str(runs[-1]))
    check("a non-repeating stretch is all ones",
          set(sim.homopolymer_lengths("ACGTACGT")) == {1})

    # and the bias must actually show up in the output
    import random

    rng2 = random.Random(7)
    model = dict(sim.ERROR_MODELS["nanopore_hac"])
    flat = "ACGT" * 250
    poly = "".join(b * 4 for b in "ACGT" * 62) + "AC"
    ef = sum(sim.apply_errors(rng2, flat, model)[2] for _ in range(40))
    ep = sum(sim.apply_errors(rng2, poly, model)[2] for _ in range(40))
    check("indels concentrate in homopolymers", ep > ef * 1.5, f"{ep} vs {ef}")


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


def build(tmp: str, platform: str, n: int, unused: int) -> None:
    import random

    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, n, unused, "RUN01", platform)
    mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
    mf.write_samplesheet(samples, os.path.join(tmp, "samplesheet.csv"))
    mf.write_runs(
        os.path.join(tmp, "runs.tsv"), "RUN01", platform, os.path.join(tmp, "run01")
    )


def run_sim(tmp: str, extra: list[str] | None = None) -> list[dict]:
    out = os.path.join(tmp, "sim")
    argv = [
        "--reference", os.path.join(tmp, "reference.tsv"),
        "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
        "--runs", os.path.join(tmp, "runs.tsv"),
        "--assay", ASSAY,
        "--outdir", out,
        "--reads-per-sample", "400",
    ] + (extra or [])
    # keep the simulator's own stdout out of the test output
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        sim.main(argv)
    with open(os.path.join(out, "RUN01.truth.tsv"), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_rates() -> None:
    print("\ninjected rates are recovered")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 24, 4)
        truth = run_sim(tmp)
        n = len(truth)
        check("reads produced", n > 5000, str(n))

        fates = Counter(r["fate"] for r in truth)
        check("every fate is a known category",
              set(fates) <= set(sim.FATES), str(set(fates) - set(sim.FATES)))

        for fate, expected in (
            ("hopped", 0.005),
            ("chimera", 0.01),
            ("carryover", 0.002),
            ("dimer", 0.02),
        ):
            obs = fates[fate] / n
            check(
                f"{fate} rate ~{expected}",
                near(obs, expected, max(0.004, expected * 0.5)),
                f"{obs:.4f}",
            )

        # every read accounted for exactly once - the invariant at source
        ids = [r["read_id"] for r in truth]
        check("read ids unique", len(set(ids)) == len(ids))
        check("fates partition the reads", sum(fates.values()) == n)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_biology() -> None:
    print("\nbiological constraints")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 24, 4)
        truth = run_sim(tmp)

        ref = {
            r["seq_id"]: r["type"]
            for r in csv.DictReader(
                open(os.path.join(tmp, "reference.tsv")), delimiter="\t"
            )
        }
        ss = {
            s["sample_id"]: s["sample_type"]
            for s in csv.DictReader(open(os.path.join(tmp, "samplesheet.csv")))
        }

        def types_in(sample_type: str, only_clean: bool = True) -> set[str]:
            return {
                ref.get(r["seq_id"], "?")
                for r in truth
                if ss.get(r["sample_id"]) == sample_type
                and (r["fate"] == "clean" if only_clean else True)
            }

        # This is the property the whole control design rests on: a mock is
        # built from plasmids, so it contains no real target sequence and
        # mock carryover can never look like a real positive.
        check(
            "mocks contain no target sequence",
            "target" not in types_in("mock"),
            str(types_in("mock")),
        )
        check("mocks do contain controls", "control" in types_in("mock"))

        for blank in ("neg_extraction", "neg_ntc"):
            seen = types_in(blank)
            check(
                f"{blank} contains only quants",
                seen <= {"quant"},
                str(seen),
            )

        check(
            "field samples carry targets",
            "target" in types_in("unknown"),
        )

        # orientation should be roughly balanced on nanopore, since demux has
        # to search both strands
        orient = Counter(r["orientation"] for r in truth)
        frac = orient["-"] / len(truth)
        check("nanopore reads are ~50% reverse", near(frac, 0.5, 0.05), f"{frac:.3f}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_hopping_lands_in_sentinels() -> None:
    print("\nhopped reads and sentinels")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 24, 4)
        truth = run_sim(tmp)

        assigned = {
            (s["index_i5"], s["index_i7"])
            for s in csv.DictReader(open(os.path.join(tmp, "samplesheet.csv")))
        }
        hopped = [r for r in truth if r["fate"] == "hopped"]
        into_sentinel = [
            r for r in hopped if (r["index_i5"], r["index_i7"]) not in assigned
        ]
        check("some hopped reads exist", len(hopped) > 10, str(len(hopped)))
        check(
            "some hopped reads land on unassigned combinations",
            len(into_sentinel) > 0,
            f"{len(into_sentinel)}/{len(hopped)}",
        )
        check(
            "hopped reads are attributed to the receiving sample, not the source",
            all(r["sample_id"] != "" for r in hopped),
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_misread_one_index() -> None:
    """--misread-i5/i7-rate: one index read as another of the plate's,
    so the read keeps the other end - the same column or the same row -
    unlike hopping, which rewrites the pair. Off unless asked for."""
    print("\none index misread")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 24, 4)
        truth = run_sim(tmp)
        check("no misreads unless asked for",
              not any(r["fate"].startswith("misread") for r in truth))
        truth = run_sim(tmp, ["--hop-rate", "0", "--misread-i5-rate", "0.02",
                              "--misread-i7-rate", "0.01"])
        n = len(truth)
        f = Counter(r["fate"] for r in truth)
        clean_ish = n * (1 - 0.02 - 0.002 - 0.005)
        check("i5 misreads at the rate asked for",
              near(f["misread_i5"] / clean_ish, 0.02, 0.006),
              f"{f['misread_i5']}/{n}")
        check("i7 misreads at the rate asked for",
              near(f["misread_i7"] / clean_ish, 0.0098, 0.004),
              f"{f['misread_i7']}/{n}")
        sheet = list(csv.DictReader(open(os.path.join(tmp,
                                                      "samplesheet.csv"))))
        by_pair = {(s["index_i5"], s["index_i7"]): s["sample_id"]
                   for s in sheet}
        i5s = {s["index_i5"] for s in sheet}
        i7s = {s["index_i7"] for s in sheet}
        mis = [r for r in truth if r["fate"].startswith("misread")]
        check("a misread index is another of the plate's own",
              all(r["index_i5"] in i5s and r["index_i7"] in i7s
                  for r in mis))
        check("and the read goes to the well with that pair, or to none",
              all(r["sample_id"] == by_pair.get(
                  (r["index_i5"], r["index_i7"]), "_unassigned")
                  for r in mis))
        check("some land on an unused combination",
              any(r["sample_id"] == "_unassigned" for r in mis))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_outputs() -> None:
    print("\noutput files")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "illumina", 12, 3)
        truth = run_sim(tmp)
        d = os.path.join(tmp, "sim", "RUN01")
        files = sorted(os.listdir(d))
        check("illumina output is pre-split per sample",
              all(f.endswith(".fastq.gz") for f in files) and len(files) >= 20,
              str(len(files)))
        check("R1 and R2 both present",
              sum("_R1_" in f for f in files) == sum("_R2_" in f for f in files))

        # FASTQ must be well formed, and read ids must match the truth table
        some = [f for f in files if "_R1_" in f][0]
        with gzip.open(os.path.join(d, some), "rt") as fh:
            lines = fh.read().splitlines()
        check("fastq has 4 lines per record", len(lines) % 4 == 0)
        check("sequence and quality lengths agree",
              all(len(lines[i + 1]) == len(lines[i + 3]) for i in range(0, len(lines), 4)))
        ids = {lines[i][1:].split()[0] for i in range(0, len(lines), 4)}
        truth_ids = {r["read_id"] for r in truth}
        check("every fastq read id is in the truth table", ids <= truth_ids,
              str(list(ids - truth_ids)[:3]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 12, 3)
        run_sim(tmp)
        pooled = os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz")
        check("nanopore output is a single pooled fastq", os.path.exists(pooled))
        with gzip.open(pooled, "rt") as fh:
            head = fh.readline()
        # Dorado writes SAM-style tags separated by TABS, and names the
        # model in RG. The fixture used to emit "model_version_id=" after
        # a space, which no parser written against real output would find.
        check("nanopore header carries dorado's RG tag",
              "\tRG:Z:" in head, head.strip()[:80])
        rg = [f for f in head.strip().split("\t") if f.startswith("RG:Z:")]
        # And the model named must be the model simulated. It used to say
        # sup regardless of the error rate applied.
        check("and names the model that was actually simulated",
              rg and "_hac@" in rg[0], rg[0][:80] if rg else "no RG tag")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_no_spike() -> None:
    print("\nunspiked mode")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 24, 4)
        truth = run_sim(tmp, ["--no-spike"])
        ref = {
            r["seq_id"]: r["type"]
            for r in csv.DictReader(
                open(os.path.join(tmp, "reference.tsv")), delimiter="\t"
            )
        }
        quant_reads = [r for r in truth if ref.get(r["seq_id"]) == "quant"]
        check("no quant reads when unspiked", len(quant_reads) == 0,
              str(len(quant_reads)))
        check("targets still present", any(ref.get(r["seq_id"]) == "target" for r in truth))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_determinism() -> None:
    print("\ndeterminism")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, "nanopore", 12, 3)
        a = run_sim(tmp)
        shutil.rmtree(os.path.join(tmp, "sim"))
        b = run_sim(tmp)
        check("same seed gives identical truth tables", a == b,
              f"{len(a)} vs {len(b)} reads")
        shutil.rmtree(os.path.join(tmp, "sim"))
        c = run_sim(tmp, ["--seed", "99"])
        check("a different seed gives different output", a != c)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_construct()
    test_error_model()
    test_rates()
    test_biology()
    test_hopping_lands_in_sentinels()
    test_misread_one_index()
    test_outputs()
    test_no_spike()
    test_determinism()

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
