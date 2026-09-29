#!/usr/bin/env python3
"""
Tests for concatemer splitting.

The measurement that matters is not how many reads get cut. It is
whether the pieces reach the samples the molecules actually came from,
so every test here is scored against known provenance: chimeras are
built from reads whose true sample is already established, and the
question asked afterwards is whether both halves went home.

The grid shape of the samplesheet is load-bearing and easy to get wrong
in a fixture. If every sample shares one i5, a head-to-tail join of two
molecules produces the pair (shared i5, second molecule's i7) - which is
just the second sample, so the join looks harmless and the test proves
nothing. The contamination only appears when i5 varies, so the fixture
below uses a grid and asserts that it did.

Run: python3 tests/test_split.py
"""

from __future__ import annotations

import contextlib
import csv
import glob
import gzip
import importlib.util
import io
import json
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


dx = _load("demux_ont", "bin/demux_ont.py")
sc = _load("split_concatemers", "bin/split_concatemers.py")
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


def quiet(fn, *a):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a)


def build(tmp, n_samples=24, unused=4, reads=120):
    """A run, its truth table, and chimeras of known provenance."""
    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, n_samples, unused, "RUN01",
                                   "nanopore")
    mf.write_reference(ref, f"{tmp}/reference.tsv")
    mf.write_samplesheet(samples, f"{tmp}/samplesheet.csv")
    mf.write_runs(f"{tmp}/runs.tsv", "RUN01", "nanopore", f"{tmp}/sim/RUN01")
    quiet(sim.main, ["--reference", f"{tmp}/reference.tsv",
                     "--samplesheet", f"{tmp}/samplesheet.csv",
                     "--runs", f"{tmp}/runs.tsv", "--assay", ASSAY,
                     "--outdir", f"{tmp}/sim", "--reads-per-sample", str(reads)])
    quiet(dx.main, ["--reads", f"{tmp}/sim/RUN01",
                    "--samplesheet", f"{tmp}/samplesheet.csv",
                    "--run-id", "RUN01", "--assay", ASSAY,
                    "--outdir", f"{tmp}/d0"])
    with open(f"{tmp}/d0/assignments.tsv", encoding="utf-8") as fh:
        truth = {r["read_id"]: r["sample_id"]
                 for r in csv.DictReader(fh, delimiter="\t") if r["sample_id"]}

    with gzip.open(f"{tmp}/sim/RUN01/RUN01.fastq.gz", "rt") as fh:
        L = fh.read().splitlines()
    recs = [(L[i][1:].split()[0], L[i + 1], L[i + 3])
            for i in range(0, len(L), 4)]
    good = [r for r in recs if r[0] in truth]

    rc, rng2, prov = dx.revcomp, random.Random(5), {}
    os.makedirs(f"{tmp}/chim", exist_ok=True)
    with gzip.open(f"{tmp}/chim/c.fastq.gz", "wt") as fh:
        for k in range(0, len(good) - 1, 2):
            a, b = good[k], good[k + 1]
            inv = rng2.random() < 0.4          # some joins head to head
            fh.write(f"@J{k}\n{a[1] + (rc(b[1]) if inv else b[1])}\n+\n"
                     f"{a[2] + (b[2][::-1] if inv else b[2])}\n")
            prov[f"J{k}"] = [truth[a[0]], truth[b[0]]]
    ss = list(csv.DictReader(open(f"{tmp}/samplesheet.csv")))
    return prov, ss


def calls_of(path):
    with open(path, encoding="utf-8") as fh:
        return {r["read_id"]: r["sample_id"]
                for r in csv.DictReader(fh, delimiter="\t")}


def test_recovers_correct_samples() -> None:
    print("\nsplitting chimeras of known provenance")
    tmp = tempfile.mkdtemp()
    try:
        prov, ss = build(tmp)
        n_i5 = len({s["index_i5"] for s in ss})
        check("the fixture is a grid, so hybrid pairs are possible at all",
              n_i5 > 1, f"only {n_i5} distinct i5 - a join cannot mislabel")

        quiet(sc.main, ["--reads", f"{tmp}/chim", "--assay", ASSAY,
                        "--outdir", f"{tmp}/split", "--quiet"])
        quiet(dx.main, ["--reads", f"{tmp}/split/reads",
                        "--samplesheet", f"{tmp}/samplesheet.csv",
                        "--run-id", "RUN01", "--assay", ASSAY,
                        "--outdir", f"{tmp}/new"])
        new = calls_of(f"{tmp}/new/assignments.tsv")

        tot = len(prov)
        right = wrong = 0
        for pid, (sa, sb) in prov.items():
            got = [new[k] for k in (f"{pid}_u1", f"{pid}_u2", pid)
                   if new.get(k)]
            if sorted(got[:2]) == sorted([sa, sb]):
                right += 1
            elif len(got) >= 2:
                wrong += 1
        check("both halves reach their own samples for the great majority",
              right > 0.90 * tot, f"{right} of {tot} ({100*right/tot:.1f}%)")
        check("and none of the recovered halves is misassigned",
              wrong == 0, f"{wrong} of {tot} went to the wrong sample")

        with open(f"{tmp}/split/split_summary.json") as fh:
            summ = json.load(fh)
        check("no sequence is lost beyond the length floor",
              summ["bases_retained_fraction"] > 0.98,
              f"{100*summ['bases_retained_fraction']:.1f}% of bases retained")
        check("the pieces outnumber the reads, as they must",
              summ["pieces_per_read"] > 1.5,
              str(summ["pieces_per_read"]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_quantifies_what_it_fixes() -> None:
    """The old behaviour, measured rather than asserted about.

    Splitting is worth its complexity only if the thing it prevents is
    real. So the same chimeras are put through the previous code path -
    no junction check, no splitting - and scored against the same
    provenance. The number that matters is reads assigned to a sample
    NEITHER molecule came from, because those are indistinguishable from
    data downstream.
    """
    print("\nwhat splitting prevents")
    tmp = tempfile.mkdtemp()
    try:
        prov, _ss = build(tmp)
        quiet(dx.main, ["--reads", f"{tmp}/chim",
                        "--samplesheet", f"{tmp}/samplesheet.csv",
                        "--run-id", "RUN01", "--assay", ASSAY,
                        "--outdir", f"{tmp}/old", "--no-chimera-check"])
        old = calls_of(f"{tmp}/old/assignments.tsv")
        tot = len(prov)
        bad = sum(1 for pid, (sa, sb) in prov.items()
                  if old.get(pid) and old[pid] not in (sa, sb))
        print(f"  ..   without the junction check, {bad:,} of {tot:,} "
              f"({100*bad/tot:.1f}%) were assigned to a\n       sample "
              "neither molecule came from")
        check("the failure being fixed is real, not hypothetical",
              bad > 0,
              "no chimera was misassigned even without the check - either "
              "the fixture cannot express the failure, or it is gone")

        # And with the check on, the same reads are refused instead.
        quiet(dx.main, ["--reads", f"{tmp}/chim",
                        "--samplesheet", f"{tmp}/samplesheet.csv",
                        "--run-id", "RUN01", "--assay", ASSAY,
                        "--outdir", f"{tmp}/guard"])
        guard = calls_of(f"{tmp}/guard/assignments.tsv")
        still = sum(1 for pid, (sa, sb) in prov.items()
                    if guard.get(pid) and guard[pid] not in (sa, sb))
        check("the guard alone removes that contamination",
              still == 0, f"{still} still misassigned with the guard on")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_clean_reads_untouched() -> None:
    """Splitting must be a no-op on a library that does not need it."""
    print("\nclean reads")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp)
        src = f"{tmp}/sim/RUN01"
        quiet(sc.main, ["--reads", src, "--assay", ASSAY,
                        "--outdir", f"{tmp}/sp", "--quiet"])
        with open(f"{tmp}/sp/split_summary.json") as fh:
            summ = json.load(fh)
        rate = summ["reads_split"] / max(1, summ["reads_in"])
        check("a single-molecule library is left alone",
              rate < 0.01, f"{100*rate:.2f}% of clean reads were cut")

        # Byte-level: the sequences that come out are the ones that went in.
        with gzip.open(f"{src}/RUN01.fastq.gz", "rt") as fh:
            a = fh.read().splitlines()
        before = {a[i][1:].split()[0]: a[i + 1] for i in range(0, len(a), 4)}
        b = []
        for part in sorted(glob.glob(f"{tmp}/sp/reads/*.fastq.gz")):
            with gzip.open(part, "rt") as fh:
                b += fh.read().splitlines()
        after = {b[i][1:].split()[0]: b[i + 1] for i in range(0, len(b), 4)}
        altered = [k for k, v in after.items()
                   if k in before and before[k] != v]
        check("and their sequence is returned unmodified",
              not altered, f"{len(altered)} reads changed, e.g. {altered[:2]}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_recovers_correct_samples()
    test_quantifies_what_it_fixes()
    test_clean_reads_untouched()
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
