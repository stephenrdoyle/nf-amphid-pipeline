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
ms = _load("merge_split", "bin/merge_split.py")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))
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


def test_basecall_provenance() -> None:
    """The basecalling model must reach the outputs, and be the right one.

    A model is not a neutral detail: it sets the error profile, which sets
    barcode recovery, which sets what reaches each sample. Two runs
    basecalled differently can differ in composition for reasons that are
    not biological - the same shape of confounding as comparing runs at
    different depths - and nothing else in these outputs would reveal it.

    Splitting rewrites read ids, so the tags are gone downstream of it.
    That makes the splitter the last stage that can record this, and makes
    "demux reports nothing" the correct behaviour rather than a gap.
    """
    print("\nbasecall provenance")
    import ampseqlib as al

    real = ("@9fc33eca-e571-4ed1-a218-6bfa69740e2c\tqs:f:20.4937\tch:i:378\t"
            "RG:Z:dd4877b8-1550-4719-b148-55edcf2e2394_dna_r10.4.1_e8.2_"
            "400bps_hac@v5.2.0\tPU:Z:FBH87365\tLB:Z:helmseq_primer_opt1_pool")
    tmp = tempfile.mkdtemp()
    try:
        p = f"{tmp}/r.fastq.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(real + "\nACGT\n+\nIIII\n")
        got = al.read_run_metadata(p)
        check("a real dorado header yields the model",
              got.get("basecall_model") == "dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
              str(got.get("basecall_model")))
        check("the run uuid is split off the model, not left joined to it",
              got.get("acquisition_id") ==
              "dd4877b8-1550-4719-b148-55edcf2e2394",
              str(got.get("acquisition_id")))
        check("the flow cell comes through", got.get("flow_cell") == "FBH87365")

        # Absent tags and absent files are answers, not errors: plenty of
        # FASTQ has neither, and this must never be the thing that fails
        # a run.
        with open(f"{tmp}/plain.fastq", "w") as fh:
            fh.write("@read_u1\nACGT\n+\nIIII\n")
        check("a header with no tags returns nothing, quietly",
              al.read_run_metadata(f"{tmp}/plain.fastq") == {})
        check("a missing file returns nothing, quietly",
              al.read_run_metadata(f"{tmp}/nope.fastq.gz") == {})

        # End to end, and specifically that the model RECORDED is the model
        # SIMULATED. The fixture used to name sup whatever error rate had
        # been applied, so this would have passed while being wrong.
        for variant in ("fast", "sup"):
            d = os.path.join(tmp, variant)
            os.makedirs(d, exist_ok=True)
            rng = random.Random(7)
            ref = mf.build_reference(rng)
            i7, i5 = mf.load_index_map(INDEX_MAP)
            mf.write_reference(ref, f"{d}/reference.tsv")
            mf.write_samplesheet(
                mf.build_samplesheet(i7, i5, 8, 2, "RUN01", "nanopore"),
                f"{d}/samplesheet.csv")
            mf.write_runs(f"{d}/runs.tsv", "RUN01", "nanopore", f"{d}/sim/RUN01")
            quiet(sim.main, ["--reference", f"{d}/reference.tsv",
                             "--samplesheet", f"{d}/samplesheet.csv",
                             "--runs", f"{d}/runs.tsv", "--assay", ASSAY,
                             "--outdir", f"{d}/sim", "--reads-per-sample", "40",
                             "--error-model", f"nanopore_{variant}"])
            quiet(sc.main, ["--reads", f"{d}/sim/RUN01", "--assay", ASSAY,
                            "--outdir", f"{d}/split", "--quiet"])
            with open(f"{d}/split/split_summary.json") as fh:
                bc = json.load(fh).get("basecalling", {})
            model = bc.get("basecall_model", "")
            check(f"a run simulated as {variant} records itself as {variant}",
                  f"_{variant}@" in model, f"recorded {model!r}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_partial_second_unit() -> None:
    """A whole molecule followed by a fragment must still be recovered.

    Found on real data, not in a fixture: the splitter reported 1.21
    pieces per read where the census had counted 1.43 adapter structures.
    The gap was reads shaped [complete A][fragment of B carrying only its
    i7 end]. Those have no SECOND opening anchor, so the original rule -
    cut before each unit that opens - never fired. The read's outer ends
    are then A's i5 and B's i7, a hybrid pair, so the junction check
    refuses the whole read and molecule A, intact and assignable, is
    discarded along with the fragment.
    """
    print("\na whole molecule followed by a fragment")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, n_samples=24, unused=4, reads=120)
        cfg = dx.load_assay(ASSAY)
        a5 = cfg["i5_inner_anchor"]
        with gzip.open(f"{tmp}/sim/RUN01/RUN01.fastq.gz", "rt") as fh:
            L = fh.read().splitlines()
        recs = [(L[i][1:].split()[0], L[i + 1], L[i + 3])
                for i in range(0, len(L), 4)]

        os.makedirs(f"{tmp}/partial", exist_ok=True)
        made = 0
        with gzip.open(f"{tmp}/partial/p.fastq.gz", "wt") as fh:
            for k in range(0, len(recs) - 1, 2):
                a, b = recs[k], recs[k + 1]
                p = b[1].find(a5)
                if p < 0:
                    continue
                # Remove B's i5 end outright, leaving only its i7 end.
                cut = p + len(a5) + 40
                frag, fq = b[1][cut:], b[2][cut:]
                if len(frag) < 120:
                    continue
                fh.write(f"@P{k}\n{a[1] + frag}\n+\n{a[2] + fq}\n")
                made += 1
        check("the fixture actually built reads of this shape", made > 50,
              f"only {made}")

        for tag, reads in (("uncut", f"{tmp}/partial"), ("split", None)):
            if tag == "split":
                quiet(sc.main, ["--reads", f"{tmp}/partial", "--assay", ASSAY,
                                "--outdir", f"{tmp}/sp", "--quiet"])
                reads = f"{tmp}/sp/reads"
            quiet(dx.main, ["--reads", reads,
                            "--samplesheet", f"{tmp}/samplesheet.csv",
                            "--run-id", "RUN01", "--assay", ASSAY,
                            "--outdir", f"{tmp}/dx_{tag}"])
        with open(f"{tmp}/dx_uncut/demux_summary.json") as fh:
            uncut = json.load(fh)
        with open(f"{tmp}/dx_split/demux_summary.json") as fh:
            split = json.load(fh)

        check("uncut, these reads are refused rather than misassigned",
              uncut["by_tier"].get("E", 0) > 0.3 * uncut["total_reads"],
              f"tier E {uncut['by_tier'].get('E', 0)} of {uncut['total_reads']}")
        check("splitting recovers the complete molecule from most of them",
              split["assigned_reads"] > 10 * max(1, uncut["assigned_reads"]),
              f"{uncut['assigned_reads']} -> {split['assigned_reads']}")
        check("and the leftover fragments are not assigned to anything",
              split["by_tier"].get("B", 0) > 0,
              "no tier B - the fragments went somewhere unexpected")
        with open(f"{tmp}/sp/split_summary.json") as fh:
            summ = json.load(fh)
        check("no bases are invented or lost by the extra cut",
              summ["bases_retained_fraction"] > 0.99,
              f"{100 * summ['bases_retained_fraction']:.2f}%")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_three_units() -> None:
    """More than two molecules in a read, in mixed orientations.

    The cut rule walks adjacent pairs of anchors, so N units should give
    N-1 cuts without anything special. What needs checking is not the
    piece count but where the pieces go: a mis-placed cut in the middle
    of a three-unit read would produce a chimeric piece that still
    demultiplexes cleanly, which is the failure that cannot be noticed
    downstream.
    """
    print("\nthree molecules in one read")
    tmp = tempfile.mkdtemp()
    try:
        build(tmp, n_samples=24, unused=4, reads=90)
        with open(f"{tmp}/d0/assignments.tsv", encoding="utf-8") as fh:
            truth = {r["read_id"]: r["sample_id"]
                     for r in csv.DictReader(fh, delimiter="\t")
                     if r["sample_id"]}
        with gzip.open(f"{tmp}/sim/RUN01/RUN01.fastq.gz", "rt") as fh:
            L = fh.read().splitlines()
        recs = [(L[i][1:].split()[0], L[i + 1], L[i + 3])
                for i in range(0, len(L), 4)]
        good = [r for r in recs if r[0] in truth]

        rc, rng, prov = dx.revcomp, random.Random(11), {}
        os.makedirs(f"{tmp}/c3", exist_ok=True)
        with gzip.open(f"{tmp}/c3/c.fastq.gz", "wt") as fh:
            for k in range(0, len(good) - 3, 3):
                parts = good[k:k + 3]
                seq = qual = ""
                for p in parts:
                    inv = rng.random() < 0.4      # mixed orientations
                    seq += rc(p[1]) if inv else p[1]
                    qual += p[2][::-1] if inv else p[2]
                fh.write(f"@T{k}\n{seq}\n+\n{qual}\n")
                prov[f"T{k}"] = [truth[p[0]] for p in parts]

        quiet(sc.main, ["--reads", f"{tmp}/c3", "--assay", ASSAY,
                        "--outdir", f"{tmp}/sp", "--quiet"])
        quiet(dx.main, ["--reads", f"{tmp}/sp/reads",
                        "--samplesheet", f"{tmp}/samplesheet.csv",
                        "--run-id", "RUN01", "--assay", ASSAY,
                        "--outdir", f"{tmp}/dx"])
        got = calls_of(f"{tmp}/dx/assignments.tsv")

        tot, recovered, wrong, whole = len(prov), 0, 0, 0
        for pid, want in prov.items():
            calls = [got[k] for k in [f"{pid}_u{j}" for j in range(1, 6)] + [pid]
                     if got.get(k)]
            recovered += len(calls)
            if sorted(calls) == sorted(want):
                whole += 1
            elif not set(calls) <= set(want):
                wrong += 1
        check("three units give roughly three pieces",
              recovered > 2.5 * tot, f"{recovered} units from {tot} reads")
        check("most reads give up all three to the right samples",
              whole > 0.80 * tot, f"{whole} of {tot} ({100*whole/tot:.1f}%)")
        # The one that matters. Losing a unit costs yield; mislabelling
        # one costs correctness, and nothing downstream could detect it.
        check("and no unit is ever assigned to a sample it did not come from",
              wrong == 0, f"{wrong} of {tot} reads had a misassigned unit")

        with open(f"{tmp}/sp/split_summary.json") as fh:
            summ = json.load(fh)
        check("two cuts invent no sequence and lose none",
              summ["bases_retained_fraction"] > 0.99,
              f"{100 * summ['bases_retained_fraction']:.3f}%")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def chunk_reads(srcs, outdir: str, n: int) -> int:
    """Pool the simulated FASTQ and cut them into n MinKNOW-style batches."""
    os.makedirs(outdir, exist_ok=True)
    lines = []
    for src in srcs:
        with gzip.open(src, "rt") as fh:
            lines += fh.read().splitlines()
    recs = [lines[i:i + 4] for i in range(0, len(lines), 4)]
    per = (len(recs) + n - 1) // n
    written = 0
    for k in range(n):
        part = recs[k * per:(k + 1) * per]
        if not part:
            continue
        with gzip.open(os.path.join(outdir, f"FAX_pass_{k}.fastq.gz"), "wt") as fh:
            for r in part:
                fh.write("\n".join(r) + "\n")
        written += 1
    return written


def test_sharded_matches_whole_run() -> None:
    """Splitting file-by-file must equal splitting everything at once.

    Nextflow fans this stage out - one task per input FASTQ, one core
    each - because it is the step that hit a 12 hour wall. That is only
    safe if the fanned-out route produces the same thing as the single
    task did, and "the same" has to be checked rather than reasoned
    about: the two routes share the cutting code but not the plumbing
    around it, which is where output naming and summary arithmetic
    live.

    A first version of this test compared two sets of empty files and
    reported that they matched, because the fixture FASTQ it used were
    truncated. Hence the non-vacuity checks below: the comparison is
    worthless unless reads went in, pieces came out, and some read was
    actually cut.
    """
    print("\nsharded splitting equals whole-run splitting")
    tmp = tempfile.mkdtemp()
    try:
        quiet(build, tmp, 24, 4, 120)
        # The concatemers live in a separate fixture file, and without
        # them this compares two runs of a no-op: the splitter would
        # pass every read through untouched and the routes would agree
        # for the wrong reason.
        inp = f"{tmp}/in"
        n_files = chunk_reads([f"{tmp}/sim/RUN01/RUN01.fastq.gz",
                               f"{tmp}/chim/c.fastq.gz"], inp, 4)
        check("the input was chunked into several files", n_files > 1,
              str(n_files))

        # -- route A: one task splits the lot --------------------------
        whole = f"{tmp}/whole"
        quiet(sc.main, ["--reads", inp, "--assay", ASSAY,
                        "--threads", "1", "--outdir", whole])
        sum_a = json.load(open(f"{whole}/split_summary.json"))

        # The comparison means nothing if nothing happened.
        check("reads actually went through", sum_a["reads_in"] > 0,
              str(sum_a["reads_in"]))
        check("pieces actually came out", sum_a["pieces_out"] > 0,
              str(sum_a["pieces_out"]))
        check("and some read was actually cut", sum_a["reads_split"] > 0,
              str(sum_a["reads_split"]))

        # -- route B: one task per file, then merge --------------------
        shard = f"{tmp}/shard"
        quiet(sc.main, ["--reads", inp, "--assay", ASSAY,
                        "--list-only", "--outdir", shard])
        listed = [l.split("\t") for l in
                  open(f"{shard}/read_files.txt").read().splitlines() if l]
        check("the listing carries an index per file",
              all(len(r) == 2 for r in listed) and len(listed) == n_files,
              str(listed[:2]))

        for idx, path in listed:
            quiet(sc.main, ["--reads", path, "--assay", ASSAY,
                            "--part-index", idx, "--threads", "1",
                            "--outdir", f"{shard}/s{idx}"])
        quiet(ms.main, ["--parts"] + sorted(glob.glob(f"{shard}/s*/parts"))
              + ["--file-report", f"{shard}/file_report.json",
                 "--expect-parts", str(n_files),
                 "--outdir", shard])
        sum_b = json.load(open(f"{shard}/split_summary.json"))

        # -- the reads ------------------------------------------------
        got_a = {os.path.basename(f): f
                 for f in glob.glob(f"{whole}/reads/*.fastq.gz")}
        got_b = {os.path.basename(f): f
                 for f in glob.glob(f"{shard}/s*/reads/*.fastq.gz")}
        check("the shards produce the same output filenames",
              set(got_a) == set(got_b),
              f"whole={sorted(got_a)} shard={sorted(got_b)}")

        same, diff = 0, []
        for name, pa in sorted(got_a.items()):
            pb = got_b.get(name)
            if not pb:
                diff.append(f"{name} missing")
                continue
            a = gzip.open(pa, "rt").read()
            b = gzip.open(pb, "rt").read()
            if a == b and a:
                same += 1
            else:
                diff.append(f"{name} differs ({len(a)} vs {len(b)} chars)")
        check("every shard's reads are identical to the whole-run reads",
              not diff and same == len(got_a), "; ".join(diff[:3]))
        check("and those comparisons were not of empty files", same > 0,
              str(same))

        # -- the summary ----------------------------------------------
        keys = sorted(set(sum_a) | set(sum_b))
        bad = [k for k in keys if sum_a.get(k) != sum_b.get(k)]
        check("the merged summary equals the whole-run summary", not bad,
              "; ".join(f"{k}: {sum_a.get(k)!r} vs {sum_b.get(k)!r}"[:70]
                        for k in bad[:3]))

        # -- the guards on the merge ----------------------------------
        try:
            ms.main(["--parts"] + sorted(glob.glob(f"{shard}/s*/parts"))[:-1]
                    + ["--expect-parts", str(n_files),
                       "--outdir", f"{tmp}/short", "--quiet"])
            check("a missing shard is refused", False, "no error raised")
        except SystemExit as exc:
            check("a missing shard is refused", "expected" in str(exc),
                  str(exc)[:70])

        parts = sorted(glob.glob(f"{shard}/s*/parts"))
        try:
            ms.main(["--parts"] + parts + [parts[0]]
                    + ["--outdir", f"{tmp}/dupe", "--quiet"])
            check("a duplicated shard is refused", False, "no error raised")
        except SystemExit as exc:
            check("a duplicated shard is refused",
                  "more than one shard" in str(exc), str(exc)[:70])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_cutting_preserves_indices() -> None:
    """A cut must not cost a molecule its index.

    The suspicion, from a real run where 51.7% of pieces reached demux
    with only one end: perhaps the cut itself was landing inside the
    index and manufacturing those single-ended pieces.

    It is not. The cut is placed at a unit's own boundary - just before
    an opening unit's outer tail, or just past a closing unit's - so on
    a clean join both pieces come back byte-identical to the molecules
    that were ligated. This checks that for every orientation, because
    inverted joins were 35.9% of the real run and the reverse-strand
    anchors are a separate code path.

    The two truncated cases are the point of the exercise: there the
    single-ended piece is single-ended because that molecule really did
    lose its index before ligation, and the OTHER molecule - intact and
    assignable - is recovered rather than thrown away with it.
    """
    print("\ncutting preserves both indices")
    cfg = sc.load_assay(ASSAY)
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    reps = 2
    span = 8 * reps
    head_len = {"5+": span + len(cfg["p5_tail"]), "7-": span + len(cfg["p7_tail"])}
    tail_len = {"7+": span + len(cfg["p7_tail"]), "5-": span + len(cfg["p5_tail"])}

    amp = "ACGT" * 60
    m1 = sim.build_molecule(amp, "GGGGTTTT", "AAAACCCC", reps, cfg)
    m2 = sim.build_molecule(amp, "CCCCAAAA", "TTTTGGGG", reps, cfg)
    rc = dx.revcomp
    whole = (m1, m2, rc(m1), rc(m2))

    def split(seq):
        return sc.split_read(seq, "I" * len(seq), a5, a7, 2, head_len, tail_len)

    def ends(s):
        f5 = a5 in s or rc(a5) in s
        f7 = a7 in s or rc(a7) in s
        return ("both" if f5 and f7 else
                "i5" if f5 else "i7" if f7 else "none")

    # A single molecule must come through untouched.
    one = split(m1)
    check("a lone molecule is not cut", len(one) == 1 and one[0][0] == m1)

    # Every orientation of a two-molecule join: two intact molecules.
    for name, seq in (("head-to-tail", m1 + m2),
                      ("tail-to-tail", m1 + rc(m2)),
                      ("head-to-head", rc(m1) + m2)):
        got = split(seq)
        check(f"{name}: cut into two", len(got) == 2, str(len(got)))
        check(f"{name}: both pieces keep both indices",
              all(ends(p) == "both" for p, _ in got),
              str([ends(p) for p, _ in got]))
        check(f"{name}: both pieces are the original molecules",
              all(p in whole for p, _ in got),
              str([len(p) for p, _ in got]))

    # Three fragments, which the first version of the cut logic mishandled.
    got = split(m1 + m2 + m1)
    check("a three-molecule concatemer gives three pieces", len(got) == 3,
          str(len(got)))
    check("and all three keep both indices",
          all(ends(p) == "both" for p, _ in got),
          str([ends(p) for p, _ in got]))

    # A truncated partner must not take the intact molecule down with it.
    got = split(m1[:-60] + m2)
    check("when the FIRST molecule lost its i7, the second is recovered whole",
          len(got) == 2 and got[1][0] == m2,
          str([ends(p) for p, _ in got]))
    check("and the fragment is honestly reported as single-ended",
          ends(got[0][0]) == "i5", ends(got[0][0]))

    got = split(m1 + m2[60:])
    check("when the SECOND molecule lost its i5, the first is recovered whole",
          len(got) == 2 and got[0][0] == m1,
          str([ends(p) for p, _ in got]))
    check("and that fragment is single-ended too",
          ends(got[1][0]) == "i7", ends(got[1][0]))


def _geometry():
    cfg = sc.load_assay(ASSAY)
    span = 8 * 2
    head_len = {"5+": span + len(cfg["p5_tail"]), "7-": span + len(cfg["p7_tail"])}
    tail_len = {"7+": span + len(cfg["p7_tail"]), "5-": span + len(cfg["p5_tail"])}
    return cfg, head_len, tail_len


def test_earlier_worse_anchor_is_not_skipped() -> None:
    """Every anchor within tolerance is found, not only the best one.

    The scan used to ask "where does this anchor fit best in the rest of
    the read" and resume past the answer. A read of [molecule A whose i7
    anchor carries one error][a fragment of B ending in a clean i7
    anchor] answers with B's anchor, and A's is skipped: one CLOSE where
    there are two, no CLOSE -> CLOSE cut, and the read reaches demux
    whole - A's i5 with B's i7.
    """
    print("\nan imperfect anchor before a perfect one")
    from ampseqlib import find_anchor_cost, find_anchor_leftmost
    cfg, head_len, tail_len = _geometry()
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]

    # The primitive first. One substitution, then an exact copy.
    bad = a7[:5] + ("A" if a7[5] != "A" else "C") + a7[6:]
    text = "T" * 30 + bad + "T" * 30 + a7 + "T" * 10
    lm = find_anchor_leftmost(text, a7, 2)
    check("the leftmost search finds the first, imperfect copy",
          lm is not None and lm[1] == 30 + len(a7) and lm[2] == 1, str(lm))
    check("and lands on its end, not a ragged prefix of an exact copy",
          find_anchor_leftmost("G" * 20 + a7 + "G" * 20, a7, 2)
          == (20, 20 + len(a7), 0),
          str(find_anchor_leftmost("G" * 20 + a7 + "G" * 20, a7, 2)))
    check("find_anchor_cost still answers with the best copy, for its "
          "other callers", find_anchor_cost(text, a7, 2)[2] == 0)
    check("both copies are scanned", len(sc.scan_all(text, a7, 2)) == 2,
          str(sc.scan_all(text, a7, 2)))

    # Then the read that went wrong.
    m1 = sim.build_molecule("ACGT" * 60, "GGGGTTTT", "AAAACCCC", 2, cfg)
    m2 = sim.build_molecule("ACGT" * 60, "CCCCAAAA", "TTTTGGGG", 2, cfg)
    p = m1.find(a7)
    check("the fixture molecule carries its i7 anchor forward", p > 0)
    m1_bad = m1[:p + 5] + ("A" if m1[p + 5] != "A" else "C") + m1[p + 6:]
    frag = m2[m2.find(a5) + len(a5) + 40:]
    seq = m1_bad + frag
    got = sc.split_read(seq, "I" * len(seq), a5, a7, 2, head_len, tail_len)
    check("a one-error anchor before an exact one still triggers the cut",
          len(got) == 2 and got[0][0] == m1_bad,
          str([len(s) for s, _ in got]) + f" vs A={len(m1_bad)}")


def test_split_accounting() -> None:
    """Every read is accounted for, and "split" means split in the output.

    Two gaps. A read all of whose pieces fell below the length floor
    vanished with no counter, so reads_in could not be reconciled with
    anything written. And reads_split and the pieces-per-read histogram
    were taken BEFORE the floor, so a read whose cut left a remainder of
    a few bases counted as split while producing one piece.
    """
    print("\nsplit read accounting")
    cfg, head_len, tail_len = _geometry()
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    m1 = sim.build_molecule("ACGT" * 60, "GGGGTTTT", "AAAACCCC", 2, cfg)
    m2 = sim.build_molecule("ACGT" * 60, "CCCCAAAA", "TTTTGGGG", 2, cfg)
    stub = m2[-80:]          # B's i7 end only, below the 100bp floor
    check("the fixture's fragment carries an anchor and is short",
          a7 in stub and len(stub) < 100)
    reads = {
        "whole": m1,                    # one piece, kept
        "joined": m1 + m2,              # two pieces, both kept
        "stub_tail": m1 + stub,         # cut, but only one piece survives
        "tiny": "ACGT" * 10,            # nothing survives at all
    }
    tmp = tempfile.mkdtemp()
    try:
        src = f"{tmp}/in.fastq"
        with open(src, "w") as fh:
            for rid, s in reads.items():
                fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
        job = {"src": src, "dest": f"{tmp}/out.fastq.gz",
               "stats": f"{tmp}/part.json", "a5": a5, "a7": a7, "edits": 2,
               "head_len": head_len, "tail_len": tail_len,
               "min_length": 100, "code": "x", "index_length": 8,
               "barcode_repeats": 2}
        st = sc.split_one_file(job)
        check("a read that left no piece is counted",
              st.get("reads_with_no_piece_kept") == 1,
              str(st.get("reads_with_no_piece_kept")))
        check("reads_split counts only reads that came out in pieces",
              st["reads_split"] == 1, str(st["reads_split"]))
        check("reads the cut rule fired on are still visible",
              st.get("reads_cut") == 2, str(st.get("reads_cut")))
        dist = {int(k): v for k, v in st["pieces_per_read_distribution"].items()}
        check("the histogram is of kept pieces and sums to reads_in",
              dist == {0: 1, 1: 2, 2: 1}
              and sum(dist.values()) == st["reads_in"], str(dist))
        check("and pieces_out is its weighted sum",
              sum(k * v for k, v in dist.items()) == st["pieces_out"])

        # Merged across shards, the new counters are summed like the old.
        # A part written without them (the workflow's stub) contributes
        # nothing rather than failing the merge.
        stub_part = {"source": "stub", "reads_in": 0, "pieces_out": 0,
                     "reads_split": 0, "pieces_dropped_below_min_length": 0,
                     "pieces_per_read_distribution": {}, "bases_in": 0,
                     "bases_out": 0}
        st2 = dict(st, source="other")
        os.makedirs(f"{tmp}/parts")
        for name, part in (("a", st), ("b", st2), ("c", stub_part)):
            with open(f"{tmp}/parts/{name}.json", "w") as fh:
                json.dump(part, fh)
        quiet(ms.main, ["--parts", f"{tmp}/parts", "--outdir", f"{tmp}/m",
                        "--expect-parts", "3", "--quiet"])
        summ = json.load(open(f"{tmp}/m/split_summary.json"))
        check("merge sums reads with no piece across shards",
              summ["reads_with_no_piece_kept"] == 2
              and summ["reads_cut"] == 4 and summ["reads_split"] == 2,
              str({k: summ.get(k) for k in ("reads_with_no_piece_kept",
                                            "reads_cut", "reads_split")}))
        check("and the merged histogram still reconciles with reads_in",
              sum(summ["pieces_per_read_distribution"].values())
              == summ["reads_in"] == 8)
        check("the cut geometry in force is reported",
              summ["settings"].get("index_length") == 8,
              str(summ["settings"]))
        try:
            quiet(ms.main, ["--parts", f"{tmp}/parts", "--outdir", f"{tmp}/m2",
                            "--expect-parts", "4", "--quiet"])
            refused = False
        except SystemExit:
            refused = True
        check("--expect-parts still refuses a missing shard", refused)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_panel_lookalikes_are_not_junctions() -> None:
    """An anchor spelled out inside an amplicon is insert, not an adapter.

    Found on a real panel: a Trichuris incognita amplicon carries the i5
    anchor at two edits, 353bp in (the context below is its own). The
    splitter cut every one of its reads there and the chimera check
    quarantined any that got through, so the species was under-counted
    with no error anywhere. Pre-flight now lists such look-alikes and both
    stages leave them alone - and only them: a real junction, including
    one beside a look-alike, is still cut and still quarantined.
    """
    print("\npanel anchor look-alikes are not junctions")
    import ampseqlib as al
    cfg = sc.load_assay(ASSAY)
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    reps, span = 2, 16
    head_len = {"5+": span + len(cfg["p5_tail"]), "7-": span + len(cfg["p7_tail"])}
    tail_len = {"7+": span + len(cfg["p7_tail"]), "5-": span + len(cfg["p5_tail"])}
    rng = random.Random(353)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    site = "GCCCCCTAGATCGACGGCAGCGGCAGTCGACTAG"   # T. incognita, as found
    insert = rand(320) + site + rand(120)
    other = rand(420)
    lk = al.find_anchor_lookalikes(
        [("T_incognita", insert), ("clean", other)],
        {"i5_inner_anchor": a5, "i7_inner_anchor": a7},
        al.LOOKALIKE_MAX_EDITS)
    check("the look-alike in the insert is found",
          any(r["seq_id"] == "T_incognita" and r["anchor"] == "i5_inner_anchor"
              and r["edits"] <= 2 for r in lk), str(lk))
    mask = al.LookalikeMask(r["context"] for r in lk)

    mol = sim.build_molecule(insert, "GGGGTTTT", "AAAACCCC", reps, cfg)
    clean = sim.build_molecule(other, "CCCCAAAA", "TTTTGGGG", reps, cfg)
    rc = dx.revcomp

    def split(seq, m=None):
        return sc.split_read(seq, "I" * len(seq), a5, a7, 2, head_len,
                             tail_len, m)

    check("without the mask the amplicon is cut in two (the bug)",
          len(split(mol)) == 2, str(len(split(mol))))
    check("with the mask it comes through whole",
          len(split(mol, mask)) == 1 and split(mol, mask)[0][0] == mol)
    check("and in the other orientation too",
          len(split(rc(mol), mask)) == 1)

    # A sequencing error or two near the site must not undo the mask.
    noisy = list(mol)
    k = mol.find(site)
    for j in (k + 3, k + 30):
        noisy[j] = "A" if noisy[j] != "A" else "C"
    check("the mask survives sequencing errors near the site",
          len(split("".join(noisy), mask)) == 1)

    # Real junctions are untouched, including one beside a look-alike.
    for name, seq in (("a clean concatemer", clean + clean),
                      ("a junction beside the look-alike", mol + clean)):
        got = split(seq, mask)
        check(f"{name} is still cut", len(got) == 2, str(len(got)))

    # The chimera check: the same distinction, one stage later.
    junction = (a5, a7, rc(a5), rc(a7))
    ins_start = len(cfg["p5_tail"]) + span + len(cfg["pcr1_overhang_fwd"])
    ins_end = len(mol) - (len(cfg["pcr1_overhang_rev"]) + span
                          + len(cfg["p7_tail"]))
    check("without the mask the chimera check quarantines it (the bug)",
          dx.has_internal_adapter(mol, ins_start, ins_end, junction, 2))
    check("with the mask it does not",
          not dx.has_internal_adapter(mol, ins_start, ins_end, junction, 2,
                                      mask))
    both = mol + clean
    check("a real internal junction is still caught, even after a "
          "masked look-alike",
          dx.has_internal_adapter(both, ins_start, len(both) - 40,
                                  junction, 2, mask))

    # End to end through the demultiplexer: assigned, not set aside.
    # build_molecule takes (amplicon, i7, i5, ...): mol carries i7
    # GGGGTTTT and i5 AAAACCCC.
    samples = [{"run_id": "R", "sample_id": "S1", "index_i5": "AAAACCCC",
                "index_i7": "GGGGTTTT"},
               {"run_id": "R", "sample_id": "S2", "index_i5": "TTTTGGGG",
                "index_i7": "CCCCAAAA"}]
    s5 = dx.BarcodeSet.build({x["index_i5"] for x in samples}, 1)
    s7 = dx.BarcodeSet.build({x["index_i7"] for x in samples}, 1)
    pairs = {(x["index_i5"], x["index_i7"]): x["sample_id"] for x in samples}
    call = dx.assign_read("r", mol, cfg, s5, s7, pairs, reps, 120, False,
                          junction_anchors=junction, lookalikes=mask)
    check("the read is assigned to its sample in tier A",
          call.tier == "A" and call.sample_id == "S1",
          f"{call.tier} {call.sample_id}")
    call0 = dx.assign_read("r", mol, cfg, s5, s7, pairs, reps, 120, False,
                           junction_anchors=junction)
    check("where without the mask it was quarantined as tier E",
          call0.tier == "E", call0.tier)


REAL_READS = ("/Users/sd21/Documents/WORK/informatics_workbook/"
              "amplicon_sequencing/nextflow_testing/run/fastq_pass")
REAL_LOOKALIKES = ("/Users/sd21/Documents/WORK/informatics_workbook/"
                   "amplicon_sequencing/nextflow_testing/results-field/"
                   "preflight/anchor_lookalikes.tsv")


def _reference_scan_hits(seq, a5, a7, edits, mask=None):
    """The splitter's anchor scan as it was before seed-and-verify, kept
    verbatim as the standard the fast scan is held to: find_anchor_leftmost
    across every remainder of the read, and a mask that aligns every
    look-alike context with find_anchor_cost."""
    from ampseqlib import find_anchor_cost, find_anchor_leftmost

    def scan_all(text, pattern, max_edits, limit=24):
        hits, off = [], 0
        while len(hits) < limit and off < len(text):
            h = find_anchor_leftmost(text[off:], pattern, max_edits)
            if h is None:
                break
            s, e, _c = h
            hits.append((off + s, off + e))
            off += max(e, s + 1)
        return hits

    def covers(read, start, end):
        for ctx in mask.contexts:
            lo = max(0, start - len(ctx))
            hi = min(len(read), end + len(ctx))
            h = find_anchor_cost(read[lo:hi], ctx, max(4, len(ctx) // 8))
            if h is not None and lo + h[0] < end and lo + h[1] > start:
                return True
        return False

    hits = []
    for kind, pat in (("5+", a5), ("7+", a7),
                      ("5-", dx.revcomp(a5)), ("7-", dx.revcomp(a7))):
        for s, e in scan_all(seq, pat, edits):
            if mask and covers(seq, s, e):
                continue
            hits.append((s, e, kind))
    hits.sort()
    kept = []
    for h in hits:
        if kept and h[0] < kept[-1][1]:
            continue
        kept.append(h)
    return kept


def _real_reads(per_file: int) -> list[str]:
    seqs = []
    for path in sorted(glob.glob(os.path.join(REAL_READS, "*.fastq.gz"))):
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seqs.append(line.strip())
                if i >= 4 * per_file:
                    break
    return seqs


def test_seeded_scan_matches_reference() -> None:
    """Seed-and-verify finds exactly the anchors the full alignment did.

    Aligning every base of every read four times over, and every
    look-alike context at every hit, was two-thirds of a laptop run. The
    splitter now aligns only near exact pieces of the pattern, which by
    the pigeonhole argument cannot miss a copy within tolerance - but a
    cut placed one base differently is a different output, so the claim
    is checked hit for hit against the old scan: same (start, end, kind)
    list, same order, with and without the look-alike mask.
    """
    print("\nseeded anchor scan matches the full alignment")
    import time
    import ampseqlib as al
    cfg = sc.load_assay(ASSAY)
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    rc = dx.revcomp
    rng = random.Random(1428)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def idx():
        return rand(8)

    def mutate(seq, n, near=None):
        # Substitutions, insertions and deletions, optionally aimed at a
        # stretch (an anchor) so that tolerance edges are exercised.
        s = list(seq)
        for _ in range(n):
            if not s:
                break
            if near is not None and near[1] > near[0]:
                j = rng.randrange(max(0, near[0] - 2),
                                  min(len(s), near[1] + 2))
            else:
                j = rng.randrange(len(s))
            op = rng.random()
            if op < 0.4:
                s[j] = rng.choice("ACGT")
            elif op < 0.7:
                s.insert(j, rng.choice("ACGT"))
            else:
                del s[j]
        return "".join(s)

    # The primitive first, against the matchers it stands in for, across
    # tolerances and on texts dense with near-copies.
    bad = []
    for _ in range(3000):
        k = rng.choice((0, 1, 2, 3))
        pat = rng.choice((a5, a7, rc(a5), rc(a7), rand(rng.randint(4, 20))))
        text = "".join(rng.choice((rand(rng.randint(0, 12)),
                                   mutate(pat, rng.randint(0, 3))))
                       for _ in range(rng.randint(0, 6)))
        search = al.SeededAnchorScan(pat, k)
        hits, off = [], 0
        while off < len(text):
            h = al.find_anchor_leftmost(text[off:], pat, k)
            if h is None:
                break
            hits.append((off + h[0], off + h[1]))
            off += max(h[1], h[0] + 1)
        if search.scan(text, 10 ** 6) != hits:
            bad.append(("scan", pat, k, text))
        if search.best(text) != al.find_anchor_cost(text, pat, k):
            bad.append(("best", pat, k, text))
    check("SeededAnchorScan.scan and .best agree with find_anchor_leftmost "
          "and find_anchor_cost on 3000 near-copy texts",
          not bad, str(bad[:3]))

    # The real look-alike, and a mask built the way pre-flight builds it.
    site = "GCCCCCTAGATCGACGGCAGCGGCAGTCGACTAG"   # T. incognita, as found
    incog = rand(320) + site + rand(120)
    lk = al.find_anchor_lookalikes(
        [("T_incognita", incog)],
        {"i5_inner_anchor": a5, "i7_inner_anchor": a7},
        al.LOOKALIKE_MAX_EDITS)
    mask = al.LookalikeMask(r["context"] for r in lk)
    check("the T. incognita look-alike is in the mask", bool(mask))

    def unit():
        insert = incog if rng.random() < 0.3 else rand(rng.randint(150, 500))
        m = sim.build_molecule(insert, idx(), idx(), 2, cfg)
        return rc(m) if rng.random() < 0.5 else m

    reads = []
    for _ in range(1500):
        kind = rng.random()
        if kind < 0.45:          # concatemers, either way round
            seq = "".join(unit() for _ in range(rng.randint(1, 4)))
        elif kind < 0.6:         # inverted joins: a unit and its reverse
            u = unit()
            seq = u + rc(u)[rng.randint(0, 60):]
        elif kind < 0.7:         # a truncated fragment on either side
            u = unit()
            seq = u[rng.randint(0, len(u) // 2):] + unit() \
                + unit()[:rng.randint(0, 200)]
        elif kind < 0.8:         # short reads, down to nothing
            seq = unit()[:rng.randint(0, 40)] if rng.random() < 0.5 \
                else mutate(rng.choice((a5, a7, rc(a5), rc(a7))), 1)
        elif kind < 0.9:         # the look-alike itself, noisy
            seq = rand(rng.randint(0, 50)) + mutate(site, rng.randint(0, 4)) \
                + rand(rng.randint(0, 50))
        else:                    # an ambiguity code: the unseeded path
            u = unit()
            j = rng.randrange(len(u))
            seq = u[:j] + rng.choice("NRY") + u[j + 1:]
        # Errors everywhere, and more aimed at the anchors.
        seq = mutate(seq, rng.randint(0, len(seq) // 60 + 1))
        for a in (a5, a7, rc(a5), rc(a7)):
            p = seq.find(a[2:10])
            if p >= 0 and rng.random() < 0.6:
                seq = mutate(seq, rng.randint(1, 3), (p - 2, p + 12))
        reads.append(seq)

    diff = []
    for seq in reads:
        for m in (None, mask):
            if sc.scan_hits(seq, a5, a7, 2, m) != \
                    _reference_scan_hits(seq, a5, a7, 2, m):
                diff.append(seq)
    n_hits = sum(len(_reference_scan_hits(s, a5, a7, 2, mask)) for s in reads)
    check(f"{len(reads)} simulated reads ({n_hits} anchor hits): identical "
          "hits with and without the mask", not diff,
          f"{len(diff)} differ, e.g. {diff[:1]}")

    # Real reads, when the reference run is on this machine.
    if not os.path.isdir(REAL_READS) or not glob.glob(
            os.path.join(REAL_READS, "*.fastq.gz")):
        print("  skip real reads: not on this machine")
        return
    real = _real_reads(1000)
    rmask = al.LookalikeMask.load(REAL_LOOKALIKES) \
        if os.path.exists(REAL_LOOKALIKES) else mask
    t0 = time.time()
    want = [_reference_scan_hits(s, a5, a7, 2, rmask) for s in real]
    t1 = time.time()
    got = [sc.scan_hits(s, a5, a7, 2, rmask) for s in real]
    t2 = time.time()
    n_diff = sum(a != b for a, b in zip(want, got))
    check(f"{len(real)} real reads ({sum(map(len, want))} anchor hits): "
          "identical hits", n_diff == 0, f"{n_diff} differ")
    print(f"       full alignment {len(real) / (t1 - t0):.0f} reads/s, "
          f"seeded {len(real) / (t2 - t1):.0f} reads/s (one core)")


DB_RELEASE = ("/Users/sd21/Documents/WORK/informatics_workbook/"
              "vespa-18s-refdb/releases/vespa18s_v0.4.0")


def test_indexed_lookalike_mask() -> None:
    """A database-sized mask answers exactly as the plain loop does.

    A panel's mask holds a handful of contexts and covers() tries each in
    turn. A database release holds ~46,000 rows, where that loop costs
    0.3 s per anchor hit, so above LOOKALIKE_INDEX_MIN contexts covers()
    first asks an index which contexts could possibly match and tries
    only those. The index is meant to be exact - it leaves out only
    contexts with no alignment in tolerance - so the claim is checked
    call for call against the loop over every context (_covers_by with
    all of them), on hits that are look-alikes as well as hits that are
    not, so a filter that dropped real candidates would show.
    """
    print("\nindexed look-alike mask matches the plain loop")
    import time
    import ampseqlib as al
    cfg = sc.load_assay(ASSAY)
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    rc = dx.revcomp
    rng = random.Random(4242)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def mutate(s, n):
        s = list(s)
        for _ in range(n):
            if not s:
                break
            j = rng.randrange(len(s))
            op = rng.random()
            if op < 0.5:
                s[j] = rng.choice("ACGT")
            elif op < 0.75:
                del s[j]
            else:
                s.insert(j, rng.choice("ACGT"))
        return "".join(s)

    # 150 inserts, each carrying an anchor near-copy at 0-3 edits: well
    # over the index threshold once both strands are counted.
    inserts = []
    for k in range(150):
        anchor = rng.choice((a5, a7))
        inserts.append((f"ins{k}", rand(rng.randint(30, 200))
                        + mutate(anchor, rng.randint(0, 3))
                        + rand(rng.randint(30, 200))))
    lk = al.find_anchor_lookalikes(
        inserts, {"i5_inner_anchor": a5, "i7_inner_anchor": a7},
        al.LOOKALIKE_MAX_EDITS)
    mask = al.LookalikeMask(r["context"] for r in lk)
    check(f"synthetic mask is indexed ({len(mask.contexts)} contexts)",
          mask._index is not None)

    def every(read, s, e):
        return mask._covers_by(range(len(mask.contexts)), read, s, e,
                               al.is_plain_acgt(read[max(0, s - mask.widest):
                                                     e + mask.widest]))

    def idx():
        return rand(8)

    calls = diff = positive = 0
    for _ in range(600):
        name, ins = rng.choice(inserts)
        unit = sim.build_molecule(mutate(ins, rng.randint(0, 12)), idx(),
                                  idx(), 2, cfg)
        read = rc(unit) if rng.random() < 0.5 else unit
        if rng.random() < 0.3:
            read = read + rc(sim.build_molecule(rand(300), idx(), idx(), 2,
                                                cfg))
        if rng.random() < 0.05:
            j = rng.randrange(len(read))
            read = read[:j] + "N" + read[j + 1:]
        spots = [(s, e) for pat in (a5, a7, rc(a5), rc(a7))
                 for s, e in sc.scan_all(read, pat, 3)]
        # and intervals that are not anchor hits at all
        spots += [(j, j + rng.randint(1, 20))
                  for j in rng.sample(range(len(read)), 3)]
        for s, e in spots:
            want = every(read, s, e)
            calls += 1
            positive += want
            if mask.covers(read, s, e) != want:
                diff += 1
    check(f"{calls} calls on noisy concatemers ({positive} covered): "
          "indexed covers() == plain loop", diff == 0 and positive > 50,
          f"{diff} differ, {positive} positive")

    # The real release, when it is on this machine: its own mask, real
    # adapter hits, and its own look-alikes inside noisy reads.
    path = os.path.join(DB_RELEASE, "anchor_lookalikes.tsv")
    if not os.path.exists(path):
        print("  skip database release: not on this machine")
        return
    t0 = time.time()
    big = al.LookalikeMask.load(path)
    print(f"       release mask: {len(big.contexts)} contexts, "
          f"built in {time.time() - t0:.1f} s")
    rows = al.read_table(path)
    rng.shuffle(rows)
    reads = []
    for r in rows[:200]:
        ctx = r["context"]
        core = mutate(ctx, rng.randint(0, 6))
        unit = sim.build_molecule(rand(rng.randint(20, 200)) + core
                                  + rand(rng.randint(20, 200)),
                                  idx(), idx(), 2, cfg)
        reads.append(rc(unit) if rng.random() < 0.5 else unit)
    # VESPA-shaped reads: release inserts behind the primer, where the
    # 18S next to the adapter genuinely matches many contexts.
    fa = os.path.join(DB_RELEASE, "refdb.fasta")
    with open(fa) as fh:
        inserts = [ln.strip() for i, ln in enumerate(fh) if i % 194 == 1]
    for ins in inserts[:150]:
        unit = sim.build_molecule("AGCAGCCGCGGTAATTCC" + ins
                                  + "GAAACTTAAAGGAATTGACGGA", idx(), idx(),
                                  2, cfg)
        reads.append(mutate(rc(unit) if rng.random() < 0.5 else unit,
                            len(unit) // 40))
    if os.path.isdir(REAL_READS):
        reads += _real_reads(100)
    calls = diff = positive = 0
    t_idx = t_all = 0.0
    for read in reads:
        for pat in (a5, a7, rc(a5), rc(a7)):
            for s, e in sc.scan_all(read, pat, 2):
                plain = al.is_plain_acgt(read[max(0, s - big.widest):
                                              e + big.widest])
                t = time.time()
                got = big.covers(read, s, e)
                t_idx += time.time() - t
                if calls < 60:
                    t = time.time()
                    want = big._covers_by(range(len(big.contexts)), read,
                                          s, e, plain)
                    t_all += time.time() - t
                else:
                    want = None
                calls += 1
                positive += got
                if want is not None and got != want:
                    diff += 1
    # The plain loop is too slow to run on every call; the first 60 are
    # compared outright, and below, for a dozen reads, every context is
    # tried on its own: any whose best alignment overlaps the hit (the
    # one thing covers() asks) must be among the index's candidates.
    miss = 0
    for read in reads[:6] + reads[-6:]:
        for pat in (a5, a7, rc(a5), rc(a7)):
            for s, e in sc.scan_all(read, pat, 2):
                ws = max(0, s - big.widest)
                text = read[ws:e + big.widest]
                plain = al.is_plain_acgt(text)
                cand = set(big._index.candidates(text, s - ws, e - ws,
                                                 big.widest, plain))
                for i in range(len(big.contexts)):
                    if i not in cand and big._covers_by([i], read, s, e,
                                                        plain):
                        miss += 1
    check(f"release mask: {calls} anchor hits ({positive} covered), first "
          "60 identical to the plain loop", diff == 0 and positive > 0,
          f"{diff} differ")
    check("release mask: no context whose best alignment overlaps the hit "
          "was left out of the candidates (12 reads, every context)",
          miss == 0, f"{miss} missed")
    print(f"       per call: indexed {1000 * t_idx / max(1, calls):.2f} ms, "
          f"plain loop {1000 * t_all / 60:.0f} ms")


def test_indexed_mask_matches_reference() -> None:
    """The indexed look-alike mask answers exactly as trying every context.

    A reference database's mask holds tens of thousands of contexts, and
    covers() used to try each one at every anchor hit. It now tries only
    those a seed-piece index names and the q-gram lemma does not rule out
    - both bounds that cannot drop a context within tolerance. Checked
    here against the definition itself, find_anchor_cost on each context
    in turn, over junctions, embedded look-alikes with edits, contexts
    carrying IUPAC codes (as database sequences do), and reads with N.
    """
    print("\nindexed look-alike mask matches trying every context")
    import ampseqlib as al
    cfg = sc.load_assay(ASSAY)
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    rng = random.Random(2071)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def mutate(seq, n):
        s = list(seq)
        for _ in range(n):
            j = rng.randrange(len(s))
            op = rng.random()
            if op < 0.5:
                s[j] = rng.choice("ACGT")
            elif op < 0.75:
                s.insert(j, rng.choice("ACGT"))
            else:
                del s[j]
        return "".join(s)

    # Contexts as find_anchor_lookalikes writes them: insert, a near-copy
    # of an anchor, insert. Some carry IUPAC codes, as database ones do.
    contexts = []
    for _ in range(400):
        mid = mutate(rng.choice((a5, a7)), rng.randint(0, 3))
        ctx = rand(24) + mid + rand(24)
        if rng.random() < 0.15:
            ctx = "".join(rng.choice("NRYS") if rng.random() < 0.03 else c
                          for c in ctx)
        contexts.append(ctx)
    mask = al.LookalikeMask(contexts)

    def reference(read, start, end):
        for ctx in mask.contexts:
            lo = max(0, start - len(ctx))
            hi = min(len(read), end + len(ctx))
            h = al.find_anchor_cost(read[lo:hi], ctx, max(4, len(ctx) // 8))
            if h is not None and lo + h[0] < end and lo + h[1] > start:
                return True
        return False

    bad, covered, n = [], 0, 0
    for _ in range(600):
        kind = rng.random()
        if kind < 0.4:      # a real junction: tail, index, anchor, insert
            read = (rand(rng.randint(0, 150)) + cfg["p5_tail"] + rand(16)
                    + a5 + rand(rng.randint(20, 150)))
            start = read.rfind(a5)
        else:               # a look-alike, with edits, somewhere in insert
            ctx = rng.choice(contexts)
            pre = rand(rng.randint(0, 120))
            read = pre + mutate(ctx.replace("N", "A"), rng.randint(0, 8)) \
                + rand(rng.randint(0, 120))
            start = len(pre) + 24
        if kind > 0.9:      # an N near the hit: the unindexed path
            j = rng.randrange(max(0, start - 30), min(len(read), start + 40))
            read = read[:j] + "N" + read[j + 1:]
        end = min(len(read), start + 14)
        got, want = mask.covers(read, start, end), reference(read, start, end)
        n += 1
        covered += want
        if got != want:
            bad.append((read, start, end, got, want))
    check(f"indexed covers() matches every-context reference on {n} hits "
          f"({covered} covered)", not bad, str(bad[:2]))
    n_amb = sum(bool(set(c) - set("ACGT")) for c in mask.contexts)
    check(f"contexts with IUPAC codes ({n_amb}) are indexed, not tried at "
          f"every hit", n_amb > 0 and not mask._index.always,
          f"{len(mask._index.always)} always tried")


def main() -> int:
    test_earlier_worse_anchor_is_not_skipped()
    test_split_accounting()
    test_basecall_provenance()
    test_partial_second_unit()
    test_three_units()
    test_recovers_correct_samples()
    test_quantifies_what_it_fixes()
    test_clean_reads_untouched()
    test_sharded_matches_whole_run()
    test_cutting_preserves_indices()
    test_panel_lookalikes_are_not_junctions()
    test_seeded_scan_matches_reference()
    test_indexed_lookalike_mask()
    test_indexed_mask_matches_reference()
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
