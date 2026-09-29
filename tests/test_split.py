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


def main() -> int:
    test_basecall_provenance()
    test_partial_second_unit()
    test_three_units()
    test_recovers_correct_samples()
    test_quantifies_what_it_fixes()
    test_clean_reads_untouched()
    test_sharded_matches_whole_run()
    test_cutting_preserves_indices()
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
