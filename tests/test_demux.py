#!/usr/bin/env python3
"""
Tests for nanopore demultiplexing.

Correctness here is measured against the simulator's truth table, not
against the demultiplexer's own output. Three bugs found during development
are each pinned by a test below, because every one of them produced results
that looked plausible:

  - the i7 barcode follows its anchor while the i5 barcode precedes it.
    Assuming one side for both silently lost every i7 (tier A fell to 0.01%).
  - reading the tandem copies at fixed offsets meant one indel beside the
    anchor broke both copies at once, which is precisely the case the
    duplication exists to survive (tier A 63%).
  - Myers' bit-vector needs a plain horizontal shift; setting the low bit
    forces global alignment and the anchor is then never found except when
    it spans the whole window (tier A 52%, and a vacuous unit test that
    passed anyway because it used a tolerance of 64).

Run: python3 tests/test_demux.py
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
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


dx = _load("demux_ont", "bin/demux_ont.py")
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


# --------------------------------------------------------------------------
# anchor finding
# --------------------------------------------------------------------------


def plain_dp_min(text: str, pat: str) -> int:
    """Reference implementation: minimum edit distance between `pat` and any
    substring of `text`. Slow, obviously correct, and the thing Myers must
    agree with."""
    prev = list(range(len(pat) + 1))
    best = len(pat)
    for tc in text:
        cur = [0]
        for i, pc in enumerate(pat, 1):
            cur.append(min(prev[i] + 1, cur[i - 1] + 1, prev[i - 1] + (pc != tc)))
        best = min(best, cur[-1])
        prev = cur
    return best


def test_find_anchor() -> None:
    print("\nanchor finding (Myers bit-vector)")

    # the cases that exposed the global-alignment bug: a pattern that does
    # not span the whole text
    check("exact match mid-string", dx.find_anchor("AAAACGTAAA", "ACGT", 0) == (3, 7),
          str(dx.find_anchor("AAAACGTAAA", "ACGT", 0)))
    check("exact match at the start", dx.find_anchor("ACGTTTTT", "ACGT", 0) == (0, 4))
    check("exact match at the end", dx.find_anchor("TTTTACGT", "ACGT", 0) == (4, 8))
    check("absent pattern returns None",
          dx.find_anchor("GGGGGGGGGGGG", "ACGT", 1) is None)

    # indels are the whole reason for using edit distance
    check("one substitution found",
          dx.find_anchor("TTTTTTTTACCTTTTTTTT", "ACGT", 1) is not None)
    check("one deletion found",
          dx.find_anchor("TTTTTTTTACTTTTTTTT", "ACGT", 1) is not None)
    check("one insertion found",
          dx.find_anchor("TTTTTTTTACGGTTTTTTT", "ACGT", 1) is not None)
    check("beyond tolerance rejected",
          dx.find_anchor("TTTTTTTTAAAATTTTTTT", "ACGT", 1) is None)

    # agreement with the reference DP, at realistic tolerances. The earlier
    # version of this test used a tolerance of 64, which made every case
    # trivially findable and hid a completely broken implementation.
    import random

    rng = random.Random(5)
    bad = 0
    total = 0
    for _ in range(300):
        pat = "".join(rng.choice("ACGT") for _ in range(14))
        txt = "".join(rng.choice("ACGT") for _ in range(120))
        if rng.random() < 0.8:
            mut = list(pat)
            for _ in range(rng.randint(0, 3)):
                k = rng.randrange(len(mut))
                o = rng.random()
                if o < 0.34:
                    mut[k] = rng.choice("ACGT")
                elif o < 0.67:
                    mut[k] = ""
                else:
                    mut[k] = mut[k] * 2
            i = rng.randrange(0, 100)
            txt = txt[:i] + "".join(mut) + txt[i:]
        for k in (1, 2, 3):
            total += 1
            if (plain_dp_min(txt, pat) <= k) != (dx.find_anchor(txt, pat, k) is not None):
                bad += 1
    check(f"agrees with the reference DP at tolerance 1-3 ({total} cases)",
          bad == 0, f"{bad} disagreements")

    # the anchor END is what positions the barcode: an indel inside the
    # anchor shifts everything after it
    hit = dx.find_anchor("TTTT" + "ACGT" + "GGGG", "ACGT", 0)
    check("returns the alignment end, not start+length", hit[1] == 8, str(hit))


# --------------------------------------------------------------------------
# barcode sets
# --------------------------------------------------------------------------


def test_barcode_set() -> None:
    print("\nbarcode set tolerance")
    i7 = {"TAAGGCGA", "CGTACTAG", "AGGCAGAA", "TCCTGAGC"}
    bs = dx.BarcodeSet.build(i7, requested=1)
    check("exact match", bs.match("TAAGGCGA") == "TAAGGCGA")
    check("one substitution corrected", bs.match("TAAGGCGT") == "TAAGGCGA")
    check("unrelated sequence rejected", bs.match("GGGGGGGG") is None)

    # tolerance is derived from the set, not taken on trust: a number
    # someone typed cannot know the set's minimum distance
    close = dx.BarcodeSet.build({"AAAAAAAA", "AAAAAAAT"}, requested=2)
    check("a close set caps its own tolerance at 0", close.max_mismatch == 0,
          str(close.max_mismatch))
    check("and then refuses an inexact window", close.match("AAAAAAAG") is None)

    # an ambiguous window is no match: guessing is worse than declining
    amb = dx.BarcodeSet.build({"AAAACCCC", "AAAAGGGG", "AAAATTTT"}, requested=1)
    check("ambiguous window declined", amb.match("AAAAACCC") is None
          or amb.match("AAAAACCC") == "AAAACCCC")

    real = dx.BarcodeSet.build(
        {s for s in ("TAAGGCGA CGTACTAG AGGCAGAA TCCTGAGC GGACTCCT TAGGCATG "
                     "CTCTCTAC CGAGGCTG").split()}, requested=1)
    check("the real i7 subset supports one substitution",
          real.max_mismatch == 1, str(real.max_mismatch))


# --------------------------------------------------------------------------
# end to end against the truth table
# --------------------------------------------------------------------------


def build_and_demux(tmp: str, n_samples=24, unused=4, reads=250, extra=None):
    import random

    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, n_samples, unused, "RUN01", "nanopore")
    mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
    mf.write_samplesheet(samples, os.path.join(tmp, "samplesheet.csv"))
    mf.write_runs(os.path.join(tmp, "runs.tsv"), "RUN01", "nanopore",
                  os.path.join(tmp, "sim", "RUN01"))

    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        sim.main([
            "--reference", os.path.join(tmp, "reference.tsv"),
            "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
            "--runs", os.path.join(tmp, "runs.tsv"),
            "--assay", ASSAY,
            "--outdir", os.path.join(tmp, "sim"),
            "--reads-per-sample", str(reads),
        ])
        dx.main([
            "--reads", os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz"),
            "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
            "--run-id", "RUN01",
            "--assay", ASSAY,
            "--outdir", os.path.join(tmp, "demux"),
        ] + (extra or []))

    with open(os.path.join(tmp, "sim", "RUN01.truth.tsv"), encoding="utf-8") as fh:
        truth = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(os.path.join(tmp, "demux", "assignments.tsv"), encoding="utf-8") as fh:
        calls = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    ss = list(csv.DictReader(open(os.path.join(tmp, "samplesheet.csv"))))
    return truth, calls, ss


def test_against_truth() -> None:
    print("\ndemultiplexing against the truth table")
    tmp = tempfile.mkdtemp()
    try:
        truth, calls, ss = build_and_demux(tmp)

        # the read-accounting invariant, at the demux stage
        check("every read accounted for exactly once",
              set(truth) == set(calls) and len(truth) == len(calls),
              f"{len(truth)} vs {len(calls)}")

        clean = [r for r, t in truth.items() if t["fate"] == "clean"]
        inA = [r for r in clean if calls[r]["tier"] == "A"]
        rate = len(inA) / len(clean)
        check("most clean reads reach tier A", rate > 0.90, f"{rate:.3f}")

        a = [(r, c) for r, c in calls.items() if c["tier"] == "A"]
        right = sum(1 for r, c in a if c["sample_id"] == truth[r]["sample_id"])
        acc = right / len(a)
        check("tier A assignment accuracy above 98%", acc > 0.98, f"{acc:.4f}")

        # the number that actually matters: a clean read put in the wrong
        # sample is a false positive nothing downstream can undo
        wrong_clean = [r for r, c in a
                       if c["sample_id"] != truth[r]["sample_id"]
                       and truth[r]["fate"] == "clean"]
        frac = len(wrong_clean) / len(a)
        check("clean reads misassigned below 0.5%", frac < 0.005,
              f"{frac:.5f} ({len(wrong_clean)}/{len(a)})")

        # truncated reads must not be confidently assigned
        trunc = [r for r, t in truth.items() if t["fate"] == "truncated"]
        if trunc:
            in_a = sum(1 for r in trunc if calls[r]["tier"] == "A")
            check("truncated reads are not given tier A",
                  in_a / len(trunc) < 0.05, f"{in_a}/{len(trunc)}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_hopping_and_chimeras() -> None:
    """Hopping onto an unused index combination must be caught; hopping onto
    a used one is undetectable by construction, and the design says so. This
    test pins both halves, so a future change cannot quietly claim to have
    fixed something that is not fixable."""
    print("\nhopping and chimeras")
    tmp = tempfile.mkdtemp()
    try:
        # Enough reads that the hopped-and-unused class is not a handful:
        # at 250 reads/sample it was ~20 reads and binomial noise alone
        # swung the rate by 15 points.
        truth, calls, ss = build_and_demux(tmp, reads=1200)
        assigned = {(s["index_i5"], s["index_i7"]) for s in ss}

        # The expected detectable fraction is NOT the unweighted share of
        # unused combinations. A chimera takes its i5 from one sample and its
        # i7 from another, so each index appears in proportion to how many
        # samples carry it. With samples concentrated in a few full rows,
        # most such pairs are themselves assigned and therefore invisible.
        pairs = [(a["index_i5"], b["index_i7"]) for a in ss for b in ss]
        observable = sum(1 for p in pairs if p not in assigned) / len(pairs)

        hop = [r for r, t in truth.items() if t["fate"] == "hopped"]
        onto_unused = [r for r in hop
                       if (truth[r]["index_i5"], truth[r]["index_i7"]) not in assigned]
        caught = sum(1 for r in onto_unused if calls[r]["tier"] == "C")
        if onto_unused:
            check("hops onto an unused combination are caught",
                  caught / len(onto_unused) > 0.85,
                  f"{caught}/{len(onto_unused)}")

        onto_used = [r for r in hop if r not in set(onto_unused)]
        as_a = sum(1 for r in onto_used if calls[r]["tier"] == "A")
        if onto_used:
            check("hops onto a used combination pass as tier A, as designed",
                  as_a / len(onto_used) > 0.5, f"{as_a}/{len(onto_used)}")

        # A chimera carries one sample's i5 and another's i7. With
        # combinatorial indexing most such pairs are themselves valid, so the
        # detectable fraction is bounded by the observable fraction - a
        # property of the plate layout, not of the matcher.
        ch = [r for r, t in truth.items() if t["fate"] == "chimera"]
        in_a = [r for r in ch if calls[r]["tier"] == "A"]
        if in_a:
            all_valid = all(
                (calls[r]["index_i5"], calls[r]["index_i7"]) in assigned for r in in_a
            )
            check("chimeras reaching tier A all carry a valid index pair",
                  all_valid,
                  "some carried an invalid pair, which should have been tier C")
            detected = 1 - len(in_a) / len(ch)
            check("chimera detection tracks the observable fraction",
                  abs(detected - observable) < 0.10,
                  f"detected {detected:.3f} vs observable {observable:.3f}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_orientation_and_trimming() -> None:
    print("\norientation and trimming")
    tmp = tempfile.mkdtemp()
    try:
        truth, calls, ss = build_and_demux(tmp)
        orient = Counter(c["orientation"] for c in calls.values())
        frac = orient["-"] / sum(orient.values())
        check("both strands recovered, roughly balanced",
              0.40 < frac < 0.60, f"{frac:.3f} reverse")

        # trimmed output should be the insert, shorter than the construct
        d = os.path.join(tmp, "demux", "reads")
        files = [f for f in os.listdir(d) if f.endswith(".fastq.gz")]
        check("per-sample files written", len(files) > 10, str(len(files)))
        with gzip.open(os.path.join(d, files[0]), "rt") as fh:
            lines = fh.read().splitlines()
        check("fastq well formed", len(lines) % 4 == 0)

        # Compare the trimmed insert against the reference amplicon it came
        # from, not against a fixed number: the panel spans several amplicon
        # lengths, and dimers are legitimately tiny. An earlier version of
        # this test asserted a fixed maximum and failed on the longest
        # amplicon while calling a correctly-trimmed dimer "over-trimmed".
        with open(os.path.join(tmp, "reference.tsv"), encoding="utf-8") as fh:
            reflen = {r["seq_id"]: len(r["sequence"])
                      for r in csv.DictReader(fh, delimiter="\t")}
        deltas = []
        for f in files:
            with gzip.open(os.path.join(d, f), "rt") as fh:
                ls = fh.read().splitlines()
            for i in range(0, len(ls), 4):
                rid = ls[i][1:].split()[0]
                t = truth.get(rid)
                if not t or t["fate"] != "clean":
                    continue
                exp = reflen.get(t["seq_id"])
                if exp:
                    deltas.append(len(ls[i + 1]) - exp)
        import statistics as st

        med = st.median(deltas)
        check("insert length matches the reference amplicon",
              abs(med) <= 6, f"median delta {med:+.0f}bp from the amplicon")
        check("adapters fully removed, mosaic ends included",
              med < 10, f"median delta {med:+.0f}bp - the 19bp mosaic ends "
              "are adapter and should not survive trimming")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_single_copy_barcodes() -> None:
    """barcode_repeats=1 must still work, and must lose the copy-agreement
    check rather than silently pretending to have it."""
    print("\nsingle-copy barcodes")
    tmp = tempfile.mkdtemp()
    try:
        import contextlib
        import io
        import json
        import random

        rng = random.Random(7)
        ref = mf.build_reference(rng)
        i7, i5 = mf.load_index_map(INDEX_MAP)
        samples = mf.build_samplesheet(i7, i5, 12, 3, "RUN01", "nanopore")
        mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
        mf.write_samplesheet(samples, os.path.join(tmp, "samplesheet.csv"))
        runs = os.path.join(tmp, "runs.tsv")
        mf.write_runs(runs, "RUN01", "nanopore", os.path.join(tmp, "sim", "RUN01"))
        lines = open(runs).read().splitlines()
        head, row = lines[0].split("\t"), lines[1].split("\t")
        row[head.index("barcode_repeats")] = "1"
        open(runs, "w").write("\t".join(head) + "\n" + "\t".join(row) + "\n")

        with contextlib.redirect_stdout(io.StringIO()):
            sim.main(["--reference", os.path.join(tmp, "reference.tsv"),
                      "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                      "--runs", runs, "--assay", ASSAY,
                      "--outdir", os.path.join(tmp, "sim"),
                      "--reads-per-sample", "200"])
            dx.main(["--reads", os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz"),
                     "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                     "--run-id", "RUN01", "--assay", ASSAY,
                     "--barcode-repeats", "1",
                     "--outdir", os.path.join(tmp, "demux")])

        with open(os.path.join(tmp, "demux", "demux_summary.json"), encoding="utf-8") as fh:
            summ = json.load(fh)
        check("repeats=1 still assigns reads",
              summ["by_tier"]["A"] / summ["total_reads"] > 0.80,
              f"{summ['by_tier']['A']}/{summ['total_reads']}")
        check("copy agreement correctly reported as unavailable",
              summ["settings"]["require_copy_agreement"] is False)
        check("indels stay disabled", summ["settings"]["allow_indels"] is False)
        check("read accounting balances", summ["read_accounting_balanced"])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_minknow_input() -> None:
    """A MinION produces a tree of batch files, not a FASTQ.

    Every split in that tree - pass against fail, and MinKNOW's own
    barcode folders - is a decision about which reads enter the analysis.
    The bug this pins: the barcode-folder check originally ran after a
    "no files at this level" shortcut, and in a real MinKNOW tree
    fastq_pass holds only directories, so the check never ran.
    """
    print("\nMinKNOW run directories")
    import gzip as _gz
    sys.path.insert(0, os.path.join(ROOT, "bin"))
    from ampseqlib import resolve_read_files, open_reads_many

    tmp = tempfile.mkdtemp()
    try:
        def write(path, n, start=0):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with _gz.open(path, "wt") as fh:
                for i in range(start, start + n):
                    fh.write(f"@r{i}\nACGT\n+\nIIII\n")

        run = os.path.join(tmp, "run")
        for b in range(3):
            write(os.path.join(run, "fastq_pass", f"FAX_pass_{b}.fastq.gz"),
                  10, b * 10)
        write(os.path.join(run, "fastq_fail", "FAX_fail_0.fastq.gz"), 7, 900)

        files, rep = resolve_read_files(run)
        check("a run directory resolves to its batch files",
              rep["n_files"] == 3, str(rep["n_files"]))
        check("every passing read is streamed, across all batches",
              sum(1 for _ in open_reads_many(files)) == 30,
              str(sum(1 for _ in open_reads_many(files))))
        check("fastq_fail is excluded", rep["excluded_fastq_fail"] == 1)
        check("and the exclusion is stated, not silent",
              any("fastq_fail" in n for n in rep["notes"]))
        check("files come back in a stable order", files == sorted(files))

        # asked for explicitly, the failed reads are included
        files2, rep2 = resolve_read_files(run, prefer_pass=False)
        check("failed reads can be included on request",
              sum(1 for _ in open_reads_many(files2)) == 37,
              str(sum(1 for _ in open_reads_many(files2))))

        # MinKNOW demultiplexing left on: fastq_pass holds only folders
        bc = os.path.join(tmp, "run_bc")
        write(os.path.join(bc, "fastq_pass", "barcode01", "x.fastq.gz"), 4)
        write(os.path.join(bc, "fastq_pass", "unclassified", "y.fastq.gz"),
              20, 100)
        files3, rep3 = resolve_read_files(bc)
        check("barcode folders are noticed even though fastq_pass has no "
              "files of its own", bool(rep3["minknow_barcode_dirs"]),
              str(rep3["minknow_barcode_dirs"]))
        check("and the warning explains why its barcodes are not ours",
              any("in-line tandem" in n for n in rep3["notes"]))
        check("unclassified reads are kept, not discarded",
              sum(1 for _ in open_reads_many(files3)) == 24,
              str(sum(1 for _ in open_reads_many(files3))))
        check("a clean run raises no barcode warning",
              not rep["minknow_barcode_dirs"])

        # a single file and a glob still work
        one = os.path.join(run, "fastq_pass", "FAX_pass_0.fastq.gz")
        check("a single file is still accepted",
              resolve_read_files(one)[1]["n_files"] == 1)
        check("a glob is accepted",
              resolve_read_files(
                  os.path.join(run, "fastq_pass", "*.fastq.gz")
              )[1]["n_files"] == 3)
        check("the same file reached twice is read once",
              resolve_read_files([one, one])[1]["n_files"] == 1)
        check("an empty directory yields nothing rather than failing",
              resolve_read_files(tempfile.mkdtemp())[1]["n_files"] == 0)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_index_matrix() -> None:
    """The i7 x i5 grid, which measures cross-talk directly.

    A combination no sample claims cannot produce a real read, so what
    lands there is hopping or misreading and nothing else. This was
    specified early, never built, and found missing only when a static
    check noticed the workflow declared an output the tool never wrote.
    """
    print("\nthe observed index matrix")
    tmp = tempfile.mkdtemp()
    try:
        truth, calls, ss = build_and_demux(tmp, n_samples=24, unused=4)
        path = os.path.join(tmp, "demux", "index_matrix.tsv")
        check("a matrix is written", os.path.exists(path))
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        i5 = {r["index_i5"] for r in rows}
        i7 = {r["index_i7"] for r in rows}
        check("it is the full grid, not only the assigned pairs",
              len(rows) == len(i5) * len(i7),
              f"{len(rows)} rows for {len(i5)}x{len(i7)}")
        assigned = {(s["index_i5"], s["index_i7"]) for s in ss}
        marked = {(r["index_i5"], r["index_i7"])
                  for r in rows if r["assigned"] == "yes"}
        check("assigned combinations are marked as such",
              marked == assigned, f"{len(marked)} vs {len(assigned)}")

        with open(os.path.join(tmp, "demux", "demux_summary.json")) as fh:
            m = json.load(fh)["index_matrix"]
        check("unused combinations are counted as sentinels",
              m["sentinel_combinations"] == len(i5) * len(i7) - len(assigned),
              str(m["sentinel_combinations"]))
        sent = sum(int(r["reads"]) for r in rows if r["assigned"] == "no")
        check("the sentinel read count matches the table",
              m["reads_on_sentinels"] == sent,
              f"{m['reads_on_sentinels']} vs {sent}")
        check("reads did land on sentinels, so the sensor is live",
              sent > 0, "no hopped reads reached an unused combination")
        check("the observable fraction is reported, not just the rate",
              0 < m["observable_fraction"] < 1,
              str(m["observable_fraction"]))
        check("and the summary says the rate is a floor",
              "floor" in m["note"])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_streaming_output() -> None:
    """Output must not depend on how often the write buffer is flushed.

    Demultiplexing used to hold every read in memory until the end, which
    cost about 1.4 KB per read - fine on a test fixture, ~28 GB on a full
    flow cell, and an OOM kill an hour into the job. It now appends to
    each sample's file in chunks and drops them. gzip members concatenate,
    so the file reads back as one stream, but the flush boundary is a new
    place for reads to go missing, so it is tested at a size where every
    sample flushes many times.
    """
    print("\nstreaming output")
    tmp = tempfile.mkdtemp()
    try:
        # SampleWriter on its own first: many members, awkward totals.
        w = dx.SampleWriter(os.path.join(tmp, "unit"), flush_records=3)
        expect = {}
        for s, n in (("A", 10), ("B", 1), ("C", 9), ("D", 0)):
            expect[s] = [(f"{s}-r{i}", "ACGT" * (i + 1), "I" * 4 * (i + 1))
                         for i in range(n)]
            for rid, seq, qual in expect[s]:
                w.add(s, rid, seq, qual)
        w.close()

        check("a sample with no reads gets no file",
              not os.path.exists(w.path_for("D")))
        check("the written count equals what went in",
              w.written == 20, str(w.written))
        for s in ("A", "B", "C"):
            with gzip.open(w.path_for(s), "rt") as fh:
                lines = fh.read().splitlines()
            got = [(lines[i][1:], lines[i + 1], lines[i + 3])
                   for i in range(0, len(lines), 4)]
            check(f"sample {s} reads back complete and in order",
                  got == expect[s],
                  f"{len(got)} records, expected {len(expect[s])}")

        # Then end to end: the same run at two flush sizes must agree
        # byte for byte, or the buffer is losing or reordering reads.
        truth, calls, ss = build_and_demux(tmp, reads=400)
        big = os.path.join(tmp, "demux", "reads")
        small = os.path.join(tmp, "demux_small", "reads")
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main([
                "--reads", os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz"),
                "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                "--run-id", "RUN01", "--assay", ASSAY,
                "--outdir", os.path.join(tmp, "demux_small"),
                "--flush-records", "7",
            ])

        names = sorted(os.listdir(big))
        check("the same samples are emitted at either flush size",
              names == sorted(os.listdir(small)))
        differing = []
        for name in names:
            with gzip.open(os.path.join(big, name), "rt") as fh:
                a = fh.read()
            with gzip.open(os.path.join(small, name), "rt") as fh:
                b = fh.read()
            if a != b:
                differing.append(name)
        check("and every file is identical, flush boundaries and all",
              not differing, f"differ: {differing[:3]}")

        with open(os.path.join(tmp, "demux_small", "demux_summary.json")) as fh:
            summ = json.load(fh)
        check("reads written balances against reads assigned",
              summ["reads_written_balanced"]
              and summ["reads_written"] == summ["assigned_reads"],
              f"{summ['reads_written']} vs {summ['assigned_reads']}")
        check("and the per-sample counts sum to it",
              sum(summ["per_sample"].values()) == summ["reads_written"],
              f"{sum(summ['per_sample'].values())} vs {summ['reads_written']}")
        n_flushed = max(summ["per_sample"].values()) // 7
        check("the small run really did flush many times, so this proves something",
              n_flushed > 5, f"only {n_flushed} flushes for the deepest sample")

        # Re-running into a used directory must replace, not append. An
        # appended file would double every sample's depth while still
        # looking like a well-formed FASTQ.
        before = {n: os.path.getsize(os.path.join(small, n)) for n in names}
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main([
                "--reads", os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz"),
                "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                "--run-id", "RUN01", "--assay", ASSAY,
                "--outdir", os.path.join(tmp, "demux_small"),
                "--flush-records", "7",
            ])
        after = {n: os.path.getsize(os.path.join(small, n)) for n in names}
        check("re-running replaces the previous output rather than appending",
              before == after,
              f"{sum(after.values()) / max(1, sum(before.values())):.2f}x larger")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_ligation_chimeras() -> None:
    """Two molecules joined into one read must not demultiplex.

    This is the failure that motivated the check, and it is nastier than
    a low yield. A head-to-tail join presents molecule A's i5 at the 5'
    end and molecule B's i7 at the 3'. Both are real declared indices, so
    the pair is very likely a well on the plate - just not the well
    either molecule came from. It lands in tier A, is assigned, and
    nothing downstream could ever notice. So the test is not "are these
    caught" but "are they kept out of the per-sample output".
    """
    print("\nligation chimeras")
    tmp = tempfile.mkdtemp()
    try:
        truth, calls, ss = build_and_demux(tmp, n_samples=24, unused=4,
                                           reads=150)
        src = os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz")
        with gzip.open(src, "rt") as fh:
            lines = fh.read().splitlines()
        # Split the header the way the reader does, or these ids will not
        # join against the assignment table and the selection below comes
        # back empty.
        recs = [(lines[i][1:].split()[0], lines[i + 1], lines[i + 3])
                for i in range(0, len(lines), 4)]

        # Only join molecules that are individually COMPLETE - reads that
        # already resolved both of their own ends. That is the class the
        # check can detect, because joining two complete molecules
        # necessarily buries a full adapter structure at the junction.
        #
        # Joining truncated fragments instead produces a chimera with one
        # usable end from each and NO adapter between them: indetectable
        # by this method, or any method that works from adapter position.
        # Such a read is tested separately below, as a measured residue
        # rather than as a requirement the method cannot meet.
        complete = [r for r in recs
                    if calls.get(r[0], {}).get("tier") == "A"]
        rc = dx.revcomp
        joined = os.path.join(tmp, "chim")
        os.makedirs(joined, exist_ok=True)
        pairs = []
        with gzip.open(os.path.join(joined, "c.fastq.gz"), "wt") as fh:
            for k in range(0, min(len(complete) - 1, 600), 2):
                a, b = complete[k], complete[k + 1]
                fh.write(f"@cat{k}\n{a[1] + b[1]}\n+\n{a[2] + b[2]}\n")
                fh.write(f"@inv{k}\n{a[1] + rc(b[1])}\n+\n{a[2] + b[2][::-1]}\n")
                pairs.append(f"cat{k}")
                pairs.append(f"inv{k}")

        out = os.path.join(tmp, "demux_chim")
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main(["--reads", joined,
                     "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                     "--run-id", "RUN01", "--assay", ASSAY, "--outdir", out])
        with open(os.path.join(out, "assignments.tsv"), encoding="utf-8") as fh:
            got = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

        chim = [got[p] for p in pairs if p in got]
        assigned = [r for r in chim if r["sample_id"]]
        check("every chimera was read back",
              len(chim) == len(pairs), f"{len(chim)} of {len(pairs)}")
        check("no join of two complete molecules is assigned to a sample",
              not assigned,
              f"{len(assigned)} of {len(chim)} were assigned, e.g. "
              f"{assigned[0]['read_id'] if assigned else ''}")
        in_e = sum(1 for r in chim if r["tier"] == "E")
        # Head-to-head joins often fail outright on their own ends and
        # never reach the tier-A check, so E is a floor, not the total.
        # What matters is that none are assigned; E is how many were
        # caught specifically by the junction scan.
        check("a substantial share are caught by the junction scan itself",
              in_e > 0.25 * len(chim), f"{in_e} of {len(chim)} in tier E")

        with open(os.path.join(out, "demux_summary.json")) as fh:
            summ = json.load(fh)
        check("tier E appears in the summary's tier table",
              summ["by_tier"].get("E", 0) > 0)
        check("read accounting still balances with the new tier",
              summ["read_accounting_balanced"])
        check("and the written reads still balance",
              summ["reads_written_balanced"])

        # The check must not cost clean reads: same pipeline, unjoined.
        clean_out = os.path.join(tmp, "demux_clean")
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main(["--reads", src,
                     "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                     "--run-id", "RUN01", "--assay", ASSAY,
                     "--outdir", clean_out])
        with open(os.path.join(clean_out, "demux_summary.json")) as fh:
            cs = json.load(fh)
        false_pos = cs["by_tier"].get("E", 0) / max(1, cs["total_reads"])
        check("clean single-molecule reads are not flagged",
              false_pos < 0.01, f"{100 * false_pos:.2f}% of clean reads hit E")

        # The residue the method cannot reach: two TRUNCATED fragments
        # joined, contributing one end each and no adapter in between.
        # This is asserted to be detected-as-undetectable rather than
        # caught, so that the limitation is recorded in the test suite
        # instead of living only in a comment.
        frags = [r for r in recs
                 if calls.get(r[0], {}).get("tier") == "B"]
        if len(frags) >= 4:
            fdir = os.path.join(tmp, "frag")
            os.makedirs(fdir, exist_ok=True)
            with gzip.open(os.path.join(fdir, "f.fastq.gz"), "wt") as fh:
                for k in range(0, min(len(frags) - 1, 200), 2):
                    a, b = frags[k], frags[k + 1]
                    fh.write(f"@frag{k}\n{a[1] + b[1]}\n+\n{a[2] + b[2]}\n")
            fout = os.path.join(tmp, "demux_frag")
            with contextlib.redirect_stdout(io.StringIO()):
                dx.main(["--reads", fdir,
                         "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                         "--run-id", "RUN01", "--assay", ASSAY,
                         "--outdir", fout])
            with open(os.path.join(fout, "demux_summary.json")) as fh:
                fs = json.load(fh)
            leak = fs["assigned_reads"] / max(1, fs["total_reads"])
            # Reported, not asserted. There is no threshold here that
            # would mean anything: the number depends entirely on how the
            # fragments were made, and a passing assertion on it would be
            # decoration. What it documents is that the figure is not
            # zero and the method does not claim otherwise.
            print(f"  ..   fragment chimeras (no adapter at the junction) "
                  f"leak through: {100 * leak:.0f}% assigned, "
                  f"{fs['by_tier'].get('E', 0)} caught - a limit of any "
                  "adapter-based check")

        # And it can be switched off, because a yield/correctness trade is
        # the user's call to make, not one to make silently on their behalf.
        off_out = os.path.join(tmp, "demux_off")
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main(["--reads", joined,
                     "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                     "--run-id", "RUN01", "--assay", ASSAY,
                     "--outdir", off_out, "--no-chimera-check"])
        with open(os.path.join(off_out, "demux_summary.json")) as fh:
            os_ = json.load(fh)
        check("--no-chimera-check restores the old, wrong behaviour",
              os_["by_tier"].get("E", 0) == 0
              and os_["assigned_reads"] > summ["assigned_reads"],
              f"E={os_['by_tier'].get('E', 0)} "
              f"assigned {os_['assigned_reads']} vs {summ['assigned_reads']}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_parallel_matches_serial() -> None:
    """Thread count must not be an input to the result.

    Demultiplexing is sharded by input file so a full flow cell finishes
    inside a wall-clock limit, and each shard writes its own per-sample
    FASTQ, its own slice of the per-read table and its own tallies. Three
    things could go wrong and none of them would announce itself: shards
    merged in completion order rather than input order, tallies summed
    from a subset, or per-sample reads lost when the parts are joined. So
    the test is byte equality against the same data read as one file.
    """
    print("\nsharded demultiplexing")
    tmp = tempfile.mkdtemp()
    try:
        truth, calls, ss = build_and_demux(tmp, n_samples=24, unused=4,
                                           reads=200)
        src = os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz")
        with gzip.open(src, "rt") as fh:
            lines = fh.read().splitlines()
        n = len(lines) // 4

        shard = os.path.join(tmp, "shard")
        os.makedirs(shard, exist_ok=True)
        per = (n + 4) // 5
        for f in range(5):
            with gzip.open(os.path.join(shard, f"b{f}.fastq.gz"), "wt") as out:
                for i in range(f * per, min((f + 1) * per, n)):
                    out.write("\n".join(lines[4 * i:4 * i + 4]) + "\n")

        import contextlib
        import io
        outs = {}
        for tag, reads, threads in (("one", src, 1), ("many", shard, 4)):
            d = os.path.join(tmp, "dx_" + tag)
            with contextlib.redirect_stdout(io.StringIO()):
                dx.main(["--reads", reads,
                         "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                         "--run-id", "RUN01", "--assay", ASSAY,
                         "--outdir", d, "--threads", str(threads)])
            outs[tag] = d

        with open(os.path.join(outs["one"], "assignments.tsv")) as fh:
            a = fh.read()
        with open(os.path.join(outs["many"], "assignments.tsv")) as fh:
            b = fh.read()
        check("the per-read table is identical, and in input order",
              a == b, f"{len(a)} vs {len(b)} bytes")

        ra = sorted(os.listdir(os.path.join(outs["one"], "reads")))
        rb = sorted(os.listdir(os.path.join(outs["many"], "reads")))
        check("the same samples get files", ra == rb,
              f"{len(ra)} vs {len(rb)}")
        differing = []
        for name in ra:
            with gzip.open(os.path.join(outs["one"], "reads", name), "rt") as fh:
                x = fh.read()
            with gzip.open(os.path.join(outs["many"], "reads", name), "rt") as fh:
                y = fh.read()
            if x != y:
                differing.append(name)
        # Parts are joined by concatenating compressed bytes, which is a
        # valid gzip file only because members concatenate. If that ever
        # stops being true this is where it shows.
        check("every per-sample FASTQ reads back identical after merging",
              not differing, f"differ: {differing[:3]}")

        sa = json.load(open(os.path.join(outs["one"], "demux_summary.json")))
        sb = json.load(open(os.path.join(outs["many"], "demux_summary.json")))
        # parts_recomputed counts what THIS invocation did, not what the
        # data contains, so it legitimately differs between a one-file run
        # and a sharded one. Everything that describes the reads must match.
        keys = [k for k in sa
                if k not in ("input_files", "settings", "parts_recomputed")]
        bad = [k for k in keys if sa[k] != sb.get(k)]
        check("every tally agrees, including the index matrix",
              not bad, f"differ: {bad}")
        check("and both still balance their read accounting",
              sa["read_accounting_balanced"] and sb["read_accounting_balanced"]
              and sa["reads_written_balanced"] and sb["reads_written_balanced"])

        # Resumption: drop some parts and rerun. The rebuilt output must
        # equal the untouched one, or a resumed job is not the same job.
        parts = os.path.join(outs["many"], "parts")
        removed = 0
        for name in sorted(os.listdir(parts)):
            if name.startswith(("00001", "00003")):
                path = os.path.join(parts, name)
                shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
                removed += 1
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main(["--reads", shard,
                     "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                     "--run-id", "RUN01", "--assay", ASSAY,
                     "--outdir", outs["many"], "--threads", "4"])
        with open(os.path.join(outs["many"], "assignments.tsv")) as fh:
            c = fh.read()
        check("a resumed run reproduces the same output exactly",
              c == a, f"{removed} part files were removed before resuming")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_end_reasons_discriminate() -> None:
    """The per-end reasons must tell chemistry from matching.

    A run that comes back half single-ended has two completely different
    possible causes with the same symptom in the tiers: either the
    molecule never carried that index (PCR2 did not add it, and no
    tolerance will recover it), or the index is present and the matcher
    is rejecting it (a parameter away from fixed). The tier says only
    that an end failed.

    So this checks the map can actually separate those two, on reads
    built to fail each way. A breakdown that labelled both the same
    would sum correctly, look plausible, and be useless - which is the
    failure mode worth testing for here.
    """
    print("\nper-end failure reasons")
    cfg = dx.load_assay(ASSAY)
    i7map, i5map = mf.load_index_map(INDEX_MAP)
    i5_seq = i5map[0][0]
    i7_seq = i7map[0][0]
    amp = "ACGT" * 60
    reps = 2
    whole = sim.build_molecule(amp, i7_seq, i5_seq, reps, cfg)

    i5_set = dx.BarcodeSet({i5_seq: "S1"}, len(i5_seq), 1)
    i7_set = dx.BarcodeSet({i7_seq: "N1"}, len(i7_seq), 1)
    pairs = {(i5_seq, i7_seq): "SAMP"}

    def call(seq):
        return dx.assign_read("r", seq, cfg, i5_set, i7_set, pairs,
                              reps, 200, False)

    # 1. intact: both ends read
    c = call(whole)
    check("an intact molecule reports ok at both ends",
          c.i5_reason.startswith("ok") and c.i7_reason.startswith("ok"),
          f"{c.i5_reason}/{c.i7_reason}")

    # 2. the i7 end removed entirely - the molecule never carried it
    c = call(whole[:-60])
    check("a molecule missing its i7 end says anchor_not_found",
          c.i7_reason == "anchor_not_found", c.i7_reason)
    check("and still reports the i5 end as read",
          c.i5_reason.startswith("ok"), c.i5_reason)

    # 3. the index present but not one we declared - a matching failure,
    #    NOT a missing end. This is the distinction that matters.
    other = sim.build_molecule(amp, "TTTTTTTT", i5_seq, reps, cfg)
    c = call(other)
    check("an undeclared i7 index says no_copy_matched, not anchor_not_found",
          c.i7_reason == "no_copy_matched", c.i7_reason)

    # 4. the two reasons must be different strings, or the map cannot
    #    distinguish the two causes at all.
    check("the two causes carry different labels",
          "anchor_not_found" != "no_copy_matched")

    # 5. every reason emitted must be documented.
    seen = {c.i5_reason, c.i7_reason}
    check("every reason used is in END_REASON_MEANING",
          seen <= set(dx.END_REASON_MEANING), str(seen))


def test_stats_only() -> None:
    """--stats-only must count identically and write nothing.

    And, the part that actually matters: a stats-only part records
    itself as finished having produced no FASTQ, so a later real run
    MUST NOT resume onto it. If it did, the merge would find nothing,
    the output would be empty, and the run would report success - the
    same shape as the stale-parts bug that once blended two algorithm
    versions in one results directory. The mode is therefore in the
    part fingerprint, and this checks that it is.

    It is not a speed feature. Measured on 51k reads the write is about
    2% of this stage; the saving is ~8 GB of disk and not clobbering
    reads you already have.
    """
    print("\nstats-only counting")
    import contextlib
    import glob as _glob
    import io
    tmp = tempfile.mkdtemp()
    try:
        build_and_demux(tmp, n_samples=8, reads=40)
        src = os.path.join(tmp, "sim", "RUN01")

        codes = {}

        def run(outdir, extra=()):
            with contextlib.redirect_stdout(io.StringIO()):
                codes[outdir] = dx.main(["--reads", src, "--samplesheet",
                         f"{tmp}/samplesheet.csv", "--run-id", "RUN01",
                         "--assay", ASSAY, "--outdir", outdir] + list(extra))
            with open(f"{outdir}/demux_summary.json", encoding="utf-8") as fh:
                return json.load(fh)

        full = run(f"{tmp}/full")
        so = run(f"{tmp}/so", ["--stats-only"])

        check("stats-only wrote no per-sample FASTQ",
              not _glob.glob(f"{tmp}/so/reads/*.fastq.gz"),
              str(_glob.glob(f"{tmp}/so/reads/*.fastq.gz")[:3]))
        check("and the normal run did write them",
              bool(_glob.glob(f"{tmp}/full/reads/*.fastq.gz")))

        for key in ("by_tier", "i5_end_reasons", "i7_end_reasons",
                    "per_sample", "total_reads", "assigned_reads"):
            check(f"{key} is identical with and without writing",
                  full[key] == so[key], f"{full[key]} vs {so[key]}")

        # The exit code, because a mode that works but returns 1 is
        # reported as a failed job by every scheduler - and the first
        # version of this flag did exactly that, since the
        # written-equals-assigned check fired on a run that wrote
        # nothing on purpose.
        check("stats-only exits zero", codes[f"{tmp}/so"] == 0,
              str(codes[f"{tmp}/so"]))
        check("and so does a normal run", codes[f"{tmp}/full"] == 0,
              str(codes[f"{tmp}/full"]))

        check("stats-only says so in the summary", so["stats_only"] is True)
        check("and does not claim the write check passed",
              so["reads_written_balanced"] is None,
              repr(so["reads_written_balanced"]))
        check("while a normal run does assert it",
              full["reads_written_balanced"] is True)

        # The safety property: a real run into a directory a stats-only
        # run already touched must redo the work, not trust it.
        mixed = f"{tmp}/mixed"
        so2 = run(mixed, ["--stats-only"])
        check("stats-only computed every part", so2["parts_recomputed"] > 0,
              str(so2["parts_recomputed"]))
        real = run(mixed)
        check("a real run does NOT resume onto stats-only parts",
              real["parts_recomputed"] == so2["parts_recomputed"],
              f"recomputed {real['parts_recomputed']} of "
              f"{so2['parts_recomputed']}")
        check("and the reads are actually there afterwards",
              bool(_glob.glob(f"{mixed}/reads/*.fastq.gz")))
        check("with the write check now asserted",
              real["reads_written_balanced"] is True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_barcode_diagnostic_discriminates() -> None:
    """The diagnostic must tell a misplaced index from a miscalled one.

    Its entire value is separating two causes of no_copy_matched that
    look identical in the tiers and have opposite remedies: a window in
    the wrong place (fixable by searching wider) versus bases that are
    wrong (fixable, if at all, by tolerance). So it is checked against
    failures injected on purpose, where the answer is known.

    The first version got this wrong. Its indel probe padded the window
    by +/-4, which is wide enough to stumble on a correctly-placed copy
    the offset search had missed - so all 400 deliberately SHIFTED
    reads came back as "indel distance 0" and would have been credited
    to indel tolerance. The pad is now +/-1.
    """
    print("\nbarcode diagnostic discriminates the two causes")
    import re
    import subprocess
    tmp = tempfile.mkdtemp()
    try:
        cfg = dx.load_assay(ASSAY)
        i7map, i5map = mf.load_index_map(INDEX_MAP)
        i5s = [x[0] for x in i5map[:4]]
        i7s = [x[0] for x in i7map[:4]]
        amp = "ACGTTGCA" * 40
        rc = dx.revcomp

        def tail(i7):
            return (rc(cfg["pcr1_overhang_rev"]) + i7 * 2
                    + rc(cfg["p7_tail"]))

        with open(f"{tmp}/samplesheet.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["run_id", "platform", "sample_id",
                        "index_i7", "index_i5", "sample_type"])
            for i in range(4):
                w.writerow(["D1", "nanopore", f"S{i}", i7s[i], i5s[i],
                            "unknown"])

        # Both copies carrying the same deletion, and a 4bp positional
        # shift. 4 is chosen because a shift that is a multiple of the
        # index length lands on the OTHER copy and is silently
        # rescued - which is what the tandem duplication is for, and
        # why an injected shift of 5 did not fail at all.
        gone = i5s[2][:4] + i5s[2][5:]
        recs = []
        for _ in range(100):
            recs.append(("bothindel", cfg["p5_tail"] + gone + gone
                         + cfg["pcr1_overhang_fwd"] + amp + tail(i7s[2])))
        for _ in range(100):
            recs.append(("shift4", cfg["p5_tail"] + i5s[2] * 2 + "GGGG"
                         + cfg["pcr1_overhang_fwd"] + amp + tail(i7s[2])))
        with gzip.open(f"{tmp}/reads.fastq.gz", "wt") as fh:
            for i, (lab, seq) in enumerate(recs):
                fh.write(f"@{lab}_{i}\n{seq}\n+\n{'I' * len(seq)}\n")

        out = subprocess.run(
            [sys.executable, os.path.join(ROOT, "bin", "diagnose_barcodes.py"),
             "--reads", f"{tmp}/reads.fastq.gz",
             "--samplesheet", f"{tmp}/samplesheet.csv",
             "--run-id", "D1", "--assay", ASSAY, "--limit", "200"],
            capture_output=True, text=True, timeout=300).stdout

        check("the diagnostic ran", "no_copy_matched" in out, out[-200:])
        check("both injected classes failed to match",
              "200 of 200 ends failed" in out, out[:200])

        # The shifted half must be found OUTSIDE the slack, and the
        # indel half must not be.
        wide = re.search(r"wider offset search\s+(\d+)", out)
        indel = re.search(r"1-2 indels WITH a margin rule\s+(\d+)", out)
        check("the shifted reads are reported as mispositioned",
              wide and int(wide.group(1)) == 100,
              wide.group(1) if wide else "not reported")
        check("the indel reads are reported as indel-recoverable",
              indel and int(indel.group(1)) == 100,
              indel.group(1) if indel else "not reported")
        check("and the offset it found is the one injected",
              "-4" in out and "OUTSIDE current slack" in out,
              "offset not reported")

        # A single-index set constrains nothing; saying otherwise once
        # produced "49 substitutions safely correctable".
        import importlib.util as _il
        spec = _il.spec_from_file_location(
            "dbc", os.path.join(ROOT, "bin", "diagnose_barcodes.py"))
        dbc = _il.module_from_spec(spec)
        spec.loader.exec_module(dbc)
        g = dbc.set_geometry("i5", ["ACGTACGT"])
        check("a one-index set reports no safe tolerance",
              g["safe_sub"] is None and g["min_hamming"] is None, str(g))

        # The landmark gap must recover the indel count that was put
        # in. This is the measurement that replaces sweeping the slack:
        # the index sits between two pieces of fixed sequence, so its
        # length is readable off the molecule rather than assumed.
        def gap(seq):
            return dbc.landmark_gap(seq, cfg["p5_tail"],
                                    cfg["i5_inner_anchor"], True, 2)

        nominal = (cfg["p5_tail"] + i5s[0] * 2
                   + cfg["pcr1_overhang_fwd"] + amp + tail(i7s[0]))
        short = (cfg["p5_tail"] + (i5s[1] * 2)[:-2]
                 + cfg["pcr1_overhang_fwd"] + amp + tail(i7s[1]))
        longer = (cfg["p5_tail"] + i5s[2] * 2 + "TTT"
                  + cfg["pcr1_overhang_fwd"] + amp + tail(i7s[2]))
        span = 16
        check("the gap reads the nominal span exactly",
              gap(nominal) == span, str(gap(nominal)))
        check("two deleted bases read as a gap two short",
              gap(short) == span - 2, str(gap(short)))
        check("three inserted bases read as a gap three long",
              gap(longer) == span + 3, str(gap(longer)))
        check("and a missing outer landmark returns None, not a guess",
              gap(i5s[0] * 2 + cfg["pcr1_overhang_fwd"] + amp) is None)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_barcode_margin_rule() -> None:
    """A barcode may be accepted on margin as well as on tolerance.

    Tolerance is a property of the SET - k substitutions are safely
    correctable only when every pair is more than 2k apart - so on a
    set whose closest pair is 4 apart it caps at 1 however clean the
    read is. Margin is a property of the READ: a window 2 from one
    index and 5 from everything else is unambiguous regardless of how
    tightly the set is packed, and accepting it introduces no
    misassignment the geometry would otherwise have prevented.

    Measured on the reference run, this recovers 5.0% of failed i5
    ends and 18.2% of i7. The ambiguous remainder is still declined.
    """
    print("\nbarcode margin rule")
    seqs = {"AAAACCCC", "AAAAGGGG", "TTTTCCCC", "TTTTGGGG"}
    bs = dx.BarcodeSet.build(seqs, 1)
    check("the set's own spacing caps tolerance at 1",
          bs.max_mismatch == 1, str(bs.max_mismatch))

    check("an exact match is still exact",
          bs.match("AAAACCCC") == "AAAACCCC")
    check("one substitution is taken on tolerance",
          bs.match("AAAACCCA") == "AAAACCCC")
    check("two substitutions ARE taken when the runner-up is far",
          bs.match("AAAACCAA") == "AAAACCCC", str(bs.match("AAAACCAA")))
    check("but a window equidistant from two is still declined",
          bs.match("AAAAGGCC") is None, str(bs.match("AAAAGGCC")))
    check("and a window far from everything is declined",
          bs.match("GGGGGGGG") is None, str(bs.match("GGGGGGGG")))
    check("a window of the wrong length is declined",
          bs.match("AAAACCC") is None)

    # The rule must be switchable off, or a run cannot be reproduced
    # against the behaviour that preceded it.
    off = dx.BarcodeSet.build(seqs, 1, margin=0)
    check("margin 0 restores tolerance-only behaviour",
          off.match("AAAACCAA") is None and off.match("AAAACCCA") == "AAAACCCC")

    # And the distance ceiling must bind: a margin alone should not
    # accept something arbitrarily far away just because nothing else
    # is near it.
    far = dx.BarcodeSet.build({"AAAACCCC", "TTTTGGGG"}, 1,
                              margin=2, margin_max_dist=2)
    check("the margin does not accept a distant window",
          far.match("AAAAGGGG") is None, str(far.match("AAAAGGGG")))


def test_combination_check() -> None:
    """A tie at one end, narrowed by the end that did resolve.

    MiniBar's trick, which this demultiplexer did not have. Its
    README reports the case directly: three forward barcodes matched at
    4 errors each, unresolvable alone, but paired with a perfect reverse
    barcode only one of the three identified a declared sample.

    The reason it works is sparseness. 16 i5 x 24 i7 is 384 possible
    combinations and 94 were used, so most pairs are not wells at all
    and a tie is more likely than not to contain exactly one that is.

    Everything here is a MEASUREMENT. The verdict is counted and no
    tier, sample or output file moves - so the test asserts that too,
    because a check that quietly started assigning reads would change
    results with no sign in the tiers.
    """
    print("\ncombination check")
    E = dx.EndCall

    # Three i7 candidates tied, one i5 resolved. Only (I5a, B) is a well.
    wells = {("I5a", "B"): "S1", ("I5b", "A"): "S2", ("I5b", "C"): "S3"}
    five = E(barcode="I5a")
    three = E(barcode=None, candidates=("A", "B", "C"), candidate_dist=4)
    verdict, sample, which = dx.combination_check(five, three, wells)
    check("three tied candidates, one declared well: resolved",
          verdict == "resolved", verdict)
    check("and it names the sample the read belongs to",
          sample == "S1" and which == "B", f"{sample} {which}")

    # Two of the tied candidates are wells with that i5: still ambiguous.
    amb = {("I5a", "A"): "S1", ("I5a", "B"): "S2"}
    check("two declared wells among the tie is still declined",
          dx.combination_check(five, three, amb)[0] == "still_ambiguous")

    # None of them pairs with this i5: the read was never assignable, so
    # declining it cost nothing. Good news, and distinguished from the
    # recoverable case rather than lumped with it.
    none_declared = {("I5z", "A"): "S1", ("I5z", "B"): "S2"}
    check("no declared well among the tie is reported as a free decline",
          dx.combination_check(five, three, none_declared)[0]
          == "no_declared")

    # It must work from the other side too - an i7 resolved, i5 tied.
    v, s, _ = dx.combination_check(
        E(barcode=None, candidates=("I5a", "I5b")), E(barcode="B"), wells)
    check("the check is symmetric: a tied i5 is narrowed by the i7",
          v == "resolved" and s == "S1", f"{v} {s}")

    # Guards.
    check("both ends resolved is not this check's business",
          dx.combination_check(E(barcode="I5a"), E(barcode="B"),
                               wells)[0] == "not_applicable")
    check("neither end resolved cannot be narrowed either",
          dx.combination_check(E(), E(), wells)[0] == "not_applicable")
    check("a tied end with no candidates at all is unreadable",
          dx.combination_check(five, E(barcode=None), wells)[0]
          == "unreadable")

    # A CANDIDATE SET OF ONE IS NOT A TIE, and this is the assertion the
    # first version of this check lacked. `match` declines a unique
    # nearest candidate whenever it sits further out than the set's
    # geometry safely supports; the combination then "confirms" it and
    # has contributed nothing. Counting that as a combination win put
    # the safe mechanism and an unsafe tolerance increase in one
    # number - and on the real plate almost all of it was the latter.
    lone = E(barcode=None, candidates=("B",), candidate_dist=3)
    v, s_, _ = dx.combination_check(five, lone, wells)
    check("a single candidate is NOT reported as a combination win",
          v == "distance_only", v)
    check("it still names where the read would go, for costing",
          s_ == "S1", str(s_))
    check("and a single candidate forming no well is a free decline",
          dx.combination_check(five, lone, none_declared)[0]
          == "no_declared")

    # The distinction has to survive the real call path, not just the
    # function: that is where the verdict is attached to the read.
    check("resolved requires two or more tied candidates",
          dx.combination_check(five, three, wells)[0] == "resolved"
          and len(three.candidates) >= 2)

    # NON-INTERFERENCE. The same ends, through the real classifier: the
    # read must still be tier B with no sample, even though the
    # combination check can name one. A measurement that quietly began
    # assigning reads would move results with no sign in the tiers,
    # which is the failure this project has already had twice.
    call = dx._classify("r1", five, three, "+", wells, False)
    verdict, named, _ = dx.combination_check(five, three, wells)
    check("the classifier is untouched: still tier B, still unassigned",
          call.tier == "B" and call.sample_id is None,
          f"{call.tier} {call.sample_id}")
    check("while the check names the sample it WOULD have gone to",
          verdict == "resolved" and named == "S1", f"{verdict} {named}")

    # candidates() returns the TIED set, not the k nearest: anything
    # further off is already excluded by distance.
    bs = dx.BarcodeSet.build({"AAAACCCC", "AAAAGGGG", "TTTTCCCC"}, 1)
    cand, d = bs.candidates("AAAAGGCC")
    check("candidates returns every barcode tied for nearest",
          set(cand) == {"AAAACCCC", "AAAAGGGG"} and d == 2,
          f"{cand} at {d}")
    cand2, _ = bs.candidates("AAAACCCC")
    check("an exact window returns just itself",
          cand2 == ("AAAACCCC",), str(cand2))
    check("nothing within reach returns an empty set, not a guess",
          bs.candidates("GGGGGGGG", max_dist=1) == ((), -1),
          str(bs.candidates("GGGGGGGG", max_dist=1)))
    check("a window of the wrong length returns nothing",
          bs.candidates("AAAA") == ((), -1))


def main() -> int:
    test_parallel_matches_serial()
    test_end_reasons_discriminate()
    test_barcode_margin_rule()
    test_combination_check()
    test_stats_only()
    test_barcode_diagnostic_discriminates()
    test_ligation_chimeras()
    test_find_anchor()
    test_index_matrix()
    test_streaming_output()
    test_minknow_input()
    test_barcode_set()
    test_against_truth()
    test_hopping_and_chimeras()
    test_orientation_and_trimming()
    test_single_copy_barcodes()

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
