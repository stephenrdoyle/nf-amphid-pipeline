#!/usr/bin/env python3
"""
Tests for PacBio demultiplexing (bin/demux_pacbio.py), the pacbio
platform adapter.

  * the matching rules of the reference script it ports (tag distance,
    clipped flanks, strictly-best tags, primers located through indels),
    each pinned by a case worked by hand;
  * every tier, on reads built for it, and the summary's accounting;
  * the PacBio fixture (tests/make_pacbio_fixture.py) read across a FASTQ
    and a BAM cell, one worker and several giving the same bytes;
  * the unaligned BAM reader;
  * lima: the barcodes and preset written for it, and its output read
    back as the samples it belongs to.

Run: python3 tests/test_demux_pacbio.py
"""

from __future__ import annotations

import contextlib
import gzip
import io
import json
import os
import random
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import demux_pacbio as dp          # noqa: E402
import demux_accounting as acc     # noqa: E402
import make_pacbio_fixture as pbf  # noqa: E402
from ampseqlib import open_reads, revcomp  # noqa: E402

PASS, FAIL = [], []

NC1 = "ACGTCTGGTTCAGGGTTGTT"
NC2 = "TTAGTTTCTTTTCCTCCGCT"
INSERT = ("GGCTAATGCAATTGCTACGTAGCATGCATTAGCTTAGCAATCGGATCGATTGCATGCAT"
          "CGATCGTAGCTAGCTAGGATCCGATCGATCGTTAGCAATGCATCG")
# Two pools' worth of 10-base symmetric tags, as S-BSST527 lays them out:
# forward and reverse agree on bases 1-9.
TAGS = [("A1", "GACATAGTCT", "GACATAGTCA"),
        ("A2", "GACGATCAGA", "GACGATCAGC"),
        ("A3", "GACGTGCTCT", "GACGTGCTCA"),
        ("A4", "GACTATAGCC", "GACTATAGCT")]


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def samples(tags=TAGS) -> list[dict]:
    return [{"run_id": "P1", "sample_id": s, "index_i5": f, "index_i7": r}
            for s, f, r in tags]


def matcher(tags=TAGS) -> dp.Matcher:
    return dp.Matcher(samples(tags), [NC1], [NC2], dict(dp.CONSTRUCT))


def read(f: str, r: str, insert: str = INSERT) -> str:
    return f + NC1 + insert + revcomp(NC2) + revcomp(r)


def quiet(fn, *a):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a)


# --------------------------------------------------------------------------


def test_rules() -> None:
    print("\nmatching rules")
    td = dp.tag_distance
    check("a tag read exactly is 0 edits", td("GACATAGTCT", "GACATAGTCT", False) == 0)
    check("bases before the tag in the flank are free",
          td("GACATAGTCT", "ACGACATAGTCT", False) == 0)
    check("a substitution is one edit", td("GACATAGTCT", "GACATAGACT", False) == 1)
    check("a deletion is one edit", td("GACATAGTCT", "ACGACATAGCT", False) == 1)
    check("the tag's end is held to the primer: a shifted tag costs",
          td("GACATAGTCT", "GACATAGTCTA", False) == 1)
    check("where the read ends in the flank, up to three missing tag "
          "bases are free", td("GACATAGTCT", "ATAGTCT", True) == 0
          and td("GACATAGTCT", "ATAGTCT", False) == 3)
    check("beyond three, the rest are charged",
          td("GACATAGTCT", "TAGTCT", True) == 1)

    m = matcher()
    end = m.end5
    check("an exact flank matches its tag", end.match("GACGATCAGA", True)[2]
          == "matched" and m.i5[end.match("GACGATCAGA", True)[0]] == "GACGATCAGA")
    k, d, why = end.match("GACGATCTCA", False)
    check("two edits from one tag and further from the rest: matched",
          why == "matched" and d == 2, f"{k} {d} {why}")
    tie = dp.TagEnd(["AAAAACCCCC", "AAAAGCCCCC"], 2)
    check("equally close to two tags: ambiguous, not either",
          tie.match("AAAATCCCCC", False)[2] == "ambiguous")
    check("further than tag_max_edits from every tag: no match",
          end.match("TTTTTTTTTT", False)[2] == "no_match")
    short = dp.TagEnd(["CTCTCTAT", "GTAAGGAG"], 2)
    check("a clipped flank holding under five tag bases names no tag "
          "(an 8-base tag at two edits matched on three)",
          short.match("TAT", True)[2] == "no_match"
          and short.match("CTCTAT", True)[2] == "matched")

    # Primers: exact, through an indel near their start, and degenerate.
    w = "GACATAGTCT" + NC1 + INSERT[:10]
    check("an exact primer is found where it is", dp.locate(w, [NC1], 3)
          == (10, 30, 0))
    w2 = "GACATAGTCT" + NC1[:2] + NC1[3:] + INSERT[:10]
    s, e, c = dp.locate(w2, [NC1], 3)
    check("a primer with a deletion: its true start, so the tag is read "
          "from the right place", (s, e, c) == (10, 29, 1), str((s, e, c)))
    check("a degenerate primer (IUPAC) matches",
          dp.locate("TT" + "AGCAGCCGCGGTAATTCC" + "GG", ["AGCAGCCGCGGTAAYTCC"], 2)
          == (2, 20, 0))


def test_tiers() -> None:
    print("\nevery tier, on reads built for it")
    m = matcher()
    f1, r1 = TAGS[0][1], TAGS[0][2]
    f2, r2 = TAGS[1][1], TAGS[1][2]
    c = m.classify(read(f1, r1))
    check("both tags of one sample: assigned, forward as read",
          c[0] == "assigned" and c[1] == "A1" and c[4] is False, str(c[:5]))
    c = m.classify(revcomp(read(f1, r1)))
    check("the same read reverse-complemented: assigned, turned",
          c[0] == "assigned" and c[1] == "A1" and c[4] is True)
    s = read(f1, r1)
    c = m.classify(s)
    check("what is written runs primer to primer, tags off",
          s[c[2]:c[3]] == NC1 + INSERT + revcomp(NC2))
    c = m.classify(read(f1, r2))
    check("tags of two samples: a pair no sample has",
          c[0] == "pair_not_declared" and c[1] is None)
    c = m.classify(read(f1, "TTTTTTTTTT"))
    check("one end a tag, the other nothing: one end only",
          c[0] == "one_end_only")
    c = m.classify(read("CCCCCCCCCC", "TTTTTTTTTT"))
    check("primers but no tag of this run: no tag", c[0] == "no_tag")
    c = m.classify(f1 + "ACGT" * 30 + revcomp(r1))
    check("no primers: no primers", c[0] == "no_primers")
    c = m.classify(f1 + NC1 + INSERT)
    check("forward primer only: no primers", c[0] == "no_primers")


def test_main(tmp: str) -> None:
    print("\nthe tool, on a library of every kind of read")
    rng = random.Random(3)
    lib = []
    want = {"A1": 0, "A2": 0, "A3": 0, "A4": 0}
    for k in range(200):
        sid, f, r = TAGS[k % 4]
        s = read(f, r)
        lib.append(s if k % 2 else revcomp(s))
        want[sid] += 1
    lib += [read(TAGS[0][1], TAGS[1][2])] * 7          # pair not declared
    lib += [read(TAGS[2][1], "AAAAAAAAAA")] * 5        # one end only
    lib += [read("CCCCCCCCCC", "AAAAAAAAAA")] * 4      # no tag
    lib += ["ACGT" * 40] * 3                           # no primers
    rng.shuffle(lib)
    fq = os.path.join(tmp, "lib.fastq.gz")
    with gzip.open(fq, "wt") as fh:
        for k, s in enumerate(lib):
            fh.write(f"@r{k}/ccs\n{s}\n+\n{'~' * len(s)}\n")
    ss = os.path.join(tmp, "samplesheet.tsv")
    with open(ss, "w") as fh:
        fh.write("run_id\tplatform\tsample_id\tindex_i7\tindex_i5\tsample_type\n")
        for sid, f, r in TAGS:
            fh.write(f"P1\tpacbio\t{sid}\t{r}\t{f}\tunknown\n")
        fh.write(f"P2\tpacbio\tother\t{TAGS[0][2]}\t{TAGS[0][1]}\tunknown\n")
    prim = os.path.join(tmp, "primers.tsv")
    with open(prim, "w") as fh:
        fh.write(f"name\tdirection\tsequence\nNC1\tfwd\t{NC1}\nNC2\trev\t{NC2}\n")
    out = os.path.join(tmp, "demux")
    code = quiet(dp.main, ["--reads", fq, "--samplesheet", ss, "--run-id", "P1",
                           "--primers", prim, "--outdir", out])
    js = json.load(open(os.path.join(out, "demux_summary.json")))
    check("runs", code == 0)
    check("every sample read to its sample, and only this run's samples",
          js["per_sample"] == want, str(js["per_sample"]))
    check("each fate counted", js["reasons"] == {
        "assigned": 200, "one_end_only": 5, "pair_not_declared": 7,
        "no_primers": 3, "no_tag": 4}, str(js["reasons"]))
    check("the tiers are nanopore's and add up to the reads",
          js["by_tier"] == {"A": 200, "B": 5, "C": 7, "D": 7}
          and js["read_accounting_balanced"] and js["total_reads"] == 219)
    check("the reads written are the reads assigned",
          js["reads_written_balanced"] is True and js["reads_written"] == 200)
    recs = list(open_reads(os.path.join(out, "reads", "A1.fastq.gz")))
    check("a sample's reads are written forward, primer to primer",
          len(recs) == 50 and all(s == NC1 + INSERT + revcomp(NC2)
                                  for _n, s, _q in recs))
    check("the undeclared pair lands on an unused cell of the matrix",
          js["index_matrix"]["reads_on_sentinels"] == 7)
    a = acc.accounting("P1", "pacbio", "single_long_accurate",
                       [{"run_id": "P1", "sample_id": s} for s, _f, _r in TAGS],
                       js, None)
    check("contract A: the shared fates, no internal_adapter, balanced",
          a["fates"] == {"assigned": 200, "one_end_only": 5, "chimeric": 7,
                         "unmatched": 7} and a["balanced"] is True
          and a["reads_in"] == 219 and a["reads_off_instrument"] == 219,
          str(a["fates"]))
    check("contract A: says amphid demultiplexed it",
          a["demultiplexed_by"] == "amphid (demux_pacbio.py)")

    bad = os.path.join(tmp, "bad.tsv")
    with open(bad, "w") as fh:
        fh.write("run_id\tplatform\tsample_id\tindex_i7\tindex_i5\tsample_type\n"
                 "P1\tpacbio\tx\tAAAA\t\tunknown\n")
    try:
        dp.run_samples(bad, "P1")
        msg = ""
    except SystemExit as exc:
        msg = str(exc)
    check("a sample without both tags is refused", "forward" in msg, msg)
    with open(bad, "w") as fh:
        fh.write("run_id\tplatform\tsample_id\tindex_i7\tindex_i5\tsample_type\n"
                 "P1\tpacbio\tx\tAAAA\tCCCC\tunknown\n"
                 "P1\tpacbio\ty\tAAAA\tCCCC\tunknown\n")
    try:
        dp.run_samples(bad, "P1")
        msg = ""
    except SystemExit as exc:
        msg = str(exc)
    check("a tag pair declared twice is refused", "more than one" in msg, msg)


def test_fixture(tmp: str) -> None:
    print("\nthe PacBio fixture: a FASTQ cell and a BAM cell, one library")
    fx = os.path.join(tmp, "fx")
    quiet(pbf.main, ["--outdir", fx])
    run = os.path.join(fx, "run")
    files = sorted(os.listdir(run))
    check("two cells, one FASTQ and one BAM",
          [f.rsplit(".", 2)[-2:] for f in files][0] == ["fastq", "gz"]
          and files[1].endswith(".bam"), str(files))
    bam = os.path.join(run, files[1])
    recs = list(dp.open_bam(bam))
    check("the BAM reader reads every record, sequence and quality",
          len(recs) > 800 and all(len(s) == len(q) and set(s) <= set("ACGT")
                                  for _n, s, q in recs))
    with open(bam, "rb") as fh:
        check("and the BAM is BGZF, ending in its EOF block",
              fh.read(4) == b"\x1f\x8b\x08\x04"
              and open(bam, "rb").read()[-28:] == pbf.BGZF_EOF)
    common = ["--reads", run, "--samplesheet", os.path.join(fx, "samplesheet.tsv"),
              "--run-id", pbf.RUN_ID, "--primers",
              os.path.join(ROOT, "tests", "data", "db", "primers.tsv")]
    outs = {}
    for t in (1, 3):
        outs[t] = os.path.join(tmp, f"fx_t{t}")
        quiet(dp.main, common + ["--threads", str(t), "--outdir", outs[t]])
    js = json.load(open(os.path.join(outs[1], "demux_summary.json")))
    truth = {}
    for line in open(os.path.join(fx, "truth.tsv")).read().splitlines()[1:]:
        sid, _src, n = line.split("\t")[:3]
        if sid:
            truth[sid] = truth.get(sid, 0) + int(n)
    got = js["per_sample"]
    check("both cells are read as one library",
          [f["format"] for f in js["input_files"]] == ["fastq", "bam"]
          and sum(f["reads"] for f in js["input_files"]) == js["total_reads"])
    check("at least 99% of every sample's reads reach it, and no sample "
          "gets more than it was given",
          all(truth[s] * 0.99 <= got[s] <= truth[s] for s in truth),
          f"{got} vs {truth}")
    same = open(os.path.join(outs[1], "demux_summary.json"), "rb").read() == \
        open(os.path.join(outs[3], "demux_summary.json"), "rb").read() and all(
        open(os.path.join(outs[1], "reads", f), "rb").read()
        == open(os.path.join(outs[3], "reads", f), "rb").read()
        for f in os.listdir(os.path.join(outs[1], "reads")))
    check("one worker or three: the same bytes", same)


def test_lima(tmp: str) -> None:
    print("\nlima: barcodes out, its output read back")

    def refused(tags) -> str:
        try:
            dp.lima_design(samples(tags))
            return ""
        except SystemExit as exc:
            return str(exc)

    # S-BSST527's layout: two tags one primer-adjacent base apart, and in
    # the same pool samples with one tag at both ends. ASYMMETRIC lost
    # every one of the latter and most of the former on the farm.
    pool = TAGS + [("A5", "GATCGTACAG", "GATCGTACAG")]
    design = dp.lima_design(samples(pool))
    check("tags one primer-adjacent base apart, with same-tag samples in "
          "the pool: SYMMETRIC", design["preset"] == "SYMMETRIC")
    check("one barcode per sample: its tag less the primer-adjacent base",
          [q for _n, q in design["barcodes"]]
          == sorted(f[:-1] for _s, f, _r in pool))
    name = {q: n for n, q in design["barcodes"]}
    check("each sample is its barcode at both ends",
          design["pairs"][(name[TAGS[0][1][:-1]],) * 2] == "A1"
          and len(design["pairs"]) == 5)
    sym = [("S1", "ACGTACGTAC", "ACGTACGTAC"), ("S2", "TGCATGCATG", "TGCATGCATG")]
    d0 = dp.lima_design(samples(sym))
    check("one tag at both ends for every sample: SYMMETRIC, whole tags",
          d0["preset"] == "SYMMETRIC"
          and [q for _n, q in d0["barcodes"]] == ["ACGTACGTAC", "TGCATGCATG"])
    combi = [("C1", "AAAACCCCGG", "TTTTGGGGAA"), ("C2", "AAAACCCCGG", "CCCCAAAATT"),
             ("C3", "GGGGTTTTCC", "TTTTGGGGAA")]
    d1 = dp.lima_design(samples(combi))
    check("distinct tags at the two ends: ASYMMETRIC, one barcode per tag",
          d1["preset"] == "ASYMMETRIC" and len(d1["barcodes"]) == 4)
    check("a pool mixing the two is refused, naming the builtin route",
          "builtin" in refused(combi + [("S1", "ACGTACGTAC", "ACGTACGTAC")]))
    check("samples whose barcodes collide once the last base is set aside "
          "are refused", "share the barcode" in refused(
              [("D1", "ACGTACGTAC", "ACGTACGTAG"),
               ("D2", "ACGTACGTAT", "ACGTACGTAT")]))
    # The read-back, on the run's own four S-BSST527-like samples.
    design = dp.lima_design(samples())
    name = {q: n for n, q in design["barcodes"]}
    n1 = name[TAGS[0][1][:-1]]
    x1 = name[TAGS[1][1][:-1]]

    d = os.path.join(tmp, "lima")
    os.makedirs(d)
    for f, k in ((f"part1.{n1}--{n1}.fastq", 6),     # A1, declared
                 (f"part1.{min(n1, x1)}--{max(n1, x1)}.fastq", 2)):  # undeclared
        with open(os.path.join(d, f), "w") as fh:
            for j in range(k):
                fh.write(f"@{f}{j}\n{NC1 + INSERT}\n+\n{'~' * (len(NC1) + len(INSERT))}\n")
    with open(os.path.join(d, "part1.lima.summary"), "w") as fh:
        fh.write("ZMWs input                    (A) : 20\n"
                 "ZMWs above all thresholds     (B) : 8 (40%)\n"
                 "ZMWs below any threshold      (C) : 12 (60%)\n\n"
                 "ZMW marginals for (C):\n"
                 "Below min length                  : 0 (0%)\n"
                 "Below min score                   : 9 (75%)\n"
                 "Below min end score               : 3 (25%)\n")
    open(os.path.join(d, "part1.lima.counts"), "w").write("IdxFirst\n")
    ss = os.path.join(tmp, "lima_ss.tsv")
    with open(ss, "w") as fh:
        fh.write("run_id\tplatform\tsample_id\tindex_i7\tindex_i5\tsample_type\n")
        for sid, f, r in TAGS:
            fh.write(f"P1\tpacbio\t{sid}\t{r}\t{f}\tunknown\n")
    prim = os.path.join(tmp, "primers.tsv")
    out = os.path.join(tmp, "lima_out")
    code = quiet(dp.main, ["--samplesheet", ss, "--run-id", "P1", "--primers",
                           prim, "--from-lima", d, "--outdir", out])
    js = json.load(open(os.path.join(out, "demux_summary.json")))
    check("read back: runs", code == 0)
    check("a declared pair's reads go to its sample",
          js["per_sample"]["A1"] == 6 and os.path.exists(
              os.path.join(out, "reads", "A1.fastq.gz")))
    check("an undeclared pair's reads are chimeric; lima's end-score "
          "rejects one end only; the rest unmatched; all of lima's input",
          js["reasons"] == {"assigned": 6, "one_end_only": 3,
                            "pair_not_declared": 2, "no_primers": 0,
                            "no_tag": 9} and js["total_reads"] == 20,
          str(js["reasons"]))
    check("lima's own summary rides along, and says it was lima",
          js["lima"]["summary"]["ZMWs input"] == 20
          and js["demultiplexed_by"].startswith("lima"))
    bc = os.path.join(tmp, "bc.fasta")
    quiet(dp.main, ["--samplesheet", ss, "--run-id", "P1", "--primers", prim,
                    "--write-lima-barcodes", bc])
    check("--write-lima-barcodes writes the FASTA and the preset",
          open(bc).read().count(">") == 4
          and open(bc + ".preset").read().strip() == "SYMMETRIC")


def main() -> int:
    tmp = tempfile.mkdtemp(prefix="demux_pacbio_")
    try:
        test_rules()
        test_tiers()
        test_main(tmp)
        test_fixture(tmp)
        test_lima(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{'=' * 62}")
    print(f"{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("\nfailures:")
        for name, detail in FAIL:
            print(f"  {name}: {detail}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
