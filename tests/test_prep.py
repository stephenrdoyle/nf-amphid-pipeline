#!/usr/bin/env python3
"""
Tests for Illumina read preparation.

The assertion that matters most is the one that is easy to skip: a merged
read must equal the *primer-trimmed* reference, not the reference. Reads
arrive with their primers removed, because those bases are synthesised
rather than observed. If the reference keeps its primers while the reads
lose theirs, every read looks 40bp short against every reference and
competitive assignment still returns a confident winner - it is simply
measuring the wrong thing, consistently, which is the failure mode least
likely to be noticed.

Two more things are checked against the simulator rather than against the
prep's own output:

  - reads with no primer must be the ones the simulator deliberately broke
    (dimer, off-target, truncated) and not clean reads
  - amplicons longer than the two reads together must land in the paired
    bucket rather than being merged on a spurious overlap

Run: python3 tests/test_prep.py
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
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


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


def qual(n: int, q: int = 37) -> str:
    return chr(33 + q) * n


def read_fasta(path: str) -> dict[str, str]:
    op = gzip.open if path.endswith(".gz") else open
    out, key = {}, None
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                key = ln[1:].split()[0]
                out[key] = ""
            elif key:
                out[key] += ln.strip()
    return out


# --------------------------------------------------------------------------
# primer matching
# --------------------------------------------------------------------------


def test_primer_match() -> None:
    print("\nprimer matching")
    ps = P.PrimerSet(
        fwd=[("F1", "ACGTACGTACGTACGTAAAA"), ("F2", "TTTTGGGGCCCCAAAATTTT")],
        rev=[("R1", "GGCCGGCCGGCCGGCCGGCC")],
    )
    body = "CATGCATGCATGCATGCATGCATG" * 6

    hit = ps.match("ACGTACGTACGTACGTAAAA" + body, "fwd")
    check("exact primer found and its end reported",
          hit == ("F1", 20), str(hit))

    # one substitution inside the primer
    hit = ps.match("ACGTACGTACGTACGTAAtA".upper() + body, "fwd")
    check("primer tolerates a substitution", hit is not None and hit[0] == "F1",
          str(hit))

    # a deletion: the end offset must move with it, or the trim is off by one
    hit = ps.match("ACGTACGTACGTACGTAAA" + body, "fwd")
    check("primer tolerates a deletion and the trim point moves",
          hit is not None and hit[1] == 19, str(hit))

    # the competing primer must not win
    hit = ps.match("TTTTGGGGCCCCAAAATTTT" + body, "fwd")
    check("the better-fitting primer wins", hit == ("F2", 20), str(hit))

    # a primer sequence in the middle of a read is not the primer
    hit = ps.match(body + "ACGTACGTACGTACGTAAAA" + body, "fwd")
    check("an internal primer match is rejected", hit is None, str(hit))

    # too many errors
    hit = ps.match("AAAAAAAAAAAAAAAAAAAA" + body, "fwd")
    check("an unrelated read start matches nothing", hit is None, str(hit))


# --------------------------------------------------------------------------
# read-through and merging
# --------------------------------------------------------------------------


def test_trim_and_merge() -> None:
    print("\nread-through trimming and merging")
    rng = random.Random(4)
    insert = "".join(rng.choice("ACGT") for _ in range(180))
    adapter = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

    s, q = P.trim_readthrough(insert + adapter, qual(len(insert) + len(adapter)),
                              adapter)
    check("adapter read-through removed exactly", s == insert,
          f"{len(s)} vs {len(insert)}")
    check("quality string trimmed with the sequence", len(q) == len(s))

    s, _ = P.trim_readthrough(insert, qual(len(insert)), adapter)
    check("a read with no read-through is left alone", s == insert)

    # merging: 250bp mates over a 300bp insert -> 200bp overlap
    amp = "".join(rng.choice("ACGT") for _ in range(300))
    r1 = amp[:250]
    r2 = P.revcomp(amp[50:])
    m = P.merge_pair(r1, qual(250), r2, qual(250))
    check("overlapping mates merge to the original insert",
          m is not None and m[0] == amp,
          "None" if m is None else f"{len(m[0])} vs {len(amp)}")

    # the mates disagree at one base; the higher quality mate should win
    bad = list(r2)
    pos = 100
    bad[pos] = "A" if bad[pos] != "A" else "C"
    q2 = list(qual(250, 37))
    q2[pos] = chr(33 + 8)
    m = P.merge_pair(r1, qual(250, 37), "".join(bad), "".join(q2))
    check("a disagreement is resolved toward the higher quality base",
          m is not None and m[0] == amp,
          "None" if m is None else f"differs at {sum(1 for a,b in zip(m[0],amp) if a!=b)}")

    # and the merged quality at a contested base must not inherit the winner's
    if m is not None:
        contested = len(amp) - 1 - pos
        check("a contested base loses confidence",
              ord(m[1][contested]) - 33 < 37,
              str(ord(m[1][contested]) - 33))

    # non-overlapping mates
    long_amp = "".join(rng.choice("ACGT") for _ in range(700))
    m = P.merge_pair(long_amp[:250], qual(250),
                     P.revcomp(long_amp[450:]), qual(250))
    check("mates that do not overlap are not merged", m is None,
          "" if m is None else f"merged to {len(m[0])}")

    # mates that overlap but disagree everywhere
    other = "".join(rng.choice("ACGT") for _ in range(300))
    m = P.merge_pair(amp[:250], qual(250), P.revcomp(other[50:]), qual(250))
    check("mates from different molecules are not merged", m is None,
          "" if m is None else "merged")


# --------------------------------------------------------------------------
# end to end against the truth table
# --------------------------------------------------------------------------


def build_run(tmp: str, n_samples=12, reads=250):
    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, n_samples, 2, "RUN01", "illumina")
    mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
    mf.write_primers(os.path.join(tmp, "primers.tsv"))
    mf.write_samplesheet(samples, os.path.join(tmp, "samplesheet.csv"))
    mf.write_runs(os.path.join(tmp, "runs.tsv"), "RUN01", "illumina",
                  os.path.join(tmp, "sim", "RUN01"))

    with contextlib.redirect_stdout(io.StringIO()):
        pf.main([
            "--reference", os.path.join(tmp, "reference.tsv"),
            "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
            "--runs", os.path.join(tmp, "runs.tsv"),
            "--primers", os.path.join(tmp, "primers.tsv"),
            "--index-map", INDEX_MAP,
            "--assay", ASSAY,
            "--outdir", os.path.join(tmp, "preflight"),
        ])
        sim.main([
            "--reference", os.path.join(tmp, "reference.tsv"),
            "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
            "--runs", os.path.join(tmp, "runs.tsv"),
            "--assay", ASSAY,
            "--outdir", os.path.join(tmp, "sim"),
            "--reads-per-sample", str(reads),
        ])
        for r1 in sorted(glob.glob(os.path.join(tmp, "sim", "RUN01",
                                                "*_R1_001.fastq.gz"))):
            sid = os.path.basename(r1).split("_S")[0]
            P.main([
                "--r1", r1, "--r2", r1.replace("_R1_", "_R2_"),
                "--sample-id", sid,
                "--primers", os.path.join(tmp, "primers.tsv"),
                "--assay", ASSAY,
                "--outdir", os.path.join(tmp, "prep"),
                "--quiet",
            ])

    with open(os.path.join(tmp, "sim", "RUN01.truth.tsv"), encoding="utf-8") as fh:
        truth = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return truth


def test_against_truth() -> None:
    print("\nprep against the truth table")
    tmp = tempfile.mkdtemp()
    try:
        truth = build_run(tmp)

        full = read_fasta(os.path.join(tmp, "preflight", "reference.fasta"))
        insert = read_fasta(os.path.join(tmp, "preflight",
                                         "reference.insert.fasta"))
        check("preflight writes a primer-trimmed panel",
              bool(insert) and len(insert) == len(full))
        trimmed = [len(full[k]) - len(insert[k]) for k in insert
                   if len(full[k]) != len(insert[k])]
        check("the trim removes both primers and nothing else",
              trimmed and set(trimmed) == {40}, str(sorted(set(trimmed))))

        # accounting, across every sample
        summaries = [json.load(open(f))
                     for f in glob.glob(os.path.join(tmp, "prep", "*.prep.json"))]
        check("every sample was prepped", len(summaries) > 1, str(len(summaries)))
        check("read accounting balances in every sample",
              all(s["accounting_balances"] for s in summaries))
        total_in = sum(s["pairs_in"] for s in summaries)
        total_out = sum(sum(s["fates"].values()) for s in summaries)
        check("no pair vanishes between input and fates",
              total_in == total_out, f"{total_in} vs {total_out}")

        # rerun one sample in-process so reads keep their identity
        cfg = P.load_assay(ASSAY)
        ps = P.load_primers(os.path.join(tmp, "primers.tsv"))
        rt1, rt2 = P.revcomp(cfg["pcr1_overhang_rev"]), P.revcomp(
            cfg["pcr1_overhang_fwd"])
        r1s = sorted(glob.glob(os.path.join(tmp, "sim", "RUN01",
                                            "*_R1_001.fastq.gz")))

        outcome = Counter()
        merged_exact = Counter()
        by_amplicon = {"merged": Counter(), "paired": Counter()}
        for r1p in r1s[:6]:
            pairs = zip(P.open_reads(r1p),
                        P.open_reads(r1p.replace("_R1_", "_R2_")))
            for (rid, s1, q1), (_rid2, s2, q2) in pairs:
                t = truth[rid]
                ref = insert.get(t["seq_id"])
                f, r = ps.match(s1, "fwd"), ps.match(s2, "rev")
                if f is None or r is None:
                    outcome[("no_primer", t["fate"])] += 1
                    continue
                a1, b1 = P.trim_readthrough(s1[f[1]:], q1[f[1]:], rt1,
                                            dict(ps.rev)[r[0]])
                a2, b2 = P.trim_readthrough(s2[r[1]:], q2[r[1]:], rt2,
                                            dict(ps.fwd)[f[0]])
                m = P.merge_pair(a1, b1, a2, b2)
                kind = "merged" if m is not None else "paired"
                outcome[(kind, t["fate"])] += 1
                if ref and t["fate"] == "clean":
                    by_amplicon[kind][len(ref)] += 1
                    if kind == "merged":
                        merged_exact["match" if m[0] == ref else "differs"] += 1
                        merged_exact["right length" if len(m[0]) == len(ref)
                                     else "wrong length"] += 1

        clean_merged = outcome[("merged", "clean")]
        clean_total = sum(v for (k, f), v in outcome.items() if f == "clean")
        check("most clean pairs merge",
              clean_merged / clean_total > 0.85,
              f"{clean_merged}/{clean_total}")

        # The assertion this suite exists for.
        acc = merged_exact["match"] / max(1, merged_exact["match"]
                                          + merged_exact["differs"])
        check("merged reads reproduce the primer-trimmed reference exactly "
              "more often than not", acc > 0.5, f"{acc:.3f}")
        lok = merged_exact["right length"] / max(
            1, merged_exact["right length"] + merged_exact["wrong length"])
        check("merged length matches the trimmed reference",
              lok > 0.97, f"{lok:.4f}")

        # no clean read should be silently discarded for want of a primer
        np_clean = outcome[("no_primer", "clean")]
        np_total = sum(v for (k, f), v in outcome.items() if k == "no_primer")
        check("reads without a primer are the deliberately broken ones",
              np_total > 0 and np_clean / np_total < 0.02,
              f"{np_clean} clean of {np_total}")

        # long amplicons cannot merge and must not pretend to
        long_ref = [len(v) for v in insert.values() if len(v) > 500]
        if long_ref and by_amplicon["paired"]:
            merged_long = sum(n for L, n in by_amplicon["merged"].items()
                              if L > 500)
            paired_long = sum(n for L, n in by_amplicon["paired"].items()
                              if L > 500)
            check("amplicons longer than the read pair stay unmerged",
                  paired_long > 0 and merged_long == 0,
                  f"{merged_long} merged, {paired_long} paired")
            short_paired = sum(n for L, n in by_amplicon["paired"].items()
                               if L < 450)
            short_merged = sum(n for L, n in by_amplicon["merged"].items()
                               if L < 450)
            check("short amplicons overwhelmingly merge",
                  short_merged / max(1, short_merged + short_paired) > 0.95,
                  f"{short_merged} merged, {short_paired} paired")

        # dereplication must not lose or invent reads
        sid = os.path.basename(r1s[0]).split("_S")[0]
        s = [x for x in summaries if x["sample_id"] == sid][0]
        fa = read_fasta(os.path.join(tmp, "prep", f"{sid}.merged.fasta.gz"))
        sizes = sum(int(h.split("size=")[1]) for h in
                    _headers(os.path.join(tmp, "prep",
                                          f"{sid}.merged.fasta.gz")))
        check("dereplicated counts sum to the merged reads",
              sizes == s["fates"]["merged"],
              f"{sizes} vs {s['fates']['merged']}")
        check("unique count matches the record count",
              len(fa) == s["unique_merged"], f"{len(fa)} vs {s['unique_merged']}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _headers(path: str):
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                yield ln.strip()


# --------------------------------------------------------------------------
# filters that the simulator is too kind to trigger
# --------------------------------------------------------------------------


def test_filters() -> None:
    print("\nquality and length filters")
    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(9)
        fwd = "ACGTACGTACGTACGTAAAA"
        rev = "GGCCGGCCGGCCGGCCGGCC"
        with open(os.path.join(tmp, "primers.tsv"), "w") as fh:
            fh.write("name\tdirection\tsequence\n")
            fh.write(f"F1\tfwd\t{fwd}\nR1\trev\t{rev}\n")

        amp = "".join(rng.choice("ACGT") for _ in range(300))
        mol = fwd + amp + P.revcomp(rev)
        r1, r2 = mol[:250], P.revcomp(mol)[:250]

        def write(path, recs):
            with gzip.open(path, "wt") as fh:
                for rid, s, q in recs:
                    fh.write(f"@{rid}\n{s}\n+\n{q}\n")

        good = ("g1", r1, qual(len(r1), 37)), ("g2", r1, qual(len(r1), 37))
        goodr = ("g1", r2, qual(len(r2), 37)), ("g2", r2, qual(len(r2), 37))
        bad = ("b1", r1, qual(len(r1), 5)),
        badr = ("b1", r2, qual(len(r2), 5)),

        write(os.path.join(tmp, "a_R1.fastq.gz"), list(good) + list(bad))
        write(os.path.join(tmp, "a_R2.fastq.gz"), list(goodr) + list(badr))

        with contextlib.redirect_stdout(io.StringIO()):
            P.main(["--r1", os.path.join(tmp, "a_R1.fastq.gz"),
                    "--r2", os.path.join(tmp, "a_R2.fastq.gz"),
                    "--sample-id", "a",
                    "--primers", os.path.join(tmp, "primers.tsv"),
                    "--outdir", os.path.join(tmp, "prep"), "--quiet"])
        s = json.load(open(os.path.join(tmp, "prep", "a.prep.json")))
        check("a low quality pair is rejected", s["fates"]["low_quality"] == 1,
              str(s["fates"]))
        check("good pairs survive the filter", s["fates"]["merged"] == 2,
              str(s["fates"]))

        # swapped mates: the same molecule, handed over the wrong way round
        write(os.path.join(tmp, "b_R1.fastq.gz"), list(goodr))
        write(os.path.join(tmp, "b_R2.fastq.gz"), list(good))
        with contextlib.redirect_stdout(io.StringIO()):
            P.main(["--r1", os.path.join(tmp, "b_R1.fastq.gz"),
                    "--r2", os.path.join(tmp, "b_R2.fastq.gz"),
                    "--sample-id", "b",
                    "--primers", os.path.join(tmp, "primers.tsv"),
                    "--outdir", os.path.join(tmp, "prep"), "--quiet"])
        s = json.load(open(os.path.join(tmp, "prep", "b.prep.json")))
        check("mates supplied the wrong way round are swapped, not dropped",
              s["fates"]["merged"] == 2, str(s["fates"]))

        # files out of step must raise rather than silently mispair
        write(os.path.join(tmp, "c_R1.fastq.gz"), list(good))
        write(os.path.join(tmp, "c_R2.fastq.gz"), list(goodr)[::-1])
        raised = False
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                P.main(["--r1", os.path.join(tmp, "c_R1.fastq.gz"),
                        "--r2", os.path.join(tmp, "c_R2.fastq.gz"),
                        "--sample-id", "c",
                        "--primers", os.path.join(tmp, "primers.tsv"),
                        "--outdir", os.path.join(tmp, "prep"), "--quiet"])
        except SystemExit:
            raised = True
        check("R1 and R2 out of step is fatal, not silent", raised)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------


def test_primer_collision() -> None:
    """Two references differing only under a primer are indistinguishable."""
    print("\nreferences that differ only under a primer")
    tmp = tempfile.mkdtemp()
    try:
        fwd = "ACGTACGTACGTACGTAAAA"
        rev = "GGCCGGCCGGCCGGCCGGCC"
        core = "".join(random.Random(2).choice("ACGT") for _ in range(300))
        a = fwd + core + P.revcomp(rev)
        # same insert, one base different inside the forward primer site
        b = "TCGTACGTACGTACGTAAAA" + core + P.revcomp(rev)

        with open(os.path.join(tmp, "primers.tsv"), "w") as fh:
            fh.write(f"name\tdirection\tsequence\nF1\tfwd\t{fwd}\nR1\trev\t{rev}\n")
        with open(os.path.join(tmp, "reference.tsv"), "w") as fh:
            fh.write("seq_id\tlocus\ttype\tpurpose\ttaxon\tparent_id\t"
                     "default_conc\tconc_units\tcorr_factor\tsequence\n")
            fh.write(f"A_h1\tITS2\ttarget\tspeciation\tSpecies_a\t\t\t\t\t{a}\n")
            fh.write(f"B_h1\tITS2\ttarget\tspeciation\tSpecies_b\t\t\t\t\t{b}\n")

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = pf.main(["--reference", os.path.join(tmp, "reference.tsv"),
                          "--primers", os.path.join(tmp, "primers.tsv"),
                          "--assay", ASSAY,
                          "--outdir", os.path.join(tmp, "pf")])
        log = out.getvalue()
        check("a pair separable only under a primer is an error",
              rc != 0 and "only under a primer" in log,
              log.strip().splitlines()[-1] if log.strip() else "no output")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# short inserts: read-through crosses the opposite primer first
# --------------------------------------------------------------------------

# The real construct, so the junction is the one a MiSeq actually reads.
OVF = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
OVR = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
FWD = "GCGTTATTTACCACGAATTGCAG"
REV = "TATGCTTAAGTTCAGCGGGT"
# What lies past each overhang (index read primer site and beyond); only
# its being there matters.
P7_TAIL = "ATCTCGTATGCCGTCTTCTGCTTG"
P5_TAIL = "GTGTAGATCTCGGTGGTCGCCGTATCATT"


def short_pair(insert: str, read_len: int = 250):
    """R1 and R2 as sequenced from fwd - insert - rc(rev), primers first."""
    r1 = (FWD + insert + P.revcomp(REV) + P.revcomp(OVR) + P7_TAIL)[:read_len]
    r2 = (REV + P.revcomp(insert) + P.revcomp(FWD) + P.revcomp(OVF)
          + P5_TAIL)[:read_len]
    return r1, r2


def test_readthrough_takes_the_opposite_primer() -> None:
    """Between a short insert and the adapter sits the opposite primer.

    Trimming at the adapter alone left rc(primer) on every read-through
    mate: the merge seed then sat on primer, the pair stayed unmerged, and
    assignment charged each mate ~20 edits for synthesised bases.
    """
    print("\nread-through removes the opposite primer too")
    rng = random.Random(31)
    rt1, rt2 = P.revcomp(OVR), P.revcomp(OVF)

    insert = "".join(rng.choice("ACGT") for _ in range(150))
    r1, r2 = short_pair(insert)
    s1, s2 = r1[len(FWD):], r2[len(REV):]
    a, qa = P.trim_readthrough(s1, qual(len(s1)), rt1, REV)
    b, qb = P.trim_readthrough(s2, qual(len(s2)), rt2, FWD)
    check("R1 read-through is cut where the insert ends",
          a == insert, f"{len(a)} vs {len(insert)}")
    check("R2 read-through is cut where the insert ends",
          b == P.revcomp(insert), f"{len(b)} vs {len(insert)}")
    m = P.merge_pair(a, qa, b, qb)
    check("a short insert merges to exactly the insert",
          m is not None and m[0] == insert,
          "None" if m is None else f"{len(m[0])} vs {len(insert)}")

    # The read stops a few bases into the adapter: the full junction is
    # not on the read, only its start.
    for past in (12, 5, 0):
        ins = "".join(rng.choice("ACGT") for _ in range(250 - len(FWD)
                                                     - len(REV) - past))
        r1, _ = short_pair(ins)
        s1 = r1[len(FWD):]
        a, _ = P.trim_readthrough(s1, qual(len(s1)), rt1, REV)
        check(f"a read ending {past}bp into the adapter loses all of the "
              "junction", a == ins, f"{len(a)} vs {len(ins)}")

    # ... or inside the primer itself, with only the primer's first bases
    ins = "".join(rng.choice("ACGT") for _ in range(250 - len(FWD) - 12))
    r1, _ = short_pair(ins)
    s1 = r1[len(FWD):]
    a, _ = P.trim_readthrough(s1, qual(len(s1)), rt1, REV)
    check("a read ending 12bp into the opposite primer loses those 12bp",
          a == ins, f"{len(a)} vs {len(ins)}")

    # a degenerate primer has to be reverse-complemented code by code
    deg = REV[:5] + "Y" + REV[6:]            # REV[5] is T, so Y covers it
    a, _ = P.trim_readthrough(short_pair(insert)[0][len(FWD):],
                              qual(250 - len(FWD)), rt1, deg)
    check("a degenerate opposite primer is still recognised",
          a == insert, f"{len(a)} vs {len(insert)}")

    # and a read that never ran through must keep every base
    trimmed = 0
    for _ in range(3000):
        s = "".join(rng.choice("ACGT") for _ in range(227))
        a, _ = P.trim_readthrough(s, qual(len(s)), rt1, REV)
        trimmed += a != s
    check("reads with no read-through are not trimmed by chance",
          trimmed == 0, f"{trimmed} of 3000 trimmed")

    # End to end through prep_sample: a short-amplicon library must merge.
    tmp = tempfile.mkdtemp()
    try:
        with open(os.path.join(tmp, "primers.tsv"), "w") as fh:
            fh.write(f"name\tdirection\tsequence\nF\tfwd\t{FWD}\n"
                     f"R\trev\t{REV}\n")
        recs1, recs2 = [], []
        for i, L in enumerate((120, 160, 190, 205, 220, 225)):
            ins = "".join(rng.choice("ACGT") for _ in range(L))
            r1, r2 = short_pair(ins)
            recs1.append((f"s{i}", r1, qual(len(r1))))
            recs2.append((f"s{i}", r2, qual(len(r2))))
        _write_fastq(os.path.join(tmp, "R1.fastq.gz"), recs1)
        _write_fastq(os.path.join(tmp, "R2.fastq.gz"), recs2)
        with contextlib.redirect_stdout(io.StringIO()):
            P.main(["--r1", os.path.join(tmp, "R1.fastq.gz"),
                    "--r2", os.path.join(tmp, "R2.fastq.gz"),
                    "--sample-id", "short", "--min-length", "50",
                    "--primers", os.path.join(tmp, "primers.tsv"),
                    "--outdir", tmp, "--quiet"])
        sj = json.load(open(os.path.join(tmp, "short.prep.json")))
        fa = read_fasta(os.path.join(tmp, "short.merged.fasta.gz"))
        lens = sorted(len(v) for v in fa.values())
        check("every short-amplicon pair merges in prep",
              sj["fates"]["merged"] == 6, str(sj["fates"]))
        check("and each merged read is exactly its insert",
              lens == [120, 160, 190, 205, 220, 225], str(lens))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _write_fastq(path, recs):
    with gzip.open(path, "wt") as fh:
        for rid, s, q in recs:
            fh.write(f"@{rid}\n{s}\n+\n{q}\n")


def test_min_overlap_is_honoured() -> None:
    """--min-overlap 20 must merge a 20bp overlap.

    The 24bp seed has to lie inside R1, so overlaps of 20-23bp were never
    found and the effective minimum was 24 whatever the option said.
    """
    print("\nshort overlaps down to --min-overlap merge")
    rng = random.Random(41)
    for ov in (20, 21, 22, 23):
        amp = "".join(rng.choice("ACGT") for _ in range(500 - ov))
        m = P.merge_pair(amp[:250], qual(250), P.revcomp(amp[250 - ov:]),
                         qual(250), min_overlap=20)
        check(f"a {ov}bp overlap merges to the amplicon",
              m is not None and m[0] == amp,
              "None" if m is None else f"{len(m[0])} vs {len(amp)}")
    amp = "".join(rng.choice("ACGT") for _ in range(481))
    m = P.merge_pair(amp[:250], qual(250), P.revcomp(amp[231:]), qual(250),
                     min_overlap=20)
    check("a 19bp overlap stays below --min-overlap 20", m is None)
    m = P.merge_pair(amp[:250], qual(250), P.revcomp(amp[231:]), qual(250),
                     min_overlap=19)
    check("and merges once --min-overlap allows it",
          m is not None and m[0] == amp)
    # one error inside a 20bp overlap is still the same molecule
    amp = "".join(rng.choice("ACGT") for _ in range(480))
    r2 = list(P.revcomp(amp[230:]))
    r2[-5] = "A" if r2[-5] != "A" else "C"
    m = P.merge_pair(amp[:250], qual(250), "".join(r2), qual(250, 20))
    check("a 20bp overlap with one error merges toward the better base",
          m is not None and m[0] == amp,
          "None" if m is None else f"{len(m[0])}")


def test_truncated_mate_is_fatal() -> None:
    """zip() stopped at the shorter file and the accounting still balanced."""
    print("\nmate files that end unequally")
    tmp = tempfile.mkdtemp()
    try:
        fwd, rev = "ACGTACGTACGTACGTAAAA", "GGCCGGCCGGCCGGCCGGCC"
        with open(os.path.join(tmp, "primers.tsv"), "w") as fh:
            fh.write(f"name\tdirection\tsequence\nF1\tfwd\t{fwd}\n"
                     f"R1\trev\t{rev}\n")
        amp = "".join(random.Random(5).choice("ACGT") for _ in range(300))
        mol = fwd + amp + P.revcomp(rev)
        r1, r2 = mol[:250], P.revcomp(mol)[:250]
        recs1 = [(f"x{i}", r1, qual(250)) for i in range(5)]
        recs2 = [(f"x{i}", r2, qual(250)) for i in range(5)]

        def run(sid, a, b):
            _write_fastq(os.path.join(tmp, f"{sid}_R1.fastq.gz"), a)
            _write_fastq(os.path.join(tmp, f"{sid}_R2.fastq.gz"), b)
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    P.main(["--r1", os.path.join(tmp, f"{sid}_R1.fastq.gz"),
                            "--r2", os.path.join(tmp, f"{sid}_R2.fastq.gz"),
                            "--sample-id", sid,
                            "--primers", os.path.join(tmp, "primers.tsv"),
                            "--outdir", tmp, "--quiet"])
            except SystemExit as e:
                return str(e)
            return None

        err = run("t2", recs1, recs2[:3])
        check("a truncated R2 is fatal", err is not None and "ran out" in err,
              str(err))
        err = run("t1", recs1[:3], recs2)
        check("a truncated R1 is fatal", err is not None and "ran out" in err,
              str(err))
        err = run("ok", recs1, recs2)
        check("complete mate files still run", err is None, str(err))
        # pre-Casava 1.8 ids carry the mate number on the id itself
        err = run("old", [(f"x{i}/1", r1, qual(250)) for i in range(5)],
                  [(f"x{i}/2", r2, qual(250)) for i in range(5)])
        check("old-style /1 and /2 mate suffixes are mates, not a mismatch",
              err is None, str(err))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_primer_match()
    test_trim_and_merge()
    test_readthrough_takes_the_opposite_primer()
    test_min_overlap_is_honoured()
    test_truncated_mate_is_fatal()
    test_filters()
    test_primer_collision()
    test_against_truth()

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
