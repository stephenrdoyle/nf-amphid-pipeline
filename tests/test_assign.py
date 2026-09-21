#!/usr/bin/env python3
"""
Tests for competitive assignment.

The bug this suite was written around is worth stating, because it produced
a perfectly plausible answer. A control plasmid differs from its target only
by an inserted tag. On an amplicon too long for the mates to merge, that tag
can fall in the gap between the two reads - and then the pair fits the
target and the control equally well. The first implementation treated a
missing tag as evidence of target, so every carryover read at that locus was
reported as a genuine detection: 97 of them, silently, with nothing in the
output suggesting a problem.

The fix is asymmetric, and the asymmetry is the point. Finding the tag
proves control. Not finding it proves nothing unless the read covered the
place the tag would be - and a tie between a target and its control is
itself proof that it did not, since an observed tag would have separated
them by its own length.

bounded_edit is cross-checked against a plain dynamic program rather than
against itself. An earlier matcher in this pipeline passed a unit test that
could not fail, and shipped a bug that halved assignment.

Run: python3 tests/test_assign.py
"""

from __future__ import annotations

import contextlib
import csv
import glob
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


A = _load("assign", "bin/assign.py")
P = _load("prep_illumina", "bin/prep_illumina.py")
pf = _load("preflight", "bin/preflight.py")
sim = _load("simulate", "bin/simulate.py")
sys.path.insert(0, os.path.join(ROOT, "bin"))
from ampseqlib import bounded_edit, revcomp  # noqa: E402

sys.path.insert(0, HERE)
import make_fixtures as mf  # noqa: E402

INDEX_MAP = os.path.join(ROOT, "assets", "index_map.tsv")
ASSAY = os.path.join(ROOT, "assets", "assay", "ampseq_helminth.json")
TAG = json.load(open(ASSAY))["control_tag"]["sequence"]

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def plain_dp(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, y in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y))
        prev = cur
    return prev[-1]


def mutate(rng, seq: str, n: int) -> str:
    s = list(seq)
    for _ in range(n):
        p = rng.randrange(len(s))
        op = rng.choice("sid")
        if op == "s":
            s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
        elif op == "i":
            s.insert(p, rng.choice("ACGT"))
        elif len(s) > 1:
            s.pop(p)
    return "".join(s)


# --------------------------------------------------------------------------
# the distance primitive
# --------------------------------------------------------------------------


def test_bounded_edit() -> None:
    print("\nbounded edit distance")
    rng = random.Random(11)
    mismatches = 0
    n = 0
    for _ in range(1500):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 60)))
        b = mutate(rng, a, rng.randrange(0, 8))
        true = plain_dp(a, b)
        for md in (0, 1, 2, 3, 6, 15, 60):
            want = true if true <= md else md + 1
            if bounded_edit(a, b, md) != want:
                mismatches += 1
            n += 1
    check("agrees with a plain dynamic program at every tolerance",
          mismatches == 0, f"{mismatches} of {n} disagreed")

    # the free-suffix mode used for the mates of an unmerged pair
    bad = 0
    for _ in range(600):
        full = "".join(rng.choice("ACGT") for _ in range(rng.randrange(40, 120)))
        cut = rng.randrange(10, len(full))
        mate = mutate(rng, full[:cut], rng.randrange(0, 4))
        d = bounded_edit(mate, full, 10, free_b_suffix=True)
        if d > 4:
            bad += 1
    check("a mate is not charged for the reference it never reached",
          bad == 0, f"{bad} mates scored as distant")

    a = "ACGT" * 20
    check("an unrelated sequence exceeds the bound",
          bounded_edit(a, "TTTT" * 20, 5) == 6)
    check("identical sequences are distance zero", bounded_edit(a, a, 5) == 0)


# --------------------------------------------------------------------------
# calling
# --------------------------------------------------------------------------


def toy_panel():
    rng = random.Random(3)
    core = "".join(rng.choice("ACGT") for _ in range(360))
    h1 = core
    h2 = mutate(random.Random(5), core, 0)
    # three substitutions, the separation real haplotypes have
    h2 = list(core)
    for p in (80, 190, 300):
        h2[p] = "A" if h2[p] != "A" else "G"
    h2 = "".join(h2)
    other = "".join(random.Random(9).choice("ACGT") for _ in range(360))

    refs = [
        A.Ref("Sp_a_h1", h1, "ITS2", "target", "Species_a"),
        A.Ref("Sp_a_h2", h2, "ITS2", "target", "Species_a"),
        A.Ref("Sp_b_h1", other, "ITS2", "target", "Species_b"),
        # control: parent with the tag inserted 30bp from the end
        A.Ref("Sp_a_h1_ctrl", h1[:330] + TAG + h1[330:], "ITS2", "control",
              "Species_a", parent_id="Sp_a_h1"),
        # and one with the tag in the middle, the pathological placement
        A.Ref("Sp_b_h1_ctrl", other[:180] + TAG + other[180:], "ITS2",
              "control", "Species_b", parent_id="Sp_b_h1"),
    ]
    return A.Panel(refs), {r.seq_id: r for r in refs}


def test_calls() -> None:
    print("\ncalling")
    panel, by_id = toy_panel()
    h1 = by_id["Sp_a_h1"].sequence
    h2 = by_id["Sp_a_h2"].sequence
    ctrl = by_id["Sp_a_h1_ctrl"].sequence

    c = A.call_sequence(h1, panel, TAG, 0.10, 2, 3)
    check("an exact read is assigned to its reference",
          c.outcome == "assigned" and c.seq_id == "Sp_a_h1", str(c))

    rng = random.Random(21)
    c = A.call_sequence(mutate(rng, h1, 2), panel, TAG, 0.10, 2, 3)
    check("two errors do not lose the haplotype",
          c.outcome == "assigned" and c.seq_id == "Sp_a_h1", str(c))

    c = A.call_sequence(h2, panel, TAG, 0.10, 2, 3)
    check("haplotypes three bases apart are separated",
          c.outcome == "assigned" and c.seq_id == "Sp_a_h2", str(c))

    # a read midway between the two haplotypes identifies neither
    mid = list(h1)
    mid[80] = h2[80]
    c = A.call_sequence("".join(mid), panel, TAG, 0.10, 2, 3)
    check("a read between two haplotypes is not forced onto one",
          c.outcome == "taxon_only" and c.taxon == "Species_a", str(c))
    check("the species survives when the haplotype does not",
          c.outcome == "taxon_only" and set(c.tied) == {"Sp_a_h1", "Sp_a_h2"},
          str(c.tied))

    c = A.call_sequence(ctrl, panel, TAG, 0.10, 2, 3)
    check("a tagged read is called control, not target",
          c.outcome == "assigned" and c.seq_id == "Sp_a_h1_ctrl", str(c))
    check("the tag is reported as found", c.tag)

    c = A.call_sequence(h1, panel, TAG, 0.10, 2, 3)
    check("an untagged read that covers the tag site is called target",
          c.outcome == "assigned" and c.seq_id == "Sp_a_h1", str(c))

    foreign = "".join(random.Random(77).choice("ACGT") for _ in range(360))
    c = A.call_sequence(foreign, panel, TAG, 0.10, 2, 3)
    check("a sequence from outside the panel matches nothing",
          c.outcome == "unmatched", str(c))


def test_tag_unobservable() -> None:
    """The bug this suite exists for."""
    print("\nwhen the tag falls between the mates")
    panel, by_id = toy_panel()
    target = by_id["Sp_b_h1"].sequence          # tag site at position 180
    control = by_id["Sp_b_h1_ctrl"].sequence

    # mates of 100bp: they see [0,100) and the last 100bp, never position 180
    for name, src in (("a target read", target), ("a control read", control)):
        r1 = src[:100]
        r2 = revcomp(src[-100:])
        c = A.call_sequence(r1, panel, TAG, 0.10, 2, 3, second=r2)
        check(f"{name} whose tag site is unobserved is not guessed at",
              c.outcome == "tag_unobservable", f"{name}: {c.outcome} {c.seq_id}")
        check(f"{name} names the pair it could not separate",
              set(c.tied) == {"Sp_b_h1", "Sp_b_h1_ctrl"}, str(c.tied))

    # the same pair, long enough to reach the tag, resolves
    r1 = control[:250]
    r2 = revcomp(control[-250:])
    c = A.call_sequence(r1, panel, TAG, 0.10, 2, 3, second=r2)
    check("mates long enough to reach the tag do resolve it",
          c.outcome == "assigned" and c.seq_id == "Sp_b_h1_ctrl", str(c))

    # and the tag is searched on BOTH mates, not only the first
    tail_tag = by_id["Sp_a_h1_ctrl"].sequence   # tag 30bp from the end
    r1 = tail_tag[:150]
    r2 = revcomp(tail_tag[-150:])
    c = A.call_sequence(r1, panel, TAG, 0.10, 2, 3, second=r2)
    check("a tag visible only to the reverse mate is still found", c.tag,
          str(c))


def test_preflight_catches_it() -> None:
    print("\npre-flight refuses a panel with an unreadable tag")
    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(31)
        fwd = "ACGTACGTACGTACGTAAAA"
        rev = "GGCCGGCCGGCCGGCCGGCC"
        core = "".join(rng.choice("ACGT") for _ in range(640))
        amp = fwd + core + revcomp(rev)
        mid = len(amp) // 2
        with open(os.path.join(tmp, "primers.tsv"), "w") as fh:
            fh.write(f"name\tdirection\tsequence\nF1\tfwd\t{fwd}\nR1\trev\t{rev}\n")
        with open(os.path.join(tmp, "reference.tsv"), "w") as fh:
            fh.write("seq_id\tlocus\ttype\tpurpose\ttaxon\tparent_id\t"
                     "default_conc\tconc_units\tcorr_factor\tsequence\n")
            fh.write(f"A_h1\tITS2\ttarget\tspeciation\tSpecies_a\t\t\t\t\t{amp}\n")
            fh.write(f"A_h1_ctrl\tITS2\tcontrol\t\tSpecies_a\tA_h1\t\t\t\t"
                     f"{amp[:mid] + TAG + amp[mid:]}\n")

        def run(read_len):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rc = pf.main([
                    "--reference", os.path.join(tmp, "reference.tsv"),
                    "--primers", os.path.join(tmp, "primers.tsv"),
                    "--assay", ASSAY,
                    "--illumina-read-length", str(read_len),
                    "--outdir", os.path.join(tmp, f"pf{read_len}"),
                ])
            return rc, out.getvalue()

        rc, log = run(250)
        check("a mid-amplicon tag on an unmergeable amplicon is an error",
              rc != 0 and "where neither read can see it" in log,
              log.strip().splitlines()[-1] if log.strip() else "")
        rc, log = run(400)
        check("the same panel passes when the reads are long enough to merge",
              rc == 0, log.strip().splitlines()[-1] if log.strip() else "")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_tag_under_a_primer() -> None:
    """A tag that primer trimming would remove.

    Inserting a long tag inside a primer site breaks the primer match, and
    pre-flight already reports that. The case that slips through is a SHORT
    tag inside a LONG primer: the insertion stays within the primer's edit
    tolerance, so the primer is still located, and the tag sits under it -
    replaced by primer sequence in every read, silently.
    """
    print("\na tag that primer trimming would remove")
    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(41)
        # 32bp primer: tolerance is 32//8 = 4 edits, so a 3bp insertion
        # inside it does not stop the primer being found
        fwd = "".join(rng.choice("ACGT") for _ in range(32))
        rev = "".join(rng.choice("ACGT") for _ in range(20))
        core = "".join(rng.choice("ACGT") for _ in range(300))
        amp = fwd + core + revcomp(rev)
        short_tag = "CGC"
        at = 16                       # halfway into the forward primer
        ctrl = amp[:at] + short_tag + amp[at:]

        with open(os.path.join(tmp, "primers.tsv"), "w") as fh:
            fh.write(f"name\tdirection\tsequence\nF1\tfwd\t{fwd}\n"
                     f"R1\trev\t{rev}\n")
        with open(os.path.join(tmp, "reference.tsv"), "w") as fh:
            fh.write("seq_id\tlocus\ttype\tpurpose\ttaxon\tparent_id\t"
                     "default_conc\tconc_units\tcorr_factor\tsequence\n")
            fh.write(f"A_h1\tITS2\ttarget\tspeciation\tSpecies_a\t\t\t\t\t{amp}\n")
            fh.write(f"A_h1_ctrl\tITS2\tcontrol\t\tSpecies_a\tA_h1\t\t\t\t{ctrl}\n")

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = pf.main([
                "--reference", os.path.join(tmp, "reference.tsv"),
                "--primers", os.path.join(tmp, "primers.tsv"),
                "--assay", ASSAY, "--expected-tag", short_tag,
                "--outdir", os.path.join(tmp, "pf"),
            ])
        log = out.getvalue()
        check("a tag under a primer binding site is an error",
              rc != 0 and "under a primer binding site" in log,
              log.strip().splitlines()[-1] if log.strip() else "no output")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


def build(tmp: str, platform: str, n_samples: int, reads: int):
    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, n_samples, 2, "RUN01", platform)
    mf.write_reference(ref, os.path.join(tmp, "reference.tsv"))
    mf.write_primers(os.path.join(tmp, "primers.tsv"))
    mf.write_samplesheet(samples, os.path.join(tmp, "samplesheet.csv"))
    mf.write_runs(os.path.join(tmp, "runs.tsv"), "RUN01", platform,
                  os.path.join(tmp, "sim", "RUN01"))
    with contextlib.redirect_stdout(io.StringIO()):
        # the simulator first: pre-flight checks that the run directory
        # named in the manifest actually exists
        sim.main([
            "--reference", os.path.join(tmp, "reference.tsv"),
            "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
            "--runs", os.path.join(tmp, "runs.tsv"),
            "--assay", ASSAY, "--outdir", os.path.join(tmp, "sim"),
            "--reads-per-sample", str(reads),
        ])
        rc = pf.main([
            "--reference", os.path.join(tmp, "reference.tsv"),
            "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
            "--runs", os.path.join(tmp, "runs.tsv"),
            "--primers", os.path.join(tmp, "primers.tsv"),
            "--index-map", INDEX_MAP, "--assay", ASSAY,
            "--outdir", os.path.join(tmp, "preflight"),
        ])
    with open(os.path.join(tmp, "sim", "RUN01.truth.tsv"), encoding="utf-8") as fh:
        truth = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    panel = A.Panel(A.load_panel(
        os.path.join(tmp, "preflight", "reference.insert.fasta"),
        os.path.join(tmp, "reference.tsv")))
    return rc, truth, panel


def test_illumina_end_to_end() -> None:
    print("\nIllumina, read to call, against the truth table")
    tmp = tempfile.mkdtemp()
    try:
        rc, truth, panel = build(tmp, "illumina", 10, 300)
        check("the fixture panel passes pre-flight", rc == 0, f"rc={rc}")
        meta = {r.seq_id: r for r in panel.refs}
        ps = P.load_primers(os.path.join(tmp, "primers.tsv"))
        cfg = P.load_assay(ASSAY)
        rt1 = revcomp(cfg["pcr1_overhang_rev"])
        rt2 = revcomp(cfg["pcr1_overhang_fwd"])

        res, kinds, tagq = Counter(), Counter(), Counter()
        for r1p in sorted(glob.glob(os.path.join(tmp, "sim", "RUN01",
                                                 "*_R1_001.fastq.gz"))):
            pairs = zip(P.open_reads(r1p),
                        P.open_reads(r1p.replace("_R1_", "_R2_")))
            for (rid, s1, q1), (_r, s2, q2) in pairs:
                t = truth[rid]
                if t["fate"] not in ("clean", "carryover"):
                    continue
                f, r = ps.match(s1, "fwd"), ps.match(s2, "rev")
                if f is None or r is None:
                    continue
                a1, b1 = P.trim_readthrough(s1[f[1]:], q1[f[1]:], rt1)
                a2, b2 = P.trim_readthrough(s2[r[1]:], q2[r[1]:], rt2)
                m = P.merge_pair(a1, b1, a2, b2)
                c = (A.call_sequence(m[0], panel, TAG, 0.10, 2, 3) if m
                     else A.call_sequence(a1, panel, TAG, 0.10, 2, 3, second=a2))
                res["total"] += 1
                truth_type = meta[t["seq_id"]].type
                tagq[(c.tag, truth_type == "control")] += 1
                if c.outcome != "assigned":
                    res[c.outcome] += 1
                elif c.seq_id == t["seq_id"]:
                    res["correct"] += 1
                else:
                    res["wrong"] += 1
                    got, want = meta[c.seq_id], meta[t["seq_id"]]
                    kinds["haplotype" if got.taxon == want.taxon
                          and got.type == want.type
                          else "target/control" if got.taxon == want.taxon
                          else "species"] += 1

        acc = res["correct"] / res["total"]
        check("nearly every clean read reaches its own reference",
              acc > 0.99, f"{acc:.4f} of {res['total']}")
        check("no read is assigned to the wrong species",
              kinds["species"] == 0, f"{kinds['species']}")
        check("no control read is reported as its target",
              kinds["target/control"] == 0, f"{kinds['target/control']}")
        check("the tag never fires on a non-control read",
              tagq[(True, False)] == 0, str(tagq[(True, False)]))
        check("the tag fires on every control read",
              tagq[(False, True)] == 0, str(tagq[(False, True)]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_nanopore_end_to_end() -> None:
    print("\nnanopore, read to call, against the truth table")
    dx = _load("demux_ont", "bin/demux_ont.py")
    tmp = tempfile.mkdtemp()
    try:
        _rc, truth, panel = build(tmp, "nanopore", 8, 300)
        meta = {r.seq_id: r for r in panel.refs}
        with contextlib.redirect_stdout(io.StringIO()):
            dx.main([
                "--reads", os.path.join(tmp, "sim", "RUN01", "RUN01.fastq.gz"),
                "--samplesheet", os.path.join(tmp, "samplesheet.csv"),
                "--run-id", "RUN01", "--assay", ASSAY,
                "--outdir", os.path.join(tmp, "demux"),
            ])
        ps = P.load_primers(os.path.join(tmp, "primers.tsv"))
        res, kinds = Counter(), Counter()
        for fq in sorted(glob.glob(os.path.join(tmp, "demux", "reads",
                                                "*.fastq.gz"))):
            for rid, seq, _q in P.open_reads(fq):
                t = truth[rid]
                if t["fate"] not in ("clean", "carryover"):
                    continue
                f = ps.match(seq, "fwd")
                r = ps.match(revcomp(seq), "rev")
                if f is None or r is None:
                    res["no_primer"] += 1
                    continue
                c = A.call_sequence(seq[f[1]: len(seq) - r[1]], panel, TAG,
                                    0.10, 2, 3)
                res["total"] += 1
                if c.outcome != "assigned":
                    res[c.outcome] += 1
                elif c.seq_id == t["seq_id"]:
                    res["correct"] += 1
                else:
                    res["wrong"] += 1
                    got, want = meta[c.seq_id], meta[t["seq_id"]]
                    kinds["haplotype" if got.taxon == want.taxon
                          and got.type == want.type
                          else "target/control" if got.taxon == want.taxon
                          else "species"] += 1

        acc = res["correct"] / max(1, res["total"])
        check("assignment survives a 1% error rate",
              acc > 0.97, f"{acc:.4f} of {res['total']}")
        check("the extra error costs haplotypes, not species",
              kinds["species"] == 0 and kinds["target/control"] == 0,
              f"species {kinds['species']}, target/control "
              f"{kinds['target/control']}")
        # what it does cost should show up as declared ambiguity
        check("what cannot be resolved is declared, not guessed",
              res["wrong"] == 0, f"{res['wrong']} wrong calls")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_outputs() -> None:
    print("\noutput tables")
    tmp = tempfile.mkdtemp()
    try:
        _rc, _truth, _panel = build(tmp, "illumina", 6, 300)
        prep = os.path.join(tmp, "prep")
        r1s = sorted(glob.glob(os.path.join(tmp, "sim", "RUN01",
                                            "*_R1_001.fastq.gz")))
        # A field sample. A blank has no reads to compose, and a mock is
        # built from control plasmids, so neither has a species composition
        # - which is itself asserted below.
        field = [f for f in r1s if os.path.basename(f).startswith(
            ("iCCH", "POS"))]
        mock = [f for f in r1s if os.path.basename(f).startswith("MOCK")]
        def run_one(path):
            one = os.path.basename(path).split("_S")[0]
            with contextlib.redirect_stdout(io.StringIO()):
                P.main(["--r1", path, "--r2", path.replace("_R1_", "_R2_"),
                        "--sample-id", one,
                        "--primers", os.path.join(tmp, "primers.tsv"),
                        "--assay", ASSAY, "--outdir", prep, "--quiet"])
                A.main(["--panel", os.path.join(tmp, "preflight",
                                                "reference.insert.fasta"),
                        "--reference", os.path.join(tmp, "reference.tsv"),
                        "--sample-id", one,
                        "--merged", os.path.join(prep, f"{one}.merged.fasta.gz"),
                        "--pairs", os.path.join(prep, f"{one}.pairs.tsv.gz"),
                        "--assay", ASSAY,
                        "--outdir", os.path.join(tmp, "assign"), "--quiet"])
            return one

        # Prevalence is low by design, so a given well may hold no target at
        # all - a real outcome, but not one that can demonstrate a
        # composition. Take the first well that has something in it.
        sid = None
        for path in field:
            one = run_one(path)
            js = json.load(open(os.path.join(tmp, "assign",
                                             f"{one}.assign.json")))
            if js["composition_denominator"] > 0:
                sid = one
                break
        check("at least one well carries target reads", sid is not None)
        if sid is None:
            return

        s = json.load(open(os.path.join(tmp, "assign", f"{sid}.assign.json")))
        check("assignment accounting balances", s["accounting_balances"],
              str(s["outcomes"]))
        check("tag and mapping do not contradict each other",
              s["tag_mapping_conflicts"] == 0,
              str(s.get("tag_vs_mapping")))

        taxa = list(csv.DictReader(
            open(os.path.join(tmp, "assign", f"{sid}.taxa.tsv")),
            delimiter="\t"))
        total = sum(float(r["rel_abundance"]) for r in taxa if r["rel_abundance"])
        check("relative abundance is a composition and sums to one",
              abs(total - 1.0) < 1e-6, f"{total:.6f}")

        counts = list(csv.DictReader(
            open(os.path.join(tmp, "assign", f"{sid}.counts.tsv")),
            delimiter="\t"))
        types = {r["type"] for r in counts if int(r["reads"]) > 0}
        check("quant and control reads are counted but kept out of the "
              "composition",
              "quant" in types
              and not any(r["taxon"].startswith("QUANT") for r in taxa),
              str(sorted(types)))
        check("every reference appears in the counts table, zeros included",
              len(counts) == len(_panel.refs), f"{len(counts)}")

        # A mock population is made of control plasmids. It is a positive
        # control for the chemistry, not an observation of species, and it
        # must not turn up in a composition table as though it were one.
        msid = run_one(mock[0])
        ms = json.load(open(os.path.join(tmp, "assign", f"{msid}.assign.json")))
        mtaxa = list(csv.DictReader(
            open(os.path.join(tmp, "assign", f"{msid}.taxa.tsv")),
            delimiter="\t"))
        by_type = ms["reads_by_reference_type"]
        check("a mock built from plasmids yields control reads",
              by_type.get("control", 0) > 0, str(by_type))
        # Not zero: index hopping from real wells puts a few target reads
        # into a mock, and pretending otherwise would hide the very effect
        # this pipeline measures. What must hold is that the mock is
        # overwhelmingly control-derived and that the trace of target is
        # left below the detection limit rather than called.
        ctrl = by_type.get("control", 0)
        denom = ms["composition_denominator"]
        check("and its target signal is trace-level beside the controls",
              denom == 0 or ctrl > 20 * denom,
              f"{denom} target vs {ctrl} control reads")
        check("the mock composition table names at most a trace taxon",
              len(mtaxa) <= 1, f"{len(mtaxa)} taxa")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_bounded_edit()
    test_calls()
    test_tag_unobservable()
    test_preflight_catches_it()
    test_tag_under_a_primer()
    test_outputs()
    test_illumina_end_to_end()
    test_nanopore_end_to_end()

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
