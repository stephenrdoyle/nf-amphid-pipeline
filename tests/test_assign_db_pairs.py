#!/usr/bin/env python3
"""
Tests for database-mode assignment on Illumina (bin/assign_db_pairs.py).

The pieces first - orienting and trimming a pair, the read-through that
makes a short insert one sequence, the denoiser that folds error copies
onto their source and keeps real variants - then the fixture's Illumina
run end to end, where every present taxon must be named and nothing
held out of the release may borrow a relative's species name.

The end-to-end test needs minimap2 (on the PATH or in ~/bin).

Run: python3 tests/test_assign_db_pairs.py
"""

from __future__ import annotations

import contextlib
import csv
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
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


adp = _load("assign_db_pairs", "bin/assign_db_pairs.py")
import ampseqlib as al          # noqa: E402
import make_db_fixture as fx    # noqa: E402

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def minimap2() -> str | None:
    for p in (shutil.which("minimap2"), os.path.expanduser("~/bin/minimap2")):
        if p and os.path.isfile(p):
            return p
    return None


def rand(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_orient() -> None:
    print("\norienting and trimming a pair")
    ps = al.PrimerSet(fwd=[fx.FWD], rev=[fx.REV])
    rng = random.Random(1)
    ins = rand(rng, 400)
    amp = fx.amplicon(ins)
    r1 = ("a", amp[:150], "I" * 150)
    r2 = ("a", al.revcomp(amp)[:150], "I" * 150)
    o = adp.orient(r1, r2, ps)
    check("R1 from the forward primer: both primers trimmed",
          o is not None and o[0][0] == ins[:132]
          and o[1][0] == al.revcomp(ins)[:128], str(o)[:80])
    o2 = adp.orient(r2, r1, ps)
    check("mates swapped (R1 from the reverse primer): put back",
          o2 is not None and o2[0][0] == o[0][0] and o2[1][0] == o[1][0])
    check("no primer pair: None",
          adp.orient(("a", rand(rng, 150), "I" * 150), r2, ps) is None)


def test_denoiser() -> None:
    print("\ndenoising")
    rng = random.Random(3)
    a = rand(rng, 260)
    b = a[:100] + ("A" if a[100] != "A" else "C") + a[101:]       # 1 off a
    c = fx.mf.mutate(rng, a, 6)                                   # a variant

    def err(s, k):
        s = list(s)
        for j in rng.sample(range(len(s)), k):
            s[j] = rng.choice([x for x in "ACGT" if x != s[j]])
        return "".join(s)

    counts = Counter({a: 400, b: 30, c: 60})
    for _ in range(40):
        counts[err(a, 1)] += 1
    for _ in range(10):
        counts[err(c, 2)] += 1
    counts[rand(rng, 260)] += 1                                   # noise
    dn = adp.Denoiser(max_diffs=3, alpha=2.0, min_size=2)
    member, lost = dn.run(counts)
    cents = set(dn.centroids)
    check("the source and a variant six bases away both survive",
          a in cents and c in cents, str(len(cents)))
    check("one base off at 30/400 - under an eighth - is an error copy",
          b not in cents and member[b] == dn.centroids.index(a))
    dn2 = adp.Denoiser(3, 2.0, 2)
    dn2.run(Counter({a: 400, b: 60}))
    check("one base off at 60/400 - over an eighth - is a sequence of its own",
          b in dn2.centroids)
    b2 = fx.mf.mutate(rng, a, 2)
    dn3 = adp.Denoiser(3, 2.0, 2)
    dn3.run(Counter({a: 400, b2: 20}))
    check("two bases off at 20/400 - over a 32nd - is kept",
          b2 in dn3.centroids)
    check("singleton error copies join their source",
          all(member[s] == dn.centroids.index(a) for s in counts
              if s != a and s != b and counts[s] == 1
              and sum(x != y for x, y in zip(s, a)) == 1))
    check("a singleton near nothing is unclustered", lost == 1, str(lost))


def test_whole_insert() -> None:
    print("\nshort inserts are kept whole")
    mm2 = minimap2()
    rng = random.Random(5)
    with tempfile.TemporaryDirectory() as tmp:
        species = fx.build_tree(random.Random(3))
        rel = os.path.join(tmp, "release")
        fx.write_release(rel, species)
        fx.write_primers(os.path.join(tmp, "primers.tsv"))
        # a 100bp insert: both mates read through into the adapter
        short = rand(rng, 100)
        amp = fx.amplicon(short)
        tail1 = al.revcomp(fx.al.load_assay(fx.ASSAY)["pcr1_overhang_rev"]) + "A" * 150
        tail2 = al.revcomp(fx.al.load_assay(fx.ASSAY)["pcr1_overhang_fwd"]) + "A" * 150
        with gzip.open(os.path.join(tmp, "r1.fq.gz"), "wt") as f1, \
                gzip.open(os.path.join(tmp, "r2.fq.gz"), "wt") as f2:
            for i in range(20):
                f1.write(f"@{i} 1\n{(amp + tail1)[:150]}\n+\n{'I' * 150}\n")
                f2.write(f"@{i} 2\n{(al.revcomp(amp) + tail2)[:150]}\n+\n{'I' * 150}\n")
        if not mm2:
            print("  SKIP: minimap2 not found")
            return
        with contextlib.redirect_stdout(io.StringIO()):
            adp.main(["--r1", os.path.join(tmp, "r1.fq.gz"),
                      "--r2", os.path.join(tmp, "r2.fq.gz"),
                      "--database", rel, "--primers",
                      os.path.join(tmp, "primers.tsv"), "--assay", fx.ASSAY,
                      "--sample-id", "X", "--outdir", tmp, "--minimap2", mm2])
        js = json.load(open(os.path.join(tmp, "X.assign_db.json")))
        rows = list(csv.DictReader(open(os.path.join(tmp, "X.clusters.tsv")),
                                   delimiter="\t"))
        check("a pair whose insert ends inside R1 is one whole sequence",
              js["reads"]["whole_insert"] == 20 and rows
              and rows[0]["consensus"] == short
              and rows[0]["origin"] == "whole_insert", str(js["reads"]))


def test_divergent_window() -> None:
    print("\na window too divergent for -x sr is aligned end to end")
    mm2 = minimap2()
    if not mm2:
        print("  SKIP: minimap2 not found")
        return
    with tempfile.TemporaryDirectory() as tmp:
        species = fx.build_tree(random.Random(3))
        rel = os.path.join(tmp, "release")
        fx.write_release(rel, species)
        fx.write_primers(os.path.join(tmp, "primers.tsv"))
        src = next(sp for sp, d in species.items() if not d["held_out"])
        ref = species[src]["seq"]
        # VESPA's NHP gregarines: a conserved start, then a stretch at a
        # third mismatched, then conserved again - a lineage ~12% from
        # its nearest record over the pair. -x sr's alignment of the
        # first window stops at the stretch.
        swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
        q = "".join(swap[b] if 40 <= i < 120 and i % 3 == 0 else b
                    for i, b in enumerate(ref))
        amp = fx.amplicon(q)
        with gzip.open(os.path.join(tmp, "r1.fq.gz"), "wt") as f1, \
                gzip.open(os.path.join(tmp, "r2.fq.gz"), "wt") as f2:
            for i in range(20):
                f1.write(f"@{i} 1\n{amp[:150]}\n+\n{'I' * 150}\n")
                f2.write(f"@{i} 2\n{al.revcomp(amp)[:150]}\n+\n{'I' * 150}\n")

        def run(clip):
            old, adp.SR_CLIP = adp.SR_CLIP, clip
            try:
                out = os.path.join(tmp, f"out{clip}")
                with contextlib.redirect_stdout(io.StringIO()):
                    adp.main(["--r1", os.path.join(tmp, "r1.fq.gz"),
                              "--r2", os.path.join(tmp, "r2.fq.gz"),
                              "--database", rel, "--primers",
                              os.path.join(tmp, "primers.tsv"), "--assay",
                              fx.ASSAY, "--sample-id", "X", "--outdir", out,
                              "--minimap2", mm2])
            finally:
                adp.SR_CLIP = old
            rows = list(csv.DictReader(open(os.path.join(out, "X.clusters.tsv")),
                                       delimiter="\t"))
            js = json.load(open(os.path.join(out, "X.assign_db.json")))
            return rows[0], js
        sr, _ = run(10.0)           # never re-mapped: -x sr alone
        row, js = run(adp.SR_CLIP)
        # the expected identity: mismatches over the two windows
        w1, w2 = len(row["consensus"].split("N" * 10)[0]), 128
        diffs = sum(a != b for a, b in zip(q[:w1], ref[:w1]))
        want = 1 - diffs / (w1 + w2)
        check("-x sr alone covers the divergent window only in part",
              float(sr["coverage"]) < 0.9, f"{sr['coverage']} {sr['identity']}")
        check("re-mapped, both windows align end to end on the right record",
              float(row["coverage"]) >= 0.99 and row["origin"] == "pair"
              and row["best_hit"] == fx._id(ref)
              and js["remapped_divergent"] == 20,
              f"{row['coverage']} {row['origin']} {js['remapped_divergent']}")
        check("identity is the mismatches over both windows, not half of them",
              abs(float(row["identity"]) - want) < 0.01,
              f"{row['identity']} vs {want:.4f}")


def test_end_to_end() -> None:
    print("\nthe fixture's Illumina run")
    mm2 = minimap2()
    if not mm2:
        print("  SKIP: minimap2 not found - database mode cannot be tested")
        FAIL.append(("minimap2 present", "not on PATH or in ~/bin"))
        return
    with tempfile.TemporaryDirectory() as tmp:
        rng = random.Random(20261028)
        species = fx.build_tree(rng)
        rel = os.path.join(tmp, "release")
        man = fx.write_release(rel, species)
        fx.write_primers(os.path.join(tmp, "primers.tsv"))
        sheet = [{"sample_id": s} for s, _t, _c in fx.SAMPLES]
        cfg = dict(fx.sim.DEFAULTS)
        cfg.update(al.load_assay(fx.ASSAY))
        fx.write_illumina(os.path.join(tmp, "run"), sheet, species, cfg,
                          random.Random(9))
        missing, wrong, ok = [], [], True
        for sid, _t, comp in fx.SAMPLES:
            with contextlib.redirect_stdout(io.StringIO()):
                adp.main(["--r1", os.path.join(tmp, "run", f"{sid}_R1.fastq.gz"),
                          "--r2", os.path.join(tmp, "run", f"{sid}_R2.fastq.gz"),
                          "--database", rel, "--primers",
                          os.path.join(tmp, "primers.tsv"), "--assay", fx.ASSAY,
                          "--sample-id", sid, "--outdir",
                          os.path.join(tmp, "out"), "--minimap2", mm2])
            rows = list(csv.DictReader(open(os.path.join(
                tmp, "out", f"{sid}.clusters.tsv")), delimiter="\t"))
            js = json.load(open(os.path.join(tmp, "out", f"{sid}.assign_db.json")))
            r = js["reads"]
            if r["reads_in"] != (r["no_primer"] + r["too_short"]
                                 + r["low_quality"] + r["clustered"]
                                 + r["unclustered"]):
                ok = False
            if js["database"]["version"] != man["version"]:
                ok = False
            for src, n in comp.items():
                if n < 5 or src.endswith("#hap2"):
                    continue
                d = species[src]
                if not d["held_out"]:
                    if not any(x["taxon"] == src and x["rank"] == "species"
                               for x in rows):
                        missing.append(f"{sid}:{src}")
                    continue
                # held out: no species call in its genus beyond what the
                # sample really holds
                genus = d["lineage"][7]
                for x in rows:
                    lin = x["lineage"].split(";")
                    if (x["rank"] == "species" and len(lin) > 7
                            and lin[7] == genus and x["taxon"] not in comp):
                        wrong.append(f"{sid}: {src} named {x['taxon']}")
                    if (rank_of(x) == "genus" and src.startswith("Teladorsagia")
                            and x["taxon"] != "Teladorsagia" and len(lin) > 6
                            and lin[6] == d["lineage"][6] and
                            x["taxon"] not in {species[s]["lineage"][7]
                                               for s in comp if s in species}):
                        wrong.append(f"{sid}: {src} named genus {x['taxon']}")
        check("every present source of 5+ pairs is named to species",
              not missing, "; ".join(missing))
        check("held-out taxa never take a relative's name at 2x150",
              not wrong, "; ".join(wrong))
        check("pair accounting closes, and the release version is recorded", ok)


def test_gap_chimera() -> None:
    print("\na chimera joined in the unread gap")
    mm2 = minimap2()
    if not mm2:
        print("  SKIP: minimap2 not found")
        return
    rng = random.Random(21)
    with tempfile.TemporaryDirectory() as tmp:
        species = fx.build_tree(random.Random(3))
        rel = os.path.join(tmp, "release")
        fx.write_release(rel, species)
        fx.write_primers(os.path.join(tmp, "primers.tsv"))
        a = fx.amplicon(species["Haemonchus_contortus"]["seq"])
        b = fx.amplicon(species["Ascaris_suum"]["seq"])
        cfg = al.load_assay(fx.ASSAY)
        t1 = al.revcomp(cfg["pcr1_overhang_rev"]) + "A" * 150
        t2 = al.revcomp(cfg["pcr1_overhang_fwd"]) + "A" * 150
        with gzip.open(os.path.join(tmp, "r1.fq.gz"), "wt") as f1, \
                gzip.open(os.path.join(tmp, "r2.fq.gz"), "wt") as f2:
            k = 0
            for left, right, n in ((a, a, 300), (b, b, 200), (a, b, 25)):
                for _ in range(n):
                    f1.write(f"@{k} 1\n{(left + t1)[:150]}\n+\n{'I' * 150}\n")
                    f2.write(f"@{k} 2\n{(al.revcomp(right) + t2)[:150]}\n+\n"
                             f"{'I' * 150}\n")
                    k += 1
        with contextlib.redirect_stdout(io.StringIO()):
            adp.main(["--r1", os.path.join(tmp, "r1.fq.gz"),
                      "--r2", os.path.join(tmp, "r2.fq.gz"),
                      "--database", rel, "--primers",
                      os.path.join(tmp, "primers.tsv"), "--assay", fx.ASSAY,
                      "--sample-id", "X", "--outdir", tmp, "--minimap2", mm2])
        rows = list(csv.DictReader(open(os.path.join(tmp, "X.clusters.tsv")),
                                   delimiter="\t"))
        js = json.load(open(os.path.join(tmp, "X.assign_db.json")))
        ch = [r for r in rows if r["origin"] == "chimera"]
        check("R1 from one species and R2 from another: a chimera joined in "
              "the gap, not named", len(ch) == 1 and ch[0]["reads"] == "25"
              and ch[0]["chimera_of"].endswith("@gap") and not ch[0]["rank"]
              and js["chimeras"]["reads"] == 25,
              str([(r["reads"], r["origin"], r["chimera_of"], r["label"])
                   for r in rows]))
        check("its parents are still named (to genus or species: at 2x150 "
              "a sister can be identical over the windows)", sum(
            1 for r in rows if r["rank"] in ("species", "genus")) == 2,
            str([(r["reads"], r["rank"], r["label"]) for r in rows]))


def rank_of(x):
    return x["rank"]


def main() -> int:
    test_orient()
    test_denoiser()
    test_whole_insert()
    test_divergent_window()
    test_end_to_end()
    test_gap_chimera()
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
