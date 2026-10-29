#!/usr/bin/env python3
"""
Tests for the VESPA benchmark's scoring (tests/bench_vespa_public.py).

The NHP pairing with the paper's samples is inferred from our results, so
if it were inferred per window set a sample could change partner between
window sets and the "difference" between them be the re-pairing. Here two
window sets are built that would pair two samples opposite ways; both
must be scored against the same partners.

Run: python3 tests/test_bench_vespa_public.py
"""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

spec = importlib.util.spec_from_file_location(
    "bench_vespa_public", os.path.join(HERE, "bench_vespa_public.py"))
bench = importlib.util.module_from_spec(spec)
sys.modules["bench_vespa_public"] = bench    # for the process pool
spec.loader.exec_module(bench)

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


PROTIST = ("Eukaryota;Amoebozoa;Evosea;Evosea_X;Archamoebea;Archamoebea_X;"
           "Entamoebidae;Entamoeba")
FUNGUS = ("Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Saccharomycotina;"
          "Saccharomycetales;Candida")
HEAD = ["sample_id", "cluster", "reads", "rank", "taxon", "lineage"]


def clusters(path: str, sid: str, protists: int, fungi: int) -> None:
    with open(path, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEAD)
        w.writerow([sid, f"{sid}_1", protists, "genus", "Entamoeba", PROTIST])
        w.writerow([sid, f"{sid}_2", fungi, "genus", "Candida", FUNGUS])


def partners(path: str) -> dict:
    with open(path) as fh:
        return {r["sample"]: r["paper_sample"]
                for r in csv.DictReader(fh, delimiter="\t")}


def test_pairing_once() -> None:
    print("\nthe NHP pairing is decided once, for every window set")
    vespa = {"groups": {"NHP1": {"Protozoa": 1.0}, "NHP2": {"Fungi": 1.0}},
             "helminths": {}, "protozoa": {}}
    with tempfile.TemporaryDirectory() as tmp:
        rep = os.path.join(tmp, "filereport.tsv")
        with open(rep, "w") as fh:
            fh.write("sample_alias\trun_accession\nNHPa\tSRR1\nNHPb\tSRR2\n")
        # 1x1: NHPa is protozoa, NHPb fungi; 2x2 the other way round
        for tag, a, b in (("w1x1", (90, 10), (10, 90)),
                          ("w2x2", (10, 90), (90, 10))):
            os.makedirs(os.path.join(tmp, tag))
            clusters(os.path.join(tmp, tag, "NHPa.clusters.tsv"), "NHPa", *a)
            clusters(os.path.join(tmp, tag, "NHPb.clusters.tsv"), "NHPb", *b)
        each = {t: {s: v[0] for s, v in bench.remap(
            ["NHPa", "NHPb"], vespa, os.path.join(tmp, t), 0.001).items()}
            for t in ("w1x1", "w2x2")}
        check("the fixture would pair differently per window set",
              each["w1x1"] != each["w2x2"], str(each))
        old = bench.vespa_tables
        bench.vespa_tables = lambda path: vespa
        try:
            def run(*extra, out=tmp):
                with contextlib.redirect_stdout(io.StringIO()):
                    bench.main(["--reads", tmp, "--filereport", rep,
                                "--source-data", "unused.xlsx",
                                "--release", tmp, "--index", "unused",
                                "--windows", "1x1,2x2", "--outdir", out,
                                "--jobs", "1", *extra])
                return [partners(os.path.join(out, f"{t}.comparison.tsv"))
                        for t in ("w1x1", "w2x2")]
            p1, p2 = run()
            check("both window sets are scored against the same partners",
                  p1 == p2 == each["w1x1"], f"{p1} {p2}")
            check("one nhp_mapping.tsv, none per window set",
                  os.path.isfile(os.path.join(tmp, "nhp_mapping.tsv"))
                  and not any(f.endswith(".nhp_mapping.tsv")
                              for f in os.listdir(tmp)))
            p1, p2 = run("--pair-on", "2x2")
            check("--pair-on picks the window set the pairing is decided on",
                  p1 == p2 == each["w2x2"], f"{p1} {p2}")
            fixed = os.path.join(tmp, "fixed.tsv")
            bench.write_mapping(fixed, {"NHPa": ("NHP2", 0.0),
                                        "NHPb": ("NHP1", 0.0)})
            p1, p2 = run("--nhp-mapping", fixed)
            check("--nhp-mapping is used as given, for every window set",
                  p1 == p2 == {"NHPa": "NHP2", "NHPb": "NHP1"}, f"{p1} {p2}")
        finally:
            bench.vespa_tables = old


def main() -> int:
    test_pairing_once()
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
