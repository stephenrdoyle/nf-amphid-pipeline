#!/usr/bin/env python3
"""
Tests for database-mode synthesis and report (bin/synthesise_db.py,
bin/report_db.py).

Built on hand-written assignment outputs rather than a pipeline run, so
each rule is tested where it can be seen to fire: a taxon in a blank
must be beaten, not merely matched; a call under the relative floor is
not a detection; a shallow well calls nothing; the off-target groups
come from the PR2 lineage; no blanks is said, not hidden; a sample that
had reads and no assignment stops the run.

Run: python3 tests/test_synthesis_db.py
"""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


syn = _load("synthesise_db", "bin/synthesise_db.py")
rep = _load("report_db", "bin/report_db.py")

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


LIN = {
    "hc": "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Haemonchus;Haemonchus_contortus",
    "sc": "Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Saccharomycotina;Saccharomycetales;Saccharomyces;Saccharomyces_cerevisiae",
    "bt": "Eukaryota;Obazoa;Opisthokonta;Metazoa;Craniata;Craniata_X;Mammalia;Bos;Bos_taurus",
    "ec": "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacterales_X;Enterobacteriaceae;Enterobacteriaceae_X;Escherichia;Escherichia_coli",
    "pl": "Eukaryota:plas;Archaeplastida;Chlorophyta;Chlorophyta_X;Trebouxiophyceae;Chlorellales;Chlorellaceae;Chlorella;Chlorella_sp.",
    "tg": "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Teladorsagia;Teladorsagia_circumcincta",
}
HITS = {"hc": "V4_hc", "sc": "V4_sc", "bt": "V4_bt", "ec": "V4_ec",
        "pl": "V4_pl", "tg": "V4_tg"}


def cluster(sid, k, key, reads, rank="species", ident="1.0000", label=None):
    lin = LIN[key].split(";")
    d = {"species": 8, "genus": 7, "family": 6, "order": 5}[rank]
    taxon = lin[d]
    return {"sample_id": sid, "cluster": f"{sid}_c{k}", "reads": reads,
            "rel_abundance": "", "rank": rank, "taxon": taxon,
            "label": label or taxon, "identity": ident, "coverage": "1.000",
            "lineage": ";".join(lin[:d + 1]), "best_hit": HITS[key],
            "best_hit_lineage": LIN[key], "near_hits": 1,
            "lca_rank": rank, "db_resolution": "species", "origin": "group",
            "merged": 0, "length": 500, "database_version": "t-1",
            "consensus": "ACGT"}


def write_case(tmp, wells, types, demux=None, primer_fit=None):
    ad = os.path.join(tmp, "assign_db")
    pf = os.path.join(tmp, "preflight")
    rel = os.path.join(tmp, "release")
    for d in (ad, pf, rel):
        os.makedirs(d, exist_ok=True)
    for sid, rows in wells.items():
        clustered = sum(r["reads"] for r in rows)
        with open(os.path.join(ad, f"{sid}.clusters.tsv"), "w") as fh:
            w = csv.DictWriter(fh, list(cluster("x", 1, "hc", 1)),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        with open(os.path.join(ad, f"{sid}.assign_db.json"), "w") as fh:
            json.dump({"sample_id": sid, "clusters": len(rows),
                       "haplotype_splits": 0,
                       "reads": {"reads_in": clustered + 5, "no_primer": 5,
                                 "too_short": 0, "unclustered": 0,
                                 "clustered": clustered},
                       "settings": {"floors": {"species": 0.995}}}, fh)
    with open(os.path.join(pf, "file_resolution.tsv"), "w") as fh:
        fh.write("run_id\tsample_id\tsanitised_id\tsample_type\n")
        for sid, t in types.items():
            fh.write(f"R1\t{sid}\t{sid}\t{t}\n")
    with open(os.path.join(pf, "preflight_summary.json"), "w") as fh:
        json.dump({"reference_database": {"minimap2": "2.31", "lookalikes": 3},
                   "assay": {"name": "vespa18s_v4"}}, fh)
    with open(os.path.join(rel, "manifest.json"), "w") as fh:
        json.dump({"name": "test_release", "version": "t-1",
                   "primers": {"forward": {"sequence": "AGC"},
                               "reverse": {"sequence": "TCC"}}}, fh)
    if primer_fit:
        with open(os.path.join(rel, "primer_fit.tsv"), "w") as fh:
            fh.write("id\tfit\n")
            for k, v in primer_fit.items():
                fh.write(f"{k}\t{v}\n")
    args = ["--assign-dir", ad, "--preflight-dir", pf, "--database", rel,
            "--run-id", "R1", "--outdir", os.path.join(tmp, "synthesis"),
            "--min-depth", "50"]
    if demux is not None:
        with open(os.path.join(tmp, "demux.json"), "w") as fh:
            json.dump({"per_sample": demux, "by_tier": {"E": 4},
                       "total_reads": 999}, fh)
        args += ["--demux-summary", os.path.join(tmp, "demux.json")]
    return args


def run(args):
    with contextlib.redirect_stdout(io.StringIO()):
        syn.main(args)
    out = args[args.index("--outdir") + 1]
    with open(os.path.join(out, "db_calls.tsv")) as fh:
        calls = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(out, "db_wells.tsv")) as fh:
        wells = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(os.path.join(out, "db_run_qc.json")) as fh:
        qc = json.load(fh)
    return calls, wells, qc


def call_of(calls, sid, taxon):
    for c in calls:
        if c["sample_id"] == sid and c["taxon"] == taxon:
            return c
    return None


def test_calls() -> None:
    print("\ncalls against blanks and floors")
    tmp = tempfile.mkdtemp()
    try:
        wells = {
            "BLANK": [cluster("BLANK", 1, "sc", 8)],
            "S1": [cluster("S1", 1, "hc", 800), cluster("S1", 2, "sc", 8),
                   cluster("S1", 3, "bt", 100), cluster("S1", 4, "ec", 60),
                   cluster("S1", 5, "pl", 30)],
            "S2": [cluster("S2", 1, "hc", 9000), cluster("S2", 2, "sc", 30),
                   cluster("S2", 3, "tg", 5, rank="family", ident="0.9600")],
            "S3": [cluster("S3", 1, "hc", 20)],
        }
        types = {"BLANK": "neg_extraction", "S1": "unknown",
                 "S2": "unknown", "S3": "unknown"}
        calls, w, qc = run(write_case(tmp, wells, types,
                                      primer_fit={"V4_hc": "near_3prime"}))
        c = call_of(calls, "S1", "Saccharomyces_cerevisiae")
        check("a taxon at the blank's own count is not detected",
              c["call"] == "below_blanks" and c["limit_reads"] == "9", str(c))
        c = call_of(calls, "S2", "Saccharomyces_cerevisiae")
        check("and is detected once it beats the blank",
              c["call"] == "detected", str(c))
        c = call_of(calls, "S2", "Trichostrongylidae")
        check("a call under the relative floor is not a detection",
              c["call"] == "below_floor" and c["rank"] == "family", str(c))
        c = call_of(calls, "S3", "Haemonchus_contortus")
        check("a shallow well calls nothing", c["call"] == "low_depth"
              and "S3" in qc["low_depth_wells"], str(c))
        check("every row carries the release version",
              all(c["database_version"] == "t-1" for c in calls))
        s1 = w["S1"]
        tot = 800 + 8 + 100 + 60 + 30
        check("off-target fractions from the PR2 lineage",
              abs(float(s1["host_fraction"]) - 100 / tot) < 1e-3
              and abs(float(s1["bacteria_archaea_fraction"]) - 60 / tot) < 1e-3
              and abs(float(s1["organelle_fraction"]) - 30 / tot) < 1e-3
              and abs(float(s1["fungi_fraction"]) - 8 / tot) < 1e-3, str(s1))
        check("reads by rank reached",
              w["S2"]["rank_family"] == "5" and w["S2"]["rank_species"] == "9030",
              str(w["S2"]))
        check("a near_3prime best hit is flagged",
              call_of(calls, "S1", "Haemonchus_contortus")["primer_fit"]
              == "near_3prime" and any("Haemonchus" in x
                                       for x in qc["near_3prime_calls"]))
        check("taxa seen in the blanks are listed",
              qc["negatives"]["taxa_in_blanks"] ==
              ["species Saccharomyces_cerevisiae"])

        # the report
        out = os.path.join(tmp, "r.html")
        with contextlib.redirect_stdout(io.StringIO()):
            rep.main(["--synthesis-dir", os.path.join(tmp, "synthesis"),
                      "--run-id", "R1", "--out", out])
        doc = open(out).read()
        check("the report says relative 18S abundance, never copies",
              "relative 18S abundance" in doc and "never organisms or copies"
              in doc)
        check("the report names the release and the blank contamination",
              "test_release" in doc and "blanks not empty" in doc)
        check("the report is self-contained",
              "<script src" not in doc and "http" not in doc.split("<body>")[1])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_no_blanks_and_cohort() -> None:
    print("\nno blanks, and an incomplete cohort")
    tmp = tempfile.mkdtemp()
    try:
        wells = {"S1": [cluster("S1", 1, "hc", 500)]}
        args = write_case(tmp, wells, {"S1": "unknown", "S2": "unknown"},
                          demux={"S1": 505, "S2": 0})
        calls, w, qc = run(args)
        c = call_of(calls, "S1", "Haemonchus_contortus")
        check("with no blanks the limit says it is an assumption",
              "assumption" in c["limit_basis"] and qc["negatives"]["n"] == 0
              and qc["negatives"]["note"], str(c))
        check("a well with no reads is reported empty, not lost",
              qc["empty_wells"] == ["S2"] and "S2" in w, str(qc["empty_wells"]))
        check("internal adapter rates are carried", qc["internal_anchors"][
            "quarantined_internal_adapter"] == 4)
        out = os.path.join(tmp, "r.html")
        with contextlib.redirect_stdout(io.StringIO()):
            rep.main(["--synthesis-dir", os.path.join(tmp, "synthesis"),
                      "--out", out])
        check("and the report says there were no blanks",
              "no blanks" in open(out).read())
        args = write_case(tmp, wells, {"S1": "unknown", "S2": "unknown"},
                          demux={"S1": 505, "S2": 300})
        try:
            run(args)
            check("a sample with reads and no assignment stops synthesis",
                  False)
        except SystemExit as e:
            check("a sample with reads and no assignment stops synthesis",
                  "cohort incomplete" in str(e), str(e))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_calls()
    test_no_blanks_and_cohort()
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
