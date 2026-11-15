#!/usr/bin/env python3
"""
Tests for database mode's layers and contracts (docs/CONTRACTS.md).

  * the protocol registry: every protocol loads, says what it does with
    every read type, and a release gets the protocol its manifest
    declares, or the one its manifest implies;
  * contract B: every protocol, on its own fixture, over every read type
    it supports, writes clusters.tsv, taxa.tsv, consensus.fasta and
    assign_db.json as the contract says - required columns in order, its
    declared extra columns after them, read fates that balance - and
    refuses the read types it does not support;
  * db_calls.tsv: synthesis over each protocol's output keeps the shared
    columns and adds a protocol's optional ones only where they exist;
  * contract A: the demux accounting folds a nanopore run's summaries in
    and says what an Illumina or predemultiplexed run cannot know, and
    synthesis reads a run through it exactly as through the summaries;
  * the predemultiplexed adapter: pre-flight resolves each sample's
    FASTQs, writes the read types, and refuses what cannot run.

The protocol runs need minimap2 (on the PATH or in ~/bin).

Run: python3 tests/test_contracts.py
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


adb = _load("assign_db", "bin/assign_db.py")
adp = _load("assign_db_pairs", "bin/assign_db_pairs.py")
sdb = _load("synthesise_db", "bin/synthesise_db.py")
acc = _load("demux_accounting", "bin/demux_accounting.py")
pf = _load("preflight", "bin/preflight.py")
nem = _load("test_nemits2", "tests/test_nemits2.py")
import ampseqlib as al            # noqa: E402
import protocols                  # noqa: E402
import make_db_fixture as fx      # noqa: E402
import make_fixtures as mf        # noqa: E402
import simulate as sim            # noqa: E402

PASS, FAIL = [], []
INDEX_MAP = os.path.join(ROOT, "assets", "index_map.tsv")


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


def quiet(fn, *a):
    with contextlib.redirect_stdout(io.StringIO()), \
            contextlib.redirect_stderr(io.StringIO()):
        return fn(*a)


# --------------------------------------------------------------------------
# fixtures: one per protocol, inputs for every read type it supports
# --------------------------------------------------------------------------


def vespa_fixture(tmp: str) -> dict:
    """The database fixture's release, three samples as nanopore reads
    and as 2x150 pairs."""
    species = fx.build_tree(random.Random(20261028))
    rel = os.path.join(tmp, "release")
    fx.write_release(rel, species)
    prim = os.path.join(tmp, "primers.tsv")
    fx.write_primers(prim)
    picked = [s for s in fx.SAMPLES if s[0] in ("BLANK-01", "S01", "S04")]
    rng = random.Random(7)
    model = sim.ERROR_MODELS["nanopore_hac"]
    single = []
    for sid, _t, comp in picked:
        path = os.path.join(tmp, f"{sid}.fastq.gz")
        with gzip.open(path, "wt") as fh:
            k = 0
            for src, n in comp.items():
                for _ in range(n):
                    m = fx.amplicon(species[src]["seq"])
                    if rng.random() < 0.5:
                        m = al.revcomp(m)
                    m, _s, _i = sim.apply_errors(rng, m, model)
                    fh.write(f"@{sid}.{k}\n{m}\n+\n{'5' * len(m)}\n")
                    k += 1
        single.append((sid, [path]))
    cfg = dict(sim.DEFAULTS)
    cfg.update(al.load_assay(fx.ASSAY))
    run = os.path.join(tmp, "run_il")
    # write_illumina writes the fixture's samples in its own order.
    fx.write_illumina(run, [{"sample_id": s} for s, _t, _c in fx.SAMPLES],
                      species, cfg, random.Random(9))
    pairs = [(sid, os.path.join(run, f"{sid}_R1.fastq.gz"),
              os.path.join(run, f"{sid}_R2.fastq.gz"),
              ["--read-length", "150"]) for sid, _t, _c in picked]
    types = {s: t for s, t, _c in picked}
    return {"release": rel, "primers": prim, "types": types,
            "single_long_noisy": single, "paired_short": pairs}


def nemits2_fixture(tmp: str) -> dict:
    """tests/test_nemits2.py's miniature ITS2 release, its nanopore
    sample and its pairs at 2x150 and 2x250 (behind a 3-base spacer)."""
    recs = nem.build_records(random.Random(20261110))
    rel = os.path.join(tmp, "release")
    nem.write_release(rel, recs)
    prim = os.path.join(tmp, "primers.tsv")
    nem.write_primers(prim)
    seqs = nem.source_seqs(recs)
    fq = os.path.join(tmp, "N.fastq.gz")
    nem.write_nanopore(fq, seqs, random.Random(5))
    pairs = []
    for L in (150, 250):
        r1, r2 = (os.path.join(tmp, f"S{L}_R{k}.fastq.gz") for k in (1, 2))
        nem.write_pairs(r1, r2, seqs, L, random.Random(L))
        pairs.append((f"S{L}", r1, r2, ["--read-length", str(L),
                                         "--primer-offset", "3"]))
    return {"release": rel, "primers": prim,
            "types": {"N": "unknown", "S150": "unknown", "S250": "unknown"},
            "single_long_noisy": [("N", [fq])], "paired_short": pairs}


FIXTURES = {"vespa_v4": vespa_fixture, "nemits2": nemits2_fixture}


def run_route(proto: str, read_type: str, fx_: dict, item, outdir: str,
              mm2: str, extra=()) -> int:
    if read_type == "paired_short":
        sid, r1, r2, more = item
        args = ["--r1", r1, "--r2", r2, "--assay", fx.ASSAY, *more]
        main = adp.main
    else:
        sid, paths = item
        args = ["--reads", *paths, "--read-type", read_type]
        main = adb.main
    return quiet(main, args + [
        "--database", fx_["release"], "--primers", fx_["primers"],
        "--sample-id", sid, "--protocol", proto, "--outdir", outdir,
        "--minimap2", mm2, *extra])


def read_tsv(path: str) -> tuple[list[str], list[dict]]:
    with open(path, encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
    return list(rdr.fieldnames or []), rows


# --------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------


def test_registry() -> None:
    print("\nthe protocol registry")
    for name in protocols.PROTOCOLS:
        p = protocols.get(name)
        check(f"{name}: loads and knows its name", p.name == name)
        covered = set(p.read_types) | set(p.not_supported)
        check(f"{name}: says what it does with every read type, once each",
              covered == set(protocols.READ_TYPES)
              and not set(p.read_types) & set(p.not_supported),
              f"{sorted(p.read_types)} / {sorted(p.not_supported)}")
        check(f"{name}: its extra columns are its own",
              not set(p.extra_columns) & set(al.CLUSTER_COLUMNS))
        for rt in p.not_supported:
            check(f"{name}: says why not {rt}",
                  rt in p.unsupported(rt) and p.not_supported[rt]
                  in p.unsupported(rt))
    # Protocols never import one another.
    for name in protocols.PROTOCOLS:
        src = open(os.path.join(ROOT, "bin", "protocols", f"{name}.py")).read()
        others = [o for o in protocols.PROTOCOLS if o != name
                  and (f"protocols.{o}" in src or f"import {o}" in src)]
        check(f"{name}: imports no other protocol", not others, str(others))
    pfor = protocols.protocol_for
    v4 = {"version": "x"}
    its2 = {"score_over": "reference span"}
    check("a manifest with no protocol and no span scoring is vespa_v4",
          pfor(v4) == "vespa_v4")
    check("a release scored over its records' span is nemits2",
          pfor(its2) == "nemits2")
    check("--db_calling_rule gap on a V4 release is span scoring (nemits2), "
          "as it always was", pfor(v4, "auto", "gap") == "nemits2")
    check("a manifest calling_rule gap is nemits2",
          pfor({"calling_rule": "gap"}) == "nemits2")
    check("the manifest's own protocol wins over inference",
          pfor({**its2, "protocol": "vespa_v4"}) == "vespa_v4")
    check("--db_protocol wins over the manifest",
          pfor({"protocol": "vespa_v4"}, "nemits2") == "nemits2")
    try:
        pfor({"protocol": "nope"})
        bad = False
    except ValueError as exc:
        bad = "nope" in str(exc) and "vespa_v4" in str(exc)
    check("an unknown protocol is refused, naming the ones there are", bad)
    check("the default organism groups are the hook's default",
          protocols.get("vespa_v4").organism_group(
              ["Eukaryota", "Obazoa", "Opisthokonta", "Metazoa", "Nematoda"],
              "Craniata") == "helminths")


def test_contract_b(fixtures: dict, mm2: str, tmp: str) -> dict:
    print("\ncontract B: every protocol on its fixture, every read type")
    outs: dict = {}
    for name, fx_ in fixtures.items():
        p = protocols.get(name)
        for rt in protocols.READ_TYPES:
            items = fx_.get(rt) or []
            if rt not in p.read_types:
                # Refused, with the protocol's reason.
                item = (fx_["single_long_noisy"][0] if rt.startswith("single")
                        else fx_["paired_short"][0])
                if rt.startswith("single"):
                    item = (item[0], item[1])
                try:
                    run_route(name, rt, fx_, item, os.path.join(tmp, "x"), mm2)
                    msg = ""
                except SystemExit as exc:
                    msg = str(exc)
                check(f"{name}: {rt} is refused with its reason",
                      p.not_supported.get(rt, "@") in msg, msg)
                continue
            check(f"{name}: has a fixture for {rt}", bool(items))
            out = os.path.join(tmp, f"{name}_{rt}")
            for item in items:
                sid = item[0]
                code = run_route(name, rt, fx_, item, out, mm2)
                head, rows = read_tsv(os.path.join(out, f"{sid}.clusters.tsv"))
                thead, trows = read_tsv(os.path.join(out, f"{sid}.taxa.tsv"))
                js = json.load(open(os.path.join(out, f"{sid}.assign_db.json")))
                fasta = open(os.path.join(out, f"{sid}.consensus.fasta")).read()
                tag = f"{name} {rt} {sid}"
                check(f"{tag}: runs", code == 0)
                check(f"{tag}: clusters.tsv has the required columns in "
                      "order, then the protocol's declared extras",
                      head == al.CLUSTER_COLUMNS + list(p.extra_columns),
                      str(head))
                check(f"{tag}: taxa.tsv has exactly its columns",
                      thead == al.TAXA_COLUMNS, str(thead))
                check(f"{tag}: assign_db.json has every required key",
                      all(k in js for k in al.ASSIGN_SUMMARY_KEYS),
                      str(sorted(set(al.ASSIGN_SUMMARY_KEYS) - set(js))))
                check(f"{tag}: assign_db.json names the protocol",
                      js["protocol"] == name, str(js.get("protocol")))
                r = js["reads"]
                check(f"{tag}: the read fates are there and balance",
                      all(k in r for k in al.ASSIGN_READ_FATES)
                      and r["reads_in"] == r["no_primer"] + r["too_short"]
                      + r.get("low_quality", 0) + r["clustered"]
                      + r["unclustered"] and r["reads_in"] > 0, str(r))
                reads = [int(x["reads"]) for x in rows]
                check(f"{tag}: clustered is the clusters' reads, most first",
                      sum(reads) == r["clustered"]
                      and reads == sorted(reads, reverse=True)
                      and js["clusters"] == len(rows))
                check(f"{tag}: relative abundances add to one",
                      abs(sum(float(x["rel_abundance"]) for x in rows) - 1)
                      < 1e-3 if rows else True)
                check(f"{tag}: taxa.tsv holds the clusters' reads",
                      sum(int(x["reads"]) for x in trows) == sum(reads))
                check(f"{tag}: one consensus record per cluster, in order",
                      [ln[1:].split()[0] for ln in fasta.splitlines()
                       if ln.startswith(">")] == [x["cluster"] for x in rows])
                ranks = al.load_taxonomy(os.path.join(
                    fx_["release"], "taxonomy.tsv")).ranks
                check(f"{tag}: every rank is one of the release's",
                      all(x["rank"] in ranks + [""] for x in rows))
                check(f"{tag}: every call names the release's version",
                      {x["database_version"] for x in rows}
                      <= {js["database"]["version"]})
                if p.extra_columns:
                    check(f"{tag}: the declared extras are filled",
                          any(x[c] for x in rows for c in p.extra_columns))
            outs[(name, rt)] = out
    return outs


def preflight_dir(tmp: str, types: dict, plate: bool = True) -> str:
    """A minimal pre-flight directory for synthesis."""
    d = os.path.join(tmp, "pf")
    os.makedirs(d, exist_ok=True)
    rows = "run_id\tsample_id\tsanitised_id\tplatform\tplate_row\tplate_col\tsample_type\n"
    for k, (sid, t) in enumerate(sorted(types.items())):
        r, c = ("A", k + 1) if plate else ("", "")
        rows += f"R1\t{sid}\t{sid}\tnanopore\t{r}\t{c}\t{t}\n"
    open(os.path.join(d, "file_resolution.tsv"), "w").write(rows)
    json.dump({}, open(os.path.join(d, "preflight_summary.json"), "w"))
    return d


def test_db_calls(fixtures: dict, outs: dict, tmp: str) -> None:
    print("\ndb_calls.tsv: synthesis over each protocol's output")
    for name, fx_ in fixtures.items():
        p = protocols.get(name)
        out = outs.get((name, "single_long_noisy"))
        if not out:
            continue
        pfd = preflight_dir(os.path.join(tmp, f"syn_{name}"), {
            sid: fx_["types"][sid] for sid, _ in fx_["single_long_noisy"]})
        syn = os.path.join(tmp, f"syn_{name}", "out")
        code = quiet(sdb.main, ["--assign-dir", out, "--preflight-dir", pfd,
                                "--database", fx_["release"], "--run-id", "R1",
                                "--protocol", name, "--outdir", syn,
                                "--min-depth", "10", "--permutations", "50"])
        head, rows = read_tsv(os.path.join(syn, "db_calls.tsv"))
        opt = [c for c in head if c not in al.DB_CALL_COLUMNS]
        check(f"{name}: synthesis runs", code == 0)
        check(f"{name}: db_calls.tsv has the shared columns in order",
              head[:len(al.DB_CALL_COLUMNS)] == al.DB_CALL_COLUMNS, str(head))
        check(f"{name}: and only declared optional ones after them, in order",
              opt == [c for c in al.DB_CALL_OPTIONAL if c in opt], str(opt))
        want = [c for c in al.DB_CALL_OPTIONAL if c in p.extra_columns]
        check(f"{name}: optional columns exactly where its clusters carry them",
              opt == want, f"{opt} vs {want}")
        qc = json.load(open(os.path.join(syn, "db_run_qc.json")))
        check(f"{name}: the protocol's hooks set the defaults",
              qc["settings"]["host_clade"] == p.host_clade
              and qc["settings"]["min_rel_abundance"] == p.min_rel_abundance)
        check(f"{name}: every call is in a known group",
              {r["group"] for r in rows} <= set(sdb.GROUPS))


def test_contract_a(tmp: str) -> None:
    print("\ncontract A: the demux accounting")
    declared = [
        {"run_id": "R1", "sample_id": "S 1", "sanitised_id": "S_1",
         "plate_row": "A", "plate_col": "1", "sample_type": "unknown"},
        {"run_id": "R1", "sample_id": "B", "sanitised_id": "B",
         "plate_row": "A", "plate_col": "2", "sample_type": "neg_extraction"},
        {"run_id": "R1", "sample_id": "X", "sanitised_id": "X",
         "plate_row": "", "plate_col": "", "sample_type": "unknown"},
        {"run_id": "R2", "sample_id": "other", "sanitised_id": "other",
         "plate_row": "B", "plate_col": "1", "sample_type": "unknown"},
    ]
    demux = {"total_reads": 100, "assigned_reads": 90,
             "by_tier": {"A": 90, "B": 6, "C": 2, "D": 1, "E": 1},
             "per_sample": {"S_1": 80, "B": 10}}
    split = {"reads_in": 85, "pieces_out": 100}
    a = acc.accounting("R1", "nanopore", "single_long_noisy", declared,
                       demux, split)
    check("nanopore: every fate counted, in the shared names, and balanced",
          a["fates"] == {"assigned": 90, "one_end_only": 6, "chimeric": 2,
                         "unmatched": 1, "internal_adapter": 1}
          and a["balanced"] is True and a["reads_in"] == 100
          and a["reads_off_instrument"] == 85, str(a["fates"]))
    check("nanopore: reads per declared sample, by declared or sanitised id",
          a["per_sample"] == {"S 1": 80, "B": 10, "X": 0}, str(a["per_sample"]))
    check("plate positions where the platform has them, this run only",
          a["plate"] == {"S 1": {"row": "A", "col": 1},
                         "B": {"row": "A", "col": 2}}, str(a["plate"]))
    check("the platform's own summaries ride along unchanged",
          a["native"] == {"demux": demux, "split": split})
    i = acc.accounting("R1", "illumina", "paired_short", declared, None, None)
    check("illumina: what the sequencer discarded is unknown, not zero",
          i["fates"] == {} and i["reads_in"] is None and i["balanced"] is None
          and set(i["per_sample"].values()) == {None})
    p = acc.accounting("R1", "predemultiplexed", "single_long_accurate",
                       declared, None, None)
    check("predemultiplexed: demultiplexed upstream, and says so",
          "upstream" in p["demultiplexed_by"] and p["fates"] == {})

    # Synthesis through the accounting is synthesis through the summaries.
    mm2 = minimap2()
    if not mm2:
        return
    fx_ = vespa_fixture(os.path.join(tmp, "a_fx"))
    out = os.path.join(tmp, "a_assign")
    for item in fx_["single_long_noisy"]:
        run_route("vespa_v4", "single_long_noisy", fx_, item, out, mm2)
    pfd = preflight_dir(os.path.join(tmp, "a_syn"), fx_["types"])
    demux = {"total_reads": 999, "assigned_reads": 900,
             "by_tier": {"A": 900, "B": 99},
             "tier_meaning": {"B": "one end only - truncated"},
             "per_sample": {s: 1 for s in fx_["types"]}}
    split = {"reads_in": 950, "pieces_out": 999}
    paths = {}
    for n, obj in (("demux", demux), ("split", split)):
        paths[n] = os.path.join(tmp, f"{n}.json")
        json.dump(obj, open(paths[n], "w"))
    a = acc.accounting("R1", "nanopore", "single_long_noisy",
                       al.read_table(os.path.join(pfd, "file_resolution.tsv")),
                       demux, split)
    paths["acc"] = os.path.join(tmp, "acc.json")
    json.dump(a, open(paths["acc"], "w"))
    common = ["--assign-dir", out, "--preflight-dir", pfd, "--database",
              fx_["release"], "--run-id", "R1", "--min-depth", "10",
              "--permutations", "50"]
    quiet(sdb.main, common + ["--demux-summary", paths["demux"],
                              "--split-summary", paths["split"],
                              "--outdir", os.path.join(tmp, "via_summaries")])
    quiet(sdb.main, common + ["--accounting", paths["acc"],
                              "--outdir", os.path.join(tmp, "via_accounting")])
    same = all(open(os.path.join(tmp, "via_summaries", f), "rb").read() ==
               open(os.path.join(tmp, "via_accounting", f), "rb").read()
               for f in ("db_calls.tsv", "db_wells.tsv", "db_prevalence.tsv",
                         "db_run_qc.json"))
    check("synthesis through the accounting is byte for byte synthesis "
          "through the platform's summaries", same)


def test_predemultiplexed(tmp: str) -> None:
    print("\nthe predemultiplexed adapter, in pre-flight")
    species = fx.build_tree(random.Random(3))
    rel = os.path.join(tmp, "release")
    fx.write_release(rel, species)
    prim = os.path.join(tmp, "primers.tsv")
    fx.write_primers(prim)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    sheet = mf.build_samplesheet(i7, i5, 4, 2, "P1", "predemultiplexed")
    ids = ["S1", "S10", "lima", "dorado"]
    for row, sid in zip(sheet, ids):
        row["sample_id"] = sid
    sheet[2]["file_prefix"] = "m84001.bc1001--bc1001"
    sheet[3]["index_i7"] = sheet[3]["index_i5"] = ""
    run = os.path.join(tmp, "run")
    os.makedirs(os.path.join(run, "dorado"))
    for f in ("S1.fastq.gz", "S10_part1.fastq.gz", "S10_part2.fastq.gz",
              "m84001.bc1001--bc1001.fastq.gz", "dorado/a.fastq",
              "dorado/b.fq.gz", "notes.txt"):
        open(os.path.join(run, f), "w").write("@r\nACGT\n+\nIIII\n")
    assay = os.path.join(ROOT, "assets", "assay", "vespa18s.json")

    def go(label, rows, read_type="", extra=(), database=True):
        d = os.path.join(tmp, label)
        os.makedirs(d, exist_ok=True)
        ss = os.path.join(d, "samplesheet.tsv")
        cols = list(rows[0]) + (["file_prefix"] if "file_prefix" not in rows[0] else [])
        with open(ss, "w") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n",
                               restval="")
            w.writeheader()
            w.writerows(rows)
        runs = os.path.join(d, "runs.tsv")
        with open(runs, "w") as fh:
            fh.write("run_id\tplatform\tpath\tread_type\tdna_input_ng\n"
                     f"P1\tpredemultiplexed\t{run}\t{read_type}\t10\n")
        args = ["--samplesheet", ss, "--runs", runs, "--index-map", INDEX_MAP,
                "--assay", assay, "--primers", prim,
                "--outdir", os.path.join(d, "out"), *extra]
        args += (["--reference-database", rel] if database else
                 ["--reference", os.path.join(ROOT, "tests", "data",
                                              "reference.tsv")])
        code = quiet(pf.main, args)
        summ = json.load(open(os.path.join(d, "out", "preflight_summary.json")))
        return code, summ, os.path.join(d, "out")

    code, summ, out = go("ok", sheet)
    errs = [e for e in summ["errors"] if "minimap2" not in e]
    check("a predemultiplexed run passes pre-flight", not errs, str(errs))
    _h, res = read_tsv(os.path.join(out, "file_resolution.tsv"))
    files = {r["sample_id"]: [os.path.basename(f) for f in
                              r["files"].split(";") if f] for r in res}
    check("one file per sample, by its id", files["S1"] == ["S1.fastq.gz"],
          str(files))
    check("chunks of one sample are all its own, and S1 does not take S10's",
          files["S10"] == ["S10_part1.fastq.gz", "S10_part2.fastq.gz"],
          str(files))
    check("a lima output found by file_prefix",
          files["lima"] == ["m84001.bc1001--bc1001.fastq.gz"], str(files))
    check("a folder per sample (Dorado, MinKNOW barcodeNN)",
          files["dorado"] == ["a.fastq", "b.fq.gz"], str(files))
    _h, rts = read_tsv(os.path.join(out, "read_types.tsv"))
    check("the run's read type is written for the workflow (contract A)",
          rts == [{"run_id": "P1", "platform": "predemultiplexed",
                   "read_type": "single_long_noisy"}], str(rts))
    check("with no read_type it says it took single_long_noisy",
          any("taken as single_long_noisy" in w for w in summ["warnings"]))
    check("a well with no indices is warned of: no plate position",
          any("no plate position" in w for w in summ["warnings"]))
    check("pre-flight records the protocol it resolved",
          summ["reference_database"].get("protocol") == "vespa_v4")

    code, summ, _ = go("missing", sheet + [dict(sheet[0], sample_id="S2")])
    check("a sample with no FASTQ is an error",
          any("'S2': no FASTQ matched" in e for e in summ["errors"]),
          str(summ["errors"]))
    code, summ, _ = go("hifi", sheet, read_type="single_long_accurate")
    check("a read type the protocol cannot take is refused before any read "
          "is touched", any("cannot take single_long_accurate" in e
                            for e in summ["errors"]), str(summ["errors"]))
    code, summ, _ = go("pairs", sheet, read_type="paired_short")
    check("pairs are an Illumina run, not a predemultiplexed one",
          any("Illumina run" in e for e in summ["errors"]), str(summ["errors"]))
    code, summ, _ = go("proto", sheet, extra=["--db-protocol", "nope"])
    check("an unknown --db_protocol is refused",
          any("protocol 'nope'" in e for e in summ["errors"]),
          str(summ["errors"]))
    code, summ, _ = go("gap", sheet, extra=["--db-protocol", "vespa_v4",
                                            "--db-calling-rule", "gap"])
    check("a calling rule the protocol does not use is refused",
          any("protocol vespa_v4 calls by floors" in e
              for e in summ["errors"]), str(summ["errors"]))
    code, summ, _ = go("panel", sheet, database=False)
    check("predemultiplexed is database mode only",
          any("database mode only" in e for e in summ["errors"]),
          str(summ["errors"]))


def main() -> int:
    mm2 = minimap2()
    test_registry()
    tmp = tempfile.mkdtemp(prefix="contracts_")
    try:
        if not mm2:
            print("\n  SKIP: minimap2 not found - protocols cannot be run")
            FAIL.append(("minimap2 present", "not on PATH or in ~/bin"))
        else:
            fixtures = {n: f(os.path.join(tmp, n)) for n, f in FIXTURES.items()}
            check("every protocol has a fixture",
                  set(fixtures) == set(protocols.PROTOCOLS))
            outs = test_contract_b(fixtures, mm2, tmp)
            test_db_calls(fixtures, outs, tmp)
        test_contract_a(os.path.join(tmp, "a"))
        test_predemultiplexed(os.path.join(tmp, "p"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\n" + "=" * 70)
    print(f"{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        for name, detail in FAIL:
            print(f"  {name}: {detail}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
