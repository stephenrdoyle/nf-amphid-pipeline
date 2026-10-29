#!/usr/bin/env python3
"""
Database mode's Illumina route against VESPA's own published results.

Not a test: a benchmark on public data, kept so it can be rerun on a new
release. VESPA (Owens et al. 2024, Nat Commun 15:402) deposited the
reads of its 12 human and 40 non-human primate faecal samples (BioProject
PRJNA944233, MiSeq 2x301) and gave its per-sample results in the paper's
Source Data file: each sample's reads by organism category (Figs 4a, 5a),
its helminth composition (4b, 5b) and its protozoa composition (4c, 5c).
The mock communities' reads were not deposited, so this measures
AGREEMENT with VESPA's DADA2 + assignTaxonomy calls, not accuracy
against a known input.

    python3 tests/bench_vespa_public.py \\
        --reads <dir of SRR*_1/_2.fastq.gz> --filereport ena_filereport.tsv \\
        --source-data 41467_2023_44521_MOESM5_ESM.xlsx \\
        --release <release dir> --index refdb.k19w10.mmi \\
        --windows 132x128,250x200 --outdir <dir>

132x128 are the release's 2x150 windows: what the same samples would
give on the 2x150 runs this pipeline is set up for. 250x200 uses most of
the 2x301 reads.

The NHP sample names in the Source Data are NOT the BioProject's: NHP11
here matches the paper's "NHP39" to within 1% in every category, and its
own namesake not at all (the human names do match). So the NHP
correspondence is inferred: each BioProject sample is paired one to one
with the paper's sample whose category profile, helminth genera and
protozoa genera fit it best (a greedy matching over all pairs, closest
first), and the pairing is written out with its distance - a pair whose
distance is large is a guess, and is flagged. The pairing is decided ONCE
per invocation - on the first window set, or --pair-on - and every window
set is scored with it; paired per window set, half the NHP samples changed
partner between 132x128 and 250x200, and the "difference" between window
sets was the re-pairing. --nhp-mapping takes a pairing from an earlier
run, so releases are scored against the same partners too.

Comparison is at GENUS: VESPA names some taxa to species and some to
genus, and this pipeline stops where its floors say; a genus is the
deepest level both can be read at for every taxon. A VESPA taxon counts
as found when any detected call in that sample has the same genus in its
lineage - including a call that stopped at genus.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import assign_db_pairs as adp   # noqa: E402
import synthesise_db as syn     # noqa: E402

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}

# VESPA's categories, from ours. VESPA's 'Protozoa' panel includes
# Blastocystis and the microsporidian Enterocytozoon; its NHP figure
# splits environmental from parasitic protozoa, pooled here.
CATEGORY = {"helminths": "Helminths", "protists": "Protozoa",
            "microsporidia": "Protozoa", "host": "Host", "fungi": "Fungi",
            "plants": "Plants + Invertebrates",
            "invertebrates": "Plants + Invertebrates",
            "bacteria_archaea": "Bacteria + Archaea",
            "organelle": "Bacteria + Archaea"}
CATS = ["Bacteria + Archaea", "Host", "Plants + Invertebrates", "Fungi",
        "Protozoa", "Helminths"]


def xlsx_sheets(path: str) -> dict[str, list[list[str]]]:
    """Every sheet of an .xlsx as rows of strings, standard library only."""
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            ss.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target") for r in
            ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for sh in wb.find("m:sheets", NS):
        tgt = rels[sh.get(f"{{{NS['r']}}}id")].lstrip("/").replace("xl/", "")
        rows = []
        for r in ET.fromstring(z.read("xl/" + tgt)).iter(f"{{{NS['m']}}}row"):
            vals = {}
            for c in r.findall("m:c", NS):
                v = c.find("m:v", NS)
                col = 0
                for ch in re.match(r"[A-Z]+", c.get("r")).group():
                    col = col * 26 + ord(ch) - 64
                if v is None:
                    continue
                vals[col - 1] = ss[int(v.text)] if c.get("t") == "s" else v.text
            if vals:
                rows.append([vals.get(i, "") for i in range(max(vals) + 1)])
        out[sh.get("name")] = rows
    return out


def vespa_tables(path: str) -> dict:
    """{'groups'|'helminths'|'protozoa': {sample: {name: share}}} from the
    Source Data's Figure 4 and Figure 5 sheets. Sample names normalised
    to the BioProject's (Human1 -> Human01)."""
    def norm(s):
        m = re.match(r"(Human|NHP)0*(\d+)$", s.strip())
        return f"{m.group(1)}{int(m.group(2)):02d}" if m and m.group(1) == "Human" \
            else (f"NHP{int(m.group(2))}" if m else s.strip())
    out = {"groups": {}, "helminths": {}, "protozoa": {}}
    sheets = xlsx_sheets(path)
    for fig in ("Figure 4", "Figure 5"):
        rows = sheets.get(fig, [])
        kind, head = None, []
        for r in rows:
            tag = (r[0] if r else "").strip()
            if tag:
                kind = ("groups" if tag.endswith("a") and len(tag) <= 3 else
                        "helminths" if "b" in tag else "protozoa")
                head = r
                continue
            if kind and len(r) > 2 and r[1].strip():
                vals = {}
                for name, v in zip(head[2:], r[2:]):
                    name = name.strip()
                    if not name or name == "TOTAL" or name in NOT_A_TAXON:
                        continue
                    try:
                        vals[name] = float(v)
                    except ValueError:
                        pass
                out[kind][norm(r[1])] = vals
    # Fig 5 splits environmental and parasitic protozoa; pool them.
    for s, v in out["groups"].items():
        if "Environmental protozoa" in v or "Parasitic protozoa" in v:
            v["Protozoa"] = v.pop("Environmental protozoa", 0) + \
                v.pop("Parasitic protozoa", 0)
    return out


# Columns of VESPA's composition tables that are bins, not taxa.
NOT_A_TAXON = {"Environmental prots", "Environmental protozoa"}


def genus_of(name: str) -> str:
    name = name.replace("Paragonimu safricanus", "Paragonimus africanus")
    return name.replace("_", " ").split()[0].strip(".")


def run_one(job):
    sid, r1, r2, a = job
    out = os.path.join(a["outdir"], a["tag"])
    if not os.path.exists(os.path.join(out, f"{sid}.clusters.tsv")):
        with contextlib.redirect_stdout(io.StringIO()):
            adp.main(["--r1", r1, "--r2", r2, "--database", a["release"],
                      "--index", a["index"], "--primers", a["primers"],
                      "--assay", a["assay"], "--sample-id", sid,
                      "--outdir", out, "--window-r1", str(a["w1"]),
                      "--window-r2", str(a["w2"]),
                      "--read-length", str(a["read_length"])])
    return sid


def ours(outdir: str, sid: str, host_clade: str, min_rel: float) -> dict:
    with open(os.path.join(outdir, f"{sid}.clusters.tsv")) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    groups, genera = Counter(), defaultdict(float)
    tot = sum(int(r["reads"]) for r in rows) or 1
    for r in rows:
        g = syn.organism_group(r["lineage"].split(";") if r["lineage"] else [],
                               host_clade)
        groups[CATEGORY.get(g, "unassigned")] += int(r["reads"])
        lin = r["lineage"].split(";")
        if len(lin) > 7 and lin[7] and int(r["reads"]) / tot >= min_rel:
            genera[lin[7]] += int(r["reads"]) / tot
    return {"groups": groups, "genera": genera, "total": tot}


def profile_distance(o, sid_v, vespa) -> float:
    """How unlike a BioProject sample (our result) and a paper sample are:
    total absolute difference in category shares, plus the share of the
    helminth and protozoa genera found by one and not the other."""
    v = vespa["groups"].get(sid_v, {})
    tot = sum(o["groups"][c] for c in CATS) or 1
    d = sum(abs(o["groups"][c] / tot - v.get(c, 0)) for c in CATS)
    for panel in ("helminths", "protozoa"):
        vg = {genus_of(n) for n, x in vespa[panel].get(sid_v, {}).items() if x > 0}
        og = set(o["genera"])
        if vg:
            d += len(vg - og) / len(vg) * 0.5
    return d


def remap(runs, vespa, out, min_rel) -> dict:
    """BioProject NHP sample -> (paper NHP sample, distance), one to one,
    closest pairs first. Human samples keep their names."""
    mine = {sid: ours(out, sid, "Craniata", min_rel) for sid in runs
            if sid.startswith("NHP")}
    theirs = [s for s in vespa["groups"] if s.startswith("NHP")]
    pairs = sorted((profile_distance(o, t, vespa), sid, t)
                   for sid, o in mine.items() for t in theirs)
    got, used = {}, set()
    for d, sid, t in pairs:
        if sid in got or t in used:
            continue
        got[sid] = (t, d)
        used.add(t)
    for sid in runs:
        got.setdefault(sid, (sid, 0.0))
    return got


def read_mapping(path: str) -> dict:
    """bioproject sample -> (paper sample, distance), as written below."""
    with open(path) as fh:
        return {r["bioproject_sample"]: (r["paper_sample"], float(r["distance"]))
                for r in csv.DictReader(fh, delimiter="\t")}


def write_mapping(path: str, got: dict) -> None:
    with open(path, "w") as fh:
        fh.write("bioproject_sample\tpaper_sample\tdistance\n")
        for sid in sorted(got, key=lambda x: (x[:3], int(re.sub(r"\D", "", x) or 0))):
            fh.write(f"{sid}\t{got[sid][0]}\t{got[sid][1]:.3f}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", required=True)
    ap.add_argument("--filereport", required=True)
    ap.add_argument("--source-data", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--primers", default=os.path.join(ROOT, "tests/data/db/primers.tsv"))
    ap.add_argument("--assay", default=os.path.join(ROOT, "assets/assay/vespa18s.json"))
    ap.add_argument("--windows", default="132x128,250x200")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--read-length", type=int, default=300,
                    help="read length of the deposited runs (2x301 MiSeq: "
                    "the release's 300x300 table)")
    ap.add_argument("--min-rel", type=float, default=0.001)
    ap.add_argument("--no-remap", action="store_true",
                    help="trust the NHP names as given")
    ap.add_argument("--pair-on", help="window set the NHP pairing is decided "
                    "on, e.g. 132x128 (default: the first of --windows)")
    ap.add_argument("--nhp-mapping", help="a pairing written by an earlier "
                    "run (nhp_mapping.tsv), used as is")
    args = ap.parse_args(argv)

    ap_remap = not args.no_remap
    with open(args.filereport) as fh:
        runs = {r["sample_alias"]: r["run_accession"]
                for r in csv.DictReader(fh, delimiter="\t")}
    vespa = vespa_tables(args.source_data)
    summary = {}
    specs = args.windows.split(",")
    for spec in specs:
        w1, w2 = (int(x) for x in spec.split("x"))
        tag = f"w{w1}x{w2}"
        a = {"outdir": args.outdir, "tag": tag, "release": args.release,
             "index": args.index, "primers": args.primers,
             "assay": args.assay, "w1": w1, "w2": w2,
             # 2x150-equivalent windows are scored against the 2x150
             # resolution table, the long ones against the reads' own.
             "read_length": 150 if (w1, w2) == (132, 128) else args.read_length}
        jobs = [(sid, os.path.join(args.reads, f"{acc}_1.fastq.gz"),
                 os.path.join(args.reads, f"{acc}_2.fastq.gz"), a)
                for sid, acc in sorted(runs.items())]
        with ProcessPoolExecutor(args.jobs) as ex:
            for sid in ex.map(run_one, jobs):
                pass
    # One pairing for every window set (see the docstring).
    name = {sid: sid for sid in runs}
    if args.nhp_mapping:
        got = read_mapping(args.nhp_mapping)
        got = {sid: got.get(sid, (sid, 0.0)) for sid in runs}
        write_mapping(os.path.join(args.outdir, "nhp_mapping.tsv"), got)
        name = {sid: v[0] for sid, v in got.items()}
    elif ap_remap:
        pair_on = args.pair_on or specs[0]
        if pair_on not in specs:
            ap.error(f"--pair-on {pair_on} is not one of --windows")
        got = remap(runs, vespa, os.path.join(args.outdir, f"w{pair_on}"),
                    args.min_rel)
        write_mapping(os.path.join(args.outdir, "nhp_mapping.tsv"), got)
        name = {sid: v[0] for sid, v in got.items()}
    for spec in specs:
        tag = "w" + spec
        out = os.path.join(args.outdir, tag)
        rows, agree = [], Counter()
        missed: dict = defaultdict(Counter)
        diffs = defaultdict(list)
        for sid in sorted(runs):
            o = ours(out, sid, "Craniata", args.min_rel)
            v = vespa["groups"].get(name[sid])
            tot = sum(o["groups"][c] for c in CATS) or 1
            row = {"sample": sid, "paper_sample": name[sid]}
            for c in CATS:
                row[f"ours {c}"] = o["groups"][c] / tot
                row[f"vespa {c}"] = (v or {}).get(c, "")
                if v is not None:
                    diffs[c].append(o["groups"][c] / tot - v.get(c, 0))
            row["ours unassigned"] = o["groups"]["unassigned"] / o["total"]
            for panel in ("helminths", "protozoa"):
                vt = vespa[panel].get(name[sid], {})
                vg = {genus_of(n) for n, x in vt.items() if x > 0}
                og = set(o["genera"])
                both = vg & og
                # Human names are the paper's own; NHP names are inferred
                # (remap), so the two cohorts are scored apart.
                cohort = "human" if sid.startswith("Human") else "nhp"
                agree[f"{cohort} {panel}: VESPA genera"] += len(vg)
                agree[f"{cohort} {panel}: also found here"] += len(both)
                for gname in vg - og:
                    missed[f"{cohort} {panel}"][gname] += 1
                row[f"{panel} VESPA only"] = ",".join(sorted(vg - og))
                row[f"{panel} found"] = ",".join(sorted(both))
            rows.append(row)
        with open(os.path.join(args.outdir, f"{tag}.comparison.tsv"), "w") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        mad = {c: sum(abs(x) for x in d) / len(d) for c, d in diffs.items() if d}
        summary[tag] = {"samples": len(rows),
                        "mean_abs_diff_in_category_share": mad,
                        "genus_agreement": dict(agree),
                        "missed_genera": {k: dict(v.most_common())
                                          for k, v in missed.items()}}
        print(tag, json.dumps(summary[tag], indent=1))
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
