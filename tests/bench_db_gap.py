#!/usr/bin/env python3
"""
The gap rule's exact-match clause (nemits2 --gap-exact) on a leave-out
simulation of nanopore ITS2 reads.

Not a test: a benchmark, kept because no nanopore NC1/NC2 run was at hand
when GAP_EXACT was taken from the release's own classifier
(its2_classify.py), whose leave-one-out scored reference sequences, not
consensus sequences of reads with errors. The risk it measures: a
consensus of a species missing from the database that is, through the
reads' errors, an exact match to a sister's record 2-3 edits from the
next species - named that sister under --gap-exact 2, at genus without.

As bench_db_floors.py, but for a span-scored release under the gap rule:
species and genera are held OUT of a copy of the release (with any shared
record naming them), reads are simulated from records held out and from
records left in - NC1, the record, the 28S stretch reads carry and
records do not, NC2 - and every sample goes through bin/assign_db.py as
the pipeline runs it, once per --gap-exact value, on the same reads. The
protocol's own call is scored:

  in_db        a record in the index          may reach species
  species_out  its species held out           genus (a congener is in)
  genus_out    its genus held out             family (a confamilial is in)

A call deeper than that, or naming anything wrong at any rank down to it,
is WRONG - unless the release itself cannot separate the source from the
species named, which is INDISTINGUISHABLE: the source record lies within
the gap rule's reach of a typical record of the named species (fewer than
nemits2.GAP_EDITS edits, each record placed inside the other), so with
both species in the database a read of either stops at genus. That takes
in the contained pairs (0 edits: 151 record pairs of 57 species pairs in
v0.2.2, flagged in contained_in / contains) and species that share an
ITS2 to an edit or two (H. sinense / H. liparis). A record far from every
other record of its own species does not count (see typical): it is more
likely misidentified than evidence the two species share an ITS2, and a
call through it is a database error. A wrong call is then a
database error or a classifier error, not a limit of the marker; both
counts are reported, and the pairs behind each.

  python3 tests/bench_db_gap.py --release <nemits2 release> \\
      --workdir /tmp/gap --samples 150 [--families Strongylidae,...]
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import csv
import io
import json
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import ampseqlib as al            # noqa: E402
import assign_db as adb           # noqa: E402
import simulate as sim            # noqa: E402
from protocols import nemits2     # noqa: E402

NC1 = "ACGTCTGGTTCAGGGTTGTT"
NC2 = "TTAGTTTCTTTTCCTCCGCT"
# The 28S stretch between the end of ITS2 and NC2 (tests/test_nemits2.py):
# reads carry it, records do not.
STRETCH = "AACCTGAGCTCAGGCGTGATTACCCGCTGAACTTAAGCATATCATTT"
ALLOWED = {"in_db": "species", "species_out": "genus", "genus_out": "family"}


def source_ok(t: dict) -> bool:
    return (t.get("species_status") == "named"
            and t.get("grade") in ("complete", "near_complete")
            and all(t.get(r) for r in ("family", "genus", "species")))


def prepare(rel: str, work: str, rng: random.Random, n_genera: int,
            n_species: int, families: set[str], mm2: str) -> dict:
    """Hold genera and species out of a copy of the release; index it."""
    seqs = dict(al.read_fasta(os.path.join(rel, "refdb.fasta")))
    man = al.load_release_manifest(rel)
    tax = al.load_taxonomy(os.path.join(rel, "taxonomy.tsv"), man)
    ok = [i for i in seqs if source_ok(tax[i]) and set(seqs[i]) <= set("ACGT")
          and (not families or tax[i]["family"] in families)]
    by_g, by_s = defaultdict(list), defaultdict(list)
    gen_in_fam, sp_in_gen = defaultdict(set), defaultdict(set)
    for i in seqs:
        t = tax[i]
        if t.get("species_status") == "named" and t.get("genus"):
            gen_in_fam[t["family"]].add(t["genus"])
            sp_in_gen[t["genus"]].add(t["species"])
    for i in ok:
        by_g[tax[i]["genus"]].append(i)
        by_s[tax[i]["species"]].append(i)
    gpool = sorted(g for g, ids in by_g.items()
                   if len(gen_in_fam[tax[ids[0]]["family"]]) >= 3)
    drop_g = set(rng.sample(gpool, min(n_genera, len(gpool))))
    spool = sorted(s for s, ids in by_s.items()
                   if tax[ids[0]]["genus"] not in drop_g
                   and len(sp_in_gen[tax[ids[0]]["genus"]]) >= 3)
    drop_s = set(rng.sample(spool, min(n_species, len(spool))))

    def held(i: str) -> bool:
        t = tax[i]
        if t.get("genus") in drop_g or t.get("species") in drop_s:
            return True
        return any(s in drop_s or s.split(" ")[0] in drop_g
                   for s in tax.shared.get(i, ()))
    out = {i for i in seqs if held(i)}
    lo = os.path.join(work, "release_leaveout")
    os.makedirs(lo, exist_ok=True)
    with open(os.path.join(lo, "refdb.fasta"), "w") as fh:
        for i, s in seqs.items():
            if i not in out:
                fh.write(f">{i}\n{s}\n")
    # The taxonomy holds only the records left in: a held-out record
    # named in a kept one's contained_in / contains is not in the
    # database, as a species missing from a real one is not.
    with open(os.path.join(rel, "taxonomy.tsv")) as src, \
            open(os.path.join(lo, "taxonomy.tsv"), "w") as fh:
        for n, line in enumerate(src):
            if n == 0 or line.split("\t", 1)[0] not in out:
                fh.write(line)
    for f in ("members.tsv", "resolution_nanopore.tsv", "manifest.json"):
        dst = os.path.join(lo, f)
        if os.path.exists(os.path.join(rel, f)) and not os.path.exists(dst):
            os.symlink(os.path.join(os.path.abspath(rel), f), dst)
    mmi = os.path.join(work, "leaveout.k19w10.mmi")
    subprocess.run([mm2, *al.INDEX_ARGS, "-d", mmi,
                    os.path.join(lo, "refdb.fasta")], check=True,
                   capture_output=True)
    pools = {
        "in_db": [i for i in ok if i not in out],
        "species_out": [i for i in ok if tax[i]["species"] in drop_s],
        "genus_out": [i for i in ok if tax[i]["genus"] in drop_g],
    }
    return {"pools": pools, "release": lo, "index": mmi,
            "held": {"genera": len(drop_g), "species": len(drop_s),
                     "records": len(out)}}


def one_sample(job: dict) -> list[dict]:
    """Simulate one sample, assign it under each --gap-exact value, and
    match each consensus to the source it came from."""
    rng = random.Random(job["seed"])
    seqs, sid, work = job["seqs"], job["sid"], job["workdir"]
    model = sim.ERROR_MODELS[job["error_model"]]
    fq = os.path.join(work, f"{sid}.fastq.gz")
    with gzip.open(fq, "wt") as fh:
        n = 0
        for t, cat in job["taxa"]:
            for _ in range(rng.choice([20, 50, 150, 400])):
                m = NC1 + seqs[t] + STRETCH + al.revcomp(NC2)
                if rng.random() < 0.5:
                    m = al.revcomp(m)
                m, _s, _i = sim.apply_errors(rng, m, model)
                fh.write(f"@{t}|{cat}|{n}\n{m}\n+\n{'5' * len(m)}\n")
                n += 1
    truth = dict(job["taxa"])
    rows = []
    for ge in job["gap_exact"]:
        out = os.path.join(work, f"out_ge{ge}")
        with contextlib.redirect_stdout(io.StringIO()):
            adb.main(["--reads", fq, "--read-type", job["read_type"],
                      "--database", job["release"], "--index", job["index"],
                      "--primers", job["primers"], "--sample-id", sid,
                      "--outdir", out, "--minimap2", job["minimap2"],
                      "--gap-exact", str(ge)])
        with open(os.path.join(out, f"{sid}.clusters.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                best = min(((al.bounded_edit(r["consensus"], seqs[t] + STRETCH,
                                             60, free_b_suffix=True), t)
                            for t in truth), default=None)
                src, cat = (best[1], truth[best[1]]) \
                    if best and best[0] <= 60 else ("", "spurious")
                rows.append({"sample": sid, "gap_exact": ge,
                             "source": src, "category": cat,
                             **{k: r[k] for k in (
                                 "reads", "rank", "taxon", "lineage", "best_hit",
                                 "identity", "edits", "gap", "origin")}})
        # a source with no consensus: its reads were unmapped or left
        # unclustered
        got = {r["source"] for r in rows if r["gap_exact"] == ge}
        for t, cat in job["taxa"]:
            if t not in got:
                rows.append({"sample": sid, "gap_exact": ge, "source": t,
                             "category": cat, "reads": "0", "rank": "",
                             "taxon": "", "lineage": "", "best_hit": "",
                             "identity": "", "edits": "", "gap": "",
                             "origin": "none"})
    os.remove(fq)
    return rows


def span_edits(a: str, b: str) -> int:
    """Edits between two records, each placed inside the other with its
    overhangs free and the better direction taken (the release's margins);
    a record contained in the other is 0."""
    return min(al.infix_distance(a, b), al.infix_distance(b, a))


_REACH: dict = {}
_TYPICAL: dict = {}


def typical(i: str, tax, seqs: dict) -> bool:
    """A record within the gap of another record of its own species, or
    its species' only record. A record that is not - an H. aduncum ITS2
    labelled H. fabri, 47 edits from every other H. fabri (v0.2.4) - is
    taken for a misidentification, and a call through it stays wrong."""
    if i not in _TYPICAL:
        sp = tax.get(i, {}).get("species", "")
        own = [k for k in seqs if k != i and tax.get(k, {}).get("species") == sp]
        _TYPICAL[i] = not sp or not own or any(
            span_edits(seqs[i], seqs[k]) < nemits2.GAP_EDITS for k in own)
    return _TYPICAL[i]


def within_gap(source: str, r: dict, tax, seqs: dict) -> bool:
    """The release cannot separate the source from what it was named:
    some typical record of the species named (or the best hit, for a call
    above species) is fewer than GAP_EDITS edits from the source record."""
    named = r["taxon"] if r["rank"] == "species" else ""
    key = (source, named or r["best_hit"])
    if key not in _REACH:
        cands = [i for i in seqs if tax.get(i, {}).get("species") == named] \
            if named else [r["best_hit"]]
        s = seqs.get(source, "")
        _REACH[key] = bool(s) and any(
            i in seqs and span_edits(s, seqs[i]) < nemits2.GAP_EDITS
            and typical(i, tax, seqs) for i in cands)
    return _REACH[key]


def outcome(r: dict, tax, ranks: list[str], seqs: dict) -> str:
    if r["category"] == "spurious":
        return "spurious"
    if r["origin"] == "none":
        return "none"
    if r["origin"] == "chimera" or not r["rank"]:
        return "higher"
    depth = ranks.index(r["rank"])
    lim = ranks.index(ALLOWED[r["category"]])
    call = r["lineage"].split(";")
    t = tax[r["source"]]
    if depth > lim or any(call[k] != t[ranks[k]]
                          for k in range(min(depth, lim) + 1)
                          if k < len(call)):
        if r["category"] != "in_db" and within_gap(r["source"], r, tax, seqs):
            return "indistinguishable"
        return "wrong"
    return "deepest" if depth == lim else "higher"


def per_source(rows: list[dict]) -> dict:
    """(sample, source) -> gap_exact -> the source's largest consensus
    (or its no-consensus row)."""
    out: dict = defaultdict(dict)
    for r in rows:
        if r["category"] == "spurious":
            continue
        k = (r["sample"], r["source"])
        cur = out[k].get(r["gap_exact"])
        if cur is None or int(r["reads"]) > int(cur["reads"]):
            out[k][r["gap_exact"]] = r
    return out


def report(rows: list[dict], tax, ranks: list[str], seqs: dict,
           args) -> None:
    src = per_source(rows)
    by = defaultdict(Counter)
    for d in src.values():
        for ge, r in d.items():
            by[(ge, r["category"])][outcome(r, tax, ranks, seqs)] += 1
    print(f"\n{args.error_model}, {args.samples} samples"
          + (f", families {args.families}" if args.families else "")
          + "\nper source, its largest consensus:"
          + "\ngap_exact  category        n  deepest   wrong  indist.  higher"
          "  no consensus")
    for (ge, cat), c in sorted(by.items()):
        n = sum(c.values())
        print(f"{ge:>9}  {cat:12s} {n:4d}  {c['deepest'] / n:7.1%} "
              f"{c['wrong'] / n:6.1%} {c['indistinguishable'] / n:7.1%} "
              f"{c['higher'] / n:7.1%} {c['none'] / n:8.1%}")
    pairs = defaultdict(Counter)
    for d in src.values():
        r = d[max(d)]
        o = outcome(r, tax, ranks, seqs)
        if o in ("wrong", "indistinguishable"):
            pairs[o][(r["category"], tax[r["source"]]["species"],
                      r["taxon"])] += 1
    for o in ("wrong", "indistinguishable"):
        print(f"\n{o} (gap_exact {max(max(d) for d in src.values())}), "
              "source species -> named:")
        for (cat, a, b), n in sorted(pairs[o].items(),
                                     key=lambda kv: (kv[0][0], -kv[1])):
            print(f"  {cat:12s} {n:4d}  {a} -> {b}")
    print("\ncalls the clause changes (gap_exact 0 -> 2):")
    changed = Counter()
    for (_s, s_), d in sorted(src.items()):
        if len(d) < 2:
            continue
        a, b = (d[k] for k in sorted(d))
        src_ = s_
        if (a["rank"], a["taxon"]) != (b["rank"], b["taxon"]):
            oa, ob = outcome(a, tax, ranks, seqs), outcome(b, tax, ranks, seqs)
            changed[(a["category"], oa, ob)] += 1
            print(f"  {a['category']:12s} {tax[src_]['species']:34s} "
                  f"{a['rank']} {a['taxon']} -> {b['rank']} {b['taxon']} "
                  f"(edits {b['edits']}, gap {b['gap']}): {oa} -> {ob}")
    print("  summary:", dict(changed))
    print("spurious clusters:", sum(1 for r in rows
                                    if r["category"] == "spurious") // 2)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--samples", type=int, default=150)
    ap.add_argument("--held-genera", type=int, default=30)
    ap.add_argument("--held-species", type=int, default=100)
    ap.add_argument("--families", default="",
                    help="draw sources from these families only (comma "
                    "separated); the whole release is still the index")
    ap.add_argument("--error-model", default="nanopore_hac",
                    choices=sorted(sim.ERROR_MODELS))
    ap.add_argument("--gap-exact", default="0,2")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--report-only", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.workdir, exist_ok=True)
    man = al.load_release_manifest(args.release)
    tax = al.load_taxonomy(os.path.join(args.release, "taxonomy.tsv"), man)
    ranks = al.ranks_of(tax)
    path = os.path.join(args.workdir, "calls.json")
    seqs = dict(al.read_fasta(os.path.join(args.release, "refdb.fasta")))
    if args.report_only:
        with open(path) as fh:
            report(json.load(fh), tax, ranks, seqs, args)
        return 0
    rng = random.Random(args.seed)
    fams = {f for f in args.families.split(",") if f}
    P = prepare(args.release, args.workdir, rng, args.held_genera,
                args.held_species, fams, args.minimap2)
    print("held out:", P["held"], {k: len(v) for k, v in P["pools"].items()},
          file=sys.stderr)
    primers = os.path.join(args.workdir, "primers.tsv")
    with open(primers, "w") as fh:
        fh.write(f"name\tdirection\tsequence\nNC1\tfwd\t{NC1}\nNC2\trev\t{NC2}\n")
    jobs = []
    for k in range(args.samples):
        taxa = [(i, "in_db") for i in rng.sample(P["pools"]["in_db"], 2)]
        taxa.append((rng.choice(P["pools"]["species_out"]), "species_out"))
        taxa.append((rng.choice(P["pools"]["genus_out"]), "genus_out"))
        jobs.append({"sid": f"G{k:04d}", "seed": rng.random(), "taxa": taxa,
                     "seqs": {t: seqs[t] for t, _c in taxa},
                     "workdir": args.workdir, "release": P["release"],
                     "index": P["index"], "primers": primers,
                     "minimap2": args.minimap2,
                     "error_model": args.error_model,
                     "read_type": "single_long_accurate"
                     if args.error_model.startswith("pacbio")
                     else "single_long_noisy",
                     "gap_exact": [int(x) for x in args.gap_exact.split(",")]})
    rows = []
    with ProcessPoolExecutor(args.workers) as ex:
        for n, got in enumerate(ex.map(one_sample, jobs), 1):
            rows.extend(got)
            if n % 25 == 0:
                print(f"  {n}/{args.samples} samples", file=sys.stderr)
    with open(path, "w") as fh:
        json.dump(rows, fh)
    report(rows, tax, ranks, seqs, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
