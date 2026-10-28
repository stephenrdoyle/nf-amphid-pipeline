#!/usr/bin/env python3
"""
Calibrate database mode's identity floors on a leave-out simulation.

Not a test: a benchmark, kept because it is where the default --db_floors
come from, and a new release or a new chemistry should be able to rerun
it. It needs a real release and minimap2, and takes a while.

The question a floor answers is "how close must a consensus be to its
best hit before we may use that hit's name at this rank?". It cannot be
answered on taxa that are in the database - there every floor looks
fine. So genera and species are held OUT of a copy of the release, reads
are simulated from them and from taxa left in, every sample goes through
bin/assign_db.py exactly as the pipeline runs it but with the floors
switched off, and each consensus's LCA rank and best-hit identity are
recorded. Floors are then applied afterwards, over a grid, to the same
consensus sequences: the trade-off between naming present taxa and
misnaming missing ones is read straight off the table.

  python3 tests/bench_db_floors.py --release <release dir> \\
      --workdir /tmp/floors --samples 150

Categories, and the deepest rank each may honestly reach:
  in_db        a record in the index          species
  species_out  its species held out           genus (a congener is in)
  genus_out    its genus held out             family (a confamilial is in)
A call deeper than that, or naming anything wrong at any rank down to
it, is WRONG.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import io
import json
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import ampseqlib as al            # noqa: E402
import assign_db as adb           # noqa: E402
import simulate as sim            # noqa: E402

ALLOWED = {"in_db": "species", "species_out": "genus", "genus_out": "family"}
FWD = "AGCAGCCGCGGTAATTCC"
REV_TEMPLATE = "TCCGTCAATTCCTTTAAGTTTC"


def full(t):
    return all(t.get(r) for r in ("family", "genus", "species"))


def prepare(rel: str, work: str, rng: random.Random, n_genera: int,
            n_species: int, mm2: str) -> dict:
    """Hold genera and species out of a copy of the release; index it."""
    seqs = dict(al.read_fasta(os.path.join(rel, "refdb.fasta")))
    tax = al.load_taxonomy(os.path.join(rel, "taxonomy.tsv"))
    ok = [i for i in seqs if full(tax[i]) and set(seqs[i]) <= set("ACGT")]
    by_g, by_s = defaultdict(list), defaultdict(list)
    gen_in_fam, sp_in_gen = defaultdict(set), defaultdict(set)
    for i in ok:
        t = tax[i]
        by_g[t["genus"]].append(i)
        by_s[t["species"]].append(i)
        gen_in_fam[t["family"]].add(t["genus"])
        sp_in_gen[t["genus"]].add(t["species"])
    # Genera whose family keeps others, so a family call is possible;
    # species whose genus keeps others, so a genus call is.
    gpool = sorted(g for g, ids in by_g.items()
                   if len(gen_in_fam[tax[ids[0]]["family"]]) >= 3)
    drop_g = set(rng.sample(gpool, min(n_genera, len(gpool))))
    spool = sorted(s for s, ids in by_s.items()
                   if tax[ids[0]]["genus"] not in drop_g
                   and len(sp_in_gen[tax[ids[0]]["genus"]]) >= 3)
    drop_s = set(rng.sample(spool, min(n_species, len(spool))))
    held = {i for g in drop_g for i in by_g[g]} | \
           {i for s in drop_s for i in by_s[s]}
    lo = os.path.join(work, "release_leaveout")
    os.makedirs(lo, exist_ok=True)
    with open(os.path.join(lo, "refdb.fasta"), "w") as fh:
        for i, s in seqs.items():
            if i not in held:
                fh.write(f">{i}\n{s}\n")
    for f in ("taxonomy.tsv", "resolution_nanopore.tsv", "manifest.json"):
        dst = os.path.join(lo, f)
        if not os.path.exists(dst):
            os.symlink(os.path.join(os.path.abspath(rel), f), dst)
    mmi = os.path.join(work, "leaveout.k19w10.mmi")
    subprocess.run([mm2, *adb.INDEX_ARGS, "-d", mmi,
                    os.path.join(lo, "refdb.fasta")], check=True,
                   capture_output=True)
    pools = {
        "in_db": [i for i in ok if i not in held],
        "species_out": [i for i in ok if tax[i]["species"] in drop_s],
        "genus_out": [i for i in ok if tax[i]["genus"] in drop_g],
    }
    return {"seqs": seqs, "tax": tax, "pools": pools, "release": lo,
            "index": mmi, "held_genera": sorted(drop_g),
            "held_species": sorted(drop_s)}


def run(args) -> int:
    rng = random.Random(args.seed)
    os.makedirs(args.workdir, exist_ok=True)
    mm2 = args.minimap2
    P = prepare(args.release, args.workdir, rng, args.held_genera,
                args.held_species, mm2)
    seqs, tax, pools = P["seqs"], P["tax"], P["pools"]
    with open(os.path.join(args.workdir, "primers.tsv"), "w") as fh:
        fh.write(f"name\tdirection\tsequence\nF\tfwd\t{FWD}\n"
                 "R\trev\tTCCGTCAATTYCTTNAASTTTC\n")
    model = dict(sim.ERROR_MODELS[args.error_model])
    rows = []
    for k in range(args.samples):
        sid = f"L{k:04d}"
        taxa = [(i, "in_db") for i in rng.sample(pools["in_db"], 2)]
        taxa.append((rng.choice(pools["species_out"]), "species_out"))
        taxa.append((rng.choice(pools["genus_out"]), "genus_out"))
        fq = os.path.join(args.workdir, f"{sid}.fastq.gz")
        with gzip.open(fq, "wt") as fh:
            n = 0
            for t, cat in taxa:
                for _ in range(rng.choice([20, 50, 150, 400])):
                    m = FWD + seqs[t] + al.revcomp(REV_TEMPLATE)
                    if rng.random() < 0.5:
                        m = al.revcomp(m)
                    m, _s, _i = sim.apply_errors(rng, m, model)
                    fh.write(f"@{t}|{cat}|{n}\n{m}\n+\n{'5' * len(m)}\n")
                    n += 1
        out = os.path.join(args.workdir, "out")
        with contextlib.redirect_stdout(io.StringIO()):
            adb.main(["--reads", fq, "--database", P["release"],
                      "--index", P["index"], "--primers",
                      os.path.join(args.workdir, "primers.tsv"),
                      "--sample-id", sid, "--outdir", out,
                      "--minimap2", mm2, "--threads", str(args.threads),
                      "--floors", ""])
        os.remove(fq)
        truth = {t: cat for t, cat in taxa}
        with open(os.path.join(out, f"{sid}.clusters.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows.append({"sample": sid, **r})
        # which source owns each cluster: re-derived from read names is
        # not possible after counting, so match the consensus to the
        # closest source sequence instead
        for r in rows:
            if r["sample"] != sid:
                continue
            best = None
            for t in truth:
                d = al.bounded_edit(r["consensus"], seqs[t], 60)
                if best is None or d < best[0]:
                    best = (d, t)
            r["source"], r["category"] = (best[1], truth[best[1]]) \
                if best and best[0] <= 60 else ("", "spurious")
            r["edits_to_source"] = best[0] if best else -1
        if (k + 1) % 10 == 0:
            print(f"  {k + 1}/{args.samples} samples", file=sys.stderr)
    with open(os.path.join(args.workdir, "clusters.json"), "w") as fh:
        json.dump(rows, fh)
    report(rows, tax, args)
    return 0


def outcome(r, tax, floors) -> str:
    if r["category"] == "spurious":
        return "spurious"
    lca = al.RANKS.index(r["lca_rank"]) if r["lca_rank"] else -1
    depth = adb.apply_floors(lca, float(r["identity"]), floors)
    lim = al.RANKS.index(ALLOWED[r["category"]])
    t, h = tax[r["source"]], tax.get(r["best_hit"], {})
    if any(t[x] != h.get(x, "") for x in al.RANKS[:min(depth, lim) + 1]) \
            or depth > lim:
        return "wrong"
    return "deepest" if depth == lim else "higher"


def report(rows, tax, args) -> None:
    grids = []
    for sp in (0.99, 0.993, 0.995, 0.997):
        for gn in (0.97, 0.975, 0.98, 0.985, 0.99):
            grids.append({"species": sp, "genus": gn, "family": 0.95})
    lines = ["species\tgenus\tcategory\tn\tdeepest\thigher\twrong"]
    for fl in grids:
        for cat in ("in_db", "species_out", "genus_out"):
            sel = [r for r in rows if r["category"] == cat]
            c = Counter(outcome(r, tax, fl) for r in sel)
            n = max(1, len(sel))
            lines.append(f"{fl['species']}\t{fl['genus']}\t{cat}\t{len(sel)}\t"
                         f"{c['deepest'] / n:.3f}\t{c['higher'] / n:.3f}\t"
                         f"{c['wrong'] / n:.3f}")
    text = "\n".join(lines)
    with open(os.path.join(args.workdir, "floors.tsv"), "w") as fh:
        fh.write(text + "\n")
    print(text)
    print(f"spurious clusters: {sum(r['category'] == 'spurious' for r in rows)}"
          f" of {len(rows)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--release", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--samples", type=int, default=150)
    ap.add_argument("--held-genera", type=int, default=100)
    ap.add_argument("--held-species", type=int, default=150)
    ap.add_argument("--error-model", default="nanopore_hac",
                    choices=sorted(sim.ERROR_MODELS))
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--report-only", action="store_true",
                    help="re-read clusters.json and print the grid")
    args = ap.parse_args(argv)
    if args.report_only:
        with open(os.path.join(args.workdir, "clusters.json")) as fh:
            rows = json.load(fh)
        tax = al.load_taxonomy(os.path.join(args.release, "taxonomy.tsv"))
        report(rows, tax, args)
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
