#!/usr/bin/env python3
"""
Build a database-mode fixture: a miniature reference database release,
and a nanopore run whose samples hold taxa in, near and missing from it.

The real release (vespa-18s-refdb, ~97,000 records) lives outside this
repository, so the tests cannot depend on it. This writes a release of a
few dozen records in exactly the release format (docs/REFERENCE_DATABASE.md
section 4.1) - manifest with hashes, taxonomy over the nine PR2 ranks, the
look-alike mask in pre-flight's columns - from a synthetic tree whose
divergences are chosen to sit either side of the identity floors:

    species     0.5% from the genus sequence, so sisters differ by ~1%,
                which is what V4 gives for real sister species
    genus       1.5% from the family sequence
    family      and above, progressively more

and a run whose samples are built to ask the questions database mode has
to answer honestly:

  * a species held OUT of the release whose genus is in it - must stop
    at genus, not borrow a sister's name;
  * a genus held out whose family is in it - must stop at family;
  * a second haplotype of a species two bases from the first, at 30% of
    its reads - must come back as its own cluster, still that species;
  * a species whose insert carries a near-copy of the i5 anchor - its
    reads must survive the splitter and demultiplexer, which only the
    release's look-alike mask makes possible;
  * host, fungal and bacterial reads, for the off-target fractions;
  * a taxon at a handful of reads, below the size a cluster needs.

  python3 tests/make_db_fixture.py              # into tests/data/db

The truth is written beside the run (truth.tsv): per sample, each source
and what database mode should call it.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import make_fixtures as mf          # noqa: E402
import simulate as sim              # noqa: E402
import ampseqlib as al              # noqa: E402

ASSAY = os.path.join(ROOT, "assets", "assay", "vespa18s.json")
INDEX_MAP = os.path.join(ROOT, "assets", "index_map.tsv")

FWD = ("AMPSEQ_VESPA18S_29F", "AGCAGCCGCGGTAATTCC")
REV = ("AMPSEQ_VESPA18S_21b8_I_R_Nmod", "TCCGTCAATTYCTTNAASTTTC")
# The reverse primer as one template actually reads, for building reads.
REV_TEMPLATE = "TCCGTCAATTCCTTTAAGTTTC"

# (lineage over the nine ranks, in release) - the tree is built from these.
# Names follow PR2's style. A trailing '*' marks a taxon held OUT of the
# release: it exists in the tree, so it has a sequence, but no record.
LINEAGES = [
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Haemonchus;Haemonchus_contortus",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Haemonchus;Haemonchus_placei*",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Haemonchus;Haemonchus_similis",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Trichostrongylus;Trichostrongylus_colubriformis",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Trichostrongylus;Trichostrongylus_axei",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Teladorsagia*;Teladorsagia_circumcincta*",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Ostertagia;Ostertagia_ostertagi",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Trichostrongylidae;Ostertagia;Ostertagia_leptospicularis",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Ascarididae;Ascaris;Ascaris_suum",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Ascarididae;Ascaris;Ascaris_lumbricoides",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Ascarididae;Toxocara;Toxocara_canis",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Ascarididae;Toxocara;Toxocara_cati",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Strongyloididae;Strongyloides;Strongyloides_papillosus",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Chromadorea;Strongyloididae;Strongyloides;Strongyloides_stercoralis",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Enoplea;Trichuridae;Trichuris;Trichuris_trichiura",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Nematoda;Enoplea;Trichuridae;Trichuris;Trichuris_ovis",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Platyhelminthes;Trematoda;Fasciolidae;Fasciola;Fasciola_hepatica",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Platyhelminthes;Trematoda;Fasciolidae;Fasciola;Fasciola_gigantica",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Chordata;Mammalia;Bovidae;Bos;Bos_taurus",
    "Eukaryota;Obazoa;Opisthokonta;Metazoa;Chordata;Mammalia;Bovidae;Ovis;Ovis_aries",
    "Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Saccharomycetes;Saccharomycetaceae;Saccharomyces;Saccharomyces_cerevisiae",
    "Eukaryota;TSAR;Alveolata;Apicomplexa;Conoidasida;Eimeriorina;Eimeriidae;Eimeria;Eimeria_tenella",
    "Eukaryota;TSAR;Alveolata;Apicomplexa;Conoidasida;Eimeriorina;Eimeriidae;Eimeria;Eimeria_acervulina",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacterales_X;Enterobacteriaceae;Enterobacteriaceae_X;Escherichia;Escherichia_coli",
]

# Fraction of positions changed between a node and its parent, by the
# rank of the node. Species at 0.5% puts sisters ~1% apart over the
# insert, the V4 resolution limit the floors were calibrated around.
DIVERGENCE = {"domain": 0.20, "supergroup": 0.12, "division": 0.08,
              "subdivision": 0.06, "class": 0.05, "order": 0.04,
              "family": 0.035, "genus": 0.015, "species": 0.005}

# The species carrying an anchor near-copy, and the one with a hidden
# second haplotype.
LOOKALIKE_SPECIES = "Trichuris_trichiura"
HAPLOTYPE_SPECIES = "Fasciola_hepatica"

# Samples: (sample_id, sample_type, {source: reads}). A source is a
# species name, a species name + '#hap2' for its second haplotype.
SAMPLES = [
    ("BLANK-01", "neg_extraction", {"Saccharomyces_cerevisiae": 12}),
    ("NTC-01", "neg_ntc", {"Escherichia_coli": 3}),
    ("S01", "unknown", {"Haemonchus_contortus": 180,
                        "Trichuris_trichiura": 80, "Bos_taurus": 40}),
    ("S02", "unknown", {"Haemonchus_placei": 150, "Ascaris_suum": 120}),
    ("S03", "unknown", {"Teladorsagia_circumcincta": 120,
                        "Strongyloides_papillosus": 120,
                        "Eimeria_tenella": 60}),
    ("S04", "unknown", {"Fasciola_hepatica": 140,
                        "Fasciola_hepatica#hap2": 60,
                        "Escherichia_coli": 50}),
    ("S05", "unknown", {"Toxocara_canis": 200, "Ovis_aries": 50,
                        "Saccharomyces_cerevisiae": 40}),
    ("S06", "unknown", {"Haemonchus_contortus": 100,
                        "Trichostrongylus_colubriformis": 100,
                        "Ostertagia_ostertagi": 60, "Ascaris_suum": 3}),
]


def _id(seq: str) -> str:
    return "V4_" + hashlib.sha1(seq.encode()).hexdigest()[:12]


def _indels(rng: random.Random, seq: str, n: int) -> str:
    s = list(seq)
    for _ in range(n):
        j = rng.randrange(20, len(s) - 20)
        if rng.random() < 0.5:
            del s[j]
        else:
            s.insert(j, rng.choice("ACGT"))
    return "".join(s)


def build_tree(rng: random.Random) -> dict:
    """Species name -> {'seq', 'lineage', 'held_out'}; lengths vary with
    indels from order level up, as V4 lengths do between lineages."""
    node_seq: dict[tuple, str] = {(): mf.random_seq(rng, 540, gc=0.48)}
    species = {}
    for line in LINEAGES:
        names = line.split(";")
        held = any(n.endswith("*") for n in names)
        names = [n.rstrip("*") for n in names]
        for depth in range(1, len(names) + 1):
            key = tuple(names[:depth])
            if key in node_seq:
                continue
            rank = al.RANKS[depth - 1]
            parent = node_seq[key[:-1]]
            child = mf.diverge(rng, parent, DIVERGENCE[rank])
            if rank in ("supergroup", "division", "class", "order"):
                child = _indels(rng, child, rng.randint(1, 4))
            node_seq[key] = child
        seq = node_seq[tuple(names)]
        sp = names[-1]
        if sp == LOOKALIKE_SPECIES:
            # A near-copy of the i5 anchor (two edits) in the middle of the
            # insert: the T. incognita situation, in a database record.
            a5 = "TCGTCGGCAGCGTC"
            site = a5[:4] + "A" + a5[5:10] + "T" + a5[11:]
            mid = len(seq) // 2
            seq = seq[:mid] + site + seq[mid + len(site):]
        species[sp] = {"seq": seq, "lineage": names, "held_out": held}
    hap = species[HAPLOTYPE_SPECIES]["seq"]
    species[HAPLOTYPE_SPECIES + "#hap2"] = {
        "seq": mf.mutate(random.Random(7), hap, 2),
        "lineage": species[HAPLOTYPE_SPECIES]["lineage"],
        "held_out": True}       # a haplotype, not its own record
    return species


def _resolution(records: list, tax: dict) -> dict:
    """Deepest rank at which each record is distinguishable: no record
    within two edits names something different at that rank."""
    out = {}
    for rid, seq in records:
        rank = "species"
        for other, oseq in records:
            if other == rid or abs(len(seq) - len(oseq)) > 2:
                continue
            if al.bounded_edit(seq, oseq, 2) <= 2:
                for r in reversed(al.RANKS):
                    if tax[rid][r] == tax[other][r]:
                        cand = r
                        break
                else:
                    cand = "domain"
                if al.RANKS.index(cand) < al.RANKS.index(rank):
                    rank = cand
        out[rid] = rank
    return out


def write_release(outdir: str, species: dict, version: str = "fixture-0.1.0") -> dict:
    """The miniature release, in the agreed format. Returns its manifest."""
    os.makedirs(outdir, exist_ok=True)
    records, tax = [], {}
    for sp, d in species.items():
        if d["held_out"]:
            continue
        rid = _id(d["seq"])
        records.append((rid, d["seq"]))
        tax[rid] = dict(zip(al.RANKS, d["lineage"]))
    records.sort()
    res = _resolution(records, tax)

    def w(name, header, rows):
        with open(os.path.join(outdir, name), "w", encoding="utf-8") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(outdir, "refdb.fasta"), "w", encoding="utf-8") as fh:
        for rid, s in records:
            fh.write(f">{rid}\n{s}\n")
    w("taxonomy.tsv", ["id"] + al.RANKS + ["deepest_rank", "members", "sources"],
      [[rid] + [tax[rid][r] for r in al.RANKS] + ["species", 1, "fixture"]
       for rid, _s in records])
    w("members.tsv", ["id", "source", "accession", "label", "role"],
      [[rid, "fixture", f"FX{k:05d}.1", ";".join(tax[rid][r] for r in al.RANKS),
        "kept"] for k, (rid, _s) in enumerate(records)])
    for name in ("resolution_nanopore.tsv", "resolution_illumina_150x150.tsv"):
        w(name, ["id", "deepest_rank", "name"],
          [[rid, res[rid], tax[rid][res[rid]]] for rid, _s in records])
    lk = al.find_anchor_lookalikes(
        records, {"i5_inner_anchor": "TCGTCGGCAGCGTC",
                  "i7_inner_anchor": "CCGAGCCCACGAGAC"},
        al.LOOKALIKE_MAX_EDITS)
    w("anchor_lookalikes.tsv", al.LOOKALIKE_COLUMNS,
      [[r[c] for c in al.LOOKALIKE_COLUMNS] for r in lk])
    w("conflicts.tsv", ["id", "decision", "members", "labels", "kept_label"], [])
    w("primer_fit.tsv", ["id", "fit", "f_edits", "r_edits"],
      [[rid, "exact", 0, 0] for rid, _s in records])
    files = {}
    for name in sorted(os.listdir(outdir)):
        if name.endswith((".tsv", ".fasta")):
            with open(os.path.join(outdir, name), "rb") as fh:
                files[name] = hashlib.sha256(fh.read()).hexdigest()
    man = {
        "name": f"vespa18s_{version}", "version": version,
        "built": "2026-10-28T00:00:00Z", "code_commit": "fixture",
        "marker": "18S rRNA V4",
        "primers": {"forward": {"name": FWD[0], "sequence": FWD[1]},
                    "reverse": {"name": REV[0], "sequence": REV[1]},
                    "max_edits_in_silico": 2},
        "insert": "between the primers; primers not included",
        "anchors": {"i5_inner_anchor": "TCGTCGGCAGCGTC",
                    "i7_inner_anchor": "CCGAGCCCACGAGAC"},
        "lookalike_max_edits": al.LOOKALIKE_MAX_EDITS,
        "counts": {"records": len(records)},
        "files": files,
    }
    with open(os.path.join(outdir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=1)
    return man


def amplicon(insert: str) -> str:
    """Primer to primer, as PCR1 leaves it."""
    return FWD[1] + insert + al.revcomp(REV_TEMPLATE)


def write_primers(path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("name\tdirection\tsequence\n")
        fh.write(f"{FWD[0]}\tfwd\t{FWD[1]}\n{REV[0]}\trev\t{REV[1]}\n")


def expected_call(species: dict, source: str) -> tuple[str, str]:
    """(rank, name) database mode should reach for a source."""
    d = species[source]
    lin = d["lineage"]
    sp = source.split("#")[0]
    if source.endswith("#hap2") or not d["held_out"]:
        return "species", sp
    # held out: the genus if any congener is in the release, else family
    genus_in = any(not v["held_out"] and v["lineage"][7] == lin[7]
                   for v in species.values())
    return ("genus", lin[7]) if genus_in else ("family", lin[6])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--outdir", default=os.path.join(ROOT, "tests/data/db"))
    ap.add_argument("--batches", type=int, default=4)
    ap.add_argument("--chimera-fraction", type=float, default=0.3)
    ap.add_argument("--error-model", default="nanopore_hac",
                    choices=sorted(sim.ERROR_MODELS))
    ap.add_argument("--seed", type=int, default=20261028)
    args = ap.parse_args(argv)

    out = args.outdir
    rng = random.Random(args.seed)
    species = build_tree(rng)
    write_release(os.path.join(out, "release"), species)

    cfg = sim.DEFAULTS.copy() if hasattr(sim, "DEFAULTS") else {}
    cfg.update(al.load_assay(ASSAY))
    i7, i5 = mf.load_index_map(INDEX_MAP)
    sheet = mf.build_samplesheet(i7, i5, len(SAMPLES), 2, "DB01", "nanopore")
    for row, (sid, stype, _c) in zip(sheet, SAMPLES):
        row["sample_id"], row["sample_type"] = sid, stype
    mf.write_samplesheet(sheet, os.path.join(out, "samplesheet.tsv"))
    write_primers(os.path.join(out, "primers.tsv"))
    # RELATIVE path, as in the nanopore fixture: resolution of relative
    # paths is part of what the run exercises.
    mf.write_runs(os.path.join(out, "runs.tsv"), "DB01", "nanopore", "run")

    model = sim.ERROR_MODELS[args.error_model]
    reads, truth = [], []
    k = 0
    for row, (sid, _t, comp) in zip(sheet, SAMPLES):
        for src, n in comp.items():
            rank, name = expected_call(species, src)
            truth.append([sid, src, n, rank, name])
            mol0 = sim.build_molecule(amplicon(species[src]["seq"]),
                                      row["index_i7"], row["index_i5"], 2, cfg)
            for _ in range(n):
                m = mol0 if rng.random() < 0.5 else al.revcomp(mol0)
                m, _s, _i = sim.apply_errors(rng, m, model)
                reads.append([f"r{k:06d}", m])
                k += 1
    rng.shuffle(reads)
    # Ligation joins, as the splitter will meet them.
    n_join = int(len(reads) * args.chimera_fraction) // 2 * 2
    pool = [[f"J{j}", reads[j][1] + reads[j + 1][1]]
            for j in range(0, n_join, 2)] + reads[n_join:]
    rng.shuffle(pool)

    fq = os.path.join(out, "run", "fastq_pass")
    os.makedirs(fq, exist_ok=True)
    for f in os.listdir(fq):
        os.remove(os.path.join(fq, f))
    per = (len(pool) + args.batches - 1) // args.batches
    variant = args.error_model.replace("nanopore_", "")
    rg = f"{sim.uuid_like(rng)}_dna_r10.4.1_e8.2_400bps_{variant}@v5.2.0"
    for b in range(args.batches):
        part = pool[b * per:(b + 1) * per]
        with gzip.open(os.path.join(fq, f"FAX00002_pass_{b:04d}.fastq.gz"),
                       "wt") as fh:
            for rid, s in part:
                q = sim.qual_string(rng, len(s), model["rate"])
                fh.write(f"@{rid}\tRG:Z:{rg}\tPU:Z:SIM00000\tLB:Z:DB01\n"
                         f"{s}\n+\n{q}\n")
    with open(os.path.join(out, "truth.tsv"), "w", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["sample_id", "source", "reads", "expected_rank",
                     "expected_name"])
        wr.writerows(truth)
    # Every source sequence, for tests that build reads of their own.
    with open(os.path.join(out, "sources.fasta"), "w", encoding="utf-8") as fh:
        for sp, d in species.items():
            fh.write(f">{sp}\n{d['seq']}\n")
    print(f"database fixture in {out}/")
    print(f"  release: {sum(1 for d in species.values() if not d['held_out'])} "
          f"records; {len(SAMPLES)} samples, {len(pool):,} reads "
          f"({n_join // 2} concatemers) in {args.batches} batch files")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(BrokenPipeError):
        sys.exit(main())
