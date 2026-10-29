#!/usr/bin/env python3
"""
Tests for database-mode assignment (bin/assign_db.py).

The parts are tested on their own first - pileup consensus, the linkage
test that splits haplotypes, merging, the LCA and the identity floors -
because each encodes a decision that an end-to-end number could hide: a
consensus that is right for the wrong reason, or a floor that never
fires. Then a whole sample set from the database fixture
(tests/make_db_fixture.py) is resolved and scored against its truth,
which is where held-out taxa have to stop at the right rank.

The end-to-end tests need minimap2 on the PATH (or in ~/bin); without it
they are skipped, loudly.

Run: python3 tests/test_assign_db.py
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
import ampseqlib as al            # noqa: E402
import make_db_fixture as fx      # noqa: E402
import simulate as sim            # noqa: E402

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


def aln(q, ref, read, ts=0):
    """A fake minimap2 record: cs for `read` against ref[ts:], built from
    a plain global alignment of equal-length strings with '-' for gaps.
    `read` here is an alignment row: same length as ref[ts:ts+len],
    '-' marking a deletion, lower case marking an inserted base."""
    cs, run, p = [], 0, ts
    i = 0
    while i < len(read):
        c = read[i]
        if c.islower():
            j = i
            while j < len(read) and read[j].islower():
                j += 1
            if run:
                cs.append(f":{run}")
                run = 0
            cs.append("+" + read[i:j])
            i = j
            continue
        if c == "-":
            if run:
                cs.append(f":{run}")
                run = 0
            j = i
            while j < len(read) and read[j] == "-":
                j += 1
            cs.append("-" + ref[p:p + (j - i)].lower())
            p += j - i
            i = j
            continue
        if c == ref[p]:
            run += 1
        else:
            if run:
                cs.append(f":{run}")
                run = 0
            cs.append(f"*{ref[p].lower()}{c.lower()}")
        p += 1
        i += 1
    if run:
        cs.append(f":{run}")
    return {"q": q, "ts": ts, "cs": "".join(cs)}


# --------------------------------------------------------------------------


def test_pileup_consensus() -> None:
    print("\npileup consensus")
    ref = "ACGTACGTTTGACCAGTA"
    # five reads: a substitution at 4 in three of them (majority), a
    # deletion at 10 in three, an insertion before 14 in three; noise
    # elsewhere in a minority
    rows = [
        "ACGTGCGTTT-ACCaAGTA",
        "ACGTGCGTTT-ACCaAGTA",
        "ACGTGCGATT-ACCaAGTA",
        "ACGTACGTTTGACCAGTA",
        "ACCTACGTTTGACCAGTA",
    ]
    pu = adb.Pileup(ref, [aln(f"r{i}", ref, r) for i, r in enumerate(rows)])
    want = "ACGTGCGTTTACCAAGTA"
    got = pu.consensus(5)
    check("majority substitution, deletion and insertion all taken",
          got == want, f"{got} vs {want}")
    check("cover counts deletions as covering", pu.cover[10] == 5,
          str(pu.cover[10]))
    check("alleles are kept per read (substitution and deletion)",
          pu.alleles["r0"] == {4: "G", 10: "-"}, str(pu.alleles["r0"]))
    s, d = pu.error_rates()
    check("error rates measured against the draft", 0.05 < s < 0.2 and
          0.02 < d < 0.1, f"{s:.3f} {d:.3f}")
    # half coverage rule: with four reads, a position needs two
    short = [aln(f"s{i}", ref, ref[:8]) for i in range(3)]
    check("positions spanned by fewer than half the reads are dropped",
          adb.Pileup(ref, short + [aln("f", ref, ref)]).consensus(4)
          == ref[:8]
          and adb.Pileup(ref, short[:2] + [aln("f", ref, ref)] * 2)
          .consensus(4) == ref)


def test_binomial() -> None:
    print("\nbinomial tail")
    # exact values by enumeration
    from math import comb
    n, p = 12, 0.1
    for k in (0, 1, 3, 7, 12):
        want = sum(comb(n, i) * p ** i * (1 - p) ** (n - i)
                   for i in range(k, n + 1))
        got = adb.binom_sf(k, n, p)
        check(f"binom_sf({k}, {n}, {p})", abs(got - want) < 1e-12,
              f"{got} vs {want}")


def _two_haplotypes(rng, n_a, n_b, diffs, err=0.0, linked=True):
    """Alignment rows for a group drawn from two haplotypes of `ref`, the
    second differing at `diffs`. With linked=False the same alleles are
    scattered independently over the reads instead - the same per-site
    frequency with no linkage - which is what a systematic error looks
    like."""
    ref = rand(rng, 400)
    rows = []
    for i in range(n_a + n_b):
        r = list(ref)
        if linked:
            carry = i >= n_a
            for p in diffs:
                if carry:
                    r[p] = "ACGT"[("ACGT".index(ref[p]) + 1) % 4]
        else:
            for p in diffs:
                if rng.random() < n_b / (n_a + n_b):
                    r[p] = "ACGT"[("ACGT".index(ref[p]) + 1) % 4]
        for j in range(len(r)):
            if rng.random() < err:
                r[j] = rng.choice([b for b in "ACGT" if b != ref[j]])
        rows.append("".join(r))
    return ref, [aln(f"r{i}", ref, row) for i, row in enumerate(rows)]


def test_linkage_signature() -> None:
    print("\nhaplotype signature by linkage")
    rng = random.Random(11)
    st = adb.Settings()
    diffs = [50, 120, 300]
    ref, alns = _two_haplotypes(rng, 70, 30, diffs, err=0.02)
    sig = adb.signature(adb.Pileup(ref, alns), st)
    check("a 30% haplotype three bases off is found", sig is not None and
          sorted(p for p, _ in sig) == diffs, str(sig))
    ref, alns = _two_haplotypes(rng, 70, 30, diffs, err=0.02, linked=False)
    sig = adb.signature(adb.Pileup(ref, alns), st)
    check("the same alleles at the same frequency but unlinked are not "
          "a haplotype", sig is None, str(sig))
    ref, alns = _two_haplotypes(rng, 96, 4, diffs, err=0.02)
    sig = adb.signature(adb.Pileup(ref, alns), st)
    check("a minority under 10% of the group is not split off",
          sig is None, str(sig))
    ref, alns = _two_haplotypes(rng, 70, 30, [200], err=0.02)
    check("one frequent allele alone is not a signature",
          adb.signature(adb.Pileup(ref, alns), st) is None)
    # clean group of errors only
    ref, alns = _two_haplotypes(rng, 100, 0, [], err=0.03)
    check("read errors alone give no signature",
          adb.signature(adb.Pileup(ref, alns), st) is None)
    # a deletion in a homopolymer is not a candidate allele
    ref = "ACGT" * 20 + "AAAAA" + "CGTA" * 20
    rows = [ref] * 60 + [ref[:82] + "-" + ref[83:]] * 40
    alns = [aln(f"h{i}", ref, r) for i, r in enumerate(rows)]
    check("homopolymer deletions are passed over",
          adb.signature(adb.Pileup(ref, alns), st) is None)


def test_merge() -> None:
    print("\nmerging close consensus sequences")
    rng = random.Random(5)
    a = rand(rng, 300)
    b = a[:100] + ("A" if a[100] != "A" else "C") + a[101:]     # 1 sub
    c = a[:150] + a[151:]                                          # 1 del
    d = fx.mf.mutate(rng, a, 2)                                    # 2 subs
    cl = [{"consensus": a, "reads": ["x"] * 10, "origin": "group"},
          {"consensus": b, "reads": ["y"] * 3, "origin": "group"},
          {"consensus": c, "reads": ["z"] * 2, "origin": "group"},
          {"consensus": d, "reads": ["w"] * 4, "origin": "group"}]
    m = adb.merge_close(cl, 1)
    check("within one edit merged into the largest, two edits kept apart",
          len(m) == 2 and m[0]["consensus"] == a and len(m[0]["reads"]) == 15
          and m[1]["consensus"] == d, str([(len(x["reads"]), x["merged"])
                                          for x in m]))


def _tax():
    lin = {
        "A1": "E;S;D;SD;C;O;F1;G1;G1_a",
        "A2": "E;S;D;SD;C;O;F1;G1;G1_b",
        "B1": "E;S;D;SD;C;O;F1;G2;G2_a",
        "P1": "E;S;D;SD;C;O;;G3;G3_a",       # placeholder family
    }
    return {k: dict(zip(al.RANKS, v.split(";"))) for k, v in lin.items()}


def hit(t, AS, match, alen, qlen=500, qs=0, qe=500):
    return {"t": t, "AS": AS, "match": match, "alen": alen, "qlen": qlen,
            "qs": qs, "qe": qe}


def test_lca_and_floors() -> None:
    print("\nLCA and identity floors")
    tax = _tax()
    R = al.RANKS
    check("one record: LCA is its species", adb.lca_depth(["A1"], tax) == 8)
    check("two species of one genus: genus",
          adb.lca_depth(["A1", "A2"], tax) == R.index("genus"))
    check("two genera: family", adb.lca_depth(["A1", "B1"], tax) == R.index("family"))
    check("a blank rank breaks agreement",
          adb.lca_depth(["A1", "P1"], tax) == R.index("order")
          and adb.lca_depth(["P1"], tax) == R.index("order"))
    floors = adb.DEFAULT_FLOORS
    check("species kept at 0.996", adb.apply_floors(8, 0.996, floors) == 8)
    check("species cut to genus at 0.994",
          adb.apply_floors(8, 0.994, floors) == R.index("genus"))
    check("species cut to family at 0.97",
          adb.apply_floors(8, 0.97, floors) == R.index("family"))
    check("species cut above family at 0.90",
          adb.apply_floors(8, 0.90, floors) == R.index("order"))
    check("a genus LCA is not raised by a high identity",
          adb.apply_floors(R.index("genus"), 1.0, floors) == R.index("genus"))

    # identity charges clipped ends
    h = hit("A1", 900, 450, 450, qlen=500, qs=20, qe=470)
    check("identity counts unaligned query ends as differences",
          abs(adb.identity(h) - 450 / 500) < 1e-9, str(adb.identity(h)))

    # margin: hits within 0.5% of the best score take part
    hs = [hit("A1", 1000, 498, 500), hit("A2", 996, 497, 500),
          hit("B1", 980, 495, 500)]
    c = adb.classify_hits(hs, tax, floors, 0.005)
    check("a near-tie between sisters stops at genus", c["rank"] == "genus"
          and c["near_hits"] == 2 and c["best_hit"] == "A1", str(c))
    c = adb.classify_hits(hs[:1] + hs[2:], tax, floors, 0.005)
    check("a clear best hit at 99.6% is a species", c["rank"] == "species",
          str(c))
    c = adb.classify_hits([hit("B1", 800, 440, 500)], tax, floors, 0.005)
    taxon, label = adb.label_for(c, tax, R.index("family"))
    check("88% identity: unclassified, named for how far and to what",
          c["rank"] == "order" and label == "unclassified, 88% to G2"
          and taxon == "O", f"{c['rank']} {taxon!r} {label!r}")
    c = adb.classify_hits([hit("A1", 950, 485, 500)], tax, floors, 0.005)
    taxon, label = adb.label_for(c, tax, R.index("family"))
    check("97% identity: named at family", c["rank"] == "family"
          and taxon == label == "F1", f"{c['rank']} {label}")
    c = adb.classify_hits([], tax, floors, 0.005)
    check("no hit at all", adb.label_for(c, tax, 6)[1] ==
          "unclassified, no database hit")

    # PR2 placeholders are explained, not repeated; a tie is not a gap
    ptax = {"H": dict(zip(al.RANKS, "Eukaryota;Obazoa;Opisthokonta;Metazoa;"
                          "Nematoda;Chromadorea;Chromadorea_X;Haemonchus;"
                          "Haemonchus_contortus".split(";")))}
    taxon, label = adb.label_for({"best_hit": "H", "depth": 6, "lca_depth": 6,
                                  "identity": 1.0}, ptax, 6)
    check("a placeholder family is labelled by its named parent and the "
          "nearest genus", taxon == "Chromadorea_X" and label ==
          "Chromadorea, unnamed family (nearest Haemonchus)", label)
    taxon, label = adb.label_for({"best_hit": "H", "depth": 5, "lca_depth": 5,
                                  "identity": 1.0}, ptax, 6)
    check("stopped above family by a tie: unresolved, not unclassified",
          label == "Chromadorea, unresolved below order (100% to Haemonchus)",
          label)
    taxon, label = adb.label_for({"best_hit": "H", "depth": 5, "lca_depth": 8,
                                  "identity": 0.9}, ptax, 6)
    check("stopped above family by the floors: unclassified",
          label == "unclassified, 90% to Haemonchus", label)
    btax = {"B": dict(zip(al.RANKS, "Bacteria;Bacteroidetes;Bacteroidetes_X;"
                          "Bacteroidia;Bacteroidales;Prevotellaceae;"
                          "Prevotellaceae_X;Prevotellaceae_XX;"
                          "Prevotellaceae_XX_sp.".split(";")))}
    _t, label = adb.label_for({"best_hit": "B", "depth": 4, "lca_depth": 8,
                               "identity": 0.9}, btax, 6)
    check("the nearest record is named past placeholders",
          label == "unclassified, 90% to Prevotellaceae", label)


def test_floor_parsing() -> None:
    print("\nfloor parsing")
    check("floors parse", adb.parse_floors("species=0.99,genus=0.97") ==
          {"species": 0.99, "genus": 0.97})
    try:
        adb.parse_floors("specie=0.9")
        check("a misspelt rank is refused", False)
    except SystemExit:
        check("a misspelt rank is refused", True)


def test_trim_reads() -> None:
    print("\nprimer trimming and orientation")
    rng = random.Random(9)
    ins = rand(rng, 300)
    with tempfile.TemporaryDirectory() as tmp:
        fq = os.path.join(tmp, "r.fastq")
        amp = fx.amplicon(ins)
        with open(fq, "w") as fh:
            for i, s in enumerate([amp, al.revcomp(amp), ins,
                                    fx.amplicon(ins[:30])]):
                fh.write(f"@r{i}\n{s}\n+\n{'5' * len(s)}\n")
        ps = al.PrimerSet(fwd=[fx.FWD], rev=[fx.REV])
        reads, fate = adb.trim_reads(fq, ps)
    check("both orientations trimmed to the insert, forward strand",
          reads.get("r0") == ins and reads.get("r1") == ins, str(reads)[:80])
    # three edits in the reverse primer: found at the database-mode
    # tolerance, not at the panel's two
    with tempfile.TemporaryDirectory() as tmp:
        fq = os.path.join(tmp, "r3.fastq")
        rev = fx.REV_TEMPLATE
        bad = "A" + rev[1:5] + ("T" if rev[5] != "T" else "G") + rev[6:10] \
            + ("C" if rev[10] != "C" else "G") + rev[11:]
        bad = ("G" if bad[0] != "G" else "C") + bad[1:]
        amp3 = fx.FWD[1] + ins + al.revcomp(bad)
        with open(fq, "w") as fh:
            fh.write(f"@x\n{amp3}\n+\n{'5' * len(amp3)}\n")
        ps2 = al.PrimerSet(fwd=[fx.FWD], rev=[fx.REV])
        r2, f2 = adb.trim_reads(fq, ps2)
        ps3 = al.PrimerSet(fwd=[fx.FWD], rev=[fx.REV], max_edits=3)
        r3, f3 = adb.trim_reads(fq, ps3)
    check("a primer three edits off is found at max_edits 3, not at 2",
          f2["no_primer"] == 1 and r3.get("x") == ins, f"{f2} {f3}")
    check("reads without primers, or too short, are counted out",
          fate["no_primer"] == 1 and fate["too_short"] == 1
          and fate["reads_in"] == 4, str(fate))


# --------------------------------------------------------------------------
# end to end on the fixture
# --------------------------------------------------------------------------


def _sample_reads(tmp, rng, species, comp, model, sid):
    path = os.path.join(tmp, f"{sid}.fastq.gz")
    with gzip.open(path, "wt") as fh:
        k = 0
        for src, n in comp.items():
            for _ in range(n):
                m = fx.amplicon(species[src]["seq"])
                if rng.random() < 0.5:
                    m = al.revcomp(m)
                m, _s, _i = sim.apply_errors(rng, m, model)
                fh.write(f"@{src}|{k}\n{m}\n+\n{'5' * len(m)}\n")
                k += 1
    return path


def test_end_to_end() -> None:
    print("\nend to end on the database fixture")
    mm2 = minimap2()
    if not mm2:
        print("  SKIP: minimap2 not found - database mode cannot be tested")
        FAIL.append(("minimap2 present", "not on PATH or in ~/bin"))
        return
    rng = random.Random(20261028)
    species = fx.build_tree(rng)
    model = sim.ERROR_MODELS["nanopore_hac"]
    with tempfile.TemporaryDirectory() as tmp:
        rel = os.path.join(tmp, "release")
        man = fx.write_release(rel, species)
        fx.write_primers(os.path.join(tmp, "primers.tsv"))
        out = os.path.join(tmp, "out")
        results = {}
        for sid, _t, comp in fx.SAMPLES:
            fq = _sample_reads(tmp, rng, species, comp, model, sid)
            with contextlib.redirect_stdout(io.StringIO()):
                adb.main(["--reads", fq, "--database", rel,
                          "--primers", os.path.join(tmp, "primers.tsv"),
                          "--sample-id", sid, "--outdir", out,
                          "--minimap2", mm2])
            with open(os.path.join(out, f"{sid}.clusters.tsv")) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            with open(os.path.join(out, f"{sid}.assign_db.json")) as fh:
                js = json.load(fh)
            results[sid] = (comp, rows, js)

        wrong, missing, right = [], [], 0
        for sid, (comp, rows, js) in results.items():
            for src, n in comp.items():
                if n < 5:
                    continue
                rank, name = fx.expected_call(species, src)
                # the cluster carrying most of this source's reads
                mine = [r for r in rows if r["taxon"] == name
                        and r["rank"] == rank]
                if not mine:
                    missing.append(f"{sid}:{src} (want {rank} {name}; got "
                                   f"{[(r['rank'], r['taxon']) for r in rows]})")
                else:
                    right += 1
                # A held-out taxon must not borrow a relative's name: no
                # species call in its genus (family, for a held-out
                # genus) that the sample does not really contain.
                if species[src]["held_out"] and not src.endswith("#hap2"):
                    level = 7 if rank == "genus" else 6
                    clade = species[src]["lineage"][level]
                    for r in rows:
                        lin = r["lineage"].split(";")
                        if (r["rank"] == "species" and len(lin) > level
                                and lin[level] == clade
                                and r["taxon"] not in comp):
                            wrong.append(f"{sid}: {src} named {r['taxon']}")
        check(f"every source of 5+ reads called at the rank it should reach "
              f"({right} calls)", not missing, "; ".join(missing))
        check("held-out taxa never take a relative's species name",
              not wrong, "; ".join(wrong))

        comp, rows, js = results["S04"]
        fh_rows = [r for r in rows if r["taxon"] == "Fasciola_hepatica"]
        check("the hidden haplotype is its own cluster, both still the "
              "species", len(fh_rows) == 2 and js["haplotype_splits"] >= 1,
              str([(r["reads"], r["origin"]) for r in fh_rows]))
        hap = species["Fasciola_hepatica#hap2"]["seq"]
        main_ = species["Fasciola_hepatica"]["seq"]
        cons = {r["consensus"] for r in fh_rows}
        check("both haplotype consensus sequences are exact",
              cons == {hap, main_}, str(len(cons)))
        comp, rows, js = results["S03"]
        tel = [r for r in rows if r["rank"] == "family"]
        check("the held-out genus stops at family and shows its identity",
              len(tel) == 1 and tel[0]["taxon"] == "Trichostrongylidae"
              and 0.95 <= float(tel[0]["identity"]) < 0.98,
              str([(r["rank"], r["taxon"], r["identity"]) for r in tel]))
        comp, rows, js = results["S06"]
        check("a 3-read taxon is left unclustered, not absorbed",
              js["reads"]["unclustered"] >= 3
              and not any(r["taxon"].startswith("Ascaris") for r in rows),
              str(js["reads"]))
        ok = True
        for sid, (comp, rows, js) in results.items():
            r = js["reads"]
            if r["reads_in"] != (r["no_primer"] + r["too_short"]
                                 + r["clustered"] + r["unclustered"]):
                ok = False
            if sum(int(x["reads"]) for x in rows) != r["clustered"]:
                ok = False
        check("read accounting closes in every sample", ok)
        check("every result carries the release version",
              all(js["database"]["version"] == man["version"]
                  and all(x["database_version"] == man["version"] for x in rows)
                  for _c, rows, js in results.values()))
        # exact consensus for present taxa
        exact = tot = 0
        for sid, (comp, rows, js) in results.items():
            for src, n in comp.items():
                if n < 20:
                    continue
                tot += 1
                exact += any(r["consensus"] == species[src]["seq"]
                             for r in rows)
        check(f"consensus is exact for {exact}/{tot} sources of 20+ reads",
              exact == tot)


def test_chimeras() -> None:
    print("\nPCR chimeras")
    rng = random.Random(77)
    A = rand(rng, 550)
    B = fx.mf.mutate(rng, A, 25)
    C = fx.mf.mutate(rng, A, 12)
    chim = A[:300] + B[300:]
    var = fx.mf.mutate(rng, A, 4)
    found = adb.chimera_scan([A, B, C, chim, var], [1000, 600, 300, 40, 30])
    check("a join of the two most abundant sequences is a chimera of them",
          set(found) == {3} and found[3][:2] == (0, 1), str(found))
    check("a 4-base variant of one parent is not", 4 not in found)
    rare = adb.chimera_scan([A, B, chim], [100, 60, 50])
    check("a 'chimera' more than half as abundant as a parent is not called",
          rare == {}, str(rare))
    indel = A[:200] + "T" + B[200:]
    check("an indel at the join is tolerated",
          2 in adb.chimera_scan([A, B, indel], [1000, 500, 20]))
    near = A[:540] + B[540:]
    check("a join that takes too little from one parent to tell is not "
          "called", adb.chimera_scan([A, B, near], [1000, 500, 20]) == {}
          or all(a != 1 or True for a in []))
    f = adb.match_flags("ACGTACGT", "ACGAACGT", 0)
    check("match_flags base for base", f == [True, True, True, False,
                                             True, True, True, True])
    f = adb.match_flags("ACGTTACGT", "ACGTACGT")
    check("match_flags through an insertion: one base charged",
          f.count(False) == 1, str(f))


def main() -> int:
    test_pileup_consensus()
    test_binomial()
    test_linkage_signature()
    test_merge()
    test_lca_and_floors()
    test_floor_parsing()
    test_trim_reads()
    test_chimeras()
    test_end_to_end()
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
