#!/usr/bin/env python3
"""
Tests for nemabiome ITS2 releases (nemabiome-its2-refdb, nemits2):
releases scored over their records' span, under the gap rule.

An ITS2 release differs from a vespa V4 one where it matters to
assignment (issue #9): its taxonomy has NCBI's eight ranks and says
whether a record with no species names two (shared) or none (unnamed);
its records stop where the 28S stretch before NC2 starts, and many are
partial, so a read always runs past its record (manifest score_over);
its Illumina windows count record bases; and its species calls follow
the gap rule - >= 99% identity and the nearest other species >= 4
edits further off - not identity floors.

The pieces first - loading, the span distances against a brute-force
alignment, ties and the gap, merging overlapping mates, primers behind
a spacer, organism groups over NCBI ranks, pre-flight's release checks -
then a miniature ITS2 release run end to end on Illumina at 2x150 and
2x250 and on nanopore, where each scenario has to come out as the gap
rule says:

  * a species 6 edits from its sister            - named to species;
  * a species 1 edit from its sister             - genus, "next +1";
  * a sequence held by a shared record           - genus, the record's
                                                    member species named;
  * a species 1 edit from an unnamed (hybrid)    - species: an unnamed
    record                                         record is no species.

The end-to-end tests need minimap2 (on the PATH or in ~/bin).

Run: python3 tests/test_nemits2.py
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
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
pf = _load("preflight", "bin/preflight.py")
sdb = _load("synthesise_db", "bin/synthesise_db.py")
import ampseqlib as al          # noqa: E402

PASS, FAIL = [], []

NC1 = "ACGTCTGGTTCAGGGTTGTT"
NC2 = "TTAGTTTCTTTTCCTCCGCT"
# The 28S stretch between the end of ITS2 and NC2: reads carry it,
# records do not.
STRETCH = "AACCTGAGCTCAGGCGTGATTACCCGCTGAACTTAAGCATATCATTT"
ADAPTER_R1 = al.revcomp("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG")
ADAPTER_R2 = al.revcomp("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
ASSAY = os.path.join(ROOT, "assets", "assay", "vespa18s.json")


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


def mutate(rng, seq, n, hi=100):
    """seq with n substitutions at distinct positions in 10..hi - inside
    R1's window, so even 2x150 pairs read every difference."""
    s = list(seq)
    for p in rng.sample(range(10, hi), n):
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


# --------------------------------------------------------------------------
# a miniature ITS2 release
# --------------------------------------------------------------------------

LINEAGE = ["Eukaryota", "Metazoa", "Nematoda", "Chromadorea", "Rhabditida"]


def build_records(rng) -> list[dict]:
    """The release's records: id, sequence, family, genus, species,
    status, and for a shared record the species its members name."""
    hc = rand(rng, 240)
    co = rand(rng, 236)
    tr = rand(rng, 244)
    oo = rand(rng, 238)
    recs = [
        ("Trichostrongylidae", "Haemonchus", "Haemonchus contortus", hc),
        ("Trichostrongylidae", "Haemonchus", "Haemonchus placei", mutate(rng, hc, 6)),
        ("Cooperiidae", "Cooperia", "Cooperia oncophora", co),
        ("Cooperiidae", "Cooperia", "Cooperia surnabada", mutate(rng, co, 1)),
        ("Trichostrongylidae", "Trichostrongylus", "", tr),
        ("Trichostrongylidae", "Trichostrongylus", "Trichostrongylus vitrinus",
         mutate(rng, tr, 12)),
        ("Trichostrongylidae", "Ostertagia", "Ostertagia ostertagi", oo),
        ("Trichostrongylidae", "Ostertagia", "", mutate(rng, oo, 1)),
    ]
    out = []
    for k, (fam, gen, sp, seq) in enumerate(recs):
        status = "named" if sp else ("shared" if gen == "Trichostrongylus"
                                     else "unnamed")
        out.append({"id": f"ITS2_{k:012x}", "seq": seq, "family": fam,
                    "genus": gen, "species": sp, "status": status,
                    "members": ["Trichostrongylus axei",
                                "Trichostrongylus colubriformis"]
                    if status == "shared" else
                    ["Ostertagia ostertagi x Ostertagia leptospicularis"]
                    if status == "unnamed" else [sp]})
    return out


def write_release(outdir: str, recs: list[dict]) -> dict:
    os.makedirs(outdir, exist_ok=True)

    def w(name, header, rows):
        with open(os.path.join(outdir, name), "w", encoding="utf-8") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(outdir, "refdb.fasta"), "w", encoding="utf-8") as fh:
        for r in recs:
            fh.write(f">{r['id']}\n{r['seq']}\n")
    w("taxonomy.tsv", ["id"] + al.NCBI_RANKS + [
        "deepest_rank", "species_status", "members", "grade", "coverage",
        "sources"],
      [[r["id"]] + LINEAGE + [r["family"], r["genus"], r["species"],
                              "species" if r["species"] else "genus",
                              r["status"], len(r["members"]), "complete",
                              "1.000", "fixture"] for r in recs])
    w("members.tsv", ["id", "source", "accession", "label", "role", "grade",
                      "coverage", "trimmed_start", "trimmed_end", "note"],
      [[r["id"], "fixture", f"FX{k:05d}.1",
        ";".join(LINEAGE + [r["family"], r["genus"], m]), "kept",
        "complete", "1.000", 0, 0, ""]
       for k, r in enumerate(recs) for m in r["members"]])
    res = [[r["id"], "species" if r["species"] else "genus",
            r["species"] or r["genus"], "complete", ""] for r in recs]
    for name in ("resolution_nanopore.tsv", "resolution_illumina_150x150.tsv",
                 "resolution_illumina_250x250.tsv"):
        w(name, ["id", "deepest_rank", "name", "grade", "margin"], res)
    w("anchor_lookalikes.tsv", al.LOOKALIKE_COLUMNS, [])
    w("conflicts.tsv", ["id", "decision", "members", "labels", "kept_label"], [])
    w("primer_fit.tsv", ["id", "fit", "f_edits", "r_edits"],
      [[r["id"], "unknown", "", ""] for r in recs])
    files = {}
    for name in sorted(os.listdir(outdir)):
        with open(os.path.join(outdir, name), "rb") as fh:
            files[name] = hashlib.sha256(fh.read()).hexdigest()
    man = {
        "name": "nemits2_fixture", "version": "fixture-its2",
        "marker": "ITS2 (nematodes)",
        "primers": {"forward": {"name": "NC1", "sequence": NC1},
                    "reverse": {"name": "NC2", "sequence": NC2}},
        "insert": "from the 3' end of NC1 to the start of the 28S stretch "
                  "before NC2",
        "score_over": "the reference's own span (a read's 28S stretch "
                      "counts neither for nor against)",
        "illumina_geometries": [
            {"read_length": 150, "windows": [130, 83],
             "file": "resolution_illumina_150x150.tsv"},
            {"read_length": 250, "windows": [230, 183],
             "file": "resolution_illumina_250x250.tsv"}],
        "anchors": {"i5_inner_anchor": "TCGTCGGCAGCGTC",
                    "i7_inner_anchor": "CCGAGCCCACGAGAC"},
        "lookalike_max_edits": al.LOOKALIKE_MAX_EDITS,
        "counts": {"records": len(recs)},
        "files": files,
    }
    with open(os.path.join(outdir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=1)
    return man


def write_primers(path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"name\tdirection\tsequence\nNC1\tfwd\t{NC1}\nNC2\trev\t{NC2}\n")


# What each scenario's reads come from, and the call the gap rule gives.
SOURCES = {
    "Haemonchus contortus": ("species", "Haemonchus contortus"),
    "Cooperia oncophora": ("genus", "Cooperia"),
    "": ("genus", "Trichostrongylus"),        # the shared record
    "Ostertagia ostertagi": ("species", "Ostertagia ostertagi"),
}


def source_seqs(recs) -> dict:
    return {r["species"]: r["seq"] for r in recs
            if r["species"] in SOURCES and r["status"] != "unnamed"}


def write_pairs(path1, path2, seqs, L, rng, n=40):
    """2xL pairs of each source, behind a 0-3 base spacer, as the
    nemabiome protocol builds them; a read past its insert runs into the
    adapter."""
    with gzip.open(path1, "wt") as f1, gzip.open(path2, "wt") as f2:
        k = 0
        for its2 in seqs.values():
            insert = its2 + STRETCH
            for _ in range(n):
                s1 = rand(rng, rng.randint(0, 3)) + NC1 + insert + \
                    al.revcomp(NC2) + ADAPTER_R1
                s2 = rand(rng, rng.randint(0, 3)) + NC2 + al.revcomp(insert) + \
                    al.revcomp(NC1) + ADAPTER_R2
                for fh, s in ((f1, s1), (f2, s2)):
                    s = (s + "A" * L)[:L]
                    fh.write(f"@p{k}\n{s}\n+\n{'I' * L}\n")
                k += 1


def write_nanopore(path, seqs, rng, n=30):
    with gzip.open(path, "wt") as fh:
        k = 0
        for its2 in seqs.values():
            for _ in range(n):
                m = NC1 + its2 + STRETCH + al.revcomp(NC2)
                m = "".join(c if rng.random() > 0.01 else rng.choice("ACGT")
                            for c in m)
                if rng.random() < 0.5:
                    m = al.revcomp(m)
                fh.write(f"@r{k}\n{m}\n+\n{'5' * len(m)}\n")
                k += 1


# --------------------------------------------------------------------------
# the pieces
# --------------------------------------------------------------------------


def test_loading() -> None:
    print("\nloading an ITS2 release")
    rng = random.Random(1)
    recs = build_records(rng)
    with tempfile.TemporaryDirectory() as tmp:
        man = write_release(tmp, recs)
        tax = al.load_taxonomy(os.path.join(tmp, "taxonomy.tsv"), man)
        check("NCBI's eight ranks are read from the header",
              tax.ranks == al.NCBI_RANKS, str(tax.ranks))
        r = tax[recs[0]["id"]]
        check("species_status and grade are kept beside the ranks",
              r.get("species_status") == "named" and r.get("grade") == "complete",
              str(r))
        check("a shared record carries the species its members name",
              tax.shared.get(recs[4]["id"]) == ["Trichostrongylus axei",
                                                "Trichostrongylus colubriformis"],
              str(tax.shared))
        check("an unnamed record's members are not taken as species",
              recs[7]["id"] not in tax.shared)
        check("a manifest's own ranks win over the known sets",
              al.ranks_in(["id", "a", "b"], {"ranks": ["a", "b"]}) == ["a", "b"])
        check("a header with neither set gives None",
              al.ranks_in(["id", "domain", "genus"]) is None)
        check("the release is scored over its records' span; V4 is not",
              al.scores_over_reference(man) and not al.scores_over_reference(
                  {"insert": "between the primers"}))
        check("the calling rule: gap for ITS2, floors for V4, explicit wins",
              al.calling_rule(man) == "gap" and al.calling_rule({}) == "floors"
              and al.calling_rule(man, "floors") == "floors"
              and al.calling_rule({"calling_rule": "gap"}) == "gap")
        check("R2's 28S stretch is taken from the geometries (150-20-83 = 47)",
              al.illumina_tail(man, 150) == 47 and al.illumina_tail(man, 251) == 47
              and al.illumina_tail({}, 150) == 0, str(al.illumina_tail(man, 150)))
    vespa = os.path.join(HERE, "data", "db", "release", "taxonomy.tsv")
    tv = al.load_taxonomy(vespa)
    check("a vespa release still loads over PR2's nine ranks",
          tv.ranks == al.RANKS and not tv.shared)


def _dp(ref, read, whole=True):
    """Brute-force anchored_distance: ref start global, read prefix free,
    read end global with ref past it free, or (whole) all of ref with the
    read past it free; best identity, longer span on a tie."""
    m, n = len(ref), len(read)
    D = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        D[i][0] = i
        for j in range(1, n + 1):
            D[i][j] = min(D[i - 1][j - 1] + (ref[i - 1] != read[j - 1]),
                          D[i - 1][j] + 1, D[i][j - 1] + 1)
    col = [D[i][n] for i in range(m + 1)]
    if whole:
        col[m] = min(D[m])
    best = None
    for i in range(1, m + 1):
        if best is None or col[i] * best[1] <= best[0] * i:
            best = (col[i], i)
    return best


def test_span_distances() -> None:
    print("\nspan distances")
    rng = random.Random(7)
    bad = []
    for _ in range(300):
        a, b = rand(rng, rng.randint(1, 50)), rand(rng, rng.randint(1, 50))
        if al.anchored_distance(a, b) != _dp(a, b):
            bad.append((a, b))
        # infix: all of a inside b
        m, n = len(a), len(b)
        D = [[0] * (n + 1) for _ in range(m + 1)]
        for i in range(1, m + 1):
            D[i][0] = i
            for j in range(1, n + 1):
                D[i][j] = min(D[i - 1][j - 1] + (a[i - 1] != b[j - 1]),
                              D[i - 1][j] + 1, D[i][j - 1] + 1)
        if al.infix_distance(a, b) != min(D[m]):
            bad.append(("infix", a, b))
    check("anchored and infix distances agree with a full alignment",
          not bad, str(bad[:2]))
    ref = rand(rng, 200)
    r2 = ref[-100:] + STRETCH
    check("R2's 28S stretch past the record's end is free",
          al.anchored_distance(ref[::-1], r2[::-1], 50) == (0, 100))
    r1 = ref + STRETCH[:20]
    check("a window longer than its record reads past it for free",
          al.anchored_distance(ref, r1, 150) == (0, 200))
    check("a partial record that starts late leaves the read before it free",
          al.anchored_distance(ref[10:], ref[:120], 90) == (0, 110))
    check("a whole consensus: the record placed inside, 28S free",
          al.infix_distance(ref, ref + STRETCH) == 0)


def test_span_call() -> None:
    print("\nties, the gap and the calling rule")
    tax = al.Taxonomy(al.NCBI_RANKS)
    lin = LINEAGE + ["Trichostrongylidae"]
    for rid, gen, sp, status in (
            ("full", "Haemonchus", "Haemonchus contortus", "named"),
            ("part", "Haemonchus", "Haemonchus placei", "named"),
            ("sis", "Haemonchus", "Haemonchus placei", "named"),
            ("hyb", "Haemonchus", "", "unnamed"),
            ("shr", "Haemonchus", "", "shared")):
        tax[rid] = dict(zip(al.NCBI_RANKS, lin + [gen, sp]),
                        species_status=status, grade="complete")
    tax.shared["shr"] = ["Haemonchus contortus", "Haemonchus placei"]
    g = adb.GAP_FLOORS
    c = adb.span_call([(0, 240, "full"), (0, 150, "part"), (8, 240, "sis")],
                      tax, "gap", g)
    check("a partial matching over a shorter span does not tie with a "
          "full-length match", c["ties"] == ["full"], str(c["ties"]))
    c = adb.span_call([(0, 240, "full"), (6, 240, "sis")], tax, "gap", g)
    check("species when the next species is >= 4 edits further off",
          c["rank"] == "species" and c["gap"] == 6, str(c))
    c = adb.span_call([(0, 240, "full"), (2, 240, "sis")], tax, "gap", g)
    check("genus when it is closer, with how close",
          c["rank"] == "genus" and adb.nearest_next(c, tax) ==
          ("Haemonchus contortus 100.0%", "Haemonchus placei +2"),
          str(adb.nearest_next(c, tax)))
    c = adb.span_call([(0, 240, "full"), (1, 240, "hyb")], tax, "gap", g)
    check("an unnamed record is not another species",
          c["rank"] == "species" and c["next"] is None, str(c))
    c = adb.span_call([(0, 240, "full"), (1, 240, "shr")], tax, "gap", g)
    check("a shared record is another species, named by its members",
          c["rank"] == "genus" and adb.nearest_next(c, tax)[1] ==
          "Haemonchus shared: contortus / placei +1",
          str(adb.nearest_next(c, tax)))
    c = adb.span_call([(3, 240, "full"), (20, 240, "sis")], tax, "gap", g)
    check("species needs 99% identity too (3 edits in 240 is 98.75%)",
          c["rank"] == "genus", str(c["rank"]))
    c = adb.span_call([(0, 240, "full"), (2, 240, "sis")], tax, "floors",
                      {"species": 0.995})
    check("the floors rule ignores the gap", c["rank"] == "species")
    c = adb.span_call([(0, 120, "full")], tax, "gap", g, min_species_span=150)
    check("a species call rests on at least --min-species-span bases",
          c["rank"] == "genus" and c["span"] == 120, str(c))
    c = adb.span_call([(0, 90, "full")], tax, "gap", g)
    check("a span under MIN_SPAN is no evidence", c["best_hit"] == "")
    hits = [{"t": f"r{i}", "AS": 100 - (i // 2), "ident": 0.99}
            for i in range(10)]
    got = adb.span_candidates(hits, 3)
    check("candidates keep every record tied with the last one taken",
          got == ["r0", "r1", "r2", "r3"], str(got))
    hits = [{"t": "full", "match": 238, "alen": 240, "ts": 0, "tlen": 240,
             "te": 240, "AS": 470},
            {"t": "part", "match": 150, "alen": 150, "ts": 0, "tlen": 150,
             "te": 150, "AS": 300},
            {"t": "clipped", "match": 200, "alen": 200, "ts": 0, "tlen": 240,
             "te": 200, "AS": 400}]
    got = adb.span_candidates(hits, 2)
    check("candidates are ranked by identity over the record's span, so a "
          "partial is not pushed out by score alone",
          got == ["part", "full"], str(got))


def test_primers_and_merging() -> None:
    print("\nprimers behind a spacer, and merging mates")
    ps = al.PrimerSet(fwd=[("NC1", NC1)], rev=[("NC2", NC2)], max_edits=3)
    rng = random.Random(3)
    body = rand(rng, 100)
    ends = [ps.match(rand(rng, k) + NC1 + body, "fwd") for k in range(4)]
    check("NC1 is found behind a 0-3 base spacer",
          [e[1] if e else None for e in ends] == [20, 21, 22, 23], str(ends))
    insert = rand(rng, 287)
    s1, s2 = insert[:228], al.revcomp(insert)[:228]
    m = adp.merge_pair(s1, "I" * 228, s2, "I" * 228)
    check("overlapping mates merge into the insert", m and m[0] == insert,
          str(len(m[0]) if m else None))
    s1e = s1[:200] + ("A" if s1[200] != "A" else "C") + s1[201:]
    m = adp.merge_pair(s1e, "I" * 228, s2, "I" * 228)
    check("an error in the overlap does not stop a merge; R1's base is kept",
          m is not None and m[0][200] == s1e[200])
    s1, s2 = insert[:127], al.revcomp(insert)[:127]
    check("mates that do not overlap are not merged",
          adp.merge_pair(s1, "I" * 127, s2, "I" * 127) is None)


def test_groups() -> None:
    print("\norganism groups over NCBI ranks")
    R = al.NCBI_RANKS
    g = sdb.organism_group
    check("a nematode is a helminth",
          g(LINEAGE + ["Trichostrongylidae"], "Craniata", R) == "helminths")
    check("a vertebrate is host",
          g(["Eukaryota", "Metazoa", "Chordata", "Mammalia"], "Chordata", R)
          == "host")
    check("a fungus is a fungus",
          g(["Eukaryota", "Fungi", "Ascomycota"], "Craniata", R) == "fungi")
    check("a call stopped at superkingdom is unassigned",
          g(["Eukaryota"], "Craniata", R) == "unassigned")
    check("PR2 lineages are read as before",
          g(["Eukaryota", "Obazoa", "Opisthokonta", "Metazoa", "Nematoda"],
            "Craniata") == "helminths")


def test_preflight() -> None:
    print("\npre-flight on an ITS2 release")
    rng = random.Random(1)
    with tempfile.TemporaryDirectory() as tmp:
        rel = os.path.join(tmp, "rel")
        write_release(rel, build_records(rng))
        diag = pf.Diagnostics()
        info = pf.check_reference_database(
            rel, {"fwd": [("NC1", NC1)], "rev": [("NC2", NC2)]}, ASSAY,
            tmp, diag)
        check("the release passes pre-flight's release checks",
              diag.ok and info.get("ranks") == al.NCBI_RANKS,
              "; ".join(diag.errors))


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


def calls_by_source(rows) -> dict:
    """(rank, taxon) of the largest cluster - every source here has one
    sequence, so its largest cluster is its call."""
    return {(r["rank"], r["taxon"]): r for r in rows}


def test_end_to_end() -> None:
    print("\nthe miniature ITS2 release end to end")
    mm2 = minimap2()
    if not mm2:
        print("  SKIP: minimap2 not found - database mode cannot be tested")
        FAIL.append(("minimap2 present", "not on PATH or in ~/bin"))
        return
    rng = random.Random(20261110)
    recs = build_records(rng)
    seqs = source_seqs(recs)
    with tempfile.TemporaryDirectory() as tmp:
        rel = os.path.join(tmp, "rel")
        man = write_release(rel, recs)
        prim = os.path.join(tmp, "primers.tsv")
        write_primers(prim)
        want = set(SOURCES.values())
        for L in (150, 250):
            r1, r2 = (os.path.join(tmp, f"S{L}_R{k}.fastq.gz") for k in (1, 2))
            write_pairs(r1, r2, seqs, L, random.Random(L))
            out = os.path.join(tmp, f"il{L}")
            with contextlib.redirect_stdout(io.StringIO()):
                adp.main(["--r1", r1, "--r2", r2, "--database", rel,
                          "--primers", prim, "--assay", ASSAY,
                          "--sample-id", "S", "--read-length", str(L),
                          "--primer-offset", "3", "--outdir", out,
                          "--minimap2", mm2])
            rows = list(csv.DictReader(open(os.path.join(out, "S.clusters.tsv")),
                                       delimiter="\t"))
            js = json.load(open(os.path.join(out, "S.assign_db.json")))
            got = calls_by_source(rows)
            check(f"2x{L}: every source is called as the gap rule says",
                  set(got) == want, str(sorted(got)))
            check(f"2x{L}: every pair is used, whatever its spacer",
                  js["reads"]["clustered"] == js["reads"]["reads_in"] == 160,
                  str(js["reads"]))
            if L == 250:
                check("2x250: overlapping mates are merged and scored whole",
                      js["reads"]["merged"] == 160, str(js["reads"]))
            else:
                check("2x150: windows of 127 + 127 read bases (R2 with its "
                      "28S stretch, less the spacer allowance)",
                      js["settings"]["read_windows"] == [127, 127],
                      str(js["settings"].get("read_windows")))
            co = got.get(("genus", "Cooperia"), {})
            check(f"2x{L}: a genus call says how close it came",
                  co.get("nearest", "").startswith("Cooperia oncophora 100")
                  and co.get("next") == "Cooperia surnabada +1",
                  f"{co.get('nearest')} | {co.get('next')}")
            tr = got.get(("genus", "Trichostrongylus"), {})
            check(f"2x{L}: a shared record is named by its members",
                  tr.get("nearest", "").startswith(
                      "Trichostrongylus axei / Trichostrongylus colubriformis"),
                  str(tr.get("nearest")))
            hc = got.get(("species", "Haemonchus contortus"), {})
            check(f"2x{L}: grade and gap are reported with the call",
                  hc.get("grade") == "complete" and hc.get("gap") == "6",
                  f"{hc.get('grade')} {hc.get('gap')}")
            check(f"2x{L}: lineages run over NCBI's ranks",
                  hc.get("lineage", "").split(";") ==
                  LINEAGE + ["Trichostrongylidae", "Haemonchus",
                             "Haemonchus contortus"], str(hc.get("lineage")))
        fq = os.path.join(tmp, "ont.fastq.gz")
        write_nanopore(fq, seqs, random.Random(5))
        out = os.path.join(tmp, "ont")
        with contextlib.redirect_stdout(io.StringIO()):
            adb.main(["--reads", fq, "--database", rel, "--primers", prim,
                      "--sample-id", "N", "--outdir", out, "--minimap2", mm2])
        rows = list(csv.DictReader(open(os.path.join(out, "N.clusters.tsv")),
                                   delimiter="\t"))
        js = json.load(open(os.path.join(out, "N.assign_db.json")))
        check("nanopore: consensus sequences carrying the 28S stretch are "
              "called as the gap rule says",
              set(calls_by_source(rows)) == want,
              str(sorted(calls_by_source(rows))))
        check("nanopore: the summary records the rule and the span scoring",
              js["settings"].get("calling_rule") == "gap"
              and js["settings"].get("score_over") == "reference span",
              str(js["settings"]))
        _ = man


def main() -> int:
    test_loading()
    test_span_distances()
    test_span_call()
    test_primers_and_merging()
    test_groups()
    test_preflight()
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
