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


adb = _load("assign_db", "bin/assign_db.py")
adp = _load("assign_db_pairs", "bin/assign_db_pairs.py")
pf = _load("preflight", "bin/preflight.py")
sdb = _load("synthesise_db", "bin/synthesise_db.py")
import ampseqlib as al          # noqa: E402
import protocols  # noqa: E402
import protocols.nemits2 as its2  # noqa: E402

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
    # Drawn last, so the records above are as they always were.
    recs.append(("Chabertiidae", "Oesophagostomum",
                 "Oesophagostomum venulosum", rand(rng, 240)))
    # An (AT)6 inside R1's window, for reads that carry (AT)14.
    recs.append(("Chabertiidae", "Oesophagostomum", "Oesophagostomum asperum",
                 "G" + rand(rng, 39) + "AT" * 6 + "G" + rand(rng, 187)))
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
    # The release's lineage column: a line within a species.
    lines = {"Ostertagia ostertagi": "goat line"}
    w("taxonomy.tsv", ["id"] + al.NCBI_RANKS + [
        "deepest_rank", "species_status", "members", "grade", "coverage",
        "sources", "lineage"],
      [[r["id"]] + LINEAGE + [r["family"], r["genus"], r["species"],
                              "species" if r["species"] else "genus",
                              r["status"], len(r["members"]), "complete",
                              "1.000", "fixture", lines.get(r["species"], "")]
       for r in recs])
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
    # Reads 3 edits (98.75%) from their species' one record, no sister
    # near: the gap rule stops at genus, the tiered rule names it,
    # flagged tentative.
    "Oesophagostomum venulosum": ("species", "Oesophagostomum venulosum"),
    # Reads whose (AT)n is 8 units (16 bases) longer than the record's: 16
    # edits base by base (93%), one by copy number - named, not held at
    # family.
    "Oesophagostomum asperum": ("species", "Oesophagostomum asperum"),
}
# The sources whose reads differ from their record, and by how much.
VARIANT = {"Oesophagostomum venulosum": 3}
EXPANDED = {"Oesophagostomum asperum": ("AT" * 6, "AT" * 14)}


def source_seqs(recs) -> dict:
    out = {}
    for r in recs:
        sp = r["species"]
        if sp not in SOURCES or r["status"] == "unnamed":
            continue
        out[sp] = (mutate(random.Random(11), r["seq"], VARIANT[sp])
                   if sp in VARIANT else r["seq"].replace(*EXPANDED[sp], 1)
                   if sp in EXPANDED else r["seq"])
    return out


def write_pairs(path1, path2, seqs, L, rng, n=40, in_28s=False):
    """2xL pairs of each source, behind a 0-3 base spacer, as the
    nemabiome protocol builds them; a read past its insert runs into the
    adapter. in_28s: every other read 2 starts inside the 28S stretch,
    0-7 bases before the motif, with no NC2 (half of PRJNA1524548)."""
    with gzip.open(path1, "wt") as f1, gzip.open(path2, "wt") as f2:
        k = 0
        for its2 in seqs.values():
            insert = its2 + STRETCH
            for _ in range(n):
                s1 = rand(rng, rng.randint(0, 3)) + NC1 + insert + \
                    al.revcomp(NC2) + ADAPTER_R1
                s2 = rand(rng, rng.randint(0, 3)) + NC2 + al.revcomp(insert) + \
                    al.revcomp(NC1) + ADAPTER_R2
                if in_28s and k % 2:
                    s2 = (al.revcomp(insert) + al.revcomp(NC1) + ADAPTER_R2)[
                        rng.randint(0, 7):]
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
              protocols.scores_over_reference(man) and not protocols.scores_over_reference(
                  {"insert": "between the primers"}))
        check("the calling rule: tiered for ITS2, floors for V4, explicit wins",
              protocols.calling_rule(man) == "tiered"
              and protocols.calling_rule({}) == "floors"
              and protocols.calling_rule(man, "floors") == "floors"
              and protocols.calling_rule(man, "gap") == "gap"
              and protocols.calling_rule({"calling_rule": "gap"}) == "gap"
              and protocols.calling_rule({"calling_rule": "tiered"}) == "tiered")
        check("a tiered calling rule, like gap, means span scoring",
              protocols.protocol_for({}, None, "tiered") == "nemits2"
              and protocols.protocol_for({"calling_rule": "tiered"}) == "nemits2")
        check("R2's 28S stretch is taken from the geometries (150-20-83 = 47), "
              "less the 7 bases before the 28S motif it is cut at",
              its2.PROTOCOL.tail(man, 150) == 40 and its2.PROTOCOL.tail(man, 251) == 40
              and its2.PROTOCOL.tail({}, 150) == 0, str(its2.PROTOCOL.tail(man, 150)))
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
    g = its2.GAP_FLOORS
    c = its2.span_call([(0, 240, "full"), (0, 150, "part"), (8, 240, "sis")],
                      tax, "gap", g)
    check("a partial matching over a shorter span does not tie with a "
          "full-length match", c["ties"] == ["full"], str(c["ties"]))
    c = its2.span_call([(0, 240, "full"), (6, 240, "sis")], tax, "gap", g)
    check("species when the next species is >= 4 edits further off",
          c["rank"] == "species" and c["gap"] == 6, str(c))
    c = its2.span_call([(0, 240, "full"), (2, 240, "sis")], tax, "gap", g,
                      gap_exact=0)
    check("genus when it is closer, with how close",
          c["rank"] == "genus" and its2.nearest_next(c, tax) ==
          ("Haemonchus contortus 100.0%", "Haemonchus placei +2"),
          str(its2.nearest_next(c, tax)))
    c = its2.span_call([(0, 240, "full"), (2, 240, "sis")], tax, "gap", g)
    check("an exact match is species with the next species 2 edits off "
          "(GAP_EXACT, as its2_classify.py)", c["rank"] == "species"
          and c["gap"] == 2, str(c["rank"]))
    c = its2.span_call([(0, 240, "full"), (1, 240, "sis")], tax, "gap", g)
    check("... but not 1 edit off", c["rank"] == "genus")
    c = its2.span_call([(1, 240, "full"), (4, 240, "sis")], tax, "gap", g)
    check("one edit from its record, the gap must still be 4",
          c["rank"] == "genus" and c["gap"] == 3, str(c))
    c = its2.span_call([(0, 240, "full"), (1, 240, "hyb")], tax, "gap", g)
    check("an unnamed record is not another species",
          c["rank"] == "species" and c["next"] is None, str(c))
    c = its2.span_call([(0, 240, "full"), (1, 240, "shr")], tax, "gap", g)
    check("a shared record is another species, named by its members",
          c["rank"] == "genus" and its2.nearest_next(c, tax)[1] ==
          "Haemonchus shared: contortus / placei +1",
          str(its2.nearest_next(c, tax)))
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "taxonomy.tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(["id"] + al.NCBI_RANKS + [
                "species_status", "contained_in", "contains"]) + "\n")
            for rid, cin, cs in (("a", "b", ""), ("b", "", "a;c"),
                                 ("c", "b", ""), ("d", "", "")):
                fh.write("\t".join([rid] + LINEAGE + [
                    "Trichostrongylidae", "Haemonchus", f"Haemonchus {rid}",
                    "named", cin, cs]) + "\n")
        got = al.load_taxonomy(path).contained
    check("contained_in and contains load as one list per record",
          got == {"a": ["b"], "b": ["a", "c"], "c": ["b"]}, str(got))
    tax.contained.update({"full": ["sis"], "sis": ["part"]})
    c = its2.span_call([(0, 240, "full"), (6, 240, "sis")], tax, "gap", g)
    check("a record the release flags as contained in another species' "
          "record is 0 edits from it, however the read's ends differ",
          c["rank"] == "genus" and c["gap"] == 0
          and its2.nearest_next(c, tax)[1] == "Haemonchus placei +0",
          str(c))
    c = its2.span_call([(0, 240, "sis"), (6, 240, "part")], tax, "gap", g)
    check("... not when the flagged record is of the same species",
          c["rank"] == "species", str(c))
    tax.contained.update({"full": ["gone"]})
    del tax.contained["sis"]
    c = its2.span_call([(0, 240, "full"), (6, 240, "sis")], tax, "gap", g)
    check("... nor when it is not in the database (held out)",
          c["rank"] == "species" and c["gap"] == 6, str(c))
    tax.contained.clear()
    c = its2.span_call([(3, 240, "full"), (20, 240, "sis")], tax, "gap", g)
    check("species needs 99% identity too (3 edits in 240 is 98.75%)",
          c["rank"] == "genus", str(c["rank"]))
    c = its2.span_call([(3, 240, "full"), (20, 240, "sis")], tax, "tiered", g)
    cols = its2.span_columns(c, tax)
    check("tiered: the same read is a tentative species at >= 98%, the gap "
          "rule's genus call beside it",
          c["rank"] == "species" and c["tentative"]
          and cols["tentative"] == "yes" and cols["gap_rule_rank"] == "genus"
          and cols["gap_rule_taxon"] == "Haemonchus", str(cols))
    check("tiered: a tentative species is labelled so, never pooled with "
          "confident calls",
          adb.label_for(c, tax, al.NCBI_RANKS.index("family"))[1]
          == "Haemonchus contortus (tentative)")
    c = its2.span_call([(3, 240, "full"), (5, 240, "sis")], tax, "tiered", g)
    check("tiered: not when the next species is closer than the gap",
          c["rank"] == "genus" and not c["tentative"], str(c))
    c = its2.span_call([(5, 240, "full"), (30, 240, "sis")], tax, "tiered", g)
    check("tiered: not under 98% (5 edits in 240 is 97.9%)",
          c["rank"] == "genus" and not c["tentative"], str(c))
    c = its2.span_call([(0, 240, "full"), (6, 240, "sis")], tax, "tiered", g)
    cols = its2.span_columns(c, tax)
    check("tiered: a call the gap rule makes is not tentative, and the "
          "gap rule's call is the same",
          c["rank"] == "species" and not c["tentative"]
          and cols["tentative"] == "" and cols["gap_rule_rank"] == "species"
          and cols["gap_rule_taxon"] == "Haemonchus contortus", str(cols))
    c = its2.span_call([(3, 240, "full"), (20, 240, "sis")], tax, "tiered",
                       g, tentative_identity=0.99)
    check("tiered: --tentative-identity 0.99 is the gap rule",
          c["rank"] == "genus" and not c["tentative"])
    c = its2.span_call([(0, 240, "full"), (2, 240, "sis")], tax, "floors",
                      {"species": 0.995})
    check("the floors rule ignores the gap", c["rank"] == "species")
    c = its2.span_call([(0, 120, "full")], tax, "gap", g, min_species_span=150)
    check("a species call rests on at least --min-species-span bases",
          c["rank"] == "genus" and c["span"] == 120, str(c))
    c = its2.span_call([(0, 90, "full")], tax, "gap", g)
    check("a span under MIN_SPAN is no evidence", c["best_hit"] == "")
    hits = [{"t": f"r{i}", "AS": 100 - (i // 2), "ident": 0.99}
            for i in range(10)]
    got = its2.span_candidates(hits, 3)
    check("candidates keep every record tied with the last one taken",
          got == ["r0", "r1", "r2", "r3"], str(got))
    hits = [{"t": "full", "match": 238, "alen": 240, "ts": 0, "tlen": 240,
             "te": 240, "AS": 470},
            {"t": "part", "match": 150, "alen": 150, "ts": 0, "tlen": 150,
             "te": 150, "AS": 300},
            {"t": "clipped", "match": 200, "alen": 200, "ts": 0, "tlen": 240,
             "te": 200, "AS": 400}]
    got = its2.span_candidates(hits, 2)
    check("candidates are ranked by identity over the record's span, so a "
          "partial is not pushed out by score alone",
          got == ["part", "full"], str(got))


def test_repeats() -> None:
    print("\nrepeat-aware scoring")
    c, w = its2.collapse("GCATATATATGC")
    check("a repeat is two copies of its unit and a marker for how many",
          c == "GCATAT" + chr(0x104) + "GC" and w == [1] * 6 + [4, 1, 1],
          repr((c, w)))
    check("a homopolymer of >= 4 bases is a repeat, of 3 is not",
          its2.collapse("GAAAAG")[0] == "GAA" + chr(0x104) + "G"
          and its2.collapse("GAAAG")[0] == "GAAAG")
    check("a repeat cut off at a window's end has no marker",
          its2.collapse("GCATATAT", open_end="end")[0] == "GCATAT")
    rng = random.Random(4)
    head, tail_ = "G" + rand(rng, 59), "G" + rand(rng, 167)
    rec = head + "AT" * 6 + tail_
    longer = head + "AT" * 14 + tail_
    sc = its2.span_scores([(longer + STRETCH, "whole")], ["r"], {"r": rec},
                          repeats=True)
    check("whole: an (AT)n 8 units longer costs one edit, not 16",
          sc[0][:2] == (1, len(rec)) and sc[0][3] == 16, str(sc))
    check("... and base by base without repeats",
          its2.span_scores([(longer + STRETCH, "whole")], ["r"], {"r": rec})[0][:2]
          == (16, len(rec)))
    rec2 = rand(rng, 120) + "TA" * 6 + "C" + rand(rng, 60)
    win = (rec2[-83:].replace("TA" * 6, "TA" * 12)) + STRETCH
    sc = its2.span_scores([(rec2[:130], "start"), (win, "end")], ["r"],
                          {"r": rec2}, tail=47, repeats=True)
    check("a pair: an expansion in R2's window costs one edit, its 28S "
          "stretch none, as in the release's scorer",
          sc[0][0] == 1 and sc[0][3] > 10, str(sc))
    tax = al.Taxonomy(al.NCBI_RANKS)
    lin = LINEAGE + ["Trichuridae", "Trichuris"]
    for rid, sp in (("a", "Trichuris trichiura"), ("b", "Trichuris incognita")):
        tax[rid] = dict(zip(al.NCBI_RANKS, lin + [sp]),
                        species_status="named", grade="complete")
    recs = {"a": rec, "b": head + "AT" * 9 + tail_}
    c = its2.span_call(its2.span_scores([(rec + STRETCH, "whole")], ["a", "b"],
                                        recs, repeats=True),
                       tax, "gap", its2.GAP_FLOORS)
    check("sister species that differ only in a repeat's length keep their "
          "margin: it is measured plain",
          c["rank"] == "species" and c["gap"] == 6, str(c))


def test_informal() -> None:
    print("\ninformal \"Genus sp.\" records, after the named ones")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "taxonomy.tsv")
        with open(path, "w") as fh:
            fh.write("\t".join(["id"] + al.NCBI_RANKS + ["species_status"]) + "\n")
            for rid, sp, st in (("n", "Necator americanus", "named"),
                                ("i", "", "unnamed"), ("h", "", "unnamed")):
                fh.write("\t".join([rid] + LINEAGE + ["Ancylostomatidae",
                                                      "Necator", sp, st]) + "\n")
        with open(os.path.join(tmp, "members.tsv"), "w") as fh:
            fh.write("id\tsource\taccession\tlabel\trole\n")
            for rid, sp in (("n", "Necator americanus"), ("i", "Necator sp. type II"),
                            ("i", "Necator sp. 1 AB-2019"),
                            ("h", "Necator americanus x Necator sp.")):
                fh.write(f"{rid}\tx\tA1\t" + ";".join(
                    LINEAGE + ["Ancylostomatidae", "Necator", sp]) + "\tkept\n")
        tax = al.load_taxonomy(path)
    check("a record whose sources are all \"Genus sp. X\" is informal; a "
          "hybrid is not", tax.informal == {"i"}, str(tax.informal))
    g = its2.GAP_FLOORS
    c = its2.named_first_call([(1, 240, "n"), (0, 240, "i")], tax, "gap", g)
    check("a read 1 edit from a named species and identical to an informal "
          "record of its genus is that species: named records first",
          c["rank"] == "species" and c["best_hit"] == "n", str(c))
    c = its2.span_call([(1, 240, "n"), (0, 240, "i")], tax, "gap", g)
    check("... scored alongside, the informal record would take it to genus",
          c["rank"] == "genus" and c["best_hit"] == "i", str(c))
    tax["n"]["release_lineage"] = "type I"
    tax["i"]["release_lineage"] = "Necator sp. type II"
    c = its2.span_call([(0, 240, "i")], tax, "gap", g)
    check("a call carries the lineage its nearest records share",
          c["release_lineage"] == "Necator sp. type II"
          and adb.label_for(c, tax, al.NCBI_RANKS.index("family"))[1]
          == "Necator (Necator sp. type II)", str(c["release_lineage"]))
    c = its2.span_call([(0, 240, "i"), (0, 240, "n")], tax, "gap", g)
    check("... and none when the tied records' lineages differ",
          c["release_lineage"] == "", str(c["release_lineage"]))
    del tax["n"]["release_lineage"], tax["i"]["release_lineage"]
    c = its2.named_first_call([(40, 240, "n"), (2, 240, "i")], tax, "gap", g)
    check("a read no named record reaches is placed by an informal record",
          c["rank"] == "genus" and c["best_hit"] == "i", str(c))
    c = its2.named_first_call([(40, 240, "n")], tax, "gap", g)
    check("... and without one stops where the named records leave it",
          c["best_hit"] == "n" and c["rank"] not in ("species", "genus"),
          str(c))


def test_indel_fold() -> None:
    print("\nfolding indel variants")
    rng = random.Random(12)
    a = rand(rng, 280)
    ins = a[:100] + "T" + a[100:]
    dl2 = a[:50] + a[51:150] + a[151:]
    other = mutate(rng, a, 6)
    cents = [(a, "", 1000), (ins, "", 50), (ins + "G", "", 200), (dl2, "", 20),
             (other, "", 30)]
    out, n = its2.fold_indel_variants(cents, 2.0)
    got = {x: c for x, _y, c in out}
    check("a whole insert one indel from a sequence 8x as abundant is "
          "folded into it, two indels at 32x",
          got.get(a) == 1070 and n == 70, str(sorted(got.values())))
    check("... one not that much rarer stays, as does a sequence 6 edits off",
          got.get(ins + "G") == 200 and got.get(other) == 30, str(got.values()))
    w1, w2 = a[:120], a[-120:]
    i1 = (w1[:40] + "A" + w1[40:])[:120]
    out, n = its2.fold_indel_variants([(w1, w2, 400), (i1, w2, 40)], 2.0)
    check("a window pair whose R1 window carries an indel (its last base "
          "pushed out) is folded too", n == 40 and len(out) == 1, str(out))


def test_primers_and_merging() -> None:
    print("\nprimers behind a spacer, and merging mates")
    ps = al.PrimerSet(fwd=[("NC1", NC1)], rev=[("NC2", NC2)], max_edits=3)
    rng = random.Random(8)
    insert = rand(rng, 240) + STRETCH
    r1 = ("x", NC1 + insert[:200], "I" * 220)
    at_nc2 = ("x", "GA" + NC2 + al.revcomp(insert)[:200], "I" * 222)
    in_28s = ("x", al.revcomp(insert)[3:203], "I" * 200)
    f = Counter()
    a = its2.PROTOCOL.orient(r1, at_nc2, ps, f)
    b = its2.PROTOCOL.orient(r1, in_28s, ps, f)
    check("read 2 at NC2 and read 2 inside 28S are both cut at the 28S "
          "motif, so they start on the same base",
          a and b and a[1][0][:60] == b[1][0][:60]
          and a[1][0].startswith(its2.R2_MOTIF) and f["read2_in_28s"] == 1,
          f"{a and a[1][0][:25]} {b and b[1][0][:25]}")
    check("... and read 1 runs through what read 2 began with",
          a[3] == NC2 and b[3] == al.revcomp(insert)[3:7 + len(its2.R2_MOTIF)],
          f"{a[3]} {b[3]}")
    check("a read 2 with neither NC2 nor the motif is not placed",
          its2.PROTOCOL.orient(r1, ("x", rand(rng, 200), "I" * 200), ps, f)
          is None)
    ps = al.PrimerSet(fwd=[("NC1", NC1)], rev=[("NC2", NC2)], max_edits=3)
    rng = random.Random(3)
    body = rand(rng, 100)
    ends = [ps.match(rand(rng, k) + NC1 + body, "fwd") for k in range(4)]
    check("NC1 is found behind a 0-3 base spacer",
          [e[1] if e else None for e in ends] == [20, 21, 22, 23], str(ends))
    insert = rand(rng, 287)
    s1, s2 = insert[:228], al.revcomp(insert)[:228]
    m = its2.merge_pair(s1, "I" * 228, s2, "I" * 228)
    check("overlapping mates merge into the insert", m and m[0] == insert,
          str(len(m[0]) if m else None))
    s1e = s1[:200] + ("A" if s1[200] != "A" else "C") + s1[201:]
    m = its2.merge_pair(s1e, "I" * 228, s2, "I" * 228)
    check("an error in the overlap does not stop a merge; R1's base is kept",
          m is not None and m[0][200] == s1e[200])
    s1, s2 = insert[:127], al.revcomp(insert)[:127]
    check("mates that do not overlap are not merged",
          its2.merge_pair(s1, "I" * 127, s2, "I" * 127) is None)


def test_groups() -> None:
    print("\norganism groups over NCBI ranks")
    R = al.NCBI_RANKS
    g = protocols.organism_group
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
            write_pairs(r1, r2, seqs, L, random.Random(L), in_28s=True)
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
            check(f"2x{L}: read 2 inside 28S (half the pairs) is used, and "
                  "both layouts denoise into one sequence per source",
                  js["reads"].get("read2_in_28s") == 120
                  and len(rows) == len(SOURCES),
                  f"{js['reads'].get('read2_in_28s')} {len(rows)}")
            check(f"2x{L}: every pair is used, whatever its spacer",
                  js["reads"]["clustered"] == js["reads"]["reads_in"] == 240,
                  str(js["reads"]))
            if L == 250:
                check("2x250: overlapping mates are merged and scored whole",
                      js["reads"]["merged"] == 240, str(js["reads"]))
            else:
                check("2x150: windows of 127 + 120 read bases (R2 from the "
                      "28S motif, less the spacer allowance)",
                      js["settings"]["read_windows"] == [127, 120],
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
            ov = got.get(("species", "Oesophagostomum venulosum"), {})
            check(f"2x{L}: a species 98.75% from its one record is named, "
                  "flagged tentative, the gap rule's genus call beside it",
                  ov.get("tentative") == "yes"
                  and ov.get("label") == "Oesophagostomum venulosum (tentative)"
                  and ov.get("gap_rule_rank") == "genus"
                  and ov.get("gap_rule_taxon") == "Oesophagostomum",
                  str({k: ov.get(k) for k in ("label", "identity", "tentative",
                                              "gap_rule_rank", "gap")}))
            oo = got.get(("species", "Ostertagia ostertagi"), {})
            check(f"2x{L}: the release's lineage of the nearest records is "
                  "reported, and said in the label",
                  oo.get("release_lineage") == "goat line"
                  and oo.get("label") == "Ostertagia ostertagi (goat line)"
                  and hc.get("release_lineage") == "",
                  f"{oo.get('release_lineage')!r} {oo.get('label')!r}")
            check(f"2x{L}: confident calls are not flagged",
                  hc.get("tentative") == "" and hc.get("gap_rule_rank") == "species")
            gap_out = os.path.join(tmp, f"il{L}gap")
            with contextlib.redirect_stdout(io.StringIO()):
                adp.main(["--r1", r1, "--r2", r2, "--database", rel,
                          "--primers", prim, "--assay", ASSAY,
                          "--sample-id", "S", "--read-length", str(L),
                          "--primer-offset", "3", "--outdir", gap_out,
                          "--minimap2", mm2, "--calling-rule", "gap"])
            gr = calls_by_source(list(csv.DictReader(open(os.path.join(
                gap_out, "S.clusters.tsv")), delimiter="\t")))
            oa = got.get(("species", "Oesophagostomum asperum"), {})
            check(f"2x{L}: an (AT)n 16 bases longer than the record costs "
                  "one edit, and the species is named",
                  oa.get("edits") == "1" and oa.get("tentative") == "",
                  str({k: oa.get(k) for k in ("identity", "edits", "span")}))
            nr_out = os.path.join(tmp, f"il{L}norep")
            with contextlib.redirect_stdout(io.StringIO()):
                adp.main(["--r1", r1, "--r2", r2, "--database", rel,
                          "--primers", prim, "--assay", ASSAY,
                          "--sample-id", "S", "--read-length", str(L),
                          "--primer-offset", "3", "--outdir", nr_out,
                          "--minimap2", mm2, "--no-repeats"])
            nr = calls_by_source(list(csv.DictReader(open(os.path.join(
                nr_out, "S.clusters.tsv")), delimiter="\t")))
            check(f"2x{L}: --no-repeats charges it base by base, and the "
                  "species is not named",
                  ("species", "Oesophagostomum asperum") not in nr,
                  str(sorted(nr)))
            check(f"2x{L}: --calling-rule gap stops it at genus",
                  ("genus", "Oesophagostomum") in gr
                  and ("species", "Oesophagostomum venulosum") not in gr,
                  str(sorted(gr)))
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
              js["settings"].get("calling_rule") == "tiered"
              and js["settings"].get("score_over") == "reference span",
              str(js["settings"]))
        _ = man


def test_plate_misassignment() -> None:
    """Index misassignment on an ITS2 plate. The model works on calls,
    not sequences, so an ITS2 release changes nothing in it: a simulated
    combinatorial plate (simulate.simulate_plate_counts), every taxon a
    nematode with an NCBI lineage, through synthesis against the
    miniature release. Every taxon is one group, helminths, so the
    grouped fit is the pooled one. The rates come back, and a blank next
    to the heavy column is not read as contamination."""
    print("\nindex misassignment on an ITS2 plate")
    sim = _load("simulate", "bin/simulate.py")
    p = sim.simulate_plate_counts(random.Random(11), 8, 12, rate_row=4e-4,
                                  rate_col=7e-4)
    names = sorted({k for cs in p["counts"].values() for k in cs})
    lin = {k: ";".join(LINEAGE + ["Trichostrongylidae", f"Genus{i}",
                                  f"Genus{i} species{i}"])
           for i, k in enumerate(names)}
    cols = ["sample_id", "cluster", "reads", "rel_abundance", "rank",
            "taxon", "label", "identity", "coverage", "lineage", "best_hit",
            "best_hit_lineage", "near_hits", "lca_rank", "db_resolution",
            "db_margin", "origin", "merged", "length", "database_version",
            "consensus"]
    with tempfile.TemporaryDirectory() as tmp:
        rel = os.path.join(tmp, "rel")
        write_release(rel, build_records(random.Random(1)))
        ad, pfd = os.path.join(tmp, "assign_db"), os.path.join(tmp, "pf")
        os.makedirs(ad)
        os.makedirs(pfd)
        for sid, cs in p["counts"].items():
            with open(os.path.join(ad, f"{sid}.clusters.tsv"), "w") as fh:
                w = csv.DictWriter(fh, cols, delimiter="\t",
                                   lineterminator="\n")
                w.writeheader()
                for i, (k, n) in enumerate(sorted(cs.items())):
                    sp = lin[k].split(";")[-1]
                    w.writerow({"sample_id": sid, "cluster": f"{sid}_c{i}",
                                "reads": n, "rank": "species", "taxon": sp,
                                "label": sp, "identity": "1.0000",
                                "coverage": "1.000", "lineage": lin[k],
                                "best_hit": f"ITS2_{k}", "near_hits": 1,
                                "lca_rank": "species", "origin": "group",
                                "merged": 0, "length": 240,
                                "database_version": "nemits2-test",
                                "consensus": "ACGT"})
            n = sum(cs.values())
            json.dump({"sample_id": sid, "clusters": len(cs),
                       "reads": {"reads_in": n, "clustered": n},
                       "settings": {}},
                      open(os.path.join(ad, f"{sid}.assign_db.json"), "w"))
        with open(os.path.join(pfd, "file_resolution.tsv"), "w") as fh:
            fh.write("run_id\tsample_id\tsanitised_id\tplate_row\t"
                     "plate_col\tsample_type\n")
            for sid, t in p["types"].items():
                r, c = p["geom"][sid]
                fh.write(f"R1\t{sid}\t{sid}\t{r}\t{c}\t{t}\n")
        json.dump({}, open(os.path.join(pfd, "preflight_summary.json"), "w"))
        out = os.path.join(tmp, "synthesis")
        with contextlib.redirect_stdout(io.StringIO()):
            sdb.main(["--assign-dir", ad, "--preflight-dir", pfd,
                      "--database", rel, "--run-id", "R1", "--outdir", out,
                      "--min-depth", "50", "--permutations", "100"])
        qc = json.load(open(os.path.join(out, "db_run_qc.json")))
        calls = list(csv.DictReader(open(os.path.join(out, "db_calls.tsv")),
                                    delimiter="\t"))
    mis = qc["index_misassignment"]
    check("every taxon is a helminth over NCBI ranks",
          {c["group"] for c in calls} == {"helminths"},
          str({c["group"] for c in calls}))
    for ax, truth in (("row", 4e-4), ("col", 7e-4)):
        a = mis[ax]
        check(f"the {ax} rate comes back within 20% on an ITS2 release",
              a["basis"] == "fitted" and abs(a["rate"] / truth - 1) <= 0.2,
              f"{a['basis']} {a['rate']} {a['reason']}")
    blank = next(s for s, t in sorted(p["types"].items())
                 if t != "unknown" and p["counts"][s].get("worm"))
    worm = lin["worm"].split(";")[-1]
    cb = next(c for c in calls if c["sample_id"] == blank
              and c["taxon"] == worm)
    check("a blank next to the heavy column is misassignment, not "
          "contamination", cb["call"] in ("misassigned", "below_blanks")
          and float(cb["misassign_expected"]) >= 0.5 * int(cb["reads"]),
          str({k: cb[k] for k in ("reads", "call", "misassign_expected")}))
    check("and the worm's limit stays at the floor in field wells",
          all(int(c["limit_reads"]) <= 2 for c in calls
              if c["taxon"] == worm and c["sample_type"] == "unknown"),
          str(sorted({c["limit_reads"] for c in calls if c["taxon"] == worm})))


def main() -> int:
    test_loading()
    test_span_distances()
    test_span_call()
    test_repeats()
    test_informal()
    test_indel_fold()
    test_primers_and_merging()
    test_groups()
    test_preflight()
    test_end_to_end()
    test_plate_misassignment()
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
