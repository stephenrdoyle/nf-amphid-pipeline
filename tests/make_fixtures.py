#!/usr/bin/env python3
"""
Build test fixtures for pre-flight: a valid reference panel, samplesheet and
run manifest using the real AMPSEQ index sets and the real control tag.

Deterministic - a fixed seed - so a failing test is reproducible.
"""

from __future__ import annotations

import argparse
import csv
import os
import random

TAG = "TGGTCAATTTCAAGAGAT"

# The real primer overhangs, so simulated constructs match the assay.
OVERHANG_F = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
OVERHANG_R = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

# NOTE: every sequence here is SYNTHETIC. These are not real helminth
# sequences and must never be used as reference data - they exist only to
# give the pipeline something structurally realistic to be tested against.
#
# What "structurally realistic" means matters. Fully random sequences make
# competitive mapping trivially separable, so a panel built from them would
# pass tests that a real panel would fail. Instead:
#
#   - targets at a locus share conserved flanks, the way ITS2 sits between
#     conserved 5.8S and 28S, and universal primers anneal in those flanks
#   - species differ in the variable core, at a realistic divergence
#   - haplotypes of one species differ by only a handful of bases, which is
#     what makes taxon-level aggregation necessary rather than decorative
#
# The result is a panel where assignment is genuinely a discrimination
# problem, and where the score-margin metric has something to measure.

# species_key, locus, taxon, core_len, divergence_from_ancestor, corr_factor
SPECIES = [
    ("Nam", "ITS2", "Necator_americanus", 260, 0.18, "1.42"),
    ("Tri", "ITS2", "Trichuris_trichiura", 275, 0.22, "0.78"),
    ("Asc", "ITS2", "Ascaris_lumbricoides", 250, 0.20, ""),
    ("Nam", "cox1", "Necator_americanus", 560, 0.12, "1.10"),
]

# species_key, locus, n_haplotypes, bases differing between haplotypes
HAPLOTYPES = {
    ("Nam", "ITS2"): (2, 3),
    ("Tri", "ITS2"): (1, 0),
    ("Asc", "ITS2"): (1, 0),
    ("Nam", "cox1"): (1, 0),
}


def random_seq(rng: random.Random, n: int, gc: float = 0.45) -> str:
    out = []
    for _ in range(n):
        if rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


def mutate(rng: random.Random, seq: str, n_changes: int) -> str:
    """Substitute exactly n positions, never to the same base."""
    s = list(seq)
    for i in rng.sample(range(len(s)), min(n_changes, len(s))):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def diverge(rng: random.Random, seq: str, fraction: float) -> str:
    return mutate(rng, seq, int(len(seq) * fraction))


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


PRIMERS: dict[str, tuple[str, str]] = {}


def write_primers(path: str) -> None:
    """The primers the fixture references were built around.

    They have to be written out, not just used internally: the Illumina
    read-prep trims them, so a fixture without them can only be tested
    against the prep's own guesses.
    """
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["name", "direction", "sequence"])
        for locus, (fwd, rev) in sorted(PRIMERS.items()):
            w.writerow([f"Sim_{locus}_F1", "fwd", fwd])
            w.writerow([f"Sim_{locus}_R1", "rev", revcomp(rev)])


def build_reference(rng: random.Random) -> list[dict]:
    rows: list[dict] = []
    loci = sorted({s[1] for s in SPECIES})

    # Per locus: universal primers sitting in conserved flanks, and an
    # ancestral core that each species diverges from.
    primers: dict[str, tuple[str, str]] = {}
    flanks: dict[str, tuple[str, str]] = {}
    ancestor: dict[str, str] = {}
    for locus in loci:
        primers[locus] = (random_seq(rng, 20), random_seq(rng, 20))
        PRIMERS[locus] = primers[locus]
        flanks[locus] = (random_seq(rng, 30), random_seq(rng, 30))
        ancestor[locus] = random_seq(rng, 600)

    for key, locus, taxon, core_len, divergence, cf in SPECIES:
        fwd, rev = primers[locus]
        f5, f3 = flanks[locus]
        core = diverge(rng, ancestor[locus][:core_len], divergence)

        n_hap, hap_diff = HAPLOTYPES[(key, locus)]
        for h in range(1, n_hap + 1):
            hap_core = core if h == 1 else mutate(rng, core, hap_diff)
            rows.append(
                dict(
                    seq_id=f"{key}_{locus}_h{h}",
                    locus=locus,
                    type="target",
                    taxon=taxon,
                    parent_id="",
                    default_conc="",
                    conc_units="",
                    corr_factor=cf,
                    sequence=fwd + f5 + hap_core + f3 + rev,
                )
            )

    # One control per target, built as parent + tag inserted near the 3'
    # end - not at the midpoint. The tag has to be readable by a single
    # mate, or a pair too long to merge cannot tell the control from its
    # target at all. Real panels put it here for that reason.
    # This is what the tag derivation must recover exactly - including the
    # cases where the insertion point happens to share bases with the tag,
    # which is what makes left-alignment necessary.
    for t in [r for r in rows if r["type"] == "target"]:
        parent = t["sequence"]
        at = max(0, len(parent) - 30)
        rows.append(
            dict(
                seq_id=f"{t['seq_id']}_ctrl",
                locus=t["locus"],
                type="control",
                taxon=t["taxon"],
                parent_id=t["seq_id"],
                default_conc="",
                conc_units="",
                corr_factor="",
                sequence=parent[:at] + TAG + parent[at:],
            )
        )

    # Quants carry the same primer sites so they amplify in the same
    # reaction and compete for reads, but an unrelated core so they cannot
    # be confused with a target.
    fwd, rev = primers["ITS2"]
    f5, f3 = flanks["ITS2"]
    for name, conc in (("QUANT_hi", "10000"), ("QUANT_lo", "100")):
        rows.append(
            dict(
                seq_id=name,
                locus="ITS2",
                type="quant",
                taxon="",
                parent_id="",
                default_conc=conc,
                conc_units="copies_per_rxn",
                corr_factor="1.0",
                sequence=fwd + f5 + random_seq(rng, 240) + f3 + rev,
            )
        )

    rows.append(
        dict(
            seq_id="pUC19_bb",
            locus="",
            type="vector",
            taxon="",
            parent_id="",
            default_conc="",
            conc_units="",
            corr_factor="",
            sequence=random_seq(rng, 800, gc=0.5),
        )
    )
    return rows


REF_COLS = [
    "seq_id",
    "locus",
    "type",
    "taxon",
    "parent_id",
    "default_conc",
    "conc_units",
    "corr_factor",
    "sequence",
]


def write_reference(rows: list[dict], path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(REF_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in REF_COLS) + "\n")


def load_index_map(path: str) -> tuple[list[tuple[str, int]], list[tuple[str, str]]]:
    i7: list[tuple[str, int]] = []
    i5: list[tuple[str, str]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "role":
                continue
            role, _name, seq, pos = parts[0], parts[1], parts[2], parts[3]
            if role == "i7":
                i7.append((seq, int(pos)))
            elif role == "i5":
                i5.append((seq, pos))
    return i7, i5


def _levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def _edit_safe(pool: list, min_edit: int = 3) -> list:
    """Drop indices that sit within min_edit of one already kept.

    Order-preserving, so the plate layout that the rest of the fixture
    builds on is unchanged apart from the omissions.
    """
    kept: list = []
    for item in pool:
        if all(_levenshtein(item[0], k[0]) >= min_edit for k in kept):
            kept.append(item)
    return kept


def build_samplesheet(
    i7: list[tuple[str, int]],
    i5: list[tuple[str, str]],
    n_samples: int,
    leave_unused: int,
    run_id: str,
    platform: str,
) -> list[dict]:
    """Fill a 96-well quadrant: odd rows (A,C,E,...) and odd columns.

    The full pool is used, including the one pair two edits apart
    (N721 vs N723). That is the real plate design, and a fixture that
    quietly used a cleaner set than the bench does would hide the
    warning pre-flight is supposed to raise.
    """
    rows_q = [r for r in i5 if (ord(r[1]) - ord("A") + 1) % 2 == 1]
    cols_q = [c for c in i7 if c[1] % 2 == 1]

    wells = [(r, c) for r in rows_q for c in cols_q]
    usable = wells[: n_samples + leave_unused]

    # Scatter the unused wells rather than taking them from the end.
    # Removing a contiguous block can complete whole rows, shrinking the
    # index pool to exactly the assigned set - which leaves no sentinel at
    # all, because nothing can hop into an index that is not in the pool.
    # Evenly spaced omissions keep every row and column represented.
    if leave_unused > 0:
        step = len(usable) / (leave_unused + 1)
        unused = {usable[min(len(usable) - 1, int(step * (k + 1)))]
                  for k in range(leave_unused)}
    else:
        unused = set()
    chosen = [w for w in usable if w not in unused][:n_samples]

    out = []
    for i, ((i5_seq, row), (i7_seq, col)) in enumerate(chosen):
        if i == 0:
            stype, sid = "neg_extraction", "BLANK-01"
        elif i == 1:
            stype, sid = "neg_ntc", "NTC-01"
        elif i == 2:
            stype, sid = "mock", "MOCK-A"
        elif i == 3:
            stype, sid = "pos_control", "POS-01"
        else:
            stype, sid = "unknown", f"iCCH-{i-3:03d}-r1"
        out.append(
            dict(
                run_id=run_id,
                platform=platform,
                sample_id=sid,
                index_i7=i7_seq,
                index_i5=i5_seq,
                sample_type=stype,
            )
        )
    return out


SS_COLS = ["run_id", "platform", "sample_id", "index_i7", "index_i5", "sample_type"]


def write_samplesheet(rows: list[dict], path: str) -> None:
    """Delimiter from the extension, the same rule the readers use.

    Writing commas into a file called .tsv produces no error anywhere -
    the reader splits on tabs, gets one wide column, and reports every
    required field as missing. So the writer has to follow the same rule
    as the reader rather than hardcode one, which is also what lets the
    fixtures move to TSV while the dozens of temporary .csv samplesheets
    inside the test suite keep working unchanged. Those are not legacy:
    they are what keeps the CSV path exercised now that the fixtures no
    longer do it.
    """
    sep = "," if path.lower().endswith(".csv") else "\t"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(sep.join(SS_COLS) + "\n")
        for r in rows:
            fh.write(sep.join(str(r[c]) for c in SS_COLS) + "\n")


def write_metadata(rows: list[dict], path: str) -> None:
    """A grouping file keyed on sample_id, the harder of the two levels.

    Keyed on the LIBRARY deliberately, because that is the level a
    subject-keyed file cannot express and so the one a fixture has to
    cover. pcr1_cleanup is the real experiment it exists for: the same
    PCR1 product split into an uncleaned arm and an AMPure-cleaned one,
    which under a subject key would give both arms the same value and a
    comparison with nothing in it.

    The negatives are annotated too. An arm's negatives are how you
    know whether that arm stayed clean, and a fixture where only the
    real samples carry a group never exercises the control column.

    One sample is left out on purpose: a library with no metadata row
    must be bucketed under '(no value)', not dropped, and a fixture in
    which every library is annotated cannot show the difference.
    """
    arms = ["none", "ampure"]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("sample_id\tpcr1_cleanup\n")
        for i, r in enumerate(rows[:-1]):
            fh.write(f"{r['sample_id']}\t{arms[i % 2]}\n")


def write_runs(path: str, run_id: str, platform: str, run_path: str) -> None:
    cols = [
        "run_id",
        "platform",
        "path",
        "barcode_repeats",
        "dna_input_ng",
        "lot_pol",
        "lot_primer",
        "lot_extract",
        "lot_plasmid",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        reps = "2" if platform == "nanopore" else ""
        fh.write(
            "\t".join(
                [run_id, platform, run_path, reps, "10", "P-2291", "IDX-0417", "QK-8830", "PL-0032"]
            )
            + "\n"
        )


def write_fake_fastqs(samples: list[dict], run_path: str) -> None:
    """Empty gzip-named placeholders: pre-flight resolves paths, it does not
    read reads, so their content is irrelevant here.

    Filenames are sanitised the way a sequencing centre would, which also
    keeps the fixture writable when a test deliberately uses a sample_id
    containing path characters.
    """
    import re

    os.makedirs(run_path, exist_ok=True)
    for i, s in enumerate(samples, 1):
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", s["sample_id"])
        for read in ("R1", "R2"):
            name = f"{safe}_S{i}_L001_{read}_001.fastq.gz"
            open(os.path.join(run_path, name), "wb").close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="tests/data")
    ap.add_argument("--index-map", default="assets/index_map.tsv")
    ap.add_argument("--n-samples", type=int, default=88)
    ap.add_argument("--leave-unused", type=int, default=8)
    ap.add_argument("--platform", default="illumina")
    ap.add_argument("--seed", type=int, default=20260918)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    ref = build_reference(rng)
    write_reference(ref, os.path.join(args.outdir, "reference.tsv"))
    write_primers(os.path.join(args.outdir, "primers.tsv"))

    i7, i5 = load_index_map(args.index_map)
    run_path = os.path.join(args.outdir, "run01")
    samples = build_samplesheet(
        i7, i5, args.n_samples, args.leave_unused, "RUN01", args.platform
    )
    write_samplesheet(samples, os.path.join(args.outdir, "samplesheet.tsv"))
    write_metadata(samples, os.path.join(args.outdir, "metadata.tsv"))
    write_runs(os.path.join(args.outdir, "runs.tsv"), "RUN01", args.platform, run_path)
    if args.platform == "illumina":
        write_fake_fastqs(samples, run_path)
    else:
        os.makedirs(run_path, exist_ok=True)

    print(f"fixtures in {args.outdir}: {len(ref)} reference entries, {len(samples)} samples")


if __name__ == "__main__":
    main()
