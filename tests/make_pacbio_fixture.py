#!/usr/bin/env python3
"""
Build the PacBio fixture: test_db's release and samples, sequenced as one
tagged PacBio library of HiFi reads.

The release is tests/data/db/release (make_db_fixture.py, committed), so
the right answer for every sample is test_db's (truth.tsv). Each read is

    [pad] [forward tag] forward primer ... insert ... rc(reverse primer) [rc(reverse tag)] [pad]

in either orientation, with HiFi errors (simulate.ERROR_MODELS
"pacbio_hifi"), and the pads 0-3 random bases, so the tag is not always
at the read's very end. The tags are i5 (forward) and i7 (reverse)
indices from the index map, two rows by five columns with two wells left
empty, so pre-flight places every sample on the plate and two tag pairs
are no sample's.

Beside the samples' reads, the library holds what a real SMRT cell does,
so every tier of the demultiplexer is exercised:

  * reads carrying another library's tags (index sequences no sample on
    this run uses): neither tag is this run's;
  * reads whose forward and reverse tags are two different samples'
    (a pair no sample has): chimeric;
  * reads cut short inside one tag: one end only;
  * primer dimers and off-target products with no primers.

Two SMRT cells: one written as FASTQ, the other as an unaligned BAM
(BGZF, as PacBio delivers hifi_reads.bam), so the run directory is read
as one library across both formats.

  python3 tests/make_pacbio_fixture.py          # into tests/data/pacbio
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import os
import random
import struct
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import make_db_fixture as dbf       # noqa: E402
import make_fixtures as mf          # noqa: E402
import simulate as sim              # noqa: E402
import ampseqlib as al              # noqa: E402

RUN_ID = "DB03"
MOVIES = ("m84001_261115_000000_s1", "m84001_261115_000000_s2")


# --------------------------------------------------------------------------
# an unaligned BAM, standard library only
# --------------------------------------------------------------------------

def _bgzf_block(data: bytes) -> bytes:
    """One BGZF block: a gzip member whose extra field holds its size."""
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = comp.compress(data) + comp.flush()
    bsize = 18 + len(body) + 8
    return (b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff"
            + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, bsize - 1)
            + body + struct.pack("<II", zlib.crc32(data) & 0xffffffff,
                                 len(data)))


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def write_bam(path: str, reads, sample: str = "") -> None:
    """reads: (name, sequence, quality) - written unaligned (flag 4), as
    PacBio's hifi_reads.bam, in BGZF blocks of up to 64 KB."""
    codes = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}
    text = f"@HD\tVN:1.6\tSO:unknown\tpb:5.0.0\n@RG\tID:{sample}\tPL:PACBIO\n"
    out = bytearray(b"BAM\x01" + struct.pack("<i", len(text)) + text.encode()
                    + struct.pack("<i", 0))
    for name, seq, qual in reads:
        nm = name.encode() + b"\x00"
        packed = bytearray()
        for k in range(0, len(seq), 2):
            hi = codes.get(seq[k], 15)
            lo = codes.get(seq[k + 1], 15) if k + 1 < len(seq) else 0
            packed.append(hi << 4 | lo)
        q = bytes(ord(c) - 33 for c in qual)
        body = (struct.pack("<iiBBHHHiiii", -1, -1, len(nm), 255, 4680, 0, 4,
                            len(seq), -1, -1, 0)
                + nm + bytes(packed) + q)
        out += struct.pack("<i", len(body)) + body
    with open(path, "wb") as fh:
        for k in range(0, len(out), 65280):
            fh.write(_bgzf_block(bytes(out[k:k + 65280])))
        fh.write(BGZF_EOF)


# --------------------------------------------------------------------------

def hifi_qual(rng: random.Random, n: int) -> str:
    """HiFi base qualities: mostly high, capped at Q93 ('~')."""
    return "".join(chr(33 + min(93, max(10, int(rng.gauss(60, 15)))))
                   for _ in range(n))


def molecule(insert: str, fwd_tag: str, rev_tag: str,
             rng: random.Random) -> str:
    pad = lambda: "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 3)))  # noqa: E731
    return (pad() + fwd_tag + dbf.amplicon(insert) + al.revcomp(rev_tag)
            + pad())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--outdir", default=os.path.join(ROOT, "tests/data/pacbio"))
    ap.add_argument("--seed", type=int, default=20261115)
    args = ap.parse_args(argv)

    out = args.outdir
    rng = random.Random(args.seed)
    # test_db's tree, from its own seed: the committed release.
    species = dbf.build_tree(random.Random(20261028))
    committed = dict(al.read_fasta(
        os.path.join(ROOT, "tests/data/db/sources.fasta")))
    if any(committed.get(k) != v["seq"] for k, v in species.items()):
        raise SystemExit("the test_db tree no longer matches "
                         "tests/data/db/sources.fasta; rerun "
                         "tests/make_db_fixture.py first")

    # test_db's samples on two rows by five columns (A and C; 1, 3, 5,
    # 7, 9), filled row by row: eight wells used, so two forward x
    # reverse tag pairs are no sample's - room for chimeric reads to
    # land, and for synthesis to measure misassignment on them.
    i7, i5 = mf.load_index_map(dbf.INDEX_MAP)
    rows = [s for s, r in i5 if r in ("A", "C")]
    cols = [s for s, c in sorted(i7, key=lambda x: x[1]) if c in (1, 3, 5, 7, 9)]
    wells = [(a, b) for a in rows for b in cols]
    sheet = [dict(r, run_id=RUN_ID, platform="pacbio", index_i5=a, index_i7=b)
             for r, (a, b) in zip(
                 al.read_table(os.path.join(ROOT, "tests/data/db/samplesheet.tsv")),
                 wells)]
    os.makedirs(out, exist_ok=True)
    mf.write_samplesheet(sheet, os.path.join(out, "samplesheet.tsv"))
    mf.write_runs(os.path.join(out, "runs.tsv"), RUN_ID, "pacbio", "run")

    model = sim.ERROR_MODELS["pacbio_hifi"]
    reads, truth = [], []
    by_sample = {r["sample_id"]: r for r in sheet}
    for sid, _t, comp in dbf.SAMPLES:
        row = by_sample[sid]
        for src, n in comp.items():
            rank, name = dbf.expected_call(species, src)
            truth.append([sid, src, n, rank, name])
            for _ in range(n):
                reads.append(molecule(species[src]["seq"], row["index_i5"],
                                      row["index_i7"], rng))
    # Another library's tags: index sequences no sample here uses.
    used5 = {r["index_i5"] for r in sheet}
    used7 = {r["index_i7"] for r in sheet}
    foreign5 = [s for s, _p in i5 if s not in used5][:3]
    foreign7 = [s for s, _p in i7 if s not in used7][:3]
    srcs = sorted(k for k in species if not k.endswith("#hap2"))
    extra = {"no_tag": 0, "pair_not_declared": 0, "one_end_only": 0,
             "no_primers": 0}
    for k in range(60):
        reads.append(molecule(species[srcs[k % len(srcs)]]["seq"],
                              foreign5[k % 3], foreign7[k % 3], rng))
        extra["no_tag"] += 1
    # A pair no sample has: one sample's forward tag, another's reverse.
    pairs = {(r["index_i5"], r["index_i7"]) for r in sheet}
    swaps = [(a["index_i5"], b["index_i7"]) for a in sheet for b in sheet
             if (a["index_i5"], b["index_i7"]) not in pairs][:5]
    for k in range(40):
        a, b = swaps[k % len(swaps)]
        reads.append(molecule(species["Haemonchus_contortus"]["seq"], a, b,
                              rng))
        extra["pair_not_declared"] += 1
    # Cut short three bases into the reverse tag.
    for k in range(30):
        row = sheet[2 + k % 4]
        reads.append(row["index_i5"]
                     + dbf.amplicon(species["Ascaris_suum"]["seq"])
                     + al.revcomp(row["index_i7"])[:3])
        extra["one_end_only"] += 1
    # Dimers and off-target products.
    for k in range(25):
        junk = "".join(rng.choice("ACGT") for _ in range(rng.randint(60, 300)))
        reads.append(sheet[0]["index_i5"] + junk)
        extra["no_primers"] += 1

    lib = []
    for k, m in enumerate(reads):
        m = m if rng.random() < 0.5 else al.revcomp(m)
        m, _s, _i = sim.apply_errors(rng, m, model)
        lib.append([m, hifi_qual(rng, len(m))])
    rng.shuffle(lib)

    run = os.path.join(out, "run")
    os.makedirs(run, exist_ok=True)
    for f in os.listdir(run):
        os.remove(os.path.join(run, f))
    half = len(lib) // 2
    cells = (lib[:half], lib[half:])
    names = [[f"{MOVIES[c]}/{k}/ccs" for k in range(len(cell))]
             for c, cell in enumerate(cells)]
    with gzip.open(os.path.join(run, f"{MOVIES[0]}.hifi_reads.fastq.gz"),
                   "wt") as fh:
        for name, (s, q) in zip(names[0], cells[0]):
            fh.write(f"@{name}\n{s}\n+\n{q}\n")
    write_bam(os.path.join(run, f"{MOVIES[1]}.hifi_reads.bam"),
              ((n, s, q) for n, (s, q) in zip(names[1], cells[1])),
              MOVIES[1])

    with open(os.path.join(out, "truth.tsv"), "w", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["sample_id", "source", "reads", "expected_rank",
                     "expected_name"])
        wr.writerows(truth)
        for k, v in extra.items():
            wr.writerow(["", f"(not a sample: {k})", v, "", ""])
    print(f"PacBio fixture in {out}/")
    print(f"  {len(dbf.SAMPLES)} samples, {len(lib):,} reads in two SMRT "
          f"cells (FASTQ and BAM); {sum(extra.values())} of them no "
          "sample's")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(BrokenPipeError):
        sys.exit(main())
