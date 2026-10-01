#!/usr/bin/env python3
"""
Build a nanopore fixture that the workflow can actually be run against.

The Illumina fixtures let `-profile test` exercise one branch of
main.nf. The nanopore branch had no fixture at all, which is how three
new processes - SPLIT_LIST, SPLIT_CONCATEMERS, MERGE_SPLIT - came to be
written, statically checked and never executed. A `-stub-run` reported
"Complete" while running exactly one process, because with no nanopore
run in the manifest and an Illumina-less file resolution, both branches
were empty and every task after pre-flight was skipped.

Two things this fixture must do that a single simulated FASTQ would
not:

  * arrive as SEVERAL batch files, the way MinKNOW writes them, because
    splitting now fans out one Nextflow task per input file and a
    one-file fixture would exercise a fan-out of one - which is the
    exact shape of the bug this pipeline has already had twice.
  * contain real concatemers, because a fixture of clean reads makes
    the splitter a no-op, and a no-op agrees with anything.

  python3 tests/make_ont_fixture.py --outdir tests/data/ont
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import io
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import importlib.util as _il      # noqa: E402
import make_fixtures as mf          # noqa: E402
import demux_ont as dx              # noqa: E402
import simulate as sim              # noqa: E402

ASSAY = os.path.join(ROOT, "assets", "assay", "ampseq_helminth.json")
INDEX_MAP = os.path.join(ROOT, "assets", "index_map.tsv")


def quiet(fn, *a):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--outdir", default=os.path.join(ROOT, "tests/data/ont"))
    ap.add_argument("--n-samples", type=int, default=8)
    ap.add_argument("--leave-unused", type=int, default=2,
                    help="wells left unassigned, so index misassignment "
                         "has a sentinel and is measurable at all")
    ap.add_argument("--reads-per-sample", type=int, default=60)
    ap.add_argument("--batches", type=int, default=4,
                    help="MinKNOW-style batch files. More than one, or the "
                         "fan-out this fixture exists to test is a fan-out "
                         "of one")
    ap.add_argument("--chimera-fraction", type=float, default=0.35,
                    help="share of reads joined head to tail, so the "
                         "splitter has something to cut")
    ap.add_argument("--seed", type=int, default=20260918)
    args = ap.parse_args()

    out = args.outdir
    reads_dir = os.path.join(out, "run", "fastq_pass")
    os.makedirs(reads_dir, exist_ok=True)
    work = os.path.join(out, ".build")
    os.makedirs(work, exist_ok=True)

    rng = random.Random(args.seed)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)

    samples = mf.build_samplesheet(i7, i5, args.n_samples, args.leave_unused,
                                   "ONT01", "nanopore")

    mf.write_reference(ref, os.path.join(out, "reference.tsv"))
    mf.write_primers(os.path.join(out, "primers.tsv"))
    mf.write_samplesheet(samples, os.path.join(out, "samplesheet.csv"))
    mf.write_metadata(samples, os.path.join(out, "metadata.tsv"))
    # The manifest points at the MinKNOW run folder, not at a file: the
    # fastq_pass / fastq_fail rules are part of what the nanopore path
    # has to get right, so the fixture has the structure they act on.
    # RELATIVE, deliberately. An absolute path here would bake this
    # machine's directory layout into a fixture that gets committed and
    # then run somewhere else - which is exactly what happened on the
    # first attempt. It is also the only way the fixture exercises
    # relative-path resolution, which both branches must agree on and
    # for a while did not.
    mf.write_runs(os.path.join(out, "runs.tsv"), "ONT01", "nanopore", "run")

    sim_runs = os.path.join(work, "runs.abs.tsv")
    mf.write_runs(sim_runs, "ONT01", "nanopore", os.path.join(out, "run"))
    quiet(sim.main, ["--reference", os.path.join(out, "reference.tsv"),
                     "--samplesheet", os.path.join(out, "samplesheet.csv"),
                     "--runs", sim_runs,
                     "--assay", ASSAY, "--outdir", work,
                     "--reads-per-sample", str(args.reads_per_sample)])

    pooled = os.path.join(work, "ONT01", "ONT01.fastq.gz")
    with gzip.open(pooled, "rt") as fh:
        lines = fh.read().splitlines()
    recs = [lines[i:i + 4] for i in range(0, len(lines), 4)]

    # Ligate a share of them head to tail. This is what the splitter is
    # for: an uncut join carries one molecule's i5 and another's i7, a
    # pair that is usually a real well and never the right one.
    rng2 = random.Random(args.seed + 1)
    n_join = int(len(recs) * args.chimera_fraction) // 2 * 2
    joined, k = [], 0
    while k + 1 < n_join:
        a, b = recs[k], recs[k + 1]
        inv = rng2.random() < 0.4          # some joins head to head
        seq = a[1] + (dx.revcomp(b[1]) if inv else b[1])
        qual = a[3] + (b[3][::-1] if inv else b[3])
        joined.append([f"@J{k}", seq, "+", qual])
        k += 2
    pool = joined + recs[k:]
    rng2.shuffle(pool)

    per = (len(pool) + args.batches - 1) // args.batches
    written = 0
    for b in range(args.batches):
        part = pool[b * per:(b + 1) * per]
        if not part:
            continue
        name = os.path.join(reads_dir, f"FAX00001_pass_{b:04d}.fastq.gz")
        with gzip.open(name, "wt") as fh:
            for r in part:
                fh.write("\n".join(r) + "\n")
        written += 1

    print(f"nanopore fixture in {out}/")
    print(f"  {len(samples)} samples, {len(pool):,} reads "
          f"({len(joined):,} concatemers) in {written} batch file(s)")
    print(f"  reads: {reads_dir}")
    if written < 2:
        print("  WARNING: one batch file only - the fan-out this fixture "
              "exists to exercise would be a fan-out of one", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
