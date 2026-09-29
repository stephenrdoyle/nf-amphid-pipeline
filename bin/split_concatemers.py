#!/usr/bin/env python3
"""Cut nanopore reads that contain more than one amplicon.

Demultiplexing assumes one read is one molecule: it reads an index at
each end of the read and pairs them. Ligation does not honour that. In a
real MinION run of this assay, 40% of reads carried more than one
adapter structure and there were 1.43 amplicon units per read.

That costs twice over. Every unit past the first is invisible, because
only the outer ends are looked at - so yield falls. And worse, a
head-to-tail join presents the first molecule's i5 and the SECOND
molecule's i7, a pair of real declared indices that very likely belongs
to some well on the plate but to neither molecule. It demultiplexes
cleanly into the wrong sample.

`demux_ont.py` refuses those reads (tier E). This recovers them instead,
by cutting the read at each junction and letting each piece be
demultiplexed on its own ends.

WHERE THE CUT GOES. Each anchor hit is labelled by what it is: an i5
anchor opens a forward unit, an i7 anchor closes one, and their reverse
complements do the same for a unit lying the other way round. A unit's
index and outer tail sit at a known distance from its anchor, so the
start of any unit after the first is a computable position - and that is
where the read is cut. Each piece then carries its own complete end
structure, and nothing is spliced or reconstructed: the cut is the only
edit.

WHAT IT CANNOT DO. A junction between two TRUNCATED fragments leaves no
adapter behind, so there is nothing to find and nothing to cut. Such a
read still presents one usable end from each of two molecules and still
demultiplexes to a well that never held it. No adapter-based method
reaches it. Leaving unused index combinations on the plate is what would
make the residue measurable.

    bin/split_concatemers.py --reads <minknow dir|fastq> \
        --assay assets/assay/ampseq_helminth.json --outdir split

Self-contained: standard library only.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    find_anchor_cost,
    load_assay,
    open_reads,
    resolve_read_files,
    revcomp,
)

# An i5 anchor opens a forward unit; an i7 anchor closes one. Reverse
# complemented, the roles swap ends: rc(i7) opens a unit lying backwards,
# rc(i5) closes it.
OPENS = ("5+", "7-")
CLOSES = ("7+", "5-")


def scan_all(text: str, pattern: str, max_edits: int, limit: int = 24):
    """Every non-overlapping occurrence, left to right."""
    hits = []
    off = 0
    while len(hits) < limit and off < len(text):
        h = find_anchor_cost(text[off:], pattern, max_edits)
        if h is None:
            break
        s, e, _c = h
        hits.append((off + s, off + e))
        off += max(e, s + 1)
    return hits


def scan_hits(seq: str, a5: str, a7: str, edits: int):
    """All four anchor orientations, in read order."""
    hits = []
    for kind, pat in (("5+", a5), ("7+", a7),
                      ("5-", revcomp(a5)), ("7-", revcomp(a7))):
        for s, e in scan_all(seq, pat, edits):
            hits.append((s, e, kind))
    hits.sort()
    # Two anchors cannot occupy the same bases. If they overlap, the
    # later one is a chance match inside the earlier - drop it rather
    # than cut on it.
    kept = []
    for h in hits:
        if kept and h[0] < kept[-1][1]:
            continue
        kept.append(h)
    return kept


def cut_positions(seq: str, hits, head_len: dict) -> list[int]:
    """Where each unit after the first begins.

    A cut goes immediately before the outer tail of every unit that
    opens after the read has already started. That keeps the preceding
    unit's own tail and index intact on its side of the cut, which is
    what lets both pieces be demultiplexed normally afterwards.
    """
    cuts: list[int] = []
    for i, (pos, _end, kind) in enumerate(hits):
        if i == 0 or kind not in OPENS:
            continue
        p = pos - head_len[kind]
        lo = (cuts[-1] if cuts else 0) + 1
        # Never cut before the previous unit's anchor: that would move
        # sequence from one unit onto the other.
        lo = max(lo, hits[i - 1][1])
        p = max(lo, min(p, len(seq) - 1))
        if p > (cuts[-1] if cuts else 0):
            cuts.append(p)
    return cuts


def split_read(seq: str, qual: str, a5: str, a7: str, edits: int,
               head_len: dict) -> list[tuple[str, str]]:
    hits = scan_hits(seq, a5, a7, edits)
    cuts = cut_positions(seq, hits, head_len)
    if not cuts:
        return [(seq, qual)]
    pieces = []
    prev = 0
    for c in cuts + [len(seq)]:
        if c - prev > 0:
            pieces.append((seq[prev:c], qual[prev:c]))
        prev = c
    return pieces


def _progress(done: int, total: int, results: list, started: float) -> None:
    """Say how far in, how fast, and how long left.

    A stage with no output is indistinguishable from a hung one, and this
    stage can run for hours. Twelve of them were spent on a run that gave
    no sign of whether it was halfway or barely started.
    """
    reads = sum(r["reads_in"] for r in results)
    dt = max(1e-6, time.time() - started)
    rate = reads / dt
    left = (total - done) * (dt / max(1, done))
    sys.stderr.write(
        f"\r  {done}/{total} files  {reads:,} reads  "
        f"{rate:,.0f} reads/s  ~{left / 60:.0f} min left    ")
    sys.stderr.flush()


def split_one_file(job: dict) -> dict:
    """Split one input FASTQ into one output FASTQ, and report on it.

    One input file is the unit of work and of resumption. A MinKNOW run
    is already delivered as dozens of batch files, so this needs no
    chunking scheme of its own: the batches are the chunks. A worker that
    finishes writes its counts beside its output, and a rerun skips every
    file whose counts are already there - so a job killed by a wall-clock
    limit resumes at the file it died on rather than at the beginning.
    """
    src, dest = job["src"], job["dest"]
    a5, a7, edits = job["a5"], job["a7"], job["edits"]
    head_len, min_length = job["head_len"], job["min_length"]

    n_reads = n_pieces = n_dropped = 0
    bases_in = bases_out = 0
    per_read: Counter = Counter()

    tmp = dest + ".partial"
    with gzip.open(tmp, "wt") as out:
        for read_id, seq, qual in open_reads(src):
            n_reads += 1
            bases_in += len(seq)
            pieces = split_read(seq, qual, a5, a7, edits, head_len)
            kept = [p for p in pieces if len(p[0]) >= min_length]
            n_dropped += len(pieces) - len(kept)
            per_read[min(len(pieces), 5)] += 1
            for k, (s, q) in enumerate(kept):
                pid = read_id if len(kept) == 1 else f"{read_id}_u{k + 1}"
                out.write(f"@{pid}\n{s}\n+\n{q}\n")
                n_pieces += 1
                bases_out += len(s)
    # Rename only once the file is complete, so an interrupted worker
    # never leaves behind a truncated output that a resume would trust.
    os.replace(tmp, dest)

    stats = {
        "source": src,
        "reads_in": n_reads,
        "pieces_out": n_pieces,
        "reads_split": sum(v for k, v in per_read.items() if k > 1),
        "pieces_dropped_below_min_length": n_dropped,
        "pieces_per_read_distribution": dict(per_read),
        "bases_in": bases_in,
        "bases_out": bases_out,
    }
    with open(job["stats"], "w", encoding="utf-8") as fh:
        json.dump(stats, fh)
    return stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--assay")
    ap.add_argument("--include-failed", action="store_true")
    ap.add_argument("--anchor-max-edits", type=int, default=2)
    # load_assay only carries through keys that exist in its DEFAULTS, and
    # these two are not among them - so reading them from the config would
    # silently use the fallback whatever the config said. Explicit, and
    # matching demux_ont's own flag.
    ap.add_argument("--barcode-repeats", type=int, default=2)
    ap.add_argument("--index-length", type=int, default=8)
    ap.add_argument("--min-length", type=int, default=100,
                    help="drop pieces shorter than this. A fragment with no "
                         "room for an insert cannot be assigned and cannot "
                         "be analysed, so keeping it only inflates counts")
    ap.add_argument("--outdir", default="split")
    ap.add_argument("--threads", type=int, default=0,
                    help="worker processes; 0 picks min(8, cores). Work is "
                         "divided by input file, and a MinKNOW run already "
                         "arrives as dozens of batch files")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    cfg = load_assay(args.assay)
    a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
    span = args.index_length * args.barcode_repeats
    head_len = {
        "5+": span + len(cfg["p5_tail"]),
        "7-": span + len(cfg["p7_tail"]),
    }

    read_files, file_report = resolve_read_files(
        args.reads, prefer_pass=not args.include_failed)
    if not read_files:
        print(f"no FASTQ found under {', '.join(args.reads)}", file=sys.stderr)
        return 1

    # One output file per input file, in a directory the next stage reads
    # whole. Merging them into a single FASTQ would mean rewriting tens of
    # gigabytes for no benefit: demultiplexing walks a directory anyway.
    reads_dir = os.path.join(args.outdir, "reads")
    stats_dir = os.path.join(args.outdir, "parts")
    os.makedirs(reads_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    jobs, done = [], []
    for i, src in enumerate(read_files):
        stem = f"{i:04d}_" + os.path.basename(src).split(".")[0]
        stats_path = os.path.join(stats_dir, stem + ".json")
        dest = os.path.join(reads_dir, stem + ".fastq.gz")
        if os.path.exists(stats_path) and os.path.exists(dest):
            with open(stats_path, encoding="utf-8") as fh:
                done.append(json.load(fh))
            continue
        jobs.append({"src": src, "dest": dest, "stats": stats_path,
                     "a5": a5, "a7": a7, "edits": args.anchor_max_edits,
                     "head_len": head_len, "min_length": args.min_length})

    threads = args.threads or min(8, os.cpu_count() or 1)
    if not args.quiet:
        print(f"{len(read_files)} input file(s), {len(done)} already done, "
              f"{len(jobs)} to do on {threads} thread(s)")

    results = list(done)
    started = time.time()
    if jobs:
        if threads > 1:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                for k, stats in enumerate(pool.map(split_one_file, jobs), 1):
                    results.append(stats)
                    if not args.quiet:
                        _progress(k, len(jobs), results, started)
        else:
            for k, job in enumerate(jobs, 1):
                results.append(split_one_file(job))
                if not args.quiet:
                    _progress(k, len(jobs), results, started)
        if not args.quiet:
            print()

    n_reads = sum(r["reads_in"] for r in results)
    n_pieces = sum(r["pieces_out"] for r in results)
    n_dropped = sum(r["pieces_dropped_below_min_length"] for r in results)
    bases_in = sum(r["bases_in"] for r in results)
    bases_out = sum(r["bases_out"] for r in results)
    per_read = Counter()
    for r in results:
        for k, v in r["pieces_per_read_distribution"].items():
            per_read[int(k)] += v

    summary = {
        "reads_in": n_reads,
        "pieces_out": n_pieces,
        "pieces_per_read": n_pieces / n_reads if n_reads else 0,
        "reads_split": sum(r["reads_split"] for r in results),
        "pieces_dropped_below_min_length": n_dropped,
        "min_length": args.min_length,
        "pieces_per_read_distribution": dict(sorted(per_read.items())),
        "bases_in": bases_in,
        "bases_out": bases_out,
        # Cutting removes nothing but the pieces below the length floor.
        # If bases go missing beyond that, the cut arithmetic is wrong and
        # the output is not a rearrangement of the input any more.
        "bases_retained_fraction": bases_out / bases_in if bases_in else 0,
        "input_files": file_report,
        "parts": len(results),
        "settings": {"anchor_max_edits": args.anchor_max_edits,
                     "threads": threads},
    }
    with open(os.path.join(args.outdir, "split_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    if not args.quiet:
        split = summary["reads_split"]
        print(f"\n{n_reads:,} reads -> {n_pieces:,} pieces "
              f"({summary['pieces_per_read']:.2f} per read)")
        print(f"  reads cut            {split:>9,}  "
              f"{100 * split / max(1, n_reads):5.1f}%")
        print(f"  pieces below {args.min_length}bp    {n_dropped:>9,}  dropped")
        print(f"  bases retained       {100 * summary['bases_retained_fraction']:>8.1f}%")
        print(f"\n{reads_dir}/  ({len(results)} file(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
