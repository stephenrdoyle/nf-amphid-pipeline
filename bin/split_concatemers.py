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
import hashlib
import json
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    find_anchor_leftmost,
    load_assay,
    open_reads,
    read_run_metadata,
    resolve_read_files,
    revcomp,
)


def code_fingerprint(**settings) -> str:
    """What produced a part, so a resumed run cannot blend two versions.

    Resumption works by skipping any input file whose counts already
    exist. That is only safe while the thing being resumed is the same
    thing: change the cut rule, rerun, and every file the old version had
    finished is kept, giving a run whose recovery depends on which parts
    happened to complete before the code changed. It looks like a
    successful resume and the summary numbers are a blend of two
    algorithms.

    So each part records a fingerprint of the code and of the settings
    that affect its output, and a part is reused only if both still
    match. Anything else is recomputed.
    """
    h = hashlib.sha256()
    here = os.path.dirname(os.path.abspath(__file__))
    for name in ("split_concatemers.py", "ampseqlib.py"):
        try:
            with open(os.path.join(here, name), "rb") as fh:
                h.update(fh.read())
        except OSError:
            h.update(name.encode())
    for k in sorted(settings):
        h.update(f"{k}={settings[k]}".encode())
    return h.hexdigest()[:16]


# An i5 anchor opens a forward unit; an i7 anchor closes one. Reverse
# complemented, the roles swap ends: rc(i7) opens a unit lying backwards,
# rc(i5) closes it.
OPENS = ("5+", "7-")
CLOSES = ("7+", "5-")


def scan_all(text: str, pattern: str, max_edits: int, limit: int = 24):
    """Every non-overlapping occurrence, left to right.

    Each search asks for the FIRST copy in the remainder, not the best
    one. Asking for the best is how this used to work, and it skipped
    anchors: a read of [unit A whose i7 anchor carries one error][a
    fragment of B ending in a clean i7 anchor] presents B's anchor as
    the best match, the search resumes past it, and A's anchor is never
    seen. With one CLOSE where there were two, the CLOSE -> CLOSE cut
    cannot fire, and demultiplexing reads A's i5 with B's i7 - the
    hybrid pair this whole stage exists to prevent. Nanopore puts an
    error in a 14bp anchor often enough that this was not a corner case.
    """
    hits = []
    off = 0
    while len(hits) < limit and off < len(text):
        h = find_anchor_leftmost(text[off:], pattern, max_edits)
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


def cut_positions(seq: str, hits, head_len: dict, tail_len: dict) -> list[int]:
    """Every boundary between one molecule and the next.

    Within a single amplicon the anchors run OPEN then CLOSE - an i5 end
    then an i7 end - so that pairing is the one adjacency that is not a
    junction. Every other adjacency is:

        CLOSE -> OPEN   two whole molecules joined; cut between them
        OPEN  -> OPEN   the first molecule lost its i7 end; cut before
                        the second one starts
        CLOSE -> CLOSE  the SECOND molecule lost its i5 end; cut after
                        the first one finishes

    That last case was missed at first, and it is not a rare one. A read
    of [complete A][fragment of B carrying only its i7 end] has no second
    opening anchor at all, so nothing marked it as needing a cut - yet its
    outer ends are A's i5 and B's i7, a hybrid pair. The junction check
    correctly refuses the whole read, and molecule A, which is intact and
    assignable, is thrown away with it. Cutting after A's tail recovers
    it; the remainder has one end and is discarded downstream as tier B,
    which is the honest outcome for a fragment.

    A cut is placed at a unit's own boundary - just before an opening
    unit's outer tail, or just after a closing unit's - so each piece
    keeps its complete end structure and can be demultiplexed normally.
    """
    cuts: list[int] = []

    def add(p: int, floor: int) -> None:
        lo = max((cuts[-1] if cuts else 0) + 1, floor)
        p = max(lo, min(p, len(seq) - 1))
        if p > (cuts[-1] if cuts else 0):
            cuts.append(p)

    for i in range(1, len(hits)):
        prev_pos, prev_end, prev_kind = hits[i - 1]
        pos, end, kind = hits[i]
        if kind in OPENS:
            # Cut before this unit's outer tail, never earlier than the
            # end of the previous anchor - that would move sequence from
            # one molecule onto the other.
            add(pos - head_len[kind], prev_end)
        elif prev_kind in CLOSES:
            # Two closes in a row: the previous unit is complete, so cut
            # once past its own tail and leave the remainder behind.
            add(prev_end + tail_len[prev_kind], prev_end)
    return cuts


def split_read(seq: str, qual: str, a5: str, a7: str, edits: int,
               head_len: dict, tail_len: dict) -> list[tuple[str, str]]:
    hits = scan_hits(seq, a5, a7, edits)
    cuts = cut_positions(seq, hits, head_len, tail_len)
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
    head_len, tail_len = job["head_len"], job["tail_len"]
    min_length = job["min_length"]

    n_reads = n_pieces = n_dropped = n_cut = n_none = 0
    bases_in = bases_out = 0
    # Counted on the pieces that SURVIVE the length floor, because that is
    # what a read became. Counting them before the floor made a read whose
    # cut left a one-base remainder (the clamp in cut_positions puts the
    # last cut at most one base from the end) count as "split" while
    # writing one piece - a cut rate the output does not show.
    per_read: Counter = Counter()
    per_read_cut: Counter = Counter()

    tmp = dest + ".partial"
    with gzip.open(tmp, "wt") as out:
        for read_id, seq, qual in open_reads(src):
            n_reads += 1
            bases_in += len(seq)
            pieces = split_read(seq, qual, a5, a7, edits, head_len, tail_len)
            kept = [p for p in pieces if len(p[0]) >= min_length]
            n_dropped += len(pieces) - len(kept)
            if len(pieces) > 1:
                n_cut += 1
            # A read every piece of which is below the floor leaves nothing
            # behind. Without its own count it vanished: reads_in could
            # not be reconciled against anything written.
            if not kept:
                n_none += 1
            per_read[min(len(kept), 5)] += 1
            per_read_cut[min(len(pieces), 5)] += 1
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
        # Reads that came out as more than one usable piece.
        "reads_split": sum(v for k, v in per_read.items() if k > 1),
        # Reads the cut rule fired on at all, whatever the floor then
        # kept - the measure of what the anchors saw, not of the output.
        "reads_cut": n_cut,
        "reads_with_no_piece_kept": n_none,
        "pieces_dropped_below_min_length": n_dropped,
        # Kept pieces per read; bin 0 is reads_with_no_piece_kept, so the
        # bins sum to reads_in.
        "pieces_per_read_distribution": dict(per_read),
        "pieces_per_read_distribution_before_min_length": dict(per_read_cut),
        "bases_in": bases_in,
        "bases_out": bases_out,
        # Recorded per part, not once for the run: when Nextflow fans
        # this stage out, the step that writes the run summary never
        # sees an original read header. Splitting rewrites read ids, so
        # if it is not captured here it cannot be recovered at all.
        "basecalling": read_run_metadata([src]),
        # The cut geometry this part used. The index length sets where
        # every cut lands, and the workflow does not pass it, so the
        # summary has to be able to say which value was in force.
        "index_length": job.get("index_length"),
        "barcode_repeats": job.get("barcode_repeats"),
        "code": job["code"],
    }
    with open(job["stats"], "w", encoding="utf-8") as fh:
        json.dump(stats, fh)
    return stats


def _one_or_all(results: list, key: str):
    vals = sorted({r.get(key) for r in results if r.get(key) is not None})
    return vals[0] if len(vals) == 1 else (vals or None)


def build_summary(results: list, file_report: dict, min_length: int,
                  anchor_max_edits: int, threads: int) -> dict:
    """Aggregate per-part statistics into the run-level summary.

    Extracted so that both routes through this stage produce the
    identical object: the single task that splits every file itself,
    and Nextflow fanning one task out per file and merging afterwards.
    Two implementations of this arithmetic would be two chances to get
    it wrong, and the ways to get it wrong are quiet ones -
    bases_retained_fraction is a ratio of sums, not a mean of ratios,
    and the pieces-per-read distribution is a histogram to be summed
    bin by bin. Either mistake yields a plausible number.
    """
    n_reads = sum(r["reads_in"] for r in results)
    n_pieces = sum(r["pieces_out"] for r in results)
    n_dropped = sum(r["pieces_dropped_below_min_length"] for r in results)
    bases_in = sum(r["bases_in"] for r in results)
    bases_out = sum(r["bases_out"] for r in results)
    per_read: Counter = Counter()
    per_read_cut: Counter = Counter()
    for r in results:
        for k, v in r["pieces_per_read_distribution"].items():
            per_read[int(k)] += v
        # .get, not [], for the counters added later: a part written
        # without them (the workflow's stub, for one) is still a part, and
        # its contribution to them is nothing rather than a KeyError.
        for k, v in r.get("pieces_per_read_distribution_before_min_length",
                          {}).items():
            per_read_cut[int(k)] += v

    # Basecalling provenance. Each part records what produced its own
    # reads; splitting rewrites read ids, so this is the last stage that
    # can see it at all. Parts should agree - they are batches of one
    # MinKNOW run - and if they do not, that is a mixed-basecall run and
    # worth saying so rather than quietly reporting whichever part was
    # read first, because the model sets the error profile and so the
    # barcode recovery.
    seen = {}
    for r in results:
        bc = r.get("basecalling") or {}
        model = bc.get("basecall_model")
        if model:
            seen.setdefault(model, bc)
    if len(seen) == 1:
        basecalling = next(iter(seen.values()))
    elif not seen:
        basecalling = {}
    else:
        basecalling = {
            "basecall_model": "MIXED",
            "models_seen": sorted(seen),
            "warning": (f"{len(seen)} basecalling models across the input "
                        "files. Error profiles differ between models, so "
                        "barcode recovery and composition are not "
                        "comparable across these files."),
        }

    return {
        "reads_in": n_reads,
        "pieces_out": n_pieces,
        "pieces_per_read": n_pieces / n_reads if n_reads else 0,
        "reads_split": sum(r["reads_split"] for r in results),
        "reads_cut": sum(r.get("reads_cut", 0) for r in results),
        "reads_with_no_piece_kept": sum(
            r.get("reads_with_no_piece_kept", 0) for r in results),
        "pieces_dropped_below_min_length": n_dropped,
        "min_length": min_length,
        "pieces_per_read_distribution": dict(sorted(per_read.items())),
        "pieces_per_read_distribution_before_min_length":
            dict(sorted(per_read_cut.items())),
        "bases_in": bases_in,
        "bases_out": bases_out,
        # Cutting removes nothing but the pieces below the length floor.
        # If bases go missing beyond that, the cut arithmetic is wrong and
        # the output is not a rearrangement of the input any more.
        "bases_retained_fraction": bases_out / bases_in if bases_in else 0,
        "input_files": file_report,
        "basecalling": basecalling,
        "parts": len(results),
        "settings": {"anchor_max_edits": anchor_max_edits,
                     "threads": threads,
                     # From the parts, so the fanned-out route reports
                     # them too. More than one value is listed rather
                     # than one picked: it would mean shards cut on
                     # different geometry.
                     "index_length": _one_or_all(results, "index_length"),
                     "barcode_repeats": _one_or_all(results,
                                                    "barcode_repeats")},
    }


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
    ap.add_argument("--part-index", type=int, default=None,
                    help="the index this input file has in the whole run. "
                         "A shard sees one file and would otherwise number "
                         "it 0000, so every shard of a fanned-out run would "
                         "name its output alike - and two MinKNOW batches "
                         "in different barcode folders can share a "
                         "basename. Passing the index from --list-only "
                         "makes a shard's output name the one the "
                         "whole-run route would have given it")
    ap.add_argument("--list-only", action="store_true",
                    help="write read_files.txt and file_report.json for the "
                         "resolved inputs, then stop. Lets a workflow build "
                         "one task per input file from this tool's own "
                         "discovery rules rather than a glob of its own")
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
    # The mirror image: how far past a closing anchor that unit's own
    # index and outer tail extend.
    tail_len = {
        "7+": span + len(cfg["p7_tail"]),
        "5-": span + len(cfg["p5_tail"]),
    }

    read_files, file_report = resolve_read_files(
        args.reads, prefer_pass=not args.include_failed)
    if not read_files:
        print(f"no FASTQ found under {', '.join(args.reads)}", file=sys.stderr)
        return 1

    # Which files constitute a nanopore run is a decision with rules -
    # fastq_pass in, fastq_fail out, pod5 directories ignored - and they
    # live in resolve_read_files(). A workflow that wants one task per
    # input file asks here rather than globbing for itself, so there is
    # one answer to "what is the input" and not two that agree until
    # somebody reorganises a run directory.
    if args.list_only:
        # index<TAB>path, so the shard that handles a file can name its
        # output exactly as the whole-run route would have.
        with open(os.path.join(args.outdir, "read_files.txt"), "w",
                  encoding="utf-8") as fh:
            for i, path in enumerate(read_files):
                fh.write(f"{i}\t{os.path.abspath(path)}\n")
        with open(os.path.join(args.outdir, "file_report.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(file_report, fh, indent=2)
        if not args.quiet:
            print(f"{len(read_files)} input file(s) listed")
        return 0

    # One output file per input file, in a directory the next stage reads
    # whole. Merging them into a single FASTQ would mean rewriting tens of
    # gigabytes for no benefit: demultiplexing walks a directory anyway.
    reads_dir = os.path.join(args.outdir, "reads")
    stats_dir = os.path.join(args.outdir, "parts")
    os.makedirs(reads_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    fingerprint = code_fingerprint(
        anchor_max_edits=args.anchor_max_edits, min_length=args.min_length,
        barcode_repeats=args.barcode_repeats, index_length=args.index_length,
        assay=args.assay or "")

    if args.part_index is not None and len(read_files) != 1:
        raise SystemExit(
            f"--part-index names the position of ONE file, but "
            f"{len(read_files)} were resolved. A shard processes a single "
            "input file.")

    jobs, done, stale = [], [], 0
    for i, src in enumerate(read_files):
        n = args.part_index if args.part_index is not None else i
        stem = f"{n:04d}_" + os.path.basename(src).split(".")[0]
        stats_path = os.path.join(stats_dir, stem + ".json")
        dest = os.path.join(reads_dir, stem + ".fastq.gz")
        if os.path.exists(stats_path) and os.path.exists(dest):
            with open(stats_path, encoding="utf-8") as fh:
                prev = json.load(fh)
            if prev.get("code") == fingerprint:
                done.append(prev)
                continue
            # Produced by different code or different settings. Redo it
            # rather than mix two algorithms in one run.
            stale += 1
        jobs.append({"src": src, "dest": dest, "stats": stats_path,
                     "a5": a5, "a7": a7, "edits": args.anchor_max_edits,
                     "head_len": head_len, "tail_len": tail_len,
                     "min_length": args.min_length, "code": fingerprint,
                     "index_length": args.index_length,
                     "barcode_repeats": args.barcode_repeats})

    threads = args.threads or min(8, os.cpu_count() or 1)
    if not args.quiet:
        print(f"{len(read_files)} input file(s), {len(done)} already done, "
              f"{len(jobs)} to do on {threads} thread(s)")
        if stale:
            print(f"  {stale} previously-finished file(s) were made by "
                  "different code or settings\n  and will be redone - a "
                  "resumed run must not blend two versions")
        # Name them. A count alone cannot distinguish "the run folder has
        # one batch" from "--reads is pointing at the wrong place", and
        # those look identical in every other way.
        for p in read_files[:3]:
            print(f"    {p}")
        if len(read_files) > 3:
            print(f"    ... and {len(read_files) - 3} more")
        for note in file_report.get("notes", []):
            print(f"  note: {note}")
        # Work is divided by input file, so the file count is a ceiling on
        # how much of the requested parallelism can actually be used. A
        # MinKNOW run arrives as dozens of batches and this never fires;
        # a run someone has already concatenated arrives as one file, and
        # then the thread count is a fiction worth saying out loud rather
        # than leaving to be inferred from the wall-clock time.
        if jobs and threads > len(jobs):
            print(f"  note: only {len(jobs)} file(s) to process, so at most "
                  f"{len(jobs)} of the {threads} threads can be used. Work "
                  "is\n        divided by input file; splitting the input "
                  "into more files would\n        use the rest.")

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

    summary = build_summary(results, file_report, args.min_length,
                            args.anchor_max_edits, threads)

    with open(os.path.join(args.outdir, "split_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    if not args.quiet:
        split = summary["reads_split"]
        n_reads = summary["reads_in"]
        n_pieces = summary["pieces_out"]
        n_dropped = summary["pieces_dropped_below_min_length"]
        print(f"\n{n_reads:,} reads -> {n_pieces:,} pieces "
              f"({summary['pieces_per_read']:.2f} per read)")
        print(f"  reads split          {split:>9,}  "
              f"{100 * split / max(1, n_reads):5.1f}%")
        print(f"  pieces below {args.min_length}bp    {n_dropped:>9,}  dropped")
        print(f"  reads with no piece  {summary['reads_with_no_piece_kept']:>9,}  "
              f"left at or above {args.min_length}bp")
        print(f"  bases retained       {100 * summary['bases_retained_fraction']:>8.1f}%")
        print(f"\n{reads_dir}/  ({len(results)} file(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
