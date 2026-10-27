#!/usr/bin/env python3
"""
Assemble the run-level split summary from per-file shards.

When Nextflow fans the splitting stage out - one task per input FASTQ,
one core each, scheduled across the farm - no single task ever sees the
whole run. This reads the per-part statistics those tasks left behind
and produces the same split_summary.json that one task splitting
everything itself would have written.

"The same" is meant literally, and is tested: the aggregation is not
reimplemented here but imported from split_concatemers, because the two
quiet ways to get it wrong both produce believable numbers.
bases_retained_fraction is a ratio of sums, and averaging the per-file
fractions instead gives a figure that looks right and weights a
thousand-read batch equally with a million-read one. The
pieces-per-read distribution is a histogram to be summed bin by bin,
and concatenating or overwriting instead gives a distribution that
still sums to something plausible.

There is nothing to do for the reads themselves. The splitter writes
one output FASTQ per input FASTQ and demultiplexing walks a directory,
so the shards ARE the output - collecting them is Nextflow's staging,
not a merge, and no bytes are rewritten.

  merge_split.py --parts split/parts --file-report file_report.json \
                 --outdir split
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from split_concatemers import build_summary  # noqa: E402


def load_parts(paths: list[str]) -> list[dict]:
    """Every per-part statistics file, in a defined order.

    Sorted by source filename rather than taken in the order the
    workflow happened to hand them over. Nothing in the arithmetic
    depends on order - they are all sums - but the summary is compared
    between runs and between the two routes through this stage, and a
    file whose key order wanders is a diff that has to be read before
    it can be dismissed.
    """
    out = []
    for path in paths:
        if os.path.isdir(path):
            found = sorted(glob.glob(os.path.join(path, "*.json")))
            if not found:
                raise SystemExit(f"no part statistics found in {path}")
            for f in found:
                with open(f, encoding="utf-8") as fh:
                    out.append(json.load(fh))
        else:
            with open(path, encoding="utf-8") as fh:
                out.append(json.load(fh))
    return sorted(out, key=lambda r: os.path.basename(r.get("source", "")))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--parts", required=True, nargs="+",
                    help="per-part JSON files, or a directory of them")
    ap.add_argument("--file-report",
                    help="file_report.json from the listing step. Without "
                         "it the summary cannot say what was excluded")
    ap.add_argument("--min-length", type=int, default=100)
    ap.add_argument("--anchor-max-edits", type=int, default=2)
    ap.add_argument("--threads", type=int, default=1,
                    help="threads each shard used; 1 when a workflow "
                         "fans this stage out by file")
    ap.add_argument("--expect-parts", type=int, default=0,
                    help="how many shards should be present. A shard that "
                         "silently went missing would shrink the read "
                         "funnel's denominator, not raise an error")
    ap.add_argument("--outdir", default="split")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    results = load_parts(args.parts)
    if not results:
        raise SystemExit("no part statistics to merge")

    # A missing shard is the failure this step exists to be exposed to.
    # Every downstream percentage is taken against reads_in, so losing
    # one part of seventy-four does not produce an error - it produces a
    # read funnel that balances against a smaller run.
    if args.expect_parts and len(results) != args.expect_parts:
        raise SystemExit(
            f"expected {args.expect_parts} shard(s), found {len(results)}.\n"
            "Refusing to write a run summary from an incomplete set: the "
            "read funnel below it would balance against the wrong total.")

    # Two shards reporting the same source file means one input was
    # processed twice, which double-counts its reads.
    sources = [r.get("source", "") for r in results]
    dupes = sorted({s for s in sources if sources.count(s) > 1})
    if dupes:
        raise SystemExit(
            f"{len(dupes)} input file(s) appear in more than one shard: "
            f"{', '.join(os.path.basename(d) for d in dupes[:5])}.\n"
            "Their reads would be counted twice.")

    file_report = {}
    if args.file_report and os.path.exists(args.file_report):
        with open(args.file_report, encoding="utf-8") as fh:
            file_report = json.load(fh)

    summary = build_summary(results, file_report, args.min_length,
                            args.anchor_max_edits, args.threads)

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "split_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    if not args.quiet:
        n, p = summary["reads_in"], summary["pieces_out"]
        print(f"{len(results)} shard(s) merged: {n:,} reads -> {p:,} pieces "
              f"({summary['pieces_per_read']:.2f} per read)")
        none = summary.get("reads_with_no_piece_kept", 0)
        if none:
            print(f"  {none:,} read(s) left no piece at or above "
                  f"{args.min_length}bp")
        bc = summary.get("basecalling", {})
        if bc.get("basecall_model") == "MIXED":
            print(f"  WARNING: {bc['warning']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
