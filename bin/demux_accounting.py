#!/usr/bin/env python3
"""
The platform layer's demux accounting (contract A, docs/CONTRACTS.md).

Every platform adapter hands the protocol layer the same three things:
each sample's reads, their read type, and this - one JSON per run saying
where the reads from the instrument went before any sample saw them, in
fates every platform shares, plus the plate position of each well where
the platform has one. Synthesis reads it, never a platform's own
summaries, so a new platform needs only its adapter.

  illumina          demultiplexed by the sequencer: nothing is known of
                    reads that reached no sample, so reads_in and the
                    fates are null and per-sample counts are left to
                    assignment (assign_db.json reads_in).
  nanopore          split at internal adapters, then demultiplexed here
                    (split_concatemers.py, demux_ont.py): every fate is
                    counted, from demux_summary.json's tiers.
  predemultiplexed  demultiplexed upstream (lima, Dorado, a sequencing
                    centre): as Illumina - nothing is known of what it
                    discarded.
  pacbio            demultiplexed here by the samplesheet's tags
                    (demux_pacbio.py, or lima read back by it): every
                    fate but internal_adapter counted, from the same
                    tiers as nanopore's - CCS has removed the adapters,
                    so a read is one molecule and there is no tier E.

The platform's own summaries ride along unchanged under "native", for
the detail a shared fate cannot hold (nanopore's end reasons, the
combination grid, the basecaller) - synthesis reads them from there.

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import read_table  # noqa: E402

CONTRACT = "amphid demux accounting"
VERSION = 1
# The shared fates, and the demux_ont tier each comes from.
FATES = {"assigned": "A", "one_end_only": "B", "chimeric": "C",
         "unmatched": "D", "internal_adapter": "E"}
DEMULTIPLEXED_BY = {
    "illumina": "the sequencer",
    "nanopore": "amphid (split_concatemers.py, demux_ont.py)",
    "predemultiplexed": "upstream of amphid",
    "pacbio": "amphid (demux_pacbio.py)",
}


def load(path: str | None) -> dict | None:
    if not path or not os.path.exists(path) or \
            os.path.basename(path).startswith("NO_"):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def accounting(run_id: str, platform: str, read_type: str,
               declared: list[dict], demux: dict | None,
               split: dict | None) -> dict:
    """The accounting of one run (module docstring)."""
    rows = [r for r in declared if not run_id or r.get("run_id") == run_id]
    plate = {}
    for r in rows:
        try:
            col = int(r.get("plate_col") or 0)
        except ValueError:
            col = 0
        if r.get("plate_row") and col:
            plate[r["sample_id"]] = {"row": r["plate_row"], "col": col}
    fates: dict = {}
    reads_in = assigned = None
    per_sample = {r["sample_id"]: None for r in rows}
    if demux:
        tiers = demux.get("by_tier") or {}
        # The fates the adapter counts: a tier it does not have (PacBio's
        # E) is not a fate of its reads, rather than zero of them.
        fates = {f: int(tiers.get(t, 0)) for f, t in FATES.items()
                 if t in tiers}
        reads_in = demux.get("total_reads")
        assigned = demux.get("assigned_reads", fates["assigned"])
        got = demux.get("per_sample") or {}
        for r in rows:
            per_sample[r["sample_id"]] = got.get(
                r["sample_id"], got.get(r.get("sanitised_id", ""), 0))
    return {
        "contract": CONTRACT,
        "version": VERSION,
        "run_id": run_id,
        "platform": platform,
        "read_type": read_type,
        # The adapter's own word where it gives one (lima, on PacBio).
        "demultiplexed_by": (demux or {}).get("demultiplexed_by")
        or DEMULTIPLEXED_BY.get(platform, platform),
        # Reads off the instrument, where the adapter saw them; on
        # nanopore that is before splitting at internal adapters.
        "reads_off_instrument": (split or {}).get("reads_in", reads_in),
        "reads_in": reads_in,
        "assigned": assigned,
        "fates": fates,
        "balanced": (sum(fates.values()) == reads_in) if fates else None,
        "per_sample": per_sample,
        "plate": plate,
        "native": {"demux": demux, "split": split},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--preflight-dir", required=True,
                    help="pre-flight's outputs: file_resolution.tsv for the "
                    "samples and plate positions, read_types.tsv for the "
                    "run's platform and read type")
    ap.add_argument("--demux-summary",
                    help="nanopore or pacbio: demux_summary.json")
    ap.add_argument("--split-summary", help="nanopore: split_summary.json")
    ap.add_argument("--out", default="demux_accounting.json")
    args = ap.parse_args(argv)

    declared = read_table(os.path.join(args.preflight_dir,
                                       "file_resolution.tsv"))
    rt_path = os.path.join(args.preflight_dir, "read_types.tsv")
    run = next((r for r in (read_table(rt_path)
                            if os.path.exists(rt_path) else [])
                if r.get("run_id") == args.run_id), {})
    platform = run.get("platform") or next(
        (r.get("platform", "") for r in declared
         if r.get("run_id") == args.run_id), "")
    acc = accounting(args.run_id, platform, run.get("read_type", ""),
                     declared, load(args.demux_summary),
                     load(args.split_summary))
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(acc, fh, indent=2)
    print(f"{args.run_id}: {platform} ({acc['read_type'] or '?'}), "
          f"{len(acc['per_sample'])} samples, "
          f"{len(acc['plate'])} with plate positions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
