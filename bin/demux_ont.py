#!/usr/bin/env python3
"""
Nanopore demultiplexing for nf-helminth-amplicon.

The nanopore runs sequence Illumina libraries, so this is in-line barcode
demultiplexing on reads carrying full Nextera adapter structure - not
kit-based ONT barcoding. Four things shape the implementation, each of them
a decision that could reasonably have gone the other way:

1. ANCHOR, NEVER SCAN. The i5 and i7 index sets are only one substitution
   apart at their closest (S502 CTCTCTAT vs N707 CTCTCTAC), so an unanchored
   search cannot tell which end of a molecule it is reading. The barcode
   window is located by the adapter sequence abutting it, and only that
   window is tested.

2. SUBSTITUTIONS ONLY, INDELS DISABLED. The index sets have minimum Hamming
   distance 4 but minimum *edit* distance 2 (i7 N721 vs N723). One
   substitution is safely correctable; one indel is not. Since indels are
   nanopore's dominant error, allowing them would misassign.

3. THE TANDEM DUPLICATION IS INDEL INSURANCE, NOT EXTRA TOLERANCE. An indel
   makes one copy unmatchable and the other carries the read; a frame slip
   produces two copies that disagree, which is what detects it. With
   barcode_repeats=1 there is no second observation, and the rotational
   collisions in these sets become live misassignment routes.

4. TIERS, NOT A SINGLE PASS. Both ends agreeing is the confident case;
   single-ended reads are assigned but flagged, so the cost of accepting
   them can be measured rather than assumed.

Self-contained: standard library only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    DEFAULTS,
    open_reads_many,
    resolve_read_files,
    find_anchor,
    hamming,
    load_assay,
    open_reads,
    read_table,
    revcomp,
)

# --------------------------------------------------------------------------
# construct defaults, overridden by the assay config
# --------------------------------------------------------------------------

# Assignment tiers. A is the confident case; D is the residue.
TIERS = ("A", "B", "C", "D")

TIER_MEANING = {
    "A": "both ends agree and match a samplesheet pair",
    "B": "one end only - the other truncated or unreadable",
    "C": "ends disagree - chimeric",
    "D": "no confident match",
}


# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# barcode matching
# --------------------------------------------------------------------------


@dataclass
class BarcodeSet:
    """The declared barcodes for one end, plus the tolerance the set itself
    can support. Tolerance is derived from the set rather than supplied:
    a number someone typed cannot know the set's minimum distance."""

    by_seq: dict[str, str] = field(default_factory=dict)  # sequence -> name
    length: int = 0
    max_mismatch: int = 1

    @classmethod
    def build(cls, seqs: set[str], requested: int = 1) -> "BarcodeSet":
        seqs = {s for s in seqs if s}
        if not seqs:
            return cls()
        lengths = {len(s) for s in seqs}
        length = max(lengths) if lengths else 0
        by_seq = {s: s for s in seqs}

        safe = requested
        ordered = sorted(seqs)
        for i, a in enumerate(ordered):
            for b in ordered[i + 1 :]:
                if len(a) != len(b):
                    continue
                # k mismatches are only safely correctable if every pair is
                # more than 2k apart.
                safe = min(safe, (hamming(a, b) - 1) // 2)
        return cls(by_seq=by_seq, length=length, max_mismatch=max(0, safe))

    def match(self, window: str) -> str | None:
        """Exact first, then the single nearest within tolerance. A tie
        between two barcodes is no match at all - guessing would be worse
        than declining."""
        if window in self.by_seq:
            return window
        if self.max_mismatch < 1 or len(window) != self.length:
            return None
        best, best_d, ties = None, self.max_mismatch + 1, 0
        for cand in self.by_seq:
            d = hamming(window, cand)
            if d < best_d:
                best, best_d, ties = cand, d, 1
            elif d == best_d:
                ties += 1
        return best if best is not None and ties == 1 else None


@dataclass
class EndCall:
    barcode: str | None = None
    copies_seen: int = 0
    copies_agreed: bool = False
    anchor_pos: int = -1
    anchor_end: int = -1
    window_start: int = -1


def call_end(
    read: str,
    anchor: str,
    bc: BarcodeSet,
    repeats: int,
    search_window: int,
    barcode_before_anchor: bool,
    slack: int = 3,
    anchor_max_edits: int = 2,
) -> EndCall:
    """Locate the anchor, read the tandem barcode copies beside it, and
    require them to agree.

    Which side matters, and the two ends differ. The construct is

        P5 - i5 x n - [i5 anchor]...amplicon...[i7 anchor] - i7 x n - rc(P7)

    so the i5 copies sit immediately BEFORE their anchor (the anchor is the
    5' end of the PCR1 forward overhang) while the i7 copies sit immediately
    AFTER theirs (the anchor is the 3' end of the reverse-complemented PCR1
    reverse overhang). Assuming one side for both silently loses every i7.

    The copies are matched independently rather than concatenated. Joining
    them pools the error budget: two mismatches across 16bp, concentrated in
    one half, can push a read onto the wrong barcode. Two independent
    observations that agree is stronger evidence than one longer noisy one.
    """
    out = EndCall()
    if not anchor or not bc.length:
        return out

    # The anchor is searched only within a window from the read end it
    # belongs to. The construct puts it at a known offset, and scanning the
    # whole read would invite a chance match elsewhere.
    window = read[:search_window]
    exact = window.find(anchor)
    if exact >= 0:
        a_start, a_end = exact, exact + len(anchor)
    else:
        hit = find_anchor(window, anchor, anchor_max_edits)
        if hit is None:
            return out
        a_start, a_end = hit
    out.anchor_pos = a_start
    out.anchor_end = a_end

    span = bc.length * repeats
    # Measured from the alignment boundary, so an indel inside the anchor
    # does not displace the barcode window.
    start = a_start - span if barcode_before_anchor else a_end
    if start + slack < 0 or start + span - slack > len(read):
        # the anchor was found but the barcode runs off the end of the read
        return out
    out.window_start = start

    # Each copy is sought at its expected offset PLUS OR MINUS a little
    # slack, rather than at a fixed position. An indel between the anchor and
    # the barcode shifts every downstream base, so a fixed window turns one
    # indel into two failed copies - which is exactly the case the tandem
    # duplication exists to survive. Slack is small and anchored, so this is
    # still a handful of tested windows rather than a scan.
    calls = []
    for k in range(repeats):
        expect = start + k * bc.length
        best, best_d = None, bc.max_mismatch + 1
        for shift in range(-slack, slack + 1):
            a = expect + shift
            if a < 0 or a + bc.length > len(read):
                continue
            hit = bc.match(read[a : a + bc.length])
            if hit is None:
                continue
            d = hamming(read[a : a + bc.length], hit)
            if d < best_d:
                best, best_d = hit, d
        calls.append(best)

    seen = [c for c in calls if c]
    out.copies_seen = len(seen)
    if not seen:
        return out

    if len(set(seen)) == 1:
        out.barcode = seen[0]
        out.copies_agreed = len(seen) == repeats
    else:
        # copies disagree: a frame slip from an upstream indel produces
        # exactly this, and guessing which copy is right would defeat the
        # point of duplicating them
        out.barcode = None
        out.copies_agreed = False
    return out


# --------------------------------------------------------------------------
# per-read assignment
# --------------------------------------------------------------------------


@dataclass
class ReadCall:
    read_id: str
    tier: str
    sample_id: str | None
    i5: str | None
    i7: str | None
    orientation: str
    insert_start: int = -1
    insert_end: int = -1


def assign_read(
    read_id: str,
    seq: str,
    cfg: dict,
    i5_set: BarcodeSet,
    i7_set: BarcodeSet,
    pair_to_sample: dict[tuple[str, str], str],
    repeats: int,
    search_window: int,
    accept_single_end: bool,
) -> ReadCall:
    """Try both orientations; keep the one that resolves more ends.

    Nanopore reads arrive in either orientation, and the construct is not
    palindromic, so the i5 end is only at the start on one strand.
    """
    best: ReadCall | None = None
    for orient, s in (("+", seq), ("-", revcomp(seq))):
        five = call_end(
            s, cfg["i5_inner_anchor"], i5_set, repeats, search_window,
            barcode_before_anchor=True,
        )
        # The i7 anchor sits near the 3' end, so search the tail - and its
        # barcode follows the anchor rather than preceding it.
        tail_from = max(0, len(s) - search_window)
        three = call_end(
            s[tail_from:], cfg["i7_inner_anchor"], i7_set, repeats,
            search_window, barcode_before_anchor=False,
        )

        resolved = (1 if five.barcode else 0) + (1 if three.barcode else 0)
        call = _classify(
            read_id, five, three, orient, pair_to_sample, accept_single_end
        )
        # Trim past the whole overhang, not just the anchor: the mosaic end
        # is shared between the two overhangs, so it cannot identify an end,
        # but it is adapter all the same and does not belong in the insert.
        five_tail = max(0, len(cfg["pcr1_overhang_fwd"]) - len(cfg["i5_inner_anchor"]))
        three_head = max(0, len(cfg["pcr1_overhang_rev"]) - len(cfg["i7_inner_anchor"]))
        if five.anchor_end >= 0:
            call.insert_start = five.anchor_end + five_tail
        if three.anchor_pos >= 0:
            call.insert_end = tail_from + three.anchor_pos - three_head

        if best is None or resolved > _resolved(best):
            best = call
            best._resolved = resolved  # type: ignore[attr-defined]
    return best  # type: ignore[return-value]


def _resolved(call: ReadCall) -> int:
    return getattr(call, "_resolved", 0)


def _classify(
    read_id: str,
    five: EndCall,
    three: EndCall,
    orient: str,
    pair_to_sample: dict[tuple[str, str], str],
    accept_single_end: bool,
) -> ReadCall:
    i5, i7 = five.barcode, three.barcode

    if i5 and i7:
        sample = pair_to_sample.get((i5, i7))
        if sample:
            return ReadCall(read_id, "A", sample, i5, i7, orient)
        # both ends read cleanly but the pair is not one that was used:
        # that is a hopped or chimeric molecule, and the combination itself
        # is the evidence
        return ReadCall(read_id, "C", None, i5, i7, orient)

    if (i5 or i7) and accept_single_end:
        # one end only. Unique resolution is required: if several samples
        # share that index, a single end cannot pick between them.
        matches = [
            smp
            for (a, b), smp in pair_to_sample.items()
            if (i5 and a == i5 and not i7) or (i7 and b == i7 and not i5)
        ]
        if len(set(matches)) == 1:
            return ReadCall(read_id, "B", matches[0], i5, i7, orient)
        return ReadCall(read_id, "D", None, i5, i7, orient)

    if i5 or i7:
        return ReadCall(read_id, "B", None, i5, i7, orient)

    # copies disagreed at one end and nothing resolved at the other
    if five.copies_seen > 1 or three.copies_seen > 1:
        return ReadCall(read_id, "C", None, None, None, orient)
    return ReadCall(read_id, "D", None, None, None, orient)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def write_index_matrix(calls, samples, path: str) -> dict:
    """The full observed i7 x i5 grid, assigned combinations and not.

    Every read carries a pair of indices, and only some pairs belong to a
    sample. Reads landing on a combination no sample claims cannot be
    anything but hopping or misreading - there is no library with that
    pair - so those cells measure the run's cross-talk directly rather
    than by inference. The assigned cells are here too, because the rate
    only means something against the depth it is a fraction of.
    """
    assigned = {(s["index_i5"], s["index_i7"]) for s in samples}
    i5s = sorted({s["index_i5"] for s in samples})
    i7s = sorted({s["index_i7"] for s in samples})
    grid: dict = {}
    for c in calls:
        if c.i5 and c.i7:
            grid[(c.i5, c.i7)] = grid.get((c.i5, c.i7), 0) + 1

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["index_i5", "index_i7", "reads", "assigned"])
        for a in i5s:
            for b in i7s:
                w.writerow([a, b, grid.get((a, b), 0),
                            "yes" if (a, b) in assigned else "no"])

    sentinel = [(a, b) for a in i5s for b in i7s if (a, b) not in assigned]
    on_sentinels = sum(grid.get(k, 0) for k in sentinel)
    on_assigned = sum(v for k, v in grid.items() if k in assigned)
    return {
        "n_i5": len(i5s), "n_i7": len(i7s),
        "combinations_possible": len(i5s) * len(i7s),
        "combinations_assigned": len(assigned),
        "sentinel_combinations": len(sentinel),
        "reads_on_sentinels": on_sentinels,
        "reads_on_assigned": on_assigned,
        # The hopping rate this run can actually see. Hopping onto an
        # assigned pair is invisible, so the true rate is higher by the
        # reciprocal of the observable fraction.
        "sentinel_rate": (on_sentinels / (on_sentinels + on_assigned)
                          if (on_sentinels + on_assigned) else 0.0),
        "observable_fraction": (len(sentinel) / (len(i5s) * len(i7s))
                                if i5s and i7s else 0.0),
        "note": (
            "Reads on a combination no sample claims can only be index "
            "cross-talk. Hopping onto an assigned pair is indistinguishable "
            "from a real read, so the observed rate is a floor: divide by "
            "the observable fraction for an estimate of the total."
        ),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Demultiplex nanopore reads carrying in-line Illumina indices."
    )
    ap.add_argument(
        "--reads", required=True, nargs="+",
        help=(
            "pooled reads: a FASTQ, a glob, or a MinKNOW run directory. A "
            "directory is walked, fastq_pass is used and fastq_fail is "
            "excluded and reported"
        ),
    )
    ap.add_argument(
        "--include-failed", action="store_true",
        help="also read fastq_fail - those reads failed basecalling QC and "
             "carry a different error profile",
    )
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--barcode-repeats", type=int, default=2)
    ap.add_argument(
        "--max-mismatch",
        type=int,
        default=1,
        help="requested substitutions per barcode copy; capped by what the "
        "index set can actually support",
    )
    ap.add_argument(
        "--search-window",
        type=int,
        default=120,
        help="how far from each read end to look for the anchor",
    )
    ap.add_argument("--accept-single-end", action="store_true",
                    help="fold tier B into the output; measure the cost first")
    ap.add_argument("--outdir", default="demux")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    cfg = load_assay(args.assay)

    samples = [
        s for s in read_table(args.samplesheet) if s["run_id"] == args.run_id
    ]
    if not samples:
        print(f"no samples for run '{args.run_id}' in the samplesheet", file=sys.stderr)
        return 1

    i5_set = BarcodeSet.build({s["index_i5"].upper() for s in samples}, args.max_mismatch)
    i7_set = BarcodeSet.build({s["index_i7"].upper() for s in samples}, args.max_mismatch)
    pair_to_sample = {
        (s["index_i5"].upper(), s["index_i7"].upper()): s["sample_id"] for s in samples
    }

    print(
        f"{len(samples)} samples | i5 set n={len(i5_set.by_seq)} "
        f"tolerance {i5_set.max_mismatch} | i7 set n={len(i7_set.by_seq)} "
        f"tolerance {i7_set.max_mismatch} | repeats {args.barcode_repeats}"
    )
    if i5_set.max_mismatch < args.max_mismatch or i7_set.max_mismatch < args.max_mismatch:
        print(
            f"  tolerance capped below the requested {args.max_mismatch}: the "
            "index sets in this run are too close to support it"
        )

    # A MinION run is a tree of batched FASTQ, not a file. Resolving it
    # here, and printing what was resolved, keeps the choice of which
    # reads entered the analysis visible rather than buried in a glob.
    read_files, file_report = resolve_read_files(
        args.reads, prefer_pass=not args.include_failed
    )
    if not read_files:
        raise SystemExit(
            f"no FASTQ found under {', '.join(args.reads)} - expected "
            ".fastq/.fq, optionally gzipped"
        )
    print(f"  {file_report['n_files']} FASTQ file(s)"
          + (f", {file_report['from_fastq_pass']} from fastq_pass"
             if file_report["from_fastq_pass"] else ""))
    for note in file_report["notes"]:
        print(f"  note: {note}")

    calls: list[ReadCall] = []
    per_sample: dict[str, list[tuple[str, str, str]]] = defaultdict(list)

    for read_id, seq, qual in open_reads_many(read_files):
        call = assign_read(
            read_id, seq, cfg, i5_set, i7_set, pair_to_sample,
            args.barcode_repeats, args.search_window, args.accept_single_end,
        )
        calls.append(call)
        if call.sample_id:
            s = seq if call.orientation == "+" else revcomp(seq)
            q = qual if call.orientation == "+" else qual[::-1]
            if 0 <= call.insert_start < call.insert_end <= len(s):
                s, q = s[call.insert_start : call.insert_end], q[call.insert_start : call.insert_end]
            per_sample[call.sample_id].append((read_id, s, q))

    # -- write per-sample reads -------------------------------------------
    reads_dir = os.path.join(args.outdir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    safe = lambda x: "".join(c if c.isalnum() or c in "._-" else "_" for c in x)
    for sample_id, reads in per_sample.items():
        with gzip.open(os.path.join(reads_dir, f"{safe(sample_id)}.fastq.gz"), "wt") as fh:
            for rid, s, q in reads:
                fh.write(f"@{rid}\n{s}\n+\n{q}\n")

    # -- per-read assignment table ----------------------------------------
    with open(os.path.join(args.outdir, "assignments.tsv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["read_id", "tier", "sample_id", "index_i5", "index_i7",
                    "orientation"])
        for c in calls:
            w.writerow([c.read_id, c.tier, c.sample_id or "", c.i5 or "",
                        c.i7 or "", c.orientation])

    # -- summary ------------------------------------------------------------
    tiers = Counter(c.tier for c in calls)
    assigned = sum(1 for c in calls if c.sample_id)
    # Every read lands in exactly one tier, and the tiers sum to the input.
    # Asserted rather than reported: if reads have vanished, the results may
    # be wrong for reasons unrelated to biology.
    total = len(calls)
    accounted = sum(tiers.values())

    summary = {
        "run_id": args.run_id,
        "total_reads": total,
        "assigned_reads": assigned,
        "by_tier": {t: tiers.get(t, 0) for t in TIERS},
        "tier_meaning": TIER_MEANING,
        "orientation": dict(Counter(c.orientation for c in calls)),
        "index_matrix": write_index_matrix(
            calls, samples, os.path.join(args.outdir, "index_matrix.tsv")),
        "input_files": file_report,
        "settings": {
            "barcode_repeats": args.barcode_repeats,
            "max_mismatch_requested": args.max_mismatch,
            "max_mismatch_i5": i5_set.max_mismatch,
            "max_mismatch_i7": i7_set.max_mismatch,
            "allow_indels": False,
            "require_copy_agreement": args.barcode_repeats >= 2,
            "accept_single_end": args.accept_single_end,
        },
        "read_accounting_balanced": accounted == total,
        "per_sample": {k: len(v) for k, v in sorted(per_sample.items())},
    }
    with open(os.path.join(args.outdir, "demux_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n{total:,} reads")
    for t in TIERS:
        n = tiers.get(t, 0)
        print(f"  tier {t}  {n:>8,}  {100*n/max(1,total):5.2f}%  {TIER_MEANING[t]}")
    print(f"  assigned {assigned:,} ({100*assigned/max(1,total):.2f}%) "
          f"to {len(per_sample)} samples")

    if accounted != total:
        print(
            f"\nread accounting does not balance: {accounted:,} of {total:,}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
