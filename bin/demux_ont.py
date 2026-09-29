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
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    DEFAULTS,
    open_reads_many,
    resolve_read_files,
    find_anchor,
    find_anchor_cost,
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
TIERS = ("A", "B", "C", "D", "E")

TIER_MEANING = {
    "A": "both ends agree and match a samplesheet pair",
    "B": "one end only - the other truncated or unreadable",
    "C": "ends disagree - chimeric",
    "D": "no confident match",
    "E": "internal adapter - the read is more than one molecule",
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
    # Why this end did or did not yield a barcode. The tier a read lands in
    # says how much evidence survived; this says what destroyed it, which is
    # the part you need to fix a run where most reads fail. Costs a string
    # assignment per end and is read by bin/diagnose_demux.py.
    reason: str = "ok"
    anchor_edits: int = -1


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
        out.reason = "no_anchor_configured"
        return out

    # The anchor is searched only within a window from the read end it
    # belongs to. The construct puts it at a known offset, and scanning the
    # whole read would invite a chance match elsewhere.
    window = read[:search_window]
    exact = window.find(anchor)
    if exact >= 0:
        a_start, a_end = exact, exact + len(anchor)
        out.anchor_edits = 0
    else:
        hit = find_anchor_cost(window, anchor, anchor_max_edits)
        if hit is None:
            out.reason = "anchor_not_found"
            return out
        a_start, a_end, out.anchor_edits = hit
    out.anchor_pos = a_start
    out.anchor_end = a_end

    span = bc.length * repeats
    # Measured from the alignment boundary, so an indel inside the anchor
    # does not displace the barcode window.
    start = a_start - span if barcode_before_anchor else a_end
    if start + slack < 0 or start + span - slack > len(read):
        # the anchor was found but the barcode runs off the end of the read
        out.reason = "barcode_off_read"
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
        out.reason = "no_copy_matched"
        return out

    if len(set(seen)) == 1:
        out.barcode = seen[0]
        out.copies_agreed = len(seen) == repeats
        out.reason = "ok" if out.copies_agreed else "ok_single_copy"
    else:
        # copies disagree: a frame slip from an upstream indel produces
        # exactly this, and guessing which copy is right would defeat the
        # point of duplicating them
        out.barcode = None
        out.copies_agreed = False
        out.reason = "copies_disagreed"
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


def has_internal_adapter(
    s: str, start: int, end: int, anchors: tuple[str, ...], max_edits: int
) -> bool:
    """Is there another molecule's end between this read's two ends?

    A read is assumed to be one molecule, and the two indices found at its
    ends are assumed to belong to the same one. Ligation does not respect
    that assumption. Two amplicons joined head to tail present the first
    one's i5 at the 5' end and the SECOND one's i7 at the 3' - a pair of
    real, declared indices that very probably belongs to some well on the
    plate, just not to either molecule. It demultiplexes cleanly into the
    wrong sample and there is nothing downstream that could notice.

    The junction is the evidence. A clean amplicon's interior is insert;
    it contains no adapter. So finding an anchor between the two outer
    anchors means the read spans more than one molecule and its ends
    cannot be trusted as a pair.

    Both strands are searched, because the join may be head to head.
    """
    if start < 0 or end < 0 or end <= start:
        return False
    interior = s[start:end]
    return any(find_anchor(interior, p, max_edits) is not None for p in anchors)


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
    junction_anchors: tuple[str, ...] = (),
    anchor_max_edits: int = 2,
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

        # Only a read that resolved both ends is at risk of a fabricated
        # pair, and only those pay for the interior scan. A read that
        # already failed has nothing to protect.
        if junction_anchors and call.tier == "A":
            if has_internal_adapter(s, call.insert_start, call.insert_end,
                                    junction_anchors, anchor_max_edits):
                call.tier = "E"
                call.sample_id = None

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


def write_index_matrix(grid, samples, path: str) -> dict:
    """The full observed i7 x i5 grid, assigned combinations and not.

    Every read carries a pair of indices, and only some pairs belong to a
    sample. Reads landing on a combination no sample claims cannot be
    anything but hopping or misreading - there is no library with that
    pair - so those cells measure the run's cross-talk directly rather
    than by inference. The assigned cells are here too, because the rate
    only means something against the depth it is a fraction of.

    `grid` is a {(i5, i7): reads} tally accumulated as reads stream past,
    not a list of calls: the grid is bounded by the index sets, a list of
    calls is bounded only by the size of the run.
    """
    assigned = {(s["index_i5"], s["index_i7"]) for s in samples}
    i5s = sorted({s["index_i5"] for s in samples})
    i7s = sorted({s["index_i7"] for s in samples})

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


def safe_name(x: str) -> str:
    """A sample id is free text; a filename is not."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in x)


class SampleWriter:
    """Per-sample FASTQ output that holds neither the run nor 384 handles.

    Reads arrive interleaved, one sample at a time being impossible, so
    there are only two obvious ways to write them out and both scale
    badly. Holding every read until the end makes memory grow with the
    run - about 1.4 KB per read, so a 20M-read flow cell wants ~28 GB and
    gets killed after an hour of work. Keeping one open gzip stream per
    sample instead costs a compressor state and a file descriptor per
    sample, on a plate that may have 384 of them.

    So each sample gets a small buffer which is appended to its file and
    dropped. gzip members concatenate, and every gzip reader handles a
    multi-member file transparently, so a file written in a thousand
    flushes reads back as one stream. Memory is then bounded by
    (samples x flush_records) rather than by the size of the run, and at
    most one file is open at a time.
    """

    def __init__(self, outdir: str, flush_records: int = 2000):
        self.dir = outdir
        os.makedirs(self.dir, exist_ok=True)
        self.flush_records = flush_records
        self._buf: dict[str, list[str]] = defaultdict(list)
        self._started: set[str] = set()
        self.counts: Counter = Counter()
        self.written = 0

    def path_for(self, sample_id: str) -> str:
        return os.path.join(self.dir, f"{safe_name(sample_id)}.fastq.gz")

    def add(self, sample_id: str, read_id: str, seq: str, qual: str) -> None:
        self._buf[sample_id].append(f"@{read_id}\n{seq}\n+\n{qual}\n")
        self.counts[sample_id] += 1
        if len(self._buf[sample_id]) >= self.flush_records:
            self._flush(sample_id)

    def _flush(self, sample_id: str) -> None:
        recs = self._buf.get(sample_id)
        if not recs:
            return
        # 'wt' the first time so a re-run truncates rather than appends to
        # the output of the run before it; 'at' after that.
        mode = "at" if sample_id in self._started else "wt"
        with gzip.open(self.path_for(sample_id), mode) as fh:
            fh.write("".join(recs))
        self._started.add(sample_id)
        self.written += len(recs)
        recs.clear()

    def close(self) -> None:
        for sample_id in list(self._buf):
            self._flush(sample_id)


def _build_sets(samplesheet: str, run_id: str, max_mismatch: int):
    """Rebuild the index sets from files. Workers get paths, not objects:
    what the workers agree on is then the input, not a pickled snapshot of
    one process's interpretation of it."""
    samples = [s for s in read_table(samplesheet) if s["run_id"] == run_id]
    i5_set = BarcodeSet.build({s["index_i5"].upper() for s in samples},
                              max_mismatch)
    i7_set = BarcodeSet.build({s["index_i7"].upper() for s in samples},
                              max_mismatch)
    pair_to_sample = {
        (s["index_i5"].upper(), s["index_i7"].upper()): s["sample_id"]
        for s in samples
    }
    return samples, i5_set, i7_set, pair_to_sample


def demux_one_file(job: dict) -> dict:
    """Demultiplex one input file into one part directory.

    Every output this produces is per-part: its own per-sample FASTQ, its
    own slice of the per-read table, its own tallies. Nothing is shared
    between workers and nothing is appended to in common, so the parts can
    be produced in any order, on any number of cores, and merged after.

    The merge is cheap because of how the pieces are shaped. gzip members
    concatenate, so per-sample reads are joined by copying bytes rather
    than by decompressing and recompressing tens of gigabytes; the tables
    are appended; the tallies are summed.
    """
    cfg = load_assay(job["assay"])
    _samples, i5_set, i7_set, pair_to_sample = _build_sets(
        job["samplesheet"], job["run_id"], job["max_mismatch"])

    junction_anchors: tuple[str, ...] = ()
    if job["chimera_check"]:
        a5, a7 = cfg["i5_inner_anchor"], cfg["i7_inner_anchor"]
        junction_anchors = (a5, a7, revcomp(a5), revcomp(a7))

    part_reads = job["part_reads"]
    if os.path.isdir(part_reads):
        for stale in os.listdir(part_reads):
            if stale.endswith(".fastq.gz"):
                os.remove(os.path.join(part_reads, stale))
    writer = SampleWriter(part_reads, flush_records=job["flush_records"])

    tiers: Counter = Counter()
    orients: Counter = Counter()
    grid: Counter = Counter()
    total = assigned = 0

    tmp_tsv = job["assign_tsv"] + ".partial"
    with open(tmp_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for read_id, seq, qual in open_reads(job["src"]):
            call = assign_read(
                read_id, seq, cfg, i5_set, i7_set, pair_to_sample,
                job["barcode_repeats"], job["search_window"],
                job["accept_single_end"],
                junction_anchors=junction_anchors,
                anchor_max_edits=job["anchor_max_edits"],
            )
            total += 1
            tiers[call.tier] += 1
            orients[call.orientation] += 1
            if call.i5 and call.i7:
                grid[(call.i5, call.i7)] += 1
            w.writerow([call.read_id, call.tier, call.sample_id or "",
                        call.i5 or "", call.i7 or "", call.orientation])
            if call.sample_id:
                assigned += 1
                s = seq if call.orientation == "+" else revcomp(seq)
                q = qual if call.orientation == "+" else qual[::-1]
                if 0 <= call.insert_start < call.insert_end <= len(s):
                    s = s[call.insert_start:call.insert_end]
                    q = q[call.insert_start:call.insert_end]
                writer.add(call.sample_id, read_id, s, q)
    writer.close()
    os.replace(tmp_tsv, job["assign_tsv"])

    stats = {
        "source": job["src"],
        "total_reads": total,
        "assigned_reads": assigned,
        "by_tier": dict(tiers),
        "orientation": dict(orients),
        # Tuple keys are not JSON; the pair is joined with a tab, which
        # cannot occur in a barcode.
        "grid": {f"{a}\t{b}": c for (a, b), c in grid.items()},
        "per_sample": dict(writer.counts),
        "reads_written": writer.written,
    }
    with open(job["stats"], "w", encoding="utf-8") as fh:
        json.dump(stats, fh)
    return stats


def merge_sample_reads(part_dirs: list[str], dest_dir: str) -> dict[str, int]:
    """Join each sample's per-part FASTQ into one file, without recompressing.

    A gzip file may be several members back to back, and every reader
    handles that transparently, so concatenating the compressed bytes is
    a valid gzip file of the concatenated content. At this scale that
    matters: the alternative is decompressing and recompressing the whole
    run to achieve nothing.
    """
    os.makedirs(dest_dir, exist_ok=True)
    for stale in os.listdir(dest_dir):
        if stale.endswith(".fastq.gz"):
            os.remove(os.path.join(dest_dir, stale))
    sizes: dict[str, int] = {}
    for part in part_dirs:
        if not os.path.isdir(part):
            continue
        for name in sorted(os.listdir(part)):
            if not name.endswith(".fastq.gz"):
                continue
            src = os.path.join(part, name)
            with open(src, "rb") as fin, \
                    open(os.path.join(dest_dir, name), "ab") as fout:
                while True:
                    chunk = fin.read(1 << 20)
                    if not chunk:
                        break
                    fout.write(chunk)
            sizes[name] = sizes.get(name, 0) + os.path.getsize(src)
    return sizes


def _progress(done: int, total: int, results: list, started: float) -> None:
    reads = sum(r["total_reads"] for r in results)
    dt = max(1e-6, time.time() - started)
    left = (total - done) * (dt / max(1, done))
    sys.stderr.write(
        f"\r  {done}/{total} files  {reads:,} reads  "
        f"{reads / dt:,.0f} reads/s  ~{left / 60:.0f} min left    ")
    sys.stderr.flush()


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
    ap.add_argument(
        "--no-chimera-check", dest="chimera_check", action="store_false",
        help="stop quarantining reads that carry an adapter between their "
             "two ends. Those reads span more than one molecule, so the "
             "index at one end need not belong to the molecule at the "
             "other - they demultiplex into a real sample that never "
             "contained them. Disabling this trades correctness for yield",
    )
    ap.add_argument("--anchor-max-edits", type=int, default=2,
                    help="edit budget for the anchor search, including the "
                         "interior scan for junctions")
    ap.add_argument(
        "--flush-records", type=int, default=2000,
        help="reads buffered per sample before appending to its file. "
             "Peak memory is roughly samples x this x read length, so it "
             "trades memory against the number of gzip members written. "
             "Lower it only to test the flush boundary",
    )
    ap.add_argument("--threads", type=int, default=0,
                    help="worker processes; 0 picks min(8, cores). Work\n                         is divided by input file and resumes per file")
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

    # Nothing below accumulates per-read state. Reads stream past once;
    # the per-read table is written as they go, the per-sample FASTQ are
    # buffered in small chunks, and everything else is a tally whose size
    # is set by the number of samples and indices, not the size of the
    # run. Peak memory is therefore flat in the read count - which is the
    # only property that makes a full flow cell safe to run.
    # Both anchors and both their reverse complements: a head-to-tail join
    # leaves the junction in the same orientation, a head-to-head join
    # leaves it inverted, and both fabricate a pair.
    # One shard per input file, exactly as the splitter divides its work.
    # The parts are independent, so they can run on any number of cores,
    # and a part already finished is skipped - which is what makes a stage
    # this long survive a wall-clock limit.
    parts_dir = os.path.join(args.outdir, "parts")
    os.makedirs(parts_dir, exist_ok=True)

    jobs, results, part_dirs, tsvs = [], [], [], []
    for i, src in enumerate(read_files):
        stem = f"{i:05d}"
        stats_path = os.path.join(parts_dir, stem + ".json")
        tsv_path = os.path.join(parts_dir, stem + ".assignments.tsv")
        part_reads = os.path.join(parts_dir, stem + ".reads")
        part_dirs.append(part_reads)
        tsvs.append(tsv_path)
        if os.path.exists(stats_path) and os.path.exists(tsv_path):
            with open(stats_path, encoding="utf-8") as fh:
                results.append(json.load(fh))
            continue
        jobs.append({
            "src": src, "stats": stats_path, "assign_tsv": tsv_path,
            "part_reads": part_reads, "assay": args.assay,
            "samplesheet": args.samplesheet, "run_id": args.run_id,
            "max_mismatch": args.max_mismatch,
            "barcode_repeats": args.barcode_repeats,
            "search_window": args.search_window,
            "accept_single_end": args.accept_single_end,
            "chimera_check": args.chimera_check,
            "anchor_max_edits": args.anchor_max_edits,
            "flush_records": max(1, args.flush_records),
        })

    threads = args.threads or min(8, os.cpu_count() or 1)
    print(f"  {len(read_files)} file(s), {len(results)} already done, "
          f"{len(jobs)} to do on {threads} thread(s)")

    started = time.time()
    if jobs:
        if threads > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                for k, stats in enumerate(pool.map(demux_one_file, jobs), 1):
                    results.append(stats)
                    _progress(k, len(jobs), results, started)
        else:
            for k, job in enumerate(jobs, 1):
                results.append(demux_one_file(job))
                _progress(k, len(jobs), results, started)
        sys.stderr.write("\n")

    # -- merge --------------------------------------------------------------
    # Parts are merged in input order, so the output is the same whatever
    # order the workers finished in and whatever thread count was used.
    reads_dir = os.path.join(args.outdir, "reads")
    merge_sample_reads(part_dirs, reads_dir)

    assign_path = os.path.join(args.outdir, "assignments.tsv")
    with open(assign_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["read_id", "tier", "sample_id", "index_i5", "index_i7",
                    "orientation"])
        for tsv in tsvs:
            if not os.path.exists(tsv):
                continue
            with open(tsv, encoding="utf-8") as part:
                for line in part:
                    fh.write(line)

    by_source = {r["source"]: r for r in results}
    ordered = [by_source[s] for s in read_files if s in by_source]
    tiers: Counter = Counter()
    orients: Counter = Counter()
    grid: Counter = Counter()
    per_sample: Counter = Counter()
    total = assigned = written = 0
    for r in ordered:
        total += r["total_reads"]
        assigned += r["assigned_reads"]
        written += r["reads_written"]
        tiers.update(r["by_tier"])
        orients.update(r["orientation"])
        per_sample.update(r["per_sample"])
        for key, c in r["grid"].items():
            a, b = key.split("\t")
            grid[(a, b)] += c

    # -- summary ------------------------------------------------------------
    # Every read lands in exactly one of the known tiers, and the tiers sum
    # to the input. Asserted rather than reported: if reads have vanished,
    # the results may be wrong for reasons unrelated to biology. Summing
    # only over TIERS rather than over everything seen also catches a tier
    # that is not one of the four, which a total-vs-total check cannot.
    accounted = sum(tiers.get(t, 0) for t in TIERS)

    summary = {
        "run_id": args.run_id,
        "total_reads": total,
        "assigned_reads": assigned,
        "by_tier": {t: tiers.get(t, 0) for t in TIERS},
        "tier_meaning": TIER_MEANING,
        "orientation": dict(orients),
        "index_matrix": write_index_matrix(
            grid, samples, os.path.join(args.outdir, "index_matrix.tsv")),
        "input_files": file_report,
        "settings": {
            "barcode_repeats": args.barcode_repeats,
            "max_mismatch_requested": args.max_mismatch,
            "max_mismatch_i5": i5_set.max_mismatch,
            "max_mismatch_i7": i7_set.max_mismatch,
            "allow_indels": False,
            "require_copy_agreement": args.barcode_repeats >= 2,
            "accept_single_end": args.accept_single_end,
            "chimera_check": args.chimera_check,
            "anchor_max_edits": args.anchor_max_edits,
        },
        "read_accounting_balanced": accounted == total,
        # A second, independent check: the reads written to per-sample
        # files must equal the reads called to a sample. The tier tally and
        # the FASTQ output are produced by different code paths, so this
        # catches a dropped buffer that the tier check cannot see.
        "reads_written": written,
        "reads_written_balanced": written == assigned,
        "per_sample": dict(sorted(per_sample.items())),
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
    if written != assigned:
        print(
            f"\n{written:,} reads were written to per-sample files but "
            f"{assigned:,} were called to a sample - the FASTQ output is "
            "incomplete and must not be analysed",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
