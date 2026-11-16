#!/usr/bin/env python3
"""
PacBio demultiplexing for amphid's database mode: the pacbio platform
adapter (docs/CONTRACTS.md, layer 1).

A PacBio amplicon read - CCS or HiFi, one read per molecule - is

    [forward tag] forward primer ... insert ... rc(reverse primer) [rc(reverse tag)]

in either orientation. The tags come from the samplesheet: index_i5 is
the forward tag, as read 5'->3' just before the forward primer, and
index_i7 the reverse tag, as read on the reverse primer's strand just
before that primer. A run is one library (one pool): its tags are the
run's samplesheet rows and no other, because a SMRT cell shared between
projects carries other libraries' tags too, and matching against every
pool's tags would hand those reads to the wrong samples (S-BSST527: a
quarter of each file was tags another pool's list declares).

Per read:

  orient   the best forward primer (<= --primer-max-edits) within the
           first --primer-window bases; if there is none the read is
           turned round and searched again. Then the reverse primer,
           reverse-complemented, within the last --primer-window bases,
           after the forward one. Either missing: no primers.
  tags     the up-to --tag-flank bases outside each primer are matched to
           the run's tags at that end by edit distance, the tag's end
           held to the primer and extra flank bases before it free; where
           the read itself ends inside the flank, up to three missing tag
           bases are free too. A match is <= --tag-max-edits and strictly
           better than every other tag at that end; a tie is ambiguous.
  call     both ends matched to a declared (forward, reverse) pair:
           assigned. Both matched, pair declared by no sample: chimeric
           (or a misread tag). One end only, or none: as tiers B and D.
  write    the read from the forward primer's start to the reverse
           primer's end, tags off, forward primer first, to
           <outdir>/reads/<sample>.fastq.gz. The primers stay on: the
           protocol's own trim (ampseqlib.iter_trimmed) cuts them, as it
           does for nanopore and lima output.

The tiers are nanopore's (demux_ont.py), so bin/demux_accounting.py
folds this summary into contract A unchanged: A assigned, B one end
only, C chimeric, D unmatched. There is no tier E: CCS removes the
SMRTbell adapters, so a read is one molecule.

This is a port of nematode-its2-refdb's work/pacbio_test/demux_pacbio.py
(S-BSST527, Halvarsson & Hoglund 2021), made fast: an exact primer is
found by str.find and only an inexact one by Myers' bit-parallel
alignment (ampseqlib.find_anchor_cost), and each distinct tag flank -
most flanks are a tag exactly - is scored once per process.

  --write-lima-barcodes FASTA   writes the run's tags as a lima barcode
                                FASTA and the pair table lima's output is
                                read back through, then stops;
  --from-lima DIR               reads lima's output (--split-named, FASTQ
                                or BAM, per input file) instead of reads,
                                and writes the same outputs. lima runs on
                                Linux only (no macOS build).

Reads: FASTQ (gzipped or not) or unaligned BAM (PacBio's hifi_reads.bam),
as files or directories of them, read as one library.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import struct
import sys
import time
from collections import Counter, defaultdict, deque
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    FASTQ_SUFFIXES, find_anchor_cost, is_plain_acgt, load_primers,
    open_reads, read_table, revcomp, sanitise_id,
)
from demux_ont import write_index_matrix  # noqa: E402

TIERS = ("A", "B", "C", "D")
TIER_MEANING = {
    "A": "both tags match a samplesheet pair",
    "B": "one tag only - the other end unreadable or not this run's",
    "C": "tags name a pair no sample has - chimeric, or a misread tag",
    "D": "no confident match - no primers, or neither tag this run's",
}
# Why a read was not assigned, finer than the tiers; they sum to the reads.
REASONS = {
    "assigned": "A",
    "one_end_only": "B",
    "pair_not_declared": "C",
    "no_primers": "D",
    "no_tag": "D",
}
REASON_MEANING = {
    "no_primers": "the forward or reverse primer was not found near its "
                  "end: a dimer, an off-target product or a truncated read",
    "no_tag": "both primers found, neither flank this run's tag: most "
              "often another library's tags on a shared SMRT cell",
    "one_end_only": "one flank matched a tag, the other matched none or "
                    "two equally",
    "pair_not_declared": "both flanks matched, but no sample has that "
                         "forward and reverse tag together",
}
END_MEANING = {
    "matched": "a tag within --tag-max-edits, strictly closest",
    "ambiguous": "two tags equally close within --tag-max-edits",
    "no_match": "no tag within --tag-max-edits",
}

# The construct's tolerances. An assay JSON may set any of them in a
# "pacbio" section; a command-line value wins over both. 40/3/2/12 are
# what the port was validated with on S-BSST527 (RSII CCS, NC1/NC2,
# 10-base tags two edits apart at their closest within a pool).
CONSTRUCT = {"primer_window": 40, "primer_max_edits": 3,
             "tag_max_edits": 2, "tag_flank": 12}
CLIPPED_FREE = 3      # tag bases free where the read ends inside the flank
MIN_TAG_BASES = 5     # and the fewest of a tag that may still name it
FLUSH_RECORDS = 2000
BAM_SUFFIXES = (".bam",)


# --------------------------------------------------------------------------
# reads: FASTQ or unaligned BAM
# --------------------------------------------------------------------------


_NIBBLE = str.maketrans("0123456789abcdef", "=ACMGRSVTWYHKDBN")


def open_bam(path: str):
    """Yield (read_id, sequence, quality) from an unaligned BAM.

    BAM is BGZF - concatenated gzip members - so gzip reads it, and the
    records are a fixed layout (SAM specification 4.2). Only the name,
    sequence and quality are read; tags are skipped. A read stored
    reverse-complemented (flag 16) is turned back."""
    with gzip.open(path, "rb") as raw:
        fh = io.BufferedReader(raw, 1 << 20)
        if fh.read(4) != b"BAM\x01":
            raise SystemExit(f"{path}: not a BAM file")
        (l_text,) = struct.unpack("<i", fh.read(4))
        fh.read(l_text)
        (n_ref,) = struct.unpack("<i", fh.read(4))
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", fh.read(4))
            fh.read(l_name + 4)
        while True:
            head = fh.read(4)
            if len(head) < 4:
                return
            (size,) = struct.unpack("<i", head)
            rec = fh.read(size)
            l_name = rec[8]
            n_cigar, flag, l_seq = struct.unpack_from("<HHi", rec, 12)
            p = 32
            name = rec[p:p + l_name - 1].decode()
            p += l_name + 4 * n_cigar
            nb = (l_seq + 1) // 2
            seq = rec[p:p + nb].hex().translate(_NIBBLE)[:l_seq]
            p += nb
            q = rec[p:p + l_seq]
            qual = ("!" * l_seq if q[:1] == b"\xff"
                    else bytes(c + 33 for c in q).decode())
            if flag & 16:
                seq, qual = revcomp(seq), qual[::-1]
            yield name, seq, qual


def read_files(paths: list[str]) -> list[str]:
    """Every FASTQ or BAM named, directories expanded (sorted), in the
    order given; anything else in a directory is ignored."""
    out = []
    for p in paths:
        if os.path.isdir(p):
            out += sorted(os.path.join(p, f) for f in os.listdir(p)
                          if f.lower().endswith(FASTQ_SUFFIXES + BAM_SUFFIXES))
        else:
            out.append(p)
    return out


def open_any(path: str):
    return (open_bam(path) if path.lower().endswith(BAM_SUFFIXES)
            else open_reads(path))


# --------------------------------------------------------------------------
# the run's tags and primers
# --------------------------------------------------------------------------


def run_samples(samplesheet: str, run_id: str) -> list[dict]:
    rows = [r for r in read_table(samplesheet) if r.get("run_id") == run_id]
    for r in rows:
        r["index_i5"] = r.get("index_i5", "").upper()
        r["index_i7"] = r.get("index_i7", "").upper()
    bad = [r["sample_id"] for r in rows
           if not r["index_i5"] or not r["index_i7"]]
    if bad:
        raise SystemExit(
            f"run '{run_id}': {', '.join(bad)} without a forward (index_i5) "
            "and a reverse (index_i7) tag; on PacBio the tags are the "
            "demultiplexing key")
    pairs = Counter((r["index_i5"], r["index_i7"]) for r in rows)
    dup = [f"{a}/{b}" for (a, b), n in pairs.items() if n > 1]
    if dup:
        raise SystemExit(f"run '{run_id}': tag pair(s) {', '.join(dup)} "
                         "declared for more than one sample")
    return rows


def construct(assay: str | None, args) -> dict:
    cfg = dict(CONSTRUCT)
    if assay and os.path.exists(assay):
        with open(assay, encoding="utf-8-sig") as fh:
            sec = json.load(fh).get("pacbio") or {}
        cfg.update({k: int(v) for k, v in sec.items() if k in CONSTRUCT})
    for k in CONSTRUCT:
        v = getattr(args, k, None)
        if v is not None:
            cfg[k] = v
    return cfg


# --------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------


def locate(window: str, primers: list[str], max_edits: int):
    """(start, end, edits) of the best-fitting primer in window, or None.

    Exact copies by str.find; otherwise Myers (find_anchor_cost) for the
    end, and the same search run backwards from that end for the true
    start - an indel inside the primer moves the start, and the tag is
    read from there."""
    best = None
    for p in primers:
        i = window.find(p) if is_plain_acgt(p) else -1
        if i >= 0:
            return i, i + len(p), 0
        hit = find_anchor_cost(window, p, max_edits)
        if hit is None or (best and hit[2] >= best[2]):
            continue
        _s, end, cost = hit
        back = find_anchor_cost(window[:end][::-1], p[::-1], cost)
        best = (end - back[1] if back else hit[0], end, cost)
    return best


def tag_distance(tag: str, flank: str, clipped: bool) -> int:
    """Edits aligning tag's end to flank's end, extra flank bases at its
    start free; where the read ends inside the flank (clipped), up to
    CLIPPED_FREE missing tag bases at the tag's start are free too.
    The reference script's tag_dist, unchanged."""
    a, b = tag[::-1], flank[::-1]
    prev = list(range(len(b) + 1))
    ends = []
    for i in range(1, len(a) + 1):
        cur = [i]
        for j in range(1, len(b) + 1):
            cur.append(min(prev[j - 1] + (a[i - 1] != b[j - 1]),
                           prev[j] + 1, cur[j - 1] + 1))
        prev = cur
        if clipped and i >= len(a) - CLIPPED_FREE:
            ends.append(prev[-1])
    return min([min(prev)] + ends)


class TagEnd:
    """One end's tags. match() is memoised by (flank, clipped): flanks
    repeat - most are a tag exactly - so each is scored once."""

    def __init__(self, tags: list[str], max_edits: int):
        self.tags = tags
        self.max = max_edits
        self.exact = {t: k for k, t in enumerate(tags)}
        self.lengths = sorted({len(t) for t in tags})
        self.uniform = len(self.lengths) == 1
        self.memo: dict = {}

    def match(self, flank: str, clipped: bool):
        """(tag index or None, edits, 'matched'|'ambiguous'|'no_match')."""
        # A flank ending in a tag exactly, every tag one length: 0 edits,
        # and no other tag can also be 0 - it would have to equal the
        # flank's last len(tag) bases, or (clipped) the whole flank.
        if self.uniform and len(flank) >= self.lengths[0]:
            k = self.exact.get(flank[-self.lengths[0]:])
            if k is not None:
                return k, 0, "matched"
        key = (flank, clipped)
        got = self.memo.get(key)
        if got is None and clipped and len(flank) < MIN_TAG_BASES:
            # Too little of a tag left to name it. The reference rule
            # cannot match under five bases of a 10-base tag at two edits
            # anyway; with 8-base tags it matched on three.
            got = (None, self.lengths[-1], "no_match")
        if got is None:
            ds = sorted((tag_distance(t, flank, clipped), k)
                        for k, t in enumerate(self.tags))
            d, k = ds[0]
            if d > self.max:
                got = (None, d, "no_match")
            elif len(ds) > 1 and ds[1][0] == d:
                got = (None, d, "ambiguous")
            else:
                got = (k, d, "matched")
            if len(self.memo) < 200000:
                self.memo[key] = got
        return got


class Matcher:
    """Everything one read needs: primers, the two tag ends, the pairs."""

    def __init__(self, samples: list[dict], fwd: list[str], rev: list[str],
                 cfg: dict):
        self.cfg = cfg
        self.fwd = fwd
        self.rev_rc = [revcomp(p) for p in rev]
        self.i5 = sorted({s["index_i5"] for s in samples})
        self.i7 = sorted({s["index_i7"] for s in samples})
        self.end5 = TagEnd(self.i5, cfg["tag_max_edits"])
        self.end7 = TagEnd(self.i7, cfg["tag_max_edits"])
        self.pair = {(self.i5.index(s["index_i5"]),
                      self.i7.index(s["index_i7"])): s["sample_id"]
                     for s in samples}

    def classify(self, seq: str):
        """(reason, sample, start, end, turned, i5 idx, i7 idx, edits5,
        edits7, why5, why7) for one read; start/end bound what is written,
        on the oriented read."""
        cfg = self.cfg
        w, pm, fl = cfg["primer_window"], cfg["primer_max_edits"], cfg["tag_flank"]
        seq = seq.upper()
        turned = False
        p1 = locate(seq[:w], self.fwd, pm)
        if p1 is None:
            seq = revcomp(seq)
            p1 = locate(seq[:w], self.fwd, pm)
            if p1 is None:
                return ("no_primers",) + (None,) * 10
            turned = True
        t0 = max(0, len(seq) - w)
        p2 = locate(seq[t0:], self.rev_rc, pm)
        if p2 is None or p2[0] + t0 <= p1[1]:
            return ("no_primers",) + (None,) * 10
        s1, e2 = p1[0], p2[1] + t0
        k5, d5, why5 = self.end5.match(seq[max(0, s1 - fl):s1], s1 < fl)
        k7, d7, why7 = self.end7.match(revcomp(seq[e2:e2 + fl]),
                                       len(seq) - e2 < fl)
        if k5 is not None and k7 is not None:
            sample = self.pair.get((k5, k7))
            reason = "assigned" if sample else "pair_not_declared"
        else:
            sample = None
            reason = ("one_end_only" if k5 is not None or k7 is not None
                      else "no_tag")
        return (reason, sample, s1, e2, turned, k5, k7, d5, d7, why5, why7)


# --------------------------------------------------------------------------
# workers
# --------------------------------------------------------------------------

_M: Matcher | None = None


def _init_worker(samples, fwd, rev, cfg):
    global _M
    _M = Matcher(samples, fwd, rev, cfg)


def _classify_batch(seqs: list[str]) -> list[tuple]:
    return [_M.classify(s) for s in seqs]


def batches(files: list[str], size: int = 2000):
    """(names, seqs, quals, file index) in batches, in file order."""
    for fi, path in enumerate(files):
        names, seqs, quals = [], [], []
        for name, seq, qual in open_any(path):
            names.append(name)
            seqs.append(seq)
            quals.append(qual)
            if len(seqs) >= size:
                yield names, seqs, quals, fi
                names, seqs, quals = [], [], []
        if seqs:
            yield names, seqs, quals, fi


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------


class SampleWriter:
    """Per-sample gzipped FASTQ, buffered and appended as gzip members
    (as demux_ont.SampleWriter), each with a zero mtime so the same reads
    give the same bytes."""

    def __init__(self, outdir: str):
        self.dir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.buf: dict = defaultdict(list)
        self.started: set = set()
        self.counts: Counter = Counter()
        self.written = 0

    def add(self, sample: str, rec: str) -> None:
        self.counts[sample] += 1
        self.buf[sample].append(rec)
        if len(self.buf[sample]) >= FLUSH_RECORDS:
            self.flush(sample)

    def flush(self, sample: str) -> None:
        recs = self.buf.get(sample)
        if not recs:
            return
        path = os.path.join(self.dir, f"{sanitise_id(sample)}.fastq.gz")
        with open(path, "ab" if sample in self.started else "wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
                fh.write("".join(recs).encode())
        self.started.add(sample)
        self.written += len(recs)
        recs.clear()

    def close(self) -> None:
        for s in list(self.buf):
            self.flush(s)


def summarise(run_id: str, samples: list[dict], reasons: Counter,
              ends: dict, orient: Counter, grid: Counter, edits: dict,
              writer: SampleWriter, files: list[dict], settings: dict,
              outdir: str, demultiplexed_by: str, native: dict | None = None
              ) -> dict:
    total = sum(reasons.values())
    tiers = Counter()
    for r, n in reasons.items():
        tiers[REASONS[r]] += n
    assigned = tiers.get("A", 0)
    per_sample = {s["sample_id"]: writer.counts.get(s["sample_id"], 0)
                  for s in samples}
    summary = {
        "run_id": run_id,
        "platform": "pacbio",
        "demultiplexed_by": demultiplexed_by,
        "total_reads": total,
        "assigned_reads": assigned,
        "by_tier": {t: tiers.get(t, 0) for t in TIERS},
        "tier_meaning": TIER_MEANING,
        "reasons": {r: reasons.get(r, 0) for r in REASONS},
        "reason_meaning": REASON_MEANING,
        # Per END, for the reads whose primers were found: which end
        # fails, and whether it failed by matching nothing or two tags.
        "forward_end": {k: ends["i5"].get(k, 0) for k in END_MEANING},
        "reverse_end": {k: ends["i7"].get(k, 0) for k in END_MEANING},
        "end_meaning": END_MEANING,
        "orientation": dict(orient),
        # Reads on every forward x reverse tag pair seen together, used or
        # not; the unused cells are tag pairs no sample has (chimeras and
        # misread tags), as on nanopore.
        "index_matrix": write_index_matrix(
            grid, samples, os.path.join(outdir, "index_matrix.tsv")),
        # Mean edits to the tag at each end, over a sample's reads.
        "tag_edits": {s: {"forward": round(edits[s][0] / n, 3),
                          "reverse": round(edits[s][1] / n, 3)}
                      for s, n in sorted(writer.counts.items())
                      if n and s in edits},
        "input_files": files,
        "settings": settings,
        "read_accounting_balanced": sum(
            tiers.get(t, 0) for t in TIERS) == total,
        "reads_written": writer.written,
        "reads_written_balanced": writer.written == assigned,
        "per_sample": per_sample,
    }
    if native is not None:
        summary["lima"] = native
    return summary


def write_reads(outdir: str):
    return SampleWriter(os.path.join(outdir, "reads"))


# --------------------------------------------------------------------------
# lima
# --------------------------------------------------------------------------


def near_symmetric(s: dict) -> bool:
    """A sample whose two tags are one, or differ only at the base next to
    the primer (S-BSST527's layout): lima cannot tell its two ends apart."""
    a, b = s["index_i5"], s["index_i7"]
    return a == b or (len(a) == len(b) and a[:-1] == b[:-1])


def lima_design(samples: list[dict]) -> dict:
    """The run's tags as lima sees them: barcodes [(name, sequence)],
    {(name, name): sample} (the pair sorted, as lima orders a pair by its
    barcodes' order in the FASTA), the --hifi-preset, and each barcode's
    (forward, reverse) tags for the matrix.

    Symmetric: every sample's two tags are one, or differ only at the
    primer-adjacent base. A tag pair one base apart is one barcode to
    lima - ASYMMETRIC sent most such reads to a same-barcode pair and
    dropped them (S-BSST527 on the farm: 15 of 59 samples lost) - so each
    sample gets one barcode, its tag less that base where any sample's
    two differ there, and lima runs SYMMETRIC. Asymmetric: no sample's
    tags are near; one barcode per distinct tag, ASYMMETRIC. A run mixing
    the two has no lima preset that keeps both: refused."""
    near = [near_symmetric(s) for s in samples]
    if all(near):
        cut = any(s["index_i5"] != s["index_i7"] for s in samples)
        seq = {s["sample_id"]: s["index_i5"][:-1] if cut else s["index_i5"]
               for s in samples}
        dup = [q for q, n in Counter(seq.values()).items() if n > 1]
        if dup:
            raise SystemExit(
                f"--pacbio_demux lima: samples share the barcode "
                f"{', '.join(dup)} once each pair's primer-adjacent base is "
                "set aside, so lima cannot tell them apart; use "
                "--pacbio_demux builtin")
        order = sorted(set(seq.values()))
        name = {q: f"t{k:03d}" for k, q in enumerate(order, 1)}
        return {"barcodes": [(name[q], q) for q in order],
                "pairs": {(name[seq[s["sample_id"]]],) * 2: s["sample_id"]
                          for s in samples},
                "preset": "SYMMETRIC",
                "ends": {name[seq[s["sample_id"]]]:
                         (s["index_i5"], s["index_i7"]) for s in samples}}
    if any(near):
        both = sorted(s["sample_id"] for s, n in zip(samples, near) if n)
        raise SystemExit(
            f"--pacbio_demux lima: {len(both)} sample(s) carry one tag at "
            f"both ends ({', '.join(both[:5])}{'...' if len(both) > 5 else ''})"
            " and the rest two different tags; lima's SYMMETRIC preset drops "
            "the second kind and ASYMMETRIC the first. Use --pacbio_demux "
            "builtin")
    order = sorted({s["index_i5"] for s in samples}
                   | {s["index_i7"] for s in samples})
    name = {q: f"t{k:03d}" for k, q in enumerate(order, 1)}
    return {"barcodes": [(name[q], q) for q in order],
            "pairs": {tuple(sorted((name[s["index_i5"]],
                                    name[s["index_i7"]]))): s["sample_id"]
                      for s in samples},
            "preset": "ASYMMETRIC",
            "ends": {name[q]: (q, q) for q in order}}


LIMA_OUT = re.compile(r"^(?P<prefix>.+)\.(?P<a>[^.]+)--(?P<b>[^.]+)"
                      r"\.(?P<ext>fastq|fq|bam)(?:\.gz)?$")


def parse_lima_summary(path: str) -> dict:
    """lima's .lima.summary as {label: count}: 'ZMWs input (A) : 1000'."""
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r"^\s*([^:(]+?)\s*(?:\([A-Z]\))?\s*:\s*(\d+)", line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


def from_lima(lima_dir: str, samples: list[dict], writer: SampleWriter):
    """Read lima's per-pair outputs back as the samples they belong to.

    Returns (reasons, grid, files, native): the fates as this script's
    reasons - a declared pair assigned, an undeclared one
    pair_not_declared, reads lima held below its end-score threshold
    one_end_only (one barcode end too weak to call), and every other read
    lima did not output no_tag - and lima's own summaries verbatim."""
    design = lima_design(samples)
    pairs, by_name = design["pairs"], design["ends"]
    reasons, grid = Counter(), Counter()
    files, summ = [], {}
    names = sorted(os.listdir(lima_dir))
    for f in names:
        if f.endswith(".lima.summary"):
            for k, v in parse_lima_summary(os.path.join(lima_dir, f)).items():
                summ[k] = summ.get(k, 0) + v
    for f in names:
        m = LIMA_OUT.match(f)
        if not m or m.group("a") not in by_name or m.group("b") not in by_name:
            continue
        a, b = m.group("a"), m.group("b")
        sample = pairs.get(tuple(sorted((a, b))))
        n = 0
        for name, seq, qual in open_any(os.path.join(lima_dir, f)):
            n += 1
            if sample:
                writer.add(sample, f"@{name}\n{seq}\n+\n{qual}\n")
        reasons["assigned" if sample else "pair_not_declared"] += n
        # The pair as forward and reverse tag: where the samplesheet says
        # which end a tag belongs to, so; a symmetric barcode stands for
        # its sample's forward tag at one end and reverse at the other.
        i5s = {s["index_i5"] for s in samples}
        qa, qb = by_name[a][0], by_name[b][1]
        if design["preset"] == "ASYMMETRIC" and qa not in i5s:
            qa, qb = by_name[b][0], by_name[a][1]
        grid[(qa, qb)] += n
        files.append({"path": f, "format": m.group("ext"), "reads": n})
    total = next((v for k, v in summ.items() if "input" in k.lower()), None)
    out = sum(reasons.values())
    if total is None:
        total = out
    weak = next((v for k, v in summ.items()
                 if "end score" in k.lower()), 0)
    weak = min(weak, max(0, total - out))
    reasons["one_end_only"] += weak
    reasons["no_tag"] += max(0, total - out - weak)
    native = {"summary": summ,
              "fate_mapping": "declared pair -> assigned; undeclared pair "
                              "-> chimeric; below lima's min end score -> "
                              "one end only; every other read lima did not "
                              "output -> unmatched"}
    return reasons, grid, files, native


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", nargs="+",
                    help="the library's FASTQ or BAM files, or directories")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--primers", required=True,
                    help="primers TSV (name, direction, sequence)")
    ap.add_argument("--assay", help="assay JSON; its optional 'pacbio' "
                    "section sets the construct's tolerances")
    ap.add_argument("--outdir", default="demux")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--primer-window", dest="primer_window", type=int,
                    help=f"bases from each end searched for its primer "
                    f"(default {CONSTRUCT['primer_window']})")
    ap.add_argument("--primer-max-edits", dest="primer_max_edits", type=int,
                    help=f"default {CONSTRUCT['primer_max_edits']}")
    ap.add_argument("--tag-max-edits", dest="tag_max_edits", type=int,
                    help=f"default {CONSTRUCT['tag_max_edits']}")
    ap.add_argument("--tag-flank", dest="tag_flank", type=int,
                    help="bases outside each primer read for its tag "
                    f"(default {CONSTRUCT['tag_flank']})")
    ap.add_argument("--write-lima-barcodes", metavar="FASTA",
                    help="write the run's tags as lima barcodes (and "
                    "FASTA.preset, the --hifi-preset), then stop")
    ap.add_argument("--from-lima", metavar="DIR",
                    help="read lima's --split-named output instead of reads")
    args = ap.parse_args(argv)
    t0 = time.time()

    samples = run_samples(args.samplesheet, args.run_id)
    if not samples:
        raise SystemExit(f"no samples for run '{args.run_id}' in the "
                         "samplesheet")
    if args.write_lima_barcodes:
        design = lima_design(samples)
        with open(args.write_lima_barcodes, "w", encoding="utf-8") as fh:
            fh.writelines(f">{n}\n{q}\n" for n, q in design["barcodes"])
        with open(args.write_lima_barcodes + ".preset", "w",
                  encoding="utf-8") as fh:
            fh.write(design["preset"] + "\n")
        return 0
    cfg = construct(args.assay, args)
    primers = load_primers(args.primers)
    fwd = [p for _n, p in primers.fwd]
    rev = [p for _n, p in primers.rev]
    if not fwd or not rev:
        raise SystemExit(f"{args.primers}: needs forward and reverse primers")
    os.makedirs(args.outdir, exist_ok=True)
    writer = write_reads(args.outdir)
    settings = {**cfg, "tags_forward": len({s["index_i5"] for s in samples}),
                "tags_reverse": len({s["index_i7"] for s in samples})}

    if args.from_lima:
        reasons, grid, files, native = from_lima(args.from_lima, samples,
                                                 writer)
        writer.close()
        summary = summarise(args.run_id, samples, reasons,
                            {"i5": {}, "i7": {}}, Counter(), grid, {},
                            writer, files,
                            {"preset": lima_design(samples)["preset"]},
                            args.outdir, "lima, read back by demux_pacbio.py",
                            native)
    else:
        if not args.reads:
            raise SystemExit("--reads is required (or --from-lima)")
        paths = read_files(args.reads)
        if not paths:
            raise SystemExit(f"no FASTQ or BAM under {' '.join(args.reads)}")
        m = Matcher(samples, fwd, rev, cfg)
        reasons, orient, grid = Counter(), Counter(), Counter()
        ends = {"i5": Counter(), "i7": Counter()}
        edits: dict = defaultdict(lambda: [0, 0])
        per_file = Counter()
        pool = (Pool(args.threads, _init_worker, (samples, fwd, rev, cfg))
                if args.threads > 1 else None)
        try:
            # imap keeps the batches' order, so each batch's calls pair
            # with the batch held back here; only sequences cross to the
            # workers.
            held = deque()

            def feed():
                for b in batches(paths):
                    held.append(b)
                    yield b[1]
            calls_of = (pool.imap(_classify_batch, feed(), chunksize=1)
                        if pool else
                        ([m.classify(q) for q in b] for b in feed()))
            for calls in calls_of:
                names, seqs, quals, fi = held.popleft()
                per_file[fi] += len(seqs)
                for name, seq, qual, c in zip(names, seqs, quals, calls):
                    (reason, sample, s1, e2, turned, k5, k7, d5, d7,
                     why5, why7) = c
                    reasons[reason] += 1
                    if reason == "no_primers":
                        continue
                    orient["reverse" if turned else "forward"] += 1
                    ends["i5"][why5] += 1
                    ends["i7"][why7] += 1
                    if k5 is not None and k7 is not None:
                        grid[(m.i5[k5], m.i7[k7])] += 1
                    if not sample:
                        continue
                    seq = seq.upper()
                    if turned:
                        seq, qual = revcomp(seq), qual[::-1]
                    writer.add(sample, f"@{name}\n{seq[s1:e2]}\n+\n"
                                       f"{qual[s1:e2]}\n")
                    edits[sample][0] += d5
                    edits[sample][1] += d7
        finally:
            if pool:
                pool.close()
                pool.join()
        writer.close()
        files = [{"path": os.path.basename(p),
                  "format": "bam" if p.lower().endswith(BAM_SUFFIXES)
                  else "fastq", "reads": per_file.get(k, 0)}
                 for k, p in enumerate(paths)]
        summary = summarise(args.run_id, samples, reasons, ends, orient,
                            grid, edits, writer, files, settings,
                            args.outdir, "amphid (demux_pacbio.py)")

    with open(os.path.join(args.outdir, "demux_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    total = summary["total_reads"]
    print(f"{args.run_id}: {total:,} reads, {summary['assigned_reads']:,} "
          f"assigned ({100 * summary['assigned_reads'] / max(1, total):.1f}%) "
          f"to {sum(1 for v in summary['per_sample'].values() if v)} of "
          f"{len(samples)} samples, {time.time() - t0:.1f} s")
    for r, n in summary["reasons"].items():
        print(f"  {r:<18} {n:>9,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
