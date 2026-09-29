#!/usr/bin/env python3
"""
Shared primitives for nf-helminth-amplicon.

Everything in here is used by at least two of the executables, and the
reason it is shared rather than copied is that the alternative has already
bitten this pipeline once: the approximate matcher had a one-character bug
that silently halved assignment, and a second copy of it would have had to
be found and fixed separately.

Standard library only, so the field profile has nothing to install.
"""

from __future__ import annotations

import csv
import glob
import gzip
import json
import os
from dataclasses import dataclass, field

# --------------------------------------------------------------------------
# construct defaults, overridden by the assay config
# --------------------------------------------------------------------------

DEFAULTS = {
    "p5_tail": "AATGATACGGCGACCACCGAGATCTACAC",
    "p7_tail": "CAAGCAGAAGACGGCATACGAGAT",
    "i5_inner_anchor": "TCGTCGGCAGCGTC",
    "i7_inner_anchor": "CCGAGCCCACGAGAC",
    # The anchors are only the DISCRIMINATING part of each overhang - the
    # part that differs between the two ends. Both overhangs also carry the
    # 19bp Nextera mosaic end, which is shared and so useless for telling the
    # ends apart, but is still adapter.
    "pcr1_overhang_fwd": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "pcr1_overhang_rev": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
}


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


# IUPAC ambiguity, as the set of concrete bases each code stands for. Real
# primers are routinely degenerate, and a matcher that treats Y as a literal
# character fails to find the primer in every sequence that carries either
# base it was designed to cover.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


# --------------------------------------------------------------------------
# approximate matching
# --------------------------------------------------------------------------


def find_anchor_cost(
    text: str, pattern: str, max_edits: int
) -> tuple[int, int, int] | None:
    """As find_anchor, but also returns the edit cost.

    Primer trimming needs the cost: several primers compete for the same
    read start, and the one that fits best is the one that was used.
    """
    m = len(pattern)
    if not m or m > 64 or not text:
        return None

    # Bits per concrete base, expanding any degenerate code in the pattern,
    # then one mask per IUPAC code the text might use - so degeneracy on
    # either side matches, and an N in either matches anything.
    base: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i, c in enumerate(pattern):
        for b in IUPAC.get(c, ""):
            base[b] |= 1 << i
    peq = {
        code: (base["A"] if "A" in bases else 0)
        | (base["C"] if "C" in bases else 0)
        | (base["G"] if "G" in bases else 0)
        | (base["T"] if "T" in bases else 0)
        for code, bases in IUPAC.items()
    }

    full = (1 << m) - 1
    top = 1 << (m - 1)
    pv, mv, score = full, 0, m
    best_end, best_cost = -1, max_edits + 1

    for j, c in enumerate(text, 1):
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) ^ pv) | eq) & full
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = (mh | ~(xv | ph)) & full
        mv = ph & xv
        if score < best_cost:
            best_cost, best_end = score, j

    if best_end < 0 or best_cost > max_edits:
        return None
    return max(0, best_end - m), best_end, best_cost


def find_anchor(text: str, pattern: str, max_edits: int) -> tuple[int, int] | None:
    """Locate `pattern` in `text` allowing indels; return (start, end) of the
    best alignment, or None.

    Edit distance rather than Hamming: a fixed-width comparison fails
    outright on an indel, and indels are nanopore's dominant error - roughly
    one read in eight is lost on a 14bp anchor that way.

    Myers' bit-parallel algorithm rather than the obvious dynamic program.
    The DP is O(window x pattern) per call and, at four calls per read across
    a hundred thousand reads, runs to hundreds of millions of Python-level
    operations. Myers is O(window) integer operations for any pattern up to
    64bp, which every anchor and primer here is.

    Returning the END of the alignment matters: an indel inside the pattern
    shifts everything after it, so a downstream offset must be measured from
    where the match actually finished, not from its nominal start plus its
    length.
    """
    hit = find_anchor_cost(text, pattern, max_edits)
    return None if hit is None else (hit[0], hit[1])


# --------------------------------------------------------------------------
# primers
# --------------------------------------------------------------------------


@dataclass
class PrimerSet:
    """The primers, kept as flat lists because a read does not know its locus
    yet - that is what assignment is for. Each read start is tried against
    every primer of the right direction and the best fit wins."""

    fwd: list[tuple[str, str]] = field(default_factory=list)
    rev: list[tuple[str, str]] = field(default_factory=list)
    max_edits: int = 2

    def match(self, seq: str, direction: str) -> tuple[str, int] | None:
        """Best primer at the start of `seq`; returns (name, end offset).

        The search window is the primer length plus a little slack, not the
        whole read: a primer sequence occurring by chance in the middle of an
        amplicon is not evidence of anything, and allowing the match to float
        would let it trim away real template.
        """
        best = None
        for name, p in self.fwd if direction == "fwd" else self.rev:
            window = seq[: len(p) + 6]
            hit = find_anchor_cost(window, p, self.max_edits)
            if hit is None:
                continue
            start, end, cost = hit
            # A primer that does not start at (or very near) the read start
            # is a coincidence, not the primer.
            if start > 4:
                continue
            if best is None or cost < best[2]:
                best = (name, end, cost)
        return None if best is None else (best[0], best[1])


def load_primers(path: str) -> PrimerSet:
    ps = PrimerSet()
    for r in read_table(path):
        d = r.get("direction", "").strip().lower()
        seq = r.get("sequence", "").strip().upper()
        if not seq:
            continue
        if d.startswith("f"):
            ps.fwd.append((r.get("name", "?"), seq))
        elif d.startswith("r"):
            ps.rev.append((r.get("name", "?"), seq))
    return ps


# --------------------------------------------------------------------------
# bounded edit distance
# --------------------------------------------------------------------------


def _strip_affixes(a: str, b: str) -> tuple[str, str]:
    """Remove the shared prefix and suffix.

    Safe for Levenshtein - an optimal alignment can always be arranged to
    match identical flanks - and it is what makes competitive scoring
    affordable on Illumina data, where a read differs from its reference in
    a handful of places and is identical either side of them.
    """
    i, n = 0, min(len(a), len(b))
    while i < n and a[i] == b[i]:
        i += 1
    a, b = a[i:], b[i:]
    j, n = 0, min(len(a), len(b))
    while j < n and a[len(a) - 1 - j] == b[len(b) - 1 - j]:
        j += 1
    return (a[: len(a) - j], b[: len(b) - j]) if j else (a, b)


def bounded_edit(a: str, b: str, max_d: int, free_b_suffix: bool = False) -> int:
    """Levenshtein distance between a and b, or max_d + 1 if it exceeds it.

    Banded: only cells within max_d of the diagonal are computed, because a
    path that strays further has already spent more than the budget. That
    turns an O(n*m) comparison into O(n*max_d), which is what makes it
    usable against a whole panel.

    With free_b_suffix, b may run past the end of a at no cost - the mode
    needed when a read covers only the start of a reference, as the mates of
    an unmerged pair do. Without it a 250bp mate would be charged for every
    base of the amplicon it never reached.
    """
    if not free_b_suffix:
        a, b = _strip_affixes(a, b)
    n, m = len(a), len(b)
    if n == 0:
        return m if (m <= max_d or free_b_suffix) else max_d + 1
    if m == 0:
        return n if n <= max_d else max_d + 1
    if not free_b_suffix and abs(n - m) > max_d:
        return max_d + 1

    inf = max_d + 1
    prev = [j if j <= max_d else inf for j in range(m + 1)]
    for i in range(1, n + 1):
        lo, hi = max(1, i - max_d), min(m, i + max_d)
        cur = [inf] * (m + 1)
        if i <= max_d:
            cur[0] = i
        best = inf
        ai = a[i - 1]
        for j in range(lo, hi + 1):
            c = prev[j - 1] + (ai != b[j - 1])
            d = prev[j] + 1
            if d < c:
                c = d
            d = cur[j - 1] + 1
            if d < c:
                c = d
            cur[j] = c
            if c < best:
                best = c
        if best > max_d:
            return max_d + 1
        prev = cur

    out = min(prev[1:]) if free_b_suffix else prev[m]
    return out if out <= max_d else max_d + 1


# --------------------------------------------------------------------------
# io
# --------------------------------------------------------------------------


def writer(fh, fieldnames=None, delimiter: str = "\t"):
    """A csv writer that ends lines the way Unix does.

    csv defaults to CRLF regardless of platform. Every table here is read
    by awk, cut, pandas and shell loops, and a trailing carriage return
    turns the last column of every row into a value with an invisible
    character on the end - which surfaces as a file that cannot be opened
    and a path that looks exactly right in the error message.
    """
    if fieldnames is None:
        return csv.writer(fh, delimiter=delimiter, lineterminator="\n")
    return csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter,
                          lineterminator="\n")


def read_table(path: str) -> list[dict[str, str]]:
    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    return [
        {k.strip(): (v or "").strip() for k, v in row.items()}
        for row in csv.DictReader(lines, delimiter=delim)
    ]


def load_assay(path: str | None) -> dict:
    cfg = dict(DEFAULTS)
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        lib = {
            k: v for k, v in data.get("library", {}).items() if not k.startswith("_")
        }
        cfg.update({k: v for k, v in lib.items() if k in DEFAULTS})
    return cfg


FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def resolve_read_files(paths, prefer_pass: bool = True) -> tuple[list, dict]:
    """Work out which FASTQ files a nanopore run actually consists of.

    A MinION does not produce a FASTQ. It produces a directory tree of
    them, written in batches while the run proceeds, split into
    fastq_pass and fastq_fail, and split again into per-barcode folders
    if MinKNOW was asked to demultiplex. Every one of those splits is a
    decision about which reads enter the analysis, so each is made
    explicitly here and reported rather than absorbed by a glob.

    Two of them matter enough to be worth stating.

    fastq_fail is excluded. Those reads failed the basecaller's own
    quality filter, and including them silently would import a different
    error distribution into the middle of a demultiplexer whose barcode
    tolerances were set for passing reads. They are counted and reported
    so the exclusion is visible.

    MinKNOW's own barcode folders are a trap for this assay. The barcodes
    here are in-line and tandem-duplicated inside the amplicon; MinKNOW
    knows nothing about them and will have been matching its own kit
    barcodes, which are not present. Its `unclassified` folder is
    therefore where the real data is, not where the failures are - so
    everything is taken and the situation is flagged.

    Returns (files, report). Files are sorted, so a run processes in the
    same order every time.
    """
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]

    found: list[str] = []
    report = {
        "inputs": [str(p) for p in paths],
        "n_files": 0,
        "from_fastq_pass": 0,
        "excluded_fastq_fail": 0,
        "minknow_barcode_dirs": [],
        "notes": [],
    }

    def is_fastq(name: str) -> bool:
        return name.lower().endswith(FASTQ_SUFFIXES)

    for p in paths:
        p = str(p)
        if os.path.isfile(p):
            found.append(p)
            continue
        if not os.path.isdir(p):
            hits = sorted(glob.glob(p))
            if not hits:
                report["notes"].append(f"nothing matched '{p}'")
            found.extend(h for h in hits if is_fastq(h))
            continue

        passed, failed, other = [], [], []
        for root, dirs, files in os.walk(p):
            dirs.sort()
            rel = os.path.relpath(root, p).lower()
            parts = set(rel.split(os.sep))
            # Note the barcode folders BEFORE deciding whether this level
            # holds files. In a MinKNOW tree fastq_pass contains only
            # directories, so a check that runs after "no files here" is a
            # check that never runs.
            for d in dirs:
                if d.lower().startswith("barcode") or d.lower() == "unclassified":
                    report["minknow_barcode_dirs"].append(
                        os.path.relpath(os.path.join(root, d), p))

            here = [os.path.join(root, f) for f in sorted(files)
                    if is_fastq(f)]
            if not here:
                continue
            if "fastq_fail" in parts:
                failed.extend(here)
            elif "fastq_pass" in parts:
                passed.extend(here)
            else:
                other.extend(here)

        if passed:
            found.extend(passed)
            report["from_fastq_pass"] += len(passed)
            if other:
                found.extend(other)
                report["notes"].append(
                    f"{len(other)} FASTQ outside fastq_pass were also used"
                )
        else:
            found.extend(other)
        if failed:
            report["excluded_fastq_fail"] += len(failed)
            if prefer_pass:
                report["notes"].append(
                    f"{len(failed)} file(s) under fastq_fail were excluded - "
                    "those reads failed the basecaller's quality filter and "
                    "carry a different error profile from the one the "
                    "barcode tolerances were set for"
                )
            else:
                found.extend(failed)

    # A run tree can present the same file by more than one route.
    seen, uniq = set(), []
    for f in found:
        key = os.path.realpath(f)
        if key not in seen:
            seen.add(key)
            uniq.append(f)
    uniq.sort()
    report["n_files"] = len(uniq)
    if report["minknow_barcode_dirs"]:
        report["notes"].append(
            "MinKNOW barcode folders are present, so its demultiplexing was "
            "on. It does not know the in-line tandem barcodes this assay "
            "uses, so its assignment is not meaningful here and its "
            "'unclassified' reads are ordinary data. All folders were read; "
            "check that MinKNOW barcoding was not also trimming"
        )
    return uniq, report


def open_reads_many(paths):
    """Stream reads from any number of FASTQ files as one sequence."""
    for p in paths:
        yield from open_reads(p)


def open_reads(path: str):
    """Yield (read_id, sequence, quality) from a FASTQ, gzipped or not."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            yield head[1:].split()[0], seq, qual


def read_run_metadata(paths) -> dict:
    """What produced these reads, taken from the first read's header.

    Dorado writes SAM-style tags after the read id, and one of them says
    which basecalling model was used:

        RG:Z:<run uuid>_dna_r10.4.1_e8.2_400bps_hac@v5.2.0

    That belongs in the outputs. A basecalling model is not a neutral
    detail - it changes the error profile, and therefore the barcode
    recovery, the anchor hit rate and ultimately how many reads reach
    each sample. Two runs basecalled with different models can differ in
    composition for reasons that have nothing to do with biology, and
    nothing downstream would reveal it. Recording the model is what makes
    that check possible later.

    Best effort by design: the tags are absent from reads written by
    anything that rewrote the headers, and absent entirely from older
    files. An empty dict is a fine answer and never an error.
    """
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        try:
            op = gzip.open if str(path).endswith(".gz") else open
            with op(path, "rt") as fh:
                head = fh.readline().rstrip("\n")
        except OSError:
            continue
        if not head.startswith("@"):
            continue
        fields = head[1:].split("\t")[1:]
        tags = {}
        for f in fields:
            parts = f.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
        if not tags:
            continue
        out = {}
        rg = tags.get("RG", "")
        if rg:
            # RG is "<36-char run uuid>_<model>". Split on the uuid width
            # rather than on the first underscore, because model names
            # contain underscores of their own.
            if len(rg) > 37 and rg[36] == "_":
                out["acquisition_id"], out["basecall_model"] = rg[:36], rg[37:]
            else:
                out["basecall_model"] = rg
        for tag, name in (("PU", "flow_cell"), ("LB", "library"),
                          ("DT", "run_started")):
            if tags.get(tag):
                out[name] = tags[tag]
        if out:
            # Deliberately NOT recording which file this came from. The
            # same reads sharded differently would then report different
            # provenance, and a run must not describe itself differently
            # according to how its input happened to be batched.
            return out
    return {}


def phred(qual: str) -> list[int]:
    return [ord(c) - 33 for c in qual]


def expected_errors(qual: str) -> float:
    """Sum of per-base error probabilities.

    A better filter than mean quality, which is an average of logarithms and
    so lets a read through on the strength of its good bases while a handful
    of terrible ones decide the outcome.
    """
    total = 0.0
    for c in qual:
        total += 10.0 ** (-(ord(c) - 33) / 10.0)
    return total
