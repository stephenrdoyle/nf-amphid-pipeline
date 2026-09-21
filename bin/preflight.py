#!/usr/bin/env python3
"""
Pre-flight validation for nf-helminth-amplicon.

Validates every input, derives everything derivable, and writes its reports
before any read is touched. A samplesheet typo should cost seconds, not a
queue slot and an afternoon.

Hard errors block the run. Warnings are reported and do not.

Self-contained by design: this runs as a Nextflow bin/ script, so it imports
nothing outside the standard library and can be tested on its own.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import find_anchor_cost  # noqa: E402
from itertools import product
from math import comb

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# The contract: fixed in code, deliberately not configurable.
#
# The logic depends on these. neg_extraction sets detection thresholds, mock
# carries an expected composition, control requires a parent, quant carries a
# concentration. Letting a config add a sixth sample type would mean every
# downstream branch has to handle a category it knows nothing about, so the
# vocabulary stays here where the code that relies on it can be seen.
# --------------------------------------------------------------------------

SEQ_TYPES = {"target", "control", "quant", "vector"}

# What a reference is FOR, as opposed to what it is. Not exclusive: a
# single-copy locus can legitimately serve both, since its copy number does
# not vary between species the way a multi-copy marker's does.
#
# This matters before any variant calling exists. rel_abundance divides by
# total target reads, so adding a second locus without saying which
# references count would count every species twice, with a denominator
# mixing two loci at different copy numbers - plausible numbers, quietly
# wrong.
PURPOSES = {"speciation", "variants"}
SAMPLE_TYPES = {"unknown", "neg_extraction", "neg_ntc", "mock", "pos_control"}
PLATFORMS = {"illumina", "nanopore"}

IUPAC = set("ACGTRYSWKMBDHVN")
ID_SAFE = re.compile(r"[^A-Za-z0-9._-]")


# --------------------------------------------------------------------------
# The assay definition: everything that makes this a particular assay, loaded
# from a JSON file so it can be swapped without touching code. Defaults below
# match assets/assay/ampseq_helminth.json and apply when no file is given.
# --------------------------------------------------------------------------

DEFAULT_ASSAY: dict = {
    "name": "builtin_default",
    "version": "0",
    "control_tag": {
        "sequence": "TGGTCAATTTCAAGAGAT",
        "max_edits": 3,
        "use_edit_distance": True,
        "anchored": True,
    },
    "indexing": {
        "barcode_repeats": 2,
        "allow_indels": False,
        "require_copy_agreement": True,
    },
    "tolerances": {
        "corr_factor_min": 0.1,
        "corr_factor_max": 10.0,
        "min_tag_length_warn": 12,
    },
    "error_profiles": {
        "illumina": 0.001,
        "nanopore_sup": 0.01,
        "nanopore_hac": 0.02,
        "nanopore_fast": 0.05,
    },
}


class Assay:
    """Assay definition with defaults filled in for anything absent, so an
    older config file keeps working when a new key is added."""

    def __init__(self, data: dict | None = None):
        self.data = data or {}

    def section(self, name: str) -> dict:
        merged = dict(DEFAULT_ASSAY.get(name, {}))
        merged.update(
            {k: v for k, v in self.data.get(name, {}).items() if not k.startswith("_")}
        )
        return merged

    @property
    def name(self) -> str:
        return self.data.get("name", DEFAULT_ASSAY["name"])

    @property
    def version(self) -> str:
        return str(self.data.get("version", DEFAULT_ASSAY["version"]))

    @property
    def error_profiles(self) -> list[tuple[str, float]]:
        profiles = self.section("error_profiles")
        return sorted(
            ((k, float(v)) for k, v in profiles.items()),
            key=lambda kv: kv[1],
        )


def load_assay(path: str | None, diag: Diagnostics) -> Assay:
    if not path:
        return Assay()
    if not os.path.exists(path):
        diag.error(f"assay config: file not found: {path}")
        return Assay()
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        diag.error(f"assay config: {path} is not valid JSON - {exc}")
        return Assay()

    if not data.get("name"):
        diag.warn(
            f"assay config {path}: no 'name' - results cannot be traced back to "
            "the assay definition that produced them"
        )
    if not data.get("version"):
        diag.warn(f"assay config {path}: no 'version'")
    return Assay(data)


# --------------------------------------------------------------------------
# diagnostics
# --------------------------------------------------------------------------


@dataclass
class Diagnostics:
    """Collects errors and warnings so every problem is reported at once,
    rather than the run dying on the first one and hiding the rest."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    # A third tier, for properties that are correct and load-bearing. A
    # deliberate design choice that a later edit could silently undo is
    # worth restating every run; filing it as a warning would teach people
    # to ignore warnings, which is worse than not saying it at all.
    notes: list[str] = field(default_factory=list)

    def error(self, msg: str) -> None:
        self.errors.append(msg)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)

    def info(self, msg: str) -> None:
        self.notes.append(msg)

    @property
    def ok(self) -> bool:
        return not self.errors

    def render(self) -> str:
        out = []
        if self.errors:
            out.append(f"ERRORS ({len(self.errors)}) - the run cannot proceed:")
            out += [f"  [E] {m}" for m in self.errors]
        if self.warnings:
            if out:
                out.append("")
            out.append(f"WARNINGS ({len(self.warnings)}):")
            out += [f"  [W] {m}" for m in self.warnings]
        if self.notes:
            if out:
                out.append("")
            out.append(f"NOTES ({len(self.notes)}) - correct, and worth knowing:")
            out += [f"  [i] {m}" for m in self.notes]
        if not out:
            out.append("All checks passed with no warnings.")
        return "\n".join(out)


# --------------------------------------------------------------------------
# sequence helpers
# --------------------------------------------------------------------------


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return seq.translate(table)[::-1]


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("hamming requires equal lengths")
    return sum(x != y for x, y in zip(a, b))


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def approx_substring_distance(pattern: str, text: str) -> int:
    """Minimum edit distance between `pattern` and any substring of `text`
    (Sellers' algorithm). Used to check a tag does not occur anywhere in the
    panel, allowing for the same tolerance the read-level search will use."""
    if not pattern:
        return 0
    prev = [0] * (len(text) + 1)  # free start anywhere in text
    for i, pc in enumerate(pattern, 1):
        cur = [i]
        for j, tc in enumerate(text, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (pc != tc)))
        prev = cur
    return min(prev)


def rotations(seq: str) -> set[str]:
    return {seq[i:] + seq[:i] for i in range(1, len(seq))}


def terminal_groups(
    seqs: dict[str, str], end: str = "5", motif: int = 18, max_mismatch: int = 2
) -> list[list[str]]:
    """Cluster sequences by shared terminal sequence.

    Entries amplified by one primer pair carry that primer at their ends, so
    they share terminal sequence however much their middles differ. Grouping
    on that separates "different species, same amplicon" - which is expected
    and where length varies freely - from "different amplicon" or "trimmed
    differently", which are the cases worth reporting.

    A small mismatch budget absorbs degenerate primer positions and the
    single-base differences real primers accumulate.
    """
    ids = sorted(seqs)
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    def term(seq):
        return seq[:motif] if end == "5" else seq[-motif:]

    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            ta, tb = term(seqs[a]), term(seqs[b])
            if len(ta) == len(tb) and hamming(ta, tb) <= max_mismatch:
                union(a, b)

    out: dict[str, list[str]] = defaultdict(list)
    for i in ids:
        out[find(i)].append(i)
    return sorted(out.values(), key=lambda g: (-len(g), g[0]))


def detection_probability(length: int, error_rate: float, max_errors: int) -> float:
    """Probability that at most `max_errors` of `length` bases are miscalled.
    A floor on real detection: allowing edits rather than only mismatches can
    only help, since an indel is absorbed rather than propagated."""
    return sum(
        comb(length, k) * error_rate**k * (1 - error_rate) ** (length - k)
        for k in range(max_errors + 1)
    )


def sanitise_id(value: str) -> str:
    """Path-safe form of an identifier. The original is always what gets
    reported; this is only what appears in directory and file names."""
    return ID_SAFE.sub("_", value.strip())


# --------------------------------------------------------------------------
# tabular input
# --------------------------------------------------------------------------


def read_table(path: str, diag: Diagnostics, label: str) -> list[dict[str, str]]:
    """Read a TSV or CSV, skipping comment lines and blank rows. Delimiter is
    taken from the extension, because sniffing guesses wrongly on files whose
    first column happens to contain no delimiter."""
    if not os.path.exists(path):
        diag.error(f"{label}: file not found: {path}")
        return []

    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]

    if not lines:
        diag.error(f"{label}: file is empty: {path}")
        return []

    # The delimiter comes from the extension, never from sniffing - a
    # sniffer guesses wrongly on a file whose first column happens to
    # contain no separator, and guessing wrongly about which column is
    # which is worse than refusing. But a .txt holding CSV is an ordinary
    # mistake, and read as TSV it yields one wide column and a complaint
    # that every required field is missing. Say what actually happened.
    other = "\t" if delim == "," else ","
    header = lines[0].rstrip("\n")
    if header.count(delim) == 0 and header.count(other) >= 2:
        names = {",": "comma", "\t": "tab"}
        diag.error(
            f"{label}: read as {names[delim]}-separated because of the "
            f"'{os.path.splitext(path)[1] or 'no'}' extension, but the "
            f"header has no {names[delim]}s and "
            f"{header.count(other)} {names[other]}s - so it is almost "
            f"certainly {names[other]}-separated. Rename it to "
            f"'{'.csv' if other == ',' else '.tsv'}', or convert the "
            "separators. The delimiter is taken from the extension "
            "deliberately: sniffing it guesses wrongly on files whose "
            "first column contains no separator, and a wrong guess "
            "mislabels columns silently"
        )
        return []

    rows = list(csv.DictReader(lines, delimiter=delim))
    if not rows:
        diag.error(f"{label}: no data rows below the header: {path}")
        return []

    cleaned = []
    for row in rows:
        cleaned.append(
            {
                (k.strip() if k else k): (v.strip() if isinstance(v, str) else "")
                for k, v in row.items()
            }
        )
    return cleaned


def require_columns(
    rows: list[dict[str, str]], needed: set[str], diag: Diagnostics, label: str
) -> bool:
    if not rows:
        return False
    present = set(rows[0].keys())
    missing = needed - present
    if missing:
        diag.error(f"{label}: missing required column(s): {', '.join(sorted(missing))}")
        return False
    return True


# --------------------------------------------------------------------------
# reference panel
# --------------------------------------------------------------------------


@dataclass
class RefEntry:
    seq_id: str
    locus: str
    type: str
    taxon: str
    parent_id: str
    default_conc: str
    conc_units: str
    corr_factor: str
    sequence: str
    purpose: set = field(default_factory=set)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def report_taxon(self) -> str:
        return self.taxon or self.seq_id


def parse_reference(path: str, diag: Diagnostics) -> list[RefEntry]:
    rows = read_table(path, diag, "reference")
    needed = {"seq_id", "locus", "type", "sequence"}
    if not require_columns(rows, needed, diag, "reference"):
        return []

    entries: list[RefEntry] = []
    for i, row in enumerate(rows, start=2):  # header is line 1
        entries.append(
            RefEntry(
                seq_id=row.get("seq_id", ""),
                locus=row.get("locus", ""),
                type=row.get("type", "").lower(),
                purpose={
                    p.strip().lower()
                    for p in row.get("purpose", "").split(",")
                    if p.strip()
                },
                taxon=row.get("taxon", ""),
                parent_id=row.get("parent_id", ""),
                default_conc=row.get("default_conc", ""),
                conc_units=row.get("conc_units", ""),
                corr_factor=row.get("corr_factor", ""),
                sequence=row.get("sequence", "").upper().replace(" ", ""),
            )
        )
        if not entries[-1].seq_id:
            diag.error(f"reference line {i}: seq_id is empty")
    return entries


def validate_reference(
    entries: list[RefEntry], diag: Diagnostics, assay: Assay | None = None
) -> None:
    assay = assay or Assay()
    tol = assay.section("tolerances")
    corr_min = float(tol["corr_factor_min"])
    corr_max = float(tol["corr_factor_max"])
    by_id = {e.seq_id: e for e in entries}

    # -- identifiers ------------------------------------------------------
    dupes = [k for k, n in Counter(e.seq_id for e in entries).items() if n > 1]
    for d in dupes:
        diag.error(f"reference: seq_id '{d}' appears more than once")

    sanitised: dict[str, list[str]] = defaultdict(list)
    for e in entries:
        sanitised[sanitise_id(e.seq_id)].append(e.seq_id)
    for token, originals in sanitised.items():
        if len(set(originals)) > 1:
            diag.error(
                f"reference: seq_ids {originals} all sanitise to '{token}' - "
                "silently merging them would be worse than refusing"
            )

    # -- identical sequences ---------------------------------------------
    # Two identical references break competitive assignment silently: no
    # alignment score can ever separate them.
    seen: dict[str, str] = {}
    for e in entries:
        if not e.sequence:
            diag.error(f"reference '{e.seq_id}': sequence is empty")
            continue
        if e.sequence in seen:
            diag.error(
                f"reference: '{e.seq_id}' and '{seen[e.sequence]}' have identical "
                "sequences - competitive assignment cannot separate them"
            )
        else:
            seen[e.sequence] = e.seq_id

    for e in entries:
        # -- type ---------------------------------------------------------
        if e.type not in SEQ_TYPES:
            diag.error(
                f"reference '{e.seq_id}': type '{e.type}' not in "
                f"{sorted(SEQ_TYPES)}"
            )

        # -- sequence characters -----------------------------------------
        bad = set(e.sequence) - IUPAC
        if bad:
            diag.error(
                f"reference '{e.seq_id}': sequence contains non-IUPAC "
                f"character(s): {''.join(sorted(bad))}"
            )
        if "-" in e.sequence or "." in e.sequence:
            diag.error(f"reference '{e.seq_id}': sequence contains gaps")

        # -- locus --------------------------------------------------------
        if e.type != "vector" and not e.locus:
            diag.error(f"reference '{e.seq_id}': locus is required for type '{e.type}'")

        # -- controls need a resolvable parent at the same locus ----------
        if e.type == "control":
            if not e.parent_id:
                diag.error(
                    f"reference '{e.seq_id}': parent_id is required for a control - "
                    "a decoy is a decoy for something specific"
                )
            elif e.parent_id not in by_id:
                diag.error(
                    f"reference '{e.seq_id}': parent_id '{e.parent_id}' does not "
                    "match any seq_id"
                )
            else:
                parent = by_id[e.parent_id]
                if parent.type != "target":
                    diag.error(
                        f"reference '{e.seq_id}': parent '{e.parent_id}' has type "
                        f"'{parent.type}', expected 'target'"
                    )
                if parent.locus != e.locus:
                    diag.error(
                        f"reference '{e.seq_id}': locus '{e.locus}' differs from "
                        f"parent '{e.parent_id}' locus '{parent.locus}'"
                    )
        elif e.parent_id:
            diag.warn(
                f"reference '{e.seq_id}': parent_id is set but type is "
                f"'{e.type}' - it is only meaningful for controls, and is ignored"
            )

        # -- quants need a concentration and its units -------------------
        if e.type == "quant":
            if not e.default_conc:
                diag.error(
                    f"reference '{e.seq_id}': default_conc is required for a quant"
                )
            else:
                try:
                    if float(e.default_conc) <= 0:
                        diag.error(
                            f"reference '{e.seq_id}': default_conc must be positive"
                        )
                except ValueError:
                    diag.error(
                        f"reference '{e.seq_id}': default_conc "
                        f"'{e.default_conc}' is not a number"
                    )
            if not e.conc_units:
                diag.error(
                    f"reference '{e.seq_id}': conc_units is required for a quant - "
                    "an implicit unit is how a copy number becomes uninterpretable"
                )
        elif e.default_conc:
            diag.warn(
                f"reference '{e.seq_id}': default_conc is set on a "
                f"'{e.type}' entry and will be ignored"
            )

        # -- correction factor -------------------------------------------
        if e.corr_factor:
            try:
                cf = float(e.corr_factor)
            except ValueError:
                diag.error(
                    f"reference '{e.seq_id}': corr_factor "
                    f"'{e.corr_factor}' is not a number"
                )
            else:
                if cf <= 0:
                    diag.error(
                        f"reference '{e.seq_id}': corr_factor must be positive"
                    )
                elif not corr_min <= cf <= corr_max:
                    diag.warn(
                        f"reference '{e.seq_id}': corr_factor {cf} is outside "
                        f"{corr_min}-{corr_max} - check it is "
                        "not a typo or an inverted value"
                    )

    quantifiable = [e for e in entries if e.type in {"target", "quant"}]
    defaulted = [e for e in quantifiable if not e.corr_factor]
    if defaulted and len(defaulted) < len(quantifiable):
        diag.warn(
            f"{len(defaulted)} of {len(quantifiable)} quantifiable entries "
            "have no corr_factor and will use the default of 1.0, while the "
            "rest carry measured values. Mixing the two puts corrected and "
            "assumed copy numbers in the same column - comparable only if "
            "the unmeasured ones really do amplify at unity: "
            + ", ".join(e.seq_id for e in defaulted[:5])
            + ("..." if len(defaulted) > 5 else "")
        )
    elif defaulted:
        diag.info(
            f"no corr_factor anywhere in the panel, so all {len(defaulted)} "
            "quantifiable entries use the default of 1.0. Copy numbers are "
            "then relative to the standards' amplification efficiency rather "
            "than corrected to it: a reference that amplifies twice as well "
            "as the standards will read twice as abundant. That is a "
            "consistent basis for comparing samples with each other, and not "
            "one for comparing species with each other"
        )

    # -- purpose -----------------------------------------------------------
    for e in entries:
        bad = e.purpose - PURPOSES
        if bad:
            diag.error(
                f"reference '{e.seq_id}': purpose {sorted(bad)} not in "
                f"{sorted(PURPOSES)}"
            )
        if e.type in {"quant", "vector"} and e.purpose:
            diag.warn(
                f"reference '{e.seq_id}': purpose is set on a '{e.type}' entry "
                "and is ignored - quants are calibration and vectors are "
                "sensors, so neither is a speciation or variant target"
            )

    # A control serves the same role as the target it mimics, so inherit
    # rather than require - and complain if someone sets it differently,
    # since a decoy that disagrees with its parent is almost certainly a
    # mistake rather than an intention.
    for e in entries:
        if e.type != "control" or e.parent_id not in by_id:
            continue
        parent_purpose = by_id[e.parent_id].purpose or {"speciation"}
        if not e.purpose:
            e.purpose = set(parent_purpose)
        elif e.purpose != parent_purpose:
            diag.warn(
                f"reference '{e.seq_id}': purpose {sorted(e.purpose)} differs "
                f"from its parent '{e.parent_id}' {sorted(parent_purpose)}. A "
                "decoy serves the same role as the target it mimics"
            )

    # Blank means speciation, so panels written before this column existed
    # keep behaving exactly as they did.
    for e in entries:
        if e.type in {"target", "control"} and not e.purpose:
            e.purpose = {"speciation"}

    # Purpose is really a property of the amplicon. One haplotype flagged
    # differently from its locus-mates is the failure mode worth catching.
    purpose_by_locus: dict[str, set] = defaultdict(set)
    for e in entries:
        if e.type in {"target", "control"} and e.locus:
            purpose_by_locus[e.locus].add(",".join(sorted(e.purpose)))
    for locus, seen in sorted(purpose_by_locus.items()):
        if len(seen) > 1:
            diag.warn(
                f"locus '{locus}' mixes purposes ({', '.join(sorted(seen))}). "
                "Purpose is normally a property of the amplicon, so check this "
                "is deliberate rather than one entry flagged by mistake"
            )

    # -- amplicon coherence within a locus --------------------------------
    # The old check here compared the length range at each locus, which
    # cannot tell genuine biological variation from a file that mixes
    # insert-only and full-amplicon entries. ITS2 really does span 350-500bp
    # across genera, so that check fired on correct panels and taught people
    # to ignore it.
    #
    # Terminal sequence is the direct test: entries from one primer pair
    # carry that primer at both ends. So cluster on the termini first, and
    # only then ask whether lengths are coherent - within a group, where a
    # large difference really does suggest one entry was trimmed differently.
    by_locus: dict[str, list[RefEntry]] = defaultdict(list)
    for e in entries:
        if e.type in {"target", "control", "quant"} and e.locus and e.sequence:
            by_locus[e.locus].append(e)

    for locus, group in sorted(by_locus.items()):
        if len(group) < 2:
            continue
        seqs = {e.seq_id: e.sequence for e in group}
        g5 = terminal_groups(seqs, "5")
        g3 = terminal_groups(seqs, "3")

        if len(g5) > 1 or len(g3) > 1:
            diag.warn(
                f"locus '{locus}' contains {len(g5)} distinct 5' and "
                f"{len(g3)} distinct 3' termini, so one locus label is "
                "covering more than one primer pair. Not an error - but "
                "`locus` is what length filtering and reporting group by, and "
                "a quant carrying one primer pair only directly controls "
                "amplification for that pair. 5' groups: "
                + "; ".join(
                    "[" + ", ".join(sorted(g)) + "]" for g in g5
                )
            )

        # An entry trimmed differently from the rest does not share their
        # terminal motif, so it falls out as a singleton here. That only
        # means something when there is a clear majority to be the odd one
        # out *of* - in a panel that legitimately uses several primer pairs,
        # singletons are expected.
        biggest = len(g5[0]) if g5 else 0
        if biggest >= 0.7 * len(group):
            for g in g5[1:]:
                if len(g) == 1:
                    diag.warn(
                        f"locus '{locus}': '{g[0]}' does not share the 5' "
                        f"terminus that {biggest} of {len(group)} entries at "
                        "this locus have. Either it is a different amplicon, "
                        "or it was trimmed differently - and mixing "
                        "insert-only with full-amplicon biases assignment "
                        "toward the longer entries"
                    )

    # -- panel sanity ----------------------------------------------------
    quants = [e for e in entries if e.type == "quant"]
    if not quants:
        diag.warn(
            "reference: no quant entries - copy-number estimation is unavailable "
            "and every negative will be an unverified negative"
        )
    elif len(quants) == 1:
        diag.warn(
            "reference: only one quant entry - a single point gives a ratio but no "
            "calibration bracket, so estimates cannot be flagged in or out of range"
        )

    if not any(e.type == "target" for e in entries):
        diag.error("reference: no target entries - there is nothing to quantify")

    if not any(e.type == "vector" for e in entries):
        diag.warn(
            "reference: no vector entry - the plasmid backbone sensor is "
            "unavailable, which is not the same as clear"
        )


# --------------------------------------------------------------------------
# control tags
# --------------------------------------------------------------------------


@dataclass
class TagResult:
    control_id: str
    parent_id: str
    locus: str
    kind: str  # "insertion" | "substitutions" | "indeterminate"
    tag: str
    start: int  # 0-based position in the control sequence
    length: int
    substitutions: list[int] = field(default_factory=list)
    note: str = ""


def derive_tag(
    control: str, parent: str, expected_tag: str = ""
) -> TagResult | None:
    """Find what distinguishes a control from its parent.

    For a pure insertion - which is how these decoys are built - the flanks
    are identical, so the longest common prefix and suffix bracket the tag.
    No alignment needed and no scoring choices to get wrong.

    One trap: the insertion point is ambiguous whenever the inserted sequence
    shares bases with the flanking sequence. If the parent ends '...TG' where
    the tag begins 'TG', the common prefix runs two bases into the tag and
    what comes back is a rotation of it - self-consistent, and detection
    against it would still work, but it would not match the tag as designed.
    So the insertion is left-aligned into its canonical leftmost position,
    the same normalisation VCF uses for indels, and then cross-checked
    against the expected sequence when one is supplied.
    """
    n = min(len(control), len(parent))

    lcp = 0
    while lcp < n and control[lcp] == parent[lcp]:
        lcp += 1

    lcs = 0
    while lcs < n - lcp and control[-1 - lcs] == parent[-1 - lcs]:
        lcs += 1

    delta = len(control) - len(parent)

    if delta > 0 and lcp + lcs >= len(parent):
        # flanks account for the whole parent: the middle of the control is
        # inserted sequence
        start = lcp
        length = delta
        if len(control[start : len(control) - lcs]) == delta:
            # left-align: slide the insertion window leftwards while doing so
            # leaves the sequence unchanged
            while start > 0 and control[start - 1] == control[start + length - 1]:
                start -= 1

            tag = control[start : start + length]
            note = ""

            if expected_tag and tag != expected_tag:
                if len(tag) == len(expected_tag) and tag in rotations(expected_tag):
                    # re-anchor onto the expected reading frame, which exists
                    # somewhere in the tandem of insertion and flank
                    offset = next(
                        i
                        for i in range(len(expected_tag))
                        if expected_tag[i:] + expected_tag[:i] == tag
                    )
                    shift = len(expected_tag) - offset
                    cand = control[start + shift : start + shift + length]
                    if cand == expected_tag:
                        start += shift
                        tag = cand
                        note = (
                            "insertion point was ambiguous because the tag shares "
                            "bases with the flanking sequence; re-anchored onto the "
                            "expected tag"
                        )
                else:
                    note = (
                        f"inserted sequence differs from the expected tag "
                        f"({expected_tag}) and is not a rotation of it"
                    )

            return TagResult(
                control_id="",
                parent_id="",
                locus="",
                kind="insertion",
                tag=tag,
                start=start,
                length=length,
                note=note,
            )

    if delta == 0:
        subs = [i for i, (a, b) in enumerate(zip(control, parent)) if a != b]
        return TagResult(
            control_id="",
            parent_id="",
            locus="",
            kind="substitutions",
            tag="",
            start=subs[0] if subs else -1,
            length=0,
            substitutions=subs,
            note=(
                "control and parent are the same length - discrimination rests on "
                f"{len(subs)} substitution(s) rather than an inserted tag"
            ),
        )

    return TagResult(
        control_id="",
        parent_id="",
        locus="",
        kind="indeterminate",
        tag="",
        start=-1,
        length=0,
        note=(
            f"length differs by {delta} bp but the flanks do not bracket a clean "
            "insertion - the construct may contain both an insertion and "
            "substitutions, or the wrong parent is assigned"
        ),
    )


def analyse_tags(
    entries: list[RefEntry],
    diag: Diagnostics,
    max_edits: int,
    expected_tag: str = "",
    assay: Assay | None = None,
) -> list[TagResult]:
    assay = assay or Assay()
    min_tag_warn = int(assay.section("tolerances")["min_tag_length_warn"])
    by_id = {e.seq_id: e for e in entries}
    results: list[TagResult] = []
    collisions: set[tuple[str, str, str, int]] = set()

    for e in entries:
        if e.type != "control" or e.parent_id not in by_id:
            continue
        parent = by_id[e.parent_id]
        res = derive_tag(e.sequence, parent.sequence, expected_tag)
        if res is None:
            continue
        res.control_id = e.seq_id
        res.parent_id = parent.seq_id
        res.locus = e.locus
        results.append(res)

        if res.kind == "insertion":
            if expected_tag and res.tag != expected_tag:
                diag.error(
                    f"control '{e.seq_id}': inserted sequence '{res.tag}' does not "
                    f"match the expected tag '{expected_tag}'. {res.note}"
                )
            # verify the derived coordinates actually hold, so a coordinate
            # bug cannot pass silently into read-level detection
            if e.sequence[res.start : res.start + res.length] != res.tag:
                diag.error(
                    f"control '{e.seq_id}': derived tag does not occur at the "
                    f"derived position {res.start} - internal inconsistency"
                )
            # A tag must not occur naturally anywhere in the panel, or the
            # read-level search will call target reads control.
            for other in entries:
                if other.seq_id in {e.seq_id, parent.seq_id}:
                    continue
                if other.type == "control":
                    continue  # controls are expected to share the tag
                d = approx_substring_distance(res.tag, other.sequence)
                if d <= max_edits:
                    # Recorded rather than reported here: controls share one
                    # tag, so a single collision would otherwise be announced
                    # once per control and bury everything else.
                    collisions.add((res.tag, other.seq_id, other.type, d))

            if res.length < min_tag_warn:
                diag.warn(
                    f"control '{e.seq_id}': tag is only {res.length} bp. Short tags "
                    "match by chance more often and detect less reliably under "
                    "nanopore error - check the predicted rates in the tag report"
                )
        elif res.kind == "substitutions":
            diag.warn(
                f"control '{e.seq_id}': {res.note}. Substitution-only decoys are "
                "far weaker discriminators than an inserted tag, especially on "
                "nanopore"
            )
            if not res.substitutions:
                diag.error(
                    f"control '{e.seq_id}' is identical to its parent "
                    f"'{parent.seq_id}' - nothing can distinguish them"
                )
        else:
            diag.warn(f"control '{e.seq_id}': {res.note}")

    # Tag collisions, reported once per affected sequence rather than once
    # per control. Severity depends on whether that sequence can amplify.
    for tag, seq_id, seq_type, d in sorted(collisions, key=lambda c: (c[3], c[1])):
        if seq_type == "vector" and d > 0:
            # A vector carries no primer sites, so it barely amplifies, and
            # its reads do not present the tag window where the anchored
            # read-level search looks. A chance match can at worst
            # misclassify the few non-specific backbone reads that appear -
            # it cannot turn target reads into control calls.
            diag.warn(
                f"the control tag has a {d}-edit match inside '{seq_id}' "
                f"(vector). Chance matches occur in a backbone of that length "
                "once indels are allowed. The risk is that a few non-specific "
                "backbone reads are misclassified, not that target reads are"
            )
        else:
            diag.error(
                f"the control tag occurs in '{seq_id}' ({seq_type}) within "
                f"{d} edit(s) of an exact match. At the read-level tolerance "
                f"of {max_edits}, reads from that sequence would be called "
                "control"
            )

    # Report how many distinct tags are in play. A shared tag means detection
    # says only "control"; the mapped locus says which one.
    tags = {r.tag for r in results if r.kind == "insertion"}
    if len(tags) == 1 and len(results) > 1:
        diag.warn(
            f"all {len(results)} controls share one tag sequence - tag detection "
            "establishes that a read is control-derived, and the mapped locus "
            "establishes which control. This is expected, and noted so the "
            "two-step assignment is not mistaken for a bug"
        )

    return results


def write_tag_report(
    results: list[TagResult], outdir: str, max_edits: int, assay: Assay | None = None
) -> None:
    assay = assay or Assay()
    profiles = assay.error_profiles
    rows = []
    for r in results:
        row = {
            "control_id": r.control_id,
            "parent_id": r.parent_id,
            "locus": r.locus,
            "kind": r.kind,
            "tag": r.tag,
            "tag_length": r.length,
            "tag_start_in_control": r.start,
            "n_substitutions": len(r.substitutions),
            "note": r.note,
        }
        for name, err in profiles:
            if r.kind == "insertion":
                p = detection_probability(r.length, err, max_edits)
                row[f"detect_{name}"] = f"{p:.5f}"
            else:
                row[f"detect_{name}"] = "NA"
        rows.append(row)

    path = os.path.join(outdir, "tag_report.tsv")
    if rows:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    else:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("# no control entries in the reference panel\n")

    # A human-readable companion: the numbers that decide whether a decoy is
    # usable on a given platform.
    lines = [
        "Control tag report",
        "==================",
        "",
        f"Read-level matching tolerance: <={max_edits} edits.",
        "Detection columns give the probability that at most that many of the",
        "tag's bases are miscalled, per platform error profile. Edit-distance",
        "matching can only do better, since an indel is absorbed rather than",
        "shifting every subsequent base.",
        "",
    ]
    for r in results:
        lines.append(f"{r.control_id}  (parent {r.parent_id}, locus {r.locus})")
        if r.kind == "insertion":
            lines.append(
                f"  inserted tag: {r.tag}  ({r.length} bp at position {r.start})"
            )
            for name, err in profiles:
                p = detection_probability(r.length, err, max_edits)
                flag = "" if p >= 0.98 else "   <-- marginal"
                lines.append(f"    {name:<16} err {err*100:>4.1f}%   {p*100:6.2f}%{flag}")
        else:
            lines.append(f"  {r.kind}: {r.note}")
        lines.append("")
    with open(os.path.join(outdir, "tag_report.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


# --------------------------------------------------------------------------
# primers: optional, and only for validating reference extent
# --------------------------------------------------------------------------
#
# Primers are not needed to assign a locus (competitive mapping gives that)
# nor to set length windows (reference length gives those). But they are the
# only thing that can check reference EXTENT - whether every entry spans
# exactly primer-to-primer, or some carry extra sequence outside the primer
# sites.
#
# That matters because an entry with extra flanking sequence offers more
# alignable bases than its competitors, biasing assignment toward it for
# reasons unrelated to which sequence is actually present. And it hides
# inside genuine biological length variation, so no length-based proxy
# finds it.


def parse_primers(path: str, diag: Diagnostics) -> dict[str, list[tuple[str, str]]]:
    rows = read_table(path, diag, "primers")
    if not require_columns(rows, {"name", "direction", "sequence"}, diag, "primers"):
        return {}
    out: dict[str, list[tuple[str, str]]] = {"fwd": [], "rev": []}
    for r in rows:
        d = r["direction"].lower()
        if d not in out:
            diag.error(
                f"primers '{r['name']}': direction '{d}' is not fwd or rev"
            )
            continue
        seq = r["sequence"].upper()
        bad = set(seq) - IUPAC
        if bad:
            diag.error(
                f"primers '{r['name']}': non-IUPAC character(s) "
                f"{''.join(sorted(bad))}"
            )
            continue
        out[d].append((r["name"], seq))
    if not out["fwd"] or not out["rev"]:
        diag.error("primers: need at least one fwd and one rev primer")
    return out


def check_reference_extent(
    entries: list[RefEntry], primers: dict, diag: Diagnostics
) -> list[dict]:
    """Locate each entry's primer pair and measure what lies outside it.

    Approximate matching, not substring search. A primer rarely matches its
    template exactly: the primer may be degenerate, and the template may
    carry a mismatch under it - which is the ordinary case, not the
    pathological one. Requiring an exact match means the primer is not found
    in precisely those entries, and the reference is then left untrimmed
    while its competitors are trimmed, handing it 40 extra alignable bases.
    The bias would fall on the divergent sequences: the ones the panel exists
    to tell apart.
    """
    if not primers.get("fwd") or not primers.get("rev"):
        return []

    def locate(seq: str, probe: str):
        hit = find_anchor_cost(seq, probe, max(1, len(probe) // 8))
        return None if hit is None else hit

    findings = []
    for e in entries:
        if e.type not in {"target", "control", "quant"} or not e.sequence:
            continue
        seq = e.sequence

        fwd, rev = [], []
        for n, p in primers["fwd"]:
            hit = locate(seq, p)
            if hit:
                fwd.append((n, hit[0], hit[1], hit[2]))
        for n, p in primers["rev"]:
            hit = locate(seq, revcomp(p))
            if hit:
                rev.append((n, hit[0], hit[1], hit[2]))

        if not fwd or not rev:
            diag.warn(
                f"reference '{e.seq_id}': no "
                + ("forward" if not fwd else "reverse")
                + " primer found in the sequence. Either it is not from this "
                "panel, or it has been trimmed to the insert - in which case "
                "it offers fewer alignable bases than the rest"
            )
            continue

        # Best fit first, then earliest (forward) or latest (reverse): a
        # primer that matches in two places should be taken where it fits
        # best, and only then where it sits at the amplicon edge.
        fwd.sort(key=lambda x: (x[3], x[1]))
        rev.sort(key=lambda x: (x[3], -x[1]))
        f_name, f_pos, f_end, _ = fwd[0]
        r_name, r_pos, r_end, _ = rev[0]
        r_len = r_end - r_pos
        findings.append(
            {
                "seq_id": e.seq_id,
                "type": e.type,
                "locus": e.locus,
                "fwd": f_name,
                "rev": r_name,
                "extra_5p": f_pos,
                "extra_3p": len(seq) - (r_pos + r_len),
                # Where the insert begins and ends: reads arrive with the
                # primers trimmed off, so the reference they compete against
                # has to be trimmed the same way or every read looks short
                # against every reference.
                "insert_start": f_end,
                "insert_end": r_pos,
            }
        )

    if not findings:
        return findings

    # Consistency is what matters, not the value itself: a panel where every
    # entry carries the same flank is coherent; a panel where some do and
    # some do not is not.
    off = [f for f in findings if f["extra_5p"] or f["extra_3p"]]
    if off and len(off) < len(findings):
        detail = "; ".join(
            f"{f['seq_id']} (+{f['extra_5p']}bp 5', +{f['extra_3p']}bp 3')"
            for f in sorted(off, key=lambda f: f["seq_id"])[:6]
        )
        diag.warn(
            f"reference extent is not uniform: {len(findings) - len(off)} of "
            f"{len(findings)} entries run exactly primer-to-primer, but "
            f"{len(off)} carry sequence outside the primer sites - {detail}. "
            "Extra flanking sequence gives those entries more alignable bases "
            "than their competitors, biasing assignment toward them. Trim to "
            "the primer boundaries so every entry uses the same convention"
        )
    return findings


def check_tag_in_insert(entries, tags, extent, diag) -> list[dict]:
    """Does the tag survive primer trimming?

    Reads arrive with their primers removed, so anything under a primer
    binding site is gone before assignment sees it. A tag overlapping that
    site is not shortened in a way anything notices - it is simply matched
    at a worse edit distance, and the tolerance absorbs the loss until it
    does not. Target-versus-control then degrades quietly with the
    tolerance, which is the worst way for a deterministic test to fail.

    Margin is reported even when it is fine, because zero margin is a
    deliberate choice that a later edit to the primers can silently undo.
    """
    bounds = {f["seq_id"]: (f["insert_start"], f["insert_end"]) for f in extent}
    tight = []
    findings = []
    for t in tags:
        if t.kind != "insertion" or t.control_id not in bounds:
            continue
        a, b = bounds[t.control_id]
        start, end = t.start, t.start + t.length
        inside = start >= a and end <= b
        margin = min(start - a, b - end)
        findings.append(
            {"control_id": t.control_id, "inside_insert": inside,
             "margin_bp": margin}
        )
        if not inside:
            lost = max(0, a - start) + max(0, end - b)
            diag.error(
                f"control '{t.control_id}': {lost}bp of the {t.length}bp tag "
                "lies under a primer binding site, so it is replaced by "
                "primer sequence in every read and never observed. Tag "
                "detection would fall back on the remaining "
                f"{t.length - lost}bp and weaken as the tolerance is "
                "approached. Move the tag inside the primers"
            )
        elif margin == 0:
            tight.append(t.control_id)

    if tight:
        diag.info(
            f"{len(tight)} control(s) place the tag flush against a primer "
            "boundary with no margin. That is the strongest position for "
            "Illumina - it is the first thing the mate reads after its "
            "primer, in that read's best bases - but it leaves no room: a "
            "primer redesign that shifts the boundary by even one base puts "
            "the tag under it. Worth re-running pre-flight whenever the "
            "primers change"
        )
    return findings


def check_tag_visibility(
    entries, tags, diag, read_length: int, min_overlap: int = 20
) -> list[dict]:
    """Can a read pair actually see the tag?

    Target and control differ by the tag and nothing else. If the amplicon
    is short enough to merge, the merged read spans the whole insert and the
    tag is always observable. If it is not, the pair sees only the two ends,
    and a tag sitting in the middle falls in the gap between the mates.

    Nothing about that failure looks like a failure. Both references still
    fit the pair equally well, and a caller that treats a missing tag as
    evidence of target will report every plasmid carryover read at that
    locus as a genuine detection. The pipeline refuses to make that call,
    which converts a silent false positive into a visible loss of data -
    but the place to fix it is the panel, so it is flagged here.

    The remedy is to put the tag within a read length of one end, which also
    makes the assay robust to a shorter run than the one it was designed
    for.
    """
    lengths = {e.seq_id: len(e.sequence) for e in entries if e.sequence}
    findings = []
    for t in tags:
        m = lengths.get(t.control_id)
        if not m or t.kind != "insertion":
            continue
        merges = m <= 2 * read_length - min_overlap
        visible = (t.start + t.length) <= read_length or t.start >= m - read_length
        findings.append(
            {
                "control_id": t.control_id,
                "amplicon_bp": m,
                "tag_start": t.start,
                "merges": merges,
                "tag_visible_unmerged": visible,
            }
        )
        if not merges and not visible:
            diag.error(
                f"control '{t.control_id}': the amplicon is {m}bp, too long "
                f"to merge at {read_length}bp reads, and the tag sits at "
                f"position {t.start} - between the two mates, where neither "
                "read can see it. This control cannot be told apart from "
                f"'{t.parent_id}' on this run configuration. Move the tag "
                f"within {read_length}bp of either end, or sequence longer"
            )
    marginal = [
        f for f in findings
        if f["merges"] and not f["tag_visible_unmerged"]
    ]
    if marginal:
        diag.warn(
            f"{len(marginal)} control(s) rely on the mates merging for the "
            f"tag to be seen at all: {', '.join(f['control_id'] for f in marginal[:4])}"
            + ("..." if len(marginal) > 4 else "")
            + f". At {read_length}bp they merge and the tag is readable, but "
            "a shorter run, or read trimming that eats into the overlap, "
            "removes target-versus-control discrimination at those loci "
            "rather than degrading it"
        )
    return findings


# --------------------------------------------------------------------------
# mock communities
# --------------------------------------------------------------------------

MOCK_UNITS = {"copies", "fraction", "ng"}


def parse_mocks(path: str, diag: Diagnostics) -> dict[str, list[dict]]:
    rows = read_table(path, diag, "mocks")
    if not require_columns(
        rows, {"sample_id", "seq_id", "amount", "amount_units", "mix"},
        diag, "mocks",
    ):
        return {}
    by_sample: dict[str, list[dict]] = {}
    for r in rows:
        sid = r.get("sample_id", "")
        if not sid:
            continue
        try:
            amount = float(r["amount"])
        except (ValueError, KeyError):
            diag.error(
                f"mocks '{sid}'/{r.get('seq_id','?')}: amount "
                f"'{r.get('amount')}' is not a number"
            )
            continue
        if amount <= 0:
            diag.error(
                f"mocks '{sid}'/{r.get('seq_id','?')}: amount must be "
                "positive. A species that was not added is simply absent "
                "from the file - a zero is a claim that it was added"
            )
            continue
        units = (r.get("amount_units") or "").lower()
        if units not in MOCK_UNITS:
            diag.error(
                f"mocks '{sid}': amount_units '{units}' is not one of "
                + ", ".join(sorted(MOCK_UNITS))
            )
            continue
        by_sample.setdefault(sid, []).append(
            {
                "seq_id": r.get("seq_id", ""),
                "amount": amount,
                "units": units,
                "mix": r.get("mix", ""),
                "replicate": r.get("replicate", ""),
                "pcr_batch": r.get("pcr_batch", ""),
            }
        )
    return by_sample


def check_mocks(
    mocks: dict, entries: list, samples: list, diag: Diagnostics
) -> dict:
    """Is this mock design capable of measuring what it is for?

    A mock series exists to measure amplification efficiency. Several ways
    of laying one out produce data that cannot, and all of them look
    reasonable on the bench - so they are caught here, while the plate is
    still a spreadsheet.
    """
    known = {e.seq_id: e for e in entries}
    quants = {e.seq_id for e in entries if e.type == "quant"}
    declared_mocks = {
        s.sample_id for s in samples if s.sample_type == "mock"
    } if samples else set()

    by_mix: dict[str, set] = {}
    batches: set = set()
    series: dict[tuple, set] = {}

    for sid, rows in sorted(mocks.items()):
        units = {r["units"] for r in rows}
        if len(units) > 1:
            diag.error(
                f"mocks '{sid}': mixes {', '.join(sorted(units))} in one "
                "sample. The rows of a well have to be on one scale or "
                "their ratios mean nothing"
            )
        for r in rows:
            if r["seq_id"] not in known:
                diag.error(
                    f"mocks '{sid}': '{r['seq_id']}' is not in the "
                    "reference panel"
                )
            batches.add(r["pcr_batch"])
        if declared_mocks and sid not in declared_mocks:
            diag.error(
                f"mocks '{sid}': not in the samplesheet as sample_type "
                "'mock'"
            )

        if "fraction" in units:
            tot = sum(r["amount"] for r in rows if r["seq_id"] not in quants)
            if abs(tot - 1.0) > 0.01:
                diag.error(
                    f"mocks '{sid}': target fractions sum to {tot:.3f}, "
                    "not 1"
                )

        targets = [r for r in rows if r["seq_id"] not in quants]
        if not any(r["seq_id"] in quants for r in rows):
            diag.warn(
                f"mocks '{sid}': no quant standard. corr_factor is defined "
                "relative to the standards, so efficiencies from this well "
                "float against each other rather than anchoring to the "
                "scale field samples are measured on"
            )
        mix = rows[0]["mix"] or sid
        by_mix.setdefault(mix, set()).add(sid)
        if len(targets) == 2:
            key = tuple(sorted(r["seq_id"] for r in targets))
            lo, hi = sorted(r["amount"] for r in targets)
            series.setdefault(key, set()).add(round(hi / lo, 4) if lo else 0)

    # A single mixture ratio cannot show whether efficiency is constant.
    for pair, ratios in sorted(series.items()):
        if len(ratios) < 3:
            diag.warn(
                f"mocks: {pair[0]} vs {pair[1]} appears at "
                f"{len(ratios)} ratio(s). At least three are needed to fit "
                "a slope, and the slope is what says whether a single "
                "corr_factor per species is valid at all. With one ratio "
                "the correction is right only at that ratio"
            )
        elif max(ratios) / max(min(ratios), 1e-9) < 20:
            diag.warn(
                f"mocks: {pair[0]} vs {pair[1]} spans a narrow range of "
                "ratios. Extending to about 99:1 also measures the "
                "minority limit of quantification at no extra cost"
            )

    single_well = [m for m, wells in by_mix.items() if len(wells) < 2]
    if single_well:
        diag.warn(
            f"{len(single_well)} mixture(s) have no replicate: "
            + ", ".join(sorted(single_well)[:4])
            + ("..." if len(single_well) > 4 else "")
            + ". Without replicates there is no confidence interval on the "
            "efficiency, and a 1.3x correction cannot be told from noise"
        )
    batches.discard("")
    if mocks and len(batches) < 2:
        diag.warn(
            "every mock was amplified in one PCR batch. Replicate wells "
            "then measure pipetting rather than amplification, and a bad "
            "plate is indistinguishable from a real efficiency difference. "
            "Spread the replicates across separate PCR occasions"
        )

    # Every species that will need correcting has to appear somewhere.
    seen = {r["seq_id"] for rows in mocks.values() for r in rows}
    missing = [
        e.seq_id for e in entries
        if e.type == "target" and e.seq_id not in seen
    ]
    if mocks and missing:
        diag.warn(
            f"{len(missing)} target(s) appear in no mock, so no efficiency "
            "can be measured for them and they will keep the default "
            "corr_factor of 1.0: "
            + ", ".join(missing[:6]) + ("..." if len(missing) > 6 else "")
        )
    return {
        "n_mock_samples": len(mocks),
        "n_mixtures": len(by_mix),
        "n_pcr_batches": len(batches),
        "pairs_with_series": {
            f"{a}|{b}": sorted(r) for (a, b), r in series.items()
        },
        "targets_without_mock": missing,
    }


# --------------------------------------------------------------------------
# combined reference FASTA
# --------------------------------------------------------------------------


def write_insert_fasta(
    entries: list[RefEntry], extent: list[dict], outdir: str, diag: Diagnostics
) -> str | None:
    """The panel trimmed to between the primers.

    Read prep removes the primers, because those bases are synthesised, not
    observed - whatever the template carried under a primer, the read
    reports the primer. Leaving them on the reference while removing them
    from the reads would make every read look truncated; leaving them on
    both would silently assert template sequence that was never measured.

    Entries with no primer match (the vector backbone, which has no primer
    sites and is there to catch read-through) are carried over whole.
    """
    if not extent:
        return None
    bounds = {f["seq_id"]: (f["insert_start"], f["insert_end"]) for f in extent}

    seen: dict[str, str] = {}
    collisions: list[tuple[str, str]] = []
    path = os.path.join(outdir, "reference.insert.fasta")
    with open(path, "w", encoding="utf-8") as fh:
        for e in entries:
            if not e.sequence:
                continue
            if e.seq_id in bounds:
                a, b = bounds[e.seq_id]
                seq = e.sequence[a:b]
            else:
                seq = e.sequence
            if e.type in {"target", "control", "quant"}:
                prior = seen.get(seq)
                if prior:
                    collisions.append((prior, e.seq_id))
                else:
                    seen[seq] = e.seq_id
            fh.write(
                f">{e.seq_id} locus={e.locus or 'NA'} type={e.type} "
                f"taxon={e.report_taxon}\n"
            )
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")

    for a, b in collisions:
        diag.error(
            f"references '{a}' and '{b}' are different sequences but become "
            "identical once the primers are trimmed, so they differ only "
            "under a primer binding site. Whatever the template carried "
            "there, the read reports the primer - these two cannot be told "
            "apart by this assay"
        )
    return path


def write_combined_fasta(entries: list[RefEntry], outdir: str) -> str:
    """One index for everything. The point of competitive mapping is that a
    read chooses between targets, decoys, quants and vector."""
    path = os.path.join(outdir, "reference.fasta")
    with open(path, "w", encoding="utf-8") as fh:
        for e in entries:
            if not e.sequence:
                continue
            fh.write(
                f">{e.seq_id} locus={e.locus or 'NA'} type={e.type} "
                f"taxon={e.report_taxon}\n"
            )
            for i in range(0, len(e.sequence), 60):
                fh.write(e.sequence[i : i + 60] + "\n")
    return path


# --------------------------------------------------------------------------
# index map, plate geometry
# --------------------------------------------------------------------------

ROWS_384 = [chr(ord("A") + i) for i in range(16)]


@dataclass
class IndexMap:
    i7: dict[str, tuple[str, int]] = field(default_factory=dict)  # seq -> (name, col)
    i5: dict[str, tuple[str, str]] = field(default_factory=dict)  # seq -> (name, row)

    def well(self, i7_seq: str, i5_seq: str) -> str | None:
        """Well label in the lab's own convention: unpadded, A1 not A01.

        Padding would sort more nicely, but the label exists to be joined
        against plate records written by hand, and A01 silently fails to
        match A1. plate_row and plate_col are emitted separately for
        anything that needs to sort or do geometry.
        """
        if i7_seq not in self.i7 or i5_seq not in self.i5:
            return None
        _, col = self.i7[i7_seq]
        _, row = self.i5[i5_seq]
        return f"{row}{col}"

    def row_col(self, i7_seq: str, i5_seq: str) -> tuple[str, int] | None:
        if i7_seq not in self.i7 or i5_seq not in self.i5:
            return None
        return self.i5[i5_seq][1], self.i7[i7_seq][1]


def parse_index_map(path: str, diag: Diagnostics) -> IndexMap:
    rows = read_table(path, diag, "index map")
    if not require_columns(
        rows, {"role", "name", "sequence", "position"}, diag, "index map"
    ):
        return IndexMap()

    im = IndexMap()
    for row in rows:
        role = row["role"].lower()
        seq = row["sequence"].upper()
        name = row["name"]
        pos = row["position"]
        if role == "i7":
            try:
                im.i7[seq] = (name, int(pos))
            except ValueError:
                diag.error(f"index map: i7 '{name}' position '{pos}' is not a column number")
        elif role == "i5":
            if pos not in ROWS_384:
                diag.error(f"index map: i5 '{name}' position '{pos}' is not a row A-P")
            else:
                im.i5[seq] = (name, pos)
        else:
            diag.error(f"index map: unknown role '{role}' (expected i5 or i7)")
    return im


def index_set_report(im: IndexMap, diag: Diagnostics) -> dict:
    """Distances within and between the index sets, which cap how much
    mismatch tolerance the demultiplexer may safely allow."""

    def within(seqs: list[str], label: str) -> dict:
        if len(seqs) < 2:
            return {}
        pairs = [(a, b) for i, a in enumerate(seqs) for b in seqs[i + 1 :]]
        mh = min(hamming(a, b) for a, b in pairs)
        ml = min(levenshtein(a, b) for a, b in pairs)
        # tolerance that still guarantees unique assignment
        safe_mm = (mh - 1) // 2
        safe_edit = (ml - 1) // 2
        if safe_edit < safe_mm:
            diag.warn(
                f"{label}: minimum edit distance is {ml} but minimum Hamming "
                f"distance is {mh}. The set tolerates {safe_mm} substitution(s) "
                f"but only {safe_edit} indel(s) - so barcode matching must "
                "disable indels and rely on the tandem duplication for "
                "indel robustness"
            )
        return {
            "n": len(seqs),
            "min_hamming": mh,
            "min_levenshtein": ml,
            "safe_substitutions": safe_mm,
            "safe_edits": safe_edit,
        }

    rep = {
        "i7": within(list(im.i7), "index set i7"),
        "i5": within(list(im.i5), "index set i5"),
    }

    # Cross-set collisions: if an i5 and an i7 are near-identical, an
    # unanchored search cannot tell which end of the molecule it is reading.
    cross = []
    for a in im.i5:
        for b in im.i7:
            d = hamming(a, b) if len(a) == len(b) else levenshtein(a, b)
            cross.append((d, im.i5[a][0], a, im.i7[b][0], b))
    if cross:
        cross.sort()
        d, n5, s5, n7, s7 = cross[0]
        rep["cross_min_distance"] = d
        rep["cross_closest"] = f"{n5} {s5} vs {n7} {s7}"
        if d <= 2:
            diag.warn(
                f"index sets are only {d} substitution(s) apart at their closest "
                f"({n5} {s5} vs {n7} {s7}) - barcode search must be anchored on "
                "adapter context, never scanned across the read, or the two ends "
                "cannot be told apart"
            )

    # Rotational collisions: an indel upstream shifts the window, and a
    # rotated barcode landing near another is a misassignment route. The
    # requirement that both tandem copies agree is what catches this.
    rot = []
    for label, table in (("i7", im.i7), ("i5", im.i5)):
        for a in table:
            for r in rotations(a):
                for b in table:
                    if b == a:
                        continue
                    if hamming(r, b) <= 1:
                        rot.append(f"{label}: {table[a][0]} rotated is 1 mismatch from {table[b][0]}")
    rep["rotational_collisions"] = sorted(set(rot))
    return rep


def barcode_strategy_report(
    idx_rep: dict, runs: dict[str, Run], diag: Diagnostics
) -> dict:
    """Recommend, per run, what the demultiplexer may safely do - and say
    plainly what a non-duplicated barcode gives up.

    The tandem duplication is not extra error tolerance; it is indel
    insurance. An indel makes one copy uninformative and the other carries
    the read, and a frame slip yields two copies that disagree, which is
    what catches it. With a single copy there is no second observation, so
    a slipped frame lands wherever it lands.
    """
    n_rot = len(idx_rep.get("rotational_collisions", []))
    safe_sub = min(
        idx_rep.get("i7", {}).get("safe_substitutions", 0),
        idx_rep.get("i5", {}).get("safe_substitutions", 0),
    )
    out: dict[str, dict] = {}

    for run in runs.values():
        reps = run.barcode_repeats
        strategy = {
            "barcode_repeats": reps,
            "max_substitutions_per_copy": safe_sub,
            "allow_indels": False,
            "require_copy_agreement": reps >= 2,
            "frame_slip_protection": "copy agreement" if reps >= 2 else "none",
        }
        out[run.run_id] = strategy

        if run.platform != "nanopore":
            # Illumina arrives pre-split; none of this applies to the
            # sequencer's own demultiplexing.
            strategy["applies"] = False
            continue
        strategy["applies"] = True

        if reps == 1:
            msg = (
                f"run '{run.run_id}' is nanopore with barcode_repeats=1, so there "
                "is no second copy to cross-check. The tandem duplication is indel "
                "insurance rather than extra tolerance: with one copy, an indel "
                "upstream shifts the read frame and nothing detects it"
            )
            if n_rot:
                msg += (
                    f". {n_rot} rotational collision(s) exist in these index sets, "
                    "and in this mode they are live misassignment routes rather "
                    "than caught ones"
                )
            diag.warn(msg)
            diag.warn(
                f"run '{run.run_id}': consider anchoring the barcode window to an "
                "exact offset from the adapter rather than a tolerant one, and "
                "weigh exact 8bp matching (about 85% retention at 2% error, high "
                "specificity) against <=1 substitution (about 99% retention, "
                "vulnerable to the rotations above)"
            )
        elif reps > 2:
            diag.warn(
                f"run '{run.run_id}': barcode_repeats={reps} is unusual - confirm "
                "the construct really carries that many copies, since the search "
                "pattern is built from this number"
            )

        if safe_sub == 0:
            diag.warn(
                f"run '{run.run_id}': the index sets do not safely tolerate even "
                "one substitution, so barcode matching must require exact copies"
            )

    return out


# --------------------------------------------------------------------------
# run manifest and samplesheet
# --------------------------------------------------------------------------


@dataclass
class Run:
    run_id: str
    platform: str
    path: str
    barcode_repeats: int = 2
    dna_input_ng: str = ""
    lots: dict[str, str] = field(default_factory=dict)


def parse_runs(path: str, diag: Diagnostics) -> dict[str, Run]:
    rows = read_table(path, diag, "run manifest")
    if not require_columns(rows, {"run_id", "platform", "path"}, diag, "run manifest"):
        return {}

    runs: dict[str, Run] = {}
    for row in rows:
        rid = row["run_id"]
        if not rid:
            diag.error("run manifest: run_id is empty")
            continue
        if rid in runs:
            diag.error(f"run manifest: run_id '{rid}' appears more than once")
            continue

        platform = row["platform"].lower()
        if platform not in PLATFORMS:
            diag.error(
                f"run manifest '{rid}': platform '{platform}' not in {sorted(PLATFORMS)}"
            )

        # 1 = a conventional single 8bp index at each end; 2 = the tandem
        # duplication. Both are dual-index; they differ in whether each index
        # carries its own internal replicate.
        reps = row.get("barcode_repeats", "") or "2"
        try:
            reps_i = int(reps)
        except ValueError:
            diag.error(f"run manifest '{rid}': barcode_repeats '{reps}' is not an integer")
            reps_i = 2
        else:
            if reps_i < 1:
                diag.error(f"run manifest '{rid}': barcode_repeats must be >= 1")

        if not os.path.exists(row["path"]):
            diag.error(f"run manifest '{rid}': path does not exist: {row['path']}")

        dna = row.get("dna_input_ng", "")
        if dna:
            try:
                if float(dna) <= 0:
                    diag.error(f"run manifest '{rid}': dna_input_ng must be positive")
            except ValueError:
                diag.error(f"run manifest '{rid}': dna_input_ng '{dna}' is not a number")
        else:
            diag.warn(
                f"run manifest '{rid}': dna_input_ng is not set - copy-number "
                "estimates can be reported per reaction but not per unit DNA"
            )

        lots = {k: v for k, v in row.items() if k.startswith("lot_") and v}
        if not lots:
            diag.warn(
                f"run manifest '{rid}': no lot_* fields recorded. Lots cannot be "
                "reconstructed later, and without them a step change in any "
                "cross-run comparison is uninterpretable"
            )

        runs[rid] = Run(
            run_id=rid,
            platform=platform,
            path=row["path"],
            barcode_repeats=reps_i,
            dna_input_ng=dna,
            lots=lots,
        )
    return runs


@dataclass
class Sample:
    run_id: str
    platform: str
    sample_id: str
    index_i7: str
    index_i5: str
    sample_type: str
    file_prefix: str = ""
    well: str | None = None
    plate_row: str | None = None
    plate_col: int | None = None
    files: list[str] = field(default_factory=list)


def parse_samplesheet(
    path: str, runs: dict[str, Run], im: IndexMap, diag: Diagnostics
) -> list[Sample]:
    rows = read_table(path, diag, "samplesheet")
    needed = {"run_id", "platform", "sample_id", "index_i7", "index_i5", "sample_type"}
    if not require_columns(rows, needed, diag, "samplesheet"):
        return []

    samples: list[Sample] = []
    for row in rows:
        s = Sample(
            run_id=row["run_id"],
            platform=row["platform"].lower(),
            sample_id=row["sample_id"],
            index_i7=row["index_i7"].upper(),
            index_i5=row["index_i5"].upper(),
            sample_type=row["sample_type"].lower(),
            file_prefix=row.get("file_prefix", ""),
        )
        if not s.sample_id:
            diag.error("samplesheet: sample_id is empty")
            continue
        samples.append(s)

    # -- identifier uniqueness within a run ------------------------------
    per_run: dict[str, list[str]] = defaultdict(list)
    for s in samples:
        per_run[s.run_id].append(s.sample_id)
    for rid, ids in per_run.items():
        for sid, n in Counter(ids).items():
            if n > 1:
                diag.error(
                    f"samplesheet: sample_id '{sid}' appears {n} times in run '{rid}'"
                )

    # -- sanitisation collisions -----------------------------------------
    tokens: dict[tuple[str, str], set[str]] = defaultdict(set)
    for s in samples:
        tokens[(s.run_id, sanitise_id(s.sample_id))].add(s.sample_id)
    for (rid, token), originals in tokens.items():
        if len(originals) > 1:
            diag.error(
                f"samplesheet run '{rid}': sample_ids {sorted(originals)} all "
                f"sanitise to '{token}' - path names would collide"
            )

    # -- platform uniformity within a run --------------------------------
    plats: dict[str, set[str]] = defaultdict(set)
    for s in samples:
        plats[s.run_id].add(s.platform)
    for rid, pset in plats.items():
        if len(pset) > 1:
            diag.error(
                f"samplesheet run '{rid}': mixed platforms {sorted(pset)}. A run is "
                "one platform; a mixed run is always a metadata error"
            )
        bad = pset - PLATFORMS
        if bad:
            diag.error(f"samplesheet run '{rid}': unknown platform(s) {sorted(bad)}")

    # -- run references resolve ------------------------------------------
    for rid in per_run:
        if rid not in runs:
            diag.error(
                f"samplesheet: run_id '{rid}' has no entry in the run manifest"
            )
    for rid in runs:
        if rid not in per_run:
            diag.warn(f"run manifest: run '{rid}' has no samples in the samplesheet")

    # -- platform agreement with the manifest ----------------------------
    for s in samples:
        run = runs.get(s.run_id)
        if run and run.platform and s.platform and run.platform != s.platform:
            diag.error(
                f"samplesheet '{s.sample_id}': platform '{s.platform}' disagrees "
                f"with run manifest '{s.run_id}' platform '{run.platform}'"
            )

    # -- sample types -----------------------------------------------------
    for s in samples:
        if s.sample_type not in SAMPLE_TYPES:
            diag.error(
                f"samplesheet '{s.sample_id}': sample_type '{s.sample_type}' not in "
                f"{sorted(SAMPLE_TYPES)}"
            )

    for rid, ids in per_run.items():
        types = {s.sample_type for s in samples if s.run_id == rid}
        if "neg_extraction" not in types:
            diag.warn(
                f"run '{rid}': no neg_extraction sample. No empirical detection "
                "threshold can be derived, so calls will fall back to the "
                "configured default and be flagged as such"
            )
        if "mock" not in types:
            diag.warn(
                f"run '{rid}': no mock sample - discrimination accuracy cannot be "
                "measured for this run"
            )

    # -- indices ----------------------------------------------------------
    for s in samples:
        for label, seq in (("index_i7", s.index_i7), ("index_i5", s.index_i5)):
            if not seq:
                if s.platform == "nanopore":
                    diag.error(
                        f"samplesheet '{s.sample_id}': {label} is required for "
                        "nanopore - the index pair is the demultiplexing key, "
                        "not metadata"
                    )
                else:
                    diag.warn(
                        f"samplesheet '{s.sample_id}': {label} is not recorded. "
                        "Illumina data arrives pre-split so the run will proceed, "
                        "but misassignment diagnostics need the index pair"
                    )
                continue
            bad = set(seq) - IUPAC
            if bad:
                diag.error(
                    f"samplesheet '{s.sample_id}': {label} contains non-IUPAC "
                    f"character(s): {''.join(sorted(bad))}"
                )
            table = im.i7 if label == "index_i7" else im.i5
            if table and seq not in table:
                diag.error(
                    f"samplesheet '{s.sample_id}': {label} '{seq}' is not in the "
                    "index map, so its plate position cannot be derived"
                )

        s.well = im.well(s.index_i7, s.index_i5)
        rc = im.row_col(s.index_i7, s.index_i5)
        if rc:
            s.plate_row, s.plate_col = rc

    # -- one sample per well ---------------------------------------------
    wells: dict[tuple[str, str], list[str]] = defaultdict(list)
    for s in samples:
        if s.well:
            wells[(s.run_id, s.well)].append(s.sample_id)
    for (rid, well), ids in wells.items():
        if len(ids) > 1:
            diag.error(
                f"run '{rid}': index pair resolves to well {well} for more than "
                f"one sample ({', '.join(ids)}) - the pair must be unique per run"
            )

    return samples


# --------------------------------------------------------------------------
# plate geometry and sentinels
# --------------------------------------------------------------------------


def quadrant_of(rows_used: set[str], cols_used: set[int]) -> str:
    """A 96-well plate occupies one of four interleaved quadrants of the 384
    map. Inferred from parity, so it works whichever quadrant was used."""
    row_idx = {ROWS_384.index(r) + 1 for r in rows_used}
    row_par = {i % 2 for i in row_idx}
    col_par = {c % 2 for c in cols_used}
    if len(row_par) > 1 or len(col_par) > 1:
        return "none (spans quadrants)"
    rp, cp = row_par.pop(), col_par.pop()
    return {
        (1, 1): "Q1 (odd rows, odd columns)",
        (1, 0): "Q2 (odd rows, even columns)",
        (0, 1): "Q3 (even rows, odd columns)",
        (0, 0): "Q4 (even rows, even columns)",
    }[(rp, cp)]


def plate_report(
    samples: list[Sample], im: IndexMap, diag: Diagnostics
) -> dict[str, dict]:
    out: dict[str, dict] = {}
    by_run: dict[str, list[Sample]] = defaultdict(list)
    for s in samples:
        by_run[s.run_id].append(s)

    for rid, group in by_run.items():
        placed = [s for s in group if s.well]
        if not placed:
            continue

        i5_used = {s.index_i5 for s in placed}
        i7_used = {s.index_i7 for s in placed}
        rows_used = {im.i5[s][1] for s in i5_used if s in im.i5}
        cols_used = {im.i7[s][1] for s in i7_used if s in im.i7}

        combos = len(i5_used) * len(i7_used)
        assigned = len({(s.index_i5, s.index_i7) for s in placed})

        # A sentinel is an unassigned COMBINATION of two indices that each
        # appear elsewhere in the run. Omitting a whole index measures
        # nothing, because nothing can hop into an index that was never
        # present in the pool.
        sentinels = [
            (im.i5[a][1], im.i7[b][1])
            for a, b in product(sorted(i5_used), sorted(i7_used))
            if (a, b) not in {(s.index_i5, s.index_i7) for s in placed}
            and a in im.i5
            and b in im.i7
        ]
        observable = len(sentinels) / combos if combos else 0.0

        quad = quadrant_of(rows_used, cols_used) if len(placed) <= 96 else "n/a (>96 samples)"

        out[rid] = {
            "n_samples": len(group),
            "n_placed": len(placed),
            "i5_used": len(i5_used),
            "i7_used": len(i7_used),
            "combinations": combos,
            "assigned": assigned,
            "sentinels": len(sentinels),
            "observable_fraction": round(observable, 4),
            "quadrant": quad,
        }

        if observable == 0:
            diag.warn(
                f"run '{rid}': every index combination in the pool is assigned, so "
                "no sentinel exists and index misassignment is UNMEASURABLE for "
                "this run - which is a different claim from measuring it to be "
                "zero. Leaving a few wells unused would recover an estimate"
            )
        # A positive finding is not a warning; it belongs in the summary.
    return out


# --------------------------------------------------------------------------
# Illumina file resolution
# --------------------------------------------------------------------------


def resolve_illumina_files(
    samples: list[Sample], runs: dict[str, Run], diag: Diagnostics
) -> None:
    """Resolve and print every file before any compute starts. A mismatch
    between samplesheet names and filenames is the commonest way this class
    of pipeline wastes a day, so zero matches and ambiguous matches both
    hard-fail rather than warn."""
    for s in samples:
        run = runs.get(s.run_id)
        if not run or run.platform != "illumina":
            continue
        if not os.path.isdir(run.path):
            continue  # already reported by the manifest check

        prefix = s.file_prefix or s.sample_id
        for read in ("R1", "R2"):
            pattern = os.path.join(run.path, f"{prefix}_S*_L*_{read}_*.fastq.gz")
            hits = sorted(glob.glob(pattern))
            if not hits:
                # try the simpler naming some centres use
                alt = os.path.join(run.path, f"{prefix}*_{read}*.fastq.gz")
                hits = sorted(glob.glob(alt))
            if not hits:
                diag.error(
                    f"samplesheet '{s.sample_id}': no {read} file matched under "
                    f"{run.path} (tried '{os.path.basename(pattern)}'). If the "
                    "sequencing centre used a different name, set file_prefix"
                )
            else:
                s.files.extend(hits)
                lanes = {re.search(r"_L(\d+)_", os.path.basename(h)) for h in hits}
                lanes = {m.group(1) for m in lanes if m}
                if len(hits) > 1 and len(lanes) < len(hits):
                    diag.error(
                        f"samplesheet '{s.sample_id}': {read} matched {len(hits)} "
                        f"files that are not distinct lanes: "
                        f"{[os.path.basename(h) for h in hits]}. An ambiguous match "
                        "is as dangerous as no match"
                    )


def write_file_resolution(samples: list[Sample], outdir: str) -> None:
    path = os.path.join(outdir, "file_resolution.tsv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            ["run_id", "sample_id", "sanitised_id", "platform", "well",
             "plate_row", "plate_col", "sample_type", "n_files", "files"]
        )
        for s in samples:
            w.writerow(
                [
                    s.run_id,
                    s.sample_id,
                    sanitise_id(s.sample_id),
                    s.platform,
                    s.well or "",
                    s.plate_row or "",
                    s.plate_col or "",
                    s.sample_type,
                    len(s.files),
                    # Full paths, not basenames. This table exists so that
                    # every downstream step uses a file pre-flight has
                    # already checked, rather than re-deriving one - and a
                    # basename cannot be opened. It stayed unnoticed until
                    # something actually tried to read from it.
                    ";".join(os.path.abspath(f) for f in s.files),
                ]
            )


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Pre-flight validation for nf-helminth-amplicon."
    )
    ap.add_argument("--reference", required=True)
    ap.add_argument(
        "--samplesheet",
        help="omit to validate the reference panel on its own, which is "
        "useful while building it",
    )
    ap.add_argument("--runs", help="run manifest TSV")
    ap.add_argument("--input-dir", help="single-run shortcut; builds a one-row manifest")
    ap.add_argument("--run-id", default="RUN01", help="run_id for --input-dir")
    ap.add_argument("--platform", choices=sorted(PLATFORMS), help="platform for --input-dir")
    ap.add_argument("--index-map")
    ap.add_argument(
        "--mocks",
        help=(
            "mock community definitions - what was put into each mock well, "
            "declared before sequencing (see assets/mocks_template.tsv)"
        ),
    )
    ap.add_argument(
        "--illumina-read-length",
        type=int,
        default=250,
        help=(
            "read length per mate, used to check that the control tag is "
            "observable on pairs that cannot merge (default: 250)"
        ),
    )
    ap.add_argument(
        "--primers",
        help="optional primer TSV (name, direction, sequence). Used only "
        "to validate that every reference entry spans exactly "
        "primer-to-primer - extent cannot be checked any other way",
    )
    ap.add_argument(
        "--assay",
        help="assay definition JSON (assets/assay/*.json). Holds the control "
        "tag, library construct, indexing design and tolerances. Built-in "
        "defaults are used when omitted",
    )
    ap.add_argument(
        "--barcode-repeats",
        type=int,
        help="tandem copies of each index in the construct; 2 is the default "
        "AMPSEQ design, 1 a conventional single index at each end. Only used "
        "when building a manifest from --input-dir; a --runs manifest sets it "
        "per run",
    )
    ap.add_argument(
        "--tag-max-edits",
        type=int,
        help="override the assay config's control_tag.max_edits",
    )
    ap.add_argument(
        "--expected-tag",
        help="override the assay config's control_tag.sequence; pass an "
        "empty string to skip the cross-check",
    )
    ap.add_argument("--outdir", default="preflight")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    diag = Diagnostics()

    # The assay definition supplies defaults; an explicit flag overrides it,
    # so a one-off run does not require editing the config.
    assay = load_assay(args.assay, diag)
    tag_cfg = assay.section("control_tag")
    tag_max_edits = (
        args.tag_max_edits
        if args.tag_max_edits is not None
        else int(tag_cfg["max_edits"])
    )
    expected_tag = (
        args.expected_tag
        if args.expected_tag is not None
        else str(tag_cfg["sequence"])
    )

    reference_only = not args.samplesheet

    # -- run manifest: explicit, or built from --input-dir ---------------
    if reference_only:
        runs = {}
    elif args.runs:
        runs = parse_runs(args.runs, diag)
    elif args.input_dir:
        if not args.platform:
            diag.error("--input-dir requires --platform")
            runs = {}
        else:
            default_reps = int(assay.section("indexing")["barcode_repeats"])
            reps = (
                args.barcode_repeats
                if args.barcode_repeats is not None
                else default_reps
            )
            manifest = os.path.join(args.outdir, "runs.derived.tsv")
            with open(manifest, "w", encoding="utf-8") as fh:
                fh.write("run_id\tplatform\tpath\tbarcode_repeats\n")
                fh.write(
                    f"{args.run_id}\t{args.platform}\t{args.input_dir}\t{reps}\n"
                )
            runs = parse_runs(manifest, diag)
    else:
        diag.error("one of --runs or --input-dir is required")
        runs = {}

    # -- reference panel --------------------------------------------------
    entries = parse_reference(args.reference, diag)
    validate_reference(entries, diag, assay)
    if entries:
        write_combined_fasta(entries, args.outdir)
    primers = parse_primers(args.primers, diag) if args.primers else {}
    extent = check_reference_extent(entries, primers, diag) if primers else []
    if entries and extent:
        write_insert_fasta(entries, extent, args.outdir, diag)
    tags = analyse_tags(entries, diag, tag_max_edits, expected_tag, assay)
    write_tag_report(tags, args.outdir, tag_max_edits, assay)
    if tags:
        check_tag_visibility(entries, tags, diag, args.illumina_read_length)
        if extent:
            check_tag_in_insert(entries, tags, extent, diag)

    # -- indices, samples, geometry, files -------------------------------
    # Skipped entirely when validating a reference panel on its own, so the
    # panel can be checked before a samplesheet exists.
    im = IndexMap()
    idx_rep: dict = {}
    barcode_strategy: dict = {}
    samples: list[Sample] = []
    plates: dict = {}
    if not reference_only:
        if not args.index_map:
            diag.error("--index-map is required unless validating the reference alone")
        else:
            im = parse_index_map(args.index_map, diag)
            idx_rep = index_set_report(im, diag)
            barcode_strategy = barcode_strategy_report(idx_rep, runs, diag)
        samples = parse_samplesheet(args.samplesheet, runs, im, diag)
        plates = plate_report(samples, im, diag)
        resolve_illumina_files(samples, runs, diag)
        write_file_resolution(samples, args.outdir)

    # -- mock communities -------------------------------------------------
    # Checked against the panel and the samplesheet, so a design that
    # cannot measure efficiency is caught while the plate is still a
    # spreadsheet rather than after the reagents are spent.
    mock_report: dict = {}
    if args.mocks:
        mocks = parse_mocks(args.mocks, diag)
        if mocks:
            mock_report = check_mocks(mocks, entries, samples, diag)

    # -- machine-readable summary ----------------------------------------
    summary = {
        "schema_version": 1,
        "mode": "reference_only" if reference_only else "full",
        "mocks": mock_report,
        "assay": {
            "name": assay.name,
            "version": assay.version,
            "source": args.assay or "builtin defaults",
            "control_tag": expected_tag,
            "tag_max_edits": tag_max_edits,
        },
        "reference": {
            "n_entries": len(entries),
            "by_type": dict(Counter(e.type for e in entries)),
            "loci": sorted({e.locus for e in entries if e.locus}),
            "n_with_corr_factor": sum(1 for e in entries if e.corr_factor),
            "n_using_default_corr_factor": sum(
                1 for e in entries
                if e.type in {"target", "quant"} and not e.corr_factor
            ),
            "extent": extent,
            "composition_denominator": sorted(
                e.seq_id
                for e in entries
                if e.type == "target" and "speciation" in e.purpose
            ),
            "variant_references": sorted(
                e.seq_id for e in entries if "variants" in e.purpose
            ),
        },
        "tags": [
            {
                "control_id": t.control_id,
                "kind": t.kind,
                "tag": t.tag,
                "length": t.length,
                "start": t.start,
            }
            for t in tags
        ],
        "index_sets": idx_rep,
        "barcode_strategy": barcode_strategy,
        "runs": {
            r.run_id: {
                "platform": r.platform,
                "barcode_repeats": r.barcode_repeats,
                "dna_input_ng": r.dna_input_ng,
                "lots": r.lots,
            }
            for r in runs.values()
        },
        "plates": plates,
        "samples": {
            "n": len(samples),
            "by_type": dict(Counter(s.sample_type for s in samples)),
        },
        "errors": diag.errors,
        "warnings": diag.warnings,
    }
    with open(os.path.join(args.outdir, "preflight_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    report = diag.render()
    with open(os.path.join(args.outdir, "validation.log"), "w", encoding="utf-8") as fh:
        fh.write(report + "\n")

    print(report)
    print(f"\nReports written to {args.outdir}/")
    if reference_only:
        print(
            "\nReference panel checked on its own. Samplesheet, indices, "
            "plate geometry and file resolution were not examined."
        )

    if not diag.ok:
        print(
            f"\nPre-flight FAILED with {len(diag.errors)} error(s). "
            "No reads were processed.",
            file=sys.stderr,
        )
        return 1
    print("\nPre-flight passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
