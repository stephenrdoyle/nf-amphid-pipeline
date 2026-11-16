"""
Database-mode protocols: everything from primer trimming to a call.

Database mode has three layers (docs/CONTRACTS.md):

  platform   gets each sample's reads - Illumina pre-split by the
             sequencer, nanopore split and demultiplexed here, or any
             platform whose per-sample FASTQs already exist
             (predemultiplexed). It hands on CONTRACT A: per-sample reads,
             their read type, a demux accounting and plate positions.
  protocol   this package. One module per database protocol, chosen by
             the release: it trims, clusters or denoises, classifies and
             writes CONTRACT B (clusters.tsv, taxa.tsv, assign_db.json).
  shared     synthesis and the report, which read contract B and never
             ask which protocol wrote it. What they need to know of a
             protocol - organism groups, the host clade, the relative
             abundance floor - is a hook on Protocol, with a default.

A protocol module defines one Protocol subclass and exposes it as
PROTOCOL. It may use anything in ampseqlib and the read-type routes'
machinery (consensus building in assign_db.py, denoising in
assign_db_pairs.py), which call it; it never imports another protocol.

Which protocol a release gets: the manifest's "protocol", else inferred
from what the manifests already say (protocol_for), so every existing
release keeps working; --db_protocol overrides either.

Standard library only.
"""

from __future__ import annotations

import importlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ampseqlib import RANKS, orient_pair, parse_floors, ranks_of  # noqa: E402

# The read types a platform adapter can hand on (contract A).
READ_TYPES = {
    "paired_short": "Illumina pairs: two short, accurate mates per molecule",
    "single_long_noisy": "one read per molecule, ~1-5% error (nanopore)",
    "single_long_accurate": "one read per molecule, <0.1% error (PacBio HiFi)",
}

# Every protocol this pipeline carries, by the name a manifest or
# --db_protocol gives. Each is bin/protocols/<name>.py.
PROTOCOLS = ("vespa_v4", "nemits2")


def scores_over_reference(manifest: dict) -> bool:
    """Whether a release's records are scored over their own span, the
    query's overhang free (manifest `score_over`; nemabiome ITS2).

    ITS2 records stop where the ~47-base 28S stretch before NC2 starts,
    and many are partial, so a read runs past its record. vespa's V4
    records are whole primer-to-primer inserts with no `score_over`, and
    there an unaligned end of the query is a real difference."""
    return "reference" in str(manifest.get("score_over", "")).lower()


# The rules that need span scoring (protocols/nemits2.py).
SPAN_RULES = ("tiered", "gap")


def calling_rule(manifest: dict, override: str | None = None) -> str:
    """'floors' (an identity floor per rank, as for V4), 'gap' (ITS2:
    species only when the next species is several edits further off) or
    'tiered' (the gap rule's calls, plus a flagged tentative species at a
    lower identity with the same gap). An explicit choice wins; then the
    manifest's `calling_rule`; then a release scored over its records'
    span takes the tiered rule, its database's recommended default;
    everything else keeps the floors."""
    if override and override != "auto":
        return override
    given = str(manifest.get("calling_rule", "")).lower()
    if given in ("floors",) + SPAN_RULES:
        return given
    return "tiered" if scores_over_reference(manifest) else "floors"


def protocol_for(manifest: dict, override: str | None = None,
                 rule_override: str | None = None) -> str:
    """The protocol a release is run with.

    --db_protocol (override) wins; then the manifest's "protocol"; then
    what releases built before the field existed already say. A release
    scored over its records' span, or called by the gap or tiered rule,
    is nemabiome ITS2's protocol - span scoring is what those rules are
    measured over, so --db_calling_rule gap on a V4 release has always
    meant span scoring too, and still does. Everything else is vespa V4.
    ValueError for a name this pipeline does not carry."""
    name = override if override and override != "auto" else \
        str(manifest.get("protocol") or "")
    if not name:
        name = "nemits2" if (scores_over_reference(manifest) or calling_rule(
            manifest, rule_override) in SPAN_RULES) else "vespa_v4"
    if name not in PROTOCOLS:
        raise ValueError(
            f"protocol '{name}' is not one this pipeline carries "
            f"({', '.join(PROTOCOLS)})")
    return name


def get(name: str) -> "Protocol":
    """The Protocol object of a module in this package."""
    if name not in PROTOCOLS:
        raise ValueError(f"protocol '{name}' is not one this pipeline "
                         f"carries ({', '.join(PROTOCOLS)})")
    return importlib.import_module(f"protocols.{name}").PROTOCOL


def resolve(manifest: dict, override: str | None = None,
            rule_override: str | None = None) -> "Protocol":
    return get(protocol_for(manifest, override, rule_override))


# --------------------------------------------------------------------------
# organism groups: the default hook
# --------------------------------------------------------------------------

HELMINTH_CLASSES = {"Nematoda", "Platyhelminthes", "Acanthocephala"}


def organism_group(lineage: list[str], host_clade: str,
                   ranks: list[str] = RANKS) -> str:
    """The reporting group of a call, from its lineage down to the rank
    it reached (synthesise_db.py's docstring). 'unassigned' when that
    lineage stops too high to tell."""
    if ranks != RANKS:
        return organism_group_named(dict(zip(ranks, lineage)), host_clade)
    lin = lineage + [""] * (len(RANKS) - len(lineage))
    dom, sup, sub, cls, order = lin[0], lin[1], lin[3], lin[4], lin[5]
    if dom in ("Bacteria", "Archaea"):
        return "bacteria_archaea"
    if dom.startswith("Eukaryota:"):
        return "organelle"
    if dom != "Eukaryota":
        return "unassigned"
    if host_clade and host_clade in lin:
        return "host"
    if cls in HELMINTH_CLASSES:
        return "helminths"
    if order == "Microsporida":
        return "microsporidia"
    if sub == "Fungi":
        return "fungi"
    if sub == "Metazoa":
        # a Metazoa call that stopped above class could still be a
        # helminth: not placed
        return "invertebrates" if cls else "unassigned"
    if sup == "Archaeplastida":
        return "plants"
    if sup == "Obazoa":
        # Choanoflagellates, apusomonads and the like are protists; an
        # Obazoa call that stopped above subdivision could be an animal.
        return "protists" if sub else "unassigned"
    return "protists" if sup else "unassigned"


def organism_group_named(lin: dict, host_clade: str) -> str:
    """organism_group() for a release whose ranks are real ones (NCBI's
    superkingdom ... species): the same groups, read by rank name."""
    top = lin.get("superkingdom") or lin.get("domain") or ""
    phylum, kingdom = lin.get("phylum", ""), lin.get("kingdom", "")
    if top in ("Bacteria", "Archaea"):
        return "bacteria_archaea"
    if top != "Eukaryota":
        return "unassigned"
    if host_clade and host_clade in lin.values():
        return "host"
    if phylum in HELMINTH_CLASSES:
        return "helminths"
    if phylum == "Microsporidia":
        return "microsporidia"
    if kingdom == "Fungi":
        return "fungi"
    if kingdom == "Metazoa":
        return "invertebrates" if phylum else "unassigned"
    if kingdom == "Viridiplantae":
        return "plants"
    return "protists" if kingdom or phylum else "unassigned"


# --------------------------------------------------------------------------
# the protocol interface
# --------------------------------------------------------------------------


class Protocol:
    """What a protocol says of itself, and the hooks the read-type routes
    and synthesis call. Every hook has a default; a protocol overrides
    what it does differently.

    The routes: assign_db.py takes single reads (it maps, builds and
    splits consensus sequences, then asks classify_consensus for the
    calls); assign_db_pairs.py takes pairs (it orients, windows and
    denoises them, asking merge_mates and read_windows on the way, then
    classify_pairs for the calls). Both write contract B with the
    protocol's extra columns."""

    name = ""
    title = ""
    # Read type -> how this protocol handles it; and read type -> why it
    # does not. Every READ_TYPES key is in exactly one of the two.
    read_types: dict = {}
    not_supported: dict = {}
    calling_rules: tuple = ("floors",)
    # Release files the protocol reads, beyond the release format's own
    # (ampseqlib.REFDB_REQUIRED).
    required_files: tuple = ()
    # Optional contract-B columns this protocol adds to clusters.tsv,
    # after the required ones (docs/CONTRACTS.md).
    extra_columns: tuple = ()
    # Pair fates beyond the shared ones, and what they mean.
    pair_fates: tuple = ()
    pair_note = ""
    # Synthesis hooks.
    host_clade = "Craniata"
    min_rel_abundance = 0.001
    # The marker, as reports name it ("relative 18S abundance").
    marker = "18S"
    # Whether synthesis writes the sequences left above family to FASTA
    # (unplaced_its2.fasta, not_its2.fasta) and reports their shares.
    unplaced_fasta = False
    # Pair chimeras (assign_db.chimera_scan, assign_db.parents_differ):
    # candidate parents per sequence, and whether a parent named only
    # above genus is judged by its best record's genus.
    chimera_max_parents = 6
    chimera_parents_by_best_hit = False

    # -- what it can run ------------------------------------------------
    def unsupported(self, read_type: str) -> str:
        """'' when this protocol handles the read type, else why not."""
        if read_type in self.read_types:
            return ""
        why = self.not_supported.get(read_type) or (
            f"'{read_type}' is not a read type "
            f"({', '.join(READ_TYPES)})")
        return f"protocol {self.name} cannot take {read_type} reads: {why}"

    def release_problems(self, release: str, manifest: dict,
                         read_types=()) -> list[tuple[str, str]]:
        """(level, message) for each thing this protocol needs and the
        release or the run lacks; level is 'error' or 'warn'."""
        out = [("error", f"protocol {self.name} needs {f}, which the "
                "release does not have")
               for f in self.required_files
               if not os.path.isfile(os.path.join(release, f))]
        out += [("error", self.unsupported(rt)) for rt in read_types
                if self.unsupported(rt)]
        return out

    # -- calling --------------------------------------------------------
    def calling(self, args, manifest: dict, tax) -> tuple[str, dict]:
        """(rule, floors) for this release and these options."""
        rule = calling_rule(manifest, args.calling_rule)
        if rule not in self.calling_rules:
            raise SystemExit(
                f"protocol {self.name} calls by "
                f"{' or '.join(self.calling_rules)}, not '{rule}'. The gap "
                "rule needs span scoring: run with --db_protocol nemits2")
        return rule, parse_floors(args.floors, ranks_of(tax))

    def calling_settings(self, rule: str, args) -> dict:
        """What assign_db.json's settings record of how calls were made,
        beyond the route's own."""
        return {}

    def column_values(self, call: dict, tax) -> dict:
        """This protocol's extra_columns for one call."""
        return {}

    def classify_consensus(self, ctx, seqs: dict) -> dict:
        """{key: call} for consensus sequences {key: sequence}.

        ctx carries mm (ampseqlib.Minimap2), index, db_fasta, tmp, tax,
        ranks, rule, floors and args. A call is a dict with at least
        depth, rank, identity, coverage, best_hit, near_hits and
        lca_depth (assign_db.label_for reads them)."""
        raise NotImplementedError(
            f"protocol {self.name} does not classify single reads")

    # -- pairs ----------------------------------------------------------
    def orient(self, r1, r2, ps, fate) -> tuple | None:
        """((R1 seq, qual), (R2 seq, qual), fwd primer name, the sequence
        R1 runs into past the insert) with primers trimmed, R1 the
        forward-primer mate; None when the pair cannot be placed. By
        default a primer at the start of each mate (ampseqlib.orient_pair),
        and the opposite primer is what R1 runs through."""
        o = orient_pair(r1, r2, ps)
        if o is None:
            return None
        return o[0], o[1], o[2], dict(ps.rev).get(o[3], "")

    def tail(self, manifest: dict, read_length: int | None = None) -> int:
        """Bases a primer-trimmed R2 reads before its record's span
        starts: 0 for records that are whole inserts."""
        return 0

    def read_windows(self, manifest: dict, read_length, w1: int, w2: int,
                     primer_offset: int = 0) -> tuple[int, int, int]:
        """(read bases kept from R1, from R2, R2's tail) for windows w1,
        w2: R2 also takes its tail, and both give up the spacer bases a
        protocol may put before each primer."""
        tail = self.tail(manifest, read_length)
        return w1 - primer_offset, w2 + tail - primer_offset, tail

    def window_settings(self, rw1: int, rw2: int, args) -> dict:
        return ({"read_windows": [rw1, rw2],
                 "primer_offset": args.primer_offset}
                if args.primer_offset else {})

    def fold_variants(self, centroids: list, args) -> tuple[list, int]:
        """(centroids, reads folded): a second denoising pass over the
        denoised pairs, centroids[k] = (R1 part, R2 part or "", reads).
        None by default."""
        return centroids, 0

    def merge_mates(self, s1: str, q1: str, s2: str, q2: str, args):
        """(insert, quality) for mates to be scored as one merged insert,
        or None to keep them as two windows."""
        return None

    def classify_pairs(self, ctx, centroids: list) -> tuple[dict, int]:
        """({k: call}, reads remapped) for denoised pairs, centroids[k] =
        (R1 window, R2 window or '' for a whole insert, reads). ctx as for
        classify_consensus, plus database and tail."""
        raise NotImplementedError(
            f"protocol {self.name} does not classify pairs")

    # -- synthesis hooks -----------------------------------------------
    def organism_group(self, lineage: list[str], host_clade: str,
                       ranks: list[str] = RANKS) -> str:
        return organism_group(lineage, host_clade, ranks)
