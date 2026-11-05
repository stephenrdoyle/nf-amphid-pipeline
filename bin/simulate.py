#!/usr/bin/env python3
"""
Read simulator for amphid.

Emits FASTQ plus a per-read truth table: which sample a read came from,
which reference sequence, and what was done to it. Almost nothing downstream
can be verified without that. A read-accounting funnel can sum perfectly
while every read sits in the wrong bucket, so the invariant has to be
testable against *intent*, not merely against itself.

Consumes the same reference panel, samplesheet, run manifest and assay
definition as the pipeline, so the two cannot drift apart: a schema change
breaks both at once rather than leaving the tests validating a stale format.

Pathologies are injected deliberately and at known rates, which makes the QC
sensors testable too - you can ask whether the pipeline recovers the 0.5%
index hopping you put there on purpose.

Self-contained: standard library only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field

# --------------------------------------------------------------------------
# construct
# --------------------------------------------------------------------------

DEFAULTS = {
    "pcr1_overhang_fwd": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "pcr1_overhang_rev": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "p5_tail": "AATGATACGGCGACCACCGAGATCTACAC",
    "p7_tail": "CAAGCAGAAGACGGCATACGAGAT",
}

# Per-base error rates, and how much of that error is indel rather than
# substitution. Nanopore's signature failure is the homopolymer indel, which
# is exactly what a Hamming-based matcher cannot absorb.
ERROR_MODELS = {
    "illumina": {"rate": 0.001, "indel_fraction": 0.02, "homopolymer_bias": 1.0},
    "nanopore_sup": {"rate": 0.01, "indel_fraction": 0.45, "homopolymer_bias": 3.0},
    "nanopore_hac": {"rate": 0.02, "indel_fraction": 0.50, "homopolymer_bias": 3.5},
    "nanopore_fast": {"rate": 0.05, "indel_fraction": 0.55, "homopolymer_bias": 4.0},
}

FATES = [
    "clean",
    "hopped",
    "misread_i5",
    "misread_i7",
    "chimera",
    "carryover",
    "dimer",
    "offtarget",
    "truncated",
]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def read_table(path: str) -> list[dict[str, str]]:
    # The delimiter rule is duplicated from ampseqlib.delimiter_for, and
    # that is deliberate. This simulator is the oracle the pipeline is
    # tested against; importing the pipeline's own parsing would mean a
    # shared bug produced agreeing right and wrong answers and no test
    # could see it. The copy is held honest by a test that asserts the
    # two agree, rather than by a comment asking the next person to
    # remember.
    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    return [
        {k.strip(): (v or "").strip() for k, v in row.items()}
        for row in csv.DictReader(lines, delimiter=delim)
    ]


@dataclass
class Ref:
    seq_id: str
    locus: str
    type: str
    taxon: str
    parent_id: str
    sequence: str
    corr_factor: float = 1.0
    default_conc: float = 0.0


@dataclass
class Samp:
    run_id: str
    platform: str
    sample_id: str
    index_i7: str
    index_i5: str
    sample_type: str


def load_inputs(args) -> tuple[list[Ref], list[Samp], dict, dict]:
    refs = []
    for r in read_table(args.reference):
        refs.append(
            Ref(
                seq_id=r["seq_id"],
                locus=r.get("locus", ""),
                type=r["type"].lower(),
                taxon=r.get("taxon", "") or r["seq_id"],
                parent_id=r.get("parent_id", ""),
                sequence=r["sequence"].upper(),
                corr_factor=float(r["corr_factor"]) if r.get("corr_factor") else 1.0,
                default_conc=float(r["default_conc"]) if r.get("default_conc") else 0.0,
            )
        )

    samples = [
        Samp(
            run_id=s["run_id"],
            platform=s["platform"].lower(),
            sample_id=s["sample_id"],
            index_i7=s["index_i7"].upper(),
            index_i5=s["index_i5"].upper(),
            sample_type=s["sample_type"].lower(),
        )
        for s in read_table(args.samplesheet)
    ]

    runs = {}
    for r in read_table(args.runs):
        runs[r["run_id"]] = {
            "platform": r["platform"].lower(),
            "path": r["path"],
            "barcode_repeats": int(r.get("barcode_repeats") or 2),
        }

    assay = dict(DEFAULTS)
    if args.assay and os.path.exists(args.assay):
        with open(args.assay, encoding="utf-8") as fh:
            data = json.load(fh)
        lib = {k: v for k, v in data.get("library", {}).items() if not k.startswith("_")}
        assay.update({k: v for k, v in lib.items() if k in DEFAULTS})
    return refs, samples, runs, assay


# --------------------------------------------------------------------------
# composition: what is actually in each sample
# --------------------------------------------------------------------------


@dataclass
class Composition:
    """What a sample contains, before any sequencing happens.

    `counts` is input copies - the thing quantification is trying to
    recover. `weights` is what the PCR actually presents to the sequencer,
    which is copies scaled by each reference's amplification efficiency.
    Keeping the two apart is the whole point: if reads were drawn straight
    from copies, the corr_factor correction would have nothing to correct
    and would test as working while doing nothing.
    """

    counts: dict[str, float] = field(default_factory=dict)
    weights: dict[str, float] = field(default_factory=dict)
    note: str = ""

    def amplify(self, by_id: dict, competition: float = 0.0) -> None:
        """Copies in, read weights out.

        Two effects, and they are different in kind. Per-reference
        efficiency is a constant: a template that primes better takes a
        fixed multiple of the reads it is due, whatever else is present.
        Competition is not - templates contend for primers and polymerase,
        so as one comes to dominate the reaction its advantage is capped
        and the observed ratio compresses toward the middle.

        `competition` is the strength of that compression: 0 reproduces
        the simple multiplicative model, and 1 flattens the mixture
        entirely. It exists so the calibrator can be tested against data
        its own model does not describe - a fitter checked only against
        the linear case will report a confident constant for a reaction
        that does not have one.
        """
        raw = {
            k: v * (by_id[k].corr_factor if k in by_id else 1.0)
            for k, v in self.counts.items()
        }
        if competition <= 0 or len(raw) < 2:
            self.weights = raw
            return
        # Compress each share toward the mixture's own even split by
        # raising it to a power below one, then renormalise. A power of
        # (1 - competition) leaves order intact while pulling extremes in.
        total = sum(raw.values()) or 1.0
        power = max(1e-3, 1.0 - competition)
        shares = {k: (v / total) ** power for k, v in raw.items()}
        norm = sum(shares.values()) or 1.0
        self.weights = {k: total * v / norm for k, v in shares.items()}


def compose(
    rng: random.Random, sample: Samp, refs: list[Ref], args
) -> Composition:
    targets = [r for r in refs if r.type == "target"]
    controls = [r for r in refs if r.type == "control"]
    quants = [r for r in refs if r.type == "quant"]

    by_id = {r.seq_id: r for r in refs}
    comp = Composition()

    # A declared mock overrides everything: its composition was decided at
    # the bench, which is the whole point of it.
    declared = (args.mock_compositions or {}).get(sample.sample_id)
    if declared:
        comp.counts = dict(declared)
        comp.note = "mock: composition declared in the mocks file"
        comp.amplify(by_id, args.competition)
        return comp

    # Quants go into every sample, including blanks - which is what makes a
    # blank interpretable rather than merely empty.
    if not args.no_spike:
        for q in quants:
            # The spike is a fixed input, so its read share varies only with
            # how much else is competing and with pipetting error.
            comp.counts[q.seq_id] = q.default_conc * rng.lognormvariate(0, 0.15)

    if sample.sample_type == "neg_ntc":
        comp.note = "no template; quants only"
        comp.amplify(by_id, args.competition)
        return comp

    if sample.sample_type == "neg_extraction":
        comp.note = "extraction blank; quants only"
        comp.amplify(by_id, args.competition)
        return comp

    if sample.sample_type == "mock":
        # Mocks are built from control plasmids, so they contain no real
        # target sequence at all. Mock carryover therefore cannot masquerade
        # as a real positive - the route is closed by construction.
        for c in controls:
            comp.counts[c.seq_id] = 10 ** rng.uniform(2, 5)
        comp.note = "mock: control plasmids only, no real targets"
        comp.amplify(by_id, args.competition)
        return comp

    if sample.sample_type == "pos_control":
        for t in targets:
            comp.counts[t.seq_id] = 10 ** rng.uniform(3, 4)
        comp.note = "positive control: all targets at known level"
        comp.amplify(by_id, args.competition)
        return comp

    # Unknown samples: a random subset infected, at log-normal burden.
    # Prevalence and burden are what a real field set looks like - most
    # samples light or negative, a few very heavy.
    for t in targets:
        if rng.random() < args.prevalence:
            comp.counts[t.seq_id] = 10 ** rng.uniform(1, 5.5)
    comp.note = "field sample"
    comp.amplify(by_id, args.competition)
    return comp


# --------------------------------------------------------------------------
# sequencing error
# --------------------------------------------------------------------------


def R2_MODEL(model: dict) -> dict:
    """Read 2 is measurably worse than read 1 - the cluster has sat through
    read 1 and the index reads by the time it is imaged. Flat-rating both
    mates would make a merge look more trustworthy than it is."""
    return {**model, "rate": model["rate"] * 2.5}


def homopolymer_lengths(seq: str) -> list[int]:
    """Run length at each position, so indels can be biased toward
    homopolymers the way nanopore's actually are."""
    out = [1] * len(seq)
    run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        out[i] = run
    return out


def apply_errors(rng: random.Random, seq: str, model: dict) -> tuple[str, int, int]:
    rate = model["rate"]
    indel_frac = model["indel_fraction"]
    hp_bias = model["homopolymer_bias"]
    runs = homopolymer_lengths(seq)

    out = []
    n_sub = n_indel = 0
    for i, base in enumerate(seq):
        local = rate * (hp_bias if runs[i] >= 3 else 1.0)
        if rng.random() < local:
            if rng.random() < indel_frac:
                n_indel += 1
                if rng.random() < 0.5:
                    continue  # deletion
                out.append(base)  # insertion: duplicate
                out.append(base)
            else:
                n_sub += 1
                out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out), n_sub, n_indel


def qual_string(rng: random.Random, n: int, rate: float) -> str:
    """Phred scores consistent with the error model, so a quality filter in
    the pipeline sees something coherent rather than a flat line."""
    import math

    q = min(40, max(6, int(-10 * math.log10(max(rate, 1e-5)))))
    return "".join(chr(33 + max(2, min(40, int(rng.gauss(q, 3))))) for _ in range(n))


# --------------------------------------------------------------------------
# construct assembly
# --------------------------------------------------------------------------


def build_molecule(
    amplicon: str, i7: str, i5: str, repeats: int, assay: dict
) -> str:
    """Top strand of the finished library molecule.

    P5 - i5 x n - PCR1 fwd overhang - amplicon - rc(PCR1 rev overhang)
       - i7 x n - rc(P7)

    The i5 primer carries its index forward and the i7 primer carries the
    reverse complement, which the top strand reverses back - so both read in
    the same orientation as the index list.
    """
    return (
        assay["p5_tail"]
        + i5 * repeats
        + assay["pcr1_overhang_fwd"]
        + amplicon
        + revcomp(assay["pcr1_overhang_rev"])
        + i7 * repeats
        + revcomp(assay["p7_tail"])
    )


# --------------------------------------------------------------------------
# simulation
# --------------------------------------------------------------------------


@dataclass
class ReadRecord:
    read_id: str
    sample_id: str
    seq_id: str
    fate: str
    i7: str
    i5: str
    n_sub: int
    n_indel: int
    length: int
    orientation: str


def weighted_choice(rng: random.Random, weights: dict[str, float]) -> str:
    total = sum(weights.values())
    x = rng.uniform(0, total)
    acc = 0.0
    for k, w in weights.items():
        acc += w
        if x <= acc:
            return k
    return next(iter(weights))


def simulate_run(
    rng: random.Random,
    run_id: str,
    run: dict,
    samples: list[Samp],
    refs: list[Ref],
    assay: dict,
    args,
    model_name: str,
) -> tuple[list[ReadRecord], dict, dict]:
    by_id = {r.seq_id: r for r in refs}
    controls = [r for r in refs if r.type == "control"]
    vector = [r for r in refs if r.type == "vector"]
    repeats = run["barcode_repeats"]
    model = ERROR_MODELS[model_name]

    # index pairs actually in the pool, and which combinations are assigned
    pool_i7 = sorted({s.index_i7 for s in samples})
    pool_i5 = sorted({s.index_i5 for s in samples})
    assigned = {(s.index_i5, s.index_i7) for s in samples}
    sentinels = [
        (a, b) for a in pool_i5 for b in pool_i7 if (a, b) not in assigned
    ]

    records: list[ReadRecord] = []
    reads_by_sample: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    # What was put in, as opposed to what came out. Quantification has to be
    # checkable against this rather than against its own consistency.
    composition: dict[str, dict[str, float]] = {}
    counter = 0

    for sample in samples:
        # The composition belongs to the SAMPLE, not to the run. Nanopore
        # here sequences the same libraries Illumina did, so a subject must
        # carry the same true burden in both - otherwise a cross-platform
        # comparison is measuring two different people and any concordance
        # figure from it is meaningless. Deriving the stream from the
        # sample id rather than from run order also keeps it stable when
        # samples are reordered or a run is re-simulated alone.
        comp = compose(
            random.Random(
                int.from_bytes(
                    hashlib.sha256(
                        f"{args.seed}:{sample.sample_id}".encode()
                    ).digest()[:8],
                    "big",
                )
            ),
            sample, refs, args,
        )
        composition[sample.sample_id] = dict(comp.counts)
        if not comp.counts:
            continue

        # Depth follows what is in the sample: a well with no target yields
        # little but dimer, which is why read-count variance on an unspiked
        # run is expected rather than alarming.
        depth = int(args.reads_per_sample * rng.lognormvariate(0, 0.35))

        for _ in range(depth):
            counter += 1
            fate = "clean"
            i7, i5 = sample.index_i7, sample.index_i5
            owner = sample.sample_id

            roll = rng.random()
            if roll < args.dimer_rate:
                fate = "dimer"
                amplicon = (
                    assay["pcr1_overhang_fwd"][-15:]
                    + revcomp(assay["pcr1_overhang_rev"][-15:])
                )
                seq_id = "_dimer"
            elif roll < args.dimer_rate + args.carryover_rate and controls:
                # Plasmid carryover: a control sequence appearing in a sample
                # that never contained one. Unambiguous, because mocks carry
                # no real target sequence.
                fate = "carryover"
                pick = rng.choice(controls)
                amplicon = pick.sequence
                seq_id = pick.seq_id
            elif (
                roll < args.dimer_rate + args.carryover_rate + args.offtarget_rate
                and vector
            ):
                fate = "offtarget"
                pick = rng.choice(vector)
                start = rng.randrange(0, max(1, len(pick.sequence) - 300))
                amplicon = pick.sequence[start : start + 300]
                seq_id = pick.seq_id
            else:
                seq_id = weighted_choice(rng, comp.weights)
                amplicon = by_id[seq_id].sequence

            # Index hopping: the pair is rewritten, so the read lands in
            # another sample entirely. Rate is flat but the harm is not -
            # the same absolute number of reads devastates a 500-read
            # negative and is noise in a 50,000-read positive.
            if rng.random() < args.hop_rate:
                fate = "hopped"
                if sentinels and rng.random() < 0.5:
                    i5, i7 = rng.choice(sentinels)
                else:
                    i5 = rng.choice(pool_i5)
                    i7 = rng.choice(pool_i7)
                owner = next(
                    (
                        s.sample_id
                        for s in samples
                        if s.index_i5 == i5 and s.index_i7 == i7
                    ),
                    "_unassigned",
                )

            # Misassignment on a combinatorial plate: ONE index is misread,
            # so the read keeps the other and lands in the same row (i7
            # misread) or the same column (i5 misread) - unlike hopping
            # above, which rewrites the pair. Drawn only when asked for, so
            # a run simulated without it is the run it always was.
            if fate != "hopped" and (args.misread_i5_rate
                                     or args.misread_i7_rate):
                for end, rate, pool in (("i5", args.misread_i5_rate, pool_i5),
                                        ("i7", args.misread_i7_rate, pool_i7)):
                    if (rate and len(pool) > 1 and fate == "clean"
                            and rng.random() < rate):
                        fate = f"misread_{end}"
                        cur = i5 if end == "i5" else i7
                        new = rng.choice([x for x in pool if x != cur])
                        if end == "i5":
                            i5 = new
                        else:
                            i7 = new
                if fate.startswith("misread_"):
                    owner = next(
                        (s.sample_id for s in samples
                         if s.index_i5 == i5 and s.index_i7 == i7),
                        "_unassigned",
                    )

            mol = build_molecule(amplicon, i7, i5, repeats, assay)

            # Ligation chimera: two molecules joined, so the ends carry
            # indices from different samples. The tier-C check exists for
            # exactly this.
            if run["platform"] == "nanopore" and rng.random() < args.chimera_rate:
                fate = "chimera"
                other = rng.choice(samples)
                mol2 = build_molecule(
                    by_id[weighted_choice(rng, comp.weights)].sequence
                    if comp.weights
                    else amplicon,
                    other.index_i7,
                    other.index_i5,
                    repeats,
                    assay,
                )
                mol = mol[: len(mol) // 2] + mol2[len(mol2) // 2 :]

            if rng.random() < args.truncation_rate:
                fate = "truncated" if fate == "clean" else fate
                keep = rng.uniform(0.3, 0.9)
                if run["platform"] == "illumina":
                    # A library molecule that lost an adapter does not
                    # cluster, so it is never seen. What truncation looks
                    # like on Illumina is a short insert - and a short insert
                    # is read *through*, into the opposite adapter.
                    cut = max(20, int(len(amplicon) * keep))
                    amplicon = (
                        amplicon[:cut] if rng.random() < 0.5
                        else amplicon[len(amplicon) - cut :]
                    )
                    mol = build_molecule(amplicon, i7, i5, repeats, assay)
                elif rng.random() < 0.5:
                    mol = mol[: int(len(mol) * keep)]
                else:
                    mol = mol[int(len(mol) * (1 - keep)) :]

            orientation = "+"
            if run["platform"] == "nanopore" and rng.random() < 0.5:
                mol = revcomp(mol)
                orientation = "-"

            rid = f"sim{counter:08d}"
            if run["platform"] == "illumina":
                # Read 1 starts at the inner edge of the forward overhang and
                # read 2 at the inner edge of the reverse one, because that is
                # where the sequencing primers anneal. The index is read
                # separately and is not in either mate - which is precisely
                # why a hopped read is indistinguishable from a real one once
                # the sequencer has written it to a sample's file.
                head = (
                    len(assay["p5_tail"])
                    + len(i5) * repeats
                    + len(assay["pcr1_overhang_fwd"])
                )
                tail = (
                    len(assay["p7_tail"])
                    + len(i7) * repeats
                    + len(assay["pcr1_overhang_rev"])
                )
                r1, s1, d1 = apply_errors(
                    rng, mol[head:][: args.read_length], model
                )
                r2, s2, d2 = apply_errors(
                    rng, revcomp(mol)[tail:][: args.read_length], R2_MODEL(model)
                )
                obs = None
                n_sub, n_indel = s1 + s2, d1 + d2
            else:
                obs, n_sub, n_indel = apply_errors(rng, mol, model)

            records.append(
                ReadRecord(
                    read_id=rid,
                    sample_id=owner,
                    seq_id=seq_id,
                    fate=fate,
                    i7=i7,
                    i5=i5,
                    n_sub=n_sub,
                    n_indel=n_indel,
                    length=len(obs) if obs is not None else len(r1) + len(r2),
                    orientation=orientation,
                )
            )
            if obs is None:
                reads_by_sample[owner].append((rid, r1, r2, seq_id, i5, i7))
            else:
                reads_by_sample[owner].append((rid, obs, seq_id))

    stats = {
        "run_id": run_id,
        "platform": run["platform"],
        "barcode_repeats": repeats,
        "error_model": model_name,
        "total_reads": len(records),
        "by_fate": dict(Counter(r.fate for r in records)),
        "sentinel_combinations": len(sentinels),
        "hopped_into_sentinels": sum(
            1
            for r in records
            if r.fate == "hopped" and (r.i5, r.i7) in set(sentinels)
        ),
    }
    return records, stats, reads_by_sample, composition


# --------------------------------------------------------------------------
# a whole plate at the level of counts: index misassignment
# --------------------------------------------------------------------------


def poisson(rng: random.Random, lam: float) -> int:
    """A Poisson draw: Knuth's product under 30, the normal above."""
    if lam <= 0:
        return 0
    if lam > 30:
        return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))
    limit, k, p = math.exp(-lam), 0, rng.random()
    while p > limit:
        k += 1
        p *= rng.random()
    return k


ROW_LETTERS = "ABCDEFGHIJKLMNOP"


def simulate_plate_counts(
    rng: random.Random, n_rows: int, n_cols: int, *, rate_row: float,
    rate_col: float, unused: int = 0, depth: float = 200000,
    blanks: int = 4, heavy_column: bool = True, min_cluster: int = 5,
    n_taxa: int = 60, grouped_by_column: bool = False,
    contaminate: dict | None = None,
) -> dict:
    """Classified reads per well and taxon on a combinatorial plate, after
    index misassignment at known per-pair rates.

    Reads, not molecules: the synthesis works on classified read counts,
    and at 384 wells of 200,000 reads a read-level simulation would take
    hours to say what this says in a second. What the plate holds is
    built to be hard for a rate estimate, as stool is:

      a host and a "kitome" taxon in every well, bacteria and fungi in
        60-100% of wells, common parasites in 10-40%, a few heavy rare
        ones in 3-8%; so most taxa next door to a well are really in it
      `blanks` blanks (an extraction blank and NTCs), shallow, holding
        only the kitome and what misassignment puts there
      with heavy_column, one column of near-pure single-taxon wells - a
        worm preparation - at three times the depth: the source that
        made a blank look contaminated on the first VESPA plate
      the last `unused` cells of the last column left empty, receiving
        misassigned reads only (the demultiplexer's unused combinations)
      with grouped_by_column, each column's samples share a taxon of
        their own - wells laid out by host or site - which is biology
        lined up with an axis, and must not be read as misassignment
      `contaminate` {well: {taxon: reads}} added to wells' true content -
        real contamination of a blank, say, on top of what leaks in

    Misassignment: each (well, taxon) receives Poisson(rate_row x its
    row-mates' reads + rate_col x its column-mates' reads). A count under
    `min_cluster` is then lost, as a consensus needs that many reads.

    Returns {geom, counts, true, types, matrix, rates}: positions,
    observed and true counts by well, sample types, and the index matrix
    with each cell's reads (every read, not only classified).
    """
    rows = ROW_LETTERS[:n_rows]
    cells = [(r, c) for c in range(1, n_cols + 1) for r in rows]
    unused_cells = set(cells[len(cells) - unused:]) if unused else set()
    used = [rc for rc in cells if rc not in unused_cells]
    names = {rc: f"{rc[0]}{rc[1]}" for rc in used}
    # blanks scattered, not in the heavy column
    heavy_col = (n_cols // 2 + 1) if heavy_column else None
    pool = [rc for rc in used if rc[1] != heavy_col]
    blank_cells = set(rng.sample(pool, min(blanks, len(pool))))
    types = {}
    for rc in used:
        if rc in blank_cells:
            types[names[rc]] = "neg_extraction" if not types or all(
                t != "neg_extraction" for t in types.values()) else "neg_ntc"
        else:
            types[names[rc]] = "unknown"
    taxa = []
    # (name, prevalence, log-share mean, log-share sd)
    taxa.append(("host", 1.0, math.log(0.15), 0.8))
    taxa.append(("kitome", 1.0, math.log(0.001), 0.3))
    for i in range(int(n_taxa * 0.5)):
        taxa.append((f"bact{i}", rng.uniform(0.6, 1.0),
                     math.log(rng.uniform(0.002, 0.05)), 1.0))
    for i in range(int(n_taxa * 0.3)):
        taxa.append((f"para{i}", rng.uniform(0.1, 0.4),
                     math.log(rng.uniform(0.005, 0.05)), 1.2))
    for i in range(max(1, int(n_taxa * 0.1))):
        taxa.append((f"rare{i}", rng.uniform(0.03, 0.08),
                     math.log(rng.uniform(0.1, 0.4)), 0.8))
    true: dict = {}
    for rc in used:
        sid = names[rc]
        if types[sid] != "unknown":
            d = depth * 0.1 * rng.lognormvariate(0, 0.5)
            true[sid] = {"kitome": poisson(rng, d * 0.5)}
            continue
        if heavy_col and rc[1] == heavy_col:
            d = depth * 3 * rng.lognormvariate(0, 0.3)
            true[sid] = {"worm": poisson(rng, d * 0.97),
                         "host": poisson(rng, d * 0.02),
                         "kitome": poisson(rng, d * 0.001)}
            continue
        d = depth * rng.lognormvariate(0, 0.5)
        shares = {}
        for name, prev, mu, sd in taxa:
            if rng.random() < prev:
                shares[name] = math.exp(rng.gauss(mu, sd))
        if grouped_by_column:
            shares[f"group{rc[1]}"] = rng.uniform(0.02, 0.1)
        tot = sum(shares.values())
        true[sid] = {k: poisson(rng, d * v / tot) for k, v in shares.items()}
    for sid, add in (contaminate or {}).items():
        for k, v in add.items():
            true[sid][k] = true[sid].get(k, 0) + v
    row_tot: dict = defaultdict(int)
    col_tot: dict = defaultdict(int)
    for rc in used:
        for k, v in true[names[rc]].items():
            row_tot[(rc[0], k)] += v
            col_tot[(rc[1], k)] += v
    keys = sorted({k for t in true.values() for k in t})
    counts: dict = {}
    for rc in used:
        sid = names[rc]
        own = true[sid]
        got = {}
        for k in keys:
            y = own.get(k, 0)
            lam = (rate_row * (row_tot[(rc[0], k)] - y)
                   + rate_col * (col_tot[(rc[1], k)] - y))
            n = y + poisson(rng, lam)
            if n >= min_cluster:
                got[k] = n
        counts[sid] = got
    i5 = {r: f"i5_{r}" for r in rows}
    i7 = {c: f"i7_{c}" for c in range(1, n_cols + 1)}
    matrix = []
    for rc in cells:
        if rc in unused_cells:
            lam = sum(rate_row * row_tot[(rc[0], k)]
                      + rate_col * col_tot[(rc[1], k)] for k in keys)
            matrix.append({"index_i5": i5[rc[0]], "index_i7": i7[rc[1]],
                           "reads": str(poisson(rng, lam)), "assigned": "no"})
        else:
            matrix.append({"index_i5": i5[rc[0]], "index_i7": i7[rc[1]],
                           "reads": str(sum(counts[names[rc]].values())),
                           "assigned": "yes"})
    return {"geom": {names[rc]: rc for rc in used}, "counts": counts,
            "true": true, "types": types, "matrix": matrix,
            "rates": {"row": rate_row, "col": rate_col}}


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------


def write_illumina(
    rng, reads_by_sample, samples, out_dir, model, read_len
) -> None:
    """Pre-split per-sample R1/R2, the way the sequencer delivers them.

    Reads that carry an index pair no sample claims go to Undetermined
    rather than being dropped. That file is not waste - it is where hopping
    into a sentinel combination becomes visible, and on Illumina it is the
    only place it *can* be seen, since a hopped read that lands on a real
    pair is written into that sample's file indistinguishably.
    """
    os.makedirs(out_dir, exist_ok=True)
    idx = {s.sample_id: i for i, s in enumerate(samples, 1)}
    safe = lambda s: "".join(c if c.isalnum() or c in "._-" else "_" for c in s)
    r2model = R2_MODEL(model)

    for sample_id, reads in reads_by_sample.items():
        if sample_id == "_unassigned":
            name, n = "Undetermined", 0
        elif sample_id.startswith("_"):
            continue
        else:
            name, n = safe(sample_id), idx.get(sample_id, 0)
        base = os.path.join(out_dir, f"{name}_S{n}_L001")
        with gzip.open(f"{base}_R1_001.fastq.gz", "wt") as r1h, gzip.open(
            f"{base}_R2_001.fastq.gz", "wt"
        ) as r2h:
            for rid, r1, r2, _seq_id, i5, i7 in reads:
                tag = f"{i7}+{i5}"
                r1h.write(
                    f"@{rid} 1:N:0:{tag}\n{r1}\n+\n"
                    f"{qual_string(rng, len(r1), model['rate'])}\n"
                )
                r2h.write(
                    f"@{rid} 2:N:0:{tag}\n{r2}\n+\n"
                    f"{qual_string(rng, len(r2), r2model['rate'])}\n"
                )


def write_nanopore(rng, reads_by_sample, out_path, model,
                   model_name="nanopore_hac", run_id="RUN01") -> None:
    """One pooled FASTQ: nanopore demultiplexing happens in the pipeline.

    Headers carry dorado's SAM-style tags, because the pipeline reads the
    basecalling model out of RG and records it as provenance. Two details
    that were previously wrong here and are worth keeping right:

    the tags are TAB separated after the read id, not space separated,
    which is what a parser written against real data expects;

    and the model named is the model actually simulated. It used to say
    sup whatever error rate had been applied, so a fixture generated as
    fast would report itself as sup - exactly the kind of false
    provenance the recording exists to prevent.
    """
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    variant = model_name.replace("nanopore_", "") or "hac"
    rg = f"{uuid_like(rng)}_dna_r10.4.1_e8.2_400bps_{variant}@v5.2.0"
    allr = [x for reads in reads_by_sample.values() for x in reads]
    rng.shuffle(allr)
    with gzip.open(out_path, "wt") as fh:
        for rid, mol, _ in allr:
            qual = qual_string(rng, len(mol), model["rate"])
            fh.write(
                f"@{rid}\tRG:Z:{rg}\tPU:Z:SIM00000\tLB:Z:{run_id}\n"
                f"{mol}\n+\n{qual}\n"
            )


def uuid_like(rng) -> str:
    """A 36-character identifier shaped like the one dorado writes, so the
    RG parser is exercised on the real layout rather than a shortcut."""
    h = "".join(rng.choice("0123456789abcdef") for _ in range(32))
    return f"{h[:8]}-{h[8:12]}-{h[12:16]}-{h[16:20]}-{h[20:]}"


def write_composition(composition: dict, run_id: str, path: str) -> None:
    """True input copies, before amplification and before sequencing.

    This is the only place the answer exists. Everything downstream is an
    estimate of it, and an estimate that can only be compared with other
    estimates is not being tested.
    """
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run_id", "sample_id", "seq_id", "input_copies"])
        for sid, counts in sorted(composition.items()):
            for seq_id, n in sorted(counts.items()):
                w.writerow([run_id, sid, seq_id, f"{n:.4f}"])


def write_truth(records: list[ReadRecord], path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "read_id", "sample_id", "seq_id", "fate",
                "index_i7", "index_i5", "n_sub", "n_indel",
                "length", "orientation",
            ]
        )
        for r in records:
            w.writerow(
                [r.read_id, r.sample_id, r.seq_id, r.fate, r.i7, r.i5,
                 r.n_sub, r.n_indel, r.length, r.orientation]
            )


# --------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Simulate reads with a truth table.")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--assay")
    ap.add_argument("--outdir", default="simulated")
    ap.add_argument("--reads-per-sample", type=int, default=2000)
    ap.add_argument("--read-length", type=int, default=250)
    ap.add_argument(
        "--error-model", choices=sorted(ERROR_MODELS), default="nanopore_hac"
    )
    ap.add_argument("--prevalence", type=float, default=0.35,
                    help="probability a field sample carries any given target")
    ap.add_argument("--hop-rate", type=float, default=0.005,
                    help="both indices rewritten at random (or onto an "
                    "unused combination)")
    ap.add_argument("--misread-i5-rate", type=float, default=0.0,
                    help="share of reads whose i5 is read as another of "
                    "the plate's i5s: the read lands in the same column")
    ap.add_argument("--misread-i7-rate", type=float, default=0.0,
                    help="share of reads whose i7 is read as another of "
                    "the plate's i7s: the read lands in the same row")
    ap.add_argument("--chimera-rate", type=float, default=0.01)
    ap.add_argument("--carryover-rate", type=float, default=0.002)
    ap.add_argument("--dimer-rate", type=float, default=0.02)
    ap.add_argument("--offtarget-rate", type=float, default=0.005)
    ap.add_argument("--truncation-rate", type=float, default=0.03)
    ap.add_argument(
        "--mocks",
        help="mock definitions; wells listed here get the declared "
             "composition instead of a simulated one",
    )
    ap.add_argument(
        "--competition", type=float, default=0.0,
        help=(
            "strength of composition-dependent PCR bias, 0-1. 0 is the "
            "simple multiplicative model; above 0 the observed ratio "
            "compresses toward the middle as one template dominates"
        ),
    )
    ap.add_argument("--no-spike", action="store_true",
                    help="omit the quant standards, to exercise the unspiked path")
    ap.add_argument("--seed", type=int, default=20260918)
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    refs, samples, runs, assay = load_inputs(args)
    os.makedirs(args.outdir, exist_ok=True)

    # Declared mock compositions, if any. These replace the simulated
    # composition for the wells they name, so a mock series can be put
    # through the whole pipeline and the calibrator checked against the
    # amounts that were actually specified.
    args.mock_compositions = {}
    if getattr(args, "mocks", None):
        for r in read_table(args.mocks):
            sid, seq = r.get("sample_id"), r.get("seq_id")
            if not sid or not seq:
                continue
            try:
                amt = float(r.get("amount", ""))
            except ValueError:
                continue
            if amt > 0:
                args.mock_compositions.setdefault(sid, {})[seq] = amt

    all_records: list[ReadRecord] = []
    all_stats = []

    for run_id, run in runs.items():
        grp = [s for s in samples if s.run_id == run_id]
        if not grp:
            continue
        # Illumina always uses its own error profile; --error-model selects
        # among the nanopore basecalling models, which is the choice that
        # actually varies.
        model_name = (
            "illumina" if run["platform"] == "illumina" else args.error_model
        )
        model = ERROR_MODELS[model_name]
        records, stats, reads_by_sample, composition = simulate_run(
            rng, run_id, run, grp, refs, assay, args, model_name
        )
        out_dir = os.path.join(args.outdir, run_id)
        if run["platform"] == "illumina":
            write_illumina(rng, reads_by_sample, grp, out_dir, model, args.read_length)
        else:
            write_nanopore(
                rng, reads_by_sample,
                os.path.join(out_dir, f"{run_id}.fastq.gz"), model,
                model_name=args.error_model, run_id=run_id,
            )
        write_truth(records, os.path.join(args.outdir, f"{run_id}.truth.tsv"))
        write_composition(
            composition, run_id,
            os.path.join(args.outdir, f"{run_id}.composition.tsv"),
        )
        all_records += records
        all_stats.append(stats)

        print(f"{run_id} ({run['platform']}): {stats['total_reads']:,} reads")
        for fate, n in sorted(stats["by_fate"].items(), key=lambda kv: -kv[1]):
            print(f"    {fate:<12} {n:>8,}  {100*n/stats['total_reads']:5.2f}%")
        print(
            f"    sentinels: {stats['sentinel_combinations']} combinations, "
            f"{stats['hopped_into_sentinels']} hopped reads landed in them"
        )

    with open(os.path.join(args.outdir, "simulation_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(
            {
                "seed": args.seed,
                "rates": {
                    "hop": args.hop_rate,
                    "misread_i5": args.misread_i5_rate,
                    "misread_i7": args.misread_i7_rate,
                    "chimera": args.chimera_rate,
                    "carryover": args.carryover_rate,
                    "dimer": args.dimer_rate,
                    "offtarget": args.offtarget_rate,
                    "truncation": args.truncation_rate,
                },
                "runs": all_stats,
            },
            fh,
            indent=2,
        )
    print(f"\nTruth tables and summary in {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
