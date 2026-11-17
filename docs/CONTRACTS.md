# Database mode: layers and contracts

Database mode (`--reference_database`, [REFERENCE_DATABASE.md](REFERENCE_DATABASE.md))
runs in three layers. Each hands the next a fixed contract, so a new
sequencing platform needs only an adapter, a new database only a protocol,
and neither touches synthesis or the report.

```
 PLATFORM (layer 1)           PROTOCOL (layer 2)                SHARED (layer 3)
 get each sample's reads      primer trimming -> calls          synthesis, report

 illumina ──────┐                                           ┌── synthesise_db.py
 nanopore ──────┤                                           │   (blank limits, floors,
 pacbio ────────┼─ CONTRACT A ─► vespa_v4 | nemits2 ─ CONTRACT B ─┤    misassignment, plate QC)
 predemultiplexed┘  reads,         (bin/protocols/,       clusters.tsv,  └── report_db.py
                    read type,      chosen by the         taxa.tsv,
                    accounting,     release)              assign_db.json
                    plate
```

Panel mode (`--reference`) is not layered this way and is unchanged.

The contracts are checked by `tests/test_contracts.py`, which runs every
protocol on its fixture, over every read type it supports, and checks
what it writes against this document and against the column lists in
`bin/ampseqlib.py` (`CLUSTER_COLUMNS`, `TAXA_COLUMNS`,
`ASSIGN_SUMMARY_KEYS`, `ASSIGN_READ_FATES`, `DB_CALL_COLUMNS`,
`DB_CALL_OPTIONAL`).

---

## Layer 1: platform adapters

| platform | reads arrive | read type | adapter |
|---|---|---|---|
| `illumina` | split per sample by the sequencer | `paired_short` | pre-flight resolves each sample's R1/R2 by `file_prefix` |
| `nanopore` | as a flow cell's FASTQs | `single_long_noisy` | `split_concatemers.py`, then `demux_ont.py` on the construct the assay JSON describes |
| `predemultiplexed` | already split per sample, by lima, Dorado or a sequencing centre | `single_long_noisy` (default) or `single_long_accurate` | pre-flight resolves each sample's FASTQ(s) by `file_prefix`, like Illumina but one read per molecule |
| `pacbio` | as a library's CCS/HiFi FASTQ or BAM | `single_long_accurate` | `demux_pacbio.py` on the samplesheet's tags (`index_i5` forward, `index_i7` reverse); reads lima already split are `predemultiplexed` |

`predemultiplexed` and `pacbio` are database mode only. Its samples are found under the
run's `path` as `<prefix>.fastq.gz` (or `.fq.gz`, `.fastq`, `.fq`),
as chunks `<prefix>.*`, `<prefix>_*` or `<prefix>-*` with a FASTQ suffix,
or as every FASTQ in a folder `<prefix>/` (Dorado and MinKNOW's
`barcodeNN/`); `<prefix>` is `file_prefix`, else `sample_id`. lima names
its output `<movie>.<bc>--<bc>`, so a lima run sets `file_prefix` to that.
Indices may be left blank on a predemultiplexed run; a well without them
has no plate position (below).

A `pacbio` run is one library: its `path` is the reads (FASTQ or
unaligned BAM, or a directory of a library's SMRT cells), and only that
run's tags are matched, since a shared SMRT cell carries other
libraries'. The adapter writes per-sample FASTQs (tags off, primers on,
forward primer first), `demux_summary.json` in nanopore's tiers A-D (no
E: CCS has removed the adapters, so a read is one molecule) and an
`index_matrix.tsv` of forward x reverse tag pairs. USAGE.md §2.11.

### Contract A

What every adapter hands the protocol layer, for each run:

1. **Per-sample reads.** One FASTQ per nanopore or PacBio sample
   (demultiplexed here); R1 and R2 per Illumina sample; one or more
   FASTQs per predemultiplexed sample, read as one. Pre-flight's
   `file_resolution.tsv` names every file a sample was resolved to.
2. **The read type**, one of:

   | read type | meaning |
   |---|---|
   | `paired_short` | two short, accurate mates per molecule (Illumina) |
   | `single_long_noisy` | one read per molecule, ~1-5% error (nanopore) |
   | `single_long_accurate` | one read per molecule, <0.1% error (PacBio HiFi) |

   Pre-flight writes it per run to `preflight/read_types.tsv`
   (`run_id`, `platform`, `read_type`); runs.tsv's optional `read_type`
   column sets it for a predemultiplexed run.
3. **The demux accounting**, `demux/demux_accounting.json`, one per run,
   written by `bin/demux_accounting.py` whatever the platform:

   | key | type | meaning |
   |---|---|---|
   | `contract`, `version` | str, int | `"amphid demux accounting"`, `1` |
   | `run_id`, `platform`, `read_type` | str | as above |
   | `demultiplexed_by` | str | the sequencer, amphid, or upstream of amphid |
   | `reads_off_instrument` | int or null | reads as they came off the instrument (nanopore: before splitting at internal adapters) |
   | `reads_in` | int or null | reads (units, after splitting) entering demultiplexing |
   | `assigned` | int or null | reads written to a sample |
   | `fates` | {fate: int} | `assigned`, `one_end_only`, `chimeric`, `unmatched`, and `internal_adapter` where the adapter looks for one (nanopore); `{}` when the adapter cannot know (Illumina, predemultiplexed) |
   | `balanced` | bool or null | whether the fates add up to `reads_in` |
   | `per_sample` | {sample_id: int or null} | reads delivered to each declared sample; null where the adapter did not count them (assignment's `reads_in` then holds the count) |
   | `plate` | {sample_id: {row, col}} | plate position of each well that has one |
   | `native` | {demux, split} | the platform's own summaries, verbatim, or null |

   Null is "not known at this layer", never zero; a fate the adapter
   does not have is left out, not zero. Synthesis reads the accounting,
   not a platform's own files: `native` carries the detail a shared fate
   cannot (nanopore's end reasons, combination grid, basecaller;
   PacBio's reasons, per-end matches, tag edits), and
   the read funnel and plate QC are built from it.
4. **Plate positions**, where the platform has them: from the index map,
   the i5 giving the row and the i7 the column (`file_resolution.tsv`
   `plate_row`, `plate_col`, and the accounting's `plate`). A PacBio
   sample has one where its forward tag is an i5 in the map and its
   reverse tag an i7. The
   misassignment model and the plate-position tests in synthesis need
   them; a well without one is left out of both rather than guessed.

---

## Layer 2: protocols

A protocol is everything from primer trimming to a call: trimming,
clustering or denoising, classification, and the files of contract B.
Each lives in `bin/protocols/<name>.py`, defines one `Protocol` subclass
(`bin/protocols/__init__.py`) and exposes it as `PROTOCOL`. A protocol
may use anything in `ampseqlib` (minimap2, primer matching, FASTQ input,
the LCA and identity floors) and is called by the read-type routes; it
never imports another protocol.

| protocol | release | classification | calling rules |
|---|---|---|---|
| `vespa_v4` | vespa-18s-refdb, 18S V4 | LCA over minimap2 hits within `--db_lca_margin`, cut back by identity floors per rank | `floors` |
| `nemits2` | nemabiome-its2-refdb, ITS2 | each candidate record scored over its own span, query overhang free; LCA of the ties | `tiered` (default), `gap`, `floors` |

### Which protocol a release gets

1. `--db_protocol`, when it is not `auto`;
2. else the release's `manifest.json` `"protocol"`;
3. else inferred from what manifests already say: a release scored over
   its records' span (`"score_over": "reference span"`), or one called by
   the gap or tiered rule (manifest `calling_rule`, or `--db_calling_rule
   gap` / `tiered`), is
   `nemits2`; everything else `vespa_v4`. So every release built before
   the field existed runs exactly as before.

Pre-flight resolves it the same way, says which in its log, and fails
when the protocol does not exist, cannot take a run's read type, does
not call by the calling rule asked for, or lacks a release file it reads
(`Protocol.required_files`, beyond the release format's own).

### Read types

| | `paired_short` | `single_long_noisy` | `single_long_accurate` |
|---|---|---|---|
| `vespa_v4` | two fixed windows, denoised, pair floors | reference-guided consensus, nanopore floors | the consensus route plus three steps (below), same floors |
| `nemits2` | overlapping mates merged, else two windows (R2 with its 28S stretch), span-scored | consensus carrying the 28S stretch, span-scored | the consensus route plus three steps (below), span-scored, gap rule |

`single_long_accurate` takes the consensus route with three steps of its
own, none of which runs for `single_long_noisy`: homopolymer-length
copies of a larger consensus are absorbed (`merge_homopolymer`), a
smaller consensus takes a larger one's run lengths and keeps its
substitutions (`correct_homopolymer`), and a consensus no read supports
exactly is dropped (`unsupported`). The first two switch off apart
(`assign_db.py --homopolymer-merge-edits 0`, `--homopolymer-correct-reach
0`); the correction keeps the read type's cap and ratio when the merge is
off. `assign_db.json` counts `homopolymer_merges` (every consensus
absorbed, by either step) and `homopolymer_corrected` (every consensus
corrected, including those that then joined another). All three were added after S-BSST527
(below); on the simulated HiFi benchmark they change no call - every row
of the floors grid is the same, and spurious clusters fall from 3 of 604
to 1 of 602. It was enabled on what was measured:

- **vespa_v4: the floors, on simulated HiFi reads.** The leave-out
  benchmark the floors come from (`tests/bench_db_floors.py`, vespa
  v0.9.0, 100 genera and 150 species held out, 150 samples) was run with
  nanopore (`nanopore_hac`) and HiFi (`pacbio_hifi`: 0.1% error, mostly
  homopolymer indels) reads. At the default floors:

  | reads | in database: deepest rank, wrong | species held out: wrong | genus held out: wrong | consensus identical to its source |
  |---|---|---|---|---|
  | nanopore_hac | 96.1%, 0% | 16.6% | 35.2% | 90.9% |
  | pacbio_hifi | 98.0%, 0% | 18.8% | 45.0% | 96.0% |

  The two runs draw different samples (one random stream serves both);
  over the source sequences both drew, their wrong-call rates are the
  same to the read - 19.6% (46 held-out species), 43.1% (51 held-out
  genera) - so the difference above is the draw, not the reads. A HiFi
  consensus is more often exact, and the floors treat it as they treat
  a nanopore one. (A held-out genus taking a confamilial's name is the
  floors' known limit on any platform, REFERENCE_DATABASE.md §11.3.)
- **nemits2: the gap and tiered rules, on S-BSST527** (Halvarsson & Höglund 2021:
  158 sheep larval cultures, NC1/NC2, PacBio RSII CCS in four pools,
  914,810 reads; nemabiome-its2-refdb v0.2.0). Compared with
  nemabiome-its2-refdb's own analysis of the same reads (its
  demultiplexer, DADA2 with PacBio settings, `classify_asvs.py`, gap
  rule):
  - 91.6% of clustered reads are in a consensus identical to one of
    that sample's ASVs (over the record's span), and 99.3% of those
    reads get the ASV's call (the 18 sequences that differ are ties a
    gap of 2-3 edits separates: species there, genus here under
    `--db_gap_edits 4`).
  - The other 8.4% are small groups whose consensus carries 1-5
    homopolymer indels: reads are grouped by their best database hit,
    and RSII CCS reads sharing an indel group together. Their calls
    are safe - 49.5% of those reads take the ASV's call, 50.1% stop at
    a higher rank consistent with it, 5 reads (one consensus) name
    another species - but they cost resolution: 94.6% of clustered
    reads are named to species against DADA2's 99.3%, almost all of the
    difference *Haemonchus contortus* left at genus. Filtering reads at
    two expected errors (DADA2's filter) changed this little (13.7% to
    12.0% of one sample's reads at genus).
  - Per sample, the composition at the rank reached differs from
    DADA2's by a median Bray-Curtis of 0.024 (154 samples with 200+
    reads; 90th percentile 0.117), and study-wide, filtered as the
    paper was (taxa >= 1% of a sample), every species' share is within
    0.4 points of DADA2's and its prevalence within one sample, but for
    *H. contortus* (60.8% against 64.5%; the rest is *Haemonchus* at
    genus, 4.1%) and two cyathostomins the 4-edit gap holds at genus
    (*Cyathostomum pateratum*, *Coronocyclus coronatus*: gap 2-3).

  - Those figures are before `single_long_accurate` absorbed
    homopolymer-length variants (`merge_homopolymer` in
    bin/assign_db.py: a consensus whose bases equal a larger one's once
    runs are collapsed, run lengths 4 or fewer bases apart in total, and
    the larger standing for 4 or more times its reads). On the same
    reads it absorbed 152 consensus sequences ("merge" below). Most of
    what then still stopped at genus was a consensus differing from the
    sample's main *H. contortus* by one or two substitutions (a real
    haplotype) and also by one or two homopolymer indels. So a smaller
    consensus is then given the run lengths of the nearest larger one
    (4x its reads, within 8 edits) when every indel between them is a
    run length - the base in a run of two or more in the sequence that
    has it - keeping its substitutions (`correct_homopolymer`; "+
    correct" below, 294 consensus sequences corrected). Where any indel
    is not a run length nothing is corrected: allowing those moved as
    many reads' consensus away from the ASVs as towards them, and lost
    304 sample-ASV pairs. Nor is a consensus that is its group's
    database record: a run length a record carries is not a read error
    (without this, one consensus identical to an ASV, TTT where a 6x
    larger neighbour had TT, was corrected away).

    Most of what was left above species was then not a sequence in the
    sample at all: a consensus of noisy reads that one database record
    fits best because their errors fall where its alleles differ (61% of
    those reads are in a consensus that is exactly two of the sample's
    ASVs spliced together). A real sequence is read exactly by about 40%
    of its reads in this run; these by 3%. DADA2 never builds them: it
    filters such reads (maxEE 2) or explains them as errors of a more
    abundant ASV. So, for `single_long_accurate`, a consensus whose
    exact reads are improbably few (binomial, p < 0.001) at a third of
    the sample's largest consensus' rate is dropped and its reads counted
    again against the rest (`unsupported`; "+ support" below: 115
    consensus sequences, 10,077 reads). The two sample-ASV pairs this
    loses are sequences DADA2 did not find in those samples at all - a
    consensus equal by chance to another sample's ASV, read exactly by
    1-3% of its reads.

    | S-BSST527 | before | merge | + correct | + support | DADA2 |
    |---|---|---|---|---|---|
    | clustered reads named to species | 94.55% | 94.90% | 96.36% | 97.07% | 99.28% |
    | *Haemonchus* at genus, study-wide (paper filter) | 4.09% | 3.75% | 2.26% | 1.57% | 0.03% |
    | *H. contortus*, study-wide | 60.81% | 61.15% | 62.63% | 63.26% | 64.49% |
    | reads naming another species than the nearest ASV | 5 | 5 | 5 | 5 | |
    | clustered reads in a consensus identical to an ASV | 91.59% | 92.05% | 92.02% | 93.01% | |
    | sample-ASV pairs with an identical consensus | 1604 | 1604 | 1611 | 1609 | |
    | per-sample Bray-Curtis, median / 90th | 0.024 / 0.117 | 0.020 / 0.104 | 0.014 / 0.078 | 0.010 / 0.064 | |
    | unclustered reads | 4,233 | 4,238 | 4,319 | 5,369 | |

    The gap rule then took the release classifier's exact-match clause
    (`GAP_EXACT`, `--db_gap_exact 2`: an exact match to its record is a
    species with the next species 2 or more edits off; before, 0.65% of
    reads were in a consensus identical to an ASV that
    `classify_asvs.py` named and the gap rule here left at genus). It is
    the protocol's rule, so it applies to every read type; it changes no
    call on the fixtures. On S-BSST527: species 97.07% -> 97.78%,
    identical sequences called differently from DADA2's 18 (2,884 reads)
    -> 2 (26 reads), per-sample Bray-Curtis median 0.010 -> 0.0098 and
    worst sample 0.317 -> 0.100, wrong-species reads 5.

    No nanopore NC1/NC2 run was at hand (nor, as far as a search of ENA
    and the literature found, is one public), so the clause was measured
    on nanopore consensus sequences by a leave-out simulation instead
    (`tests/bench_db_gap.py`: v0.2.0 with 100 species and 30 genera held
    out, or, for the strongylid families, 60 and 8; 300 samples each;
    `--gap-exact` 0 and 2 on the same reads). Every call it changed was a
    species in the database named to species instead of genus; no
    held-out species or genus was called differently, so wrong calls did
    not move. A held-out species named a sister whose record it contains
    exactly (identical wherever both have sequence: 54 species pairs in
    v0.2.0, mostly plant parasites and ascaridoids, none in livestock
    strongylids) is counted apart, as indistinguishable - with both in
    the database the gap rule stops such a pair at genus:

    | sources | in database: species, 0 -> 2 | species held out: wrong (+ indistinguishable) | genus held out: wrong |
    |---|---|---|---|
    | whole release, nanopore_hac | 82.0% -> 90.8% | 3.7% (+1.3%), unchanged | 0.3% -> 0.3% |
    | strongylids, nanopore_hac | 93.3% -> 97.2% | 5.7% -> 5.7% | 2.3% -> 2.3% |
    | strongylids, nanopore_sup | 93.5% -> 97.3% | 5.7% -> 5.7% | 2.3% -> 2.3% |

    From nemits2 v0.2.2 the release flags those pairs itself
    (`contained_in` / `contains` in `taxonomy.tsv`: 151 record pairs of
    57 species pairs), and the gap rule takes a flagged record as 0
    edits from the other species, as the release's margin does: a read
    matching the shorter record and differing from the longer only where
    the shorter has no sequence is no evidence for the shorter's
    species. Only records in the database count (the benchmark's
    leave-out copy now holds only the records left in, so a held-out
    record cannot be seen through a flag). The whole-release run on
    v0.2.2, with and without the flags on the same reads, gives the
    table above unchanged: of 2,640 calls, 26 differ, 24 of them already
    at genus (the gap 1 -> 0) and 2, both from one 5-read consensus of
    *Aphelenchoides besseyi* reads named *A. oryzae*, from a wrong
    species to genus. The fixtures carry no flags, so no fixture output
    changes. On S-BSST527 (four pools, `-profile field`) v0.2.2 with the
    flags gives every one of the 2,621 consensus sequences the call
    v0.2.0 gave it; the synthesis tables differ only in the release's
    version: no consensus there ties with a flagged record.

    **Indistinguishable, widened (2026-11-17).** Exact containment
    missed pairs the release cannot separate either. *H. sinense* and
    *H. liparis* share an ITS2 to an edit (31 of v0.2.4's 54 wrong
    calls). The benchmark now counts a held-out species' call as
    indistinguishable when its source record is fewer than `GAP_EDITS`
    edits from a typical record of the species named. A typical record
    is one within the gap of another record of its own species, or that
    species' only record. A record far from the rest of its species is
    taken as misidentified, so a call through it stays wrong. One case
    still slips through: a group of misidentified records, being close
    to each other, passes as typical. That is v0.2.4's short croaker
    lineage under "*H. fabri*", which accounts for 2 calls.

    Rescored from the same calls (whole release, `nanopore_hac`, 300
    samples, `--gap-exact` 2; v0.2.5 holds out a different draw):

    | species held out | wrong | indistinguishable |
    |---|---|---|
    | v0.2.4, exact containment | 18.0% | 0.0% |
    | v0.2.4, within the gap | 8.0% | 10.0% |
    | v0.2.5, exact containment | 19.0% | 0.7% |
    | v0.2.5, within the gap | 5.7% | 14.0% |

    On v0.2.4 the wrong calls are mostly through two records v0.2.5
    drops: MW370746 (an *H. aduncum* ITS2 labelled *H. fabri*, 8 calls)
    and MW413285 (*H. zhoushanense*, 3 edits from *H. amoyense*, 7
    calls). On v0.2.5 they are:
    - 7 *H. sinense* sources beyond the gap of every *H. liparis*
      record;
    - 7 through PV138334, an *H. urticae* record 5 edits from the
      species' other record and 1 from *H. cruciferae*;
    - 3 single calls (*P. vulnus* -> *P. penetrans*, one *H. elachista*
      -> *H. oryzae*, *Heterorhabditis taysearae* -> Steinernematidae).

    The genus held out has one wrong call (0.3%): *Orloffia bisonis*
    named *Teladorsagia*.

    `single_long_noisy` is unchanged, byte for byte, on every fixture.

    **Under the tiered rule (2026-11-16).** The figures above are the gap
    rule's. Since the tiered rule became the default for a release
    scored over its records' span, the four pools were run again
    (`-profile field`, v0.2.2) with every change since: the tiered rule,
    repeat-aware scoring, informal records after the named ones, the
    release's lineage. The ASVs were classified again with the release's
    current `classify_asvs.py` on the same release. The gap rule's call,
    which each tiered call carries beside it, reproduces the figures
    above (97.79% to species, Bray-Curtis median 0.0098, the same 2
    identical sequences called differently), so the other changes moved
    nothing here and what differs is the rule:

    | S-BSST527, v0.2.2 | gap rule | tiered |
    |---|---|---|
    | clustered reads named to species (DADA2) | 97.79% (99.28%) | 98.44% (99.32%) |
    | *Haemonchus* at genus, study-wide, paper filter (DADA2) | 1.57% (0.03%) | 0.91% (0.00%) |
    | *H. contortus*, study-wide (DADA2) | 63.25% (64.49%) | 63.90% (64.53%) |
    | identical sequences called differently from DADA2's | 2 (26 reads) | 1 (22 reads) |
    | per-sample Bray-Curtis, median / 90th / max | 0.0098 / 0.052 / 0.100 | 0.0047 / 0.038 / 0.075 |
    | taxa at >= 1% of a sample in amphid only | 63 | 42 |
    | reads naming another species than the nearest ASV | 5 | 5 |

    The tentative tier names 2,848 reads (0.65% of clustered) to
    species, mostly *H. contortus* consensus sequences 98-99% to their
    records; the one identical sequence still called differently is a
    *Cylicostephanus goldi* consensus (22 reads) that `classify_asvs.py`
    leaves at genus. Every other species' study-wide share is within 0.4
    points of DADA2's under either rule.

    **Grouping tied reads (2026-11-16).** Run on v0.2.4, which drops 9
    records none of them sheep parasites, the same reads gave slightly
    worse figures (98.31% to species, *Haemonchus* at genus 1.01%), and
    475 consensus sequences (15,347 reads, 3.5%) were built differently.
    About 8% of these reads score exactly as well on several records,
    and minimap2 breaks such a tie by the index's order: dropping
    unrelated records moved 133 of one sample's 3,037 reads to another
    group. For `single_long_accurate` a read is now grouped under the
    record, among those tied for its best score, that most of the
    sample's reads tie on (`group_ties`; minimap2 reports up to 50 hits
    per read to see them; `assign_db.py --group-ties 0` turns it off,
    and `assign_db.json` records `group_ties`). The ASVs keep the
    classification above; amphid runs `assign_db.py` alone on the
    pipeline's demultiplexed reads, which reproduces the pipeline's
    output:

    | S-BSST527, tiered | primary hit, v0.2.2 | v0.2.4 | most-tied record, v0.2.2 | v0.2.4 |
    |---|---|---|---|---|
    | clustered reads named to species | 98.44% | 98.31% | 98.90% | 98.90% |
    | *Haemonchus* at genus, study-wide, paper filter | 0.91% | 1.01% | 0.50% | 0.50% |
    | *H. contortus*, study-wide | 63.90% | 63.80% | 64.24% | 64.25% |
    | per-sample Bray-Curtis, median / 90th | 0.0050 / 0.038 | 0.0063 / 0.038 | 0.0044 / 0.025 | 0.0043 / 0.029 |
    | reads naming another species than the nearest ASV | 5 | 5 | 5 | 5 |
    | reads on a v0.2.2 consensus v0.2.4 does not build | | 3.50% | | 0.45% |

    It also gathers reads that share a CCS homopolymer error - which tie
    on records of their own - back with their sequence's other reads.
    It cost a *Cyathostomum* haplotype, 13% of HRG14-084 and kept by
    DADA2 at genus: tied with a commoner record, it joined that record's
    group and was counted as *C. catinatum* (that sample's Bray-Curtis
    0.017 -> 0.116). Grouping by the whole set of tied records kept it
    but rebuilt the error groups (98.16% to species, *Haemonchus* at
    genus 1.13%, 2.4x the time), and giving only sets of 5+ reads their
    own group did no better (97.97%, 1.30%, 3.1x). Most-tied costs 1.2x
    the assignment time.

    **Splitting haplotypes above a third of a group (2026-11-17).** The
    haplotype was in its group at 48%, two adjacent substitutions, 143
    of 144 carriers shared - and the split's linkage test never linked
    it: co-occurrence 3x what independence gives cannot happen for
    alleles above a third of the reads. Two substitutions are now also
    linked when independence would give that many co-carriers with a
    hypergeometric chance below `hap_p` (1e-6). Deletions are left to
    the ratio: low-quality CCS reads share error-prone positions, and
    linking their deletions split *H. contortus* groups into consensus
    sequences six deletions off the sample's, matching no ASV (17 such,
    1,333 reads, before deletions were excluded). On S-BSST527:

    | S-BSST527, tiered, v0.2.4 | most-tied | + split |
    |---|---|---|
    | haplotype splits | 100 | 109 |
    | clustered reads named to species | 98.90% | 98.84% |
    | *Haemonchus* at genus, study-wide, paper filter | 0.50% | 0.50% |
    | per-sample Bray-Curtis, median / 90th / max | 0.0043 / 0.029 / 0.116 | 0.0043 / 0.026 / 0.072 |
    | consensus gained identical to an ASV of the sample | | 6 (802 reads) |
    | consensus gained 2+ edits from every ASV of the sample | | 4 (177 reads) |
    | reads naming another species than the nearest ASV | 5 | 5 |

    The six gained are haplotypes DADA2 finds too (HRG14-084's among
    them; *Teladorsagia*, *Cylicocyclus* and *Trichostrongylus* in
    others), each differing from its group's other consensus by
    substitutions alone; none was lost. Species calls fall 0.06 points
    because some of them are named only to genus, as DADA2 names them.
    Release-to-release churn stays 0.50% (v0.2.2 against v0.2.4). The
    nanopore leave-out simulation (`tests/bench_db_gap.py`, v0.2.4, 150
    samples) is byte-identical, though its samples carry no hidden
    haplotypes to find.
    `single_long_noisy` is unchanged, byte for byte.

  HiFi from Sequel II or Revio is more accurate than RSII CCS, so the
  resolution cost should be smaller there; it has not been measured on
  such a run.

The routes that call a protocol are `bin/assign_db.py` (single reads:
trim, map, build and split consensus sequences, then
`Protocol.classify_consensus`) and `bin/assign_db_pairs.py` (pairs:
orient, window - `Protocol.read_windows`, `Protocol.merge_mates` -
denoise, then `Protocol.classify_pairs`). Both take `--protocol`.

---

## Contract B: what a protocol writes

Per sample, in `assign_db/`:

### `<sample>.clusters.tsv`

One row per sequence (consensus, or denoised pair), most reads first.
These columns, in this order, from every protocol:

| column | meaning |
|---|---|
| `sample_id` | the sanitised sample id |
| `cluster` | `<sample>_c<n>` (consensus) or `<sample>_p<n>` (pair), n from 1 by reads |
| `reads` | reads (pairs) counted to it |
| `rel_abundance` | reads over the sample's clustered reads, 6 decimals; blank when none |
| `rank` | the rank the call reached; blank when none (or a chimera) |
| `taxon` | the name at that rank |
| `label` | what the report says: the taxon, or "unclassified, NN% to X", "X, unresolved below R", "chimera of A and B", "chimera across the gap: R1 NN% to X, R2 NN% to Y" (nemits2 window pairs) |
| `identity` | identity to the best record, 4 decimals |
| `coverage` | share of the query the best record's alignment covers, 3 decimals |
| `lineage` | `;`-joined names down to the rank reached |
| `best_hit` | the best record's id |
| `best_hit_lineage` | its full lineage |
| `near_hits` | records tied with it |
| `lca_rank` | the rank the tied records agree to, before floors |
| `db_resolution` | the deepest rank the release resolves the best record to on this platform |
| `db_margin` | a species call's edits to the nearest other species, "large", or blank |
| `origin` | `group`, `split…` (haplotype), `pair`, `pair_split`, `whole_insert` or `chimera` |
| `merged` | consensus sequences merged into this one |
| `chimera_of` | `A+B@breakpoint` for a chimera, with "within genus, kept" for one that is kept; `R1:<record>+R2:<record>@gap` for a join in the gap found window by window (nemits2) |
| `length` | sequence length |
| `database_version` | the release's version |
| `consensus` | the sequence (a pair's gap as ten Ns) |

Then the protocol's declared optional columns, `Protocol.extra_columns`,
in its order. Only `nemits2` declares any:

| column | meaning |
|---|---|
| `span` | record bases the best record was scored over |
| `edits` | edits to it over that span |
| `gap` | edits further the nearest other species is (blank unless the ties name one species) |
| `nearest` | the tied species and their identity ("Haemonchus contortus 99.6%") |
| `next` | the nearest other species and how many edits further ("Haemonchus placei +2"), or "none close" |
| `grade` | the best record's grade (complete, partial…) |
| `tentative` | `yes` when the species was named only by the tiered rule's lower floor (its label is "<species> (tentative)") |
| `gap_rule_rank` | the gap rule's own call, beside a tiered one: its rank… |
| `gap_rule_taxon` | …and taxon (the same as `rank` / `taxon` unless the call is tentative) |
| `release_lineage` | the release's `lineage` that every tied best record carries - a line within a species ("goat line") or a genus-level record's informal name ("Necator sp. type II") - else blank; a species or genus call says it in its label too |

### `<sample>.taxa.tsv`

Exactly: `sample_id`, `rank`, `taxon`, `label`, `reads`, `rel_abundance`,
`clusters`, `best_identity`, `lineage`, `database_version`. One row per
(rank, taxon, label), most reads first; reads add up to the clusters'.

### `<sample>.consensus.fasta`

One record per clusters.tsv row, in its order:
`>cluster reads=N rank=R taxon=T identity=I`.

### `<sample>.assign_db.json`

At least: `sample_id`, `mode` (`"database"`), `protocol` (the
protocol that made the calls: `vespa_v4`, `nemits2`), `database` {version,
name, path}, `minimap2`, `elapsed_s`, `reads`, `reads_note`, `clusters`
(rows in clusters.tsv), `haplotype_splits`, `chimeras` {sequences, reads,
within_genus_kept, and `across_gap` {sequences, reads} - of those
counted, the joins in the gap found window by window - when there are
any}, `settings`, `rank_reached`.

`reads` holds at least `reads_in`, `no_primer`, `too_short`, `clustered`
and `unclustered`, and they balance:

    reads_in = no_primer + too_short + low_quality + clustered + unclustered

(`low_quality` is 0 where a route does not filter on quality). Other
fates are breakdowns within these and say so in `reads_note`
(`unmapped`, `whole_insert`, `merged`, `read2_in_28s`). `clustered` is the sum of
clusters.tsv's `reads`.

`settings` records what the protocol was run with; a protocol adds its
own keys (`nemits2`: `score_over`, `calling_rule`, `gap_edits`,
`gap_exact`, `tentative_identity`, `repeat_aware`, `informal_records`,
`min_span`,
`min_species_span`, `span_slack`, `read_windows`, `primer_offset`).

### `synthesis/db_calls.tsv`

The shared layer's own output, for the report and for anyone joining
across runs: `run_id`, `sample_id`, `sample_type`, `rank`, `taxon`,
`label`, `reads`, `rel_abundance`, `call`, `limit_reads`,
`limit_basis`, `misassign_expected`, `misassign_excess`, `identity`,
`clusters`, `lineage`, `group`, `group_reads`, `group_share`,
`primer_fit`, `margin`, `rank_meaning`, `database_version`; then
`nearest`, `next`, `grade`, `tentative`, `gap_rule_rank`,
`gap_rule_taxon` and `release_lineage` where the clusters carry them (a span-scored
protocol), from each call's largest cluster.

### `synthesis/unplaced_its2.fasta`, `synthesis/not_its2.fasta`

Written by a protocol that declares `unplaced_fasta` (`nemits2`): what
stayed above family, in the two groups the release's notes for amphid
(section 7) ask for. Unplaced ITS2 is a cluster with a database hit
whose call stopped below the family floor ("unclassified, NN% to X");
not ITS2-like a cluster with no hit; chimeras are neither. One record
per sequence across the run, most reads first:
`>unplaced_its2_NNNN reads=N libraries=K nearest=Genus_species
identity=I next=Genus_species+E` (`not_its2_NNNN` carries only reads and
libraries). `db_wells.tsv` then adds `unplaced_its2_fraction` and
`not_its2_fraction` before `database_version`, and `db_run_qc.json`
`marker` and `unplaced` {group: sequences, reads, share, file,
meaning}. Other protocols write neither file and none of these fields.

---

## Layer 3: shared

Synthesis (`synthesise_db.py`: blank limits, the absolute and relative
abundance floors, the misassignment model, plate QC, the read funnel) and
the report (`report_db.py`) read contracts A and B only. What they need to
know of a protocol is a hook on `Protocol`, with a default:

| hook | default | used for |
|---|---|---|
| `organism_group(lineage, host_clade, ranks)` | PR2 groups for PR2's nine ranks, by rank name for NCBI's eight | every call's reporting group |
| `host_clade` | `Craniata` | `--db_host_clade` when not given |
| `min_rel_abundance` | `0.001` | `--db_min_rel_abundance` when not given |

Both protocols take the defaults today. Synthesis resolves the protocol
as assignment does (`--protocol`, `--calling-rule`).

---

## Adding one

**A protocol.** `bin/protocols/<name>.py` with a `Protocol` subclass:
its read types (and why not the others), calling rules, required files,
`classify_consensus` and/or `classify_pairs`, any extra columns and the
hooks it changes. Add its name to `PROTOCOLS`, a fixture to
`tests/test_contracts.py`, and how releases declare it to this document.

**A platform.** An adapter that delivers contract A: per-sample FASTQs
that pre-flight can name, a read type, and the run's accounting (fates in
the shared names, its own summary under `native`). Add it to
`PLATFORMS` and `PLATFORM_READ_TYPES` in pre-flight and to the routing in
`main.nf`.
