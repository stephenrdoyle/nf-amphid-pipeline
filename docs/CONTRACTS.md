# Database mode: layers and contracts

Database mode (`--reference_database`, [REFERENCE_DATABASE.md](REFERENCE_DATABASE.md))
runs in three layers. Each hands the next a fixed contract, so a new
sequencing platform needs only an adapter, a new database only a protocol,
and neither touches synthesis or the report.

```
 PLATFORM (layer 1)           PROTOCOL (layer 2)                SHARED (layer 3)
 get each sample's reads      primer trimming -> calls          synthesis, report

 illumina ──────┐                                           ┌── synthesise_db.py
 nanopore ──────┼─ CONTRACT A ─► vespa_v4 | nemits2 ─ CONTRACT B ─┤   (blank limits, floors,
 predemultiplexed┘  reads,         (bin/protocols/,       clusters.tsv,  │    misassignment, plate QC)
 (pacbio: next)     read type,      chosen by the         taxa.tsv,      └── report_db.py
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
| `pacbio` | (next) | `single_long_accurate` | a built-in symmetric-tag demultiplexer, lima on Linux |

`predemultiplexed` is database mode only. Its samples are found under the
run's `path` as `<prefix>.fastq.gz` (or `.fq.gz`, `.fastq`, `.fq`),
as chunks `<prefix>.*`, `<prefix>_*` or `<prefix>-*` with a FASTQ suffix,
or as every FASTQ in a folder `<prefix>/` (Dorado and MinKNOW's
`barcodeNN/`); `<prefix>` is `file_prefix`, else `sample_id`. lima names
its output `<movie>.<bc>--<bc>`, so a lima run sets `file_prefix` to that.
Indices may be left blank on a predemultiplexed run; a well without them
has no plate position (below).

### Contract A

What every adapter hands the protocol layer, for each run:

1. **Per-sample reads.** One FASTQ per nanopore sample (demultiplexed
   here); R1 and R2 per Illumina sample; one or more FASTQs per
   predemultiplexed sample, read as one. Pre-flight's
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
   | `fates` | {fate: int} | `assigned`, `one_end_only`, `chimeric`, `unmatched`, `internal_adapter`; `{}` when the adapter cannot know (Illumina, predemultiplexed) |
   | `balanced` | bool or null | whether the fates add up to `reads_in` |
   | `per_sample` | {sample_id: int or null} | reads delivered to each declared sample; null where the adapter did not count them (assignment's `reads_in` then holds the count) |
   | `plate` | {sample_id: {row, col}} | plate position of each well that has one |
   | `native` | {demux, split} | the platform's own summaries, verbatim, or null |

   Null is "not known at this layer", never zero. Synthesis reads the
   accounting, not a platform's own files: `native` carries the detail a
   shared fate cannot (nanopore's end reasons, combination grid,
   basecaller), and the read funnel and plate QC are built from it.
4. **Plate positions**, where the platform has them: from the index map,
   the i5 giving the row and the i7 the column (`file_resolution.tsv`
   `plate_row`, `plate_col`, and the accounting's `plate`). The
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
| `nemits2` | nemabiome-its2-refdb, ITS2 | each candidate record scored over its own span, query overhang free; LCA of the ties | `gap` (default), `floors` |

### Which protocol a release gets

1. `--db_protocol`, when it is not `auto`;
2. else the release's `manifest.json` `"protocol"`;
3. else inferred from what manifests already say: a release scored over
   its records' span (`"score_over": "reference span"`), or one called by
   the gap rule (manifest `calling_rule`, or `--db_calling_rule gap`), is
   `nemits2`; everything else `vespa_v4`. So every release built before
   the field existed runs exactly as before.

Pre-flight resolves it the same way, says which in its log, and fails
when the protocol does not exist, cannot take a run's read type, does
not call by the calling rule asked for, or lacks a release file it reads
(`Protocol.required_files`, beyond the release format's own).

### Read types

| | `paired_short` | `single_long_noisy` | `single_long_accurate` |
|---|---|---|---|
| `vespa_v4` | two fixed windows, denoised, pair floors | reference-guided consensus, nanopore floors | not yet: floors not validated on PacBio |
| `nemits2` | overlapping mates merged, else two windows (R2 with its 28S stretch), span-scored | consensus carrying the 28S stretch, span-scored | not yet: gap rule not validated on PacBio |

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
| `label` | what the report says: the taxon, or "unclassified, NN% to X", "X, unresolved below R", "chimera of A and B" |
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
| `chimera_of` | `A+B@breakpoint` for a chimera, with "within genus, kept" for one that is kept |
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
within_genus_kept}, `settings`, `rank_reached`.

`reads` holds at least `reads_in`, `no_primer`, `too_short`, `clustered`
and `unclustered`, and they balance:

    reads_in = no_primer + too_short + low_quality + clustered + unclustered

(`low_quality` is 0 where a route does not filter on quality). Other
fates are breakdowns within these and say so in `reads_note`
(`unmapped`, `whole_insert`, `merged`). `clustered` is the sum of
clusters.tsv's `reads`.

`settings` records what the protocol was run with; a protocol adds its
own keys (`nemits2`: `score_over`, `calling_rule`, `gap_edits`,
`min_span`, `min_species_span`, `span_slack`, `read_windows`,
`primer_offset`).

### `synthesis/db_calls.tsv`

The shared layer's own output, for the report and for anyone joining
across runs: `run_id`, `sample_id`, `sample_type`, `rank`, `taxon`,
`label`, `reads`, `rel_abundance`, `call`, `limit_reads`,
`limit_basis`, `misassign_expected`, `misassign_excess`, `identity`,
`clusters`, `lineage`, `group`, `group_reads`, `group_share`,
`primer_fit`, `margin`, `rank_meaning`, `database_version`; then
`nearest`, `next` and `grade` where the clusters carry them (a
span-scored protocol), from each call's largest cluster.

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
