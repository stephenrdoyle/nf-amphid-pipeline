# Design note: assignment against a reference database

**Status: implemented end to end on nanopore and Illumina** -
assignment, synthesis and report (`--reference_database`, 2026-10-28;
§11 records what was decided and why).

Today the pipeline is given a reference panel - a short, curated list of
targets - and every read is assigned competitively to whichever target it
fits best. The answer can only ever name something on the panel.

This note proposes a second mode, `--reference_database`, in which the
reference is a large taxonomically-annotated database (PR2) and the
question changes from "which of these targets is present?" to "what
eukaryotes are in this sample?". The driving use is the VESPA 18S V4
protocol (Owens et al. 2024, *Nat Commun* 15:402,
doi:10.1038/s41467-023-44521-3), run on both Illumina and nanopore.

---

## 1. Fixed inputs and constraints

| | |
|---|---|
| Primers | `AMPSEQ_VESPA18S_29F` `AGCAGCCGCGGTAATTCC` / `AMPSEQ_VESPA18S_21b8_I_R_Nmod` `TCCGTCAATTYCTTNAASTTTC` (locus-specific parts; both carry the standard Nextera tails). Our reverse primer differs from VESPA's published 21b8R (`TCAATTYCTTIAASTTTC`, Supplementary Table 1): a 5′ `TCCG` extension that matches the template, and N in place of deoxyinosine |
| Library | Two-step Nextera, the same construct the pipeline already handles |
| Database | a release of vespa-18s-refdb (§4): PR2 v5.1.0 backbone plus genome, GenBank and local layers, nine PR2 ranks |
| Illumina | 2x150, fixed for now |
| Nanopore | the existing ligation route, Illumina libraries sequenced on MinION |
| Controls | **no control plasmids and no quantification standards** |
| Amplicon length | not fixed: V4 length varies between lineages |

## 2. What VESPA itself does, and what we take from it

VESPA ran on MiSeq, inferred sequence variants with DADA2, joined read
pairs with `justConcatenate = TRUE` (the pairs did not overlap even at
2x250), and classified with `assignTaxonomy` against PR2 5.1.0. Its
resolution claim - 98.3% of known endosymbiont amplicons resolvable to
species - is from in-silico amplicons collapsed at 100% identity, on
known endosymbionts only.

We take the primers, the database, and the idea of measuring resolution
in silico before trusting it. We do not take the classifier: see §5.

## 3. What does not change

**Nanopore splitting and demultiplexing are primer-agnostic.**
`split_concatemers.py` and `demux_ont.py` locate units by the inner
Nextera anchors (`i5_inner_anchor`, `i7_inner_anchor` in the assay JSON)
and treat everything between them as insert. VESPA's first-round primers
carry the same Nextera tails. Nothing in either script reads the amplicon
primers or assumes an amplicon length beyond the 100bp `min_unit_length`
floor.

**The one panel dependency there is the look-alike mask.** Pre-flight's
`check_anchor_lookalikes` scans the panel for near-copies of the anchors
and writes `anchor_lookalikes.tsv`, which both the splitter and the
demultiplexer use to avoid cutting or quarantining those reads. It exists
because a *Trichuris incognita* amplicon carried the i5 anchor at two
edits and silently lost most of its reads. In database mode the mask
comes with the database release (§4), and an organism absent from the
database cannot be protected by it - so the rate of internal anchor hits
should be reported per run.

*Implemented:* pre-flight copies the release's `anchor_lookalikes.tsv`
to `preflight/anchor_lookalikes.tsv` unchanged, so both stages read it
from the path they always have. A release mask is ~46,000 rows (31,228
distinct contexts over both strands) against a panel's handful, and the
mask's original test - every context aligned at every anchor hit - cost
0.3 s per hit. It is now indexed (§11.1); answers are identical.

## 4. The database is built elsewhere

Building the reference - in-silico PCR over PR2, genome and GenBank
layers, curation - is its own project with its own inputs and release
cycle: the **vespa-18s-refdb** repository, beside this one. Its design
note holds the build and the measurements. What this pipeline needs to
know from them:

- **Inserts are long.** Median 546bp between the primers; cestodes reach
  830-970bp. At 2x150, 99.9% of inserts leave an unread gap between the
  reads - the Illumina route cannot merge (§5).
- **Resolution is a property of the release.** Helminth species-level
  resolution is roughly 78-89% over the whole insert (nanopore) and
  65-78% over the 2x150 windows, best case. The release carries a table
  for each, and the report states the rank each call actually reached.
- **The look-alike mask matters.** Nearly a thousand unique amplicons in
  PR2 carry an adapter anchor within the splitter's tolerance, 50 of them
  helminth. The release ships the mask in pre-flight's own format.
- **The database has gaps that will not close soon** (*Teladorsagia*,
  *Cooperia*, *Ancylostoma*, *Onchocerca* among them), which is why §4.2
  below is part of the design, not an afterthought.

### 4.1 Release format (the contract between the two repositories)

A release is one directory, named for its build version, holding:

| file | contents |
|---|---|
| `refdb.fasta` | one record per unique insert, primers trimmed; the header is the record id only |
| `taxonomy.tsv` | id, the nine PR2 ranks (merged to the ranks all members share; placeholder ranks left blank), the source layers contributing, and the number of member sequences |
| `members.tsv` | every source sequence collapsed into each id: accession or genome and contig, layer, original label |
| `resolution_nanopore.tsv` | per id, the deepest rank at which it is distinguishable over the whole insert |
| `resolution_illumina_<R1>x<R2>.tsv` | the same over the two read windows for that read length |
| `anchor_lookalikes.tsv` | adapter-anchor near-copies, in exactly the columns the pipeline's pre-flight writes (`seq_id anchor strand position edits context`), so the pipeline reads it unchanged |
| `conflicts.tsv` | identical sequences whose labels disagree at a high rank, and what the build did with each |
| `primer_fit.tsv` | per id, how well the primers fit its best-matching template: forward and reverse edits, mismatch positions counted from each primer's 3′ end, the template at each site, and a class - `exact`, `away_3prime`, `near_3prime` (a mismatch in the last 5 bases), `indel`, `edits_only` (a genome-only record: counts known, sites not), `unknown` (from vespa-18s-refdb v0.4.0) |
| `manifest.json` | build version, primers, PR2 version, WBPS release, version of every other layer, build date and code commit |

The pipeline records `manifest.json`'s build version in every result.

### 4.2 Gaps must be visible, not hidden

Some of what the primers amplify will not be in any database. So:

- **Identity thresholds per rank.** A read is assigned to species only
  if it is close enough; otherwise it stops at genus, family or above.
  Nothing is forced onto the nearest entry.
- **Unassigned and poorly-matched reads are clustered** per run and
  reported as clusters with a consensus, a read count and the best hit -
  "unclassified, 92% to *Trichostrongylus*", not a false
  *Trichostrongylus* call.
- **Confirmed clusters go back to vespa-18s-refdb** as local
  sequences, so each run improves the reference instead of repeating
  the same misreport.

## 5. Assignment

One rule for both platforms: **map to the trimmed database, then resolve
to the lowest common ancestor.** A rank is reported only when every
equally good hit agrees at it and identity clears that rank's threshold.
Using `assignTaxonomy` on Illumina and mapping on nanopore would mix
method with chemistry in every cross-platform comparison.
`assignTaxonomy` may be run on Illumina as an optional cross-check
against VESPA's own output.

**Illumina.** At 2x150 the pairs cannot overlap for any V4 insert of
realistic length, so merging is not an option - and dropping pairs that
fail to merge would bias silently against long-V4 lineages. Database mode
gets a non-merging route: trim the primer from each read, filter each
read, dereplicate pairs, map pairs. The insert between the windows is
never observed; the Illumina resolution table says what that costs.
Poly-G tails must be trimmed if any run is on a two-colour instrument.
The existing merge-based route is unchanged for panel mode.

**Nanopore.** The whole insert is in one read. ~~Per-read mapping with
LCA first~~ - **decided against on the benchmark (§11.2): nanopore
classifies per-sample consensus sequences, not reads.** Per-read error
does exceed the divergence between sister taxa, and the result is worse
than stopping at genus: reads from species *missing* from the database
are confidently named as a present sister.

**Thresholds** per rank: species 0.995, genus 0.99, family 0.95 on the
best hit's identity, from leave-out simulations at consensus error rates
(§11.3). The genus floor is the weak one.

## 6. No controls

- Output is **relative 18S abundance**, never organisms or copies. 18S
  copy number varies widely between taxa (VESPA's own observation, and
  why they built cloned mocks), so read share is not organism share.
- The control-tag check is switched off.
- Negative controls (extraction blanks, no-template) become the main
  defence against contamination and cross-talk; database mode should
  require them.
- Calls are judged against the negative controls and a relative-abundance
  floor, not a fixed depth alone.
- ~~Open: whether pre-flight and assignment accept a run with no control
  wells.~~ *Resolved:* in database mode pre-flight skips every panel and
  tag check (there is no panel), refuses `--spiked` and `--mocks` by
  name, and **warns** - does not refuse - when the plate has no negative
  controls: refusing to analyse reads already on disk helps nobody. The
  assay `assets/assay/vespa18s.json` has an empty control tag. Whether a
  missing blank should become an error is open (§10).

## 7. Report

- A taxon table per sample, at the **rank actually reached** for each
  call, with the database version and primers recorded.
- **Predicted amplification**, from the release's `primer_fit.tsv`: a
  call on a record classed `near_3prime`, or the absence of a taxon whose
  records all are, says that the primers were predicted to fit poorly -
  so a low or missing signal may be the PCR, not the sample. It matters
  most for flatworms: about two thirds of their records carry a mismatch
  4 bases from the forward primer's 3′ end (vespa-18s-refdb, docs §6.4).
- Off-target fractions per sample - bacteria and archaea, host, fungi -
  which PR2's taxonomy gives directly and which matches how VESPA
  reports.
- Cross-platform and cross-run comparisons made at the rank both sides
  can reach for that taxon, read from the two resolution tables.

## 8. Validation

- **Illumina:** VESPA's public reads (SRA PRJNA944233), including the
  Equimolar and Log EukMix mocks, whose composition is in VESPA
  Supplementary Table 4. Score against both the mock truth and VESPA's
  published results.
- **Nanopore:** no public equivalent. Needs either a characterised
  sample or VESPA-style libraries run on MinION.
- **Simulator:** draw inserts from the trimmed database, deliberately
  including sister taxa as decoys, and score precision and recall
  separately at each rank, so it is visible where resolution runs out.

## 9. Order of work

1. ~~A first database release from vespa-18s-refdb (§4).~~ v0.4.0.
2. ~~Illumina pair mapping with LCA~~ - built (§11.5). Validated against
   VESPA's own published per-sample results; the mocks' reads were not
   deposited.
3. ~~Nanopore per-read mapping with LCA~~ - benchmarked and replaced by
   4 (§11.2).
4. ~~Nanopore consensus~~ - built (§11), validated on simulation; still
   needs a characterised nanopore sample.
5. ~~Synthesis and report for database mode (§§6-7)~~ - built (§11.4).

## 10. Open questions

- ~~Aligner~~ - minimap2 (§11). ~~Identity thresholds per rank~~ - see
  §11.3; the genus floor still needs work.
- ~~Pre-flight and assignment behaviour with no control wells~~ - §6.
- Should a database-mode plate with no negative control be refused
  rather than warned about?
- The splitter's look-alike test costs ~2.4 ms per anchor hit with a
  release mask (two hits per read), against nothing with a panel's.
  Measured on VESPA-shaped reads, the splitter's anchor scan runs at
  ~195 reads/s per core with the v0.4.0 mask against ~3,300 without:
  on a 3.2M-read run about 4.5 CPU-hours, ~35 minutes on an 8-core
  laptop.
  Restricting the release mask to look-alikes within 2 edits (1,038
  rows, against 45,301 at 3) would cut that ~40-fold, but 3 edits is
  pre-flight's deliberate margin for sequencing error.
- Unclustered reads (groups below `--db_min_cluster_reads`) are counted
  but not classified. Whether to give them a coarse per-read call - at
  family or above only, where read error does not reach - is open.
- No container image carries both python3 and minimap2, so the two
  database-mode processes declare none.
- The consensus cannot extend past its starting reference's ends: a
  taxon whose insert is longer at an end than every database record
  near it loses those bases (soft-clipped). Inserts are primer-to-primer
  on both sides, so this needs an indel right at an end.

## 11. Implementation, nanopore (2026-10-28)

### 11.1 Aligner and index

**minimap2** (2.31, built from source for Apple Silicon): the pipeline's
first external dependency, used in database mode only; panel mode stays
standard-library Python. Indexed with `-x map-ont -k 19 -w 10` - the same
accuracy as the map-ont default k=15 on the leave-out benchmark, 1.6x
faster. The index builds in 1-2 s and ~120 MB, so `DB_INDEX` builds it
once per run rather than shipping it in the release. Pre-flight checks
that minimap2 runs and records its version.

**The look-alike mask index** (`ampseqlib._ContextIndex`). An exact
filter in front of the unchanged per-context test, by the pigeonhole
principle: cut a context of tolerance k into k + 2 pieces; an alignment
within k edits leaves two pieces exact, on read diagonals within k of
each other, and - since covers() only asks about alignments overlapping
the anchor hit - on diagonals between `start - m - 2k` and `end + k`.
Pieces with ambiguity codes are indexed under every spelling (PR2 keeps
its Ns; ~1,200 contexts carry one), and so are read pieces with one.
Diagonal pairs are found with sorted packed integers so the inner loops
run in C. 0.3 s → ~2.4 ms per hit on both panel reads and VESPA-shaped
reads; the overlap bound matters on the latter, where the 18S just past
the primer genuinely matches hundreds of insert-start contexts. Tested
call for call against the plain loop (`tests/test_split.py`).

### 11.2 Classify consensus sequences, not reads

Benchmarked with minimap2 on reads simulated from the release with 40
genera and 60 species held out of the index:

- **Single reads:** an LCA over equally best hits named 79-85% of reads
  from *missing* species or genera at species level. No identity floor
  can fix it: at 3% read error a present species' reads sit at ~97%
  identity and a missing sister's at ~96%, because sister taxa differ
  by ~1% over V4. The distributions overlap at 1% read error too.
- **Consensus (~0.2% error):** a species floor of ~0.995 kept 90% of
  correct species calls while wrong calls on missing taxa fell to 5-8%.

So `bin/assign_db.py`, per sample:

1. trims reads at the primers (orienting them; the release's records
   are primer-trimmed inserts);
2. maps them to the release and groups them by primary hit;
3. builds a reference-guided consensus per group of ≥5 reads: pile the
   reads on the hit via minimap2 `cs` tags, majority base per position
   (insertions where most spanning reads carry one), re-align to the
   draft, repeat until stable (≤4 rounds). Groups over 300 reads are
   subsampled for the consensus only;
4. **splits hidden haplotypes by linkage** (Stephen Doyle's idea): a
   candidate allele must be ≥10% of the reads at its position and
   binomially improbable under the error rate *measured on that group*,
   not a homopolymer deletion; two alleles are linked when ≥3 reads
   carry both and that is ≥3x their chance co-occurrence; the largest
   connected set of linked alleles (≥2) is the signature, reads carrying
   at least half of it form a second group, and both are resolved again
   (at most two levels deep);
5. merges consensus sequences within one edit;
6. counts **every** read by mapping it to the final consensus set, so
   reads from groups too small to have a consensus are still counted -
   unless their own database hit fits better by more than 10 alignment
   score points, which marks a taxon too rare here for a consensus; those
   are left *unclustered* rather than absorbed into a relative;
7. classifies each consensus: LCA over database hits within 0.5% of the
   best alignment score (`--db_lca_margin`), cut back to the deepest rank
   whose identity floor the best hit clears (`--db_floors`).

On simulated samples (prototype): splitting raised hidden-haplotype
recovery from 0% to 37% and exact consensus for present taxa from 87% to
98%, at about three false splits per 30 samples; merging cut spurious
consensus sequences from 33 to 5-10. On the pipeline fixture
(`tests/make_db_fixture.py`) every source is recovered at its expected
rank, with exact consensus.

**Identity** is matches over the alignment block plus any query left
unaligned at either end, not minimap2's block-only ratio: both sequences
are primer-to-primer, so an unaligned end is a real difference.

**Below family** a consensus is reported as `unclassified, NN% to X`
(X the best hit's genus), with its rank reached, lineage and identity
still in the row - a gap in the database reads as a gap.

**Output** per sample (`assign_db/`): `clusters.tsv` (one row per
consensus: reads, relative abundance, rank, taxon, label, identity,
coverage, lineage, best hit and its full lineage, near-hit count, LCA
rank before floors, the best hit's resolution rank from
`resolution_nanopore.tsv`, origin, consensus), `taxa.tsv` (reads per
taxon at the rank reached; a species call and a genus call in one genus
stay separate rows), `consensus.fasta`, and `assign_db.json` (read
accounting, settings, release and minimap2 versions). Every row carries
the release version.

### 11.3 Floors

Prototype values, from the leave-out benchmark at consensus error rates:
species 0.995, genus 0.98, family 0.95. The genus floor is too
permissive: 17-37% of consensus sequences from missing genera were still
named a wrong genus. `tests/bench_db_floors.py` reruns the calibration
on a larger leave-out simulation through `bin/assign_db.py` itself, with
floors off, and reports a grid of floor values.

**Larger leave-out run (2026-10-28).** v0.4.0 with 100 genera and 150
species held out, 150 samples of four sources each (two present, one
missing species, one missing genus; 20-400 reads each, `nanopore_hac`
errors), 649 consensus sequences of which one matched no source.
Share of consensus sequences per outcome (deepest = the deepest rank it
may honestly reach; wrong = deeper than that, or a wrong name):

| floors (species / genus) | present: deepest / wrong | missing species: deepest / higher / wrong | missing genus: deepest / higher / wrong |
|---|---|---|---|
| 0.990 / 0.98 | 0.98 / 0 | 0.16 / 0.62 / 0.22 | 0.57 / 0.16 / 0.27 |
| 0.995 / 0.98 (prototype) | 0.98 / 0 | 0.24 / 0.62 / **0.15** | 0.57 / 0.16 / **0.27** |
| **0.995 / 0.99** (default) | 0.98 / 0 | 0.18 / 0.68 / 0.14 | 0.64 / 0.16 / 0.20 |
| 0.997 / 0.99 | 0.98 / 0 | 0.24 / 0.68 / 0.09 | 0.64 / 0.16 / 0.20 |

(full grid in the benchmark's `floors.tsv`). Read with care:

- **Present taxa are untouched by every floor tried**, because simulated
  consensus is exact (identity 1.000). Real consensus will carry some
  residual error, which is what a floor of 0.997 (~1.6 differences over
  550bp) would start to cost; that needs a characterised nanopore run.
- **Raising the genus floor to 0.99 is free here and cuts wrong calls on
  missing genera from 27% to 20%** - and moves them to *family* (deepest
  0.57 → 0.64), which is the honest answer.
- What remains wrong does not respond to floors: missing taxa with a
  relative >99% identical over V4. That is the resolution limit, and the
  identity column beside each call is how it stays visible.

**Decided (2026-10-28): genus raised to 0.99.** Species stays at 0.995
until a real nanopore run shows how much residual error consensus
carries; 0.997 would cut wrong calls on missing species from 15% to 9%.

### 11.4 Synthesis and report

Database mode has its own pair, `bin/synthesise_db.py` and
`bin/report_db.py` (sharing the panel report's style), because the panel
synthesis is built on standards and plasmids this protocol does not have.

- **Unit.** Relative 18S abundance - each call's share of its well's
  classified (clustered) reads - and the words say so on every page.
- **Rank.** Calls stay at the rank they reached. A species call and a
  genus call in one genus are separate rows; nothing is rolled up.
- **Calls** (`synthesis/db_calls.tsv`): the panel's blank-derived limit,
  per (rank, taxon): more reads than any blank showed (absolute), more
  than the largest share of any blank at least `--min_depth` deep
  (proportional), and at least `--absolute_floor`; then at least
  `--db_min_rel_abundance` (0.1%) of the well. A well under `--min_depth`
  classified reads calls nothing (`low_depth`). Without blanks the limit
  is the fixed floor and every row says that is an assumption.
- **Organism groups**, reported the way VESPA reports them (Owens et al.
  2024, Figs 3-5): each well's classified reads by group first, then each
  *target* group on its own with shares recalculated within the group, a
  presence table per taxon and sample, and prevalence and richness. The
  group comes from the PR2 lineage down to the rank the call *reached* -
  a call that stopped at domain is `unassigned`, not guessed - first
  match wins: bacteria and archaea (domain); organelles (PR2's
  `Eukaryota:plas`, `:mito`, `:nucl`, `:apic`, `:chro`); host
  (`--db_host_clade`, Craniata - vertebrates - by default); **helminths**
  (class Nematoda, Platyhelminthes or Acanthocephala); **microsporidia**
  (order Microsporida); fungi (subdivision Fungi); other invertebrates
  (any other Metazoa); plants (supergroup Archaeplastida); **protists**
  (any other eukaryote). The bold three are the target groups.
  `db_calls.tsv` carries each call's group and its share within the
  group; `db_wells.tsv` reads, fractions and target-group richness per
  well; `db_prevalence.tsv` samples detected per taxon.
- **Predicted amplification:** a detected call whose best hit is classed
  `near_3prime` in `primer_fit.tsv` is flagged in the table and listed in
  the report.
- **Internal adapters:** `db_run_qc.json` carries the splitter's cut count
  and the demultiplexer's internal-adapter quarantine (§3), with the
  mask's size, and the report puts them before the results.
- **Cohort:** a sample with demultiplexed reads but no assignment stops
  synthesis; a well with no reads is reported empty.
- The report puts the release, the blanks, depth, internal adapters and
  primer-fit warnings first; then all reads by group; then helminths,
  protists and microsporidia, each with its within-group composition,
  presence table and prevalence; then off-target, naming depth and every
  call.

Release note for vespa-18s-refdb: v0.4.0's `taxonomy.tsv` has one row
whose domain is the literal `domain` and nine with no domain.

### 11.5 Illumina

`bin/assign_db_pairs.py`, as `ASSIGN_DB_PAIRS`. Its outputs have the
nanopore route's shape, so synthesis and the report read either platform.

- **Not merged.** At 2x150 no V4 insert overlaps, and dropping pairs
  that fail to merge would bias against long-V4 lineages. Each mate is
  primer-trimmed (either mate may carry the forward primer) and cut to a
  fixed window of insert - the release's `illumina_windows`, 132 and 128
  for 2x150, or `--db_window_r1`/`--db_window_r2` - so error copies line
  up base for base. A pair whose insert ends inside R1 (read-through,
  found with prep_illumina's junction search) is kept whole as one
  sequence, so short-V4 lineages are not lost either.
- **Denoised, not consensus-built.** Illumina's error is in the counting,
  not the bases: UNOISE-style, a sequence d mismatches from a more
  abundant one and at most 1/2^(2d+1) of its abundance is its error copy;
  sequences seen once join the nearest survivor within the radius or
  stay unclustered. The radius is 1.5% of the sequence, at least 3
  mismatches: a fixed 3 left 28% of VESPA's 2x300 pairs unclustered,
  against 9% at 2x150. Candidates are found exactly by pigeonhole pieces,
  and mismatches are counted with integer XOR (same answer, ~30x faster).
- **Leftovers.** What no sequence took within the radius first joins one
  within twice it (noisy copies of abundant taxa - 2x300 R2 tails), then a
  rescue pass clusters the rest among themselves: the leftover with the
  most others within the radius founds a cluster of them, its sequence
  their per-position majority. A rare taxon whose error-prone reads never
  repeat exactly (VESPA's *Endolimax*) otherwise vanished whole. On the
  largest VESPA sample (224k pairs) unclustered pairs fell from 41k to
  10k and the sample takes ~5 minutes.
- **Classified as a pair.** Each window is mapped (minimap2 `-x sr`,
  secondary hits within half the best score); a record counts when both
  windows hit it in a pair's shape (R1 forward from its start, R2 reverse
  from its end), scored by the two alignment scores together; identity is
  over both windows. When no record carries both, each window is
  classified alone and the call is the LCA of the two, at the lower
  identity (`origin` = `pair_split`).
- **Floors of its own** (`--db_floors_illumina`, species 0.997, genus
  0.995, family 0.95). The bases that separate sister taxa are mostly in
  the gap a 2x150 pair never reads, so the nanopore floors let a missing
  species take its sister's name. Leave-out calibration at 2x150 (same
  design as §11.3, `tests/bench_db_floors.py --illumina 150`, 613
  sequences):

  | floors (species / genus) | present: deepest | missing species: wrong | missing genus: wrong |
  |---|---|---|---|
  | 0.995 / 0.99 (nanopore's) | 0.83 | 0.16 | 0.18 |
  | **0.997 / 0.995** | 0.81 | **0.09** | **0.17** |

  At 250x200 windows (2x301 reads, `--illumina 301 --windows 250x200`)
  the same floors are again best: present taxa 0.95 at the deepest rank,
  wrong on missing species 0.12 and on missing genera 0.18 - so one set
  of Illumina floors serves both read lengths.

  17-19% of present taxa stop above species whatever the floor: the
  windows do not separate them from a relative in the database - the
  2x150 resolution limit, as the release's `resolution_illumina_150x150`
  table predicts.
- **SRA-Lite quality.** Public SRA data often keep two quality values
  (Q30 pass, Q3 failed read); the expected-error filter then only drops
  failed reads, and the summary says the filter was blind.
- **LCA and unlabelled records.** A record with no taxonomy at all (v0.4.0
  has nine, unresolved conflicts) is left out of the LCA rather than
  vetoing it - on both platforms.

Release notes for vespa-18s-refdb found on the way: `V4_ef59fbd391ef`
(PR2 `QJVE01000361.1.1814_U`, labelled *Arthrobacter stackebrandtii*)
is human 18S - a bacterial assembly carrying host contamination - and
ties with *Homo sapiens*, so human reads stop at no rank.
