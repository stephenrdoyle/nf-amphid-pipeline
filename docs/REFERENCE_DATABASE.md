# Design note: assignment against a reference database

**Status: proposal. Nothing here is implemented.**

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

**Nanopore.** The whole insert is in one read. Per-read mapping with LCA
first. Per-read error is likely to exceed the divergence between some
sister taxa, so reads from close relatives will often stop at genus or
family. If the mock data show that this costs too much, add a
cluster-and-consensus step per sample and classify consensus sequences.

**Thresholds** per rank are not chosen here. They should come from the
database preparation's view of how far apart neighbouring taxa actually
are over V4.

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
- Open: whether pre-flight and assignment accept a run with no control
  wells. `--spiked` off already gives presence/absence without
  quantification standards; the tagged-plasmid requirement still needs
  checking.

## 7. Report

- A taxon table per sample, at the **rank actually reached** for each
  call, with the database version and primers recorded.
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

1. A first database release from vespa-18s-refdb (§4).
2. Illumina pair mapping with LCA, validated on the VESPA mocks.
3. Nanopore per-read mapping with LCA, validated on a mock run.
4. Nanopore consensus, only if step 3 shows it is needed.

## 10. Open questions

- Aligner, and identity thresholds per rank (from the release's
  resolution tables).
- Pre-flight and assignment behaviour with no control wells.
