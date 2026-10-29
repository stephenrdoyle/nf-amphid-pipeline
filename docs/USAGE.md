# amphid — usage

Targeted amplicon panel analysis for Illumina and nanopore. One
invocation handles one sequencing run on either platform: the `platform`
column of your inputs decides which branch it takes. To analyse several
runs, invoke the pipeline once per run (see §2.3).

Every column list and default below was read out of `bin/preflight.py`
and `nextflow.config` rather than written from memory. If the two ever
disagree, the code is right and this file is stale.

---

## 1. Quick start

```bash
NF=/nfs/users/nfs_s/sd21/lustre_link/software/nf-amphid-pipeline

nextflow run $NF -resume -bg \
    --reference   reference.tsv \
    --samplesheet samplesheet.tsv \
    --runs        runs.tsv \
    --primers     primers.tsv \
    --spiked \
    --outdir      results
```

**`--spiked` only if the quantification standards went into every
sample.** Leave it off for an unspiked run: calls are then presence /
absence with no copy numbers. See §3.8.

**LSF is the default, not opt-in.** A bare `nextflow run` submits to the
farm. `-profile local` is the explicit opt-out, for a laptop or the head
node. This is deliberate: the alternative — LSF hidden behind a profile —
means every forgotten flag silently runs the whole cohort on the head
node, which happened twice on the first day this workflow was used.

`-bg` detaches the Nextflow head process, which otherwise has to survive
the whole run. Without it, run inside `screen` or `tmux`: an ssh
disconnect kills the head process and orphans the LSF jobs.

`-resume` picks up where a killed or failed run stopped. Splitting,
demultiplexing and assignment are all resumable at task granularity.

### 1.1 In the field: a laptop, no farm, no network

Every tool in `bin/` is standard-library Python 3, so a laptop needs only
**Python 3**, **Java 17+** and **Nextflow** - no conda, no containers, no
plugins. Set it up **while you still have internet**: Nextflow downloads
its own runtime the first time it runs, and cannot do so offline.

Before leaving:

```bash
nextflow -version                                   # fetches the runtime; needs network
nextflow run . -profile test,field                  # needs the fixture reads, see the README
nextflow run . -profile test_ont,field
```

In the field:

```bash
export NXF_OFFLINE=true
export NXF_ANSI_LOG=false
nextflow run . -profile field \
    --reference reference.tsv --samplesheet samplesheet.tsv \
    --runs runs.tsv --primers primers.tsv --spiked --outdir results
```

`NXF_OFFLINE=true` stops Nextflow checking for updates, which offline
otherwise waits for a timeout. `NXF_ANSI_LOG=false` swaps the live,
redrawing progress display for plain lines - worth setting on the farm
too. The live display redraws everything printed during a run, and once
pre-flight's warnings plus the progress table are taller than the
terminal, the warnings are printed again on every refresh. Plain output
shows them once, at the start, and reads the same in a log file or after
an ssh drop. Put both lines in `~/.zprofile` so every terminal has them. `-profile field` runs on the laptop's own
cores with no farm time limits, deletes `work/` once a run succeeds, and
hard-links results instead of copying them, so a run is not stored
twice. Keep `--outdir` on the same disk as the folder you launch from:
hard links cannot cross disks.

**Database mode needs minimap2** (`--reference_database`, §2.9) - the
pipeline's one external program, used in that mode only; panel runs need
nothing beyond Python. On an Apple Silicon laptop, build it from source
before leaving (it needs Xcode's command-line tools, `xcode-select
--install`, and nothing else):

```bash
git clone https://github.com/lh3/minimap2 && cd minimap2
make arm_neon=1 aarch64=1                           # Apple Silicon; plain `make` on Intel
mkdir -p ~/bin && cp minimap2 ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zprofile  # then open a new terminal
minimap2 --version                                  # 2.31 is what database mode was built against
```

Pre-flight checks that `minimap2` runs before any read is touched, and
records its version. Copy the database release folder onto the laptop
too, whole: pre-flight checks every file against the hashes in its
`manifest.json`, so a half-copied release is refused rather than used.
Then, before leaving:

```bash
python3 tests/make_db_fixture.py
nextflow run . -profile test_db,field
```

Plan disk for a split nanopore run larger than the run itself while it
is in progress. Splitting is the slow stage on a laptop; it divides by
MinKNOW batch file, so all cores are used. Basecalling quality matters
more in the field than anything here: on simulated plates, FAST rather
than HAC basecalling cost about five points of correctly resolved
barcode ends and doubled the wrong-sample calls.

---

## 2. Input files

Four files are required for a normal run: a reference panel, a
samplesheet, a run manifest, and primers. Everything else has a default
or is optional.

Tables are read by header name, not by column order. Extra columns are
ignored, so you can keep working notes in them.

### 2.1 `--reference` — the panel (TSV, **required**)

What the reads are competitively mapped against. Targets, control
plasmids, quantification spike-ins and vector backbone all live in one
file, because they compete for the same reads and the assignment is only
meaningful if they are searched together.

| column | required | meaning |
|---|---|---|
| `seq_id` | yes | unique identifier for this sequence |
| `locus` | yes | which marker, e.g. `ITS2`, `cox1` |
| `type` | yes | `target`, `control`, `quant`, or `vector` |
| `sequence` | yes | the sequence itself; IUPAC degeneracy allowed (`ACGTRYSWKMBDHVN`) |
| `taxon` | no | species name for reporting; falls back to `seq_id` |
| `purpose` | no | `speciation` or `variants`. Not exclusive |
| `parent_id` | no | links a haplotype to its parent sequence |
| `default_conc` | no | input concentration, for quantification |
| `conc_units` | no | units for `default_conc` |
| `corr_factor` | no | amplification efficiency correction; defaults to 1.0 |

On `corr_factor`: pre-flight warns if some quantifiable entries carry a
measured value and others do not. Mixing them puts corrected and assumed
copy numbers in the same column, which is only comparable if the
unmeasured ones really do amplify at unity.

### 2.2 `--samplesheet` — one row per library (TSV or CSV, **required**)

A samplesheet row is a *library*, and it is frozen at prep time.

**TSV is the default; CSV still works.** Every table in this pipeline
picks its delimiter from the file extension — `.csv` means comma,
anything else means tab — so the two are interchangeable per file and a
sheet already written as CSV needs no conversion. The rule lives in one
place, `ampseqlib.delimiter_for`, because a reader and a writer that
disagree about a delimiter do not raise: they produce one wide column
whose header matches nothing, and the error that comes back names every
required field as missing.

| column | required | meaning |
|---|---|---|
| `run_id` | yes | must match a `run_id` in the run manifest |
| `platform` | yes | `illumina` or `nanopore` |
| `sample_id` | yes | unique within the run |
| `index_i7` | yes | i7 index sequence or name |
| `index_i5` | yes | i5 index sequence or name |
| `sample_type` | yes | `unknown`, `neg_extraction`, `neg_ntc`, `mock`, `pos_control` |
| `subject_id` | no | who the library came from. Groups replicates with no metadata file needed; also the join key for a subject-keyed `--metadata`. Defaults to `sample_id` |
| `file_prefix` | no | overrides filename matching for Illumina file resolution |
| `expected_taxa` | no | on a `pos_control` row: what it should contain, `;`-separated target taxa or seq_ids |

The `sample_type` vocabulary is **not** configurable. The logic depends
on it:

| `sample_type` | what it is | what the pipeline does with it |
|---|---|---|
| `unknown` | a field sample | called, quantified, summarised |
| `mock` | a known mixture of the tagged control plasmids | processed like a field sample, with its plasmids called as its content (§3.7) |
| `pos_control` | a deliberately positive sample | checked for detection: every `expected_taxa` entry must be detected, or - without that column - any target. A failure is raised before the results. Kept out of field-sample summaries and plate QC |
| `neg_extraction` | an extraction blank | sets the detection limits; reported apart from NTCs, since it carries the kit and handling |
| `neg_ntc` | a no-template PCR control | sets the detection limits; reported apart, since it carries the PCR reagents and carryover |

Any negative carrying more than `--absolute_floor` reads of a reference is
named in run QC and the report: a negative can never be called detected,
since the limits are built from the negatives, so contamination would
otherwise show only as a quietly raised limit. When the samplesheet
declares any controls, the report gains a **Controls** section with each
positive's result and the negatives by type. Pre-flight checks that every
`expected_taxa` entry is a target in the panel.

Leaving a few index combinations unused is worth doing. If every
combination in the pool is assigned, no sentinel exists and index
misassignment becomes *unmeasurable* — which is a different claim from
measuring it to be zero. Pre-flight warns when this happens.

### 2.3 `--runs` — one row per sequencing run (TSV, **required**)

Where the reads are. On nanopore the reads belong to the *run*, not to
any one sample, which is why this file exists separately.

**One run per invocation.** Pre-flight refuses a manifest with more
than one row. Synthesis, detection limits, plate QC and the report are
not yet separated by run, so several runs would be pooled into one
result - blanks from one plate setting limits on another - and named
after the first. Split the manifest, run each separately, and pair the
results afterwards with `--comparison`.

| column | required | meaning |
|---|---|---|
| `run_id` | yes | matches the samplesheet |
| `platform` | yes | `illumina` or `nanopore` |
| `path` | yes | Illumina: directory of per-sample FASTQ. Nanopore: MinKNOW run folder, or a pooled FASTQ |
| `barcode_repeats` | no | tandem copies of each index (1 or 2). Defaults to the assay config |
| `dna_input_ng` | no | DNA input, for the quantification chain |
| `lot_*` | no | any column starting `lot_` is recorded as reagent provenance |

**A relative `path` is relative to this file**, not to the directory you
launched from. Pre-flight resolves it and writes the absolute version to
`preflight/runs_resolved.tsv`, which is what the workflow then uses.

The `platform` column must say the same thing here and in the
samplesheet. If they disagree, one branch sees no rows and the other
sees no run — and the run now fails with `no sample reached assignment`
rather than completing having done nothing.

For nanopore, `fastq_pass` is used and `fastq_fail` is excluded. Those
reads failed the basecaller's own quality filter, and including them
silently would import a different error distribution into a
demultiplexer whose barcode tolerances were set for passing reads. The
exclusion is counted and reported rather than hidden. `--include_failed`
overrides this.

### 2.4 `--primers` — PCR1 primers (TSV, required)

| column | required | meaning |
|---|---|---|
| `name` | yes | primer name |
| `direction` | yes | `fwd` or `rev` |
| `sequence` | yes | primer sequence, IUPAC allowed |

Used to trim primer sequence off the panel, producing
`preflight/reference.insert.fasta`, which is what every read is
assigned against on both platforms. The workflow refuses to start
without `--primers`, and pre-flight fails if the primers are located in
no reference entry - in both cases there would be no trimmed panel, and
the run would otherwise skip assignment for every sample.

### 2.5 `--metadata` — grouping metadata (TSV, optional)

Kept deliberately *out* of the samplesheet. A samplesheet row is a
library, is frozen at prep, and travels to the sequencing centre and
into every results directory; metadata keeps changing and should not
ride along by default. Pre-flight validates the join rather than
assuming it.

Every field in this file becomes a grouping in `group_summary.tsv`
(§5), so this is how you ask "how many reads per arm".

**The key column declares the level, and there must be exactly one:**

| key column | the row describes | applies to |
|---|---|---|
| `subject_id` | a subject or isolate | every library made from it |
| `sample_id` | one library | that library alone |

| column | required | meaning |
|---|---|---|
| one of the two above | yes | the join key |
| anything else | no | your own fields, each of which becomes a grouping |

Which one you want is decided by the variable, not by preference. A
country or a population is a property of the subject. A **bench
treatment is not**: testing uncleaned against AMPure-cleaned PCR1
splits one subject's product into two libraries that differ only in
the treatment, so keyed on `subject_id` both arms would get the same
value and the comparison would come out empty with nothing saying why.
That file has to be keyed on `sample_id`.

```tsv
sample_id	pcr1_cleanup
iCCH-001-r1	none
iCCH-001-r2	ampure
```

A file with both key columns is refused: it gives a library two
sources for one field and no rule for which wins.

Replicates need no file at all. Several libraries declaring the same
`subject_id` in the samplesheet are already a grouping, and are
summarised as one whenever at least one subject has more than one
library.

A library with no metadata row is grouped under `(no value)`, never
dropped — so the library counts within a field always add up to the
plate, and a small group is distinguishable from a broken join.

**Refused as a hard error**, not a warning: `latitude`, `longitude`,
`lat`, `lon`, `lng`, `gps`, `coordinates`, `coord`, `easting`,
`northing`, `postcode`, `zip`, `address`, `household_gps`,
`cluster_gps`. Precise location plus infection status is re-identifying
for participants who cannot consent to that, and a results directory is
copied around freely.

### 2.6 `--mocks` — mock community composition (TSV, optional)

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | which well holds this mock |
| `seq_id` | yes | which reference entry |
| `amount` | yes | how much went in |
| `amount_units` | yes | units for `amount` |
| `mix` | yes | which mix this row belongs to |

Supplying this enables `CALIBRATE_MOCKS`, which recovers per-reference
amplification efficiencies and checks them for nonlinearity.

### 2.7 `--index_map` — index primer map (TSV, has a default)

Defaults to `assets/index_map.tsv`. Columns: `role`, `name`, `sequence`,
`position`. This is what lets plate position be derived from the index
pair (i5 → row, i7 → column), which in turn makes the plate-geometry QC
sensors possible.

### 2.8 `--assay` — assay definition (JSON, has a default)

Defaults to `assets/assay/ampseq_helminth.json`. Holds the control tag,
library construct, indexing design and tolerances — everything that makes
this *a particular assay* rather than code. Pointing the pipeline at a
different panel or construct is a file swap, not a code change. Read by
both Nextflow and `preflight.py`, so the values are defined once.

### 2.9 `--reference_database` — a database release instead of a panel (directory)

For protocols with no panel, no control plasmids and no quantification
standards - the VESPA 18S V4 protocol is the one it was built for. Give a
release directory built by the vespa-18s-refdb repository *instead of*
`--reference`, with `--assay assets/assay/vespa18s.json` and the
protocol's primers. The design, and why each step is there, is in
[REFERENCE_DATABASE.md](REFERENCE_DATABASE.md).

```bash
nextflow run . -profile field \
    --reference_database ~/refdb/vespa18s_v0.4.0 \
    --assay assets/assay/vespa18s.json --primers vespa_primers.tsv \
    --samplesheet samplesheet.tsv --runs runs.tsv --outdir results
```

What changes:

- **Pre-flight** validates the release instead of a panel: every file
  present and unchanged since the build, its primers equal to
  `--primers`, its adapter anchors equal to the assay's. The release's
  look-alike mask becomes `preflight/anchor_lookalikes.tsv`, which the
  splitter and demultiplexer read exactly as they read a panel's.
  `--spiked`, `--mocks` and Illumina runs are refused in this mode, and a
  plate without negative controls is warned about: with no plasmids or
  standards, blanks are the only measure of contamination left.
- **Assignment** is per-sample consensus, not per read: reads are grouped
  by their best database hit, a consensus is built for each group, hidden
  second haplotypes are split off by linkage, and each consensus is
  named only down to the rank its identity supports. Output is in
  `assign_db/` (§5).
- **Synthesis and the report are database mode's own** (§5): relative
  18S abundance at the rank each call reached, judged against the run's
  blanks and `--db_min_rel_abundance`, with off-target fractions
  (bacteria and archaea, organelles, fungi, host) per well and the
  release version on every row. No copies, no plasmid checks.

Nanopore or Illumina; needs minimap2 (§1.1). On Illumina the pairs are
not merged - at 2x150 no V4 insert overlaps - but cut to a fixed window
of insert from each mate (`--db_window_r1`, `--db_window_r2`), denoised
the way DADA2 or UNOISE would, and each surviving pair is classified on
both windows together. A pair whose insert is shorter than a read is
kept whole.

---

## 3. Parameters

### 3.1 Inputs

| parameter | default | meaning |
|---|---|---|
| `--reference` | — | reference panel TSV (required) |
| `--samplesheet` | — | samplesheet CSV (required) |
| `--runs` | — | run manifest TSV; one of this or `--input_dir` |
| `--input_dir` | — | single-run shortcut; requires `--platform` |
| `--platform` | — | `illumina` or `nanopore`; only with `--input_dir` |
| `--run_id` | `RUN01` | run id when using `--input_dir` |
| `--primers` | — | PCR1 primers TSV; **required** |
| `--metadata` | — | grouping metadata TSV, keyed on `subject_id` **or** `sample_id` |
| `--mocks` | — | mock composition TSV |
| `--comparison` | — | a previous run's `calls.tsv`, to pair against this one |
| `--index_map` | `${projectDir}/assets/index_map.tsv` | index primer map |
| `--min_index_edit` | `null` | fail if a nanopore run's indices are closer than this in **edit** distance. Unset, the shortfall is only reported |
| `--assay` | `${projectDir}/assets/assay/ampseq_helminth.json` | assay definition |

### 3.2 Assay overrides

`null` means "defer to the assay config", which is not the same as an
empty value.

| parameter | default | meaning |
|---|---|---|
| `--barcode_repeats` | `null` | tandem index copies: 2 is the AMPSEQ design, 1 a conventional single index |
| `--tag_max_edits` | `null` | control-tag matching tolerance |
| `--expected_tag` | `null` | control tag sequence. An explicitly empty value means "skip the cross-check" |

On `barcode_repeats`: the difference between 1 and 2 is whether each
index carries its own replicate, and therefore whether a frame slip from
an indel can be detected at all.

### 3.3 Read preparation (Illumina)

| parameter | default | meaning |
|---|---|---|
| `--min_length` | `100` | minimum read length after trimming |
| `--min_overlap` | `20` | minimum overlap to merge a pair |
| `--max_ee` | `2.0` | maximum expected errors |
| `--illumina_read_length` | `250` | used by the panel geometry checks |

### 3.4 Concatemer splitting (nanopore)

| parameter | default | meaning |
|---|---|---|
| `--anchor_max_edits` | `2` | tolerance when locating adapter junctions |
| `--min_unit_length` | `100` | pieces shorter than this cannot carry an insert |

Ligation during library prep joins amplicons into single reads. Left
whole, such a read carries one molecule's i5 and another's i7 — a pair of
real indices belonging to a well that held neither. Cutting must happen
*before* demultiplexing, because once that pair has been looked up there
is no evidence left that anything was wrong.

On the real MinION run this affected 41.6% of reads.

### 3.5 Demultiplexing (nanopore)

| parameter | default | meaning |
|---|---|---|
| `--max_mismatch` | `1` | barcode matching tolerance |
| `--include_failed` | `false` | also use `fastq_fail` reads |

### 3.6 Assignment

| parameter | default | meaning |
|---|---|---|
| `--max_dist_frac` | `0.10` | maximum edit distance as a fraction of length |
| `--min_margin` | `2` | how far the best hit must beat the second |

### 3.6.1 Database mode

| parameter | default | meaning |
|---|---|---|
| `--reference_database` | — | a vespa-18s-refdb release directory, in place of `--reference` (§2.9) |
| `--db_min_cluster_reads` | `5` | reads a group needs for a consensus of its own; smaller groups' reads are still counted |
| `--db_lca_margin` | `0.005` | database hits within this fraction of the best alignment score take part in the LCA |
| `--db_floors` | `species=0.995,genus=0.99,family=0.95` | identity floor per rank; a consensus is named at a rank only when its best hit is this close. Below family it is reported as `unclassified, NN% to X` |
| `--db_min_rel_abundance` | `0.001` | a call must also be this share of its well's classified reads; `--min_depth` and `--absolute_floor` apply as in panel mode |
| `--db_primer_max_edits` | `3` | edits allowed when finding each primer in database mode (panel mode uses 2). Off-target bacterial 16S differs from the forward primer at its 3' end, often by three; at two it is dropped as `no_primer` rather than counted as bacteria |
| `--db_host_clade` | `Craniata` | a name anywhere in a best hit's PR2 lineage that marks host reads, for the off-target fractions (Craniata: vertebrates) |
| `--db_floors_illumina` | `species=0.997,genus=0.995,family=0.95` | the same for Illumina pairs, over both windows together: stricter, because the windows miss most of what separates sister taxa |
| `--db_window_r1` | `null` | Illumina database mode: insert bases kept from R1, so error copies line up base for base; `null` takes the release's `illumina_windows` (132 for 2x150) |
| `--db_window_r2` | `null` | the same for R2 (128 for 2x150); set both for longer reads, e.g. 250 and 200 for 2x300 |

### 3.7 Control plasmids

| parameter | default | meaning |
|---|---|---|
| `--plasmid_warn` | `0.05` | plasmid share above which a non-mock well is flagged |

The tagged control plasmids are handled **per well, by `sample_type`**.

- **In a `mock` well they are the content.** A mock is assembled from
  the plasmids, so it is processed like any sample and its plasmids are
  called, counted and - on a `--spiked` run - quantified, exactly as
  targets are. Its composition chart shows them. Its near-total plasmid
  share is expected, so it is never flagged as contaminated.
- **In every other well they are contamination**, recognised by the
  tag. They are never called as a finding. Each well's `plasmid_share`
  (plasmid reads over all panel amplicon reads) is in `calls.tsv`, and on
  a spiked run so is `plasmid_copies` - the contaminating plasmid in
  copies, on the same scale as the targets. Wells above `--plasmid_warn`
  are listed in run QC and the report, which also give the run-wide
  share across non-mock wells.

A plasmid titration, or any plate whose subject is the plasmids, types
those wells as `mock`. The run-wide switch this replaces could only get
one of the two right on a plate holding both mocks and field samples: off, every mock came out empty and flagged;
on, carryover in a field sample was called as a detection.

**What this means for counts built from `calls.tsv`.** It holds one row
per (sample, *called* reference): targets in every well, plus the control
plasmids in mock wells. `group_summary.tsv` counts everything in it,
with no filter on `reference_type` - so a group containing mocks pools
plasmid calls (the mock's content) with target calls. Split by
`reference_type` or `sample_type` in `calls.tsv` when that matters.

**Flags on the command line.** From Nextflow 26.04 every command-line
parameter arrives as a string, so `--include_failed false` or
`--allow_missing false` is the *string* `"false"` - truthy in Groovy,
and it switches the option **on**. Give a flag bare to turn it on and
leave it off to turn it off; to set one explicitly, use
`nextflow.config` or a `-params-file`, which preserve types. `--spiked`
is compared as text, so `--spiked false` does mean unspiked and
`--spiked yes` is refused.

### 3.8 Calling and run QC

| parameter | default | meaning |
|---|---|---|
| `--spiked` | `false` | a flag: the quantification standards were added to every sample. Off, the run is unspiked |
| `--min_depth` | `100` | below this a well is `assay_failed`, not negative |
| `--absolute_floor` | `2` | a detection needs **more** reads than this, however clean the blanks |
| `--permutations` | `2000` | permutations for the calibration test |
| `--slope_tolerance` | `0.15` | deviation from a slope of 1 that earns a well a note; outside 0.5–2 copy numbers are withheld regardless |
| `--min_composition_reads` | `100` | minimum reads to report composition |
| `--replicate_suffix` | `null` | regex stripped from `sample_id` to pair replicates across runs, e.g. `'-r[0-9]+$'` |

**`--spiked` decides what the run can report.** Without it - the
default - the run is unspiked: calls are presence/absence, no copy
numbers are given, every negative is unverified, and any standard reads
that turn up are reported as carryover rather than fitted. With it,
every well is calibrated against its standards, copy numbers are given,
and a well missing its standards has failed (`assay_failed`), however
many wells that is. The report follows suit: an unspiked report states
once that the run was not spiked and leaves out quantification.

Whether the standards went in is declared, not inferred. Inferring it -
spiked if a quarter of wells recovered a standard - fails both ways: a
run whose spikes mostly failed reads as unspiked, so real failures pass
as unverified negatives, and carryover across a quarter of the plate
reads as a spike. That rule survives only as a cross-check: when the
reads disagree with the declaration, run QC and the report say so.
`--spiked` with no `quant` entries in the reference is refused at
pre-flight. The shell drivers take `SPIKED=1` for the same thing.

`--min_depth` is worth checking against your actual depths. A fixture or
a shallow run with wells below it will come back entirely
`assay_failed`, which is correct behaviour and looks alarming.

### 3.9 Output

| parameter | default | meaning |
|---|---|---|
| `--outdir` | `./results` | results directory |
| `--publish_dir_mode` | `copy` | `copy`, `symlink`, `link`, … |

Note that `copy` does **not** remove files from a previous run. Re-reading
a stale `--outdir` after a failed run is a good way to measure last
week's results; delete it, or use a fresh one.

---

## 4. Profiles

| profile | what it does |
|---|---|
| *(none)* | **LSF**. This is the default |
| `local` | run on this machine; explicit opt-out from LSF |
| `stub` | placeholder tasks, for proving the DAG connects |
| `test` | Illumina fixtures, end to end |
| `test_ont` | nanopore fixtures, end to end, with a real 4-way split fan-out |
| `test_db` | database mode on nanopore: a miniature release, held-out taxa, a hidden haplotype (`python3 tests/make_db_fixture.py` first) |
| `test_db_il` | database mode on Illumina 2x150: the same release and samples as `test_db`, as read pairs |
| `farm` | no-op, kept so existing commands still work |
| `singularity` / `docker` / `conda` | enable that engine |
| `debug` | echo hostnames, keep work directories |

Combine with commas: `-profile test_ont,local`.

Every tool here is standard-library Python 3 with no third-party
imports, so nothing needs installing and no container is required. The
modules still declare one, used only if an engine is enabled — a
reproducibility choice, not a dependency.

---

## 5. Outputs

```
<outdir>/
  preflight/
    validation.log              every check, passed and failed
    file_resolution.tsv         each sample resolved to its files
    sample_metadata.tsv         the grouping join, one row per library
    runs_resolved.tsv           the run manifest with absolute paths
    reference.fasta             the panel
    reference.insert.fasta      the panel, primers trimmed
    anchor_lookalikes.tsv       adapter-anchor near-copies inside panel
                                sequences; the nanopore splitter and
                                chimera check leave these sites alone
    tag_report.tsv / .txt       control tag derivation
    preflight_summary.json
  split/                        nanopore only
    reads/                      amplicon units, one file per input file
    split_summary.json          read funnel starts here
    parts/                      per-shard statistics
    file_report.json            which input files were used and excluded
  demux/                        nanopore only
    reads/                      one FASTQ per sample
    assignments.tsv             per-read outcome
    demux_summary.json          the tier breakdown
    index_matrix.tsv            observed i7 x i5 grid, unclaimed pairs marked
  prep/                         Illumina only
    <sample>.prep.json
  assign/                       panel mode
    <sample>.counts.tsv         reads per reference entry
    <sample>.taxa.tsv
    <sample>.assign.json
  assign_db/                    database mode (--reference_database)
    <sample>.clusters.tsv       one row per consensus: reads, rank reached,
                                taxon, identity, best hit, consensus
    <sample>.taxa.tsv           reads per taxon at the rank reached
    <sample>.consensus.fasta
    <sample>.assign_db.json     read accounting, settings, release version
  synthesis/                    in database mode:
    db_calls.tsv                one row per sample x taxon at the rank
                                reached: share, call, limit, identity,
                                organism group, share within the group
    db_wells.tsv                per well: read funnel, reads and fractions
                                per organism group, target-group richness,
                                reads by rank
    db_prevalence.tsv           per taxon: samples detected of samples tested
    db_run_qc.json              release, blanks, internal adapter rates,
                                the run's read funnel
  synthesis/
    calls.tsv                   the results table
    calibration.tsv             per-well internal control
    read_fates.tsv              the read funnel, + per-well fit distance
    group_summary.tsv           the run consolidated by each grouping field
    run_qc.json                 run-level QC
  calibration/                  only with --mocks
  comparison/                   only with --comparison
  <run_id>_report.html          self-contained
  pipeline_info/                timeline, report, trace
```

### Reading `assign_db/<sample>.clusters.tsv` (database mode)

One row per consensus sequence. `rank` is the deepest rank the call
reached: the lowest common ancestor of every database hit within
`--db_lca_margin` of the best, cut back to the deepest rank whose
identity floor (`--db_floors`) the best hit clears. `identity` is that
best hit's, counting unaligned ends as differences - read it beside the
rank: a species call at 0.996 and one at 1.000 are not the same claim,
and a genus call at 0.993 says a near relative is in the database but
this sequence is not. Below family a row is labelled `unclassified, NN%
to X`, which is what a gap in the database looks like.

| column | meaning |
|---|---|
| `reads`, `rel_abundance` | reads counted to this consensus, and their share of the sample's clustered reads. **Relative 18S abundance**, never organisms or copies: 18S copy number varies widely between taxa |
| `rank`, `taxon`, `label` | the rank reached, the name there, and the name or `unclassified, NN% to X` |
| `lineage` | the best hit's names down to the rank reached |
| `identity`, `coverage` | the best hit's identity and the share of the consensus it aligns |
| `best_hit`, `best_hit_lineage` | the closest record and its full lineage, even below the rank reached |
| `near_hits`, `lca_rank` | how many records shared the LCA, and where the LCA alone stopped before the floors |
| `db_resolution` | the deepest rank at which the best hit is distinguishable over the whole insert, from the release's `resolution_nanopore.tsv` |
| `origin` | `group`, or `split…` when the consensus was a hidden haplotype split from its group by linkage |
| `database_version` | the release every call was made against |

`<sample>.assign_db.json` holds the read accounting: `reads_in` =
`no_primer` + `too_short` + `clustered` + `unclustered`. Unclustered reads
belong to taxa with too few reads here for a consensus
(`--db_min_cluster_reads`) and are not classified - one read is not
enough to say what it is.

### Reading `synthesis/db_calls.tsv` (database mode)

One row per sample and taxon at the rank reached. `call` is `detected`
(beats the blanks and is at least `--db_min_rel_abundance` of the well),
`below_blanks`, `below_floor`, or `low_depth` (the well has fewer than
`--min_depth` classified reads). `limit_reads` is the fewest reads that
beat the blanks and `limit_basis` says how it was set; with no blanks on
the plate it is the fixed floor, and says so. `group` is the organism
group, as VESPA reports them - helminths, protists and microsporidia (the
targets), host, fungi, plants, invertebrates, bacteria_archaea, organelle,
or unassigned - and `group_share` the call's share within its group, so
a helminth composition is not flattened by a well that is mostly host.
The groups are defined in docs/REFERENCE_DATABASE.md §11.4. And
`primer_fit` carries the release's prediction for its best hit -
`near_3prime` means its share may be held down by the PCR.
`db_wells.tsv` gives each well's funnel, reads and fractions per group,
target-group richness and reads by rank reached; `db_prevalence.tsv`
how many samples each taxon was detected in.

### Reading `calls.tsv`

Calls are three-state, not two:

| call | meaning |
|---|---|
| `detected` | above the detection limit |
| `not_detected` | below it, **and the assay demonstrably worked** |
| `assay_failed` | the well cannot support either claim |

The distinction matters: a negative result is only a negative result if
the internal control says the well amplified. `assay_failed` is not a
missing value to be dropped — it is the honest answer for a well whose
spike-in did not come back.

Copy numbers are a second, separate claim. They need a usable standard
curve, so the `internal_control` column says whether the well had one:

| `internal_control` | meaning |
|---|---|
| `spiked` | both standards recovered, slope plausible — copies reported |
| `slope_implausible` | both standards recovered, but the fitted slope is outside 0.5–2 (or negative). Calls stand; copies are withheld, because inverting that curve would invent them |
| `unspiked` | no standard present — calls stand, but a negative is unverified |
| `spike_failed` | standards present but wrong — the well is `assay_failed` |

### Reading the fit distance

`read_fates.tsv` carries four columns per well; `run_qc.json` carries
the pooled histogram under `assignment_distance`; the report draws it
under **How well the reads fit**.

| column | meaning |
|---|---|
| `assigned_dist_median` | the error mode — the platform, not the assay |
| `assigned_dist_max` | furthest assigned read; one read, so noisy |
| `assigned_dist_over_trough` | assigned reads beyond the run's trough |
| `assigned_dist_over_trough_frac` | the same as a share of the well |

**The trough is found, not chosen.** Synthesis walks the run's own
pooled histogram out from the error mode; if it dips and then rises
again by at least 2x, the dip is the boundary. A monotone histogram has
no trough, so `trough` is null, the two columns are empty for the whole
plate, and that is the clean answer — no threshold had to be invented
to produce it.

A 95th percentile shipped here first and was useless. The population it
existed to find was **0.70%** of assigned reads on a plate spiked with
1.5% chimeras, so it sat above p95 *and* p99, and the column read `2`
in every well, spiked and clean alike. A quantile cannot resolve a tail
smaller than its own complement. The trough share separates the same
wells 1.38% against 0.00%.

These describe reads that **were assigned** — it is not a loss
category. Read the shape, not the median:

- **one mode at the platform's error rate** — a clean run
- **a second mode further out** — reads that are not errors

The second case is what a read built from two panel members looks like.
It sits at roughly half the divergence between its two parents, is
assigned to one of them, and **passes the margin rule** — because the
read genuinely is closer to one reference, the substitutions separating
two references never being evenly spread. Measured on a two-reference
panel 5% divergent, a 50/50 chimera came back `assigned` at distance 7
with a margin of 4, against a `min_margin` of 2. So it is counted as a
confident species call and appears in no outcome tally anywhere. This
is the only place it is visible.

No threshold is applied and nothing is filtered. Where the second mode
sits depends on which references recombined and on the platform's error
rate, so a cutoff set here would be invented rather than measured.
Sorting `read_fates.tsv` on `assigned_dist_over_trough_frac` finds the
odd well; the pooled histogram finds the artifact. Nothing is filtered
on either.

### Reading `group_summary.tsv`

One row per (`field`, `value`): the run consolidated by each grouping
field in turn. Present whenever there is a grouping to write — a
`--metadata` file, or simply a subject with more than one library.

Three things about it decide how to read it:

**Groups are never pooled.** Where the groups are the arms of a
designed comparison, summing them answers a question nobody asked and
destroys the one that was. Each arm gets its own row, side by side.

**Controls are excluded from every metric**: negatives are counted in
`n_controls` and positives in `n_positives`. A median depth with an NTC
in the denominator moves when the number of NTCs does, and a
deliberately positive well inflates a detection rate - either would make
two runs incomparable for a reason invisible in the table. What the negatives say gets its own
column instead: `deepest_negative`, the read count of the deepest
negative in that group — which is how you tell whether an arm stayed
clean.

**Every library appears.** One with no metadata row sits under
`(no value)`. So `n_libraries` summed over a field's rows equals the
plate, and a small group is distinguishable from a broken join.

| column | meaning |
|---|---|
| `field`, `value` | the grouping and one of its values |
| `n_libraries` | all libraries in the group, controls included |
| `n_controls` | negatives, excluded from everything below |
| `n_positives` | positive controls, likewise excluded |
| `n_analysed` | the denominator for every metric that follows |
| `n_subjects`, `n_empty` | distinct subjects; libraries with zero reads |
| `reads_assigned` | reads reaching assignment — "how many reads per group" |
| `median_depth`, `min_depth`, `max_depth` | the spread, not just the centre |
| `reads_<fate>` | the eight read fates, summed over the group |
| `target_fraction` | on-target share — the usable fraction |
| `plasmid_share` | median plasmid carryover |
| `detected`, `not_detected`, `assay_failed` | target calls |
| `detection_rate` | detected / (detected + not detected) |
| `wells_<state>` | internal-control states |
| `deepest_negative` | the deepest negative control in this group |
| `note` | why a field was skipped, when one was |

There is **no test statistic here on purpose**. With two arms and a
handful of libraries each, a p-value would be read as an answer when
it is barely evidence. A field with more distinct values than
libraries-worth of meaning — above 50 — is skipped with a note, on the
grounds that it is an identifier rather than a category.

---

## 6. Checking a run finished honestly

Green output is not the same as a measured result. This project has had
three false passes with exactly that shape, so check the numbers:

```bash
# every declared sample reached the results?
cut -f2 results/synthesis/calls.tsv | tail -n +2 | sort -u | wc -l

python3 -c "import json;q=json.load(open('results/synthesis/run_qc.json'));
print('complete:', q['cohort_complete'], 'samples:', q['n_samples'],
      'absent:', q['samples_absent'])"

# did the split fan out, and does the funnel start where it should?
python3 -c "import json;s=json.load(open('results/split/split_summary.json'));
print('parts:', s['parts'], 'reads_in:', s['reads_in'],
      'pieces/read:', round(s['pieces_per_read'],3),
      'bases retained:', round(s['bases_retained_fraction'],4))"
```

`cohort_complete: false` with a non-empty `samples_absent` means wells
declared in the samplesheet did not reach synthesis. Synthesis refuses
outright when those wells had reads — reporting on part of a plate is the
one failure mode that produces a confident wrong answer rather than an
error.

`bases_retained_fraction` should be ~1.0. Splitting removes nothing but
pieces below the length floor, so if bases go missing beyond that, the
cut arithmetic is wrong and the output is no longer a rearrangement of
the input.

---

## 7. Troubleshooting

**`no sample reached assignment`** — both platform branches are empty.
The `platform` column disagrees between the samplesheet and the run
manifest, or names something other than `illumina`/`nanopore`.

**`run '<id>': reads path does not exist`** — the message prints what the
manifest said and what it resolved to. A relative path resolves against
the manifest, not your working directory.

**`cohort incomplete: N samples declared, M reached synthesis`** — tasks
were lost. On Nextflow this is usually a fan-out capped by a single-item
queue channel; `tests/test_nextflow.py::test_fanout_is_not_capped` checks
for that shape statically. `--allow-missing` overrides, and the QC then
records that the plate was partial.

**`expected N shard(s), found M`** — a split shard did not arrive.
Refusing here is deliberate: a missing part produces a read funnel that
balances perfectly against a smaller run than the one sequenced.

**A warning about index edit distance.** On nanopore, pre-flight
reports when the indices a run uses are closer than 3 edits apart,
because below that no indel can be corrected and indels are the
dominant error. It is a warning, not a failure: an existing plate
cannot be changed, and refusing to analyse reads already sequenced
helps nobody. `--min_index_edit 3` makes it a hard requirement, which
is what you want when designing a *new* plate.

`bin/design_indices.py` finds the largest subset of a primer pool
meeting a distance floor. Often the fix is dropping one primer rather
than ordering new ones — in the reference pool a single pair
(`N721 TACGCTGC` vs `N723 TAGCGCTC`, 2 edits) is the only offender, and
dropping either leaves 16 x 23 = 368 usable wells at edit distance 3.

**Most reads are tier B (one end only)** — read
`demux/demux_summary.json`, fields `i5_end_reasons` and `i7_end_reasons`.
They count per END, and each map sums to `total_reads`. The two causes
look identical in the tiers and have opposite remedies:

| dominant reason | what it means | what fixes it |
|---|---|---|
| `anchor_not_found` | the molecule never carried that index — PCR2 did not add it | the bench, not the software. No tolerance recovers a missing index |
| `no_copy_matched` | the index IS there but matched nothing declared | check the index set, or `--max_mismatch` |
| `copies_disagreed` | the tandem copies read as different indices | the frame slip the duplication exists to catch |
| `barcode_off_read` | anchor found, barcode runs past the read end | truncated reads |
| `ok_single_copy` | read from one copy; the replicate could not confirm it | usually benign |

On the reference MinION run tier B was 51.7%, and the same 52.3% was
present *before* concatemer splitting — so splitting did not cause it.

**Everything is `assay_failed`** — check `--min_depth` against your
actual depths in `calls.tsv`.

**"I passed that flag and it did nothing."** Read what the task actually
ran, rather than what you typed:

```bash
grep -l -- "--spiked" work/*/*/.command.sh
```

`.command.sh` is the rendered script for one task, so it settles whether
a parameter reached the tool — upstream of any assumption about Nextflow
parsing, your config, or the module.

**Make sure it is the right task's.** A failed attempt leaves its work
directory behind looking entirely authoritative, and the hash printed in
an error message is the *failed* one. After a fix and a `-resume` the
successful task has a different hash, and reading the old directory will
tell you the opposite of the truth. `grep -l` across `work/*/*/` finds
the task that has the flag rather than the one you happen to know the
path of; `nextflow log` lists them properly.

This is worth reaching for early. Twice in one session it resolved in
seconds a question that inference had got wrong in both directions.

**A `-stub-run` reports `Complete` having run almost nothing** — read the
task table, not the last line. Every process should show a count.

---

## 8. Re-analysing a run you already demultiplexed

`--stats-only` recomputes the tiers, the per-end reasons and the index
matrix **without** writing per-sample FASTQ:

```bash
python3 bin/demux_ont.py --reads results/split/reads \
    --samplesheet samplesheet.tsv --run-id RUN01 \
    --assay assets/assay/ampseq_helminth.json --threads 16 \
    --stats-only --outdir demux-reasons
```

Every count is identical to a normal run — tiers, end reasons,
`per_sample`, `assigned_reads`. What you save is ~8 GB of disk and a
second copy of reads you already have.

**It saves little time.** Measured on 51k reads: 15.2s writing, 16.3s
not. The write is about 2% of this stage and the per-read scanning is
the rest. The flag was added on the assumption that writing dominated,
and the measurement disproved it; the honest reason to use it is disk
and not clobbering existing output.

`stats_only` is part of the per-part fingerprint, so a stats-only run
cannot be resumed onto by a real one. Without that, a real run would
find the parts marked finished, find no FASTQ to merge, and report
success on an empty output. `reads_written_balanced` comes back `null`
rather than `true`, because a check that did not run is not a check
that passed.

## 9. The shell drivers

```bash
bin/run_illumina.sh  results reference.tsv samplesheet.tsv runs.tsv primers.tsv
THREADS=16 bin/run_nanopore.sh results reference.tsv samplesheet.tsv runs.tsv primers.tsv
```

Predate the Nextflow layer and are kept because every stage in them has
been executed on real data. They are single-platform, take positional
arguments, and reimplement in shell what `-resume` does properly. The
nanopore driver takes `THREADS`, `SPIKED`, `PLASMID_WARN`,
`MIN_DEPTH`, `PERMUTATIONS`, `METADATA` and `RUN_ID` from the
environment.

Use them until the workflow has been shaken down on real data, then
stop. The longer both exist the further they drift, and they have
drifted once already.
