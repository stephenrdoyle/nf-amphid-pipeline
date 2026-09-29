# nf-helminth-amplicon — usage

Targeted amplicon panel analysis for Illumina and nanopore. One
invocation handles either platform, or both at once: the `platform`
column of your inputs decides which branch each run takes.

Every column list and default below was read out of `bin/preflight.py`
and `nextflow.config` rather than written from memory. If the two ever
disagree, the code is right and this file is stale.

---

## 1. Quick start

```bash
NF=/nfs/users/nfs_s/sd21/lustre_link/software/nf-helminth-amplicon

nextflow run $NF -resume -bg \
    --reference   reference.tsv \
    --samplesheet samplesheet.csv \
    --runs        runs.tsv \
    --primers     primers.tsv \
    --outdir      results
```

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

### 2.2 `--samplesheet` — one row per library (CSV, **required**)

A samplesheet row is a *library*, and it is frozen at prep time.

| column | required | meaning |
|---|---|---|
| `run_id` | yes | must match a `run_id` in the run manifest |
| `platform` | yes | `illumina` or `nanopore` |
| `sample_id` | yes | unique within the run |
| `index_i7` | yes | i7 index sequence or name |
| `index_i5` | yes | i5 index sequence or name |
| `sample_type` | yes | `unknown`, `neg_extraction`, `neg_ntc`, `mock`, `pos_control` |
| `subject_id` | no | joins to the metadata file; one subject may be several libraries |
| `file_prefix` | no | overrides filename matching for Illumina file resolution |

The `sample_type` vocabulary is **not** configurable. The logic depends
on it: blanks set detection limits, NTCs measure contamination, mocks
calibrate efficiency, and positives confirm the assay ran.

Leaving a few index combinations unused is worth doing. If every
combination in the pool is assigned, no sentinel exists and index
misassignment becomes *unmeasurable* — which is a different claim from
measuring it to be zero. Pre-flight warns when this happens.

### 2.3 `--runs` — one row per sequencing run (TSV, **required**)

Where the reads are. On nanopore the reads belong to the *run*, not to
any one sample, which is why this file exists separately.

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

### 2.4 `--primers` — PCR1 primers (TSV, optional but effectively required)

| column | required | meaning |
|---|---|---|
| `name` | yes | primer name |
| `direction` | yes | `fwd` or `rev` |
| `sequence` | yes | primer sequence, IUPAC allowed |

Used to trim primer sequence off the panel, producing
`preflight/reference.insert.fasta`. The nanopore path **requires** this:
without it there is no primer-trimmed panel and the run stops.

### 2.5 `--metadata` — subject metadata (TSV, optional)

Kept deliberately *out* of the samplesheet. A samplesheet row is a
library and is frozen at prep; a metadata row is a subject and keeps
changing; and one subject is usually several libraries. Pre-flight
validates the join rather than assuming it.

| column | required | meaning |
|---|---|---|
| `subject_id` | yes | joins to `subject_id` in the samplesheet |
| anything else | no | your own fields |

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
| `--primers` | — | PCR1 primers TSV |
| `--metadata` | — | subject metadata TSV |
| `--mocks` | — | mock composition TSV |
| `--comparison` | — | a previous run's `calls.tsv`, to pair against this one |
| `--index_map` | `${projectDir}/assets/index_map.tsv` | index primer map |
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

### 3.7 Control plasmids

| parameter | default | meaning |
|---|---|---|
| `--controls_as_targets` | `false` | give control plasmids detection calls too |
| `--plasmid_warn` | `0.05` | plasmid share above which a well is flagged |

By default control plasmids are contamination and their share of each
well is reported. Set `--controls_as_targets` on a run whose *subject* is
the plasmids — a titration, or a mock-only pool. **The proportions are
emitted either way**; the flag changes calling only, because a plasmid in
a titration is the measurement and the same plasmid in a faecal sample is
contamination.

### 3.8 Calling and run QC

| parameter | default | meaning |
|---|---|---|
| `--min_depth` | `100` | below this a well is `assay_failed`, not negative |
| `--absolute_floor` | `2` | reads below this are never a detection |
| `--permutations` | `2000` | permutations for the calibration test |
| `--slope_tolerance` | `0.15` | tolerated deviation in the standard curve slope |
| `--min_composition_reads` | `100` | minimum reads to report composition |
| `--replicate_suffix` | `null` | regex stripped from `sample_id` to pair replicates across runs, e.g. `'-r[0-9]+$'` |

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
    runs_resolved.tsv           the run manifest with absolute paths
    reference.fasta             the panel
    reference.insert.fasta      the panel, primers trimmed
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
  assign/
    <sample>.counts.tsv         reads per reference entry
    <sample>.taxa.tsv
    <sample>.assign.json
  synthesis/
    calls.tsv                   the results table
    calibration.tsv             per-well internal control
    read_fates.tsv              the read funnel
    run_qc.json                 run-level QC
  calibration/                  only with --mocks
  comparison/                   only with --comparison
  <run_id>_report.html          self-contained
  pipeline_info/                timeline, report, trace
```

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

**A `-stub-run` reports `Complete` having run almost nothing** — read the
task table, not the last line. Every process should show a count.

---

## 8. Re-analysing a run you already demultiplexed

`--stats-only` recomputes the tiers, the per-end reasons and the index
matrix **without** writing per-sample FASTQ:

```bash
python3 bin/demux_ont.py --reads results/split/reads \
    --samplesheet samplesheet.csv --run-id RUN01 \
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
bin/run_illumina.sh  results reference.tsv samplesheet.csv runs.tsv primers.tsv
THREADS=16 bin/run_nanopore.sh results reference.tsv samplesheet.csv runs.tsv primers.tsv
```

Predate the Nextflow layer and are kept because every stage in them has
been executed on real data. They are single-platform, take positional
arguments, and reimplement in shell what `-resume` does properly. The
nanopore driver takes `THREADS`, `CONTROLS_AS_TARGETS`, `PLASMID_WARN`,
`MIN_DEPTH`, `PERMUTATIONS` and `RUN_ID` from the environment.

Use them until the workflow has been shaken down on real data, then
stop. The longer both exist the further they drift, and they have
drifted once already.
