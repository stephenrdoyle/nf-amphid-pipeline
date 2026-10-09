# nf-helminth-amplicon

Targeted amplicon panel analysis for Illumina and nanopore, built around
competitive mapping against a defined reference panel with spiked
quantification standards and tagged control plasmids.

**Status: all stages implemented and tested; the Nextflow layer has not
yet been run on real data.**

Every analysis stage — pre-flight, concatemer splitting, demultiplexing,
read preparation, assignment, synthesis, reporting, mock calibration and
cross-run comparison — exists, is covered by tests against a simulator
with a truth table, and has been executed on real nanopore data through
the shell drivers in `bin/`.

**Both branches now run end to end under Nextflow.** `-profile
test,local` takes the Illumina fixtures through, 8 samples in and 8
called; `-profile test_ont,local` does the same for nanopore, with
splitting fanned out to one task per MinKNOW batch file.

The nanopore run has been checked against the shell driver on the same
fixture, and the two agree exactly: identical per-sample depths,
identical per-target read counts, and identical calls for all 8
samples across 5 targets. That is the sharded Nextflow route and the
threaded single-task route reaching the same answer through different
orchestration.

Getting there needed one correction worth recording. The first
comparison showed every call differing, which looked like the two
routes disagreeing. They did not: `-profile test_ont` sets
`min_depth = 10` while the driver hardcoded the default of 100, and
the fixture's wells are 39-96 reads deep, so one route called
everything `assay_failed` and the other called it properly. The reads
and depths had matched all along. The driver now takes `MIN_DEPTH` and
`PERMUTATIONS` from the environment so the two can be compared like
for like - a driver that hardcodes what the workflow exposes cannot be
used to check the workflow.

The claim that used to stand here — that a `-stub-run` proved the DAG
connects — was wrong, and worth recording because of how it was wrong.
The stub `file_resolution.tsv` carried one row saying `nanopore`, while
the test manifest said `illumina`. So the Illumina branch saw no
matching row, the nanopore branch saw no matching run, every process
after pre-flight was skipped, and Nextflow printed `Complete`. One
process had executed, and the green output read as a passing test of
the whole workflow.

Correcting that row to `illumina` fixed the symptom and left the cause:
**a stub does not read its inputs**, so a hardcoded `file_resolution`
*decides* the platform instead of reflecting it, and
`-profile test_ont -stub-run` still ran the Illumina branch with an
invented sample called `STUB_SAMPLE`. A workflow whose routing IS its
input data cannot be stub-tested by a stub that invents that data.

So `PREFLIGHT` and `SPLIT_LIST` now run **for real** under `-stub-run`.
Both are standard-library Python that touch no reads and finish in
about a second — between them they read a samplesheet and list a
directory. They are the two processes whose real execution is cheaper
than a convincing fake, and faking them is what made the stub run
meaningless: the first invents the platform column that chooses the
branch, the second invents the file list whose length IS the fan-out
under test. Every other process still stubs normally.

A related lesson, from the fix that followed. The nanopore branch of
main.nf used to re-read the run manifest itself and resolve relative
paths in Groovy - duplicating what preflight.py already did. That
duplicate was wrong twice in succession: first it resolved against the
launch directory, then the corrected version used an `isAbsolute()`
guard that could never fire, because Nextflow's `file()` has already
made the path absolute by the time you ask. Both shipped, because
Groovy cannot be executed in the environment this is developed in.

The third attempt deleted the duplicate instead of fixing it.
Pre-flight now writes `runs_resolved.tsv` with absolute paths and the
workflow reads that, which restores the rule main.nf already stated
and the nanopore branch was quietly breaking: *a path pre-flight
blessed is the only path the rest of the pipeline sees*. The
resolution happens once, in Python, where it has tests. If logic
cannot be executed, the answer is to move it somewhere that can be
rather than to write it more carefully.

Three of this project's false passes have now had the same shape: a
comparison of two empty files that reported "identical", a sample count
re-read from a stale results directory, and a stub run that executed
one task. In each case the output looked like success and measured
nothing. The guards added since — non-vacuity assertions in the tests,
a cohort completeness check in synthesis, and `.ifEmpty { error ... }`
on the assignment channel so a run where no sample is processed fails
instead of completing — exist because reading the green line is not
the same as checking the number.

Nextflow cannot be installed in the environment this was developed in,
so the workflow is also checked statically — flags against each tool's
real argparse, `.out` references against emits, params against the
config, declared outputs against what the tools write, call arity
against each process's input block, distinct placeholders for optional
inputs, stub outputs that are parsed being structurally valid, and every
tool in `bin/` being executable with a shebang.

Six of those checks exist because the corresponding bug got through. The
first execution of the workflow found, in order: two optional inputs
colliding on the name `NO_FILE`; a CPU request larger than the host; a
memory request that `withLabel` precedence silently restored after a
profile had lowered it; an empty stub file that `splitCsv` could not
parse; a relative manifest path resolving into Nextflow's staging
directory; and nine tools without the executable bit, which only matters
because Nextflow invokes them as commands rather than as arguments to
`python3`. None was findable without running it, and none was about the
science.

## Why pre-flight comes first

All input validation, reference derivation and geometry inference happen
before a single read is touched. A samplesheet typo should cost seconds, not
a queue slot and an afternoon. It also exercises the input contract before
anything expensive depends on it.

## Quick start

```bash
# generate test fixtures (real index sets, real control tag)
python3 tests/make_fixtures.py --outdir tests/data

# run pre-flight directly - no Nextflow needed
python3 bin/preflight.py \
    --reference tests/data/reference.tsv \
    --samplesheet tests/data/samplesheet.tsv \
    --runs tests/data/runs.tsv \
    --index-map assets/index_map.tsv \
    --outdir preflight

# nanopore fixtures - several MinKNOW batch files, real concatemers
python3 tests/make_ont_fixture.py

# or through Nextflow. The two profiles exercise the two branches of
# main.nf; a failure should say which.
nextflow run . -profile test,local          # Illumina branch
nextflow run . -profile test_ont,local      # nanopore branch
```

All logic lives in `bin/preflight.py` as standard-library Python with no
imports beyond it, so it can be run and tested without Nextflow. The Nextflow
layer is orchestration only.

## Running a real run

Two ways, and they are not equals.

**Nextflow is the intended interface.** One invocation handles a run of
either platform, or both at once: `main.nf` reads the `platform` column
of the samplesheet, routes Illumina libraries and nanopore runs down
their own subworkflows, and merges them back into one assignment step.
That is the design, and the switch is built.

```bash
nextflow run . -profile test,local          # fixtures, end to end
nextflow run . -profile farm \
    --reference reference.tsv --samplesheet samplesheet.tsv \
    --runs runs.tsv --primers primers.tsv --outdir results
```

**The shell drivers are a stopgap** from before the Nextflow layer
existed, kept because every stage in them has been executed on real
data and the workflow has not. They are single-platform, take positional
arguments, and reimplement in shell what `-resume` does properly.

```bash
bin/run_illumina.sh  results reference.tsv samplesheet.tsv runs.tsv primers.tsv
THREADS=16 bin/run_nanopore.sh results reference.tsv samplesheet.tsv runs.tsv primers.tsv
```

Use the drivers until the workflow has been shaken down on real data,
then stop: the longer both exist the further they drift, and they have
drifted once already.

### Choosing the driver

The `platform` column has to say the same thing in **both** the run
manifest and the samplesheet, and it decides which driver applies. The two
paths are not interchangeable:

| | Illumina | nanopore |
|---|---|---|
| `path` in `runs.tsv` | directory of per-sample FASTQ | MinKNOW run folder, or a pooled FASTQ |
| demultiplexing | done by the sequencing centre | done here, from in-line indices |
| `files` in `file_resolution.tsv` | two paths per sample | **empty, by design** |

That last row is worth knowing before it surprises you. On the nanopore
path no read belongs to a sample until the barcodes have been read, so
pre-flight resolves nothing per sample and the `files` column is blank.
An empty `files` column with `platform = illumina` is a different thing
entirely — that is an error, and pre-flight raises it as one.

### Nanopore resource use

Demultiplexing is one serial pass over the pool, with **peak memory flat
in the read count** at about 30 MB. Throughput is ~3,400 reads/s on one
core with the junction check off, and depends on the tier-A fraction
with it on, since only reads that resolved both ends pay for the
interior scan: ~1,100 reads/s on a clean run where nearly everything is
tier A, ~2,200 on a run where a quarter is. Measured on 121,560
simulated reads of mean length 494bp: 38.8s off, 112.9s on, 28 MB either
way.

| reads | time (clean run) | peak RAM |
|---|---|---|
| 1M | ~15 min | ~30 MB |
| 10M | ~2h 30m | ~30 MB |
| 20M | ~5h | ~30 MB |

Those are the pessimistic figures — a run with a low tier-A fraction is
faster, because fewer reads reach the check. It is once per run, not per
sample.

The flatness is deliberate and was not free. The first version held every
read in memory until the end — about 1.4 KB per read, so a 20M-read flow
cell wanted ~28 GB and would have been OOM-killed an hour into the job.
Per-sample FASTQ are now written in buffered chunks as reads stream past
(`--flush-records`, default 2000), which means each file is a
multi-member gzip. That is valid and every reader handles it; the only
visible cost is about 0.2% in compressed size.

Two invariants are checked before the output is declared usable: every
read lands in exactly one of the four tiers and the tiers sum to the
input, and the reads written to per-sample files equal the reads called
to a sample. Either failing is a non-zero exit, not a warning.

### When demultiplexing goes badly

The tiers say how much evidence survived, not what destroyed it:

| tier | meaning | assigned? |
|---|---|---|
| A | both ends read, and the i5×i7 pair is one the samplesheet claims | yes |
| B | one end only — the other truncated or unreadable | only with `--accept-single-end` |
| C | the ends disagree, or the pair belongs to no sample | never |
| D | no confident match | never |
| E | an adapter sits *between* the two ends — the read is more than one molecule | never |

Tier E exists because of a failure mode that costs correctness rather
than yield. Two amplicons ligated head to tail present the first one's
i5 at the 5' end and the **second one's** i7 at the 3'. Both are real,
declared indices, so the pair is very likely a well on the plate — just
not the well either molecule came from. Without the check it lands in
tier A, is assigned, and nothing downstream could notice. The junction is
the evidence: a clean amplicon's interior is insert and contains no
adapter, so an anchor between the two outer anchors means the ends
cannot be trusted as a pair.

The check has a limit worth stating. It detects joins between *complete*
molecules, because those bury a full adapter structure at the junction. A
chimera of two **truncated fragments** contributes one usable end from
each and leaves no adapter between them — undetectable this way, or by
any method working from adapter position. `--no-chimera-check` disables
the quarantine and restores the old behaviour; that is a correctness
trade and belongs to whoever is analysing the run, not to the tool.

### Splitting concatemers

Quarantining those reads protects the results but throws away the data.
`split_concatemers.py` runs **before** demultiplexing and cuts each read
at its junctions, so every amplicon is demultiplexed on its own two
ends. The cut is the only edit — nothing is spliced or reconstructed —
and on a library without concatemers it is a no-op.

Measured against known provenance: chimeras built from reads whose true
sample was already established, on a 3×12 index grid.

| | units recovered | sent to a sample neither molecule came from |
|---|---|---|
| uncut, no check | 21% | **19.8%** |
| uncut, check on | 0% | 0% (refused instead) |
| split, then demultiplexed | 98% | 0% |

That middle row is the safeguard alone: correct, but empty-handed. The
last row is why splitting is worth the extra stage.

The contamination rate scales with how full the index grid is, because a
hybrid pair only misleads when it happens to be a combination some well
actually uses. At 3×12 with 24 of 36 combinations declared it is 19.8%;
a 96-well plate using 94 of 96 combinations would be worse. Leaving a
few wells unused is what turns this from invisible into measurable.

Splitting is the slowest stage in the pipeline — roughly 750 reads/s per
core, since every base of every read is scanned for adapter. On a full
flow cell that is hours, so it runs in parallel and resumes:

```bash
bin/split_concatemers.py --reads <minknow dir> --threads 16 --outdir split
```

Demultiplexing is sharded the same way, and for the same reason:

```bash
THREADS=16 bin/run_nanopore.sh results reference.tsv samplesheet.tsv runs.tsv primers.tsv
```

Work is divided by **input file**, which needs no chunking scheme
because a MinKNOW run already arrives as dozens of batch files. Each
worker writes its own outputs and its own tallies; a rerun skips every
input already finished, so a job killed by a wall-clock limit resumes
where it stopped rather than at the beginning.

Merging is cheap because of how the parts are shaped. gzip members
concatenate, so per-sample reads are joined by copying compressed bytes
rather than decompressing and recompressing the run; the per-read tables
are appended in input order; the tallies are summed. Parts are always
merged in input order, so **the thread count is not an input to the
result** — that is asserted byte for byte in the tests, for the per-read
table, every per-sample FASTQ, and every tally including the index
matrix, along with a resumed run reproducing an uninterrupted one.

Measured on 4 cores: splitting 3.3x, demultiplexing 3.7x.

A large tier B means the two ends are failing separately; a large tier D
means whole reads are failing. Neither says why, so:

```bash
bin/diagnose_demux.py --reads <minknow dir> --samplesheet samplesheet.tsv \
    --run-id RUN01 --assay assets/assay/ampseq_helminth.json
```

Read-only, seconds on a 20k-read subsample. For each end separately it
reports the step that lost it — anchor not found, anchor found but the
barcode window runs off the read, no barcode copy within tolerance,
copies disagreeing — checks whether the anchor sits where the construct
predicts, and sweeps the search window and anchor edit budget so the
value of loosening either is a measurement rather than an argument.

## Inputs

Four, of which one is optional and one is derivable.

### `--reference` — the panel

Sequence identity only; concentrations for quants, correction factors where
measured. Targets, control decoys, quantification standards and vector
backbone all live in one file, because competitive mapping requires a single
combined index — the point is that a read chooses between them.

```
seq_id  locus  type  taxon  parent_id  default_conc  conc_units  corr_factor  sequence
```

`type` is one of `target`, `control`, `quant`, `vector`. A control requires
`parent_id` resolving to a target at the same locus. A quant requires
`default_conc` and `conc_units`. `taxon` defaults to `seq_id` when blank, so
haplotype aggregation can be adopted later without a schema change.

### `--samplesheet` — what was sequenced

```
run_id,platform,sample_id,index_i7,index_i5,sample_type[,file_prefix]
```

`sample_type` is one of `unknown`, `neg_extraction`, `neg_ntc`, `mock`,
`pos_control`, and each drives specific behaviour rather than being
descriptive. `plate` and `well` are **not** columns — they are derived from
the index pair (see below).

### `--runs` — run-level properties

```
run_id  platform  path  barcode_repeats  dna_input_ng  lot_*
```

`--input_dir` is a shortcut that builds a single-row manifest, so the common
single-run case stays a one-liner and the multi-run case needs no new code
path. `lot_*` fields are unused by the analysis; they exist because lots
cannot be reconstructed later, and without them a step change in any
cross-run comparison is uninterpretable.

### `--assay` — the assay definition

`assets/assay/ampseq_helminth.json` holds everything that makes this a
particular assay: the control tag, the PCR1 overhangs and adapter anchors,
the index primer set and whether indices are tandem-duplicated, the
quantification design, and the tolerances tuned around them. Pointing the
pipeline at a different panel or construct is a file swap.

Both Nextflow and `preflight.py` read it, so values are defined once rather
than duplicated between a config and a script. `name` and `version` are
recorded in `preflight_summary.json`, so a result can be traced back to the
assay definition that produced it. A sparse config is fine — anything absent
falls back to built-in defaults, so adding a key later doesn't break existing
files.

Command-line flags override it for one-off runs (`--tag_max_edits`,
`--expected_tag`, `--barcode_repeats`); unset means "use the config".

**Three things are deliberately not configurable**: the `sample_type`
vocabulary, the `seq_type` vocabulary, and the structure of assignment. The
logic depends on them — `neg_extraction` sets thresholds, `mock` carries an
expected composition, `control` requires a parent. Allowing a sixth sample
type would mean every downstream branch has to handle a category it knows
nothing about. Configurability isn't free, and this is where it would cost
most.

### `--primers` — optional, for validating reference extent

A separate file from the assay config, because primers travel with the
*panel* rather than the construct: change the targets and the primers may
change, while the overhangs, tag and index design stay put.

```
name	direction	sequence
Helm_ITS2_F1	fwd	GCGTTATTTACCACGAATTGCAG
Helm_ITS2_R1	rev	TATGCTTAAGTTCAGCGGGT
```

It does exactly one job: check that every reference entry spans exactly
primer-to-primer. Primers aren't needed to assign a locus (competitive
mapping gives that) or to set length windows (reference length gives those) —
but they're the only way to detect **mixed reference extent**, where some
entries carry sequence outside their primer sites and so offer more alignable
bases than their competitors.

That failure hides inside genuine length variation. ITS2 legitimately spans
350–500 bp across genera, so no length-based check can separate "different
species" from "trimmed differently" — which is why the length-range check
that used to live here was removed rather than tuned. Pre-flight now reports,
per entry, which primer pair was found and how many bases lie outside it.

### `assets/index_map.tsv` — index primers and plate geometry

The AMPSEQ_IDXx2 primers encode plate position: i5 gives the row (A–P), i7
the column (1–24), which is a 384-well map. So `well` is a function of the
index pair, and recording it separately would only create an opportunity for
the two to disagree.

## What pre-flight produces

| Output | Contents |
|---|---|
| `reference.fasta` | combined index for competitive mapping |
| `tag_report.tsv/.txt` | control tags, positions, and predicted detection rate per platform |
| `file_resolution.tsv` | every sample with its matched FASTQ files |
| `preflight_summary.json` | machine-readable per-run record, schema-versioned |
| `validation.log` | errors and warnings |

Hard errors block the run; warnings do not. Everything is reported at once
rather than dying on the first problem and hiding the rest.

## Two findings worth knowing

**The index sets are Hamming-robust but not indel-robust.** Minimum Hamming
distance is 4 in both sets, so one substitution is safely correctable — but
minimum *edit* distance in the i7 set is only 2 (N721 vs N723). Since indels
are nanopore's dominant error mode, barcode matching must disable indels and
rely on the tandem duplication for indel robustness. Pre-flight computes this
and warns rather than trusting a supplied tolerance.

**Cross-set distance is 1.** i5 `S502 CTCTCTAT` and i7 `N707 CTCTCTAC` differ
by a single base, so barcode search must be anchored on adapter context and
never scanned across the read — otherwise the two ends of a molecule cannot
be told apart.

## Tandem versus single indices

`barcode_repeats` is 2 by default (the AMPSEQ tandem design) and can be set
to 1 for a conventional single index at each end. Both are dual-index; the
difference is whether each index carries its own internal replicate.

That replicate is **indel insurance, not extra error tolerance**. An indel
makes one copy uninformative and the other carries the read, and a frame slip
produces two copies that disagree — which is what detects it. With one copy
there is no second observation, so a slipped frame lands wherever it lands,
and the four rotational collisions in these index sets become live
misassignment routes rather than caught ones. Pre-flight says so explicitly
when it sees `barcode_repeats=1` on a nanopore run, and the recommended
matching strategy per run is written into `preflight_summary.json`.

None of it applies to Illumina, where the sequencer demultiplexes.

## Control tag detection

Control plasmids carry an inserted tag absent from their targets
(`TGGTCAATTTCAAGAGAT`). Target-versus-control is therefore a deterministic
presence test, not a score comparison — which matters, because an inserted
tag makes the control longer and a control-derived read aligned against the
target can still score highly enough to be misassigned.

Pre-flight derives each tag by aligning the control to its parent. One trap
handled explicitly: the insertion point is ambiguous whenever the tag shares
bases with the flanking sequence, so the insertion is left-aligned into its
canonical position (the same normalisation VCF uses for indels) and
cross-checked against `--expected_tag`. Without that, the derived tag comes
back as a *rotation* of the real one — self-consistent, and detection would
still work, but it would not match the tag as designed.

## Simulated data

```bash
python3 tests/make_fixtures.py --outdir example \
    --platform nanopore --n-samples 88 --leave-unused 8

python3 bin/simulate.py \
    --reference   example/reference.tsv \
    --samplesheet example/samplesheet.tsv \
    --runs        example/runs.tsv \
    --assay       assets/assay/ampseq_helminth.json \
    --outdir      example/reads \
    --reads-per-sample 1500 --error-model nanopore_hac
```

`bin/simulate.py` consumes the *same* reference panel, samplesheet, run
manifest and assay definition as the pipeline, so the two cannot drift apart
— a schema change breaks both at once rather than leaving the tests
validating a stale format.

It emits FASTQ (pre-split per sample for Illumina, one pooled file for
nanopore) plus a **per-read truth table**: which sample a read came from,
which reference sequence, and what was done to it. That is what makes the
read-accounting invariant testable against *intent* — a funnel can sum
perfectly while every read sits in the wrong bucket.

Pathologies are injected at known rates, so the QC sensors can be checked
against a known answer rather than inspected and believed:

| Flag | Default | What it creates |
|---|---|---|
| `--hop-rate` | 0.005 | index hopping, some onto unassigned combinations |
| `--chimera-rate` | 0.01 | ligation chimeras with mismatched ends |
| `--carryover-rate` | 0.002 | control plasmid appearing in a field sample |
| `--dimer-rate` | 0.02 | primer dimer |
| `--offtarget-rate` | 0.005 | vector-derived off-target product |
| `--truncation-rate` | 0.03 | partial reads, some not spanning the tag |
| `--no-spike` | off | omit the quants, to exercise the unverified-negative path |

Error models (`--error-model`) cover `illumina`, `nanopore_sup`,
`nanopore_hac` and `nanopore_fast`, with indels biased toward homopolymers
the way nanopore's actually are — which is what justifies edit-distance
matching over Hamming.

The fixture sequences are synthetic but structurally realistic: targets at a
locus share conserved flanks with the primers annealing in them, species
differ in the variable core, and haplotypes of one species differ by three
bases in 360. Fully random sequences would make competitive mapping trivially
separable, and a panel built from them would pass tests a real panel fails.

## Tests

```bash
python3 tests/test_preflight.py    # 143 assertions
python3 tests/test_simulate.py     #  43 assertions
python3 tests/test_demux.py        #  57 assertions
python3 tests/test_prep.py         #  32 assertions
python3 tests/test_assign.py       #  39 assertions
python3 tests/test_synthesis.py    #  65 assertions, ~90s
python3 tests/test_report.py       #  53 assertions
python3 tests/test_nextflow.py     #  14 static checks on the workflow
```

The pre-flight suite: unit tests for the assay config, tag derivation (including the ambiguous-boundary
case), distance helpers, plate geometry and sanitisation; then end-to-end runs
where each case breaks exactly one thing and asserts the specific complaint.
A validator that silently passes bad input is worse than no validator.

## Profiles

| Profile | Use |
|---|---|
| `local` | development, optimisation, and near-term nanopore runs |
| `farm` | routine Illumina work — LSF, Singularity |
| `test` | smoke test against the generated fixtures |

## Next increments

2. ~~Read simulator with a per-read truth table~~ — done
3. ~~Nanopore demultiplexing~~ — anchored tandem barcodes, four assignment tiers
4. ~~Illumina read path~~ — orient, trim, merge, dereplicate
5. ~~Shared assignment~~ — competitive mapping, tag test, declared ambiguity
6. ~~Run-level synthesis~~ — calibration, limits, three-state calling, QC
7. ~~Reporting~~ — self-contained run report, cross-run and paired-platform
   comparison

The simulator comes second deliberately: almost nothing downstream can be
verified without known truth, and the read-accounting invariant needs to be
testable against intent rather than merely self-consistent.

## Read extent: why the panel is trimmed twice

Pre-flight writes two FASTAs, and they are not interchangeable.

`reference.fasta` is the panel as supplied — full amplicons, primer sites
included. `reference.insert.fasta` is the same panel trimmed to between the
primers.

Reads are assigned against the second one. The primer bases in a read are
synthesised, not observed: whatever the template carried under a primer, the
read reports the primer. Keeping those bases would assert sequence that was
never measured, and — because two references differing only under a primer
would then look different — would promise a discrimination the assay cannot
deliver. Pre-flight now refuses a panel containing such a pair.

Trimming only one side is worse than trimming neither. Reads without primers
scored against references with them look uniformly ~40bp short, and
competitive mapping still returns a confident winner; it is simply measuring
the wrong thing, consistently.

Primers are located by approximate match, not substring search. A primer
rarely matches its template exactly — the primer may be degenerate and the
template may carry a mismatch under it. With exact matching the primer goes
unfound in exactly those entries, which are then left untrimmed while their
competitors are trimmed. The 40 extra alignable bases would land on the most
divergent sequences: the ones the panel exists to tell apart.

## Merging, and what it costs at 2x150

An insert merges when `2 x read_length - insert_length` clears the minimum
overlap. For the current panel, trimmed inserts run 305-477bp:

| Run | Inserts that merge | Smallest overlap |
|---|---|---|
| 2x300 | 22 of 22 | 123bp |
| 2x250 | 22 of 22 | **23bp** |
| 2x150 | none | — |

The 2x250 figure is marginal rather than comfortable. The longest insert
clears the 20bp minimum by three bases, so any quality trimming at the read
ends pushes it into the unmerged path — not as a clean switch, but
gradually, and only for the longest amplicons in the panel. That is a
composition bias that would arrive looking like biology.

Unmerged pairs are not discarded. They are kept as pairs and assigned with
both mates required to agree — the same concordance requirement the nanopore
demultiplexer applies to its two ends. It is a weaker measurement than a
merged read, because the middle of the amplicon is never observed, so the
run configuration is worth choosing deliberately rather than inheriting.

## Where the control tag has to sit

A control plasmid differs from its target by an inserted 18bp tag and
nothing else. That makes target-versus-control a presence test rather than
a comparison — but only if a read can actually see the tag.

On an amplicon too long for the mates to merge, the pair observes the two
ends and nothing in between. A tag in the middle falls in that gap. Both
references then fit the pair equally well, and nothing about the result
looks wrong: the read gets a confident best hit, to the target.

The panel is designed so this cannot happen: the tag sits flush against the
reverse primer in every control, making it the last 18bp of the insert and
therefore the *first* thing R2 reads after its own primer — in that read's
best bases, at any Illumina read length. That is the strongest available
placement and it was chosen deliberately.

The pipeline still has to handle the case, because the test fixture did not
follow that convention: it placed the tag at the midpoint, and on a 678bp
amplicon too long to merge, the first implementation treated a missing tag
as evidence of target. That turned **97 plasmid carryover reads into
genuine detections**, silently — the exact failure the tagged-decoy design
exists to prevent. The fixture has been corrected to match the real panel,
and the pathological placement now has a test of its own.

The rule is asymmetric now, and the asymmetry is the point:

| Tag found | Coverage over the tag site | Call |
|---|---|---|
| yes | — | control, deterministically |
| no | tag site observed | target |
| no | tag site not observed | `tag_unobservable` — declared, not guessed |

The third row is detectable without knowing where the tag is. If a read
ties between a target and its control, the tag site was necessarily
unobserved, because an observed tag would have separated them by its own
length. The tie is the evidence.

Pre-flight checks this before any sequencing happens. Given a read length,
it errors on any control whose amplicon cannot merge *and* whose tag lies
beyond a read length from either end, and names the fix: move the tag
toward an end, or sequence longer.

Three properties of the current panel are checked every run, because each
is deliberate and each could be undone by an edit elsewhere:

| Property | Status | Why it matters |
|---|---|---|
| tag readable by one mate | all 11 controls | survives 2x250, 2x150, 2x100 |
| tag inside the insert | all 11, margin 0 | anything under a primer is replaced by primer sequence |
| tag clear of the *other* reverse primers | all 11 | the panel uses R1, R2 and R3 across species |

The zero margin is the point of the design, not a defect — but it means a
primer redesign that shifts a boundary by one base puts the tag underneath
it. Pre-flight now reports this as a **note** rather than a warning: the
panel is correct, and saying so every run is how it stays correct. Filing
load-bearing facts as warnings only teaches people to ignore warnings.

## What assignment refuses to do

Assignment reports four things it cannot resolve rather than picking:

| Outcome | Meaning |
|---|---|
| `assigned` | one reference fits best by a clear margin |
| `taxon_only` | haplotype unresolved, but every candidate is one taxon — the species call stands |
| `ambiguous` | references from different taxa fit equally well |
| `tag_unobservable` | target and control tied, with the tag between the mates |
| `unmatched` | nothing in the panel within the distance limit |

Measured against the simulator's per-read truth: **zero misassignments** on
either platform. Illumina resolves essentially everything; at nanopore's
1% error rate the loss appears as `taxon_only` — haplotypes collapsing into
their species — and never as a wrong species or a control read reported as
its target.

Composition is computed over `speciation` targets only, so `rel_abundance`
sums to one. Quant spike-ins and control plasmids are counted, and counted
separately: a mock population built from plasmids yields control reads and
an empty composition table, which is what it should yield.

## Counts to calls

Synthesis is where the pipeline is most tempted to assert more than it
measured, so each of the three decisions carries its evidence.

### Quantification

Two spike-ins, two logs apart, go into every sample post-extraction at a
fixed dose per unit of faecal DNA. Read counts are compositional — a share
of the library, not an amount — and the standards are what convert one back
into the other. The unit is copies per unit DNA, because a fixed mass of
DNA is what was actually spiked.

Two points also make the calibration checkable. A linear assay fits a slope
of 1 through them; the fitted slope is reported for every well, because a
slope well below 1 means the high standard is saturating and every copy
number on that plate is compressed toward the middle.

Measured against the simulator's composition truth — input copies, recorded
before amplification, invisible to everything downstream:

| | |
|---|---|
| median recovery | 1.01x |
| within twofold of input | 96.8% |
| sensitivity / specificity | 0.855 / 0.993 |
| false detections in blanks and mocks | 0 |

The misses are at the bottom of the range: 11–15 true copies against a low
standard of 100. That is the assay's floor, not a defect in the code.

### Amplification correction, and the cost of not having one

`corr_factor` was being loaded and never applied, so the correction had
never actually been exercised. The simulator now applies efficiency between
copies and reads, which makes the correction testable by running the same
data twice:

| Reference | corr_factor | Uncorrected | Corrected |
|---|---|---|---|
| Ascaris | 1.80 | **1.87x too high** | 1.04x |
| Trichuris | 0.55 | **0.53x too low** | 0.96x |
| Necator ITS2 | 1.00 | 1.02x | 1.02x |

The bias is exactly the amplification factor, which is what it should be.

**`corr_factor` defaults to 1.0.** A blank column no longer means "no
correction" as a separate code path — dividing by 1 is what that already
did arithmetically — so every reference is on the same footing. Two columns
travel with each result: `corr_factor`, the number used, and
`corr_factor_source`, either `panel` or `default`. Measured-and-found-to-be-
unity and assumed-to-be-unity are the same number backed by different
evidence, and the distinction survives in the output rather than in the
arithmetic.

What the default costs is stated plainly by pre-flight and in `run_qc.json`:
with everything at 1.0, copy numbers are **relative to** the standards'
amplification efficiency rather than corrected to it. That is a sound basis
for comparing the same species across samples, which is what most of this
work needs. It is not a basis for comparing species with each other — the
table above is the size of the error when efficiencies genuinely differ.

Pre-flight distinguishes two situations. A panel with no `corr_factor`
anywhere gets a note. A panel where *some* entries carry measured values and
others do not gets a warning, because that puts corrected and assumed copy
numbers side by side in the same column.

### Detection

The limit is measured from the run's own negative controls, which carry the
same hopping, carryover and reagent contamination as everything else. Two
models are computed and the stricter applies: an absolute floor, because
hopped reads arrive at a rate set by the run rather than by the receiving
well, and a proportional one, because carryover scales with depth. The
call records which bound was active.

### The third state

`not_detected` is a claim about the sample. It is only available when the
well could have produced a positive, and the spike-ins are what establish
that. A well whose standards are missing is `assay_failed`, never
`not_detected`.

Distinguishing "the spike failed" from "there was no spike" is impossible
within one well — both give zero standard reads — but easy across a run. If
most wells recovered their standards, the few that did not have failed; if
none did, the protocol did not include a spike and every negative on the
run is reported as unverified.

## Plate-position QC, and what it is worth

Every well holds an unrelated person's sample, so well position should
explain nothing. Where it does, the cause is the plate — a multichannel
channel, a thermal gradient, evaporation at the edges — and not biology.
Row, column and edge effects are tested by permuting well labels, which is
a more honest null than any parametric one and needs no distributional
assumption about read counts.

A p-value is only worth reporting if it means what it says, so the test was
calibrated against simulated null data rather than assumed correct:

| | nominal | observed |
|---|---|---|
| false positives | 1% | 0.5% |
| | 5% | 4.5% |
| | 10% | 11.5% |

Power, at 88 wells, is **60%** against a twofold effect in two of eight
rows. So it will catch a gross plate problem and will miss a subtle one —
worth knowing before treating a clean result as evidence of a clean plate.
The p-values are unadjusted, and `run_qc.json` says so: with one test per
taxon per axis, roughly one in twenty falls below 0.05 by chance.

## Two runs of the same samples

The only external check this assay has is measuring something twice — a
replicate on another plate, or the same library on both platforms.
Everything else is the assay agreeing with itself.

`compare_runs.py` pairs by `sample_id`, with `--strip` taking a regex to
remove a replicate suffix and recover the subject. A subject that appears
twice on one side is reported as ambiguous rather than resolved by taking
whichever came first.

**Raw agreement is not the headline.** An assay whose answer is usually
`not_detected` agrees with itself easily: two runs that each say no 95% of
the time agree 90% of the time while being no better than chance. Cohen's
kappa is reported beside it. In the test suite, a constructed case with 90%
agreement comes out at kappa −0.05 — slightly *worse* than the marginal
rates predict, and describing that as "90% concordant" would be a fair
summary of the number and a false one about the assay.

Copy-number agreement is reported as a spread, not a correlation. A
correlation across four logs of burden looks excellent while individual
samples are out by tenfold, because the range does the work.

On a paired fixture — 24 samples, same libraries, Illumina and nanopore:

| | |
|---|---|
| paired observations | 120 |
| raw agreement | 0.925 |
| Cohen's kappa | 0.831 |
| both detected, within twofold | 97.1% |
| median log10 copy ratio | +0.011 |
| found by one platform only | 9, all nanopore |
| of those, found by the deeper run | 7 of 9 |

That last row is the point of recording it: most of the disagreement is a
difference in sampling depth, not between platforms.

### A bug this caught in the simulator

The first paired fixture returned **kappa 0.04** — no agreement beyond
chance between two runs of the *same libraries*. The comparison was right;
the fixture was not. The simulator drew a fresh composition for each run,
so a sample was a different sample on each platform.

Composition is now derived from the sample id, so a subject carries the
same true burden in every run it appears in — which is what "nanopore
sequences the same Illumina libraries" means. `test_report.py` asserts that
property directly rather than trusting it.

## The report

One self-contained HTML file per run. Self-contained is a requirement
rather than a preference: it has to open on a laptop in a field site with
no network, so there are no CDN links, no web fonts and no JavaScript.
Charts are inline SVG drawn from the data, and the test suite fails the
build if any external reference or script appears.

The report is ordered so that **what would invalidate the results comes
before the results**. Failed wells, unverified negatives, calibration
slopes away from unity, depth spread, plate-position structure and the
corr_factor basis are all in an opening section, because a run can fail any
of them and still produce a full set of confident-looking numbers.

Plate maps are drawn as plates. Well position is the one QC axis a table
cannot show: a dead row, a cold column, evaporation down one edge are
obvious in the shape and invisible in a sorted list.

Nothing in the report computes a new number — every figure is read from the
tables the pipeline already wrote, so it cannot disagree with them, and the
tests check the headline counts against `calls.tsv` directly.

### The composition figure: three panels, one axis

Depth, then where the reads went, then composition. The order is the
argument. The two panels
belong together, and separating them is how a composition gets over-read: a
stacked bar looks exactly the same whether it rests on eighty reads or
eighty thousand.

Two things the figure does beyond the conventional version:

**Bars are hatched when the composition is poorly determined.** The
threshold is on *target* reads, not total depth, because those are not the
same thing — a deep well can hold almost no target, and it is the target
count that the proportions are divided by. A sample with four target reads
would otherwise render as a confident 100%.

**Samples with no target reads are drawn empty**, not as an arbitrary
split. Zero over zero is not a composition.

Samples are ordered by well rather than by name, so a column effect appears
as neighbouring bars rather than scattered through the figure. Colours are
Okabe-Ito, which stays readable under the common forms of colour vision
deficiency; taxa beyond the palette are pooled as "other" in grey, and the
caption names them.

Composition is computed over `speciation` targets only, so a second locus
carried for variant calling never enters the denominator.

## Amplification efficiency, measured rather than assumed

`corr_factor` converts a read count into an amount. It can only be
measured from mixtures whose composition was decided before sequencing,
and `bin/calibrate_mocks.py` is where that happens — but the first
question it asks is whether a single number per species is the right
shape of answer at all.

For any pair of species, if bias is a fixed property of each template:

```
log(observed A/B) = log(input A/B) + log(eA/eB)
```

a line of slope 1 whose intercept is the correction. If the templates
compete instead, a dominant species is capped and the observed ratio
compresses toward the middle: the slope falls below 1, and no constant is
right anywhere except the one mixture it was fitted at. Real samples are
mostly one dominant species plus a minority, so a constant fitted at 50/50
fails exactly where it is needed — and fails by **inflating the minority**.

So the slope is tested first, and a `corr_factor` is emitted only when the
data support one.

| Simulated bias | Fitted slope | Verdict | corr_factor |
|---|---|---|---|
| multiplicative | 1.00 | constant | 1.75 [1.68–1.85] and 0.54 [0.52–0.57], against true 1.80 and 0.55 |
| competition 0.5 | **0.49** [0.48–0.51] | composition dependent | **none emitted** |

The fitted slope recovers the strength of the competition, not merely its
presence. Intervals come from resampling wells, so a correction that
cannot be told from noise reports as such.

### The plate has to be right before it is run

`assets/mocks_template.tsv` declares what went into each mock well, and
pre-flight validates the design while it is still a spreadsheet. On a
deliberately poor layout it raises all four of:

- only one mixture ratio — no slope can be fitted, so the correction is
  right only at that ratio
- no quant standard in the mock — efficiencies float instead of anchoring
  to the scale field samples are measured on
- every well in one PCR batch — replicates then measure pipetting rather
  than amplification
- targets appearing in no mock — those keep the default of 1.0

A 99:1 mixture also measures the minority limit of quantification, so the
extreme ratios pay for themselves twice.

**What this does not measure.** Templates are mixed before PCR, so this is
amplification bias alone. Applied to eggs in stool it corrects one step of
three; extraction and lysis are untouched.

## Nanopore input: a run is a directory, not a file

MinKNOW writes batched FASTQ into a tree, so `--reads` takes a file, a
glob, or a run directory. Every split in that tree is a decision about
which reads enter the analysis, so each is made explicitly and reported.

`fastq_fail` is excluded and counted. Those reads failed the basecaller's
own filter and carry a different error profile from the one the barcode
tolerances were set for; including them silently would change what the
demultiplexer is doing without saying so. `--include-failed` overrides it.

MinKNOW's own barcode folders are flagged. The barcodes in this assay are
in-line and tandem-duplicated inside the amplicon — MinKNOW knows nothing
about them and will have been matching kit barcodes that are not present,
so its `unclassified` folder holds the real data rather than the failures.
Everything is read, and the situation is reported.

That check originally never ran: it sat behind a "no files at this level"
shortcut, and in a real MinKNOW tree `fastq_pass` contains only
directories. It has a test now.

### Histograms that can be checked

The card figures in the report reconcile exactly with `calls.tsv` and
`run_qc.json` — that was verified by parsing the rendered page and
recomputing every number from the tables. The histograms did not, and four
defects were fixed:

**No y-axis.** A bar height was a shape, not a number. There was no way to
read a count off the figure or check it against anything. The axis is now
drawn, gridlines are labelled in samples, and the tallest bar carries its
count.

**Few distinct values were binned anyway.** Twenty-one calibration slopes
taking only two values — twenty at 1.00 and one at 0.85 — were spread
across twenty-six bins and drew as two bars at opposite ends of a wide
axis. That reads as a bimodal distribution with real spread; the truth is
"all identical except one". Eight or fewer distinct values are now drawn as
themselves, one labelled bar each.

**An invented axis.** When every value was identical the range collapsed
and was stretched to the next whole number, showing a spread that is not in
the data. It now says so in words.

**Mixed scales on a log axis.** A zero cannot be plotted logarithmically,
and the code left such values linear, dropping them into an arbitrary bin.
They are now excluded and counted in the caption.

Axis labels also read as numbers rather than `2.051e+04`.

## The read funnel

The per-sample figures start at *reads that reached assignment*, which is
already net of everything read preparation discarded. That hides one
specific failure: a well that produced plenty of data and lost it to
primer dimer draws the same short bar as a well that was simply loaded
light. The first needs the reaction fixed, the second needs the pool
rebalanced.

`synthesise.py` builds a run-level funnel from the sequencer to the
species call, naming the loss at each step. Both platforms are covered,
each starting from its own first count — `--prep-dir` for Illumina,
`--demux-summary` for nanopore:

```
read pairs from the sequencer      91,471  100.0%
passed read preparation            86,553   94.6%
    lost to no_primer               4,918
reached assignment                 86,553   94.6%
mapped to the panel                86,553   94.6%
mapped to a target                 81,787   89.4%
    lost to quant                   4,571
    lost to control                   195
```

Two checks come with it.

**A cross-boundary check**, on whichever path produced the reads. Both
platforms have the same seam: a step hands reads to assignment, and
nothing else in the pipeline compares the two counts, so a file dropped
between them would leave every stage balancing internally while data
quietly went missing. On Illumina that seam is preparation; on nanopore it
is demultiplexing, and tier A out of demux must equal the reads that
reached assignment.

**Per-well outliers.** A run median hides the one plate position that
failed, so any well losing more than half its reads before assignment is
named individually, with where they went. Those wells look merely shallow
in the depth panel; the funnel says they are not.

Both are tested by construction rather than by inspection: a healthy run
must raise neither, a run with one dimer-drowned well must name exactly
that well, and a run where reads vanish between tools must report the
discrepancy. Where no prep directory is supplied the funnel starts at
assignment and reports the boundary check as *not made* rather than as
passed.

### The nanopore funnel

```
reads in the run                    8,109  100.0%
demultiplexed to a sample           7,732   95.3%
    lost to tier B (one end only)         343
    lost to tier C (ends disagree)         32
    lost to tier D (no confident match)     2
reached assignment                  7,732   95.3%
mapped to the panel                 7,298   90.0%
    lost to no_primer                     417
mapped to a target                  2,978   36.7%
```

The tiers are named rather than pooled as "unassigned", because they point
at different problems: B is a read whose second barcode was unreadable, C
is a chimera, D is neither.

This path was silently broken when first written. The demux summary keys
its tiers under `by_tier`; the code looked for `tiers`, found nothing, and
skipped the demultiplexing stage without complaint — the funnel simply
began at assignment and nothing said a step was missing. It was found by
running the nanopore path end to end rather than by reading the code, and
there is now a test that constructs a demux summary and asserts the stage
appears.

## The Nextflow layer, and what is not yet verified about it

The workflow is written and **has never been executed**. Nextflow could not
be installed where this was developed, so `main.nf`, the two platform
subworkflows and the eight modules are the one part of the pipeline with
no runtime evidence behind them.

That is not an excuse for shipping them unchecked. Most of what goes wrong
in an unrun workflow is decidable without running it — a module passing
`--foo` to a tool whose flag is `--bar`, an `out.thing` nothing emits, an
undeclared `params.x` that silently reads as null and drops the flag from
the command line. `tests/test_nextflow.py` checks all of those by parsing
the `.nf` files and asking each tool for its real `--help`, so the list
cannot drift out of step by hand.

It found four things on its first run, and one of them was not a
workflow bug at all:

- sixteen parameters used but never declared — each would have read as
  null and quietly removed its flag, so the job would have *succeeded*
  with different settings than intended
- two references to pre-flight outputs that did not exist
- an output declared for a file the tool never wrote — which turned out to
  be the **i7 x i5 index matrix, specified early and never implemented**.
  It exists now: the full grid of observed index combinations, with
  unclaimed pairs marked, because a combination no sample owns cannot
  produce a real read and so measures cross-talk directly
- and a bug in the checker itself: `nf_files()` skipped any path
  containing `"work"`, which silently excluded `subworkflows/` from every
  check

### Splitting fans out across the farm

Concatemer splitting is the stage that walks every base of every read,
and on the first real MinION run it hit a twelve hour wall clock limit
without finishing. It is now one Nextflow task per input FASTQ rather
than one task with sixteen threads.

The reason it is per file and single core is that **the work divides by
input file**. `split_one_file()` is the unit; `--threads` only ever
sized a pool over those file-level jobs, and the tool says so itself
when given more threads than files:

```
note: only N file(s) to process, so at most N of the 16 threads can be used.
```

A shard holding one file therefore has nothing for a second core to do.
Asking for sixteen anyway would reserve them and leave them idle - on a
74 file run, over a thousand cores held for nothing - and it would also
be *slower*, because LSF places a 1-core job into whatever gap exists
now and a 16-core job only when a node frees up.

The stage is three processes:

| process | what it does | why separate |
|---|---|---|
| `SPLIT_LIST` | resolves the run to a list of FASTQ, one index and path per line | which files constitute a run is a decision with rules (`fastq_pass` in, `fastq_fail` out, pod5 ignored, barcode folders walked) and they live in `resolve_read_files()`. Nextflow asks the tool rather than globbing, so there is one answer and not two that agree until a run directory is reorganised |
| `SPLIT_CONCATEMERS` | one task per file, `cpus 1`, no `--threads` | the unit of work |
| `MERGE_SPLIT` | sums the per-shard statistics into `split_summary.json` | every read-funnel percentage downstream is taken against its `reads_in` |

**There is no merge of the reads.** The splitter writes one output FASTQ
per input FASTQ and demultiplexing walks a directory, so the shards
*are* the output - gathering them is Nextflow staging them side by
side, not a rewrite of tens of gigabytes.

Three things that had to be got right:

- **Output naming.** A shard sees one file and would number it `0000`,
  so every shard of a fanned-out run names its output alike - and two
  MinKNOW batches in different barcode folders can share a basename.
  `--part-index`, carried from `SPLIT_LIST`, gives a shard the name the
  whole-run route would have given it.
- **One aggregation, not two.** `build_summary()` is shared by both
  routes rather than reimplemented in the merger. The two quiet ways to
  get it wrong both yield believable numbers:
  `bases_retained_fraction` is a ratio of sums, and averaging per-file
  fractions instead weights a thousand-read batch equally with a
  million-read one; the pieces-per-read distribution is a histogram to
  be summed bin by bin.
- **A missing shard is fatal.** `--expect-parts` is checked, because
  losing one part of seventy-four raises no error on its own - it
  produces a read funnel that balances perfectly against a smaller run
  than the one that was sequenced. A duplicated shard is refused too,
  since it would double-count a file's reads.

`tests/test_split.py::test_sharded_matches_whole_run` runs both routes
on the same input and requires the output FASTQ to be identical file by
file and the summaries to be equal. The first version of that test
compared two sets of *empty* files and passed, because the fixture it
used was truncated - so it now also asserts that reads went in, pieces
came out, and at least one read was actually cut. A comparison of two
no-ops agrees for the wrong reason.

Demultiplexing is deliberately **not** fanned out this way. A sample's
reads are scattered across every batch file, so the gather would be
keyed by sample rather than by file - a different and larger change.

The driver script `bin/run_nanopore.sh` is unaffected and still uses
the tool's own threading; the two routes are tested against each other.

### The fan-out trap

The first fixture run produced **one** sample where the samplesheet has
eight, and said nothing about it. The log line was:

```
PREP_ILLUMINA (BLANK-01) [  0%] 0 of 1
```

`0 of 1` reads like a small test. It was seven wells disappearing.

Nextflow pairs a process's inputs **positionally and stops at the
shortest**. `Channel.fromPath` on a single file gives a *queue* channel
holding one item, and so does every process output. Pair a queue of
eight samples with a queue of one reference file and you get one task —
correctly, by Nextflow's own semantics, with no warning. Synthesis then
reports confidently on the sample that survived.

The fix is `.first()`, which converts a queue channel to a **value**
channel. A value channel is inexhaustible: it is re-read for every item
it is paired with and imposes no length on the pairing. In `main.nf`
every single-file input is converted at the point of creation rather
than at each call site, so a consumer added later cannot reintroduce it.

**A static check was not enough on its own.** The deeper problem was
that nothing in the pipeline ever asked whether the whole plate had
arrived. Synthesis reads `file_resolution.tsv` — the cohort the run
declared — and globs `*.counts.tsv` — the cohort that turned up — and
had never compared the two. It does now, before computing anything,
and a shortfall is fatal:

```
cohort incomplete: 8 samples declared for run 'RUN01', 1 reached synthesis.
  lost: MOCK-A, NTC-01, POS-01, iCCH-001-r1, ...
These samples have reads but no *.counts.tsv, so assignment never ran for them.
Refusing to report on a partial cohort.
```

The one legitimate shortfall is a nanopore well that received no reads:
an empty well demultiplexes to nothing and has no FASTQ to assign. That
case is distinguishable, because the demultiplexing summary counts
reads per sample — a well it never saw is empty, a well it wrote reads
for that then vanished is a lost task. Only the second is fatal; the
first is recorded in `run_qc.json` as `samples_absent`, with
`cohort_complete: false`, rather than left as a silent difference
between the samplesheet and the report. `--allow-missing` overrides,
and the QC still records that the plate was partial.

This guard protects the driver scripts as well as Nextflow, and would
have caught the original bug at the first fixture run.

The static check remains, because catching it before submission beats
catching it an hour in —
`test_fanout_is_not_capped` walks every call of a process whose first
input is a `tuple val(meta)`, follows the requirement outward through
subworkflow `take:` parameters to their callers, and rejects any scalar
argument that is not value-shaped. Re-introducing either of the two
original instances makes it fail.

### What the checker still cannot tell you

Whether a `collect()` blocks, whether the DAG deadlocks, whether a
`join` silently drops a key, whether `stageAs` lands files where a
process expects them. Those need Nextflow. On the first real run,
verify in this order:

1. `nextflow run . -profile test,local -stub-run` — every process has a
   stub, so this exercises the DAG's shape without touching data. A
   deadlock or a channel arity error shows up here in seconds.
2. `-profile test,local` on the fixtures — real data, one machine.
   **Count the samples before reading any result**: `cut -f2
   results-test/synthesis/calls.tsv | tail -n +2 | sort -u | wc -l`
   must equal the samplesheet, and `ls results-test/assign/*.counts.tsv
   | wc -l` must agree with it. A number that is merely plausible is
   not a pass.
3. A single real plate on `-profile farm`, checking that the per-sample
   tasks fan out rather than running serially.

The two places most likely to need adjustment are the `stageAs` patterns
in `SYNTHESISE` (which expects `assign/` and `prep/` directories assembled
from collected channels) and the `transpose()` in the nanopore subworkflow
that fans one demultiplexing task out into per-sample assignment.
