#!/usr/bin/env bash
# Run the nanopore path end to end, without Nextflow.
#
# The companion to run_illumina.sh. The shape differs because the two
# platforms differ: Illumina arrives already demultiplexed, one pair of
# files per sample, so pre-flight can resolve files per sample and the
# work fans out immediately. A MinION run arrives as one pool with the
# indices still in the reads, so nothing is per-sample until the barcodes
# have been read. Demultiplexing is therefore a single serial stage, and
# only after it does the per-sample loop exist.
#
#   bin/run_nanopore.sh <outdir> <reference> <samplesheet> <runs> <primers>
#
# Resumable throughout. Splitting and demultiplexing are divided by input
# file and each records which files it finished, so a rerun picks up at
# the file it died on; a sample whose counts exist is skipped. A job
# killed by a wall-clock limit is simply started again.
#
#   THREADS=16 bin/run_nanopore.sh ...
#
# MIN_DEPTH and PERMUTATIONS pass through to synthesis. They exist here
# because the driver otherwise hardcodes what Nextflow exposes as
# params, and the two then cannot be compared like for like: run the
# same fixture both ways and every call differs, not because the
# pipelines disagree but because one used min_depth 10 and the other
# 100. That happened, and for a moment it looked like a real
# discrepancy between the sharded and threaded routes.
#
# sets the worker count for both parallel stages. They are the two that
# read every base of every read, and on a full flow cell they are hours
# of work - the rest is minutes.

set -euo pipefail

OUT=${1:?outdir}
REF=${2:?reference table}
SS=${3:?samplesheet}
RUNS=${4:?run manifest}
PRIMERS=${5:?primers}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ASSAY=${ASSAY:-$HERE/../assets/assay/ampseq_helminth.json}
IDXMAP=${IDXMAP:-$HERE/../assets/index_map.tsv}

# Per-stage wall clock. A driver that prints only headers cannot answer
# "where did the seventeen hours go", and on this pipeline the answer
# has been different every time - splitting, demultiplexing and the
# per-sample loop each dominated at some point.
RUN_T0=$(date +%s)
STAGE_T0=$RUN_T0
stage() {
    local now; now=$(date +%s)
    if [ -n "${STAGE_NAME:-}" ]; then
        printf '   [%s took %dm %ds]\n' "$STAGE_NAME" \
            $(( (now - STAGE_T0) / 60 )) $(( (now - STAGE_T0) % 60 ))
    fi
    STAGE_NAME=$1
    STAGE_T0=$now
}
finish() {
    stage ""
    local now; now=$(date +%s)
    printf '\nTotal: %dh %dm\n' $(( (now - RUN_T0) / 3600 )) \
        $(( ((now - RUN_T0) % 3600) / 60 ))
}

mkdir -p "$OUT"/{preflight,split,demux,assign,synthesis}

stage "pre-flight"
echo "== pre-flight =================================================="
python3 "$HERE/preflight.py" \
    --reference "$REF" --samplesheet "$SS" --runs "$RUNS" \
    --primers "$PRIMERS" --index-map "$IDXMAP" --assay "$ASSAY" \
    --outdir "$OUT/preflight"

RES="$OUT/preflight/file_resolution.tsv"
PANEL="$OUT/preflight/reference.insert.fasta"
[ -s "$PANEL" ] || { echo "no primer-trimmed panel - was --primers given?" >&2; exit 1; }

# The run manifest, not the samplesheet, says where the reads are: on this
# platform the reads belong to the run, not to any one sample. Take the
# first nanopore run in the manifest unless RUN_ID names another.
IFS=$'\t' read -r RUN_ID READS REPEATS < <(
    awk -F'\t' -v OFS='\t' -v want="${RUN_ID:-}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) col[$i] = i
            next
        }
        col["platform"] && $(col["platform"]) != "nanopore" { next }
        want != "" && $(col["run_id"]) != want { next }
        {
            print $(col["run_id"]), $(col["path"]),
                  (col["barcode_repeats"] ? $(col["barcode_repeats"]) : "")
            exit
        }
    ' "$RUNS"
) || true

if [ -z "${RUN_ID:-}" ]; then
    echo "no nanopore run found in $RUNS." >&2
    echo "The platform column must say 'nanopore' - if this plate went on" >&2
    echo "the MinION, both the run manifest and the samplesheet need" >&2
    echo "changing, not just one." >&2
    exit 1
fi
# A relative path in the manifest is relative to the MANIFEST, not to
# wherever this script was invoked from. Pre-flight has always resolved
# it that way; this driver and the nanopore branch of main.nf both read
# the manifest directly and did not, so the same runs.tsv worked in one
# place and not another. Third occurrence of this bug, so it is now
# fixed in all three.
case "$READS" in
    /*) ;;
    *)  READS="$(cd "$(dirname "$RUNS")" && pwd)/$READS" ;;
esac
[ -e "$READS" ] || {
    echo "run $RUN_ID: reads path does not exist: $READS" >&2
    echo "A relative path is taken relative to the run manifest" >&2
    echo "($(cd "$(dirname "$RUNS")" && pwd)), not to \$PWD." >&2
    exit 1
}

echo
stage "split"
echo "== split concatemers ($RUN_ID) ================================="
echo "   reads: $READS"
# Ligation joins amplicons into single reads. Left whole, the outer ends
# of such a read are indices from two different molecules, and the pair
# they form is usually a real well - so the read demultiplexes into a
# sample that never contained it. Cutting first is what makes the index
# pair mean something again. On a library with no concatemers this is a
# no-op and the reads pass through unchanged.
# Resumable per input file: a MinKNOW run arrives as dozens of batch
# files, each is split independently, and a rerun skips the ones already
# finished. THREADS sets the worker count.
python3 "$HERE/split_concatemers.py" \
    --reads "$READS" --assay "$ASSAY" \
    ${REPEATS:+--barcode-repeats "$REPEATS"} \
    ${THREADS:+--threads "$THREADS"} \
    --outdir "$OUT/split"

echo
stage "demultiplex"
echo "== demultiplex ($RUN_ID) ======================================="
# Also sharded by input file and resumable per file, so this no longer
# has to be redone wholesale after a timeout.
python3 "$HERE/demux_ont.py" \
    --reads "$OUT/split/reads" \
    --samplesheet "$SS" \
    --run-id "$RUN_ID" \
    --assay "$ASSAY" \
    ${REPEATS:+--barcode-repeats "$REPEATS"} \
    ${THREADS:+--threads "$THREADS"} \
    --outdir "$OUT/demux"

# If demultiplexing recomputed anything, the per-sample FASTQ have
# changed and any counts already sitting in assign/ were derived from
# the old ones. Resuming past that point would mix two vintages in one
# results directory, which is exactly how a run ends up 88% one
# algorithm and 12% another without saying so.
RECOMPUTED=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1])).get('parts_recomputed',0))" "$OUT/demux/demux_summary.json" 2>/dev/null || echo 0)
if [ "${RECOMPUTED:-0}" -gt 0 ]; then
    n_stale=$(ls "$OUT"/assign/*.counts.tsv 2>/dev/null | wc -l)
    if [ "$n_stale" -gt 0 ]; then
        echo "   demultiplexing recomputed $RECOMPUTED part(s); discarding"
        echo "   $n_stale stale assignment result(s) derived from the old reads"
        rm -f "$OUT"/assign/*.counts.tsv "$OUT"/assign/*.assign.json
    fi
fi

# Demultiplexing names its files after the sanitised sample id, because a
# sample id is free text and a filename is not. Map back through pre-flight
# rather than guessing, so the id that reaches synthesis is the id the
# samplesheet uses and the join cannot silently miss.
unsanitise() {
    awk -F'\t' -v want="$1" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $(col["sanitised_id"]) == want { print $(col["sample_id"]); exit }
    ' "$RES"
}

shopt -s nullglob
FQS=("$OUT"/demux/reads/*.fastq.gz)
shopt -u nullglob
total=${#FQS[@]}
if [ "$total" -eq 0 ]; then
    echo "demultiplexing produced no per-sample reads." >&2
    echo "See $OUT/demux/demux_summary.json - the tier breakdown says why." >&2
    exit 1
fi

echo
stage "assign"
echo "== assign ($total samples) ====================================="
# Samples are independent, so this fans out too. Each command is built
# with %q so a sample id containing a space or a quote cannot break the
# line it is written into, and a sample whose counts already exist is
# never queued - which is what keeps this resumable alongside the rest.
JOBS="$OUT/assign/jobs.sh"
: > "$JOBS"
todo=0
for fq in "${FQS[@]}"; do
    base=$(basename "$fq" .fastq.gz)
    sample_id=$(unsanitise "$base")
    [ -n "$sample_id" ] || sample_id=$base
    [ -s "$OUT/assign/$sample_id.counts.tsv" ] && continue
    printf 'python3 %q/assign.py --panel %q --reference %q --sample-id %q --reads %q --primers %q --assay %q --outdir %q --quiet\n' \
        "$HERE" "$PANEL" "$REF" "$sample_id" "$fq" "$PRIMERS" "$ASSAY" \
        "$OUT/assign" >> "$JOBS"
    todo=$((todo + 1))
done
echo "  $total samples, $((total - todo)) already done, $todo to do on ${THREADS:-1} thread(s)"
if [ "$todo" -gt 0 ]; then
    xargs -a "$JOBS" -d '\n' -P "${THREADS:-1}" -I CMD bash -c CMD
fi

echo
stage "synthesis"
echo "== run-level synthesis ========================================="
# No --prep-dir: there is no merge or pair step on this platform. The read
# funnel is fed by the demultiplexing summary instead, which is where the
# reads that never reached a sample are accounted for.
python3 "$HERE/synthesise.py" \
    --assign-dir "$OUT/assign" \
    --demux-summary "$OUT/demux/demux_summary.json" \
    --split-summary "$OUT/split/split_summary.json" \
    --samplesheet "$SS" --reference "$REF" \
    --preflight-dir "$OUT/preflight" --run-id "$RUN_ID" \
    ${CONTROLS_AS_TARGETS:+--controls-as-targets} \
    ${PLASMID_WARN:+--plasmid-warn "$PLASMID_WARN"} \
    ${MIN_DEPTH:+--min-depth "$MIN_DEPTH"} \
    ${PERMUTATIONS:+--permutations "$PERMUTATIONS"} \
    --outdir "$OUT/synthesis"

stage "report"
echo "== report ======================================================"
python3 "$HERE/report.py" \
    --synthesis-dir "$OUT/synthesis" --run-id "$RUN_ID" \
    --out "$OUT/${RUN_ID}_report.html"

finish
echo
echo "Report:     $OUT/${RUN_ID}_report.html"
echo "Calls:      $OUT/synthesis/calls.tsv"
echo "Read fates: $OUT/synthesis/read_fates.tsv"
echo "Demux:      $OUT/demux/demux_summary.json"
