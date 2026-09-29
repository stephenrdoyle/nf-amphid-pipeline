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

mkdir -p "$OUT"/{preflight,split,demux,assign,synthesis}

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
[ -e "$READS" ] || { echo "run $RUN_ID: reads path does not exist: $READS" >&2; exit 1; }

echo
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

n=0
echo
echo "== assign ($total samples) ====================================="
for fq in "${FQS[@]}"; do
    n=$((n + 1))
    base=$(basename "$fq" .fastq.gz)
    sample_id=$(unsanitise "$base")
    [ -n "$sample_id" ] || sample_id=$base

    if [ ! -s "$OUT/assign/$sample_id.counts.tsv" ]; then
        python3 "$HERE/assign.py" \
            --panel "$PANEL" --reference "$REF" --sample-id "$sample_id" \
            --reads "$fq" --primers "$PRIMERS" \
            --assay "$ASSAY" --outdir "$OUT/assign" --quiet
    fi
    printf '\r  [%d/%d] %-40s' "$n" "$total" "$sample_id"
done
echo

echo
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
    --outdir "$OUT/synthesis"

echo "== report ======================================================"
python3 "$HERE/report.py" \
    --synthesis-dir "$OUT/synthesis" --run-id "$RUN_ID" \
    --out "$OUT/${RUN_ID}_report.html"

echo
echo "Report:     $OUT/${RUN_ID}_report.html"
echo "Calls:      $OUT/synthesis/calls.tsv"
echo "Read fates: $OUT/synthesis/read_fates.tsv"
echo "Demux:      $OUT/demux/demux_summary.json"
