#!/usr/bin/env bash
# Run the Illumina path end to end, without Nextflow.
#
# Every stage here has been executed and tested; the Nextflow layer has
# not. Slower and serial, but it will finish and the numbers will be
# right. Use it until `nextflow run` has been shaken down.
#
#   bin/run_illumina.sh <outdir> <reference> <samplesheet> <runs> <primers>
#
# Resumable: a sample whose outputs already exist is skipped, so a run
# killed by a wall-clock limit can simply be started again.

set -euo pipefail

OUT=${1:?outdir}
REF=${2:?reference table}
SS=${3:?samplesheet}
RUNS=${4:?run manifest}
PRIMERS=${5:?primers}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ASSAY=${ASSAY:-$HERE/../assets/assay/ampseq_helminth.json}
IDXMAP=${IDXMAP:-$HERE/../assets/index_map.tsv}
RUN_ID=${RUN_ID:-$(awk 'NR==2{print $1}' "$RUNS")}

mkdir -p "$OUT"/{preflight,prep,assign,synthesis}

echo "== pre-flight =================================================="
python3 "$HERE/preflight.py" \
    --reference "$REF" --samplesheet "$SS" --runs "$RUNS" \
    --primers "$PRIMERS" --index-map "$IDXMAP" --assay "$ASSAY" \
    --outdir "$OUT/preflight"

# Pre-flight resolves every sample to its files and checks they exist, so
# the loop below reads its output rather than re-deriving paths. A file
# pre-flight blessed is the only file the rest of the run sees.
RES="$OUT/preflight/file_resolution.tsv"
PANEL="$OUT/preflight/reference.insert.fasta"
[ -s "$PANEL" ] || { echo "no primer-trimmed panel - was --primers given?" >&2; exit 1; }

n=0
total=$(( $(wc -l < "$RES") - 1 ))
echo "== prepare and assign ($total samples) ========================="
while IFS=$'\t' read -r run_id sample_id rest; do
    [ "$run_id" = "run_id" ] && continue
    n=$((n + 1))
    # The last field is the semicolon-joined list of resolved paths.
    files=${rest##*$'\t'}
    R1=$(echo "$files" | tr ';' '\n' | grep -E '_R1[._]' | head -1)
    R2=$(echo "$files" | tr ';' '\n' | grep -E '_R2[._]' | head -1)
    if [ -z "$R1" ] || [ -z "$R2" ]; then
        echo "  [$n/$total] $sample_id - no mate pair resolved, skipped" >&2
        continue
    fi

    if [ ! -s "$OUT/assign/$sample_id.counts.tsv" ]; then
        python3 "$HERE/prep_illumina.py" \
            --r1 "$R1" --r2 "$R2" --sample-id "$sample_id" \
            --primers "$PRIMERS" --assay "$ASSAY" \
            --outdir "$OUT/prep" --quiet
        python3 "$HERE/assign.py" \
            --panel "$PANEL" --reference "$REF" --sample-id "$sample_id" \
            --merged "$OUT/prep/$sample_id.merged.fasta.gz" \
            --pairs "$OUT/prep/$sample_id.pairs.tsv.gz" \
            --assay "$ASSAY" --outdir "$OUT/assign" --quiet
    fi
    printf '\r  [%d/%d] %-40s' "$n" "$total" "$sample_id"
done < "$RES"
echo

echo "== run-level synthesis ========================================="
python3 "$HERE/synthesise.py" \
    --assign-dir "$OUT/assign" --prep-dir "$OUT/prep" \
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
