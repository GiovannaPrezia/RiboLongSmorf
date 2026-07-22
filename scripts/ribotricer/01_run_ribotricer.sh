#!/bin/bash

set -eu

CONFIG_FILE="${1:-config.yaml}"

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

get_yaml () {
python - "$CONFIG_FILE" "$1" <<'PY'
import sys, yaml

with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)

value = data
for part in sys.argv[2].split("."):
    value = value[part]

print(value)
PY
}

PROJECT_ROOT=$(get_yaml "project_root")
THREADS=$(get_yaml "threads")

GTF="$PROJECT_ROOT/08_annotation/gencode.v45.annotation.gtf"
GENOME_FA="$PROJECT_ROOT/09_genome/GRCh38.primary_assembly.genome.fa"

ALIGN_DIR="$PROJECT_ROOT/05_alignment/ribo_seq"
RIBOTRICER_DIR="$PROJECT_ROOT/10_Ribotricer"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$RIBOTRICER_DIR" "$LOG_DIR"

CANDIDATE_PREFIX="$RIBOTRICER_DIR/gencode_v45"
CANDIDATE_ORFS="${CANDIDATE_PREFIX}_candidate_orfs.tsv"

echo "======================================"
echo "Running Ribotricer"
echo "======================================"

[[ -f "$GTF" ]] || {
    echo "ERROR: GTF not found: $GTF"
    exit 1
}

[[ -f "$GENOME_FA" ]] || {
    echo "ERROR: Genome FASTA not found: $GENOME_FA"
    exit 1
}

if [[ ! -f "$CANDIDATE_ORFS" ]]; then
    echo "Preparing candidate ORFs..."

    ribotricer prepare-orfs \
        --gtf "$GTF" \
        --fasta "$GENOME_FA" \
        --prefix "$CANDIDATE_PREFIX" \
        > "$LOG_DIR/ribotricer_prepare_orfs.log" 2>&1

    echo "Candidate ORFs generated: $CANDIDATE_ORFS"
else
    echo "Candidate ORFs already exist: $CANDIDATE_ORFS"
fi

for MODE in all_lengths 28_36; do

    BAM_DIR="$ALIGN_DIR/$MODE"
    OUT_MODE_DIR="$RIBOTRICER_DIR/$MODE"

    mkdir -p "$OUT_MODE_DIR"

    if [[ ! -d "$BAM_DIR" ]]; then
        echo "Skipping mode $MODE: BAM directory not found."
        continue
    fi

    for BAM in "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam; do

        [[ -f "$BAM" ]] || continue

       SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
OUT_PREFIX="$OUT_MODE_DIR/$SAMPLE"

RIBOTRICER_OUTPUT="${OUT_PREFIX}_translating_ORFs.tsv"

if [[ ! -f "$RIBOTRICER_OUTPUT" ]]; then

    echo "Ribotricer detect-orfs: $SAMPLE"

    ribotricer detect-orfs \
        --ribotricer_index "$CANDIDATE_ORFS" \
        --bam "$BAM" \
        --prefix "$OUT_PREFIX" \
        > "$LOG_DIR/${SAMPLE}.ribotricer.log" 2>&1

else

    echo "Ribotricer output already exists, skipping detect-orfs: $SAMPLE"

fi

if [[ -f "$RIBOTRICER_OUTPUT" ]]; then

    echo "Processing Ribotricer results: $SAMPLE"

python "$PIPELINE_DIR/scripts/ribotricer/02_process_ribotricer_results.py" \
    "$RIBOTRICER_OUTPUT" \
    "$CANDIDATE_ORFS" \
    "$GTF" \
    "$OUT_PREFIX" \
    "$SAMPLE"

else

    echo "WARNING: Ribotricer output not found for $SAMPLE: $RIBOTRICER_OUTPUT"

fi

    done

    HC_OUTPUT="$OUT_MODE_DIR/${MODE}_high_confidence_candidates.tsv"

    echo "Generating high-confidence candidates for mode: $MODE"

    python "$PIPELINE_DIR/scripts/ribotricer/03_high_confidence_candidates.py" \
        "$OUT_MODE_DIR" \
        "$HC_OUTPUT"

done

echo "Ribotricer finished."
