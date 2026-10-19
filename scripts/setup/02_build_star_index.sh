#!/bin/bash

set -eu

CONFIG_FILE="${1:-config.yaml}"

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

GENOME_FA="$PROJECT_ROOT/09_genome/GRCh38.primary_assembly.genome.fa"
GTF="$PROJECT_ROOT/08_annotation/gencode.v50.annotation.gtf"
STAR_INDEX="$PROJECT_ROOT/09_genome/hg38_star_index"

mkdir -p "$STAR_INDEX"

if [[ -f "$STAR_INDEX/Genome" ]]; then
    echo "STAR index already exists: $STAR_INDEX"
    exit 0
fi

if [[ ! -f "$GENOME_FA" ]]; then
    echo "ERROR: Genome FASTA not found: $GENOME_FA"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF not found: $GTF"
    exit 1
fi

echo "Building STAR genome index (this may take several minutes)..."

STAR \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 49

echo "STAR index generated in: $STAR_INDEX"
