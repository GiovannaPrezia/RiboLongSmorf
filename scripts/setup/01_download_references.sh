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

GENOME_DIR="$PROJECT_ROOT/09_genome"
ANNOTATION_DIR="$PROJECT_ROOT/08_annotation"

mkdir -p "$GENOME_DIR" "$ANNOTATION_DIR"

GENOME_FA="$GENOME_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$ANNOTATION_DIR/gencode.v50.annotation.gtf"

echo "Downloading GENCODE v45 genome FASTA..."

if [[ ! -f "$GENOME_FA" ]]; then
    wget -O "$GENOME_FA.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"

    gunzip "$GENOME_FA.gz"
else
    echo "Genome FASTA already exists: $GENOME_FA"
fi

echo "Downloading GENCODE v45 annotation GTF..."

if [[ ! -f "$GTF" ]]; then
    wget -O "$GTF.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"

    gunzip "$GTF.gz"
else
    echo "GTF already exists: $GTF"
fi

echo "References downloaded successfully."
