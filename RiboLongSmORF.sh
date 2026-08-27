#!/bin/bash
set -e

# ==========================================================
# CHECK CONDA
# ==========================================================

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found in PATH."
    echo "Please install Miniconda or Anaconda first."
    exit 1
fi

CONFIG_FILE="${1:-config.yaml}"
ENV_NAME="ribolongsmorf_env"

# ==========================================================
# CHECK FILES
# ==========================================================

[[ -f "$CONFIG_FILE" ]] || {
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
}

[[ -f "environment.yml" ]] || {
    echo "ERROR: environment.yml not found."
    exit 1
}

[[ -f "ribolongsmorf_pipe.sh" ]] || {
    echo "ERROR: ribolongsmorf_pipe.sh not found."
    exit 1
}

[[ -f "scripts/setup/01_download_references.sh" ]] || {
    echo "ERROR: 01_download_references.sh not found."
    exit 1
}

[[ -f "scripts/setup/02_build_star_index.sh" ]] || {
    echo "ERROR: 02_build_star_index.sh not found."
    exit 1
}

[[ -f "scripts/ribotricer/01_run_ribotricer.sh" ]] || {
    echo "ERROR: 01_run_ribotricer.sh not found."

    echo ""
    echo "Please verify:"
    echo "scripts/ribotricer/01_run_ribotricer.sh"

    exit 1
}

[[ -f "scripts/ribotricer/02_process_ribotricer_results.py" ]] || {
    echo "ERROR: 02_process_ribotricer_results.py not found."

    echo ""
    echo "Please verify:"
    echo "scripts/ribotricer/02_process_ribotricer_results.py"

    exit 1
}

[[ -f "scripts/ribotricer/03_high_confidence_candidates.py" ]] || {
    echo "ERROR: 03_high_confidence_candidates.py not found."

    echo ""
    echo "Please verify:"
    echo "scripts/ribotricer/03_high_confidence_candidates.py"

    exit 1
}

[[ -f "scripts/qc_plots/09_qc_master_table.R" ]] || {
    echo "ERROR: 09_qc_master_table.R not found."

    echo ""
    echo "Please verify:"
    echo "scripts/qc_plots/09_qc_master_table.R"

    exit 1
}

# ==========================================================
# HEADER
# ==========================================================

cat << EOF

╔════════════════════════════════════════════════════════════╗
║                      RiboLongSmORF                        ║
║                                                          ║
║     Ribo-seq Processing & lncRNA-smORF Discovery         ║
║                                                          ║
║     Automated pipeline for translational profiling       ║
╚════════════════════════════════════════════════════════════╝

EOF

# ==========================================================
# CONDA ENVIRONMENT
# ==========================================================

echo ""
echo "[1/4] Checking Conda environment..."
echo ""

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then

    echo "✓ Environment '$ENV_NAME' found."

else

    echo "Environment '$ENV_NAME' not found."
    echo "Creating environment from environment.yml..."
    echo ""

    conda env create -f environment.yml

    echo ""
    echo "✓ Environment successfully created."

fi

# ==========================================================
# PROJECT PATHS
# ==========================================================

PROJECT_ROOT="$(
    conda run --no-capture-output -n "$ENV_NAME" \
        python -c \
        'import sys, yaml; print(yaml.safe_load(open(sys.argv[1]))["project_root"])' \
        "$CONFIG_FILE"
)"

GENOME_FA="$PROJECT_ROOT/09_genome/GRCh38.primary_assembly.genome.fa"
GTF="$PROJECT_ROOT/08_annotation/gencode.v45.annotation.gtf"
STAR_INDEX="$PROJECT_ROOT/09_genome/hg38_star_index"

# ==========================================================
# REFERENCES
# ==========================================================

echo ""
echo "[2/4] Checking reference files..."
echo ""

if [[ -s "$GENOME_FA" && -s "$GTF" ]]; then
    echo "✓ Genome FASTA already exists: $(readlink -f "$GENOME_FA")"
    echo "✓ GTF already exists: $(readlink -f "$GTF")"
else
    echo "One or more reference files are missing."
    echo "Downloading reference files..."

    conda run --no-capture-output -n "$ENV_NAME" \
        bash scripts/setup/01_download_references.sh "$CONFIG_FILE"

    echo ""
    echo "✓ Reference setup completed."
fi

# ==========================================================
# STAR INDEX
# ==========================================================

echo ""
echo "[3/4] Checking STAR genome index..."
echo ""

if [[ -s "$STAR_INDEX/Genome" &&
      -s "$STAR_INDEX/SA" &&
      -s "$STAR_INDEX/SAindex" &&
      -s "$STAR_INDEX/genomeParameters.txt" ]]; then

    echo "✓ Complete STAR index already exists:"
    echo "  $(readlink -f "$STAR_INDEX")"

else
    echo "STAR index was not found or is incomplete."
    echo "Building STAR genome index..."
    echo "This step can take 30-90+ minutes."
    echo ""

    conda run --no-capture-output -n "$ENV_NAME" \
        bash scripts/setup/02_build_star_index.sh "$CONFIG_FILE"

    echo ""
    echo "✓ STAR index ready."
fi

# ==========================================================
# PIPELINE
# ==========================================================

echo ""
echo "[4/4] Starting RiboLongSmORF pipeline..."
echo ""

conda run --no-capture-output -n "$ENV_NAME" \
    bash ribolongsmorf_pipe.sh "$CONFIG_FILE"

echo ""
echo "✓ RiboLongSmORF execution finished."
echo ""
