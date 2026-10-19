#!/bin/bash

set -eu

CONFIG_FILE="${1:-config.yaml}"

[[ -f "$CONFIG_FILE" ]] || {
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
}

get_yaml () {
python - "$CONFIG_FILE" "$1" <<'PY'
import sys, yaml

config_file = sys.argv[1]
key = sys.argv[2]

with open(config_file) as f:
    data = yaml.safe_load(f)

value = data
for part in key.split("."):
    value = value[part]

print(value)
PY
}

PROJECT_NAME=$(get_yaml "project_name")
PROJECT_DESCRIPTION=$(get_yaml "project_description")
BASE_DIR=$(get_yaml "project_root")
THREADS=$(get_yaml "threads")

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$PIPELINE_DIR/scripts"

# Fixed reference paths
GENOME_FA="$BASE_DIR/09_genome/GRCh38.primary_assembly.genome.fa"
GTF="$BASE_DIR/08_annotation/gencode.v50.annotation.gtf"
STAR_INDEX="$BASE_DIR/09_genome/hg38_star_index"

RNA_DICT_FASTA="$PIPELINE_DIR/resources/rnas_dictionary/rnas_dictionary_human.fa"
RNA_DICT="$PIPELINE_DIR/resources/rnas_dictionary/indexes/bowtie1/rnas_dictionary_human"

mapfile -t SIZE_MODES < <(python - "$CONFIG_FILE" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)
for mode in data["size_modes"]:
    print(mode)
PY
)

mapfile -t SAMPLES < <(python - "$CONFIG_FILE" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)
for sample in data["samples"]:
    print(sample["sample_name"])
PY
)

mapfile -t RUNS < <(python - "$CONFIG_FILE" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)
for sample in data["samples"]:
    print(sample["run_id"])
PY
)

declare -A RUN_MAP

while read -r RUN SAMPLE; do
    RUN_MAP["$RUN"]="$SAMPLE"
done < <(python - "$CONFIG_FILE" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)
for sample in data["samples"]:
    print(sample["run_id"], sample["sample_name"])
PY
)

SRA_DIR="$BASE_DIR/01_SRAs/"
RAW_DIR="$BASE_DIR/02_fastq/"
FASTQC_RAW_DIR="$RAW_DIR/fastqc_raw"

TRIM_DIR="$BASE_DIR/03_trimmed/"
FASTQC_TRIM_DIR="$TRIM_DIR/fastqc_trimmed"

CLEAN_DIR="$BASE_DIR/04_cleaned/"
FASTQC_CLEAN_DIR="$CLEAN_DIR/fastqc_cleaned"

ALIGN_DIR="$BASE_DIR/05_alignment/"
STAR_QC_DIR="$BASE_DIR/06_star_qc/"

COUNT_DIR="$BASE_DIR/07_counts/"
RIBOTRICER_DIR="$BASE_DIR/10_Ribotricer"
MULTIQC_DIR="$BASE_DIR/11_MultiQC"

FIG_DIR="$BASE_DIR/12_QC_Figures"

REPORT_DIR="$BASE_DIR/13_Report"
REPORT_TABLE_DIR="$REPORT_DIR/tables"
REPORT_PDF_DIR="$REPORT_DIR/pdf"
REPORT_HTML_DIR="$REPORT_DIR/html"

LOG_DIR="$BASE_DIR/logs"
QC_DIR="$BASE_DIR/QC_tables"

mkdir -p \
"$SRA_DIR" \
"$RAW_DIR" "$FASTQC_RAW_DIR" \
"$TRIM_DIR" "$FASTQC_TRIM_DIR" \
"$CLEAN_DIR" "$FASTQC_CLEAN_DIR" \
"$ALIGN_DIR" "$STAR_QC_DIR" \
"$COUNT_DIR" "$RIBOTRICER_DIR" "$MULTIQC_DIR" \
"$LOG_DIR" "$QC_DIR" \
"$FIG_DIR" \
"$REPORT_DIR" "$REPORT_TABLE_DIR" "$REPORT_PDF_DIR" "$REPORT_HTML_DIR"

[[ -f "$GTF" ]] || { echo "ERROR: GTF file not found: $GTF"; exit 1; }
[[ -f "$GENOME_FA" ]] || { echo "ERROR: Genome FASTA not found: $GENOME_FA"; exit 1; }
[[ -f "$STAR_INDEX/Genome" ]] || { echo "ERROR: STAR index not found: $STAR_INDEX"; exit 1; }
[[ -f "${RNA_DICT}.1.ebwt" ]] || { echo "ERROR: Bowtie contaminant index not found: ${RNA_DICT}.1.ebwt"; exit 1; }

echo ""
echo "============================================================"
echo "                    RiboLongSmORF"
echo "============================================================"
echo ""
printf "%-13s: %s\n" "Project" "$PROJECT_NAME"
printf "%-13s: %s\n" "Description" "$PROJECT_DESCRIPTION"
printf "%-13s: %s\n" "Threads" "$THREADS"
printf "%-13s: %s\n" "Size Modes" "$(printf "%s | " "${SIZE_MODES[@]}" | sed 's/ | $//')"
echo ""
echo "============================================================"
echo ""
echo "1 - Run complete pipeline"
echo "2 - Run step-by-step mode"
echo ""

read -p "Select option: " MODE

if [[ "$MODE" == "1" ]]; then
    PIPELINE_MODE="continuous"
    MODULE="12"
elif [[ "$MODE" == "2" ]]; then
    PIPELINE_MODE="interactive"

    echo ""
    echo "Available modules:"
    echo "0 - Download SRA file and convert to FASTQ"
    echo "1 - Raw FastQC"
    echo "2 - Cutadapt trimming"
    echo "3 - Post-trimming QC"
    echo "4 - noN/noPolyG filtering + all_lengths/25_40"
    echo "5 - Bowtie contaminant removal"
    echo "6 - FastQC on cleaned reads"
    echo "7 - STAR alignment"
    echo "8 - featureCounts"
    echo "9  - MultiQC + QC Figures"
    echo "10 - Ribotricer lncRNA-smORF discovery"
    echo "11 - Final QC report"
    echo "12 - Complete pipeline"
    echo ""
    

    read -p "Select module: " MODULE
else
    echo "Invalid option."
    exit 1
fi

MASTER_LOG="$LOG_DIR/${PROJECT_NAME}_pipeline_master.log"
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")

echo "" | tee "$MASTER_LOG"
echo "============================================================" | tee -a "$MASTER_LOG"
printf "%-13s: %s\n" "Pipeline Mode" "$PIPELINE_MODE" | tee -a "$MASTER_LOG"
printf "%-13s: %s\n" "Module" "$MODULE" | tee -a "$MASTER_LOG"
printf "%-13s: %s\n" "Started" "$START_TIME" | tee -a "$MASTER_LOG"
echo "============================================================" | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

CONTINUE_ALL=false

pause_step () {

    if [[ "$PIPELINE_MODE" == "interactive" ]]; then
        echo ""
        while true; do
            read -p "Continue (c), switch to continuous mode (s), pause for manual adjustment (r), or quit (q)? " ANSWER
            case "$ANSWER" in
                c|C) break ;;
                s|S)
    PIPELINE_MODE="continuous"
    CONTINUE_ALL=true
    echo "Continuous mode activated."
    break
    ;;
                r|R) echo "Pipeline paused for manual adjustment."; exit 0 ;;
                q|Q) echo "Pipeline stopped."; exit 0 ;;
                *) echo "Invalid answer. Use c, s, r, or q." ;;
            esac
        done
    fi
}

run_sra_conversion () {
    echo "Downloading SRA files and converting to FASTQ..." | tee -a "$MASTER_LOG"

    for RUN in "${RUNS[@]}"; do
        SAMPLE="${RUN_MAP[$RUN]}"
        SRA_FILE="$SRA_DIR/$RUN/$RUN.sra"
        OUT_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"

        echo "Run: $RUN → $SAMPLE" | tee -a "$MASTER_LOG"

        if [[ -f "$OUT_FASTQ" ]]; then
            echo "FASTQ already exists, skipping: $OUT_FASTQ" | tee -a "$MASTER_LOG"
            continue
        fi

        if [[ ! -f "$SRA_FILE" ]]; then
            echo "SRA file not found. Downloading with prefetch..." | tee -a "$MASTER_LOG"

            prefetch "$RUN" \
                --output-directory "$SRA_DIR" \
                > "$LOG_DIR/${RUN}.prefetch.log" 2>&1
        else
            echo "SRA file already exists: $SRA_FILE" | tee -a "$MASTER_LOG"
        fi

        [[ -f "$SRA_FILE" ]] || {
            echo "ERROR: SRA file still not found after prefetch: $SRA_FILE" | tee -a "$MASTER_LOG"
            exit 1
        }

        echo "Converting SRA to FASTQ..." | tee -a "$MASTER_LOG"

        fasterq-dump "$SRA_FILE" \
            -O "$RAW_DIR" \
            -e "$THREADS" \
            > "$LOG_DIR/${RUN}.fasterq_dump.log" 2>&1

        if [[ -f "$RAW_DIR/${RUN}.fastq" ]]; then
            mv "$RAW_DIR/${RUN}.fastq" "$RAW_DIR/${SAMPLE}.fastq"
            pigz -p "$THREADS" "$RAW_DIR/${SAMPLE}.fastq"
        else
            echo "ERROR: FASTQ was not generated for $RUN" | tee -a "$MASTER_LOG"
            exit 1
        fi

        seqkit stats "$OUT_FASTQ" | tee "$QC_DIR/${SAMPLE}.seqkit_raw_stats.txt"
    done
}

run_fastqc_raw () {
    SAMPLE="$1"
    RAW_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"

    FASTQC_HTML="$FASTQC_RAW_DIR/${SAMPLE}_fastqc.html"
    FASTQC_ZIP="$FASTQC_RAW_DIR/${SAMPLE}_fastqc.zip"

    if [[ -s "$FASTQC_HTML" && -s "$FASTQC_ZIP" ]]; then
        echo "Raw FastQC already completed, skipping: $SAMPLE" \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Raw FastQC: $SAMPLE" | tee -a "$MASTER_LOG"

        fastqc "$RAW_FASTQ" \
            -o "$FASTQC_RAW_DIR" \
            -t "$THREADS" \
            > "$LOG_DIR/${SAMPLE}.fastqc_raw.log" 2>&1
    fi
}

run_cutadapt_interactive() {
    SAMPLE="$1"

    RAW_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"
    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"
    REPORT="$LOG_DIR/${SAMPLE}.pre_cutadapt_qc.txt"
    CUTADAPT_LOG="$LOG_DIR/${SAMPLE}.cutadapt.log"

    if [[ -s "$TRIM_FASTQ" && -s "$CUTADAPT_LOG" ]] \
        && gzip -t "$TRIM_FASTQ" 2>/dev/null; then
        echo "Cutadapt already completed, skipping: $SAMPLE" \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "======================================"
        echo "Pre-QC before Cutadapt"
        echo "Sample: $SAMPLE"
        echo "======================================"

        zcat "$RAW_FASTQ" | awk '
            NR%4==2 {
                seq=$0
                n_reads++
                len=length(seq)

                if (min_len == "" || len < min_len) min_len=len
                if (len > max_len) max_len=len
                if (seq ~ /N/) n_with_N++
                if (seq ~ /AAAAAAAAAA/) n_polyA++

                if (n_reads == 100000) exit
            }

            END {
                print "Reads analyzed:", n_reads
                print "Observed length:", min_len "-" max_len " nt"
                print "Reads with N:", n_with_N "/" n_reads, "(" n_with_N/n_reads*100 "%)"
                print "Reads with polyA >=10:", n_polyA "/" n_reads, "(" n_polyA/n_reads*100 "%)"
            }
        ' | tee "$REPORT"

        echo "" | tee -a "$REPORT"
        echo "Top sequences:" | tee -a "$REPORT"

        zcat "$RAW_FASTQ" \
            | awk 'NR%4==2 {print $0}' \
            | head -n 100000 \
            | sort | uniq -c | sort -nr | head -n 10 \
            | awk '{print NR". "$2" ("$1" reads)"}' \
            | tee -a "$REPORT"

        echo "" | tee -a "$REPORT"
        echo "Most common first 3 nt:" | tee -a "$REPORT"

        zcat "$RAW_FASTQ" \
            | awk 'NR%4==2 {print substr($0,1,3)}' \
            | head -n 100000 \
            | sort | uniq -c | sort -nr | head -n 10 \
            | awk -v total=100000 \
                '{printf "%s: %.2f%% (%s reads)\n", $2, $1/total*100, $1}' \
            | tee -a "$REPORT"

        echo ""
        read -p "Remove first 3 nt? [y/n]: " REMOVE_3NT
        read -p "Remove polyA A{10}? [y/n]: " REMOVE_POLYA
        read -p "Minimum length? [20]: " MIN_LEN

        MIN_LEN=${MIN_LEN:-17}
        CUTADAPT_ARGS=()

        if [[ "$REMOVE_3NT" == "y" || "$REMOVE_3NT" == "Y" ||
              "$REMOVE_3NT" == "s" || "$REMOVE_3NT" == "S" ]]; then
            CUTADAPT_ARGS+=("-u" "3")
        fi

        if [[ "$REMOVE_POLYA" == "y" || "$REMOVE_POLYA" == "Y" ||
              "$REMOVE_POLYA" == "s" || "$REMOVE_POLYA" == "S" ]]; then
            CUTADAPT_ARGS+=("-a" "A{10}")
        fi

        CUTADAPT_ARGS+=("-m" "$MIN_LEN")

        cutadapt \
            "${CUTADAPT_ARGS[@]}" \
            -o "$TRIM_FASTQ" \
            "$RAW_FASTQ" \
            > "$CUTADAPT_LOG"

        tail -n 20 "$CUTADAPT_LOG"
    fi
}


run_cutadapt_default() {
    SAMPLE="$1"

    RAW_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"
    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"
    CUTADAPT_LOG="$LOG_DIR/${SAMPLE}.cutadapt.log"

    if [[ -s "$TRIM_FASTQ" && -s "$CUTADAPT_LOG" ]] \
        && gzip -t "$TRIM_FASTQ" 2>/dev/null; then
        echo "Cutadapt already completed, skipping: $SAMPLE" \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Cutadapt trimming: $SAMPLE" | tee -a "$MASTER_LOG"

        cutadapt \
            -u 3 \
            -a "A{10}" \
            -m 20 \
            -o "$TRIM_FASTQ" \
            "$RAW_FASTQ" \
            > "$CUTADAPT_LOG"

        tail -n 20 "$CUTADAPT_LOG"
    fi
}


run_cutadapt() {
    SAMPLE="$1"

    if [[ "$PIPELINE_MODE" == "interactive" ]]; then
        run_cutadapt_interactive "$SAMPLE"
    else
        run_cutadapt_default "$SAMPLE"
    fi
}


run_qc_trimmed() {
    SAMPLE="$1"

    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"
    LENGTH_REPORT="$QC_DIR/${SAMPLE}.length_distribution_trimmed.txt"
    TOP_REPORT="$QC_DIR/${SAMPLE}.top_sequences_trimmed.txt"
    FASTQC_HTML="$FASTQC_TRIM_DIR/${SAMPLE}.trim_fastqc.html"
    FASTQC_ZIP="$FASTQC_TRIM_DIR/${SAMPLE}.trim_fastqc.zip"

    if [[ -s "$LENGTH_REPORT" && -s "$TOP_REPORT" &&
          -s "$FASTQC_HTML" && -s "$FASTQC_ZIP" ]]; then
        echo "Post-trimming QC already completed, skipping: $SAMPLE" \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Post-trimming QC: $SAMPLE" | tee -a "$MASTER_LOG"

        if [[ ! -s "$LENGTH_REPORT" || ! -s "$TOP_REPORT" ]]; then
            zcat "$TRIM_FASTQ" \
                | head -n 4000000 \
                | awk 'NR%4==2 {print length($0)}' \
                | sort -n | uniq -c > "$LENGTH_REPORT"

            zcat "$TRIM_FASTQ" \
                | head -n 4000000 \
                | awk 'NR%4==2 {count[$0]++} END {for(seq in count) print count[seq], seq}' \
                | sort -nr | head -n 30 > "$TOP_REPORT"
        fi

        if [[ ! -s "$FASTQC_HTML" || ! -s "$FASTQC_ZIP" ]]; then
            fastqc "$TRIM_FASTQ" \
                -o "$FASTQC_TRIM_DIR" \
                -t "$THREADS" \
                > "$LOG_DIR/${SAMPLE}.fastqc_trimmed.log" 2>&1
        fi
    fi
}


run_filter() {
    SAMPLE="$1"

    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"
    FILTERED_ALL="$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.all_lengths.fastq.gz"
    FILTERED_RPF="$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.25_40.fastq.gz"

    if [[ -s "$FILTERED_ALL" && -s "$FILTERED_RPF" ]] \
        && gzip -t "$FILTERED_ALL" 2>/dev/null \
        && gzip -t "$FILTERED_RPF" 2>/dev/null; then
        echo "Filtering already completed, skipping: $SAMPLE" \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Filtering noN/noPolyG: $SAMPLE" | tee -a "$MASTER_LOG"

        zcat "$TRIM_FASTQ" \
            | awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(s !~ /N/ && s !~ /^G+$/) print h"\n"s"\n"p"\n"q}' \
            | gzip > "$FILTERED_ALL"

        zcat "$TRIM_FASTQ" \
            | awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(s !~ /N/ && s !~ /^G+$/ && length(s)>=25 && length(s)<=40) print h"\n"s"\n"p"\n"q}' \
            | gzip > "$FILTERED_RPF"
    fi
}

run_bowtie () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        INPUT_FASTQ="$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.${SIZE_MODE}.fastq.gz"
        MODE_CLEAN_DIR="$CLEAN_DIR/$SIZE_MODE"
        mkdir -p "$MODE_CLEAN_DIR"

        CLEAN_FASTQ="$MODE_CLEAN_DIR/${SAMPLE}.${SIZE_MODE}.clean.fastq"
        CLEAN_FASTQ_GZ="${CLEAN_FASTQ}.gz"
        BOWTIE_LOG="$LOG_DIR/${SAMPLE}.${SIZE_MODE}.bowtie.log"

        if [[ -s "$CLEAN_FASTQ_GZ" && -s "$BOWTIE_LOG" ]] \
            && gzip -t "$CLEAN_FASTQ_GZ" 2>/dev/null; then
            echo "Bowtie already completed, skipping: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"
            continue
        else
            echo "Bowtie contaminant removal: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"

            bowtie \
                -p "$THREADS" \
                -v 1 \
                --un "$CLEAN_FASTQ" \
                "$RNA_DICT" \
                "$INPUT_FASTQ" \
                /dev/null \
                2> "$BOWTIE_LOG"

            gzip -f "$CLEAN_FASTQ"
        fi
    done
}

run_fastqc_clean () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        CLEAN_FASTQ="$CLEAN_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}.clean.fastq.gz"
        FASTQC_HTML="$FASTQC_CLEAN_DIR/${SAMPLE}.${SIZE_MODE}.clean_fastqc.html"
        FASTQC_ZIP="$FASTQC_CLEAN_DIR/${SAMPLE}.${SIZE_MODE}.clean_fastqc.zip"
        FASTQC_LOG="$LOG_DIR/${SAMPLE}.${SIZE_MODE}.fastqc_clean.log"

        if [[ -s "$FASTQC_HTML" && -s "$FASTQC_ZIP" ]]; then
            echo "FastQC cleaned already completed, skipping: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"
            continue
        else
            echo "FastQC on cleaned reads: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"

            fastqc "$CLEAN_FASTQ" \
                -o "$FASTQC_CLEAN_DIR" \
                -t "$THREADS" \
                > "$FASTQC_LOG" 2>&1
        fi
    done
}


run_star () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        MODE_ALIGN_DIR="$ALIGN_DIR/$SIZE_MODE"
        mkdir -p "$MODE_ALIGN_DIR"

        CLEAN_FASTQ="$CLEAN_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}.clean.fastq.gz"

        BAM="$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Aligned.sortedByCoord.out.bam"
        BAI="${BAM}.bai"
        STAR_LOG="$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Log.final.out"
        IDXSTATS="$STAR_QC_DIR/${SAMPLE}.${SIZE_MODE}.idxstats.txt"
        STAR_QC_LOG="$STAR_QC_DIR/${SAMPLE}.${SIZE_MODE}_Log.final.out"

        if [[ -s "$BAM" && -s "$BAI" &&
              -s "$STAR_LOG" && -s "$IDXSTATS" &&
              -s "$STAR_QC_LOG" ]] \
            && samtools quickcheck -v "$BAM" 2>/dev/null; then

            echo "STAR already completed, skipping: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"
            continue
        else
            echo "STAR alignment: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"

            STAR \
                --runThreadN "$THREADS" \
                --genomeDir "$STAR_INDEX" \
                --readFilesIn "$CLEAN_FASTQ" \
                --readFilesCommand zcat \
                --outFileNamePrefix "$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NH HI AS nM MD \
                --seedSearchStartLmaxOverLread 0.5 \
                --outFilterMismatchNmax 2 \
                --outMultimapperOrder Random \
                --outFilterMultimapNmax 20 \
                --outSAMmultNmax 1

            samtools index "$BAM"

            samtools idxstats "$BAM" \
                > "$IDXSTATS"

            cp "$STAR_LOG" "$STAR_QC_LOG"
        fi
    done
}
run_featurecounts () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        MODE_ALIGN_DIR="$ALIGN_DIR/$SIZE_MODE"
        MODE_COUNT_DIR="$COUNT_DIR/$SIZE_MODE"
        mkdir -p "$MODE_COUNT_DIR"

        BAM="$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Aligned.sortedByCoord.out.bam"
        COUNTS="$MODE_COUNT_DIR/${SAMPLE}.${SIZE_MODE}.CDS_counts.txt"
        SUMMARY="${COUNTS}.summary"
        FC_LOG="$LOG_DIR/${SAMPLE}.${SIZE_MODE}.featureCounts.log"

        if [[ -s "$COUNTS" && -s "$SUMMARY" && -s "$FC_LOG" ]] \
            && grep -q "^Geneid" "$COUNTS" \
            && grep -q "^Assigned" "$SUMMARY"; then
            echo "featureCounts already completed, skipping: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"
            continue
        else
            echo "featureCounts: $SAMPLE - $SIZE_MODE" \
                | tee -a "$MASTER_LOG"

            featureCounts \
                -T "$THREADS" \
                -t exon \
                -g gene_id \
                -a "$GTF" \
                -o "$COUNTS" \
                "$BAM" \
                2> "$FC_LOG"
        fi
    done
}

run_report () {
    REPORT_HTML="$MULTIQC_DIR/${PROJECT_NAME}_MultiQC.html"
    REPORT_TSV="$MULTIQC_DIR/${PROJECT_NAME}_QC_report.tsv"
    REPORT_COPY="$REPORT_TABLE_DIR/${PROJECT_NAME}_QC_report.tsv"

    if [[ -s "$REPORT_HTML" && -s "$REPORT_TSV" && -s "$REPORT_COPY" ]]; then
        echo "MultiQC and QC report already exist, skipping report generation." \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Generating MultiQC and final QC report..." | tee -a "$MASTER_LOG"

        multiqc "$BASE_DIR" \
            -o "$MULTIQC_DIR" \
            -n "${PROJECT_NAME}_MultiQC.html" \
            -f

        echo -e "Sample\tRun\tMode\tInput_STAR\tUnique_Reads\tUnique_%\tAvg_Mapped_Length\tToo_Many_Loci_%\tToo_Short_%\tCDS_Assigned\tNoFeatures\tAmbiguity" \
            > "$REPORT_TSV"

        for RUN in "${RUNS[@]}"; do
            SAMPLE="${RUN_MAP[$RUN]}"

            for SIZE_MODE in "${SIZE_MODES[@]}"; do
                STAR_LOG="$ALIGN_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}_Log.final.out"
                FC_SUMMARY="$COUNT_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}.CDS_counts.txt.summary"

                [[ -f "$STAR_LOG" && -f "$FC_SUMMARY" ]] || continue

                INPUT_STAR=$(grep "Number of input reads" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
                UNIQUE_READS=$(grep "Uniquely mapped reads number" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
                UNIQUE_PCT=$(grep "Uniquely mapped reads %" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
                AVG_LEN=$(grep "Average mapped length" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
                TOO_MANY=$(grep "% of reads mapped to too many loci" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
                TOO_SHORT=$(grep "% of reads unmapped: too short" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')

                ASSIGNED=$(grep "^Assigned" "$FC_SUMMARY" | awk '{print $2}')
                NOFEATURES=$(grep "^Unassigned_NoFeatures" "$FC_SUMMARY" | awk '{print $2}')
                AMBIGUITY=$(grep "^Unassigned_Ambiguity" "$FC_SUMMARY" | awk '{print $2}')

                echo -e "${SAMPLE}\t${RUN}\t${SIZE_MODE}\t${INPUT_STAR}\t${UNIQUE_READS}\t${UNIQUE_PCT}\t${AVG_LEN}\t${TOO_MANY}\t${TOO_SHORT}\t${ASSIGNED}\t${NOFEATURES}\t${AMBIGUITY}" \
                    >> "$REPORT_TSV"
            done
        done

        cp "$REPORT_TSV" "$REPORT_COPY"
    fi
}


run_riboseq_qc_figures () {
    echo "Generating Ribo-seq QC figures..." | tee -a "$MASTER_LOG"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        BAM_DIR="$ALIGN_DIR/$SIZE_MODE"
        COUNT_MODE_DIR="$COUNT_DIR/$SIZE_MODE"
        OUTDIR="$FIG_DIR/$SIZE_MODE"

        mkdir -p "$OUTDIR"

        FIGURES=(
            "$OUTDIR/${PROJECT_NAME}_read_length_bam.pdf"
            "$OUTDIR/${PROJECT_NAME}_psite_region.pdf"
            "$OUTDIR/${PROJECT_NAME}_psite_metagene.pdf"
            "$OUTDIR/${PROJECT_NAME}_pca_featurecounts.pdf"
            "$OUTDIR/${PROJECT_NAME}_frame_preference.pdf"
            "$OUTDIR/${PROJECT_NAME}_periodicity.pdf"
        )

        ALL_EXIST=true
        for FIG in "${FIGURES[@]}"; do
            [[ -s "$FIG" ]] || ALL_EXIST=false
        done

        if [[ "$ALL_EXIST" == true ]]; then
            echo "QC figures already exist, skipping: $SIZE_MODE" \
                | tee -a "$MASTER_LOG"
            continue
        else
            echo "Generating QC figures: $SIZE_MODE" \
                | tee -a "$MASTER_LOG"

            Rscript "$SCRIPT_DIR/qc_plots/01_read_length_bam.R" "$BAM_DIR" "$OUTDIR" "$PROJECT_NAME"
            Rscript "$SCRIPT_DIR/qc_plots/02_psite_region.R" "$BAM_DIR" "$GTF" "$OUTDIR" "$PROJECT_NAME" 12
            Rscript "$SCRIPT_DIR/qc_plots/03_psite_metagene.R" "$BAM_DIR" "$GTF" "$OUTDIR" "$PROJECT_NAME" 12
            Rscript "$SCRIPT_DIR/qc_plots/04_pca_featurecounts.R" "$COUNT_MODE_DIR" "$OUTDIR" "$PROJECT_NAME"
            Rscript "$SCRIPT_DIR/qc_plots/05_frame_preference.R" "$BAM_DIR" "$GTF" "$OUTDIR" "$PROJECT_NAME" 12
            Rscript "$SCRIPT_DIR/qc_plots/06_periodicity.R" "$BAM_DIR" "$GTF" "$OUTDIR" "$PROJECT_NAME" 12
        fi
    done

    ALIGNMENT_SUMMARY="$FIG_DIR/${PROJECT_NAME}_alignment_summary.pdf"
    CONTAMINANT_SUMMARY="$FIG_DIR/${PROJECT_NAME}_contaminant_summary.pdf"

    if [[ -s "$ALIGNMENT_SUMMARY" ]]; then
        echo "Alignment summary already exists, skipping."
    else
        Rscript "$SCRIPT_DIR/qc_plots/07_alignment_summary.R" "$STAR_QC_DIR" "$FIG_DIR" "$PROJECT_NAME"
    fi

    if [[ -s "$CONTAMINANT_SUMMARY" ]]; then
        echo "Contaminant summary already exists, skipping."
    else
        Rscript "$SCRIPT_DIR/qc_plots/08_contaminant_summary.R" "$LOG_DIR" "$FIG_DIR" "$PROJECT_NAME"
    fi
}


run_ribotricer () {
    RIBOTRICER_MARKER="$RIBOTRICER_DIR/.ribotricer_completed"

    if [[ -s "$RIBOTRICER_MARKER" ]]; then
        echo "Ribotricer already completed, skipping." \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Running Ribotricer lncRNA-smORF discovery..." \
            | tee -a "$MASTER_LOG"

        if bash "$SCRIPT_DIR/ribotricer/01_run_ribotricer.sh" "$CONFIG_FILE"; then
            touch "$RIBOTRICER_MARKER"
            echo "Ribotricer completed successfully." \
                | tee -a "$MASTER_LOG"
        else
            echo "ERROR: Ribotricer failed; completion marker not created." \
                | tee -a "$MASTER_LOG"
            return 1
        fi
    fi
}


run_final_report () {
    FINAL_REPORT="$REPORT_TABLE_DIR/${PROJECT_NAME}_QC_master_table.tsv"

    if [[ -s "$FINAL_REPORT" ]]; then
        echo "Final QC report already exists, skipping." \
            | tee -a "$MASTER_LOG"
        return 0
    else
        echo "Generating final project report..." \
            | tee -a "$MASTER_LOG"

        Rscript "$SCRIPT_DIR/qc_plots/09_qc_master_table.R" \
            "$BASE_DIR" \
            "$PROJECT_NAME"
    fi
}

if [[ "$MODULE" == "0" || "$MODULE" == "12" ]]; then
    run_sra_conversion
    pause_step
fi

for SAMPLE in "${SAMPLES[@]}"; do
    echo "======================================" | tee -a "$MASTER_LOG"
    echo "Sample: $SAMPLE" | tee -a "$MASTER_LOG"
    echo "======================================" | tee -a "$MASTER_LOG"

    ACTIVE_MODULE="$MODULE"

    if [[ "$CONTINUE_ALL" == true ]]; then
        ACTIVE_MODULE="12"
    fi

    case "$ACTIVE_MODULE" in
        1)
            run_fastqc_raw "$SAMPLE"
            pause_step
            ;;
        2)
            run_cutadapt "$SAMPLE"
            pause_step
            ;;
        3)
            run_qc_trimmed "$SAMPLE"
            pause_step
            ;;
        4)
            run_filter "$SAMPLE"
            pause_step
            ;;
        5)
            run_bowtie "$SAMPLE"
            pause_step
            ;;
        6)
            run_fastqc_clean "$SAMPLE"
            pause_step
            ;;
        7)
            run_star "$SAMPLE"
            pause_step
            ;;
        8)
            run_featurecounts "$SAMPLE"
            pause_step
            ;;
        12)
            run_fastqc_raw "$SAMPLE"
            pause_step

            run_cutadapt "$SAMPLE"
            pause_step

            run_qc_trimmed "$SAMPLE"
            pause_step

            run_filter "$SAMPLE"
            pause_step

            run_bowtie "$SAMPLE"
            pause_step

            run_fastqc_clean "$SAMPLE"
            pause_step

            run_star "$SAMPLE"
            pause_step

            run_featurecounts "$SAMPLE"
            pause_step
            ;;
        0|9|10|11)
            ;;
        *)
            echo "Invalid module."
            exit 1
            ;;
    esac
done

FINAL_MODULE="$MODULE"

if [[ "$CONTINUE_ALL" == true ]]; then
    FINAL_MODULE="12"
fi

if [[ "$FINAL_MODULE" == "9" || "$FINAL_MODULE" == "12" ]]; then
    run_report
    run_riboseq_qc_figures
fi

if [[ "$FINAL_MODULE" == "10" || "$FINAL_MODULE" == "12" ]]; then
    run_ribotricer
fi

if [[ "$FINAL_MODULE" == "11" || "$FINAL_MODULE" == "12" ]]; then
    run_final_report
fi

echo "Finished at: $(date)" | tee -a "$MASTER_LOG"
