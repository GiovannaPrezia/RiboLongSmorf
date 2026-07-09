#!/usr/bin/env python3

import sys
import re
import numpy as np
import pandas as pd


def parse_gtf_attributes(attr):
    result = {}

    for item in str(attr).strip().split(";"):
        item = item.strip()

        if not item or " " not in item:
            continue

        key, value = item.split(" ", 1)
        result[key] = value.replace('"', "").strip()

    return result


def load_gtf_annotation(gtf_file):
    gtf_cols = [
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute"
    ]

    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        names=gtf_cols,
        low_memory=False
    )

    gtf = gtf[gtf["feature"] == "transcript"].copy()

    records = []

    for _, row in gtf.iterrows():
        attrs = parse_gtf_attributes(row["attribute"])

        records.append({
            "transcript_id": attrs.get("transcript_id"),
            "gtf_gene_id": attrs.get("gene_id"),
            "gtf_gene_name": attrs.get("gene_name"),
            "gtf_gene_type": attrs.get("gene_type", attrs.get("gene_biotype")),
            "gtf_transcript_type": attrs.get(
                "transcript_type",
                attrs.get("transcript_biotype")
            )
        })

    annotation = pd.DataFrame(records).drop_duplicates()

    return annotation


def find_column(columns, candidates):
    for col in candidates:
        if col in columns:
            return col
    return None


def extract_transcript_id(df):
    tx_col = find_column(
        df.columns,
        ["transcript_id", "transcript", "transcript_ids", "transcript_name"]
    )

    if tx_col is not None:
        df["transcript_id"] = df[tx_col].astype(str)
        return df

    id_col = find_column(
        df.columns,
        ["ORF_ID", "orf_id", "orf_name", "id"]
    )

    if id_col is None:
        raise ValueError(
            "Could not find transcript_id or ORF ID column. "
            f"Columns found: {list(df.columns)}"
        )

    df["transcript_id"] = df[id_col].astype(str).apply(
        lambda x: re.search(r"(ENST[0-9]+(?:\.[0-9]+)?)", x).group(1)
        if re.search(r"(ENST[0-9]+(?:\.[0-9]+)?)", x)
        else None
    )

    return df


def add_orf_length(df):
    length_col = find_column(
        df.columns,
        [
            "ORF_length_nt",
            "orf_length_nt",
            "length_nt",
            "length",
            "ORF_length",
            "orf_length",
            "orf_len"
        ]
    )

    if length_col is None:
        raise ValueError(
            "Could not identify ORF length column. "
            f"Columns found: {list(df.columns)}"
        )

    df[length_col] = pd.to_numeric(df[length_col], errors="coerce")

    # Ribotricer 'length' is nucleotide length.
    if length_col in ["ORF_length_aa", "orf_length_aa", "length_aa"]:
        df["ORF_length_aa"] = df[length_col].round(0).astype("Int64")
        df["ORF_length_nt"] = df["ORF_length_aa"] * 3
    else:
        df["ORF_length_nt"] = df[length_col].round(0).astype("Int64")
        df["ORF_length_aa"] = (df["ORF_length_nt"] / 3).round(0).astype("Int64")

    return df


def add_clean_annotation(df):
    """
    Keep Ribotricer annotation when available, but add GTF-derived fields
    with explicit names. Then create final clean columns:
    gene_id, gene_name, gene_type, transcript_type
    without _x/_y suffixes.
    """

    if "gtf_gene_id" in df.columns:
        df["gene_id"] = df["gene_id"].where(
            df.get("gene_id").notna(),
            df["gtf_gene_id"]
        ) if "gene_id" in df.columns else df["gtf_gene_id"]

    if "gtf_gene_name" in df.columns:
        df["gene_name"] = df["gene_name"].where(
            df.get("gene_name").notna(),
            df["gtf_gene_name"]
        ) if "gene_name" in df.columns else df["gtf_gene_name"]

    if "gtf_gene_type" in df.columns:
        df["gene_type"] = df["gene_type"].where(
            df.get("gene_type").notna(),
            df["gtf_gene_type"]
        ) if "gene_type" in df.columns else df["gtf_gene_type"]

    if "gtf_transcript_type" in df.columns:
        df["transcript_type"] = df["transcript_type"].where(
            df.get("transcript_type").notna(),
            df["gtf_transcript_type"]
        ) if "transcript_type" in df.columns else df["gtf_transcript_type"]

    return df


def add_ranking(df):
    phase_col = find_column(
        df.columns,
        [
            "phase_score",
            "phaseScore",
            "periodicity_score",
            "phase_score_valid_codons"
        ]
    )

    read_col = find_column(
        df.columns,
        [
            "read_count",
            "reads",
            "count",
            "RPF_count",
            "valid_codons"
        ]
    )

    if phase_col is None:
        df["phase_component"] = 0
    else:
        df["phase_component"] = pd.to_numeric(
            df[phase_col],
            errors="coerce"
        ).fillna(0)

    if read_col is None:
        df["read_component"] = 0
    else:
        df["read_component"] = np.log10(
            pd.to_numeric(df[read_col], errors="coerce").fillna(0) + 1
        )

    df["ranking_score"] = (
        df["phase_component"] +
        df["read_component"]
    )

    return df.sort_values("ranking_score", ascending=False)


def main():
    if len(sys.argv) != 5:
        sys.exit(
            "Usage: python 02_process_ribotricer_results.py "
            "<ribotricer_orfs.tsv> <gencode.gtf> <output_prefix> <sample_name>"
        )

    ribotricer_file = sys.argv[1]
    gtf_file = sys.argv[2]
    output_prefix = sys.argv[3]
    sample_name = sys.argv[4]

    print(f"Loading Ribotricer results: {ribotricer_file}")

    df = pd.read_csv(ribotricer_file, sep="\t")
    df["sample"] = sample_name

    df = extract_transcript_id(df)
    df = add_orf_length(df)

    smorfs = df[
        (df["ORF_length_aa"] >= 20) &
        (df["ORF_length_aa"] <= 150)
    ].copy()

    print(f"Total ORFs       : {len(df):,}")
    print(f"smORFs 20-150 aa : {len(smorfs):,}")

    smorfs_file = f"{output_prefix}_smorfs_20_150aa.tsv"
    smorfs.to_csv(smorfs_file, sep="\t", index=False)

    print(f"Saved: {smorfs_file}")

    annotation = load_gtf_annotation(gtf_file)

    annotated = smorfs.merge(
        annotation,
        on="transcript_id",
        how="left"
    )

    annotated = add_clean_annotation(annotated)

    annotated_file = f"{output_prefix}_smorfs_annotated.tsv"
    annotated.to_csv(annotated_file, sep="\t", index=False)

    print(f"Saved: {annotated_file}")

    lnc_types = {
        "lncRNA",
        "lincRNA",
        "antisense",
        "processed_transcript",
        "sense_intronic",
        "sense_overlapping",
        "macro_lncRNA",
        "3prime_overlapping_ncRNA",
        "non_coding"
    }

    if "gene_type" not in annotated.columns:
        raise ValueError(
            "gene_type column not found after annotation. "
            f"Columns found: {list(annotated.columns)}"
        )

    print("Using biotype column for lncRNA filtering: gene_type")

    lncrna_smorfs = annotated[
        annotated["gene_type"].isin(lnc_types)
    ].copy()

    lncrna_file = f"{output_prefix}_lncrna_smorfs.tsv"
    lncrna_smorfs.to_csv(lncrna_file, sep="\t", index=False)

    print(f"lncRNA-smORFs    : {len(lncrna_smorfs):,}")
    print(f"Saved: {lncrna_file}")

    ranked = add_ranking(lncrna_smorfs)

    ranked_file = f"{output_prefix}_ranked_lncrna_smorfs.tsv"
    ranked.to_csv(ranked_file, sep="\t", index=False)

    print(f"Saved: {ranked_file}")


if __name__ == "__main__":
    main()
