#!/usr/bin/env python3

from __future__ import annotations

import bisect
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


MIN_AA = 20
MAX_AA = 150

LNC_GENE_TYPES = {
    "lncrna",
    "lincrna",
    "antisense",
    "processed_transcript",
    "sense_intronic",
    "sense_overlapping",
    "macro_lncrna",
    "3prime_overlapping_ncrna",
    "non_coding",
}

COORD_RE = re.compile(r"(\d+)-(\d+)")
TRANSCRIPT_RE = re.compile(r"(ENST[0-9]+(?:\.[0-9]+)?)")


def clean_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return (
        str(value)
        .strip()
        .lstrip("\ufeff")
        .strip()
        .strip('"')
        .strip("'")
        .strip()
    )


def normalize_type(value: Any) -> str:
    return (
        clean_text(value)
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )


def strip_version(identifier: Any) -> str:
    value = clean_text(identifier)
    return value.split(".", 1)[0] if value else ""


def parse_gtf_attributes(attribute: Any) -> dict[str, str]:
    result: dict[str, str] = {}

    for item in clean_text(attribute).split(";"):
        item = item.strip()

        if not item:
            continue

        if " " in item:
            key, value = item.split(" ", 1)
        elif "=" in item:
            key, value = item.split("=", 1)
        else:
            continue

        result[key.strip()] = value.strip().strip('"')

    return result


def find_column(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    available = {clean_text(column): column for column in columns}

    for candidate in candidates:
        if candidate in available:
            return available[candidate]

    lower_available = {
        clean_text(column).lower(): column
        for column in columns
    }

    for candidate in candidates:
        match = lower_available.get(candidate.lower())
        if match is not None:
            return match

    return None


def clean_dataframe_headers(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [clean_text(column) for column in df.columns]
    return df


def load_gtf(gtf_file: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    gtf_columns = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        names=gtf_columns,
        dtype={
            "seqname": "string",
            "feature": "string",
            "start": "Int64",
            "end": "Int64",
            "strand": "string",
            "attribute": "string",
        },
        low_memory=False,
    )

    gene_records: list[dict[str, Any]] = []
    transcript_records: list[dict[str, Any]] = []
    cds_records: list[tuple[str, str, int, int, str, str]] = []

    for row in gtf.itertuples(index=False):
        attrs = parse_gtf_attributes(row.attribute)

        gene_id = attrs.get("gene_id", "")
        gene_id_base = strip_version(gene_id)
        gene_name = attrs.get("gene_name", gene_id_base)
        gene_type = (
            attrs.get("gene_type")
            or attrs.get("gene_biotype")
            or ""
        )

        if row.feature == "gene" and gene_id_base:
            gene_records.append(
                {
                    "gene_id_base": gene_id_base,
                    "gtf_gene_id": gene_id,
                    "gtf_gene_name": gene_name,
                    "gtf_gene_type": gene_type,
                    "gtf_gene_chrom": clean_text(row.seqname),
                    "gtf_gene_strand": clean_text(row.strand),
                }
            )

        elif row.feature == "transcript":
            transcript_id = attrs.get("transcript_id", "")
            if transcript_id:
                transcript_records.append(
                    {
                        "transcript_id_base": strip_version(transcript_id),
                        "gtf_transcript_id": transcript_id,
                        "transcript_gene_id_base": gene_id_base,
                        "gtf_transcript_type": (
                            attrs.get("transcript_type")
                            or attrs.get("transcript_biotype")
                            or ""
                        ),
                    }
                )

        elif (
            row.feature == "CDS"
            and normalize_type(gene_type) == "protein_coding"
        ):
            cds_records.append(
                (
                    clean_text(row.seqname),
                    clean_text(row.strand),
                    int(row.start),
                    int(row.end),
                    gene_id_base,
                    gene_name,
                )
            )

    genes = (
        pd.DataFrame(gene_records)
        .drop_duplicates(subset=["gene_id_base"])
    )

    transcripts = (
        pd.DataFrame(transcript_records)
        .drop_duplicates(subset=["transcript_id_base"])
    )

    cds_index = build_interval_index(cds_records)

    return genes, transcripts, cds_index


def build_interval_index(
    records: Iterable[tuple[str, str, int, int, str, str]]
) -> dict[tuple[str, str], dict[str, Any]]:
    grouped: dict[
        tuple[str, str],
        list[tuple[int, int, str, str]]
    ] = defaultdict(list)

    for chrom, strand, start, end, gene_id, gene_name in records:
        grouped[(chrom, strand)].append(
            (int(start), int(end), gene_id, gene_name)
        )

    index: dict[tuple[str, str], dict[str, Any]] = {}

    for key, intervals in grouped.items():
        intervals.sort(key=lambda item: (item[0], item[1]))

        starts = [item[0] for item in intervals]
        prefix_max_end: list[int] = []
        current_max = -1

        for _, end, _, _ in intervals:
            current_max = max(current_max, end)
            prefix_max_end.append(current_max)

        index[key] = {
            "intervals": intervals,
            "starts": starts,
            "prefix_max_end": prefix_max_end,
        }

    return index


def query_overlaps(
    interval_index: dict,
    chrom: str,
    strand: str,
    query_start: int,
    query_end: int,
) -> list[tuple[int, int, str, str]]:
    data = interval_index.get((chrom, strand))

    if data is None:
        return []

    intervals = data["intervals"]
    starts = data["starts"]
    prefix_max_end = data["prefix_max_end"]

    index = bisect.bisect_right(starts, query_end) - 1
    overlaps: list[tuple[int, int, str, str]] = []

    while index >= 0:
        if prefix_max_end[index] < query_start:
            break

        start, end, gene_id, gene_name = intervals[index]

        if end >= query_start and start <= query_end:
            overlaps.append((start, end, gene_id, gene_name))

        index -= 1

    return overlaps


def union_length(intervals: Iterable[tuple[int, int]]) -> int:
    sorted_intervals = sorted(intervals)

    if not sorted_intervals:
        return 0

    total = 0
    current_start, current_end = sorted_intervals[0]

    for start, end in sorted_intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end

    total += current_end - current_start + 1

    return total


def parse_coordinate_blocks(value: Any) -> tuple[tuple[int, int], ...]:
    pairs = COORD_RE.findall(clean_text(value))

    if not pairs:
        return tuple()

    blocks = sorted(
        (min(int(start), int(end)), max(int(start), int(end)))
        for start, end in pairs
    )

    merged: list[list[int]] = []

    for start, end in blocks:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    return tuple((start, end) for start, end in merged)


def blocks_to_string(blocks: Iterable[tuple[int, int]]) -> str:
    return ",".join(f"{start}-{end}" for start, end in blocks)


def blocks_length(blocks: Iterable[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in blocks)


def extract_transcript_id(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    transcript_column = find_column(
        df.columns,
        [
            "transcript_id",
            "transcript",
            "transcript_ids",
            "transcript_name",
        ],
    )

    if transcript_column is not None:
        df["transcript_id"] = df[transcript_column].map(clean_text)
    else:
        id_column = find_column(
            df.columns,
            ["ORF_ID", "orf_id", "orf_name", "id"],
        )

        if id_column is None:
            raise ValueError(
                "Could not find transcript_id or ORF ID column. "
                f"Columns found: {list(df.columns)}"
            )

        df["transcript_id"] = (
            df[id_column]
            .astype(str)
            .str.extract(TRANSCRIPT_RE, expand=False)
        )

    df["transcript_id_base"] = df["transcript_id"].map(strip_version)

    return df


def add_orf_length(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    aa_column = find_column(
        df.columns,
        ["ORF_length_aa", "orf_length_aa", "length_aa"],
    )

    nt_column = find_column(
        df.columns,
        [
            "ORF_length_nt",
            "orf_length_nt",
            "length_nt",
            "length",
            "ORF_length",
            "orf_length",
            "orf_len",
        ],
    )

    if aa_column is None and nt_column is None:
        raise ValueError(
            "Could not identify ORF length column. "
            f"Columns found: {list(df.columns)}"
        )

    if aa_column is not None:
        aa = pd.to_numeric(df[aa_column], errors="coerce")
        df["ORF_length_aa"] = aa.round().astype("Int64")
        df["ORF_length_nt"] = (df["ORF_length_aa"] * 3).astype("Int64")
    else:
        nt = pd.to_numeric(df[nt_column], errors="coerce")
        df["ORF_length_nt"] = nt.round().astype("Int64")
        df["ORF_length_aa"] = np.floor(
            df["ORF_length_nt"] / 3
        ).astype("Int64")

    return df


def add_genomic_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    chrom_column = find_column(
        df.columns,
        ["chrom", "chr", "chromosome", "seqname", "contig"],
    )
    strand_column = find_column(
        df.columns,
        ["strand", "orf_strand"],
    )
    coordinate_column = find_column(
        df.columns,
        [
            "coordinate",
            "coordinates",
            "genomic_coordinates",
            "ORF_coordinate",
            "orf_coordinate",
            "ORF_coordinates",
            "orf_coordinates",
        ],
    )

    missing = []

    if chrom_column is None:
        missing.append("chrom")
    if strand_column is None:
        missing.append("strand")
    if coordinate_column is None:
        missing.append("coordinate/genomic_coordinates")

    if missing:
        raise ValueError(
            "Columns required to construct genomic_orf_id were not found: "
            + ", ".join(missing)
            + f". Columns found: {list(df.columns)}"
        )

    df["chrom"] = df[chrom_column].map(clean_text)
    df["strand"] = df[strand_column].map(clean_text)
    df["coordinate_original"] = df[coordinate_column].map(clean_text)
    df["_blocks"] = df["coordinate_original"].map(parse_coordinate_blocks)

    invalid = df["_blocks"].map(len).eq(0)

    if invalid.any():
        examples = (
            df.loc[invalid, "coordinate_original"]
            .drop_duplicates()
            .head(10)
            .tolist()
        )
        raise ValueError(
            f"{int(invalid.sum())} rows have invalid coordinates. "
            f"Examples: {examples}"
        )

    df["coordinate"] = df["_blocks"].map(blocks_to_string)
    df["coordinate_length_nt"] = df["_blocks"].map(blocks_length)

    df["genomic_orf_id"] = (
        df["chrom"]
        + "|"
        + df["strand"]
        + "|"
        + df["coordinate"]
    )

    return df


def annotate_from_gtf(
    df: pd.DataFrame,
    genes: pd.DataFrame,
    transcripts: pd.DataFrame,
) -> pd.DataFrame:
    df = df.copy()

    original_columns = set(df.columns)

    for column in [
        "gene_id",
        "gene_name",
        "gene_type",
        "transcript_type",
    ]:
        if column in df.columns:
            df[f"ribotricer_{column}"] = df[column]

    if not transcripts.empty:
        df = df.merge(
            transcripts,
            on="transcript_id_base",
            how="left",
            validate="many_to_one",
        )

    if "gene_id" in original_columns:
        df["gene_id_base"] = df["gene_id"].map(strip_version)
    elif "ribotricer_gene_id" in df.columns:
        df["gene_id_base"] = df["ribotricer_gene_id"].map(strip_version)
    else:
        df["gene_id_base"] = ""

    if "transcript_gene_id_base" in df.columns:
        missing_gene = df["gene_id_base"].eq("")
        df.loc[missing_gene, "gene_id_base"] = df.loc[
            missing_gene,
            "transcript_gene_id_base",
        ]

    if not genes.empty:
        df = df.merge(
            genes,
            on="gene_id_base",
            how="left",
            validate="many_to_one",
        )

    df["gene_id"] = df.get(
        "gtf_gene_id",
        pd.Series("", index=df.index),
    ).fillna("")

    if "ribotricer_gene_id" in df.columns:
        missing = df["gene_id"].eq("")
        df.loc[missing, "gene_id"] = df.loc[
            missing,
            "ribotricer_gene_id",
        ].map(clean_text)

    df["gene_name"] = df.get(
        "gtf_gene_name",
        pd.Series("", index=df.index),
    ).fillna("")

    if "ribotricer_gene_name" in df.columns:
        missing = df["gene_name"].eq("")
        df.loc[missing, "gene_name"] = df.loc[
            missing,
            "ribotricer_gene_name",
        ].map(clean_text)

    df["gene_type"] = df.get(
        "gtf_gene_type",
        pd.Series("", index=df.index),
    ).fillna("")

    if "ribotricer_gene_type" in df.columns:
        missing = df["gene_type"].eq("")
        df.loc[missing, "gene_type"] = df.loc[
            missing,
            "ribotricer_gene_type",
        ].map(clean_text)

    df["transcript_type"] = df.get(
        "gtf_transcript_type",
        pd.Series("", index=df.index),
    ).fillna("")

    if "ribotricer_transcript_type" in df.columns:
        missing = df["transcript_type"].eq("")
        df.loc[missing, "transcript_type"] = df.loc[
            missing,
            "ribotricer_transcript_type",
        ].map(clean_text)

    df["gene_type_normalized"] = df["gene_type"].map(normalize_type)

    return df


def filter_translating_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    status_column = find_column(
        df.columns,
        ["status", "translation_status", "translating"],
    )

    if status_column is None:
        return df.copy(), 0

    status = df[status_column].map(normalize_type)

    accepted = status.isin(
        {"translating", "true", "1", "yes", "pass", "passed"}
    )

    removed = int((~accepted).sum())

    return df.loc[accepted].copy(), removed


def get_numeric_column(
    df: pd.DataFrame,
    candidates: Iterable[str],
    output_name: str,
    default: float = 0.0,
) -> pd.DataFrame:
    df = df.copy()

    column = find_column(df.columns, candidates)

    if column is None:
        df[output_name] = default
    else:
        df[output_name] = pd.to_numeric(
            df[column],
            errors="coerce",
        ).fillna(default)

    return df


def unique_join(values: Iterable[Any]) -> str:
    cleaned = {
        clean_text(value)
        for value in values
        if clean_text(value)
    }

    return ";".join(sorted(cleaned))


def first_nonempty(values: Iterable[Any]) -> str:
    for value in values:
        cleaned = clean_text(value)
        if cleaned:
            return cleaned
    return ""


def collapse_genomic_orfs(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df = get_numeric_column(
        df,
        ["read_count", "reads", "count", "RPF_count"],
        "read_count_numeric",
    )

    df = get_numeric_column(
        df,
        [
            "phase_score",
            "phaseScore",
            "periodicity_score",
            "phase_score_valid_codons",
        ],
        "phase_score_numeric",
    )

    df = get_numeric_column(
        df,
        ["valid_codons"],
        "valid_codons_numeric",
    )

    df = get_numeric_column(
        df,
        ["read_density"],
        "read_density_numeric",
    )

    if "ORF_ID" not in df.columns:
        id_column = find_column(
            df.columns,
            ["orf_id", "orf_name", "id"],
        )
        if id_column is not None:
            df["ORF_ID"] = df[id_column].map(clean_text)
        else:
            df["ORF_ID"] = ""

    group_columns = ["genomic_orf_id", "chrom", "strand", "coordinate"]

    collapsed = (
        df.groupby(group_columns, dropna=False, sort=False)
        .agg(
            sample=("sample", first_nonempty),
            ORF_IDs=("ORF_ID", unique_join),
            transcript_ids=("transcript_id", unique_join),
            gene_ids=("gene_id", unique_join),
            gene_names=("gene_name", unique_join),
            gene_types=("gene_type", unique_join),
            transcript_types=("transcript_type", unique_join),
            ORF_types=(
                find_column(df.columns, ["ORF_type", "orf_type"])
                or "ORF_ID",
                unique_join,
            ),
            start_codons=(
                find_column(df.columns, ["start_codon", "start_codons"])
                or "ORF_ID",
                unique_join,
            ),
            ORF_length_nt_min=("ORF_length_nt", "min"),
            ORF_length_nt_max=("ORF_length_nt", "max"),
            ORF_length_aa_min=("ORF_length_aa", "min"),
            ORF_length_aa_max=("ORF_length_aa", "max"),
            coordinate_length_nt=("coordinate_length_nt", "first"),
            read_count=("read_count_numeric", "max"),
            phase_score=("phase_score_numeric", "max"),
            valid_codons=("valid_codons_numeric", "max"),
            read_density=("read_density_numeric", "max"),
            source_rows=("genomic_orf_id", "size"),
            n_ORF_IDs=("ORF_ID", lambda values: len({
                clean_text(value)
                for value in values
                if clean_text(value)
            })),
            n_transcripts=("transcript_id", lambda values: len({
                clean_text(value)
                for value in values
                if clean_text(value)
            })),
            n_gene_ids=("gene_id", lambda values: len({
                clean_text(value)
                for value in values
                if clean_text(value)
            })),
            n_gene_types=("gene_type", lambda values: len({
                normalize_type(value)
                for value in values
                if normalize_type(value)
            })),
        )
        .reset_index()
    )

    collapsed["length_conflict"] = (
        collapsed["ORF_length_nt_min"]
        .ne(collapsed["ORF_length_nt_max"])
        | collapsed["ORF_length_aa_min"]
        .ne(collapsed["ORF_length_aa_max"])
    )

    collapsed["ORF_length_nt"] = collapsed["ORF_length_nt_max"].astype("Int64")
    collapsed["ORF_length_aa"] = collapsed["ORF_length_aa_max"].astype("Int64")

    collapsed["gene_type_set_normalized"] = collapsed["gene_types"].map(
        lambda value: ";".join(
            sorted({
                normalize_type(item)
                for item in clean_text(value).split(";")
                if normalize_type(item)
            })
        )
    )

    collapsed["all_genes_are_lncRNA"] = collapsed[
        "gene_type_set_normalized"
    ].map(
        lambda value: (
            bool(value)
            and set(value.split(";")).issubset(LNC_GENE_TYPES)
        )
    )

    collapsed["mixed_gene_types"] = collapsed[
        "gene_type_set_normalized"
    ].map(
        lambda value: (
            bool(value)
            and bool(set(value.split(";")).intersection(LNC_GENE_TYPES))
            and not set(value.split(";")).issubset(LNC_GENE_TYPES)
        )
    )

    return collapsed


def annotate_cds_overlap(
    df: pd.DataFrame,
    cds_index: dict,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []

    for row in df.itertuples(index=False):
        blocks = parse_coordinate_blocks(row.coordinate)

        same_strand_intersections: list[tuple[int, int]] = []
        opposite_strand_intersections: list[tuple[int, int]] = []

        same_strand_genes: set[str] = set()
        same_strand_gene_ids: set[str] = set()
        opposite_strand_genes: set[str] = set()
        opposite_strand_gene_ids: set[str] = set()

        opposite_strand = "-" if row.strand == "+" else "+"

        for block_start, block_end in blocks:
            for start, end, gene_id, gene_name in query_overlaps(
                cds_index,
                row.chrom,
                row.strand,
                block_start,
                block_end,
            ):
                overlap_start = max(block_start, start)
                overlap_end = min(block_end, end)

                if overlap_start <= overlap_end:
                    same_strand_intersections.append(
                        (overlap_start, overlap_end)
                    )
                    same_strand_gene_ids.add(gene_id)
                    same_strand_genes.add(gene_name)

            for start, end, gene_id, gene_name in query_overlaps(
                cds_index,
                row.chrom,
                opposite_strand,
                block_start,
                block_end,
            ):
                overlap_start = max(block_start, start)
                overlap_end = min(block_end, end)

                if overlap_start <= overlap_end:
                    opposite_strand_intersections.append(
                        (overlap_start, overlap_end)
                    )
                    opposite_strand_gene_ids.add(gene_id)
                    opposite_strand_genes.add(gene_name)

        same_nt = union_length(same_strand_intersections)
        opposite_nt = union_length(opposite_strand_intersections)
        orf_length = int(row.ORF_length_nt or 0)

        if same_nt == 0:
            overlap_class = "NO_SAME_STRAND_PC_CDS"
        elif orf_length > 0 and same_nt >= orf_length:
            overlap_class = "FULL_SAME_STRAND_PC_CDS"
        else:
            overlap_class = "PARTIAL_SAME_STRAND_PC_CDS"

        records.append(
            {
                "genomic_orf_id": row.genomic_orf_id,
                "same_strand_pc_cds_overlap_nt": same_nt,
                "same_strand_pc_cds_overlap_fraction": (
                    same_nt / orf_length
                    if orf_length > 0
                    else np.nan
                ),
                "same_strand_pc_cds_overlap_class": overlap_class,
                "same_strand_pc_gene_ids": ";".join(
                    sorted(same_strand_gene_ids)
                ),
                "same_strand_pc_genes": ";".join(
                    sorted(same_strand_genes)
                ),
                "opposite_strand_pc_cds_overlap_nt": opposite_nt,
                "opposite_strand_pc_cds_overlap_fraction": (
                    opposite_nt / orf_length
                    if orf_length > 0
                    else np.nan
                ),
                "opposite_strand_pc_gene_ids": ";".join(
                    sorted(opposite_strand_gene_ids)
                ),
                "opposite_strand_pc_genes": ";".join(
                    sorted(opposite_strand_genes)
                ),
            }
        )

    overlap_df = pd.DataFrame(records)

    return df.merge(
        overlap_df,
        on="genomic_orf_id",
        how="left",
        validate="one_to_one",
    )


def add_ranking(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["phase_component"] = pd.to_numeric(
        df["phase_score"],
        errors="coerce",
    ).fillna(0)

    df["read_component"] = np.log10(
        pd.to_numeric(
            df["read_count"],
            errors="coerce",
        ).fillna(0)
        + 1
    )

    df["ranking_score"] = (
        df["phase_component"]
        + df["read_component"]
    )

    return df.sort_values(
        [
            "ranking_score",
            "phase_score",
            "read_count",
            "genomic_orf_id",
        ],
        ascending=[False, False, False, True],
    )


def save_tsv(df: pd.DataFrame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"Saved: {path}")


def main() -> None:
    if len(sys.argv) != 5:
        sys.exit(
            "Usage: python 02_process_ribotricer_results.py "
            "<ribotricer_translating_ORFs.tsv> "
            "<gencode.gtf> "
            "<output_prefix> "
            "<sample_name>"
        )

    ribotricer_file = sys.argv[1]
    gtf_file = sys.argv[2]
    output_prefix = sys.argv[3]
    sample_name = sys.argv[4]

    print("============================================================")
    print("RiboLongSmORF — processing one Ribotricer sample")
    print(f"Input : {ribotricer_file}")
    print(f"Sample: {sample_name}")
    print("============================================================")

    raw = pd.read_csv(
        ribotricer_file,
        sep="\t",
        low_memory=False,
    )
    raw = clean_dataframe_headers(raw)
    raw["sample"] = sample_name

    total_input_rows = len(raw)

    translating, removed_nontranslating = filter_translating_rows(raw)

    translating = extract_transcript_id(translating)
    translating = add_orf_length(translating)
    translating = add_genomic_coordinates(translating)

    genes, transcripts, cds_index = load_gtf(gtf_file)

    annotated_rows = annotate_from_gtf(
        translating,
        genes,
        transcripts,
    )

    collapsed = collapse_genomic_orfs(annotated_rows)

    smorfs = collapsed[
        collapsed["ORF_length_aa"].between(
            MIN_AA,
            MAX_AA,
            inclusive="both",
        )
        & ~collapsed["length_conflict"]
    ].copy()

    lncrna_smorfs = smorfs[
        smorfs["all_genes_are_lncRNA"]
        & ~smorfs["mixed_gene_types"]
    ].copy()

    lncrna_smorfs = annotate_cds_overlap(
        lncrna_smorfs,
        cds_index,
    )

    strict = lncrna_smorfs[
        lncrna_smorfs[
            "same_strand_pc_cds_overlap_nt"
        ].eq(0)
    ].copy()

    excluded_pc_cds = lncrna_smorfs[
        lncrna_smorfs[
            "same_strand_pc_cds_overlap_nt"
        ].gt(0)
    ].copy()

    ranked = add_ranking(strict)

    duplicate_audit = collapsed[
        collapsed["source_rows"].gt(1)
        | collapsed["n_ORF_IDs"].gt(1)
        | collapsed["n_transcripts"].gt(1)
        | collapsed["length_conflict"]
        | collapsed["n_gene_ids"].gt(1)
        | collapsed["n_gene_types"].gt(1)
    ].copy()

    gene_type_ambiguity = smorfs[
        ~smorfs["all_genes_are_lncRNA"]
        | smorfs["mixed_gene_types"]
    ].copy()

    save_tsv(
        annotated_rows,
        f"{output_prefix}_translated_annotated_rows.tsv",
    )
    save_tsv(
        collapsed,
        f"{output_prefix}_translated_genomic_orfs.tsv",
    )
    save_tsv(
        smorfs,
        f"{output_prefix}_smorfs_20_150aa.tsv",
    )
    save_tsv(
        lncrna_smorfs,
        f"{output_prefix}_lncrna_smorfs_before_cds_filter.tsv",
    )
    save_tsv(
        strict,
        f"{output_prefix}_lncrna_smorfs_strict.tsv",
    )
    save_tsv(
        excluded_pc_cds,
        f"{output_prefix}_lncrna_smorfs_excluded_pc_cds.tsv",
    )
    save_tsv(
        ranked,
        f"{output_prefix}_ranked_lncrna_smorfs.tsv",
    )
    save_tsv(
        duplicate_audit,
        f"{output_prefix}_duplicate_audit.tsv",
    )
    save_tsv(
        gene_type_ambiguity,
        f"{output_prefix}_gene_type_ambiguity.tsv",
    )

    summary = pd.DataFrame(
        [
            {
                "stage": "input_rows",
                "n": total_input_rows,
                "unit": "rows",
            },
            {
                "stage": "translating_rows",
                "n": len(translating),
                "unit": "rows",
            },
            {
                "stage": "removed_nontranslating_rows",
                "n": removed_nontranslating,
                "unit": "rows",
            },
            {
                "stage": "unique_genomic_orfs",
                "n": len(collapsed),
                "unit": "genomic_orf_id",
            },
            {
                "stage": "smorfs_20_150aa",
                "n": len(smorfs),
                "unit": "genomic_orf_id",
            },
            {
                "stage": "lncrna_smorfs_gene_level",
                "n": len(lncrna_smorfs),
                "unit": "genomic_orf_id",
            },
            {
                "stage": "excluded_same_strand_pc_cds",
                "n": len(excluded_pc_cds),
                "unit": "genomic_orf_id",
            },
            {
                "stage": "strict_lncrna_smorfs",
                "n": len(strict),
                "unit": "genomic_orf_id",
            },
            {
                "stage": "duplicate_or_alias_audit",
                "n": len(duplicate_audit),
                "unit": "genomic_orf_id",
            },
            {
                "stage": "gene_type_ambiguity",
                "n": len(gene_type_ambiguity),
                "unit": "genomic_orf_id",
            },
        ]
    )

    save_tsv(
        summary,
        f"{output_prefix}_processing_summary.tsv",
    )

    print("============================================================")
    print(f"Input rows                         : {total_input_rows:,}")
    print(f"Translating rows                   : {len(translating):,}")
    print(f"Unique genomic ORFs                : {len(collapsed):,}")
    print(f"smORFs {MIN_AA}-{MAX_AA} aa                 : {len(smorfs):,}")
    print(f"lncRNA-smORFs, gene-level filter   : {len(lncrna_smorfs):,}")
    print(f"Excluded by same-strand PC CDS     : {len(excluded_pc_cds):,}")
    print(f"Strict lncRNA-smORFs               : {len(strict):,}")
    print("============================================================")


if __name__ == "__main__":
    main()
