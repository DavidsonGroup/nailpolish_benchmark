#!/usr/bin/env python3
"""
Summarise alignment statistics from a SAM file, one row per (read, group) pair.

USAGE
-----
    summarise_sam.py <input.sam>

INPUT READ FORMAT
-----------------
Each read must carry a GG:Z: tag (type Z, string) listing one or more
comma-separated group names, e.g.:

    GG:Z:cellA,cellB,cellC

The script emits one output row per group the read belongs to.
Reads that lack the GG tag are silently skipped.
Secondary and supplementary alignments are always skipped.

OUTPUT DATAFRAME
----------------
A DataFrame indexed by (query_name, group) with raw alignment columns:
    mapped, mapped_quality, length, alignment_score, ambiguous_bases,
    M, I, D, N, S, H, P, =, X, B, edit_dist
"""

import argparse

import numpy as np
import pandas as pd
import pysam

CIGAR_COLUMNS = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B", "edit_dist"]

CHUNK_SIZE = 2**20


def get_tag_or_na(read, tag):
    if read.has_tag(tag):
        return read.get_tag(tag)
    return np.nan


def read_sam(in_file):
    """Yield one record dict per (read, group) pair."""
    samfile = pysam.AlignmentFile(in_file, "r")
    for read in samfile.fetch():
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("GG"):
            raise ValueError(f"Read {read.query_name} lacks required GG tag, {read}")

        groups = read.get_tag("GG").split(",")
        mapped = read.is_mapped
        cigar_data = dict(zip(CIGAR_COLUMNS, read.get_cigar_stats()[0]))

        base = {
            "query_name": read.query_name,
            "mapped": mapped,
            "mapped_quality": read.mapping_quality,
            "length": read.query_length,
            "chrom": read.reference_name if mapped else None,
            "pos": read.reference_start if mapped else None,
            "alignment_score": read.get_tag("AS") if mapped else None,
            "ambiguous_bases": read.get_tag("nn") if mapped else None,
            **cigar_data,
        }

        for group in groups:
            yield {**base, "group": group}


def process_chunk(records):
    df = pd.DataFrame.from_records(records)
    df = df.astype({"ambiguous_bases": "Int64"})
    df = df.set_index(["query_name", "group"])
    return df


def read_to_df(in_file):
    chunks = []
    chunk = []

    for record in read_sam(in_file):
        chunk.append(record)

        if len(chunk) >= CHUNK_SIZE:
            chunks.append(process_chunk(chunk))
            chunk.clear()

    if chunk:
        chunks.append(process_chunk(chunk))

    return pd.concat(chunks) if chunks else pd.DataFrame()

def _summarise_group(group_df, name):
    parameters = [
        "pct_equality_over_mapped_bases",
        "pct_equality_over_all_bases",
        "error_rate",
        "bases_mapped",
        "length",
        "pct_mapped_over_total",
        "alignment_score",
    ]
    results = {
        "group": name,
        "count": len(group_df),
        "pct_of_reads_mapping": group_df["mapped"].sum() / len(group_df) * 100,
    }
    for param in parameters:
        results[param + "_mean"] = group_df[param].mean()
        results[param + "_median"] = group_df[param].median()
    return pd.Series(results)


def summarise_df(df):
    df = df.copy()

    # NM = #mismatches + #I + #D + #ambiguous_bases (minimap2 convention)
    df["mismatches"] = df["edit_dist"] - df["I"] - df["D"]
    df["indel_rate"] = (df["I"] + df["D"]) / df["length"]
    df["bases_mapped"] = df["M"] + df["I"] + df["="] + df["X"]
    df["bases_mapped"] = df["bases_mapped"].replace(0, np.nan)
    df["error_rate"] = df["edit_dist"] / df["bases_mapped"]

    df["pct_equality_over_mapped_bases"] = df["="] / df["bases_mapped"] * 100
    df["pct_equality_over_all_bases"] = df["="] / df["length"] * 100
    df["pct_mapped_over_total"] = df["bases_mapped"] / df["length"]

    rows = [
        _summarise_group(group_df, name)
        for name, group_df in df.groupby(level="group")
    ]
    return pd.DataFrame(rows).set_index("group")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sam", help="Input SAM file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("--parquet", required=True, help="Output parquet file")
    args = parser.parse_args()

    df = read_to_df(args.sam)
    df.to_parquet(args.parquet, index=True)
    print(f"Saved to {args.parquet}")

    summary = summarise_df(df)
    summary.to_csv(args.output, index=True)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
