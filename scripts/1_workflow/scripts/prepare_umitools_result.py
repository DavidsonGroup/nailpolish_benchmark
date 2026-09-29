#!/usr/bin/env python3
"""
Prepare a umi_tools-processed BAM for summarise_bam.py by adding GG:Z group tags.

Reads a BAM produced by umi_tools dedup and a TSV produced by umi_tools group,
and writes a new BAM where each read carries a GG:Z tag:
  - simplex    (final_umi_count == 1)
  - duplicate  (final_umi_count > 1)

Reads not present in the TSV are omitted from the output.
"""

import argparse

import pandas as pd
import pysam


def prepare_umitools_bam(input_path, tsv_path, output_path):
    df = (
        pd.read_csv(tsv_path, sep="\t", usecols=["read_id", "final_umi_count"])
        .drop_duplicates(subset="read_id", keep="first")
        .set_index("read_id")
    )

    mode_r = "r" if input_path.endswith(".sam") else "rb"

    with pysam.AlignmentFile(input_path, mode_r) as infile, \
         pysam.AlignmentFile(output_path, "w", header=infile.header) as outfile:

        for read in infile:
            if read.is_secondary or read.is_supplementary:
                continue

            query_name = read.query_name
            if query_name not in df.index:
                continue

            count = df.at[query_name, "final_umi_count"]
            tag = "simplex" if count == 1 else "duplicate"
            read.set_tag("GG", tag, value_type="Z")
            outfile.write(read)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", help="Input BAM file from umi_tools dedup")
    parser.add_argument("tsv", help="TSV file from umi_tools group")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    args = parser.parse_args()

    prepare_umitools_bam(args.input, args.tsv, args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
