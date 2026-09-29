#!/usr/bin/env python3
"""
Prepare a nailpolish-aligned BAM for summarise_bam.py by adding GG:Z group tags.

Reads a BAM/SAM with nailpolish tags (nT, nI, nC) and writes a new BAM where
each relevant read carries a GG:Z tag identifying its group:
  - np_simplex        (nT == "simplex")
  - np_consensus      (nT == "consensus")
  - np_best_original            (original with lowest error rate in its cluster)
  - np_worst_original           (original with highest error rate in its cluster)
  - np_random_original          (randomly selected original from its cluster, seed=42)
  - np_multiread_consensus        } only for clusters with more than one original
  - np_multiread_best_original   } read (checked via nL tag for consensus,
  - np_multiread_worst_original  } len(originals) > 1 for original reads)
  - np_multiread_random_original }

A single read may satisfy multiple criteria; it receives all applicable labels
as a comma-separated GG value. Reads not matching any group are omitted.

Clusters are assumed to be contiguous in the input BAM.
"""

import argparse
import random

import pysam

CIGAR_COLUMNS = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B", "edit_dist"]
rng = random.Random(42)


def get_tag_or_none(read, tag):
    if read.has_tag(tag):
        return read.get_tag(tag)
    return None


def edit_dist_rate(read):
    cigar = dict(zip(CIGAR_COLUMNS, read.get_cigar_stats()[0]))
    denom = cigar["M"] + cigar["I"] + cigar["="] + cigar["X"]
    try:
        return cigar["edit_dist"] / denom
    except ZeroDivisionError:
        return float("nan")


def flush_originals(originals, outfile):
    if not originals:
        return

    best = min(originals, key=lambda x: x[1])
    worst = max(originals, key=lambda x: x[1])
    random_pick = rng.choice(originals)

    # groups to be selected
    groups: dict[int, tuple] = {}  # id(read) -> (read, [group_labels])

    for label, item in [
        ("np_best_original", best),
        ("np_worst_original", worst),
        ("np_random_original", random_pick),
    ]:
        rid = id(item[0])
        if rid not in groups:
            groups[rid] = (item[0], [])
        groups[rid][1].append(label)

    if len(originals) > 1:
        for label, item in [
            ("np_multiread_best_original", best),
            ("np_multiread_worst_original", worst),
            ("np_multiread_random_original", random_pick),
        ]:
            groups[id(item[0])][1].append(label)

    for read, labels in groups.values():
        read.set_tag("GG", ",".join(labels), value_type="Z")
        outfile.write(read)


def prepare_nailpolish_bam(input_path, output_path):
    mode_r = "r" if input_path.endswith(".sam") else "rb"

    with pysam.AlignmentFile(input_path, mode_r) as infile, \
         pysam.AlignmentFile(output_path, "w", header=infile.header) as outfile:

        current_cluster_id = None
        current_group_originals = []

        for read in infile:
            if read.is_secondary or read.is_supplementary:
                continue

            read_type = get_tag_or_none(read, "nT")

            if read_type == "simplex":
                read.set_tag("GG", "np_simplex", value_type="Z")
                outfile.write(read)
            elif read_type == "consensus":
                cluster_length = get_tag_or_none(read, "nL")
                gg = "np_consensus,np_multiread_consensus" if cluster_length and cluster_length > 1 else "np_consensus"
                read.set_tag("GG", gg, value_type="Z")
                outfile.write(read)
            elif read_type == "original":
                cluster_id = (get_tag_or_none(read, "nI"), get_tag_or_none(read, "nC"))

                if cluster_id != current_cluster_id:
                    flush_originals(current_group_originals, outfile)

                    # start new cluster
                    current_group_originals = []
                    current_cluster_id = cluster_id

                current_group_originals.append((read, edit_dist_rate(read)))

        flush_originals(current_group_originals, outfile)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", help="Input BAM/SAM file")
    parser.add_argument("-o", "--output", required=True, help="Output SAM file")
    args = parser.parse_args()

    prepare_nailpolish_bam(args.input, args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
