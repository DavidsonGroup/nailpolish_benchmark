#!/usr/bin/env python3
import pysam
import sys
from collections import defaultdict

INTERSECT_BED = sys.argv[1]
INPUT_BAM     = sys.argv[2]
OUTPUT_BAM    = sys.argv[3]

GTF_ATTR_COL  = 20
OVERLAP_COL   = 21

def parse_gene_id(attr_string: str) -> str:
    for field in attr_string.split(";"):
        field = field.strip()
        if field.startswith("gene_id"):
            return field.split('"')[1]
    raise ValueError(f"gene_id not found in attribute string: {attr_string}")

# read_name -> (best_overlap, gene_id)
best: dict[str, tuple[int, str]] = {}

with open(INTERSECT_BED) as fh:
    for line in fh:
        cols = line.rstrip("\n").split("\t")
        read_name   = cols[3]
        gene        = parse_gene_id(cols[GTF_ATTR_COL])
        overlap     = int(cols[OVERLAP_COL])

        if gene == ".":
            continue

        prev = best.get(read_name)
        if prev is None or overlap > prev[0]:
            best[read_name] = (overlap, gene)

with pysam.AlignmentFile(INPUT_BAM, "rb") as bam_in, \
     pysam.AlignmentFile(OUTPUT_BAM, "wb", header=bam_in.header) as bam_out:
    for read in bam_in:
        hit = best.get(read.query_name)
        if hit:
            read.set_tag("XT", hit[1], value_type="Z")
            read.set_tag("XS", "Assigned", value_type="Z")
        else:
            read.set_tag("XT", "NA", value_type="Z")
            read.set_tag("XS", "Unassigned_NoFeatures", value_type="Z")
        bam_out.write(read)