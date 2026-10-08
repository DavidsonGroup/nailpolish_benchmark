#!/usr/bin/env python3
# /// script
# dependencies = ["tqdm"]
# ///

# the purpose of this script is to generate a dataset with the correct UB, CB
# tags added.

# You can download H146_universal_clusters.csv and H146_bulk_ONT.fastq from
# the LongBench S3 bucket.

import sys
from tqdm import tqdm

BARCODES_FILE = "/vast/projects/davidson_longread/cheng.o/proj/25.06.17-np-benchmark/snakemake_bulk/datasets_bulk/longbench_barcodes/H146_universal_clusters.csv"
FASTQ = "/vast/projects/davidson_longread/cheng.o/proj/25.06.17-np-benchmark/snakemake_bulk/datasets_bulk/longbench_fastq/H146_bulk_ONT.fastq"
OUTPUT_FASTQ = "longbench_bulk_H146_raw.fastq"

table = {}
with open(BARCODES_FILE) as f:
    next(f)  # skip header
    for line in f:
        rid, ub = line.strip().split(";")
        table[rid] = ub

max_len = max(len(ub) for ub in table.values())

count_with_bc = 0
count_without_bc = 0
skip = False
with open(FASTQ) as fin, open(OUTPUT_FASTQ, "w") as fout:
    pbar = tqdm(unit="reads")
    for idx, line in enumerate(fin):
        if idx % 4 == 0:  # header line
            line = line.rstrip()
            rid = line.split("\t")[0][1:]

            if rid in table:
                ub = table[rid].ljust(max_len, "T")
                fout.write(f"@{rid}\tCB:Z:ATCGATCGATCG\tUB:Z:{ub}\n")
                count_with_bc += 1
                skip = False
            else:
                count_without_bc += 1
                skip = True
            pbar.update(1)
        elif not skip:
            fout.write(line)

print(f"Reads with barcode: {count_with_bc}")
print(f"Reads without barcode: {count_without_bc}")