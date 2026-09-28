#!/usr/bin/env python3

import argparse
import csv
import json
import sys


def headers(path):
    with open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                yield line.rstrip("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastq")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()

    out = open(args.output, "w", newline="")
    w = csv.writer(out)
    w.writerow(["read_name", "alignment_idx", "new_nodes", "sequence_len", "valid_nodes"])

    for header in headers(args.fastq):
        fields = header[1:].split("\t")
        name = fields[0]
        tag = next((f for f in fields[1:] if f.startswith("nA:Z:")), None)
        if tag is None:
            continue
        for idx, obj in enumerate(json.loads(tag[5:])):
            w.writerow([name, idx, obj["new_nodes"], obj["sequence_len"], obj["valid_nodes"]])

    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()