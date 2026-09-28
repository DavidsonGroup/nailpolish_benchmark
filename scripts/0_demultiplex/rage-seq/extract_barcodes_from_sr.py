#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "pysam",
# ]
# ///

import pysam
import argparse
from pathlib import Path
import sys

def extract_barcodes(bam_path: str, mm_max: int = 0) -> None:
    index = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            index += 1
            if read.is_unmapped:
                continue
            if not read.has_tag("CB") or not read.has_tag("MM"):
                continue
            if read.get_tag("MM") > mm_max:
                continue

            cb = read.get_tag("CB")
            # strip GEM group suffix (e.g. -1)
            cb = cb.rsplit("-", 1)[0]

            print(cb)

            if index % 1000000 == 0:
                print("Count", index, file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract cell barcodes from a 10x BAM file"
    )
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument(
        "--mm", type=int, default=0,
        help="Maximum MM (barcode mismatch) value to allow (default: 0)"
    )
    args = parser.parse_args()

    try:
        extract_barcodes(args.bam, args.mm)
    except BrokenPipeError:
        pass
