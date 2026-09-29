#!/usr/bin/env bash
# Usage: prepare_for_sicelore.sh <input.bam> <output.bam>

module load samtools

# must convert extended CIGAR back to M for sicelore
# and add CS tag (=SEQ)
samtools view -h --sanitize cigarx "$1" \
  | awk 'BEGIN{OFS="\t"} /^@/{print; next} {
      if($10 != "*") $(NF+1) = "CS:Z:" $10
      print
  }' \
  | samtools view -b -o "$2"