#!/bin/bash

#SBATCH --job-name=2_remove_remnants
#SBATCH --partition=regular
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --mail-user=ntfy-wehi-cheng-o@ntfy.sh
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

set -euo pipefail

module load cutadapt/

input_fastq="temp/longbench_bulk_H146_trimmed.fastq"
intermediate_fastq="temp/longbench_bulk_H146_remove_adapters.fastq"
output_fastq="longbench_bulk_H146_cleaned.fastq"

# -O: Requires ≥18 bp overlap; -e: Allows 10% errors; -m: Keeps reads ≥100 nt (-m) after trimming, -j: Uses 16 threads   
# Step 1: fixed adapters
cutadapt -j 16 \
  --revcomp -O 18 -e 0.1 -m 100 \
  -g ^CTTGCGGGCGGCGGACTCTCCTCT \
  -a CTTGCGGGCGGCGGACTCTCCTCT$ \
  -o "$intermediate_fastq" \
  "$input_fastq"

# Step 2: poly-A / homopolymers
cutadapt -j 16 \
  --revcomp -m 100 \
  --poly-a \
  -a A{10} -a T{10} \
  -g ^A{10} -g ^T{10} \
  -n 2 \
  -o "$output_fastq" \
  "$intermediate_fastq"
