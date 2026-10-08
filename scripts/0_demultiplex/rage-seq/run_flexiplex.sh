#!/bin/bash
# Usage: sbatch slurm-job-script
# Prepared By: Oliver Cheng
#              cheng.o@wehi.edu.au

#SBATCH --job-name=rageseq-demultiplex

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=32000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=8:00:00

#SBATCH --mail-user=ntfy-wehi-cheng-o@ntfy.sh
#SBATCH --mail-type=ALL

#SBATCH --output=logs/%j.out

zcat ERR3273784.fastq.gz | flexiplex -d 10x3v2 -k barcodes_corrected.txt -p 16 -e 4 > rage-seq.fastq

