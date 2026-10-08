#!/bin/bash
# Usage: sbatch slurm-job-script
# Prepared By: Oliver Cheng
#              cheng.o@wehi.edu.au

#SBATCH --job-name=demultiplex-longbench
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=32000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=16:00:00
#SBATCH --mail-user=ntfy-wehi-cheng-o@ntfy.sh
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%j.out

zcat SC_ONT.fastq.gz | \
        flexiplex -d 10x3v3 -p 32 -n LB_SC_ONT -k all_barcodes.txt > lb_sc_ont.fastq

zcat SC_PB.fastq.gz | \
        flexiplex -d 10x3v3 -p 32 -n LB_SC_PB -k all_barcodes.txt > lb_sc_pb.fastq
