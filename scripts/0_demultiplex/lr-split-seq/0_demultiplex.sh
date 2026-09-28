#!/bin/bash
# Usage: sbatch slurm-job-script
# Prepared By: Oliver Cheng
#              cheng.o@wehi.edu.au

#SBATCH --job-name=lr-split-seq-demultiplex
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=32000
#SBATCH --time=00:10:00
#SBATCH --mail-user=ntfy-wehi-cheng-o@ntfy.sh
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%j.out

INPUT="SRR13948564.fastq"

FLEXIPLEX="./flexiplex-v1.02.6"
LINKER1="GTGGCCGATGTTTCGCATCGGCGTACGACT"
LINKER2="ATCCACGTGCTTGAGACTGTGG"


# cat $INPUT | \
# 	$FLEXIPLEX -u "??????????" -b "????????" -x $LINKER1 -k BCs_R3.txt -f 2 -e 1 -p 4 | \
# 	$FLEXIPLEX -b "????????" -x $LINKER2 -k BCs_R2.txt -f 2 -e 1 -p 4 | \
# 	sed "/[@,+]/! s/^/987654321/g" | \
# 	$FLEXIPLEX -x "987654321" -b "????????" -k BCs_R1.txt -f 0 -e 1 -p 4 \
# 	> SRR13948564.DMX.fastq

cat $INPUT | \
	$FLEXIPLEX -u "??????????" -b "????????" -x $LINKER1 -k BCs_R3.txt -f 4 -e 2 -p 4 | \
	sed "/[@,+]/! s/^/987654321/g" | \
	$FLEXIPLEX -x "987654321" -b "????????" -x $LINKER2 -k BCs_R2.txt -f 2 -e 1 -p 4 | \
	sed "/[@,+]/! s/^/987654321/g" | \
	$FLEXIPLEX -x "987654321" -b "????????" -k BCs_R1.txt -f 0 -e 2 -p 4 \
	> SRR13948564.DMX.fastq