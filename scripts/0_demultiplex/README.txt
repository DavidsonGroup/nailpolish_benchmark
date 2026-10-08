This folder contains demultiplexing scripts to reproduce the input .fastq files
provided to the pipeline.

Accessions for each dataset can be found in the Supplementary Tables.

LR-split-seq, LongBench SC ONT, LongBench SC PB, LongBench Bulk H146,
Visium, DNM1L, RAGE-seq:
    See the corresponding folder in this directory.


# hiPSC (gridION, gridION_Q20, promethION) & scmixology2:
Each file was demultiplexed using Flexiplex as follows with barcodes derived from
the long-read data:

    flexiplex -f 0 -p $THREADS $IN_FASTQ > raw_barcodes.txt
    flexiplex-filter --whitelist 3M-february-2018.txt \
                    --outfile filtered_barcodes.txt \
                    flexiplex_barcodes_counts.txt
    flexiplex $IN_FASTQ -k filtered_barcodes.txt -p $THREADS > $OUT_FASTQ

