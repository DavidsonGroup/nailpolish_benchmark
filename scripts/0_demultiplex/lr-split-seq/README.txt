Barcodes are from https://github.com/fairliereese/LR-splitpipe/tree/master/barcodes:
- curl https://raw.githubusercontent.com/fairliereese/LR-splitpipe/refs/heads/master/barcodes/bc_data_v2.csv | cut -f 2 -d',' | tail -n +2 > BCs_R1.txt
- curl https://raw.githubusercontent.com/fairliereese/LR-splitpipe/refs/heads/master/barcodes/bc_data_v1.csv | cut -f 2 -d',' | tail -n +2 > BCs_R2.txt
- curl https://raw.githubusercontent.com/fairliereese/LR-splitpipe/refs/heads/master/barcodes/bc_data_v1.csv | cut -f 2 -d',' | tail -n +2 > BCs_R3.txt

Steps:
- Download SRR13948564.fastq e.g. using SRA-toolkit
- Run ./0_demultiplex.sh
- Run ./1_extract_cb_ub.sh
