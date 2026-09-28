#!/bin/bash

awk 'NR%4==1 || NR%4==3 {
    match($1, /^[@+]([^_]+)_#\?_([^#]+)#/, arr)
    print $1 "\tCB:Z:" arr[1] "\tUB:Z:" arr[2]
} NR%4==2 || NR%4==0 { print }' /vast/projects/davidson_longread/cheng.o/data/Max_2nd/hs/barcode08_merged/barcode08_merged_fp_with_barcodes.fastq > barcode01_with_cb_ub.fastq
