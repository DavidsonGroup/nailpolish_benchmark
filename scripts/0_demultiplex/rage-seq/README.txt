1. Download tumour sapmles from https://www.ebi.ac.uk/ena/browser/view/PRJEB28878
2. Execute `./extract_barcodes_from_sr.py tumour_possorted_genome_bam.bam > all_barcodes_mm_0.txt`
3. Create barcodes count file `sort all_barcodes_mm_0.txt | uniq -c | awk '{print $2"\t"$1}' | sort -t$'\t' -k2 -rn > all_barcodes_mm_0_sorted.txt`
4. `uvx flexiplex-filter all_barcodes_mm_0_sorted.txt --max-rank 0 --outfile barcodes_corrected.txt`
   (or pipx or `pip install flexiplex-filter && flexiplex-filter`)
5. Run ./run_flexiplex.sh
