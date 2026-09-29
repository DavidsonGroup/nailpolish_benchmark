Run ./extract_nA_tags.py on the consensus .fastq output from 1_workflow/:

For each sample, run:
  ./extract_nA_tags.py results/nailpolish/<sample>/consensus.fastq --output <sample>.csv
and place them in data/fdd_stats