This directory contains the data used by the plotting code
in the plots/ folder. You must first generate the data yourself
using the preprocessing code in scripts/.

Pre-generated data included in this repo
- benchmarks/         results/_benchmarks/* from scripts/1_workflow
- csv/                results/_results/*.csv from scripts/1_workflow
- read_counts.csv     see below

Not included - must run yourself and move/symlink in
- fusions/            output from scripts/3_fusions
- fdd_stats/          output from scripts/2_fdd_stats


Format of read_counts.csv:
This spreadsheet's data is derived from the Nailpolish stdout logs:
- sample: sample name
- pre_sin, pre_dup_groups, pre_dup_reads:                    from results/nailpolish/<sample>/summary_stdout.txt
- post_sin, post_dup_groups, post_dup_reads, post_filtered:  from results/nailpolish/<sample>/np.log