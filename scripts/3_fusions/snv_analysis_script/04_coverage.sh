#!/bin/bash
#SBATCH --job-name=np_snv_04_coverage
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=script_output/%x_%A_%a.out
#SBATCH --error=script_output/%x_%A_%a.err

# Stage 4: per-site depth at the Sanger truth positions, one array task per arm

set -euo pipefail
module purge
module load samtools/1.23.1

BASE=/vast/projects/davidson_longread/tan.j/nailpolish_analysis
OUT=$BASE/output/SNV

# Must match 03_prepare_refs.sh, which writes this BED
TRUTH_BED=$OUT/ref/truth_sanger.bed

ARMS=(original nailpolish)
ARM=${ARMS[$SLURM_ARRAY_TASK_ID]}
BAM=$OUT/GenomeAlignment/${ARM}.sorted.bam
DEPTH_DIR=$OUT/depth
mkdir -p "$DEPTH_DIR"

[[ -f $BAM ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[[ -s $TRUTH_BED ]] || { echo "ERROR: run 03_prepare_refs.sh first" >&2; exit 1; }

# -a reports zero-depth positions too, so uncovered truth sites appear as depth 0 rather than vanishing from the table 
DEPTH_TSV=$DEPTH_DIR/truth_depth_sanger_${ARM}.tsv
samtools depth -a -b "$TRUTH_BED" "$BAM" \
  | awk -v a="$ARM" 'BEGIN{OFS="\t"; print "CHROM","POS","depth","arm"} {print $1,$2,$3,a}' \
  > "$DEPTH_TSV"

n_rows=$(( $(wc -l < "$DEPTH_TSV") - 1 ))
n_span=$(awk '{n += $3 - $2} END{print n+0}' "$TRUTH_BED")
echo "[$ARM] truth-site depth rows: $n_rows (BED spans $n_span bp)"
awk 'NR>1{n++; if($3>=10) c++} END{printf "[covered >=10x: %d / %d]\n", c+0, n+0}' "$DEPTH_TSV"

# GATE: `samtools depth -a` silently omits contigs with no reads at all and a truth set whose positions are 
# missing from this file joins to NA downstream rather than to zero, which reads as a plausible small number 
# instead of failing. Assert the row count matches the BED span so that loss is caught here.
if [[ "$n_rows" -ne "$n_span" ]]; then
    echo "FATAL: $n_rows depth rows for a $n_span bp truth BED. samtools depth dropped" >&2
    echo "       positions (see METHODS.md 6.9); the depth-conditioned analysis would be" >&2
    echo "       silently wrong. Stopping." >&2
    exit 1
fi
