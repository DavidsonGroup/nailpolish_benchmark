#!/bin/bash
#SBATCH --job-name=np_snv_01_minimap2
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=script_output/%x_%A_%a.out
#SBATCH --error=script_output/%x_%A_%a.err

# Stage 1: splice-aware genome alignment, one array task per arm

set -euo pipefail
module purge
module load minimap2/2.28 samtools/1.23.1

BASE=/vast/projects/davidson_longread/tan.j/nailpolish_analysis
RAW=$BASE/raw_data
REF=$BASE/ref/GRCh38.primary_assembly.genome.add.spikein.fa
OUT=$BASE/output/SNV

ARMS=(original nailpolish)
declare -A FQ=(
  [original]=$RAW/lb_bulk_h146_original.fastq
  [nailpolish]=$RAW/lb_bulk_h146_nailpolish.fastq
)

ARM=${ARMS[$SLURM_ARRAY_TASK_ID]}
FQ_IN=${FQ[$ARM]}
BAM=$OUT/GenomeAlignment/${ARM}.sorted.bam

SORT_TMP_DIR=/vast/scratch/users/tan.j/np_snp_sort_${ARM}_${SLURM_JOB_ID}
mkdir -p "$SORT_TMP_DIR"
trap 'rm -rf "$SORT_TMP_DIR"' EXIT

echo "[$ARM] fastq   : $FQ_IN"
echo "[$ARM] ref     : $REF"
echo "[$ARM] out bam : $BAM"
echo "[$ARM] sort tmp: $SORT_TMP_DIR ($(df -h --output=avail "$(dirname "$SORT_TMP_DIR")" | tail -1 | tr -d ' ') free)"
[[ -f $FQ_IN ]] || { echo "ERROR: input FASTQ not found: $FQ_IN" >&2; exit 1; }

minimap2 -ax splice:hq -y --MD -t "$SLURM_CPUS_PER_TASK" "$REF" "$FQ_IN" \
  | samtools sort -@ 8 -m 4G -T "$SORT_TMP_DIR/sort" -o "$BAM" -
samtools index -@ 8 "$BAM"
samtools flagstat -@ 8 "$BAM" > "$OUT/GenomeAlignment/${ARM}.flagstat"

echo
echo "=== [$ARM] flagstat ==="
cat "$OUT/GenomeAlignment/${ARM}.flagstat"

# Gate: -y must have carried the comment tags through
n_tag=$(samtools view "$BAM" 2>/dev/null | head -100000 | grep -c 'nT:Z:' || true)
echo
echo "[$ARM] alignments carrying nT:Z: in first 100k records: $n_tag"
if [[ "$ARM" == "nailpolish" && "$n_tag" -eq 0 ]]; then
    echo "ERROR: nailpolish BAM has no nT:Z: tags - minimap2 -y did not take effect." >&2
    exit 1
fi

samtools quickcheck -v "$BAM" && echo "[$ARM] quickcheck OK"
