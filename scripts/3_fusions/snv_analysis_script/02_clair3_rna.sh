#!/bin/bash
#SBATCH --job-name=np_snv_02_clair3_rna
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=0-1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=script_output/%x_%A_%a.out
#SBATCH --error=script_output/%x_%A_%a.err

# Stage 2: genome-wide Clair3-RNA, one array task per arm

set -euo pipefail
module purge
module load apptainer/1.4.1

BASE=/vast/projects/davidson_longread/tan.j/nailpolish_analysis
REF=$BASE/ref/GRCh38.primary_assembly.genome.add.spikein.fa
OUT=$BASE/output/SNV

SIF=/vast/scratch/users/tan.j/tools/clair3-rna_v0.2.2.sif
PLATFORM=ont_r10_dorado_cdna

ARMS=(original nailpolish)

export PYTHONNOUSERSITE=1
export APPTAINER_CACHEDIR=/vast/scratch/users/tan.j/.apptainer_cache
export APPTAINER_TMPDIR=/vast/scratch/users/tan.j/.apptainer_tmp

ARM=${ARMS[$SLURM_ARRAY_TASK_ID]}
BAM=$OUT/GenomeAlignment/${ARM}.sorted.bam
ODIR=$OUT/clair3_rna/${ARM}
mkdir -p "$ODIR"

[[ -f $BAM ]] || { echo "ERROR: BAM not found: $BAM (run 01_minimap2 first)" >&2; exit 1; }

apptainer exec -B "$OUT","$(dirname "$REF")" "$SIF" \
  env CONDA_PREFIX=/opt/conda/envs/clair3_rna PYTHONNOUSERSITE=1 \
  /opt/bin/run_clair3_rna \
    --bam_fn      "$BAM" \
    --ref_fn      "$REF" \
    --output_dir  "$ODIR" \
    --platform    "$PLATFORM" \
    --threads     "$SLURM_CPUS_PER_TASK" \
    --include_all_ctgs \
    --enable_phasing_model

echo
echo "=== [$ARM] outputs ==="
ls -la "$ODIR"/*.vcf.gz

module load bcftools/1.23
for v in "$ODIR"/output.vcf.gz "$ODIR"/output_enable_phasing.vcf.gz; do
    if [[ -f $v ]]; then
        bcftools index -f -t "$v"
        echo "$v: $(bcftools index -n "$v") records"
    fi
done
