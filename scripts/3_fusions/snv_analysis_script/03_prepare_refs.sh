#!/bin/bash
#SBATCH --job-name=np_snv_03_prepare_refs
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=script_output/%x_%J.out
#SBATCH --error=script_output/%x_%J.err

# Stage 3: build exons.bed, render the Sanger WES truth document, normalise its VCF

# Two outputs are required downstream:
#   ref/H146_sanger_wes_truth_table.tsv  - read directly by the figure document
#   ref/truth_sanger.bed                 - read by 04_coverage.sh, which makes the depth files

set -euo pipefail
module purge
module load micromamba bedtools/2.31.1 R/4.5.3 pandoc/3.2
BASE=/vast/projects/davidson_longread/tan.j/nailpolish_analysis
SNP_DIR=${SNP_DIR:-/vast/projects/davidson_longread/tan.j/nailpolish_analysis/script/snv_analysis_script}

REF=$BASE/ref/GRCh38.primary_assembly.genome.add.spikein.fa
OUT=$BASE/output/SNV

GTF=/vast/scratch/users/tan.j/LINDTIE/ref/gencode.v49.primary_assembly.annotation.gtf
ENV=/vast/scratch/users/tan.j/condaenvs/clair3rna_analysis_env
export PYTHONNOUSERSITE=1

# Sanger WES CaVEMan, 2414 SNVs
TRUTH_RMD=03b_h146_sanger_wes_truth_vcf_processing.Rmd
TRUTH_STEM=H146_SangerWES
TRUTH_BED=$OUT/ref/truth_sanger.bed
EXPECTED_TRUTH_N=2414

eval "$(micromamba shell hook --shell=bash)"
micromamba activate "$ENV"

D=$SNP_DIR
echo "truth set: $TRUTH_STEM, expecting $EXPECTED_TRUTH_N records"

# exons.bed is the ONLY region BED still built here: 03b reads it to compute the `exonic`
# column, which is the figure's unthresholded-panel denominator. Merged as well as sorted,
# so overlapping transcript exons are not counted twice. LongBench's contig regex silently
# drops chr21/chr22 (it matches 20 but not 21/22); fixed to 2[0-2] here.
CTG_RE='^chr([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]'
awk '$3=="exon"' "$GTF" | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' \
  | grep -E "$CTG_RE" | bedtools sort -i - | bedtools merge -i - > "$OUT/ref/exons.bed"

# ---------- truth VCF ----------
export SNP_OUT="$OUT"
Rscript -e "rmarkdown::render('$D/$TRUTH_RMD',
                              output_dir = '$OUT/comparison/report',
                              knit_root_dir = '$D')"

RAW_VCF=$OUT/ref/${TRUTH_STEM}.raw.vcf
[[ -s $RAW_VCF ]] || { echo "ERROR: $RAW_VCF not produced by $TRUTH_RMD" >&2; exit 1; }

bcftools sort -Oz -o "$OUT/ref/${TRUTH_STEM}.sorted.vcf.gz" "$RAW_VCF"
bcftools index -f -t "$OUT/ref/${TRUTH_STEM}.sorted.vcf.gz"

# GATE: -c w WARNS about REF/reference mismatches without rewriting them
bcftools norm -f "$REF" -m -any -c w \
    -Oz -o "$OUT/ref/${TRUTH_STEM}.vcf.gz" "$OUT/ref/${TRUTH_STEM}.sorted.vcf.gz" \
    2> "$OUT/ref/${TRUTH_STEM}.norm.log"
bcftools index -f -t "$OUT/ref/${TRUTH_STEM}.vcf.gz"

n_mismatch=$(grep -ci "REF_MISMATCH\|does not match the reference" "$OUT/ref/${TRUTH_STEM}.norm.log" || true)
n_rec=$(bcftools index -n "$OUT/ref/${TRUTH_STEM}.vcf.gz")
echo "truth records: $n_rec ; REF mismatches: $n_mismatch"
cat "$OUT/ref/${TRUTH_STEM}.norm.log"
if [[ ${n_mismatch:-0} -gt 5 ]]; then
    echo "FATAL: $n_mismatch REF mismatches - truth table is not on this build. Stopping." >&2
    exit 1
fi


if [[ "$n_rec" -ne "$EXPECTED_TRUTH_N" ]]; then
    echo "FATAL: $n_rec records after norm, expected $EXPECTED_TRUTH_N" >&2
    exit 1
fi

bcftools query -f '%CHROM\t%POS0\t%END\n' "$OUT/ref/${TRUTH_STEM}.vcf.gz" \
  | sort -k1,1 -k2,2n > "$TRUTH_BED"

echo
echo "=== reference artefacts ==="
wc -l "$OUT/ref/exons.bed" "$TRUTH_BED"
