#!/bin/bash

#SBATCH --job-name=run_JAFFAL_longbench_ONT_cDNA_H146_pre_nailpolish
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=script_output/%x_%J.out
#SBATCH --error=script_output/%x_%J.err

set -euo pipefail

module purge
module load micromamba

eval "$(micromamba shell hook --shell=bash)"
micromamba activate /vast/projects/LINDTIE/condaenvs/jaffal_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

input_fastq="/vast/projects/davidson_longread/tan.j/nailpolish_analysis/raw_data/lb_bulk_h146_original.fastq"
output_dir="/vast/projects/davidson_longread/tan.j/nailpolish_analysis/output/JAFFAL/lb_bulk_h146_original"

mkdir -p "$output_dir"

if [[ ! -f "$input_fastq" ]]; then
    echo "ERROR: input FASTQ not found: $input_fastq"
    exit 1
fi

cd "$output_dir"

/vast/projects/LINDTIE/tools/JAFFA-version-2.5/tools/bin/bpipe run \
/vast/projects/LINDTIE/tools/JAFFA-version-2.5/JAFFAL.groovy \
"$input_fastq"
