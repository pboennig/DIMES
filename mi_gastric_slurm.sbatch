#!/bin/bash
#
#SBATCH --job-name=mi_gastric_array
#SBATCH --mail-user=pboennig@stanford.edu
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --array=1-4 # job array index

module load r
echo ${SLURM_ARRAY_TASK_ID}
sed -n "${SLURM_ARRAY_TASK_ID}p" gastric_cell_lines.txt
sed -n "${SLURM_ARRAY_TASK_ID}p" gastric_cell_lines.txt | xargs -n1 Rscript load_scripts/load_gastric_MI_only.R
