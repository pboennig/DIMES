#!/bin/bash
#
#SBATCH --job-name=mi_gastric_duplicate
#SBATCH --time=12:00:00
#SBATCH --out=output/duplicate/dup_%a.out
#SBATCH --array=1-2 # job array index

module load r
echo ${SLURM_ARRAY_TASK_ID}
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/duplicate_gastric_cell_lines.txt
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/duplicate_gastric_cell_lines.txt | xargs -n1 srun Rscript load_scripts/load_duplicate_gastric.R
