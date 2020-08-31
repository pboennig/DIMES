#!/bin/bash
#
#SBATCH --job-name=mi_gastric_duplicate
#SBATCH --time=12:00:00
#SBATCH --array=1-2 # job array index
#SBATCH --out=output/dup_%a.out

module load r
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/duplicate_gastric_cell_lines.txt | xargs -n1 Rscript load_scripts/load_duplicate_gastric.R 
