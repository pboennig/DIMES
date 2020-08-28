#!/bin/bash
#
#SBATCH --job-name=single_gastric
#SBATCH --time=12:00:00
#SBATCH --out=output/single/single_%a.out
#SBATCH --array=1-5 # job array index

module load r
echo ${SLURM_ARRAY_TASK_ID}
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/gastric_temp.txt
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/gastric_temp.txt | xargs -n1 Rscript load_scripts/load_gastric.R
