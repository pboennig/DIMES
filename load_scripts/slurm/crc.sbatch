#!/bin/bash
#
#SBATCH --job-name=single_gastric
#SBATCH --time=12:00:00
#SBATCH --out=output/single/crc_%a.out
#SBATCH --array=1-6 # job array index

module load r
echo ${SLURM_ARRAY_TASK_ID}
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/crc_cell_lines.txt
sed -n "${SLURM_ARRAY_TASK_ID}p" load_scripts/crc_cell_lines.txt | xargs -n1 Rscript load_scripts/load_crc.R
