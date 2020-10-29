#!/bin/bash
#
#SBATCH --job-name=HGC27
#SBATCH --partition=bigmem 
#SBATCH --qos=bigmem 
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-2 # job array index
#SBATCH --out=output/HGC-27_%a.out

module load r
echo "load_gastric_HGC-27_${SLURM_ARRAY_TASK_ID}.R"
Rscript load_scripts/load_gastric_HGC-27_${SLURM_ARRAY_TASK_ID}.R 
