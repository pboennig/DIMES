#!/bin/bash
#SBATCH --job-name=melanoma
#SBATCH --partition=bigmem 
#SBATCH --qos=bigmem 
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --out=output/melanoma.out

module load r
Rscript load_scripts/load_melanoma.R
