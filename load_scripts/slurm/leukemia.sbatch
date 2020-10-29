#!/bin/bash
#
#SBATCH --job-name=leuk
#SBATCH --time=12:00:00
#SBATCH --out=output/leukemia.out

module load r
Rscript load_scripts/load_leukemia.R
