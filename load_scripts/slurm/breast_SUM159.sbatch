#!/bin/bash
#
#SBATCH --job-name=SUM159
#SBATCH --time=12:00:00
#SBATCH --out=output/SUM159.out

module load r
Rscript load_scripts/load_breast_SUM159.R
