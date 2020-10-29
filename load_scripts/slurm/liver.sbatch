#!/bin/bash
#
#SBATCH --job-name=liver
#SBATCH --partition=bigmem 
#SBATCH --qos=bigmem 
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --out=output/liver.out

module load r
Rscript load_scripts/load_liver.R
