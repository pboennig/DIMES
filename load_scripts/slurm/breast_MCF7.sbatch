#!/bin/bash
#
#SBATCH --job-name=MCF7
#SBATCH --partition=bigmem 
#SBATCH --qos=bigmem 
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --out=output/MCF7.out

module load r
Rscript load_scripts/load_breast_MCF7.R
