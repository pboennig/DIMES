#!/bin/bash
#SBATCH -output=output/dimes_fix/crc.out

module load r
Rscript load_scripts/load_crc.R A549
