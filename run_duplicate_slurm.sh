#!/bin/bash
#
#SBATCH --job-name=mi_gastric_duplicate

cat load_scripts/duplicate_gastric_cell_lines.txt | xargs -n1 srun Rscript load_scripts/load_duplicate_gastric_MI_only.R
