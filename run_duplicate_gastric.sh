#!/bin/bash

cat load_scripts/duplicate_gastric_cell_lines.txt | xargs -n1 Rscript load_scripts/load_duplicate_gastric.R
