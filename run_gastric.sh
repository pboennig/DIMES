#!/bin/bash

cat load_scripts/gastric_cell_lines.txt | xargs -n1 Rscript load_scripts/load_gastric.R
