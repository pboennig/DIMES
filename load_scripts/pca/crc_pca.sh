#!/bin/bash
cat load_scripts/crc_cell_lines.txt | xargs -n1 Rscript load_scripts/load_crc_PCA.R
