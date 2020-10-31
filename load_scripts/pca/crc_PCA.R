#!/usr/bin/env Rscript

library(Seurat)
source('dimes.R')
library(purrr)
library(glue)
library(pryr)
library(here)

args = commandArgs(trailingOnly = TRUE)
cell_line <- args[1]
print(cell_line)
exprs <- read.csv("load_scripts/cell_line/GSE81861/GSE81861_Cell_Line_COUNT.csv", row.names=1)
# function to turn gene ids into ensembl gene IDs (ENSG)
get_gene_name <- function(id) {
  sub("[^_]*_[^_]*_", "", id, perl=TRUE) %>% sub("\\.(.*)", "", ., perl=TRUE)
}
short_genes <- sapply(rownames(exprs), get_gene_name)
repeats <- names(short_genes[duplicated(short_genes)])
exprs <- exprs[!(rownames(exprs) %in% repeats), ] 
rownames(exprs) <- short_genes[!(names(short_genes) %in% repeats)]

cells <- grep(paste(cell_line, "_", sep=""), colnames(exprs))
crc <- CreateSeuratObject(counts=exprs[,cells], min.cells=3, min.features=200)
crc <- NormalizeData(object = crc, normalization.method = "LogNormalize", scale.factor = 1e4)
crc <- FindVariableFeatures(crc)
crc <- ScaleData(object=crc, features=rownames(crc))

aucs <- ppi_comp(crc, run.PCA = T)
write.csv(aucs, here("results", "aucs", "pca", glue("{cell_line}_PCA.csv")))
