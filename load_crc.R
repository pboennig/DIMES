library(Seurat)
library(GEOquery)
source('dimes.R')

geo <- getGEO("GSE81861")
data <- read.csv("crc/GSE81861_Cell_Line_COUNT.csv")
rownames(exprs) <- genes
colnames(exprs) <- barcodes
liver <- CreateSeuratObject(counts=exprs, min.cells=3, min.features=200)
liver <- NormalizeData(object = liver, normalization.method = "LogNormalize", scale.factor = 1e4)
liver <- FindVariableFeatures(liver)
liver <- ScaleData(object=liver, features=rownames(liver))
#only want the first 2 samples (cell lines)
liver_cells <- colnames(liver)[grep("1|2", colnames(liver))]

crc_aucs <- ppi_comp(liver, liver_cells)