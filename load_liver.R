#loading GEO database GSE103691
library(Seurat)
library(GEOquery)
source('dimes.R')

#10X data
exprs <- readMM("cell_line/GSE103867/GSE103867_matrix.mtx")
barcodes <- read.table("cell_line/GSE103867/GSE103867_barcodes.tsv")$V1
genes <- read.table("cell_line/GSE103867/GSE103867_genes.tsv")$V1
rownames(exprs) <- genes
colnames(exprs) <- barcodes

# only first two datasets (cell lines)
huh1 <- barcodes[grep("1", barcodes)]
huh7 <- barcodes[grep("2", barcodes )]
exprs <- exprs[,union(huh1,huh7)]

make_obj <- function(cells, exprs) {
  liver <- CreateSeuratObject(counts=exprs[,cells], min.cells=3, min.features=200)
  liver <- NormalizeData(object = liver, normalization.method = "LogNormalize", scale.factor = 1e4)
  liver <- FindVariableFeatures(liver)
  ScaleData(object=liver, features=rownames(liver))
}

objects <- sapply(X=list(huh1, huh7), FUN=make_obj, exprs=exprs)
huh1_aucs <- ppi_comp(objects[[1]], huh1)
huh7_aucs <- ppi_comp(objects[[2]], huh7)
write.csv(huh1_aucs, file="aucs/HUH1.csv")
write.csv(huh7_aucs, file="aucs/HUH7.csv")
rm(list=ls())
