#loading GEO database GSE103691
library(Seurat)
library(GEOquery)
exprs <- readMM("GSE103867/GSE103867_matrix.mtx")
barcodes <- read.table("GSE103867/GSE103867_barcodes.tsv")$V1
genes <- read.table("GSE103867/GSE103867_genes.tsv")$V2
rownames(exprs) <- genes
colnames(exprs) <- barcodes
liver <- CreateSeuratObject(counts=exprs, min.cells=3, min.features=200)

