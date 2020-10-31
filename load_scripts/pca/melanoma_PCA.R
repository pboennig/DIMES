library(Seurat)
source('dimes.R')
library(purrr)
library(glue)
library(here)

exprs <- read.csv("load_scripts/cell_line/GSE108383/GSE108383_Melanoma_fluidigm.txt", sep='\t')
make_objs <- function(cell_line, exprs) { 
  cells <- grep(cell_line, colnames(exprs))
  melanoma <- CreateSeuratObject(counts=exprs[,cells], min.cells=3, min.features=200)
  melanoma <- NormalizeData(object = melanoma, normalization.method = "LogNormalize", scale.factor = 1e4)
  melanoma <- FindVariableFeatures(melanoma)
  ScaleData(object=melanoma, features=rownames(melanoma))
}
cell_lines <- c("A375", "451Lu")
objs <- sapply(X=cell_lines, FUN=make_objs, exprs=exprs)
aucs <- lapply(X=objs, FUN=ppi_comp, run.PCA=T)
map2(aucs, names(aucs), 
     function(auc, cell_line) { write.csv(auc, here("results","aucs", "pca", glue("{cell_line}_PCA.csv")))})
