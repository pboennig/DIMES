library(Seurat)
source('dimes.R')

exprs <- read.csv("cell_line/GSE108383/GSE108383_Melanoma_fluidigm.txt", sep='\t')
make_objs <- function(cell_line, exprs) { 
  cells <- grep(cell_line, colnames(exprs))
  melanoma <- CreateSeuratObject(counts=exprs[,cells], min.cells=3, min.features=200)
  melanoma <- NormalizeData(object = melanoma, normalization.method = "LogNormalize", scale.factor = 1e4)
  melanoma <- FindVariableFeatures(melanoma)
  ScaleData(object=melanoma, features=rownames(melanoma))
}
cell_lines <- c("A375", "451Lu")
objs <- sapply(X=cell_lines, FUN=make_objs, exprs=exprs)
aucs <- sapply(X=objs, FUN=ppi_comp)
for (i in 1:length(cell_lines)) {
  j = 2 * i - 1
  df <- data.frame(metric=aucs[[j]],auc=aucs[[j+1]])
  file_name <- paste("aucs/", cell_lines[i], ".csv", sep="")
  write.csv(df, file=file_name)
}