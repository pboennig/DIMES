library(Seurat)
source('dimes.R')

exprs <- read.csv("cell_line/GSE81861/GSE81861_Cell_Line_COUNT.csv", row.names=1)
# function to turn gene ids into ensembl gene IDs (ENSG)
get_gene_name <- function(id) {
  sub("[^_]*_[^_]*_", "", id, perl=TRUE) %>% sub("\\.(.*)", "", ., perl=TRUE)
}
short_genes <- sapply(rownames(exprs), get_gene_name)
repeats <- names(short_genes[duplicated(short_genes)])
exprs <- exprs[!(rownames(exprs) %in% repeats), ] 
rownames(exprs) <- short_genes[!(names(short_genes) %in% repeats)]
make_objs <- function(cell_line, exprs) { 
  cells <- grep(paste(cell_line, "_", sep=""), colnames(exprs))
  crc <- CreateSeuratObject(counts=exprs[,cells], min.cells=3, min.features=200)
  crc <- NormalizeData(object = crc, normalization.method = "LogNormalize", scale.factor = 1e4)
  crc <- FindVariableFeatures(crc)
  ScaleData(object=crc, features=rownames(crc))
}
cell_lines <- c("A549", "H1437", "HCT116", "IMR90", "GM12878", "H1")
objs <- sapply(X=cell_lines, FUN=make_objs, exprs=exprs)
aucs <- sapply(X=objs, FUN=ppi_comp) 
for (i in 1:length(cell_lines)) {
  j = 2 * i - 1
  df <- data.frame(metric=aucs[[j]],auc=aucs[[j+1]])
  file_name <- paste("aucs/", cell_lines[i], ".csv", sep="")
  write.csv(df, file=file_name)
}