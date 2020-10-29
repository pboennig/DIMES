library(Seurat)
source("dimes.R")
library(glue)
library(here)
library(purrr)
library(stringr)

name <- "GSM4238679_HGC-27" 
gsm_number <- str_match(name, "GSM(.*?)_")[,2] %>% as.integer
cell_line_name <- gsub("^[^_]*_", "", name)

GSE <- "GSE142750"

read_MM_triplet <- function(basename) {
  exprs <- readMM(here("load_scripts", "cell_line", GSE, glue("{basename}_matrix.mtx")))
  barcodes <- read.table(here("load_scripts", "cell_line", GSE, glue("{basename}_barcodes.tsv")))$V1
  genes <- read.table(here("load_scripts", "cell_line", GSE, glue("{basename}_genes.tsv")))$V1
  rownames(exprs) <- genes
  colnames(exprs) <- barcodes
  gastric <- CreateSeuratObject(counts=exprs, min.cells=3, min.features=200)
  gastric <- NormalizeData(object = gastric, normalization.method = "LogNormalize", scale.factor = 1e4)
  gastric <- FindVariableFeatures(gastric)
  ScaleData(object=gastric, features=rownames(gastric))
}

replicateObjects <- sapply(c(paste("GSM", gsm_number, "_", cell_line_name,"_1", sep=""), paste("GSM", gsm_number+1, "_", cell_line_name,"_2", sep="")), read_MM_triplet)
anchorSet <- FindIntegrationAnchors(object.list = replicateObjects, scale = F)
allGenes <- intersect(rownames(replicateObjects[1]), rownames(replicateObjects[2]))
mergedObject <- IntegrateData(anchorSet, normalization.method="SCT", features.to.integrate=allGenes)

metrics <- c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 'MI')
aucs <- ppi_comp(mergedObject, metrics=metrics)
write.csv(aucs, here("results", "aucs",glue("{cell_line_name}_1.csv")))
