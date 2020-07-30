library(Seurat)
source("dimes.R")
library(glue)
library(here)
library(purrr)
basename <- commandArgs(trailingOnly = TRUE)[1]
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

metrics <- c('pearson', 'spearman', 'kendall', 'bicor', 'binomial','cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
                                                          'weighted_rank', 'hamming', 'rho_p', 'phi_s')
obj <- read_MM_triplet(basename)
aucs <- ppi_comp(obj, metrics=metrics)
cell_line_name <- gsub("^[^_]*_", "", basename)
write.csv(aucs, here("aucs",glue("{cell_line_name}.csv")))
