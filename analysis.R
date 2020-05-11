knitr::opts_chunk$set(
    cache = TRUE,
    cache.lazy = FALSE,
    tidy = TRUE
)
library(Seurat)
library(readxl)
library(dplyr)
library(tibble)

source("mi_scrna.R")
load("./adenoallraw.RData")
type <- "tumor"
lung <- CreateSeuratObject(adeno.all$TUMOR, min.cells = 3, min.features = 0, project = "nsclc")

metadata <- read_excel("./MetaData.xlsx", sheet = 1)
strsplit.ind2 <- seq(from=2, by=2, length.out = nrow(metadata))
metadata <- metadata %>%
    mutate(piece = unlist(strsplit(Patient_piece, split="-"))[strsplit.ind2]) %>%
    filter(CellType == type)

lung <- AddMetaData(object = lung, metadata = metadata)

lung <- NormalizeData(object = lung, normalization.method = "LogNormalize", scale.factor = 1e4)
all.genes <- rownames(x = lung)
lung <- ScaleData(object = lung, features = all.genes)
cells <- metadata[metadata$'PatientNumber MS' %in% c(3:4),]$cell
metrics <- c('pearson', 'spearman', 'kendall', 'bicor', 'zi_kendall', 'binomial', 'MI', 
  'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
  'weighted_rank', 'hamming')

graphs <- build_gene_graphs(lung, cells, metrics)

for (metric in metrics) {
    save_gene_clusters(graphs$metric, metric)
    save_network(graphs$metric, metric)
}
