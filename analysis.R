knitr::opts_chunk$set(
    cache = TRUE,
    cache.lazy = FALSE,
    tidy = TRUE
)
library(Seurat)
library(readxl)
library(dplyr)
library(tibble)

library(wesanderson) #colors
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
lung <- FindVariableFeatures(lung)

genes <- head(VariableFeatures(lung), 70) # top 70 most variable genes
counts <- GetAssayData(object=lung,slot="counts")[genes, ] %>% as.matrix # rows = genes, columns = cells

comb.cells <- metadata[metadata$'PatientNumber MS' %in% c(3:4) &
                       metadata$cell %in% colnames(counts),]$cell

rm(list=setdiff(ls(), c("comb.cells", "counts")))
source("mi_scrna.R")
print("Creating distance matrices...")
metrics <- c('pearson', 'spearman', 'kendall', 'bicor', 'zi_kendall', 'binomial', 'MI', 
  'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
  'weighted_rank', 'hamming')

for (metric in metrics) {
    graph <- dist_matrix(cells=comb.cells,counts=counts,metric=metric) %>%
                as_tbl_graph(directed=FALSE)
    save_gene_clusters(graph, metric)
    save_network(graph, metric)
}

