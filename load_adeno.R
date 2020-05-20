library(Seurat)
library(readxl)
load("./adenoallraw.RData")
type <- "tumor"
lung <- CreateSeuratObject(adeno.all$TUMOR, min.cells = 3, min.features = 200, project = "nsclc")
metadata <- read_excel("./MetaData.xlsx", sheet = 1)
strsplit.ind2 <- seq(from=2, by=2, length.out = nrow(metadata))
metadata <- metadata %>%
    mutate(piece = unlist(strsplit(Patient_piece, split="-"))[strsplit.ind2]) %>%
    filter(CellType == type)
lung <- AddMetaData(object = lung, metadata = metadata)
lung <- NormalizeData(object = lung, normalization.method = "LogNormalize", scale.factor = 1e4)
lung <- FindVariableFeatures(lung)
all.genes <- rownames(x = lung)
lung <- ScaleData(object = lung, features = all.genes)
cells <- metadata[metadata$'PatientNumber MS' %in% c(3:4),]$cell
rm(list = c("metadata", "strsplit.ind2", "all.genes"))
