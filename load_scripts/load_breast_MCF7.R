source("dimes.R")
library(Seurat)
library(GEOquery)
library(here)

GSE114459.counts <- read.csv("load_scripts/cell_line/GSE114459/GSE114459_Polyak.csv")
GSE114459.geo <- getGEO('GSE114459')
pdata <- pData(GSE114459.geo$GSE114459_series_matrix.txt.gz)

GSE114459.parental <- GSE114459.counts[,grep("Parental",colnames(GSE114459.counts))]
GSE114459.obj <- CreateSeuratObject(GSE114459.parental, 
                                    min.cells = 3, min.features = 200)
GSE114459.obj[["percent.mt"]] <- PercentageFeatureSet(GSE114459.obj, 
                                                          pattern = "^MT-")
VlnPlot(GSE114459.obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

GSE114459.obj <- subset(GSE114459.obj, percent.mt < 10)

GSE114459.obj <- NormalizeData(GSE114459.obj)
GSE114459.obj <- ScaleData(GSE114459.obj)

GSE114459.obj <- FindVariableFeatures(GSE114459.obj) #necessary for ppi_comp
gse114459.aucs <- ppi_comp(GSE114459.obj, colnames(GSE114459.obj))
write.csv(gse114459.aucs,here("results", "aucs", "MCF7.csv"))
