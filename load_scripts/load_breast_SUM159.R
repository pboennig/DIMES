source("dimes.R")
library(Seurat)
library(here)

GSE131984.counts <- read.csv("load_scripts/cell_line/GSE131984/GSE131984_scRNA_counts.txt",
                             sep="\t")
colnames(GSE131984.counts) <- gsub("\\.", "-", colnames(GSE131984.counts))
  
GSE131984.meta <- read.csv("load_scripts/cell_line/GSE131984/GSE131984_scRNA_metadata.txt",
                           sep="\t")
GSE131984.meta.pre <- GSE131984.meta[which(GSE131984.meta$orig.ident == "Pre"),]
GSE131984.pre.counts <- GSE131984.counts[,which(colnames(GSE131984.counts) %in% 
                                                  rownames(GSE131984.meta.pre))]
GSE131984.obj <- CreateSeuratObject(GSE131984.counts, 
                                    min.cells = 3, min.features = 200)
GSE131984.obj <- AddMetaData(object = GSE131984.obj, metadata = GSE131984.meta) 

GSE131984.pre.obj <- subset(GSE131984.obj, orig.ident == "Pre")
GSE131984.pre.obj[["percent.mt"]] <- PercentageFeatureSet(GSE131984.pre.obj, 
                                                          pattern = "^MT-")
GSE131984.pre.obj <- NormalizeData(GSE131984.pre.obj)
GSE131984.pre.obj <- ScaleData(GSE131984.pre.obj)
GSE131984.pre.obj <- FindVariableFeatures(GSE131984.pre.obj)
GSE131984.aucs <- ppi_comp(GSE131984.pre.obj, colnames(GSE131984.pre.obj))
write.csv(GSE131984.aucs,here("results", "aucs", "SUM159.csv"))
