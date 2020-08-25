source("dimes.R")
library(GEOquery)
library(Seurat)

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
write.csv(gse114459.aucs,file="aucs/MCF7.csv")

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
write.csv(GSE131984.aucs, file="aucs/SUM159.csv")
