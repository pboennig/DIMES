# taken from: https://gist.github.com/stevenworthington/3178163
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})

packages <- c("tidygraph",
"ggplot2", "ggraph", "purrr",
"Matrix", "dismay", "pryr",
"wesanderson", "tibble", "dplyr", "WGCNA",
"STRINGdb", "igraph", "pROC", "glue")
ipak(packages)

library(tidygraph)
library(ggplot2)
library(ggraph)
library(purrr)
library(Matrix)
library(dismay)
library(pryr)
library(wesanderson) #colors
library(tibble)
library(dplyr)
library(WGCNA)
library(STRINGdb)
library(igraph)
library(pROC)
library(glue)

#' Calculates the AUCs for different *metrics* on a particular Seurat *obj* by comparing
#' to STRINGdb network
ppi_comp <- function(obj, cells=colnames(obj), num_genes= nrow(obj), metrics=c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 'MI',
                                                          'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
                                                          'weighted_rank', 'hamming', 'rho_p', 'phi_s'), run.PCA = FALSE, num.PCs = 30) {
  if (run.PCA) {
	# rho and phi require non-negative values, which we can't guarantee in PCA space
	# jaccard returns 0 if all entries are nonzero, which is typically true in PCA space
	metrics <- metrics[!(metrics %in% c('rho_p', 'phi_s', 'jaccard'))]
  }
  if (num_genes < nrow(obj)) {
    genes <- data.frame(gene=head(VariableFeatures(obj), num_genes)) 
  } else {
    genes <-data.frame(gene=VariableFeatures(obj))
  }
  string_db <- STRINGdb$new(version="11", score_threshold=50, input_directory="STRINGdb/", species=9606)
  mapped <- string_db$map(genes, "gene", removeUnmappedRows=TRUE)
  graphs <- build_gene_graphs(obj, cells, mapped, metrics, run.PCA=run.PCA, num.PCs=num.PCs)
  # filter out genes that didn't appear in assayData 
  print("DIMES: Graphs built")  
  # mutinfo functions remove '-' characters
  if ('MI' %in% metrics) {
    pears_names <- graphs$pearson %>% pull(name)
    graphs$MI <- graphs$MI %>% activate(nodes) %>% mutate(name = pears_names)
  }
  
  string_subset <- string_db$get_subnetwork(mapped$STRING_id)
  string_subset <- delete.edges(string_subset, which(E(string_subset)$experiments < 400)) %>%
                   as_adjacency_matrix %>% as.matrix
  string_subset <- string_subset / 2# adjacency matrix has 2 if edge exists, want it to be 0/1
  scores <- lapply(graphs, calc_score, mapped, string_subset)
  rm(graphs) # no longer need graphs (they take up a bunch of memory)
  print("DIMES: Scores calculated")  
  
  string_subset <- string_subset[rownames(scores[[1]]), ]
  string_subset <- string_subset[, colnames(scores[[1]])]
  names(scores) <- metrics
  scores <- lapply(scores, as.vector)
  status <- as.vector(string_subset) # score matrix w/ labels in status column
  aucs <- sapply(X=scores, FUN=function(pred) auc(response=status, predictor=pred, direction="<")[[1]])
  
  res <- data.frame(metric=metrics, auc=aucs)
  rownames(res) <- NULL
  print("DIMES: Cleaning up...")
  rm(list=c("string_subset", "scores", "mapped", "string_db")) 
  gc()
  res
}

#' Calculate scores for AUC calculation based on correspondence between genes and STRING_ids
calc_score <- function(g, mapped, string_subset) {
  ig <- as.igraph(g) 
  ig <- set.vertex.attribute(ig, "name", value=mapped[which(mapped$gene==V(ig)$name),]$STRING_id)
  score <- get.adjacency(ig, attr="weight")
  mapped_string_ids <- intersect(rownames(string_subset), rownames(score))
  correspondence <- match(mapped_string_ids, rownames(score))
  score <- score[correspondence,]
  score[,correspondence]
}


#' Given seurat object *obj*, cell list *cells* and *metrics*, return *graphs*,
#' where graphs$m is a tbl_graph composed using distance metric m from metrics
build_gene_graphs <- function(obj, cells, mapped, metrics, run.PCA, num.PCs) {
    counts <- GetAssayData(object=obj,slot="counts",assay="RNA")
    scale_assay <- "RNA"
    if ("integrated" %in% Assays(obj)) {
      scale_assay <- "integrated"
    }
    log_counts <- GetAssayData(object=obj,slot="scale.data",assay=scale_assay)
    mapped <- mapped[mapped$gene %in% rownames(counts),]
    counts <- counts[mapped$gene,] %>% as.matrix
    log_counts <- log_counts[mapped$gene,] %>% as.matrix
     
    cells <- intersect(cells, colnames(counts)) # only take cells that are actually in our counts
    counts <- counts[, cells]
    log_counts <- log_counts[, cells]
    if (run.PCA) {
	counts <- prcomp(counts)$x[,1:num.PCs]
	log_counts <- prcomp(log_counts)$x[,1:num.PCs]
    }
    print(dim(counts))
    count_mat <- t(counts) # dismay expects as transpose of Seurat default
    log_mat <- t(log_counts) # dismay expects as transpose of Seurat default
    graphs <- lapply(metrics, dist_matrix, count_mat, log_mat)
    names(graphs) <- metrics
    graphs
}

#' Wrapper around dismay, switch on metrics that need raw counts
dist_matrix <-function(metric, count_mat, log_mat, eps=.001) {
    if (metric %in% c('phi_s', 'rho_p')) {
        sim_mat <- dismay(count_mat - min(count_mat) + eps, metric) 
    } else {
        sim_mat <- dismay(log_mat, metric) 
    }
    sim_mat <- sim_mat - min(sim_mat) + eps # make all positive with eps for numerical stability
    d_mat <- 1 / (sim_mat) # convert to distance matrix
    diag(d_mat) <- 0 # distance from gene to itself = 0
    print(glue("DIMES: built graph for {metric}"))
    as_tbl_graph(d_mat)
}

