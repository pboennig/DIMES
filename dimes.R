source("ipak.R")
## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})
packages <- c("tidygraph",
"ggplot2",
"ggraph",
"purrr",
"Matrix",
"dismay",
"pryr",
"wesanderson",
"tibble",
"dplyr",
"WGCNA",
"STRINGdb",
"igraph",
"plotROC",
"glue")
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
library(plotROC)
library(glue)

#' Calculates the AUCs for different *metrics* on a particular Seurat *obj* by comparing
#' to STRINGdb network
ppi_comp <- function(obj, cells=colnames(obj), num_genes= nrow(obj), metrics=c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 'MI',
                                                          'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
                                                          'weighted_rank', 'hamming', 'rho_p', 'phi_s')) {
  if (num_genes < nrow(obj)) {
    genes <- data.frame(gene=head(VariableFeatures(obj), num_genes)) 
  } else {
    genes <-data.frame(gene=VariableFeatures(obj))
  }
  string_db <- STRINGdb$new(version="11", score_threshold=50, input_directory="STRINGdb/", species=9606)
  mapped <- string_db$map(genes, "gene", removeUnmappedRows=TRUE)
  graphs <- build_gene_graphs(obj, cells, mapped, metrics)
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
  print("DIMES: Scores calculated")  
  
  string_subset <- string_subset[rownames(scores[[1]]), ]
  string_subset <- string_subset[, colnames(scores[[1]])]
  names(scores) <- metrics
  s <- data.frame(lapply(scores, as.vector))
  s$status <- as.vector(string_subset) # score matrix w/ labels in status column
  longtest <- melt_roc(s, d = "status", m = metrics) # merge into long data (all metrics together in one column)
  p <- ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(labels=FALSE, increasing=FALSE) + style_roc() # plot ROCs
  aucs <- data.frame(metric=metrics, auc=calc_auc(p)$AUC)
  print("DIMES: Cleaning up...")
  rm(list=c("string_subset", "scores", "mapped", "string_db")) 
  gc()
  aucs
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
build_gene_graphs <- function(obj, cells, mapped, metrics) {
    counts <- GetAssayData(object=obj,slot="counts")
    log_counts <- GetAssayData(object=obj,slot="scale.data")
    mapped <- mapped[mapped$gene %in% rownames(counts),]
    counts <- counts[mapped$gene,] %>% as.matrix
    log_counts <- log_counts[mapped$gene,] %>% as.matrix
    
    cells <- intersect(cells, colnames(counts)) # only take cells that are actually in our counts
    count_mat <- t(counts[, cells]) # dismay expects as transpose of Seurat default
    log_mat <- t(log_counts[, cells]) # dismay expects as transpose of Seurat default
    graphs <- lapply(metrics, dist_matrix, count_mat, log_mat)
    names(graphs) <- metrics
    graphs
}

#' Wrapper around dismay, switch on metrics that need raw counts
dist_matrix <-function(metric, count_mat, log_mat, eps=.001) {
    if (metric %in% c('phi_s', 'rho_p')) {
        sim_mat <- dismay(count_mat, metric) 
    } else {
        sim_mat <- dismay(log_mat, metric) 
    }
    sim_mat <- sim_mat - min(sim_mat) + eps # make all positive with eps for numerical stability
    d_mat <- 1 / (sim_mat) # convert to distance matrix
    diag(d_mat) <- 0 # distance from gene to itself = 0
    print(glue("DIMES: built graph for {metric}"))
    as_tbl_graph(d_mat)
}

