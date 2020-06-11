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

#' Calculates the AUCs for different *metrics* on a particular Seurat *obj* by comparing
#' to STRINGdb network
ppi_comp <- function(obj, cells=colnames(obj), num_genes=200, metrics=c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 'MI',
                                                          'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
                                                          'weighted_rank', 'hamming', 'rho_p', 'phi_s')) {
  genes <- data.frame(gene=head(VariableFeatures(obj), num_genes)) # top 70 most variable genes
  string_db <- STRINGdb$new(version="11", score_threshold=50, input_directory=".", species=9606)
  mapped <- string_db$map(genes, "gene", removeUnmappedRows=TRUE)
  graphs <- build_gene_graphs(obj, cells, mapped$gene, metrics)
  
  # mutinfo functions remove '-' characters
  if ('MI' %in% metrics) {
    pears_names <- graphs$pearson %>% pull(name)
    graphs$MI <- graphs$MI %>% activate(nodes) %>% mutate(name = pears_names)
  }
  
  string_subset <- string_db$get_subnetwork(mapped$STRING_id)
  string_subset <- delete.edges(string_subset, which(E(string_subset)$experiments < 400)) %>%
                   as_adjacency_matrix
  string_subset <- string_subset / 2# adjacency matrix has 2 if edge exists, want it to be 0/1
  
  scores <- lapply(graphs, calc_score, mapped, string_subset)
  
  names(scores) <- metrics
  s <- data.frame(scores)
  s$status <- as.vector(string_subset) # score matrix w/ labels in status column
  longtest <- melt_roc(s, d = "status", m = metrics) # merge into long data (all metrics together in one column)
  p <- ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(labels=FALSE, increasing=FALSE) + style_roc() # plot ROCs
  aucs <- data.frame(metric=metrics, auc=calc_auc(p)$AUC)
<<<<<<< HEAD:dimes.R
  rm(list=setdiff(ls(), "aucs"))
=======
  aucs <- aucs[order(-aucs$auc),] #calculate AUCs and sort by AUC value
>>>>>>> 166dbc54a8b84aabb0b8e68bc47e7ffaf4dbfc31:dimess.R
  aucs
}

#' Calculate scores for AUC calculation
calc_score <- function(g, mapped, string_subset) {
    g <- as.igraph(g)
    g <- set.vertex.attribute(g, "name", value=mapped[which(mapped$gene==V(g)$name),]$STRING_id)
    score <- get.adjacency(g, attr="weight") # higher score -> closer distance
    correspondence <- match(rownames(string_subset), rownames(score))
    score <- score[correspondence,]
    score[,correspondence] %>% as.vector 
}

#' Given seurat object *obj*, cell list *cells* and *metrics*, return *graphs*,
#' where graphs$m is a tbl_graph composed using distance metric m from metrics
build_gene_graphs <- function(obj, cells, genes, metrics) {
    counts <- GetAssayData(object=obj,slot="counts")[genes, ] %>% as.matrix # rows = genes, columns = cells
    log_counts <- GetAssayData(object=obj,slot="scale.data")[genes, ] %>% as.matrix # rows = genes, columns = cells
    cells <- cells[cells %in% colnames(counts)] # only take cells that are actually in our counts
    count_mat <- t(as.matrix(counts[, cells])) # dismay expects as transpose of Seurat default
    log_mat <- t(as.matrix(log_counts[, cells])) # dismay expects as transpose of Seurat default
    graphs <- lapply(metrics, dist_matrix, count_mat, log_mat)
    names(graphs) <- metrics
    graphs
}

#' Wrapper around dismay with inversion to convert to distance matrix.
dist_matrix <-function(metric, count_mat, log_mat, eps=.001) {
    if (metric %in% c('phi_s', 'rho_p')) {
        sim_mat <- dismay(count_mat, metric) 
    } else {
        sim_mat <- dismay(log_mat, metric) 
    }
    sim_mat <- sim_mat - min(sim_mat) + eps # make all positive with eps for numerical stability
    d_mat <- 1 / (sim_mat) # convert to distance matrix
    diag(d_mat) <- 0 # distance from gene to itself = 0
    as_tbl_graph(d_mat)
}

# 2. Functions for saving and visualizing networks ---------------------------------------------------
#' Save picture of network built from *x* colored by authority (Kleinberg Centrality) at *filename*.png
plot_network <- function(x, filename, layout_alg="drl") {
     p <- mutate(x, centrality = centrality_authority()) %>%
        ggraph(layout = layout_alg) +
            scale_color_gradientn(colors = wes_palette("Zissou1", 21, type = "continuous")) +
            geom_edge_link(aes(alpha = weight)) +
            scale_edge_alpha(range = c(.1, .5)) +
            geom_node_point(aes(color = centrality), size = 5) + 
            geom_node_text(aes(label = name), repel = TRUE, color = "firebrick3") +
            labs(edge_alpha = "distance") + 
            ggtitle(paste(filename, "expression network"))
    print(p)
}

#' Plot picture of genes colored by cluster membership from *x* at *filename* clusters.png
plot_gene_clusters <- function(x, filename, layout_alg="drl") {
    p <- mutate(x, group = as.factor(group_louvain())) %>% 
        ggraph(layout = layout_alg) +
            geom_node_point(aes(color = group), size = 5) + 
            geom_node_text(aes(label = name), repel = TRUE) + 
            ggtitle(paste(filename, "clusters"))
    print(p)
}

# 3. Functions for listing and ranking genes ----------------------------------
#' Get top *num* most central genes from network *x* as defined by Kleinberg
#' centrality
most_central <- function(x, num) {
    g <- mutate(x, centrality = centrality_authority()) %>%
        top_n(n = num, wt = centrality) %>%
        as_tibble %>%
        select("name")
}

#' Get top *num* genes as ranked by hierarchical clustering
most_hclust <- function(x, num) {
    g <- mutate(x, n_rank_hclust = node_rank_hclust()) %>%
        top_n(n = num, wt = n_rank_hclust) %>%
        as_tibble %>%
        select("name")
}

most_tsp <- function(x, num, counts) {
    g <- mutate(x, n_rank_tsp = node_rank_traveller()) %>%
        top_n(n = num, wt = n_rank_tsp) %>%
        as_tibble %>%
        select("name")
}

#' Print list of chosen genes as oncogene or tumor supressor
print_list <- function(name, genes, counts) {
    print(name)
    med <- mean(counts)
    for (i in 1:10) {
        gene <- genes[i]
        if (mean(counts[gene, ]) > med) {
            type <- "Oncogene"
        } else {
            type <- "Tumor Supressor"
        }
        print(paste(i, gene, type))
    }
}
