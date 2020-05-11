library(tidygraph)
library(ggplot2)
library(ggraph)
library(purrr)
library(Matrix)
library(dismay)
library(pryr)
library(wesanderson) #colors

#' Given seurat object *obj*, cell list *cells* and *metrics*, return *graphs*,
#' where graphs$m is a tbl_graph composed using distance metric m from metrics
build_gene_graphs <- function(obj, cells, metrics, num_genes=70) {
    obj <- FindVariableFeatures(obj)
    genes <- head(VariableFeatures(obj), num_genes) # top 70 most variable genes
    counts <- GetAssayData(object=obj,slot="counts")[genes, ] %>% as.matrix # rows = genes, columns = cells
    cells <- cells[cells %in% colnames(counts)] # only take cells that are actually in our counts
    mat <- t(as.matrix(counts[, cells])) # dismay expects as transpose of Seurat default
    graphs <- vector(mode="list",length=length(metrics))
    names(graphs) <- metrics
    for (metric in metrics) {
        graphs$metric <- dist_matrix(mat, metric) %>% as_tbl_graph(directed=FALSE)
    }
    graphs
}

#' Core function to build a given distance matrix between genes (the rownames of *counts*) using *metric*. 
#' Wrapper around dismay with inversion to convert to distance matrix.
dist_matrix <-function(mat, metric, eps=.01) {
    sim_mat <- dismay(mat, metric) 
    sim_mat <- sim_mat - min(sim_mat) + eps # make all positive with eps for numerical stability
    d_mat <- 1 / (sim_mat) # convert to distance matrix
    diag(d_mat) <- 0 # distance from gene to itself = 0
    d_mat
}

# 2. Functions for saving and visualizing networks ---------------------------------------------------
#' Save picture of network built from *x* colored by authority (Kleinberg Centrality) at *filename*.png
save_network <- function(x, filename, layout_alg="drl") {
    mutate(x, centrality = centrality_authority()) %>%
        ggraph(layout = layout_alg) +
            scale_color_gradientn(colors = wes_palette("Zissou1", 21, type = "continuous")) +
            geom_edge_link(aes(alpha = weight)) +
            scale_edge_alpha(range = c(.1, .5)) +
            geom_node_point(aes(color = centrality), size = 5) + 
            geom_node_text(aes(label = name), repel = TRUE, color = "firebrick3") +
            labs(edge_alpha = "distance") +
            ggtitle(paste(filename, "expression network"))

    ggsave(paste("plots/",filename, ".png",sep=""), width = 11, height = 8.5)
}

#' Save picture of genes colored by cluster membership from *x* at *filename* clusters.png
save_gene_clusters <- function(x, filename) {
    mutate(x, group = as.factor(group_louvain())) %>% 
        ggraph(layout = layout_alg) +
            geom_node_point(aes(color = group), size = 5) + 
            geom_node_text(aes(label = name), repel = TRUE) +
            ggtitle(paste(filename, "clusters"))

    ggsave(paste("plots/",filename, " clusters.png",sep=""), width = 11, height = 8.5)
}

# 3. Functions for listing and ranking genes ----------------------------------
#' Get top *num* most central genes from network *x* as defined by Kleinberg
#' centrality
most_central <- function(x, num, counts) {
    g <- mutate(x, centrality = centrality_authority()) %>%
        top_n(n = num, wt = centrality) %>%
        as_tibble %>%
        select("name") %>%
        mutate(mean_expr = sapply(pull(., name), function(x) return(mean(counts[x, ])))) %>%
        mutate(status = if_else(mean_expr > mean(counts), "oncogene", "TS"))
}

#' Get top *num* genes as ranked by hierarchical clustering
most_hclust <- function(x, num, counts) {
    g <- mutate(x, n_rank_hclust = node_rank_hclust()) %>%
        top_n(n = num, wt = n_rank_hclust) %>%
        as_tibble %>%
        select("name") %>%
        mutate(mean_expr = sapply(pull(., name), function(x) return(mean(counts[x, ])))) %>%
        mutate(status = if_else(mean_expr > mean(counts), "oncogene", "TS"))
}

most_tsp <- function(x, num, counts) {
    g <- mutate(x, n_rank_tsp = node_rank_traveller()) %>%
        top_n(n = num, wt = n_rank_tsp) %>%
        as_tibble %>%
        select("name") %>%
        mutate(mean_expr = sapply(pull(., name), function(x) return(mean(counts[x, ])))) %>%
        mutate(status = if_else(mean_expr > mean(counts), "oncogene", "TS"))
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
