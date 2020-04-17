library(tidygraph)
library(ggplot2)
library(ggraph)
library(ggdendro)
library(purrr)
library(Matrix)

mmnorm <- function(x, mx, mn) {
    exp((x - mn) / (mx - mn))
}

#' Core function to build a given distance matrix between genes (the rownames of *counts*) using *metric*. 
#' Necessary for all other functions. Can optionally run PCA to reduce dimensionality before calculating
#' similarity metrics. Assumes *counts* is rows of genes with columns of cells. 
#' TODO: user is able to define the metric, want a default metric (NMI) and maybe both flag
#' Look into PCA limits/warning message
build_dist_matrix <- function(cells, counts, num_PCs=20, metric, run_PCA = FALSE, gene_mat = TRUE) {
    FUN <- match.fun(metric)
    genes <- rownames(counts)
    
    # builds similarity matrix via metric
    if (run_PCA) {
        pca <- counts[,cells] %>% prcomp 
        pca$x <- round(pca$x, 0) 
        d <- pca$x[,1:num_PCs]
    } else {
        d <- counts[, cells]
    }
    # reduce dimensionality and discretize
    d_mat <- matrix(nrow=length(genes), ncol=length(genes), dimnames = list(genes, genes))
    for (i in c(1:length(genes))) {
        x <- genes[i]
        for (j in c(i:length(genes))) {
            y <- genes[j]    
            d_mat[x, y] <- 1 - FUN(d[x,], d[y,])
            d_mat[y, x] <- d_mat[x, y] # distance metric should be symmetric
        }
        d_mat[x, x] <- 0
    }
    d_mat[] <- vapply(d_mat, partial(mmnorm, mx = max(d_mat), mn = min(d_mat)), numeric(1))
    d_mat[] 
}

# 2. Functions for saving and visualizing networks ---------------------------------------------------

#' Save picture of network built from *x* colored by authority (Kleinberg Centrality) at *filename*.png
#' TODO: option on build_dist_matrix
save_network <- function(x, filename) {
    mutate(x, centrality = centrality_authority()) %>%
        ggraph(layout = "kk") +
            scale_color_gradientn(colors = wes_palette("Zissou1", 21, type = "continuous")) +
            geom_edge_link(aes(alpha = weight)) +
            scale_edge_alpha(range = c(.1, .5)) +
            geom_node_point(aes(color = centrality), size = 12) + 
            geom_node_text(aes(label = name), fontface = "bold") +
            theme_graph(base_family = "Arial") +
            labs(edge_alpha = "distance") +
            ggtitle(paste(filename, "expression network"))

    ggsave(paste("plots/",filename, ".png",sep=""), width = 11, height = 8.5)
}

#' Save picture of genes colored by cluster membership from *x* at *filename* clusters.png
#' TODO: also option on build_dist_matrix
save_gene_clusters <- function(x, filename) {
    mutate(x, group = as.factor(group_louvain())) %>% 
        ggraph(layout = "kk") +
            geom_node_point(aes(color = group), size = 12) + 
            geom_node_text(aes(label = name), fontface = "bold") +
            theme_graph() +
            ggtitle(paste(filename, "clusters"))

    ggsave(paste("plots/",filename, " clusters.png",sep=""), width = 11, height = 8.5)
}

#' Wrapper function to create cluster and network plots using both mutual information and correlation
create_clusters <- function(dataset_name, mi.graph, cor.graph) {
    save_gene_clusters(mi.graph, filename = paste(dataset_name, "MI"))
    save_network(mi.graph, filename = paste(dataset_name, "MI"))
    save_gene_clusters(cor.graph, filename = paste(dataset_name, "corr"))
    save_network(mi.graph, filename = paste(dataset_name, "corr"))
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
