library(dplyr)
library(STRINGdb)
library(igraph)
source("dimess.R")

source("load_adeno.R") #load adeno data

metrics <- c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 'MI',
  'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
  'weighted_rank', 'hamming')
graphs <- build_gene_graphs(lung, cells, metrics)

string_db <- STRINGdb$new(version="11")
mapped <- graphs[[1]] %>% 
            activate(nodes) %>% 
            as.data.frame %>%
            string_db$map("name", removeUnmappedRows = TRUE)

string_subset <- string_db$get_subnetwork(mapped$STRING_id) %>% as_adjacency_matrix
calc_acc <- function(m, graphs, mapped, string_subset) {
  num_connections <- sum(string_subset / 2) # how many PPI connections
  adj_mat <- as.igraph(graphs[[m]]) %>%
          simplify %>% # remove self-edges
          delete.vertices(setdiff(V(.)$name, mapped$name)) %>%
          set.vertex.attribute("name", value = mapped[which(mapped$name == V(.)$name),]$STRING_id) %>%
          delete.edges(E(.)[E(.)$weight > sort(E(.)$weight)[num_connections]]) %>%
          as_adjacency_matrix
  1 - norm(adj_mat - string_subset / 2) / num_connections
}
accs <- lapply(metrics, calc_acc, graphs, mapped, string_subset)
names(accs) <- metrics