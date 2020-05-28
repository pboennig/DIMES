library(dplyr)
library(STRINGdb)
library(igraph)
library(plotROC)
source("dimes.R")

source("load_adeno.R") #load adeno data

metrics <- c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 'MI',
  'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan',
  'weighted_rank', 'hamming', 'rho_p', 'phi_s')
graphs <- build_gene_graphs(lung, cells, metrics)
# mutinfo functions remove '-' characters
pears_names <- graphs$pearson %>% pull(name)
graphs$MI <- graphs$MI %>% activate(nodes) %>% mutate(name = pears_names)

string_db <- STRINGdb$new(version="11", score_threshold=50, input_directory=".", species=9606)

verts <- graphs[[1]] %>% activate(nodes) %>% as.data.frame # get vertices names
mapped <- string_db$map(verts, "name", removeUnmappedRows = TRUE) # map to string_ids
to_remove <- verts$name[which(!(verts$name %in% mapped$name))] # remove unmapped vertices
mapped <- mapped[mapped$name %in% intersect(verts$name, mapped$name),] 
mapped <- mapped[!duplicated(mapped$name),] # some datasets have duplicates
string_subset <- string_db$get_subnetwork(mapped$STRING_id)
string_subset <- delete.edges(string_subset, which(E(string_subset)$experiments < 400)) %>%
                 as_adjacency_matrix
string_subset <- string_subset / 2# adjacency matrix has 2 if edge exists, want it to be 0/1

calc_score <- function(g, mapped, string_subset, to_remove) {
  g <- as.igraph(g) %>%
    delete.vertices(to_remove)
  g <- set.vertex.attribute(g, "name", value=mapped[which(mapped$name==V(g)$name),]$STRING_id)
  score <- get.adjacency(g, attr="weight") # higher score -> closer distance
  correspondence <- match(rownames(string_subset), rownames(score))
  score <- score[correspondence,]
  score[,correspondence] %>% as.vector 
}
scores <- lapply(graphs, calc_score, mapped, string_subset, to_remove)

names(scores) <- metrics
s <- data.frame(scores)
s$status <- as.vector(string_subset) # score matrix w/ labels in status column
longtest <- melt_roc(s, d = "status", m = metrics) # merge into long data (all metrics together in one column)
p <- ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(labels=FALSE, increasing=FALSE) + style_roc() # plot ROCs
aucs <- data.frame(metric=metrics, auc=calc_auc(p)$AUC)
aucs <- aucs[order(-aucs$auc),] #calculate AUCs and sort by AUC value
