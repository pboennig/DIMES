library("gt")
library("tidyverse")
library("glue")
library("purrr")
library("ggplot2")

make_tibbles <- function(cell_lines) {
  tbl <- read_csv(glue("aucs/{cell_lines[1]}.csv"), col_types = "_cd")
  colnames(tbl)[2]  <- cell_lines[1]
  if (length(cell_lines) > 1) {
    for (i in 2:length(cell_lines)) {
      tbl[, cell_lines[i]] <- read.csv(glue("aucs/{cell_lines[i]}.csv", col_types="_cd"))$auc
    }
  }
  tbl
}

save_tables <- function(tbl, cancer_type) {
  gt(data = tbl) %>% 
    tab_header(title = cancer_type) %>%
    gtsave(glue("tables/{cancer_type}.png"))
}

plot_long_tbl <- function(long_tbl, title, colorCol) {
  colorCol <- sym(colorCol)
  ggplot(long_tbl, aes(x=metric, y = AUC)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color= !! colorCol), width=.2) + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    ggtitle(glue("{title} AUCS"))
}

plot_cancer <- function(tbl, cancer_type) {
  long_tbl <- tbl %>% pivot_longer(cols=-metric, names_to="CellLine", values_to="AUC")
  long_tbl$metric <- reorder(long_tbl$metric, -long_tbl$AUC, median)
  plot_long_tbl(long_tbl, cancer_type, "CellLine")
  ggsave(glue("plots/{cancer_type} plot.png"))
}

overall_plot <- function(tibbles) {
  long_tbls <- lapply(tibbles, pivot_longer, cols=-metric, names_to="CellLine", values_to="AUC")
  for (i in 1:length(long_tbls)) {
    long_tbls[[i]]$CancerType <- rep(names(long_tbls)[i], nrow(long_tbls[[i]]))
  }
  combined <- bind_rows(long_tbls)
  combined$metric <- reorder(combined$metric, -combined$AUC, median)
  plot_long_tbl(combined, "All cancers by cell line", "CellLine")
  ggsave(glue("plots/combined plot CellLine.png"))
  plot_long_tbl(combined, "All cancers by type", "CancerType")
  ggsave(glue("plots/combined plot CancerType.png"))
}

summarize_AUC <- function(tbl) {
  tbl %>%
    group_by(metric) %>%
    summarise(CancerType = CancerType,
              mean_AUC = mean(AUC),
              sd_AUC = sd(AUC),
              n_AUC = n(),
              SE_AUC = sd(AUC) / sqrt(n()),
              .groups = "drop") %>%
    distinct
}
  

mean_cancer_type <- function(tibbles, numCellLines) {
  dodge_w <- .5
  long_tbls <- lapply(tibbles, pivot_longer, cols=-metric, names_to="CellLine", values_to="AUC")
  for (i in 1:length(long_tbls)) {
    long_tbls[[i]]$CancerType <- rep(names(long_tbls)[i], nrow(long_tbls[[i]]))
  }
  combined <- lapply(long_tbls, summarize_AUC)
  combined <- bind_rows(combined)
  combined$metric <- reorder(combined$metric, -combined$mean_AUC, median)
  combined$NumberCellLines <- numCellLines[combined$CancerType]
  ggplot(combined, aes(x=metric, y = mean_AUC, color = CancerType)) + 
    geom_point(position = position_dodge(width=dodge_w), aes(size=NumberCellLines)) +
    scale_size(range = c(1,3), name = "Number of cell lines", breaks=c(1,2,6,9)) +
    theme_classic() +
    geom_errorbar(
      aes(ymin = mean_AUC-SE_AUC, ymax = mean_AUC + SE_AUC),
      width = 0.1,
      linetype = "solid",
      position=position_dodge(width=dodge_w)) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
  ggsave(glue("plots/combined plot mean CancerType.png"), width = 10, height = 6)
}

cellLineMap <- list()
cellLineMap[["Breast"]] <- c("MCF7", "SUM159")
cellLineMap[["Liver"]] <- c("HUH1", "HUH7")
cellLineMap[["CRC"]] <- c("A549", "H1437", "HCT116", "IMR90", "GM12878", "H1")
cellLineMap[["Leukemia"]] <- c("K652")
cellLineMap[["Melanoma"]] <- c("A375", "451Lu")
cellLineMap[["Gastric"]] <- c("HGC-27", "KATOIII", "MKN-45", "NCI-N87", "NUGC-4",
                              "SNU-16", "SNU-601",  "SNU-638",  "SNU-668")
numCellLines <- sapply(cellLineMap, length)
tibbles <- lapply(cellLineMap, make_tibbles)
#map2(tibbles, names(cellLineMap), save_tables)
#map2(tibbles, names(cellLineMap), plot_cancer)
#overall_plot(tibbles)
mean_cancer_type(tibbles, numCellLines)
