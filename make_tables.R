library("gt")
library("dplyr")
library("purrr")
library("glue")
library("htmltools")
library("webshot")

make_tables <- function(cancer_type, cell_lines) {
  df <- get_df(cell_lines[1])
  colnames(df)[2] <- cell_lines[1]
  if (length(cell_lines) > 1) {
    for (i in 2:length(cell_lines)) {
      df[, cell_lines[i]] <- get_df(cell_lines[i])$auc
    }
  }
  gt(df) %>% 
    tab_header(title = cancer_type) %>% 
    gtsave(glue("aucs/{cancer_type}.png"))
}

cellLineMap <- list()
cellLineMap[["Breast"]] <- c("MCF7", "SUM159")
cellLineMap[["Liver"]] <- c("HUH1", "HUH7")
cellLineMap[["CRC"]] <- c("A549", "H1437", "HCT116", "IMR90", "GM12878", "H1")
cellLineMap[["Leukemia"]] <- c("K652")
map2(names(cellLineMap), cellLineMap, make_tables)
