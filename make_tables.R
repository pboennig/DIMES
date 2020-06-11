library("formattable")
library("dplyr")
library("htmltools")
library("webshot")

export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

get_df <- function(cell_line) {
  df <- paste("aucs/", cell_line, ".csv", sep="") %>% read.csv
  df[3] <- round(df[3], 2)
  df[,-1]
}
make_tables <- function(cell_lines) {
  df <- get_df(cell_lines[1])
  colnames(df)[2] <- cell_lines[1]
  if (length(cell_lines) > 1) {
    for (i in 2:length(cell_lines)) {
      df[, cell_lines[i]] <- get_df(cell_lines[i])$auc
    }
  }
  alignment <- c("l", rep("c", length(cell_lines)))
  formattable(df, align=alignment)
}

breast <- c("MCF7", "SUM159")
export_formattable(make_tables(breast), "aucs/breast.png", width="25%")
liver <- c("HUH1", "HUH7")
export_formattable(make_tables(liver), "aucs/liver.png", width="25%")
crc <- c("A549", "H1437", "HCT116", "IMR90", "GM12878", "H1")
export_formattable(make_tables(crc), "aucs/crc.png", width="50%")
