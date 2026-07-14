library(dplyr)

# Small placeholder list of ubiquitous lab contaminants (keratins, trypsin,
# serum albumin). This is NOT a CRAPome substitute - real projects should
# supply their own CRAPome export (or equivalent) via load_contaminant_list().
.default_contaminants <- c(
  "KRT1", "KRT2", "KRT9", "KRT10", "KRT14", "KRT16",
  "PRSS1", "PRSS2", "ALB"
)

#' Load a contaminant list from a user-supplied CSV/TSV (one gene symbol per
#' row, in `gene_col`). Falls back to a small built-in placeholder list if
#' `path` is NULL.
load_contaminant_list <- function(path = NULL, gene_col = "Gene") {
  if (is.null(path)) {
    warning("No contaminant list supplied - using a small built-in placeholder ",
            "(keratins/trypsin/albumin only). Supply a real CRAPome export via ",
            "`path` for production use.")
    return(.default_contaminants)
  }
  tbl <- if (grepl("\\.tsv$", path)) {
    read.delim(path, stringsAsFactors = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE)
  }
  if (!gene_col %in% colnames(tbl)) {
    stop("Column '", gene_col, "' not found in ", path)
  }
  unique(toupper(trimws(tbl[[gene_col]])))
}

#' Remove rows whose gene symbol matches the contaminant list. Returns the
#' filtered data frame; prints how many rows were removed.
filter_contaminants <- function(df, contaminants, gene_col = "Gene") {
  is_contam <- toupper(df[[gene_col]]) %in% toupper(contaminants)
  message(sum(is_contam), " of ", nrow(df), " proteins removed as contaminants.")
  df[!is_contam, , drop = FALSE]
}
