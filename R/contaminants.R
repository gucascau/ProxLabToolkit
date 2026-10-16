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

#' Load a full CRAPome matrix export (e.g. Datasets/crapome_reprint-apms.org.csv:
#' one row per gene, a few annotation columns, then one spectral-count column
#' per control experiment named CC<n>) and summarise it per gene, the way the
#' CRAPome web tool does: CRAPome_n_found (control runs with SC > 0),
#' CRAPome_Freq (n_found / total runs), CRAPome_AvgSC (mean SC over the runs
#' it was found in) and CRAPome_MaxSC. Gene symbols are upper-cased; rows
#' sharing a symbol are collapsed by taking the max SC per experiment.
#' Unlike load_contaminant_list(), this keeps the frequencies rather than
#' treating every listed gene as a contaminant - almost every expressed
#' protein appears in some CRAPome control.
load_crapome <- function(path, gene_col = "geneSymbol", count_pattern = "^CC[0-9]+$") {
  tbl <- if (grepl("\\.tsv$", path)) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!gene_col %in% colnames(tbl)) {
    stop("Column '", gene_col, "' not found in ", path)
  }
  count_cols <- grep(count_pattern, colnames(tbl), value = TRUE)
  if (length(count_cols) == 0) {
    stop("No control-experiment columns matching '", count_pattern, "' in ", path)
  }
  mat <- as.matrix(tbl[, count_cols])
  mat[is.na(mat)] <- 0
  genes <- toupper(trimws(tbl[[gene_col]]))
  keep <- !is.na(genes) & genes != ""
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]

  dup <- genes %in% genes[duplicated(genes)]
  if (any(dup)) {
    collapsed <- lapply(split(which(dup), genes[dup]),
                        function(i) apply(mat[i, , drop = FALSE], 2, max))
    mat <- rbind(mat[!dup, , drop = FALSE], do.call(rbind, collapsed))
    genes <- c(genes[!dup], names(collapsed))
  }

  n_found <- rowSums(mat > 0)
  out <- data.frame(
    Gene = genes,
    CRAPome_n_found = n_found,
    CRAPome_Freq = n_found / ncol(mat),
    CRAPome_AvgSC = ifelse(n_found > 0, rowSums(mat) / pmax(n_found, 1), 0),
    CRAPome_MaxSC = apply(mat, 1, max),
    stringsAsFactors = FALSE, row.names = NULL
  )
  message("CRAPome: ", nrow(out), " genes across ", ncol(mat), " control experiments.")
  out
}

#' Join CRAPome summary columns (from load_crapome()) onto a result table by
#' gene symbol, plus a logical CRAPome_Flag (CRAPome_Freq >= freq_cutoff).
#' Genes absent from CRAPome get n_found/Freq/SC = 0 and are never flagged.
#' Any existing CRAPome_* columns are replaced, so re-annotating is safe.
annotate_crapome <- function(df, crapome, gene_col = "Gene", freq_cutoff = 0.2) {
  df <- df[, !grepl("^CRAPome_", colnames(df)), drop = FALSE]
  idx <- match(toupper(df[[gene_col]]), crapome$Gene)
  for (col in c("CRAPome_n_found", "CRAPome_Freq", "CRAPome_AvgSC", "CRAPome_MaxSC")) {
    vals <- crapome[[col]][idx]
    vals[is.na(vals)] <- 0
    df[[col]] <- vals
  }
  df$CRAPome_Flag <- df$CRAPome_Freq >= freq_cutoff
  df
}
