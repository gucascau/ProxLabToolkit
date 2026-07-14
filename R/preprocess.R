library(limma)

#' Select one comparison's MaxLFQ.Intensity columns out of a loaded FragPipe
#' protein table, zero -> NA, and drop proteins with more than spec$max_na
#' missing values across the comparison's samples. Returns a numeric matrix
#' (genes x samples), columns ordered group1 then group2.
extract_intensity_matrix <- function(pro_file_update, spec) {
  sample_cols <- c(spec$group1_cols, spec$group2_cols)
  maxlfq_cols <- paste0(sample_cols, ".MaxLFQ.Intensity")
  missing <- setdiff(maxlfq_cols, colnames(pro_file_update))
  if (length(missing) > 0) {
    stop("Missing MaxLFQ.Intensity columns: ", paste(missing, collapse = ", "))
  }

  intens <- as.data.frame(pro_file_update[, maxlfq_cols])
  rownames(intens) <- pro_file_update$Gene
  intens[intens == 0] <- NA
  colnames(intens) <- gsub("\\.MaxLFQ\\.Intensity$", "", colnames(intens))

  keep <- rowSums(is.na(intens)) <= spec$max_na
  as.matrix(intens[keep, sample_cols, drop = FALSE])
}

#' log2-transform then quantile-normalize (limma::normalizeBetweenArrays).
#' Returns both the log2-only matrix (still has NAs; proDA wants this one)
#' and the quantile-normalized matrix (still has NAs; feeds imputation).
log2_quantile_normalize <- function(mat) {
  intens_tran <- log2(mat)
  intens_tran_normalized <- normalizeBetweenArrays(intens_tran, method = "quantile")
  list(intens_tran = intens_tran, intens_tran_normalized = intens_tran_normalized)
}
