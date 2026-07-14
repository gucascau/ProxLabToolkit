library(DEP)
library(dplyr)
library(tibble)

#' Run the DEP package's own native workflow (make_se -> normalize_vsn ->
#' impute(knn) -> test_diff -> add_rejections -> get_results) on the RAW
#' (not log2/quantile-normalized) intensity matrix - DEP's vsn normalization
#' does its own variance-stabilizing transform.
#'
#' NOTE (ported as-is from the original scripts, not "fixed"): normalize_vsn()
#' is applied to the unfiltered SE, not the filter_missval() output - the
#' missing-value filter is computed but never actually feeds into what gets
#' normalized/imputed downstream. This matches the original pipeline's
#' behavior; flagging it here in case it's worth revisiting.
#'
#' Result column names from DEP::get_results() are contrast-dependent
#' (e.g. "<levelA>_vs_<levelB>_ratio") and the level order depends on R's
#' default alphabetical factor ordering, not the order groups were supplied
#' in - so those columns are detected by suffix rather than hard-coded, and
#' also exposed as standardized logFC/P.Value/adj.P.Val columns.
run_dep_package <- function(intens_raw, spec, imputation_rowmax = 0.9) {
  n1 <- length(spec$group1_cols)
  n2 <- length(spec$group2_cols)

  deg_intens <- as.data.frame(intens_raw) %>% rownames_to_column("ID")
  deg_intens$name <- deg_intens$ID
  lfq_cols <- 2:(1 + n1 + n2)

  experimental_design <- data.frame(
    label = colnames(deg_intens)[lfq_cols],
    condition = factor(c(rep(spec$group1_label, n1), rep(spec$group2_label, n2))),
    replicate = c(seq_len(n1), seq_len(n2))
  )

  data_se <- make_se(deg_intens, lfq_cols, experimental_design)
  data_norm <- normalize_vsn(data_se)
  data_imp_knn <- DEP::impute(data_norm, fun = "knn", rowmax = imputation_rowmax)
  data_diff <- DEP::test_diff(data_imp_knn, type = "all")
  dep <- add_rejections(data_diff, alpha = 0.05, lfc = 0)
  res <- get_results(dep)

  ratio_col <- grep("_ratio$", colnames(res), value = TRUE)
  pval_col <- grep("_p\\.val$", colnames(res), value = TRUE)
  padj_col <- grep("_p\\.adj$", colnames(res), value = TRUE)
  if (length(ratio_col) != 1 || length(pval_col) != 1) {
    stop("Expected exactly one contrast from DEP::get_results(); found ratio cols: ",
         paste(ratio_col, collapse = ", "))
  }

  res$logFC <- res[[ratio_col]]
  res$P.Value <- res[[pval_col]]
  res$adj.P.Val <- if (length(padj_col) == 1) res[[padj_col]] else NA_real_
  res %>% rename(Gene = name)
}
