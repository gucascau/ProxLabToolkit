library(proDA)
library(dplyr)

#' proDA models left-censored missingness directly rather than imputing, so
#' it takes the log2-only matrix (intens_tran from log2_quantile_normalize) -
#' NOT quantile-normalized or imputed - and does its own median_normalization.
run_proda <- function(intens_tran, spec) {
  n1 <- length(spec$group1_cols)
  n2 <- length(spec$group2_cols)
  mat <- intens_tran[, c(spec$group1_cols, spec$group2_cols), drop = FALSE]

  normalized_abundance_matrix <- median_normalization(mat)
  sample_info_df <- data.frame(
    name = colnames(normalized_abundance_matrix),
    condition = factor(c(rep(spec$group1_label, n1), rep(spec$group2_label, n2))),
    replicate = c(seq_len(n1), seq_len(n2))
  )

  fit <- proDA(normalized_abundance_matrix, design = ~condition,
               col_data = sample_info_df, reference_level = spec$group2_label)

  coef_name <- paste0("condition", spec$group1_label)
  test_res <- test_diff(fit, coef_name, sort_by = "pval")

  test_res$logFC <- test_res$diff
  test_res$P.Value <- test_res$pval
  test_res$adj.P.Val <- if ("adj_pval" %in% colnames(test_res)) test_res$adj_pval else NA_real_
  test_res %>% rename(Gene = name)
}
