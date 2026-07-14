library(dplyr)
library(purrr)
library(metap)

#' Join every method's result table into one wide table, prefixing each
#' method's columns with its name (e.g. "Limma_logFC") and keyed by a single
#' `ID` column (Gene symbol - every run_* wrapper in this toolkit already
#' standardizes on a "Gene" column). Ports PvalueIntegration.Rmd's join logic
#' (inner_join across all methods - a protein must be callable by every
#' method to enter the consensus).
join_method_results <- function(results_list, id_col = "Gene") {
  renamed <- imap(results_list, function(df, method_name) {
    # some methods (e.g. DEP's get_results()) keep their own "ID" column
    # distinct from the Gene/name column we join on - drop it first so it
    # doesn't collide with the "ID" key column created below.
    if ("ID" %in% colnames(df) && id_col != "ID") {
      df <- df %>% select(-ID)
    }
    df <- df %>% rename(ID = all_of(id_col))
    other_cols <- setdiff(colnames(df), "ID")
    df %>% rename_with(~ paste0(method_name, "_", .x), all_of(other_cols))
  })
  reduce(renamed, function(x, y) inner_join(x, y, by = "ID"))
}

#' Pick, per protein, the logFC with the largest absolute value among the
#' given columns (ties/NA/Inf ignored). Ports PvalueIntegration.Rmd lines 126-133.
add_max_logfc <- function(combined, logfc_cols) {
  combined %>%
    rowwise() %>%
    mutate(max_logFC = {
      vals <- c_across(all_of(logfc_cols))
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) NA_real_ else vals[which.max(abs(vals))]
    }) %>%
    ungroup()
}

#' Combine p-values across methods four ways (Fisher, Tippett, Stouffer, an
#' empirical permutation test), each BH-adjusted. Ports PvalueIntegration.Rmd
#' lines 150-221 verbatim (same formulas/seed/n_perm default).
add_combined_pvalues <- function(combined, pvalue_cols, n_perm = 100000, seed = 123) {
  pmat <- combined[, pvalue_cols, drop = FALSE]

  combined$FinalP_empirical <- apply(pmat, 1, function(x) {
    pvals <- x[!is.na(x)]
    observed <- sum(-2 * log(pvals))
    set.seed(seed)
    perm_stats <- replicate(n_perm, sum(-2 * log(runif(length(pvals)))))
    mean(perm_stats >= observed)
  })
  combined$FinalP_Fisher <- apply(pmat, 1, function(x) sumlog(x[!is.na(x)])$p)
  combined$FinalP_Tippett <- apply(pmat, 1, function(x) 1 - (1 - min(x[!is.na(x)]))^length(x[!is.na(x)]))
  combined$FinalP_Stouffer <- apply(pmat, 1, function(x) sumz(x[!is.na(x)])$p)

  combined$FinaladjP_Fisher <- p.adjust(combined$FinalP_Fisher, method = "BH")
  combined$FinaladjP_Tippett <- p.adjust(combined$FinalP_Tippett, method = "BH")
  combined$FinaladjP_Stouffer <- p.adjust(combined$FinalP_Stouffer, method = "BH")
  combined$FinaladjP_empirical <- p.adjust(combined$FinalP_empirical, method = "BH")

  combined %>% arrange(FinalP_Fisher)
}

#' Full consensus step: join all methods, pick max_logFC, combine p-values.
#' Every run_* wrapper in this toolkit standardizes on logFC/P.Value columns,
#' so by default logfc_cols/pvalue_cols are just "<method>_logFC" /
#' "<method>_P.Value" for every method in results_list - override to match
#' the original scripts' choice of excluding a method from the combined
#' statistic (the original excluded Wilcoxon/Annova from p-value combination
#' for reasons not stated in the code; this default includes every method
#' passed in instead, since Ttest/Annova are no longer first-class methods
#' here - worth a deliberate decision rather than silently inheriting that).
combine_pvalues <- function(results_list, id_col = "Gene",
                             logfc_cols = paste0(names(results_list), "_logFC"),
                             pvalue_cols = paste0(names(results_list), "_P.Value"),
                             n_perm = 100000, seed = 123) {
  combined <- join_method_results(results_list, id_col = id_col)
  combined <- add_max_logfc(combined, logfc_cols)
  add_combined_pvalues(combined, pvalue_cols, n_perm = n_perm, seed = seed)
}

#' The p<=0.05 & abs(logFC)>0.25 threshold used throughout the project.
call_deps <- function(combined, p_col = "FinalP_Fisher", logfc_col = "max_logFC",
                        p_cutoff = 0.05, logfc_cutoff = 0.25) {
  combined %>% filter(.data[[p_col]] <= p_cutoff & abs(.data[[logfc_col]]) > logfc_cutoff)
}
