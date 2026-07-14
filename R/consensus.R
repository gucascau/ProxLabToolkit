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
#'
#' Caveat: all four assume the input p-values are independent tests. Here
#' they mostly aren't - Limma/DEqMS/Wilcoxon/ROTS/DEP/proDA are computed on
#' the same underlying intensity matrix, so their p-values are correlated.
#' Under positive correlation, Fisher's method in particular is
#' anti-conservative (inflates the false-positive rate). See
#' `add_corrected_fisher()` for a correlation-adjusted alternative to
#' `FinalP_Fisher`, and `add_vote_counts()`/`call_deps_vote()` for an
#' approach that doesn't assume independence at all.
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

#' Correlation-corrected Fisher combination (Kost & McDermott 2002): adjusts
#' plain Fisher (`FinalP_Fisher` above) for the fact that the methods being
#' combined share the same underlying data and so aren't independent tests.
#' The average pairwise correlation among methods' -2*log(p) statistics is
#' estimated once, globally, across all proteins (a simplification - true
#' correlation likely varies by protein/method-pair, but a single global
#' estimate is the standard practical approach and avoids needing a
#' per-protein correlation, which isn't estimable from one row of data).
#' That correlation then rescales Fisher's chi-square statistic per row via
#' Kost & McDermott's cubic approximation (valid for average correlation in
#' [0, 0.9], which the estimate is clamped to).
add_corrected_fisher <- function(combined, pvalue_cols) {
  pmat <- as.matrix(combined[, pvalue_cols, drop = FALSE])
  stat_mat <- -2 * log(pmat)

  cor_mat <- suppressWarnings(cor(stat_mat, use = "pairwise.complete.obs"))
  r <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
  if (!is.finite(r)) r <- 0
  r <- max(0, min(r, 0.9))

  combined$FinalP_FisherCorrected <- apply(pmat, 1, function(x) {
    x <- x[!is.na(x)]
    k <- length(x)
    if (k == 0) return(NA_real_)
    if (k == 1) return(x)
    stat <- sum(-2 * log(x))
    cov_ij <- r * (3.263 + 0.710 * r + 0.027 * r^2)
    var_t <- 4 * k + 2 * choose(k, 2) * cov_ij
    e_t <- 2 * k
    f <- 2 * e_t^2 / var_t
    c_scale <- var_t / (2 * e_t)
    pchisq(stat / c_scale, df = f, lower.tail = FALSE)
  })
  combined$FinaladjP_FisherCorrected <- p.adjust(combined$FinalP_FisherCorrected, method = "BH")
  combined
}

#' Vote-counting alternative to p-value combination: per protein, count how
#' many methods it was testable in (`n_methods_tested`) and how many of
#' those calls clear that method's own p_cutoff/logfc_cutoff
#' (`n_methods_hit`). Unlike Fisher/Tippett/Stouffer, this makes no
#' independence assumption about the methods, at the cost of losing a single
#' interpretable combined p-value/effect size. See `call_deps_vote()`.
add_vote_counts <- function(combined, logfc_cols, pvalue_cols,
                             p_cutoff = 0.05, logfc_cutoff = 0.25) {
  stopifnot(length(logfc_cols) == length(pvalue_cols))
  hits <- Map(function(lfc_col, p_col) {
    !is.na(combined[[lfc_col]]) & !is.na(combined[[p_col]]) &
      abs(combined[[lfc_col]]) > logfc_cutoff & combined[[p_col]] <= p_cutoff
  }, logfc_cols, pvalue_cols)
  tested <- lapply(pvalue_cols, function(p_col) !is.na(combined[[p_col]]))

  combined$n_methods_tested <- Reduce(`+`, tested)
  combined$n_methods_hit <- Reduce(`+`, hits)
  combined
}

#' Full consensus step: join all methods, pick max_logFC, combine p-values
#' (plain + correlation-corrected Fisher), and tally per-protein vote counts.
#' Every run_* wrapper in this toolkit standardizes on logFC/P.Value columns,
#' so by default logfc_cols/pvalue_cols are just "<method>_logFC" /
#' "<method>_P.Value" for every method in results_list - override to match
#' the original scripts' choice of excluding a method from the combined
#' statistic (the original excluded Wilcoxon/Annova from p-value combination
#' for reasons not stated in the code; this default includes every method
#' passed in instead, since Ttest/Annova are no longer first-class methods
#' here - worth a deliberate decision rather than silently inheriting that).
#' `p_cutoff`/`logfc_cutoff` are only used for the vote-count columns
#' (n_methods_tested/n_methods_hit) - they don't affect the p-value columns.
combine_pvalues <- function(results_list, id_col = "Gene",
                             logfc_cols = paste0(names(results_list), "_logFC"),
                             pvalue_cols = paste0(names(results_list), "_P.Value"),
                             n_perm = 100000, seed = 123,
                             p_cutoff = 0.05, logfc_cutoff = 0.25) {
  combined <- join_method_results(results_list, id_col = id_col)
  combined <- add_max_logfc(combined, logfc_cols)
  combined <- add_combined_pvalues(combined, pvalue_cols, n_perm = n_perm, seed = seed)
  combined <- add_corrected_fisher(combined, pvalue_cols)
  add_vote_counts(combined, logfc_cols, pvalue_cols,
                   p_cutoff = p_cutoff, logfc_cutoff = logfc_cutoff)
}

#' The p<=0.05 & abs(logFC)>0.25 threshold used throughout the project.
#' Defaults to plain Fisher for continuity with prior results; pass
#' p_col = "FinaladjP_FisherCorrected" for the correlation-adjusted version
#' (recommended - see add_corrected_fisher()). Optionally also pass
#' `min_methods` to additionally require at least that many individual
#' methods to have hit their own thresholds (n_methods_hit, from
#' add_vote_counts()/combine_pvalues()) - combines the p-value-combination
#' and vote-counting approaches instead of choosing one. Note n_methods_hit
#' reflects whatever p_cutoff/logfc_cutoff were passed to combine_pvalues()
#' (default 0.05/0.25) - pass matching values there if you change them here.
call_deps <- function(combined, p_col = "FinalP_Fisher", logfc_col = "max_logFC",
                        p_cutoff = 0.05, logfc_cutoff = 0.25, min_methods = NULL) {
  out <- combined %>% filter(.data[[p_col]] <= p_cutoff & abs(.data[[logfc_col]]) > logfc_cutoff)
  if (!is.null(min_methods)) {
    out <- out %>% filter(n_methods_hit >= min_methods)
  }
  out
}

#' Vote-counting alternative to call_deps(): call a protein a DEP if it
#' clears its own p_cutoff/logfc_cutoff (already applied by
#' add_vote_counts()/combine_pvalues() into n_methods_hit) in at least
#' `min_methods` of the methods that tested it. Doesn't assume the methods
#' are independent tests, unlike call_deps()'s p-value-combination approach.
call_deps_vote <- function(combined, min_methods = 3) {
  combined %>% filter(n_methods_hit >= min_methods)
}
