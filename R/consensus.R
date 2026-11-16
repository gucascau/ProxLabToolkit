library(dplyr)
library(purrr)
library(metap)

#' Join every method's result table into one wide table, prefixing each
#' method's columns with its name (e.g. "Limma_logFC") and keyed by a single
#' `ID` column (Gene symbol - every run_* wrapper in this toolkit already
#' standardizes on a "Gene" column).
#'
#' The first method (Limma, run on the full filtered + imputed intensity
#' matrix) defines the protein set; the others are left-joined, so a protein
#' a method could not test gets NA there and is simply not counted for that
#' method downstream (every later step skips NA p-values/logFCs). The original
#' PvalueIntegration.Rmd used inner_join across all methods, which silently
#' dropped every protein MSstats cannot model - i.e. proteins detected only in
#' the bait group, often the bait itself and its strongest proximal partners.
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
  reduce(renamed, function(x, y) left_join(x, y, by = "ID"))
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
#' they mostly aren't - Limma/DEqMS/ROTS/DEP/proDA are computed on
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

#' Correlation-corrected Fisher combination. The methods combined here share
#' the same underlying data, so their p-values are not independent and plain
#' Fisher (`FinalP_Fisher` above) is anti-conservative. Fisher's statistic
#' T = sum(-2 log p_m) is instead referred to a scaled chi-square c * chi2_f
#' whose mean (2k) and variance match those of T under dependence:
#'   Var(T) = 4k + 2 * sum_{m<l} Cov(X_m, X_l),  X_m = -2 log p_m,
#'   c = Var(T) / (4k),  f = 8k^2 / Var(T),  p = P(chi2_f >= T / c).
#'
#' method = "robust" (default): each method's p-values are turned into normal scores u = qnorm(p), which are
#' standard normal under the null, and the correlation rho of each pair of methods is estimated with the
#' Gnanadesikan-Kettenring estimator on median/MAD-standardized scores (robust_cor()). Because it uses the bulk of
#' proteins, strongly enriched true interactors (extreme scores in every method) do not inflate it. The covariance of
#' -2 log p follows from Kost & McDermott's (2002) cubic, cov = 3.263 rho + 0.710 rho^2 + 0.027 rho^3, which is
#' defined for exactly this input (correlation of the underlying normal statistics). rho is bounded to [0, max_rho].
#'
#' method = "brown" (the empirical Brown's method, Poole et al. 2016):
#' Cov(X_m, X_l) is estimated directly as the sample covariance of the -2 log p
#' statistics of each pair of methods, across all proteins of this comparison
#' (pairwise-complete), with Var(X_m) fixed at its null value of 4. A row tested
#' by only some methods uses the covariances of those methods. Each covariance
#' is bounded to [0, 4] (between independence and identical tests), so Var(T)
#' is never below its value under independence and, in the upper tail used for
#' calling, the corrected p-value is never smaller than plain Fisher's. Including true
#' interactors in the estimate raises the covariances; in comparisons with many
#' true interactors they reach the bound (tests treated as identical), so this
#' option is strongly conservative exactly where there is signal.
#'
#' method = "kost": the earlier version (Kost & McDermott 2002 cubic, applied
#' to the mean Pearson correlation of the -2 log p statistics). The cubic
#' expects the correlation of the underlying normal statistics, so applied to
#' -2 log p it under-corrects; kept only to reproduce earlier results.
add_corrected_fisher <- function(combined, pvalue_cols, method = c("robust", "brown", "kost"),
                                 max_rho = 0.95) {
  method <- match.arg(method)
  pmat <- as.matrix(combined[, pvalue_cols, drop = FALSE])
  stat_mat <- -2 * log(pmat)
  k_all <- ncol(stat_mat)

  if (method == "robust") {
    # normal scores of the p-values (standard normal under the null), robust pairwise correlation from the bulk of
    # proteins, then Kost & McDermott's cubic, which maps the correlation of the underlying normal statistics to the
    # covariance of -2 log p
    u <- qnorm(pmin(pmax(pmat, 1e-300), 1 - 1e-12))
    rho <- matrix(0, k_all, k_all)
    for (a in seq_len(k_all - 1)) for (b in (a + 1):k_all) {
      r <- robust_cor(u[, a], u[, b])
      rho[a, b] <- rho[b, a] <- r
    }
    rho <- pmin(pmax(rho, 0), max_rho)
    cov_mat <- rho * (3.263 + 0.710 * rho + 0.027 * rho^2)
    diag(cov_mat) <- 4
  } else if (method == "brown") {
    fin <- stat_mat
    fin[!is.finite(fin)] <- NA
    cov_mat <- suppressWarnings(cov(fin, use = "pairwise.complete.obs"))
    cov_mat[!is.finite(cov_mat)] <- 0
    cov_mat <- pmin(pmax(cov_mat, 0), 4)
    diag(cov_mat) <- 4
  } else {
    cor_mat <- suppressWarnings(cor(stat_mat, use = "pairwise.complete.obs"))
    r <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
    if (!is.finite(r)) r <- 0
    r <- max(0, min(r, 0.9))
    cov_mat <- matrix(r * (3.263 + 0.710 * r + 0.027 * r^2), k_all, k_all)
    diag(cov_mat) <- 4
  }

  ok <- !is.na(pmat)
  combined$FinalP_FisherCorrected <- vapply(seq_len(nrow(pmat)), function(i) {
    idx <- which(ok[i, ])
    k <- length(idx)
    if (k == 0) return(NA_real_)
    if (k == 1) return(pmat[i, idx])
    stat <- sum(stat_mat[i, idx])
    var_t <- sum(cov_mat[idx, idx])          # = 4k + 2 * sum of pairwise covariances
    e_t <- 2 * k
    f <- 2 * e_t^2 / var_t
    c_scale <- var_t / (2 * e_t)
    pchisq(stat / c_scale, df = f, lower.tail = FALSE)
  }, numeric(1))
  combined$FinaladjP_FisherCorrected <- p.adjust(combined$FinalP_FisherCorrected, method = "BH")
  attr(combined, "fisher_correction") <- list(method = method, cov = cov_mat)
  combined
}

#' Vote-counting alternative to p-value combination: per protein, count how
#' many methods it was testable in (`n_methods_tested`) and how many of
#' those calls clear that method's own p_cutoff/logfc_cutoff
#' (`n_methods_hit`); with direction = "up" (default) only enrichment over the
#' control counts as a hit. Unlike Fisher/Tippett/Stouffer, this makes no
#' independence assumption about the methods, at the cost of losing a single
#' interpretable combined p-value/effect size. See `call_deps_vote()`.
add_vote_counts <- function(combined, logfc_cols, pvalue_cols,
                             p_cutoff = 0.05, logfc_cutoff = 0.25,
                             direction = c("up", "both")) {
  direction <- match.arg(direction)
  stopifnot(length(logfc_cols) == length(pvalue_cols))
  hits <- Map(function(lfc_col, p_col) {
    !is.na(combined[[lfc_col]]) & !is.na(combined[[p_col]]) &
      passes_logfc(combined[[lfc_col]], logfc_cutoff, direction) & combined[[p_col]] <= p_cutoff
  }, logfc_cols, pvalue_cols)
  tested <- lapply(pvalue_cols, function(p_col) !is.na(combined[[p_col]]))

  combined$n_methods_tested <- as.integer(Reduce(`+`, tested))
  combined$n_methods_hit <- as.integer(Reduce(`+`, hits))
  combined
}

#' Full consensus step: join all methods, pick max_logFC, combine p-values
#' (plain + correlation-corrected Fisher), and tally per-protein vote counts.
#' Every run_* wrapper in this toolkit standardizes on logFC/P.Value columns,
#' so by default logfc_cols/pvalue_cols are just "<method>_logFC" /
#' "<method>_P.Value" for every method in results_list - override to match
#' the original scripts' choice of excluding a method from the combined
#' statistic. Wilcoxon is no longer run at all: with 3 vs 3 replicates its
#' smallest attainable two-sided p-value is 2/choose(6, 3) = 0.1, so it could
#' never reach significance and only diluted the combined statistic.
#' `p_cutoff`/`logfc_cutoff` are only used for the vote-count columns
#' (n_methods_tested/n_methods_hit) - they don't affect the p-value columns.
combine_pvalues <- function(results_list, id_col = "Gene",
                             logfc_cols = paste0(names(results_list), "_logFC"),
                             pvalue_cols = paste0(names(results_list), "_P.Value"),
                             n_perm = 100000, seed = 123,
                             p_cutoff = 0.05, logfc_cutoff = 0.25,
                             direction = "up", correction = "robust") {
  combined <- join_method_results(results_list, id_col = id_col)
  combined <- add_max_logfc(combined, logfc_cols)
  combined <- add_combined_pvalues(combined, pvalue_cols, n_perm = n_perm, seed = seed)
  combined <- add_corrected_fisher(combined, pvalue_cols, method = correction)
  add_vote_counts(combined, logfc_cols, pvalue_cols,
                   p_cutoff = p_cutoff, logfc_cutoff = logfc_cutoff, direction = direction)
}

#' Default call: adjusted p <= 0.05 & logFC > 0.25, the convention of other proteomics tools. A more sensitive
#' setting is p_cutoff = 0.10 (options$dep_p_cutoff / options$intensity_p_cutoff in run_pl_interactome()): on
#' the project's spike-ins it raised sensitivity from 37% to 46% at more null false calls (15% vs 9% of
#' control-vs-control splits with any call). With
#' direction = "up" (default) only proteins enriched over the control are
#' called, as a proximity interactome should be; direction = "both" restores
#' the earlier abs(logFC) behavior (also calls depleted proteins).
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
                        p_cutoff = 0.05, logfc_cutoff = 0.25, min_methods = NULL,
                        direction = c("up", "both")) {
  direction <- match.arg(direction)
  out <- combined %>%
    filter(.data[[p_col]] <= p_cutoff & passes_logfc(.data[[logfc_col]], logfc_cutoff, direction))
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

#' logFC filter shared by the vote counts, call_deps() and fuse_branches():
#' "up" keeps enrichment over the control only, "both" uses abs(logFC).
passes_logfc <- function(lfc, cutoff, direction = "up") {
  if (direction == "up") !is.na(lfc) & lfc > cutoff else !is.na(lfc) & abs(lfc) > cutoff
}

#' Robust correlation of two vectors (pairwise-complete): Gnanadesikan-Kettenring estimator with the MAD as scale,
#' r = (s(x+y)^2 - s(x-y)^2) / (s(x+y)^2 + s(x-y)^2) on median/MAD-standardized x and y. Resistant to a minority of
#' outlying points (breakdown 50%); returns 0 when fewer than 10 complete pairs or a zero scale.
robust_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(0)
  x <- x[ok]
  y <- y[ok]
  sx <- mad(x)
  sy <- mad(y)
  if (sx == 0 || sy == 0) return(0)
  x <- (x - median(x)) / sx
  y <- (y - median(y)) / sy
  a <- mad(x + y)^2
  b <- mad(x - y)^2
  if (a + b == 0) return(0)
  (a - b) / (a + b)
}
