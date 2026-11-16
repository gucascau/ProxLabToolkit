library(dplyr)

#' Convert a vector of raw p-values into an empirical-Bayes posterior
#' probability of "true" (1 - local FDR), putting the intensity branch's
#' corrected-Fisher p-value onto the same [0,1] probability scale as the count
#' branch's P_count, so the two can be combined below:
#'   pi0 = #{p > 0.5} / (0.5 n), bounded to [1/n, 1]  (Storey, single lambda = 0.5)
#'   lfdr = qvalue::lfdr(p, pi0)
#' A fixed lambda is used instead of qvalue()'s default smoother over lambda in
#' [0.05, 0.95]: the corrected-Fisher p-values match the null only in mean and
#' variance, so they are sparse near 1, and extrapolating pi0 towards lambda = 1
#' then badly underestimates it (e.g. pi0 = 0.14 for a spike-in comparison that
#' is ~98% null); the smoother also fails outright when no p-value reaches 0.95.
#' Falls back to 1 - BH-adjusted p only if qvalue is unavailable or lfdr()
#' fails. Non-finite p-values get NA. The estimator used and pi0 are stored in
#' attr(, "pi0_method") and attr(, "pi0").
to_posterior <- function(pvalues) {
  post <- rep(NA_real_, length(pvalues))
  ok <- is.finite(pvalues)
  p <- pvalues[ok]
  if (length(p) == 0) return(post)
  pi0 <- min(1, max(mean(p > 0.5) / 0.5, 1 / length(p)))
  lf <- NULL
  if (requireNamespace("qvalue", quietly = TRUE)) {
    lf <- tryCatch(qvalue::lfdr(p, pi0 = pi0), error = function(e) NULL)
    if (is.null(lf)) warning("qvalue::lfdr() failed; falling back to 1 - BH-adjusted p-value.")
  } else {
    warning("qvalue package not installed; falling back to 1 - BH-adjusted p-value.")
  }
  how <- "storey_lambda0.5"
  if (is.null(lf)) { lf <- p.adjust(p, method = "BH"); how <- "1-BH"; pi0 <- NA_real_ }
  post[ok] <- 1 - lf
  attr(post, "pi0_method") <- how
  attr(post, "pi0") <- pi0
  post
}

#' Fuse the count branch (run_count_branch()) and intensity branch
#' (combine_pvalues()'s output) into one table: a single Fused_Score for
#' ranking (naive-Bayes / product-of-experts combination of the two
#' branches' posterior probabilities - odds multiply, assuming the two
#' measurement modalities are reasonably independent given true-interactor
#' status, the same independence assumption the intensity consensus already
#' relies on across its own 7 sub-methods) plus a categorical Evidence tag
#' (Both / Count_only / Intensity_only / Neither) derived from each branch's
#' own existing threshold, for interpretability.
#' The intensity posterior is computed from the raw (unadjusted) corrected
#' Fisher p-values (intensity_raw_p_col): the local FDR of qvalue() is defined
#' on raw p-values, and BH-adjusted values are not uniform under the null.
#' The intensity threshold for Evidence still uses the BH-adjusted column
#' (intensity_p_col), as call_deps() does.
#' Defaults intensity_p_col to the correlation-corrected Fisher p-value
#' (FinaladjP_FisherCorrected, see add_corrected_fisher() in consensus.R),
#' matching call_deps()'s recommended default - the plain FinalP_Fisher is
#' anti-conservative (see consensus.R), which previously let genes like
#' KMT2C (significant only by spectral count, not by any individual
#' intensity method) get mislabeled "Both" off a borderline plain-Fisher p.
fuse_branches <- function(count_results, intensity_combined,
                            count_cutoff = 0.75,
                            intensity_p_col = "FinaladjP_FisherCorrected",
                            intensity_raw_p_col = "FinalP_FisherCorrected",
                            intensity_logfc_col = "max_logFC",
                            intensity_p_cutoff = 0.05,
                            intensity_logfc_cutoff = 0.25,
                            intensity_min_methods = 3,
                            direction = c("up", "both")) {
  direction <- match.arg(direction)
  merged <- inner_join(count_results, intensity_combined, by = c("Gene" = "ID"))

  merged$count_posterior <- merged$P_count
  post <- to_posterior(merged[[intensity_raw_p_col]])
  merged$intensity_posterior <- as.numeric(post)
  merged$intensity_pi0_method <- attr(post, "pi0_method")
  merged$intensity_pi0 <- attr(post, "pi0")

  eps <- 1e-6
  clip <- function(p) pmin(pmax(p, eps), 1 - eps)
  odds <- function(p) clip(p) / (1 - clip(p))

  fused_odds <- odds(merged$count_posterior) * odds(merged$intensity_posterior)
  merged$Fused_Score <- fused_odds / (1 + fused_odds)

  count_positive <- merged$count_posterior >= count_cutoff
  # Same three criteria as call_deps() in run_pl_dea(), so a protein tagged
  # intensity-positive here is exactly an intensity-consensus DEP. Without the
  # vote criterion, a protein passing p/logFC but not the vote was tagged
  # "Both" yet fell out of the final list (neither a DEP nor Count_only).
  votes_ok <- if (is.null(intensity_min_methods)) TRUE else
    !is.na(merged$n_methods_hit) & merged$n_methods_hit >= intensity_min_methods
  intensity_positive <- !is.na(merged[[intensity_p_col]]) &
    merged[[intensity_p_col]] <= intensity_p_cutoff &
    passes_logfc(merged[[intensity_logfc_col]], intensity_logfc_cutoff, direction) & votes_ok

  merged$Evidence <- case_when(
    count_positive & intensity_positive ~ "Both",
    count_positive & !intensity_positive ~ "Count_only",
    !count_positive & intensity_positive ~ "Intensity_only",
    TRUE ~ "Neither"
  )

  merged %>% arrange(desc(Fused_Score))
}
