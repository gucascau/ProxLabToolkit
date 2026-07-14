library(dplyr)

#' Convert a vector of p-values into an empirical-Bayes posterior probability
#' of "true" (1 - local FDR), putting the intensity branch's Fisher-combined
#' p-value onto the same [0,1] probability scale as the count branch's
#' P_count, so the two can be combined below. Prefers qvalue::lfdr(); falls
#' back to a cruder 1 - BH-adjusted-p approximation if qvalue isn't
#' installed or its pi0 estimation fails (e.g. when a large fraction of
#' proteins are genuinely non-null, which qvalue's method can struggle with).
to_posterior <- function(pvalues) {
  if (requireNamespace("qvalue", quietly = TRUE)) {
    post <- tryCatch(1 - qvalue::qvalue(pvalues)$lfdr, error = function(e) NULL)
    if (!is.null(post)) return(post)
    warning("qvalue::qvalue() failed; falling back to 1 - BH-adjusted p-value.")
  } else {
    warning("qvalue package not installed; falling back to 1 - BH-adjusted p-value.")
  }
  1 - p.adjust(pvalues, method = "BH")
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
#' Defaults intensity_p_col to the correlation-corrected Fisher p-value
#' (FinaladjP_FisherCorrected, see add_corrected_fisher() in consensus.R),
#' matching call_deps()'s recommended default - the plain FinalP_Fisher is
#' anti-conservative (see consensus.R), which previously let genes like
#' KMT2C (significant only by spectral count, not by any individual
#' intensity method) get mislabeled "Both" off a borderline plain-Fisher p.
fuse_branches <- function(count_results, intensity_combined,
                            count_cutoff = 0.75,
                            intensity_p_col = "FinaladjP_FisherCorrected",
                            intensity_logfc_col = "max_logFC",
                            intensity_p_cutoff = 0.05,
                            intensity_logfc_cutoff = 0.25) {
  merged <- inner_join(count_results, intensity_combined, by = c("Gene" = "ID"))

  merged$count_posterior <- merged$P_count
  merged$intensity_posterior <- to_posterior(merged[[intensity_p_col]])

  eps <- 1e-6
  clip <- function(p) pmin(pmax(p, eps), 1 - eps)
  odds <- function(p) clip(p) / (1 - clip(p))

  fused_odds <- odds(merged$count_posterior) * odds(merged$intensity_posterior)
  merged$Fused_Score <- fused_odds / (1 + fused_odds)

  count_positive <- merged$count_posterior >= count_cutoff
  intensity_positive <- merged[[intensity_p_col]] <= intensity_p_cutoff &
    abs(merged[[intensity_logfc_col]]) > intensity_logfc_cutoff

  merged$Evidence <- case_when(
    count_positive & intensity_positive ~ "Both",
    count_positive & !intensity_positive ~ "Count_only",
    !count_positive & intensity_positive ~ "Intensity_only",
    TRUE ~ "Neither"
  )

  merged %>% arrange(desc(Fused_Score))
}
