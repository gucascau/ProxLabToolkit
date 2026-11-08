library(imputeLCMD)
library(impute)

#' TRUE for proteins that are entirely missing in at least one whole
#' condition group (structural/MNAR, likely below the limit of detection),
#' as opposed to scattered/random missingness (MAR).
classify_mnar <- function(mat, cond) {
  apply(mat, 1, function(x) any(tapply(x, cond, function(y) all(is.na(y)))))
}

#' Mixed MAR/MNAR imputation: knn for proteins with scattered missingness,
#' MinProb (imputeLCMD) for proteins missing an entire condition group.
#'
#' Deliberately NOT DEP::impute(fun = "mixed", ...): that convenience wrapper
#' (via MsCoreUtils::impute_mixed -> impute_MinProb) transposes the matrix
#' before calling imputeLCMD::impute.MinProb() on only the MNAR subset. With a
#' modest MNAR class size that scrambles the "samples vs features" orientation
#' MinProb expects, so its quantile/sd reference ends up computed across
#' proteins-within-a-sample instead of replicates-within-a-protein - verified
#' against real data to produce unstable, sometimes implausibly high
#' "low-abundance" draws. Calling imputeLCMD::impute.MinProb() directly on the
#' full matrix (correct orientation, stable quantile/sd from the whole dynamic
#' range) avoids that.
#'
#' `rule` (default "partial_mnar"):
#'   "mixed"        knn for scattered gaps, MinProb for proteins missing an
#'                  entire group (the original rule).
#'   "partial_mnar" as "mixed", plus: when a protein is complete in one group
#'                  and missing in at least half (but not all) of the other
#'                  group's runs, that group's gaps are treated as below the
#'                  detection limit (MinProb) rather than missing at random.
#'                  knn otherwise fills such gaps from neighbours detected in
#'                  those runs, which shrinks the fold change of proteins that
#'                  are consistently present with the bait but only sporadically
#'                  in controls. On semi-synthetic spike-ins (40 null splits +
#'                  40 spike-ins) this raised sensitivity at matched error by
#'                  ~10 % with no change in FDR or null-split calls
#'                  (Scripts/Benchmark/impute_variants.R, rule "ocmaj"); on all
#'                  260 null + 260 spike-in splits +6 % sensitivity at matched
#'                  error and FDR 21.8 -> 19.7 % at the default call. Default
#'                  since 2026-11-08; use rule = "mixed" to reproduce earlier runs.
impute_mixed <- function(mat, cond, q = 0.01, tune_sigma = 1, seed = 1000000,
                         rule = c("partial_mnar", "mixed")) {
  rule <- match.arg(rule)
  mnar <- classify_mnar(mat, cond)

  # colmax = 1: impute.knn() otherwise stops when any sample is > 80 % missing,
  # which happens for sparse controls after filtering on the bait (values for
  # samples under the limit are unchanged)
  mar_knn <- impute::impute.knn(mat, colmax = 1)$data

  set.seed(seed)
  mnar_minprob <- imputeLCMD::impute.MinProb(mat, q = q, tune.sigma = tune_sigma)

  out <- mat
  na_idx <- is.na(out)
  out[na_idx & !mnar] <- mar_knn[na_idx & !mnar]
  out[na_idx & mnar] <- mnar_minprob[na_idx & mnar]

  if (rule == "partial_mnar") {
    grp <- unique(cond)
    nmiss <- sapply(grp, function(g) rowSums(na_idx[, cond == g, drop = FALSE]))
    if (is.null(dim(nmiss))) nmiss <- matrix(nmiss, nrow = 1)
    for (j in seq_along(grp)) {
      cols <- which(cond == grp[j])
      other_complete <- rowSums(nmiss[, -j, drop = FALSE]) == 0
      flag <- !mnar & other_complete & nmiss[, j] > 0 & nmiss[, j] >= length(cols) / 2
      sub <- na_idx
      sub[!flag, ] <- FALSE
      sub[, -cols] <- FALSE
      out[sub] <- mnar_minprob[sub]
    }
  }
  out
}
