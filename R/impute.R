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
impute_mixed <- function(mat, cond, q = 0.01, tune_sigma = 1, seed = 1000000) {
  mnar <- classify_mnar(mat, cond)

  mar_knn <- impute::impute.knn(mat)$data

  set.seed(seed)
  mnar_minprob <- imputeLCMD::impute.MinProb(mat, q = q, tune.sigma = tune_sigma)

  out <- mat
  na_idx <- is.na(out)
  out[na_idx & !mnar] <- mar_knn[na_idx & !mnar]
  out[na_idx & mnar] <- mnar_minprob[na_idx & mnar]
  out
}
