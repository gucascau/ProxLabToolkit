#' Extract a raw spectral/PSM count matrix (preys x samples) for one
#' comparison spec. Unlike extract_intensity_matrix(), zero is a real
#' observation here (a prey genuinely wasn't seen), not a missing value, so
#' there is no 0->NA step or NA-count filter.
extract_count_matrix <- function(pro_file_update, spec, count_type = "Total.Spectral.Count") {
  sample_cols <- c(spec$group1_cols, spec$group2_cols)
  count_cols <- paste0(sample_cols, ".", count_type)
  missing <- setdiff(count_cols, colnames(pro_file_update))
  if (length(missing) > 0) {
    stop("Missing count columns: ", paste(missing, collapse = ", "))
  }
  mat <- as.matrix(pro_file_update[, count_cols])
  mat[is.na(mat)] <- 0
  rownames(mat) <- pro_file_update$Gene
  colnames(mat) <- sample_cols
  mat
}

#' Original, SAINT-*inspired* two-component negative-binomial mixture, fit by
#' EM jointly across all preys (the core SAINT idea: borrow strength across
#' the whole dataset to learn what "background" and "true interactor" counts
#' look like, rather than testing each prey in isolation). This is NOT a port
#' of the published SAINTexpress algorithm - it's a from-scratch, simplified
#' reimplementation of the same underlying philosophy, documented here so it
#' is never mistaken for the validated original tool's output:
#'
#' For each prey i, let expected_bait_sum_i be what the summed bait-replicate
#' count would be if this prey behaved exactly like its own matched-control
#' replicates (i.e. FC = 1, no enrichment), after scaling for each
#' replicate's total library size. Bait counts are then modeled as a mixture
#' of two negative-binomial components:
#'   - "background": mean = expected_bait_sum_i            (FC fixed at 1)
#'   - "true interactor": mean = expected_bait_sum_i * FC1  (FC1 > 1, a single
#'     shared enrichment level, EM-estimated across all preys)
#' with dispersions for each component also EM-estimated (weighted
#' method-of-moments in the M-step, rather than refitting a GLM every
#' iteration, for speed/stability).
#'
#' Returns, per prey: P_count (posterior probability of the "true
#' interactor" component - the SAINT-style confidence score) and a point
#' log2 fold-change for ranking ties/reporting.
fit_saint_style_mixture <- function(bait_counts, control_counts,
                                      max_iter = 200, tol = 1e-4,
                                      pseudocount = 0.5, min_dispersion = 1e-4,
                                      pi_init = 0.1, fc1_init = 4) {
  control_sum <- rowSums(control_counts)
  bait_sum <- rowSums(bait_counts)

  sf_bait <- colSums(bait_counts)
  sf_control <- colSums(control_counts)
  mean_libsize <- mean(c(sf_bait, sf_control))
  sf_bait <- sf_bait / mean_libsize
  sf_control <- sf_control / mean_libsize

  expected_bait_sum <- (control_sum + pseudocount) / sum(sf_control) * sum(sf_bait)

  nb_dens <- function(x, mu, phi) {
    dnbinom(x, size = 1 / phi, mu = pmax(mu, 1e-6))
  }
  weighted_dispersion <- function(x, mu, w) {
    wsum <- sum(w)
    wmean <- sum(w * mu) / wsum
    wvar <- sum(w * (x - mu)^2) / wsum
    max(min_dispersion, (wvar - wmean) / wmean^2)
  }

  pi_true <- pi_init
  fc1 <- fc1_init
  phi0 <- max(min_dispersion, (var(bait_sum) - mean(expected_bait_sum)) / mean(expected_bait_sum)^2)
  phi1 <- phi0
  resp <- rep(pi_init, length(bait_sum))

  for (iter in seq_len(max_iter)) {
    d0 <- nb_dens(bait_sum, expected_bait_sum, phi0)
    d1 <- nb_dens(bait_sum, expected_bait_sum * fc1, phi1)
    numerator <- pi_true * d1
    denom <- numerator + (1 - pi_true) * d0
    denom[denom <= 0 | is.na(denom)] <- .Machine$double.eps
    resp <- numerator / denom
    resp[is.na(resp)] <- 0

    pi_new <- min(max(mean(resp), 1e-4), 0.9)

    if (sum(resp) > 1e-6) {
      fc1_new <- max(sum(resp * bait_sum) / sum(resp * expected_bait_sum), 1.01)
    } else {
      fc1_new <- fc1
    }

    phi0_new <- weighted_dispersion(bait_sum, expected_bait_sum, 1 - resp)
    phi1_new <- weighted_dispersion(bait_sum, expected_bait_sum * fc1_new, resp)

    delta <- abs(pi_new - pi_true) + abs(fc1_new - fc1)
    pi_true <- pi_new
    fc1 <- fc1_new
    phi0 <- phi0_new
    phi1 <- phi1_new
    if (delta < tol) break
  }

  log2fc <- log2((bait_sum + pseudocount) / (expected_bait_sum + pseudocount))
  list(P_count = resp, count_logFC = log2fc, fc1 = fc1, pi_true = pi_true,
       phi0 = phi0, phi1 = phi1, iterations = iter)
}

#' Orchestrates extract_count_matrix() + fit_saint_style_mixture() for one
#' comparison spec, returning a standardized (Gene, P_count, count_logFC) table.
run_count_branch <- function(pro_file_update, spec, count_type = "Total.Spectral.Count", ...) {
  mat <- extract_count_matrix(pro_file_update, spec, count_type = count_type)
  bait_counts <- mat[, spec$group1_cols, drop = FALSE]
  control_counts <- mat[, spec$group2_cols, drop = FALSE]
  fit <- fit_saint_style_mixture(bait_counts, control_counts, ...)
  data.frame(Gene = rownames(mat), P_count = fit$P_count, count_logFC = fit$count_logFC,
             row.names = NULL)
}
