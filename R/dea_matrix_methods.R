library(limma)
library(DEqMS)
library(dplyr)
library(tibble)

# Shared design/contrast fit reused by run_limma() and run_deqms() so the two
# stay consistent (DEqMS is just Limma's fit + a PSM-count variance correction).
.fit_limma <- function(dat_log_exp, spec) {
  cond <- factor(c(rep(spec$group1_label, length(spec$group1_cols)),
                    rep(spec$group2_label, length(spec$group2_cols))))
  design <- model.matrix(~0 + cond)
  colnames(design) <- gsub("^cond", "", colnames(design))
  contrast <- makeContrasts(contrasts = paste0(spec$group1_label, "-", spec$group2_label), levels = design)
  fit1 <- lmFit(dat_log_exp, design)
  fit2 <- contrasts.fit(fit1, contrasts = contrast)
  eBayes(fit2)
}

run_limma <- function(dat_log_exp, spec) {
  fit3 <- .fit_limma(dat_log_exp, spec)
  topTable(fit3, coef = 1, adjust = "fdr", number = Inf) %>%
    rownames_to_column("Gene")
}

#' DEqMS needs a minimum-PSM-count-per-protein covariate (see load_psm_counts()
#' / min_psm_for_comparison()) as a named vector keyed by Gene, matching
#' rownames(dat_log_exp). Missing genes get a PSM count of 0 (+1 pseudocount,
#' as the original pipeline does, since spectraCounteBayes needs count >= 1).
run_deqms <- function(dat_log_exp, spec, min_psm_by_gene) {
  fit3 <- .fit_limma(dat_log_exp, spec)
  count <- min_psm_by_gene[rownames(fit3$coefficients)]
  count[is.na(count)] <- 0
  fit3$count <- count + 1
  fit4 <- spectraCounteBayes(fit3)
  # outputResult() keeps limma's ordinary statistics in t/P.Value/adj.P.Val and
  # DEqMS's PSM-count-moderated ones in sca.*. Expose DEqMS's as the standard
  # columns so the consensus uses DEqMS itself rather than a copy of limma.
  outputResult(fit4, coef_col = 1) %>%
    rownames_to_column("Gene") %>%
    rename(limma_t = t, limma_P.Value = P.Value, limma_adj.P.Val = adj.P.Val) %>%
    mutate(t = sca.t, P.Value = sca.P.Value, adj.P.Val = sca.adj.pval)
}

run_rots <- function(dat_log_exp, spec, B = 10000, K = 500, seed = 1234) {
  n1 <- length(spec$group1_cols)
  n2 <- length(spec$group2_cols)
  groups <- c(rep(0, n1), rep(1, n2))
  results <- ROTS::ROTS(data = dat_log_exp, groups = groups, B = B, K = K, seed = seed)
  data.frame(Gene = rownames(results$data), logFC = results$logfc,
             P.Value = results$pvalue, adj.pval = results$FDR)
}
