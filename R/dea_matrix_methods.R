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
  outputResult(fit4, coef_col = 1) %>% rownames_to_column("Gene")
}

run_wilcoxon <- function(dat_log_exp, spec) {
  n1 <- length(spec$group1_cols)
  n2 <- length(spec$group2_cols)
  g1 <- seq_len(n1)
  g2 <- (n1 + 1):(n1 + n2)
  pval <- apply(dat_log_exp, 1, function(x) wilcox.test(x[g1], x[g2])$p.value)
  logFC <- rowMeans(dat_log_exp[, g1, drop = FALSE]) - rowMeans(dat_log_exp[, g2, drop = FALSE])
  data.frame(Gene = rownames(dat_log_exp), logFC = logFC, P.Value = pval,
             adj.pval = p.adjust(pval, method = "fdr"))
}

run_rots <- function(dat_log_exp, spec, B = 10000, K = 500, seed = 1234) {
  n1 <- length(spec$group1_cols)
  n2 <- length(spec$group2_cols)
  groups <- c(rep(0, n1), rep(1, n2))
  results <- ROTS::ROTS(data = dat_log_exp, groups = groups, B = B, K = K, seed = seed)
  data.frame(Gene = rownames(results$data), logFC = results$logfc,
             P.Value = results$pvalue, adj.pval = results$FDR)
}
