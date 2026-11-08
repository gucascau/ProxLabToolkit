# Core unit tests (base R; run by R CMD check, or: Rscript tests/test_core.R from the package root).
# Uses the installed package when available, otherwise sources R/.
if (requireNamespace("ProxLabToolkit", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ProxLabToolkit))
} else {
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
}
set.seed(1)

## 1  impute_mixed(): partial_mnar sends gaps below the detection limit when the other group is complete
cond <- rep(c("bait", "ctrl"), each = 3)
mat <- matrix(rnorm(60 * 6, 25, 1), 60, 6, dimnames = list(paste0("P", 1:60), paste0(cond, 1:3)))
mat["P1", 4:5] <- NA          # complete in bait, 2/3 missing in control -> partial_mnar target
mat["P2", 4:6] <- NA          # whole control missing -> MNAR in both rules
mat["P3", 1] <- NA            # one scattered gap -> kNN in both rules
pm <- impute_mixed(mat, cond, rule = "partial_mnar")
mx <- impute_mixed(mat, cond, rule = "mixed")
stopifnot(!anyNA(pm), !anyNA(mx))
low <- quantile(mat, 0.10, na.rm = TRUE)
med <- median(mat, na.rm = TRUE)
# MinProb draws (random, around the 1st percentile) sit in the low tail; kNN values sit near the middle
stopifnot(mean(pm["P1", 4:5]) < low, all(pm["P1", 4:5] < med))
stopifnot(mean(mx["P1", 4:5]) > low)                     # the old rule used kNN here
stopifnot(mean(pm["P2", 4:6]) < low, mean(mx["P2", 4:6]) < low)
stopifnot(isTRUE(all.equal(pm["P3", ], mx["P3", ])))    # scattered gaps unaffected by the rule
stopifnot(classify_mnar(mat, cond)[["P2"]], !classify_mnar(mat, cond)[["P1"]])

## 2  add_corrected_fisher(): calibrated under the null even when methods are correlated
n <- 4000
z <- rnorm(n)
pm3 <- sapply(1:3, function(i) pnorm(0.9 * z + sqrt(1 - 0.81) * rnorm(n), lower.tail = FALSE))
df <- setNames(as.data.frame(pm3), paste0("M", 1:3, "_P.Value"))
out <- add_corrected_fisher(df, names(df))
plain <- apply(pm3, 1, function(p) pchisq(-2 * sum(log(p)), 6, lower.tail = FALSE))
fpr_corr <- mean(out$FinalP_FisherCorrected <= 0.05)
fpr_plain <- mean(plain <= 0.05)
stopifnot(fpr_plain > 0.08)                              # plain Fisher is anti-conservative here
stopifnot(fpr_corr < 0.075, fpr_corr > 0.02)             # corrected Fisher stays near the nominal 5%
stopifnot(all(abs(out$FinaladjP_FisherCorrected - p.adjust(out$FinalP_FisherCorrected, "BH")) < 1e-12))

## 3  call_deps(): thresholds and method-count filter
comb <- data.frame(Gene = c("A", "B", "C", "D"), FinaladjP_FisherCorrected = c(0.01, 0.01, 0.2, 0.01),
                   max_logFC = c(1, 0.1, 2, -1), n_methods_hit = c(4, 4, 4, 2))
d <- call_deps(comb, p_col = "FinaladjP_FisherCorrected", logfc_cutoff = 0.25, min_methods = 3)
stopifnot(identical(d$Gene, "A"))
stopifnot(identical(call_deps(comb, p_col = "FinaladjP_FisherCorrected")$Gene, c("A", "D")))

## 4  comparison_spec(): structure
sp <- comparison_spec("x_vs_y", "x", paste0("x_", 1:3), "y", paste0("y_", 1:3))
stopifnot(sp$name == "x_vs_y", length(sp$group1_cols) == 3, sp$max_na == 3)

cat("ProxLabToolkit core tests: all passed\n")
