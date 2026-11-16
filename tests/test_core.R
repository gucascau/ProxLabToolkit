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
stopifnot(attr(out, "fisher_correction")$method == "robust")

## 2b Brown's method: calibrated with stronger, unequal correlation and methods missing at random,
##    and never more liberal than plain Fisher in the rejection region
z <- rnorm(n)
rho <- c(0.95, 0.9, 0.6, 0.3)
pm4 <- sapply(rho, function(r) pnorm(r * z + sqrt(1 - r^2) * rnorm(n), lower.tail = FALSE))
pm4[sample(length(pm4), 0.1 * length(pm4))] <- NA
df4 <- setNames(as.data.frame(pm4), paste0("M", 1:4, "_P.Value"))
b4 <- add_corrected_fisher(df4, names(df4), method = "brown")
r4 <- add_corrected_fisher(df4, names(df4))                      # robust (default)
k4 <- add_corrected_fisher(df4, names(df4), method = "kost")
plain4 <- apply(pm4, 1, function(p) { p <- p[!is.na(p)]; pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE) })
stopifnot(mean(b4$FinalP_FisherCorrected <= 0.05) < 0.07, mean(r4$FinalP_FisherCorrected <= 0.05) < 0.07)
tail <- !is.na(plain4) & plain4 <= 0.05
stopifnot(all(b4$FinalP_FisherCorrected[tail] >= plain4[tail] - 1e-12))   # in the rejection region, never below plain Fisher
stopifnot(all(r4$FinalP_FisherCorrected[tail] >= plain4[tail] - 1e-12))
stopifnot(mean(b4$FinalP_FisherCorrected <= 0.05) <= mean(k4$FinalP_FisherCorrected <= 0.05))

## 2c true interactors must not inflate the robust dependence estimate (they do inflate Brown's)
zs <- rnorm(n)
pn <- sapply(c(0.6, 0.6, 0.6), function(r) pnorm(r * zs + sqrt(1 - r^2) * rnorm(n), lower.tail = FALSE))
sig <- 1:300                                                       # 7.5% true interactors, strong in every method
pn[sig, ] <- pnorm(rnorm(300 * 3, 4), lower.tail = FALSE)
dn <- setNames(as.data.frame(pn), paste0("M", 1:3, "_P.Value"))
cr <- attr(add_corrected_fisher(dn, names(dn)), "fisher_correction")$cov
cb <- attr(add_corrected_fisher(dn, names(dn), method = "brown"), "fisher_correction")$cov
r0 <- 0.36                                                         # true correlation between methods (0.6 x 0.6)
cnull <- r0 * (3.263 + 0.710 * r0 + 0.027 * r0^2)                 # covariance implied by the true null correlation
stopifnot(abs(mean(cr[upper.tri(cr)]) - cnull) < 0.35, mean(cb[upper.tri(cb)]) > cnull + 0.5)
stopifnot(abs(robust_cor(1:100 + rnorm(100), 1:100 + rnorm(100)) - 1) < 0.05, robust_cor(1:5, 1:5) == 0)

## 3  call_deps(): thresholds and method-count filter
comb <- data.frame(Gene = c("A", "B", "C", "D"), FinaladjP_FisherCorrected = c(0.01, 0.01, 0.2, 0.01),
                   max_logFC = c(1, 0.1, 2, -1), n_methods_hit = c(4, 4, 4, 2))
d <- call_deps(comb, p_col = "FinaladjP_FisherCorrected", logfc_cutoff = 0.25, min_methods = 3)
stopifnot(identical(d$Gene, "A"))
stopifnot(identical(call_deps(comb, p_col = "FinaladjP_FisherCorrected")$Gene, "A"))          # enriched only (default)
stopifnot(identical(call_deps(comb, p_col = "FinaladjP_FisherCorrected", direction = "both")$Gene, c("A", "D")))

## 3b add_vote_counts(): depletion is not a hit by default
v <- data.frame(M1_logFC = c(1, -1), M1_P.Value = c(0.01, 0.01))
stopifnot(identical(add_vote_counts(v, "M1_logFC", "M1_P.Value")$n_methods_hit, c(1L, 0L)))
stopifnot(identical(add_vote_counts(v, "M1_logFC", "M1_P.Value", direction = "both")$n_methods_hit, c(1L, 1L)))

## 3c fuse_branches(): intensity posterior comes from raw p-values, not BH-adjusted ones
ic <- data.frame(ID = paste0("G", 1:500), FinalP_FisherCorrected = c(runif(450), rbeta(50, 0.05, 1)),
                 max_logFC = 1, n_methods_hit = 4)
ic$FinaladjP_FisherCorrected <- p.adjust(ic$FinalP_FisherCorrected, "BH")
cr <- data.frame(Gene = ic$ID, P_count = 0.5)
fz <- fuse_branches(cr, ic)
fz <- fz[match(ic$ID, fz$Gene), ]
stopifnot(isTRUE(all.equal(fz$intensity_posterior, as.numeric(to_posterior(ic$FinalP_FisherCorrected)))))

## 3d to_posterior(): strongly enriched comparison (no p above 0.86, where qvalue's default smoother fails)
##     works, keeps NAs, and stays a probability; a uniform null gives pi0 near 1 and low posteriors
pe <- c(runif(100, 0, 0.86), rbeta(900, 0.1, 5), NA)
po <- to_posterior(pe)
stopifnot(attr(po, "pi0_method") == "storey_lambda0.5", is.na(po[1001]),
          all(po[-1001] >= 0 & po[-1001] <= 1), cor(po[-1001], pe[-1001], method = "spearman") < -0.9)

## 4  comparison_spec(): structure
sp <- comparison_spec("x_vs_y", "x", paste0("x_", 1:3), "y", paste0("y_", 1:3))
stopifnot(sp$name == "x_vs_y", length(sp$group1_cols) == 3, sp$max_na == 3)


## 5  to_posterior(): p-values that are uniform except sparse near 1 (as corrected Fisher gives) keep pi0 near 1
pu <- c(runif(2400, 0, 0.85), runif(50, 0.85, 1))
stopifnot(attr(to_posterior(pu), "pi0") > 0.8, median(to_posterior(pu)) < 0.3)

cat("ProxLabToolkit core tests: all passed\n")
