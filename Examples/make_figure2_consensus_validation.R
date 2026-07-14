# Figure 2 - the consensus-combination statistical fix, shown on real data:
# (A) the 7 DEA methods' p-values are highly correlated, not independent;
# (B) plain Fisher combination over-calls DEPs relative to the
#     correlation-corrected version and vote-counting;
# (C) a p-value QQ-plot makes the resulting inflation visible directly.
#
# Reads the already-generated negtKTSwt_R467W_vs_Mock comparison outputs -
# run Examples/example_R467W_vs_Mock.R first if these don't exist yet.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(reshape2)
})

out_dir  <- "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Results/NegKTSR467WIndvsMockInd/PLToolkitDemo"
fig_dir  <- "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/scripts/PLToolkit/Figs"
combined <- read.csv(file.path(out_dir, "negtKTSwt_R467W_vs_Mock__FinalCombinedPvalue_IntegratedMethods.csv"))

# validated palette (see dataviz skill references/palette.md)
col_blue   <- "#2a78d6"  # categorical slot 1 - plain Fisher
col_aqua   <- "#1baf7a"  # categorical slot 2 - corrected Fisher
col_yellow <- "#eda100"  # categorical slot 3 - vote-counting
col_muted  <- "#898781"
col_grid   <- "#e1e0d9"
seq_blue   <- c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95")  # steps 100/200/350/450/600

# ---- Panel A: method correlation heatmap ---------------------------------
pcols <- c("Limma_P.Value", "Wilcoxon_P.Value", "ROTS_P.Value", "DEP_P.Value",
           "proDA_P.Value", "DEqMS_P.Value", "MSstats_P.Value")
method_labels <- c("Limma", "Wilcoxon", "ROTS", "DEP", "proDA", "DEqMS", "MSstats")

stat_mat <- -2 * log(as.matrix(combined[, pcols]))
colnames(stat_mat) <- method_labels
cor_mat <- cor(stat_mat, use = "pairwise.complete.obs")

cor_df <- melt(cor_mat, varnames = c("Method1", "Method2"), value.name = "r")
cor_df$Method1 <- factor(cor_df$Method1, levels = method_labels)
cor_df$Method2 <- factor(cor_df$Method2, levels = rev(method_labels))

panelA <- ggplot(cor_df, aes(Method1, Method2, fill = r)) +
  geom_tile(color = "#fcfcfb", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.6, color = "#0b0b0b") +
  scale_fill_gradientn(colors = seq_blue, limits = c(0.5, 1),
                        name = "Pearson r\n(-2log p)") +
  coord_fixed() +
  labs(title = "A. DEA methods are correlated, not independent",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

# ---- Panel B: DEP counts under 3 calling strategies -----------------------
n_plain     <- sum(combined$FinalP_Fisher <= 0.05 & abs(combined$max_logFC) > 0.25, na.rm = TRUE)
n_corrected <- sum(combined$FinaladjP_FisherCorrected <= 0.05 & abs(combined$max_logFC) > 0.25, na.rm = TRUE)
n_vote      <- sum(combined$n_methods_hit >= 3, na.rm = TRUE)

bar_df <- data.frame(
  Method = factor(c("Plain Fisher", "Corrected Fisher", "Vote (>=3 methods)"),
                   levels = c("Plain Fisher", "Corrected Fisher", "Vote (>=3 methods)")),
  n = c(n_plain, n_corrected, n_vote)
)

panelB <- ggplot(bar_df, aes(Method, n, fill = Method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.4, color = "#0b0b0b") +
  scale_fill_manual(values = c("Plain Fisher" = col_blue,
                                "Corrected Fisher" = col_aqua,
                                "Vote (>=3 methods)" = col_yellow), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B. Plain Fisher over-calls DEPs", x = NULL, y = "Proteins called") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = col_grid),
        plot.title = element_text(size = 10, face = "bold"))

# ---- Panel C: p-value QQ-plot (inflation vs. uniform null) ---------------
qq_from <- function(p, label) {
  p <- sort(p[!is.na(p) & p > 0])
  n <- length(p)
  data.frame(expected = -log10((seq_len(n) - 0.5) / n), observed = -log10(p), Method = label)
}
qq_df <- rbind(
  qq_from(combined$FinalP_Fisher, "Plain Fisher"),
  qq_from(combined$FinalP_FisherCorrected, "Corrected Fisher")
)
qq_df$Method <- factor(qq_df$Method, levels = c("Plain Fisher", "Corrected Fisher"))

panelC <- ggplot(qq_df, aes(expected, observed, color = Method)) +
  geom_abline(slope = 1, intercept = 0, color = col_muted, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 0.9, alpha = 0.5) +
  scale_color_manual(values = c("Plain Fisher" = col_blue, "Corrected Fisher" = col_aqua)) +
  labs(title = "C. Plain Fisher inflates small p-values",
       x = expression(-log[10](expected)), y = expression(-log[10](observed))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

# ---- Compose and save ------------------------------------------------------
fig2 <- (panelA | panelB | panelC) + plot_layout(widths = c(1.3, 1, 1))

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "Figure2_ConsensusValidation.pdf"), fig2, width = 13, height = 4.5)
ggsave(file.path(fig_dir, "Figure2_ConsensusValidation.png"), fig2, width = 13, height = 4.5, dpi = 300)

cat("Plain Fisher DEPs:", n_plain, " | Corrected Fisher DEPs:", n_corrected,
    " | Vote(>=3) DEPs:", n_vote, "\n")
cat("Saved to", file.path(fig_dir, "Figure2_ConsensusValidation.{pdf,png}"), "\n")
