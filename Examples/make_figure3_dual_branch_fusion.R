# Figure 3 - the dual-branch fusion adds independent, orthogonal
# corroboration on top of the intensity consensus, rather than just
# re-analyzing the same data twice:
# (A) the count branch (SAINT-style spectral-count mixture) is a sharply
#     bimodal, high-specificity filter - most proteins sit near "background",
#     with a small, cleanly separated "true interactor" cluster;
# (B) the intensity branch's posterior (qvalue::lfdr on the Fisher-combined
#     p-value) is also bimodal, independently of the count branch;
# (C) fusing both branches recovers the bait itself (WT1) as the strongest
#     possible positive control, and flags a small High-confidence tier
#     within the broader intensity-consensus DEP list.
#
# Reads the already-generated negtKTSwt_R467W_vs_Mock comparison outputs -
# run Examples/example_fused_R467W_vs_Mock.R first if these don't exist yet.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

out_dir <- "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Results/NegKTSR467WIndvsMockInd/PLToolkitDemo"
fig_dir <- "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/scripts/PLToolkit/Figs"
fused   <- read.csv(file.path(out_dir, "negtKTSwt_R467W_vs_Mock__FusedInteractomeScore.csv"))

# validated palette (see dataviz skill references/palette.md)
col_blue   <- "#2a78d6"  # categorical slot 1 - Both
col_aqua   <- "#1baf7a"  # categorical slot 2 - Count_only
col_yellow <- "#eda100"  # categorical slot 3 - Intensity_only
col_neither <- "#c3c2b7" # muted background ink - Neither (majority, not a series)
col_muted  <- "#898781"
count_cutoff <- 0.75

fused$Evidence <- factor(fused$Evidence,
                          levels = c("Both", "Count_only", "Intensity_only", "Neither"))
ev_colors <- c(Both = col_blue, Count_only = col_aqua,
               Intensity_only = col_yellow, Neither = col_neither)

# intensity_posterior has exact 0s (1088 of them) - log10(x) is undefined
# there, so add a small pseudocount for the log-x histogram only (Panel B).
ip_eps <- 1e-4
fused$intensity_posterior_log <- pmax(fused$intensity_posterior, ip_eps)

wt1 <- fused[fused$Gene == "WT1", ]

# ---- Panel A: P_count is sharply bimodal (background vs true interactor) --
wt1_bin_height <- sum(fused$P_count >= 0.95)  # WT1's histogram bin (rightmost)

panelA <- ggplot(fused, aes(P_count)) +
  geom_histogram(bins = 60, fill = col_blue, color = NA) +
  scale_y_log10() +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = count_cutoff, linetype = "dashed", color = col_muted, linewidth = 0.5) +
  annotate("text", x = count_cutoff, y = 2000, label = "count_cutoff = 0.75",
           hjust = 1.05, size = 3, color = col_muted) +
  geom_segment(data = wt1, aes(x = P_count, xend = P_count,
                                y = wt1_bin_height * 3, yend = wt1_bin_height * 1.15),
               arrow = arrow(length = unit(0.1, "cm")), color = "#0b0b0b", linewidth = 0.5) +
  annotate("text", x = wt1$P_count, y = wt1_bin_height * 6, label = "WT1 (bait)",
           hjust = 1.05, size = 3.2, color = "#0b0b0b") +
  labs(title = "A. Count branch: a sharp background / true-interactor split",
       x = "P_count (posterior, SAINT-style mixture)", y = "Proteins (log scale)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(t = 15, r = 10, b = 5, l = 5),
        plot.title = element_text(size = 10, face = "bold"))

# ---- Panel B: intensity_posterior is also bimodal, independently ----------
wt1_ip_bin_height <- sum(fused$intensity_posterior >= 0.95)

panelB <- ggplot(fused, aes(intensity_posterior_log)) +
  geom_histogram(bins = 60, fill = col_yellow, color = NA) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(clip = "off") +
  geom_segment(data = wt1, aes(x = intensity_posterior_log, xend = intensity_posterior_log,
                                y = wt1_ip_bin_height * 3, yend = wt1_ip_bin_height * 1.15),
               arrow = arrow(length = unit(0.1, "cm")), color = "#0b0b0b", linewidth = 0.5) +
  annotate("text", x = wt1$intensity_posterior_log, y = wt1_ip_bin_height * 6, label = "WT1 (bait)",
           hjust = 1.05, size = 3.2, color = "#0b0b0b") +
  labs(title = "B. Intensity branch is also bimodal, independently",
       x = paste0("Intensity-branch posterior (log scale, +", ip_eps, " pseudocount)"),
       y = "Proteins (log scale)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(t = 15, r = 10, b = 5, l = 5),
        plot.title = element_text(size = 10, face = "bold"))

# ---- Panel C: fused evidence scatter, WT1 as positive control -------------
panelC <- ggplot() +
  geom_jitter(data = fused[fused$Evidence == "Neither", ],
              aes(intensity_posterior, count_posterior, color = Evidence),
              width = 0.015, height = 0.015, size = 0.8, alpha = 0.25) +
  geom_jitter(data = fused[fused$Evidence != "Neither", ],
              aes(intensity_posterior, count_posterior, color = Evidence),
              width = 0.015, height = 0.015, size = 1.4, alpha = 0.75) +
  geom_point(data = wt1, aes(intensity_posterior, count_posterior),
             shape = 21, color = "#0b0b0b", fill = col_blue, size = 3, stroke = 0.8) +
  geom_text_repel(data = wt1, aes(intensity_posterior, count_posterior, label = "WT1 (bait)"),
                   nudge_y = -0.12, size = 3.2, color = "#0b0b0b",
                   segment.color = col_muted, min.segment.length = 0) +
  scale_color_manual(values = ev_colors, drop = FALSE) +
  labs(title = "C. Fusion recovers the bait and tiers the consensus DEPs",
       x = "Intensity-branch posterior", y = "Count-branch posterior (P_count)",
       color = "Evidence") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 10, face = "bold"))

# ---- Compose and save ------------------------------------------------------
fig3 <- (panelA | panelB | panelC) + plot_layout(widths = c(1, 1, 1.3))

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "Figure3_DualBranchFusion.pdf"), fig3, width = 14.5, height = 4.5)
ggsave(file.path(fig_dir, "Figure3_DualBranchFusion.png"), fig3, width = 14.5, height = 4.5, dpi = 300)

cat("Evidence breakdown:\n")
print(table(fused$Evidence))
cat("\nWT1 row:\n")
print(wt1[, c("Gene", "P_count", "intensity_posterior", "Fused_Score", "Evidence")])
cat("\nSaved to", file.path(fig_dir, "Figure3_DualBranchFusion.{pdf,png}"), "\n")
