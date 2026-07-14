# Demo: full dual-branch (count + intensity) fused interactome scoring for
# the NegKTS R467W Ind vs Mock Ind comparison.

toolkit_dir <- "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/scripts/PLToolkit"

for (f in list.files(file.path(toolkit_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

spec <- comparison_spec(
  name = "negtKTSwt_R467W_vs_Mock",
  group1_name = "negtKTSwt_R467W_Ind", group1_cols = paste0("negtKTSwt_R467W_Ind_", 1:3),
  group2_name = "Mock_Ind", group2_cols = paste0("Mock_Ind_", 1:3),
  max_na = 3
)

paths <- list(
  combined_protein_path = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Preproceed/combined_protein.csv",
  psm_dir = "/Users/XXW004/Library/CloudStorage/OneDrive-NationwideChildren'sHospital/NinaMannLab/WT1/Preproceed/",
  msstats_path = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Preproceed/MSstats.csv",
  contaminant_list_path = NULL,
  out_dir = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Results/NegKTSR467WIndvsMockInd/PLToolkitDemo/"
)

out <- run_pl_interactome(spec, paths)

cat("Count branch rows:", nrow(out$count_results), "\n")
cat("Fused table rows:", nrow(out$fused), "\n")
cat("\nEvidence tag breakdown:\n")
print(table(out$fused$Evidence))

cat("\nWT1 row (bait self-recovery check):\n")
print(out$fused[out$fused$Gene == "WT1", c("Gene", "P_count", "count_logFC",
                                             "intensity_posterior", "max_logFC",
                                             "Fused_Score", "Evidence")])

cat("\nTop 15 by Fused_Score:\n")
print(head(out$fused[, c("Gene", "P_count", "intensity_posterior", "Fused_Score", "Evidence")], 15))

cat("\nA few Count_only examples (if any):\n")
print(head(out$fused[out$fused$Evidence == "Count_only",
                       c("Gene", "P_count", "intensity_posterior", "Fused_Score", "Evidence")], 5))

cat("\nA few Intensity_only examples (if any):\n")
print(head(out$fused[out$fused$Evidence == "Intensity_only",
                       c("Gene", "P_count", "intensity_posterior", "Fused_Score", "Evidence")], 5))
