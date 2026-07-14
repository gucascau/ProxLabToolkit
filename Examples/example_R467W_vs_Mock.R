# Demo: reproduce the NegKTS R467W Ind vs Mock Ind comparison through the
# PLToolkit, for side-by-side comparison against the original monolithic Rmds.

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

out <- run_pl_dea(spec, paths)

cat("Proteins after filtering:", nrow(out$intens_raw), "\n")
cat("Methods run:", paste(names(out$results), collapse = ", "), "\n")
cat("Consensus table rows:", nrow(out$combined), "\n")
print(head(out$combined))

# setting up the adjusted P values and max_logFC for the DEPs. 
deps <- call_deps(out$combined)
cat("Consensus DEPs (Fisher p<=0.05 & |max_logFC|>0.25):", nrow(deps), "\n")
