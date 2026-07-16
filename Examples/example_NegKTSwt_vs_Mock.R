# Demo: NegKTSwt (wild-type, KTS-negative WT1) Ind vs Mock Ind comparison
# through the ProxLabToolkit. Uses wild-type WT1, unlike example_R467W_vs_Mock.R
# (a pathogenic/deficient variant) - preferred as the toolkit's flagship
# example since it reflects normal WT1 function rather than a disease variant.

toolkit_dir <- "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/scripts/ProxLabToolkit"

for (f in list.files(file.path(toolkit_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

spec <- comparison_spec(
  name = "NegKTSwt_vs_Mock",
  group1_name = "NegKTSwt_Ind", group1_cols = paste0("NegKTSwt_Ind_", 1:3),
  group2_name = "Mock_Ind", group2_cols = paste0("Mock_Ind_", 1:3),
  max_na = 3
)

paths <- list(
  combined_protein_path = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Preproceed/combined_protein.csv",
  psm_dir = "/Users/XXW004/Library/CloudStorage/OneDrive-NationwideChildren'sHospital/NinaMannLab/WT1/Preproceed/",
  msstats_path = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Preproceed/MSstats.csv",
  combined_peptide_path = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Preproceed/combined_peptide.tsv",
  contaminant_list_path = NULL,
  out_dir = "/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Results/NegKTSwtvsMockInd/ProxLabToolkitDemo/"
)

out <- run_pl_dea(spec, paths)

cat("Proteins after filtering:", nrow(out$intens_raw), "\n")
cat("Methods run:", paste(names(out$results), collapse = ", "), "\n")
cat("Consensus table rows:", nrow(out$combined), "\n")
print(head(out$combined))

# Consensus DEPs: correlation-corrected Fisher adjusted p-value + logFC cutoff
# + require at least 3 of the individual methods to also hit their own thresholds.
deps <- call_deps(out$combined, p_col = "FinaladjP_FisherCorrected",
                   logfc_cutoff = 0.25, min_methods = 3)
cat("Consensus DEPs (adj. corrected-Fisher p<=0.05 & |max_logFC|>0.25 & >=3 methods):",
    nrow(deps), "\n")
