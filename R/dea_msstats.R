library(MSstats)
library(dplyr)

#' MSstats works off FragPipe's separate MSstats.csv (feature-level) export,
#' not combined_protein.csv, and does its own quantification/normalization -
#' so it's given the raw MSstats table plus the comparison spec directly.
#' protein_gene_map should be a 2-column data frame (Protein, Gene) - e.g.
#' `pro_file_update[, c("Protein", "Gene")]` from load_fragpipe_protein() -
#' used to attach gene symbols to MSstats' raw UniProt-style protein IDs.
run_msstats <- function(msstats_raw, spec, protein_gene_map) {
  sample_cols <- c(spec$group1_cols, spec$group2_cols)
  raw <- msstats_raw %>%
    mutate(FinalSample = paste0(Condition, "_", BioReplicate)) %>%
    filter(FinalSample %in% sample_cols)

  processed <- dataProcess(raw, logTrans = 2)

  contrast_matrix <- matrix(c(-1, 1), nrow = 1)
  colnames(contrast_matrix) <- c(spec$group2_name, spec$group1_name)
  rownames(contrast_matrix) <- paste0(spec$group1_name, "_vs_", spec$group2_name)

  model <- groupComparison(contrast.matrix = contrast_matrix, processed)
  results <- model$ComparisonResult %>% arrange(pvalue)

  results <- results %>% inner_join(protein_gene_map, by = "Protein")
  results$logFC <- results$log2FC
  results$P.Value <- results$pvalue
  results$adj.P.Val <- if ("adj.pvalue" %in% colnames(results)) results$adj.pvalue else NA_real_
  results
}
