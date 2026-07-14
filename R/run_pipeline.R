library(dplyr)

#' Run the full proximity-labeling DEP pipeline for one comparison spec:
#' load FragPipe protein table -> (optional) contaminant filter -> extract
#' intensity matrix -> log2/quantile normalize -> mixed MAR/MNAR impute ->
#' Limma/DEqMS/Wilcoxon/ROTS/DEP/proDA/(optional MSstats) -> consensus.
#'
#' `paths`: combined_protein_path (required), psm_dir (optional - skips DEqMS
#' if omitted), msstats_path (optional - skips MSstats if omitted),
#' contaminant_list_path (optional), combined_peptide_path (optional - adds a
#' min_peptides_per_sample QC column to the consensus/DEPs tables, flagging
#' proteins whose MaxLFQ intensity rests on very few peptides in some sample;
#' purely informational, doesn't filter anything), out_dir (required).
#' `options`: organism, rots_B/rots_K/rots_seed, impute_seed, contaminants
#' (a pre-loaded character vector, takes precedence over contaminant_list_path),
#' dep_p_col/dep_p_cutoff/dep_logfc_cutoff/dep_min_methods (passed straight to
#' call_deps() for the final DEPs list - see its docs in consensus.R).
run_pl_dea <- function(spec, paths, options = list()) {
  opts <- modifyList(list(
    organism = "Homo sapiens",
    rots_B = 10000, rots_K = 500, rots_seed = 1234,
    impute_seed = 1000000,
    dep_p_col = "FinaladjP_FisherCorrected", dep_p_cutoff = 0.05,
    dep_logfc_cutoff = 0.25, dep_min_methods = 3
  ), options)

  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

  pro <- load_fragpipe_protein(paths$combined_protein_path, organism = opts$organism)

  contaminants <- opts$contaminants
  if (is.null(contaminants) && !is.null(paths$contaminant_list_path)) {
    contaminants <- load_contaminant_list(paths$contaminant_list_path)
  }
  if (!is.null(contaminants)) {
    pro <- filter_contaminants(pro, contaminants)
  }

  intens_raw <- extract_intensity_matrix(pro, spec)
  norm <- log2_quantile_normalize(intens_raw)
  cond <- c(rep(spec$group1_label, length(spec$group1_cols)),
            rep(spec$group2_label, length(spec$group2_cols)))
  dat_log_exp <- impute_mixed(norm$intens_tran_normalized, cond, seed = opts$impute_seed)

  write.csv(data.frame(Gene = rownames(dat_log_exp), dat_log_exp, check.names = FALSE),
            file.path(paths$out_dir, paste0(spec$name, "__normalized_imputed_matrix.csv")),
            row.names = FALSE)

  gene_map <- pro[, c("Protein", "Gene")]

  results <- list()
  results$Limma <- run_limma(dat_log_exp, spec)
  results$Wilcoxon <- run_wilcoxon(dat_log_exp, spec)
  results$ROTS <- run_rots(dat_log_exp, spec, B = opts$rots_B, K = opts$rots_K, seed = opts$rots_seed)
  results$DEP <- run_dep_package(intens_raw, spec)
  results$proDA <- run_proda(norm$intens_tran, spec)

  if (!is.null(paths$psm_dir)) {
    sample_labels <- c(spec$group1_cols, spec$group2_cols)
    psm <- load_psm_counts(paths$psm_dir, sample_labels)
    psm <- min_psm_for_comparison(psm, sample_labels)
    psm <- psm %>% left_join(gene_map, by = "Protein")
    min_psm_by_gene <- setNames(psm$min_psm, psm$Gene)
    results$DEqMS <- run_deqms(dat_log_exp, spec, min_psm_by_gene)
  } else {
    message("paths$psm_dir not supplied - skipping DEqMS (needs per-sample psm.tsv PSM counts).")
  }

  if (!is.null(paths$msstats_path)) {
    msstats_raw <- load_msstats_raw(paths$msstats_path)
    results$MSstats <- run_msstats(msstats_raw, spec, gene_map)
  } else {
    message("paths$msstats_path not supplied - skipping MSstats.")
  }

  for (method in names(results)) {
    method_dir <- file.path(paths$out_dir, method)
    dir.create(method_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(results[[method]],
              file.path(method_dir, paste0(spec$name, "__", method, "_results.csv")),
              row.names = FALSE)
  }

  combined <- combine_pvalues(results)

  if (!is.null(paths$combined_peptide_path)) {
    sample_labels <- c(spec$group1_cols, spec$group2_cols)
    pep_support <- load_peptide_support(paths$combined_peptide_path, sample_labels)
    pep_support <- min_peptide_support_for_comparison(pep_support, sample_labels)
    combined <- combined %>%
      left_join(pep_support[, c("Gene", "min_peptides_per_sample")], by = c("ID" = "Gene"))
  } else {
    message("paths$combined_peptide_path not supplied - skipping peptide-support QC column.")
  }

  write.csv(combined,
            file.path(paths$out_dir, paste0(spec$name, "__FinalCombinedPvalue_IntegratedMethods.csv")),
            row.names = FALSE)

  deps <- call_deps(combined, p_col = opts$dep_p_col, p_cutoff = opts$dep_p_cutoff,
                     logfc_cutoff = opts$dep_logfc_cutoff, min_methods = opts$dep_min_methods)
  write.csv(deps, file.path(paths$out_dir, paste0(spec$name, "__DEPs.csv")),
            row.names = FALSE)

  list(results = results, combined = combined, deps = deps, dat_log_exp = dat_log_exp,
       intens_raw = intens_raw, protein_table = pro)
}

#' Full dual-branch interactome pipeline: runs the intensity branch
#' (run_pl_dea()) and the count branch (run_count_branch(), SAINT-style EM
#' mixture on spectral counts vs the same matched control) and fuses them
#' (fuse_branches()) into one bait-vs-control confidence score plus an
#' interpretable Count/Intensity/Both/Neither evidence tag. Also re-writes
#' run_pl_dea()'s DEPs.csv: starts from the intensity-consensus DEPs
#' (call_deps() on the intensity branch alone), then unions in any
#' Evidence == "Count_only" proteins from the fused table - genes the
#' spectral-count branch confidently calls (count_posterior >= count_cutoff)
#' but that never cleared the intensity-consensus threshold (e.g. because
#' MaxLFQ quantification is unreliable for them - see fuse_branches() docs
#' and the KMT2C example). Joins the count/intensity posteriors, Fused_Score,
#' and Evidence columns onto every row, plus a `Confidence` column: "High"
#' (Evidence == "Both"), "Count_only" (added purely on count-branch
#' evidence), "Standard" (intensity-consensus DEP without count
#' corroboration - typically Evidence == "Intensity_only", since with the
#' default aligned cutoffs an intensity-consensus DEP will always also pass
#' fuse_branches()'s own, looser intensity_positive check).
run_pl_interactome <- function(spec, paths, options = list()) {
  dea <- run_pl_dea(spec, paths, options)

  count_opts <- options$count_model %||% list()
  count_results <- do.call(run_count_branch, c(list(pro_file_update = dea$protein_table, spec = spec), count_opts))

  fused <- fuse_branches(count_results, dea$combined,
                          count_cutoff = options$count_cutoff %||% 0.75,
                          intensity_p_cutoff = options$intensity_p_cutoff %||% 0.05,
                          intensity_logfc_cutoff = options$intensity_logfc_cutoff %||% 0.25)

  write.csv(fused, file.path(paths$out_dir, paste0(spec$name, "__FusedInteractomeScore.csv")),
            row.names = FALSE)

  count_only_genes <- fused$Gene[fused$Evidence == "Count_only" & !(fused$Gene %in% dea$deps$ID)]
  extra_deps <- dea$combined %>% filter(ID %in% count_only_genes)
  deps_all <- bind_rows(dea$deps, extra_deps)

  fused_cols <- fused[, c("Gene", "count_posterior", "intensity_posterior", "Fused_Score", "Evidence")]
  deps_all <- deps_all %>% left_join(fused_cols, by = c("ID" = "Gene"))
  deps_all$Confidence <- case_when(
    !is.na(deps_all$Evidence) & deps_all$Evidence == "Both" ~ "High",
    !is.na(deps_all$Evidence) & deps_all$Evidence == "Count_only" ~ "Count_only",
    TRUE ~ "Standard"
  )
  dea$deps <- deps_all %>% arrange(desc(Fused_Score))

  write.csv(dea$deps, file.path(paths$out_dir, paste0(spec$name, "__DEPs.csv")),
            row.names = FALSE)

  c(dea, list(count_results = count_results, fused = fused))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
