library(dplyr)

#' Run the full proximity-labeling DEP pipeline for one comparison spec:
#' load FragPipe protein table -> (optional) contaminant filter -> extract
#' intensity matrix -> log2/quantile normalize -> mixed MAR/MNAR impute ->
#' Limma/DEqMS/Wilcoxon/ROTS/DEP/proDA/(optional MSstats) -> consensus.
#'
#' `paths`: combined_protein_path (required), psm_dir (optional - skips DEqMS
#' if omitted), msstats_path (optional - skips MSstats if omitted),
#' contaminant_list_path (optional), out_dir (required).
#' `options`: organism, rots_B/rots_K/rots_seed, impute_seed, contaminants
#' (a pre-loaded character vector, takes precedence over contaminant_list_path).
run_pl_dea <- function(spec, paths, options = list()) {
  opts <- modifyList(list(
    organism = "Homo sapiens",
    rots_B = 10000, rots_K = 500, rots_seed = 1234,
    impute_seed = 1000000
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
  write.csv(combined,
            file.path(paths$out_dir, paste0(spec$name, "__FinalCombinedPvalue_IntegratedMethods.csv")),
            row.names = FALSE)

  list(results = results, combined = combined, dat_log_exp = dat_log_exp,
       intens_raw = intens_raw, protein_table = pro)
}

#' Full dual-branch interactome pipeline: runs the intensity branch
#' (run_pl_dea(), unchanged) and the count branch (run_count_branch(), SAINT-
#' style EM mixture on spectral counts vs the same matched control) and fuses
#' them (fuse_branches()) into one bait-vs-control confidence score plus an
#' interpretable Count/Intensity/Both/Neither evidence tag. Does not change
#' run_pl_dea()'s own behavior/output - this is an additive wrapper around it.
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

  c(dea, list(count_results = count_results, fused = fused))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
