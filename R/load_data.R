library(dplyr)
library(data.table)
library(purrr)
library(tidyr)

#' Build a comparison spec used throughout the toolkit.
#' group1_label/group2_label are sanitized (make.names) versions of the group
#' names, used internally for design matrices/contrasts so every DEA method
#' uses one consistent pair of labels instead of the ad hoc per-method labels
#' the original scripts used (e.g. "negtKTSwtR467W" vs "negtKTSwt_R467W").
comparison_spec <- function(name,
                             group1_name, group1_cols,
                             group2_name, group2_cols,
                             group1_label = make.names(group1_name),
                             group2_label = make.names(group2_name),
                             max_na = 3) {
  stopifnot(length(group1_cols) > 0, length(group2_cols) > 0)
  list(
    name = name,
    group1_name = group1_name, group1_cols = group1_cols, group1_label = group1_label,
    group2_name = group2_name, group2_cols = group2_cols, group2_label = group2_label,
    max_na = max_na
  )
}

#' Read a FragPipe combined_protein.csv, restrict to one organism, and
#' collapse duplicate gene symbols to one row (the one with the highest
#' total MaxLFQ intensity), falling back to Entry.Name when Gene is blank.
load_fragpipe_protein <- function(path, organism = "Homo sapiens") {
  pro_file <- read.csv(path, header = TRUE)
  maxlfq_cols <- grep("\\.MaxLFQ\\.Intensity$", colnames(pro_file), value = TRUE)
  if (length(maxlfq_cols) == 0) {
    stop("No '.MaxLFQ.Intensity' columns found in ", path)
  }

  pro_file_update <- pro_file %>% filter(Organism == organism)
  pro_file_update <- pro_file_update %>%
    group_by(Gene) %>%
    slice_max(order_by = rowSums(across(all_of(maxlfq_cols))), n = 1, with_ties = FALSE) %>%
    ungroup()

  pro_file_update$Gene <- ifelse(is.na(pro_file_update$Gene) | pro_file_update$Gene == "",
                                  pro_file_update$Entry.Name, pro_file_update$Gene)
  pro_file_update$Gene <- gsub("_HUMAN", "", pro_file_update$Gene)

  attr(pro_file_update, "maxlfq_cols") <- maxlfq_cols
  pro_file_update
}

#' Count distinct peptides per protein per sample from each sample's psm.tsv
#' (FragPipe per-sample output), for use as DEqMS's PSM-count covariate.
#' psm_dir/<sample>/psm.tsv is expected for every name in sample_labels.
load_psm_counts <- function(psm_dir, sample_labels) {
  process_one <- function(sample_name) {
    psm <- fread(file.path(psm_dir, sample_name, "psm.tsv"))
    psm <- psm[!grepl("Contaminant|Reverse", Protein), ]
    psm_expanded <- separate_rows(psm, Protein, sep = ";")
    name <- paste0(sample_name, "_psm")
    psm_expanded %>%
      group_by(Protein) %>%
      summarise(!!name := n(), .groups = "drop")
  }

  psm_counts <- map(sample_labels, process_one)
  psm_expanded <- Reduce(function(x, y) merge(x, y, by = "Protein", all = TRUE), psm_counts)
  psm_expanded[is.na(psm_expanded)] <- 0
  psm_expanded
}

#' Minimum PSM count across a set of samples (e.g. one comparison's two
#' groups), keyed by raw FragPipe Protein ID.
min_psm_for_comparison <- function(psm_expanded, sample_cols) {
  cols <- paste0(sample_cols, "_psm")
  missing <- setdiff(cols, colnames(psm_expanded))
  if (length(missing) > 0) {
    stop("Missing PSM columns: ", paste(missing, collapse = ", "))
  }
  psm_expanded$min_psm <- apply(psm_expanded[, cols, drop = FALSE], 1, function(x) min(as.numeric(x)))
  psm_expanded
}

#' Thin wrapper around reading a FragPipe MSstats.csv export.
load_msstats_raw <- function(path) {
  readr::read_csv(path, na = c("", "NA", "0"), show_col_types = FALSE)
}
