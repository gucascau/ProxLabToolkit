# Proximity-Labeling DEP Toolkit

## Introduction

Proximity-labeling (PL) proteomics is a technique for mapping a protein's
local interaction neighborhood in living cells. A promiscuous labeling
enzyme — a biotin ligase (BioID, TurboID) or engineered peroxidase (APEX) —
is fused to a bait protein of interest and expressed in cells. Upon
activation (biotin addition, or biotin-phenol + H2O2 for APEX), the enzyme
covalently tags proteins within a short radius (roughly 10-20 nm) with
biotin. Biotinylated proteins are then enriched with streptavidin beads,
digested, and identified/quantified by mass spectrometry, yielding a
snapshot of the bait's proximal interactome — including transient and
weak interactions that co-IP/AP-MS approaches often miss — rather than a
strict binary interaction list. Because labeling is spatial and
activity-driven rather than affinity-driven, the resulting quantitative
data has analytical characteristics (background labeling, bait/enzyme
activity variability, non-random missingness) that differ from standard
AP-MS data and shape most of this toolkit's design choices, described next.

### Bioinformatic limitations of proximity-labeling (BioID/TurboID/APEX) data

Proximity-labeling experiments generate protein interaction data with several
analytical challenges that distinguish them from standard AP-MS/expression
proteomics, and that motivate most of the design choices in this toolkit:

- **Non-specific, stochastic biotinylation.** Labeling is diffusion-driven and
  not stoichiometric, so signal reflects a mix of true proximal interactors and
  background labeling of abundant bystander proteins — a stringent statistical
  threshold and control comparison (e.g. mock/untagged bait) are required
  rather than presence/absence calls.
- **Streptavidin pulldown contaminants.** Endogenously biotinylated proteins
  (carboxylases), keratins, and other common streptavidin-bead binders show up
  in every pulldown regardless of bait; these need explicit CRAPome-style
  filtering (`load_contaminant_list()` / `filter_contaminants()`) rather than
  being left for the statistical model to reject.
- **MNAR-dominated missingness.** Proteins near the labeling/detection floor
  are more likely to be missing precisely because their true abundance is low,
  not at random — standard imputation (mean/kNN alone) biases these toward the
  bulk distribution. This toolkit classifies MNAR vs MAR per protein and
  imputes them differently (`classify_mnar()` / `impute_mixed()`).
- **Variable bait expression and enzyme activity across replicates.** Batch-
  and replicate-level differences in labeling efficiency (not just true
  biology) show up as intensity shifts, so normalization
  (`log2_quantile_normalize()`) is applied before any comparison.
- **No single differential-abundance method is robust on its own.** Moderated
  t-tests, rank-based tests, and count-based tests each make different
  assumptions and are differently sensitive to the noise patterns above, so
  any one method alone under- or over-calls hits. This toolkit runs multiple
  methods (Limma, DEqMS, proDA, ROTS, DEP, MSstats, Wilcoxon) and combines
  their p-values (`combine_pvalues()`) into a consensus call, reducing
  method-specific artifacts at the cost of requiring all methods to be run
  consistently on the same input.
- **Small replicate numbers.** These experiments are typically run with few
  biological replicates, limiting statistical power — a reason consensus
  scoring across methods is preferred here over relying on strict correction
  within a single test.

## Project Description

Variants in the transcription factor WT1 (Wilms tumor 1) cause severe,
progressive glomerulopathy marked by heavy proteinuria and progression to
end-stage kidney disease (ESKD). WT1 is expressed as two major isoforms,
WT1(+)KTS and WT1(-)KTS, distinguished by the inclusion or exclusion of
three amino acids (KTS) between the third and fourth zinc fingers; the two
isoforms have both overlapping and distinct functions. WT1 is known to
regulate the expression of many essential podocyte proteins, but its
transcriptional co-activators and co-repressors remain unidentified.

This project uses proximity-dependent labeling to define the WT1
interactome and to resolve the distinct and overlapping roles of the
WT1(+)KTS and WT1(-)KTS isoforms — the toolkit in this repository
implements the downstream differential-enrichment analysis (DEA) pipeline
for that data.

<img src="Figs/ProximalLabling.png" width="1500"/>

## Directory layout

```
WT1/scripts/PLToolkit/
  R/
    load_data.R          # FragPipe protein loader + gene dedup, PSM-count loader, MSstats raw loader
    contaminants.R        # contaminant list loader + filter step
    preprocess.R           # intensity matrix extraction (0->NA, column select/reorder, NA-count filter), log2 + quantile normalize
    impute.R                # MNAR classification + validated mixed KNN/MinProb imputation
    dea_matrix_methods.R    # run_limma, run_deqms, run_wilcoxon, run_rots (share the same imputed matrix + cond vector)
    dea_dep_package.R       # run_dep_package (wraps DEP's own make_se -> filter_missval -> normalize_vsn -> impute -> test_diff pipeline)
    dea_proda.R             # run_proda
    dea_msstats.R           # run_msstats
    consensus.R              # combine_pvalues (Fisher/Tippett/Stouffer/empirical-permutation, ported from PvalueIntegration.Rmd) + call_deps threshold helper
    run_pipeline.R            # run_pl_dea(spec, paths, options) orchestrator; sources all files above
  example_R467W_vs_Mock.R     # demo: runs the full toolkit for the one comparison we already validated, for side-by-side comparison against the existing Rmd's output
  PLAN.md                     # this file
```

## Comparison spec (replaces hard-coded sample vectors)

A plain list, built by a small constructor for validation — no new file format
invented since none existed to reuse:

```r
spec <- comparison_spec(
  name        = "negtKTSwt_R467W_vs_Mock",
  group1_name = "negtKTSwt_R467W_Ind", group1_cols = paste0("negtKTSwt_R467W_Ind_", 1:3),
  group2_name = "Mock_Ind",            group2_cols = paste0("Mock_Ind_", 1:3),
  max_na      = 3
)
```
`run_pl_dea()` takes this plus a `paths` list (combined_protein.csv, psm.tsv root,
MSstats.csv, contaminant list path (optional), output dir).

`group1_label`/`group2_label` (auto-derived via `make.names()` unless overridden)
are the single consistent internal labels used for every method's design
matrix/contrast, replacing the original scripts' inconsistent per-method ad hoc
labels (e.g. DEqMS/Limma/ROTS used `"negtKTSwtR467W"` while proDA used
`"negtKTSwt_R467W"` for what was really the same group).

## Key function contracts

- `load_fragpipe_protein(path, organism = "Homo sapiens")` — generalizes the
  existing block (main Rmd lines ~835-864): filter organism, dedupe by `Gene`
  via `slice_max(rowSums(MaxLFQ cols))`, fall back to `Entry.Name` for blank
  genes, strip `_HUMAN`.
- `load_contaminant_list(path = NULL)` / `filter_contaminants(df, gene_col, contaminants)`
  — new. Accepts a user-supplied CRAPome/contaminant CSV (gene symbol column);
  ships a small built-in fallback (keratins, trypsin, albumin) clearly documented
  as a placeholder — real projects should supply their own CRAPome export. Prints
  how many rows were removed.
- `extract_intensity_matrix(df, spec)` — generalizes the 0->NA / column
  reorder / `sum(is.na(x)) < max_na` filter block (main Rmd lines ~880-892).
- `log2_quantile_normalize(mat)` — log2 + `limma::normalizeBetweenArrays`.
- `classify_mnar(mat, cond)` / `impute_mixed(mat, cond, q=0.01, tune_sigma=1, seed)`
  — the already-validated fix: knn for MAR rows, `imputeLCMD::impute.MinProb()`
  on the **full** matrix for MNAR rows (not DEP's `fun="mixed"` wrapper — confirmed
  broken for small MNAR subsets, see Context).
- `run_limma`, `run_deqms` (needs `load_psm_counts()`, generalizing the
  per-sample `psm.tsv` reading block, main Rmd lines ~569-632), `run_wilcoxon`,
  `run_rots` — all take `(dat_log_exp, cond, ...)`, return a standardized data
  frame keyed by Gene.
- `run_dep_package(intens_raw, spec)` — ports the DEP-package-native workflow
  (main Rmd lines ~1429-1520: `make_se` -> `filter_missval` -> `normalize_vsn`
  -> `DEP::impute(fun="knn")` -> `test_diff` -> `add_rejections` -> `get_results`).
  Column names off `get_results()` are condition-dependent (`"<g1>_vs_<g2>_ratio"`);
  detect them with `grep("_ratio$"/"_p\\.val$")` rather than hard-coding.
- `run_proda(intens_tran, spec)` — ports main Rmd lines ~1679-1727
  (`median_normalization` -> `proDA()` -> `test_diff`), operating on the
  log2-only (not quantile-normalized) matrix as the original does.
- `run_msstats(msstats_csv_path, spec)` — ports main Rmd lines ~1561-1641
  (`dataProcess` -> contrast matrix -> `groupComparison`), a separate raw input
  from `combined_protein.csv`.
- `combine_pvalues(results_list, methods, id_col)` — faithfully ports
  `PvalueIntegration.Rmd` lines 90-227: `inner_join` all method tables by Gene,
  `max_logFC` = largest-abs-value across a configurable set of per-method logFC
  columns, then Fisher (`metap::sumlog`)/Tippett/Stouffer (`metap::sumz`)/
  empirical-permutation (100k perms) p-value combination, each BH-adjusted.
  The existing script's default combination set is Limma+DEqMS+proDA+ROTS+Ttest+
  DEP+MSstats (Ttest included, Wilcoxon/Annova excluded) — inconsistent with the
  project's stated rationale for dropping Ttest. The ported function defaults to
  the same set for continuity but takes it as a parameter; worth the user's
  attention when they next revisit the consensus method choice.
- `call_deps(combined, p_col, logfc_col, p_cutoff=0.05, logfc_cutoff=0.25)` —
  the `p<=0.05 & abs(logFC)>0.25` threshold used everywhere in the project.
- `run_pl_dea(spec, paths, options)` — orchestrates all of the above in order,
  writes each method's CSV to `paths$out_dir/<Method>/...` (same naming
  convention as today, so `PvalueIntegration`/`FigureGeneration.Rmd`-style
  downstream consumption still works unmodified), returns the consensus table.

## Explicitly out of scope (confirmed with user)

- No R package / DESCRIPTION / roxygen — plain sourced `.R` files.
- No SAINT-style scoring — CRAPome-style filtering only.
- Figure generation (`FigureGeneration.Rmd`'s volcano/heatmap/GO-enrichment logic)
  is not being ported in this pass — a natural follow-up, not built now.
- T-test / ANOVA are not given first-class wrapper functions (project's own
  stated rationale excludes them); only noted where they still feed the
  existing consensus p-value combination.

## Verification

Ran `example_R467W_vs_Mock.R` against the real data
(`/Users/XXW004/Documents/Projects/MannNina/Project/WT1/Preproceed/combined_protein.csv`)
end-to-end and confirmed:
1. Protein counts match the already-validated manual run (2119 proteins after
   the `<4`-NA filter, ~120 classified MNAR).
2. Imputed values for known MNAR proteins (e.g. `ABRAXAS2`, `AGO3`) stay low
   (~22, near the population floor) — same spot-check already done manually.
3. Each of the 7 method functions runs without error and returns a non-empty,
   correctly-keyed results table.
4. `combine_pvalues()` output has the same column shape as
   `FinalCombinedPvalue_IntegratedMethods.csv` and a comparable (not
   necessarily identical) DEP count under the same 0.05/0.25 thresholds.
