## R/multiomics/multiomics_helpers.R
## Pure data-processing logic for the Multi-Omics module. Every multi-omics
## statistic here (DIABLO fits, SNF fusion/clustering, gene<->CpG concordance,
## pathway enrichment, the leakage-safe nested-CV benchmark) was already
## computed once by Research_05_multiomics_sexstratified's own numbered
## script pipeline (data_preparation/01-06, analyses 01-07's scripts 07-18),
## including two data-leakage bugs the pipeline's own AUDIT.md documents
## finding and fixing. This module's job is to load those already-computed
## tables (instant - flat CSV reads) and, where the underlying raw values
## support it, relabel confidence against live-adjustable thresholds - never
## to re-run DIABLO/SNF/CV or invent a new combined statistic. Mirrors
## crossomics_biomarkerconv_helpers.R's cx_bc_relabel() idiom exactly.
##
## Every loader is fail-soft (list(ok, df, error)) - never errors, never
## fabricates a result when a file is missing.

## ---------------------------------------------------------------------------
## Paths (MULTI_DATA_ROOT/MULTI_DATA_AVAILABLE/MULTI_TABLE_REGISTRY are
## defined in global.R, alongside the analogous CX_* constants)
## ---------------------------------------------------------------------------

## Generic "read one CSV off disk" loader, same contract as
## load_default_cx_table() (global.R) - by full path rather than by registry
## label, so both the Dataset tab (which reads via the registry) and every
## other sub-module (which reads a specific cell's file directly) share one
## implementation.
multi_read_table <- function(path) {
  if (is.null(path) || !MULTI_DATA_AVAILABLE || !file.exists(path)) {
    return(list(ok = FALSE, df = NULL, error = sprintf("Not available in this deployment (%s).", if (is.null(path)) "no path" else basename(path))))
  }
  df <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read table:", conditionMessage(df))))
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "This table has no rows."))
  list(ok = TRUE, df = df, error = NULL)
}

## Same as multi_read_table() but resolves a MULTI_TABLE_REGISTRY label first
## - the Dataset tab's own convenience wrapper, matching load_default_cx_table().
multi_read_registry_table <- function(label) {
  path <- MULTI_TABLE_REGISTRY[[label]]
  multi_read_table(path)
}

## ---------------------------------------------------------------------------
## Analysis "cells" - the six sex x drug/outcome cohorts the pipeline ran
## (AUDIT.md's "pipeline as executed" diagram): four drug x sex response
## cells (DIABLO+SNF, both omics) plus two drug-pooled response cells
## (DIABLO only, no SNF was run drug-pooled). Kept as one lookup so every
## sub-module offers the same six choices in the same order.
## ---------------------------------------------------------------------------

MULTI_CELLS <- list(
  list(key = "female_Adalimumab", label = "Female - Adalimumab (response)", sex = "female", drug = "Adalimumab", question = "female-Adalimumab", has_snf = TRUE),
  list(key = "male_Adalimumab",   label = "Male - Adalimumab (response)",   sex = "male",   drug = "Adalimumab", question = "male-Adalimumab",   has_snf = TRUE),
  list(key = "female_Etanercept", label = "Female - Etanercept (response)", sex = "female", drug = "Etanercept", question = "female-Etanercept", has_snf = TRUE),
  list(key = "male_Etanercept",   label = "Male - Etanercept (response)",   sex = "male",   drug = "Etanercept", question = "male-Etanercept",   has_snf = TRUE),
  list(key = "female_response",   label = "Female - drug-pooled (response)", sex = "female", drug = NA_character_, question = "female-response", has_snf = FALSE),
  list(key = "male_response",     label = "Male - drug-pooled (response)",   sex = "male",   drug = NA_character_, question = "male-response",   has_snf = FALSE)
)
MULTI_CELL_CHOICES <- setNames(vapply(MULTI_CELLS, function(c) c$key, character(1)), vapply(MULTI_CELLS, function(c) c$label, character(1)))
multi_cell_by_key <- function(key) Find(function(c) identical(c$key, key), MULTI_CELLS)

## The patient/sample matching table encodes sex as single-letter codes
## ("f"/"m"), while every analysis-result table (DIABLO/SNF/concordance/...)
## spells it out ("female"/"male") - normalize before comparing the two, so
## a cell filter against the matching table doesn't silently match zero rows.
multi_norm_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("f", "female"), "female", ifelse(x %in% c("m", "male"), "male", x))
}

## Table36 (master six-part summary) mixes two different row shapes in one
## table: sex-specific integrated-model rows (sex = female/male) and
## drug-pooled single-omics *baseline* rows (sex = "pooled") - a plain
## sex+drug filter would wrongly mix cells together (both response outcomes
## for a sex share the same sex value), so cell selection instead matches
## Table36's own descriptive `question` column, and the single-omics
## baseline is fetched separately by drug only.
multi_table36_for_cell <- function(df, cell) {
  if (is.null(df) || !"question" %in% colnames(df)) return(NULL)
  df[df$question %in% cell$question, , drop = FALSE]
}
multi_table36_single_omics_baseline <- function(df, cell) {
  if (is.null(df) || is.na(cell$drug) || !all(c("sex", "drug") %in% colnames(df))) return(NULL)
  df[df$sex %in% "pooled" & df$drug %in% cell$drug, , drop = FALSE]
}

## Subsets a table to one cell by its `sex`/`drug` columns, when present.
## `drug = NA` (the drug-pooled response cells) matches rows where the table
## simply has no `drug` column at all, or leaves it untouched otherwise -
## never drops rows on a column that isn't there.
multi_filter_cell <- function(df, sex = NULL, drug = NULL) {
  if (is.null(df)) return(df)
  out <- df
  if (!is.null(sex) && "sex" %in% colnames(out)) out <- out[multi_norm_sex(out$sex) %in% multi_norm_sex(sex), , drop = FALSE]
  if (!is.null(drug) && !is.na(drug) && "drug" %in% colnames(out)) out <- out[out$drug %in% drug, , drop = FALSE]
  out
}

## ---------------------------------------------------------------------------
## Sample harmonization summary (spec: "Transcriptomics: 80 samples /
## Methylomics: 65 samples / Shared samples: 58" worked example) - computed
## directly from the pipeline's own patient/sample matching table
## (metadata/patient_sample_matching_table.csv), never silently merged.
## ---------------------------------------------------------------------------

multi_sample_harmonization <- function(matching_df) {
  need_cols <- c("RNA_available_PBMC", "methylation_available")
  if (is.null(matching_df) || !all(need_cols %in% colnames(matching_df))) {
    return(list(ok = FALSE, error = "Patient/sample matching table is missing the expected availability columns."))
  }
  rna_ok  <- matching_df$RNA_available_PBMC %in% c(TRUE, "TRUE", "Yes", "yes", 1)
  meth_ok <- matching_df$methylation_available %in% c(TRUE, "TRUE", "Yes", "yes", 1)
  list(
    ok = TRUE,
    n_total = nrow(matching_df),
    n_rna = sum(rna_ok),
    n_meth = sum(meth_ok),
    n_matched = sum(rna_ok & meth_ok),
    n_rna_only = sum(rna_ok & !meth_ok),
    n_meth_only = sum(!rna_ok & meth_ok),
    n_neither = sum(!rna_ok & !meth_ok),
    by_sex = tryCatch(as.data.frame(table(sex = multi_norm_sex(matching_df$sex[rna_ok & meth_ok]))), error = function(e) NULL),
    by_cell = do.call(rbind, lapply(MULTI_CELLS, function(cl) {
      sub <- matching_df[multi_norm_sex(matching_df$sex) %in% multi_norm_sex(cl$sex), , drop = FALSE]
      if (!is.na(cl$drug) && "treatment" %in% colnames(sub)) sub <- sub[sub$treatment %in% cl$drug, , drop = FALSE]
      r_ok <- sub$RNA_available_PBMC %in% c(TRUE, "TRUE", "Yes", "yes", 1)
      m_ok <- sub$methylation_available %in% c(TRUE, "TRUE", "Yes", "yes", 1)
      data.frame(cell = cl$label, n_total = nrow(sub), n_rna = sum(r_ok), n_methylation = sum(m_ok), n_matched = sum(r_ok & m_ok))
    }))
  )
}

## ---------------------------------------------------------------------------
## Biomarker candidate relabeling (Table40/Table44b already carry
## performance_ci_lo/performance_ci_hi/biomarker_status per the source
## pipeline's own nested-CV benchmark - this only recomputes a *display*
## confidence tier from those retained raw values against user-adjustable
## thresholds, exactly like cx_bc_relabel(); it never computes a new p-value
## or combines evidence across rows into a fabricated joint statistic).
## ---------------------------------------------------------------------------

MULTI_BIOMARKER_DEFAULT_PARAMS <- list(min_auroc = 0.6, require_excludes_chance = TRUE)

multi_biomarker_relabel <- function(df, params = MULTI_BIOMARKER_DEFAULT_PARAMS) {
  p <- utils::modifyList(MULTI_BIOMARKER_DEFAULT_PARAMS, params %||% list())
  out <- df
  excludes_chance <- !is.na(out$performance_ci_lo) & out$performance_ci_lo > 0.5
  meets_auroc <- !is.na(out$performance_auroc) & out$performance_auroc >= p$min_auroc
  out$display_confidence <- ifelse(
    meets_auroc & (!p$require_excludes_chance | excludes_chance), "High confidence (panel-level)",
    ifelse(meets_auroc, "Moderate (AUROC meets threshold, CI includes chance)", "Below threshold")
  )
  out$panel_excludes_chance <- excludes_chance
  out
}

## ---------------------------------------------------------------------------
## Gene<->CpG concordance - small tally of the pipeline's own
## `biological_pattern`/`region` classification (Table42/Table45), for a
## quick bar chart. No new classification logic - purely a count of an
## already-computed categorical column.
## ---------------------------------------------------------------------------

multi_concordance_pattern_tally <- function(df) {
  if (is.null(df) || !"biological_pattern" %in% colnames(df)) return(NULL)
  as.data.frame(table(pattern = df$biological_pattern), responseName = "n")
}

## ---------------------------------------------------------------------------
## Cross-cell performance comparison (Table36's own "integrated vs
## single-omics, per cell" rows) - answers spec's "show whether integration
## actually improves performance, don't assume it does" requirement directly
## from the pipeline's own numbers, no new statistical test performed here.
## ---------------------------------------------------------------------------

multi_performance_by_omics <- function(table36_df) {
  if (is.null(table36_df)) return(NULL)
  table36_df[order(table36_df$sex, table36_df$drug, table36_df$omics), ,
             drop = FALSE]
}

## ---------------------------------------------------------------------------
## DIABLO variance explained per component, per omics block - reads the
## saved block.splsda fit's own `prop_expl_var` slot directly (real
## mixOmics::block.splsda output, computed once by the pipeline) via plain
## readRDS() + list-indexing; never re-fit, never calls a mixOmics function,
## so library(mixOmics) is not required for this module to show it.
## ---------------------------------------------------------------------------

MULTI_BLOCK_LABELS <- c(expression = "Transcriptomics", methylation = "Methylomics")

multi_diablo_fit <- function(cell) {
  path <- MULTI_DIABLO_FIT_REGISTRY[[cell$key]]
  if (is.null(path) || !MULTI_DATA_AVAILABLE || !file.exists(path)) {
    return(list(ok = FALSE, fit = NULL, error = "No saved DIABLO fit for this cell in this deployment."))
  }
  fit <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(fit, "error")) return(list(ok = FALSE, fit = NULL, error = paste("Could not read the saved fit:", conditionMessage(fit))))
  list(ok = TRUE, fit = fit, error = NULL)
}

## Tidies fit$prop_expl_var (a named list of named numeric vectors, one per
## block including the internal "Y" outcome block) into one data.frame,
## restricted to the real omics blocks (excludes "Y", which is DIABLO's own
## outcome-indicator block, not an omics layer).
multi_diablo_variance_df <- function(fit) {
  pv <- fit$prop_expl_var
  if (is.null(pv)) return(NULL)
  blocks <- setdiff(names(pv), "Y")
  if (length(blocks) == 0) return(NULL)
  do.call(rbind, lapply(blocks, function(b) {
    v <- pv[[b]]
    data.frame(block = MULTI_BLOCK_LABELS[[b]] %||% b, component = factor(names(v), levels = names(v)), variance_explained = as.numeric(v))
  }))
}

## ---------------------------------------------------------------------------
## Dynamic QC scorecard (Overview tab) - every status computed from real
## session/deployment state, never hardcoded to "pass". `multi_results` is
## the shared reactiveValues every sub-module publishes its own loaded state
## into (see each mod_multi_*.R's trailing observe() block).
## ---------------------------------------------------------------------------

multi_qc_scorecard <- function(multi_results) {
  r <- multi_results %||% list()
  overview <- r$overview
  item <- function(label, status, detail) list(label = label, status = status, detail = detail)
  list(
    item("Data availability", if (MULTI_DATA_AVAILABLE) "pass" else "fail",
         if (MULTI_DATA_AVAILABLE) "Research_05_multiomics_sexstratified is available in this deployment." else "The multi-omics pipeline output folder is not available."),
    item("Sample harmonization", if (!is.null(overview$harmonization) && isTRUE(overview$harmonization$ok)) "pass" else "warn",
         if (!is.null(overview$harmonization) && isTRUE(overview$harmonization$ok)) sprintf("%s of %s patients matched across both omics layers.", format(overview$harmonization$n_matched, big.mark = ","), format(overview$harmonization$n_total, big.mark = ",")) else "Load cohort tables on the Overview tab to compute this."),
    item("Model performance (honesty check)",
         if (!is.null(overview$summary36)) {
           n_ex <- sum(overview$summary36$excludes_chance %in% TRUE)
           if (n_ex == 0) "warn" else "pass"
         } else "warn",
         if (!is.null(overview$summary36)) sprintf("%d of %d method x cell results exclude chance performance.", sum(overview$summary36$excludes_chance %in% TRUE), nrow(overview$summary36)) else "Load cohort tables on the Overview tab to compute this."),
    item("Integration cell loaded", if (!is.null(r$integration)) "pass" else "warn",
         if (!is.null(r$integration)) sprintf("Loaded: %s", r$integration$cell$label) else "No cell loaded yet on the Integration tab."),
    item("Patient stratification loaded", if (!is.null(r$stratification)) "pass" else "warn",
         if (!is.null(r$stratification)) sprintf("Loaded: %s SNF clusters", r$stratification$drug) else "No SNF clusters loaded yet on the Stratification tab."),
    item("Joint biomarker candidates loaded", if (!is.null(r$biomarker)) "pass" else "warn",
         if (!is.null(r$biomarker)) sprintf("%d candidate features loaded.", length(unique(r$biomarker$df$feature))) else "No candidates loaded yet on the Biomarker Discovery tab."),
    item("Pathway enrichment loaded", if (!is.null(r$pathway)) "pass" else "warn",
         if (!is.null(r$pathway)) sprintf("%d enriched terms loaded.", nrow(r$pathway$df)) else "No pathway table loaded yet on the Pathway tab.")
  )
}

## ---------------------------------------------------------------------------
## Analysis Summary table (Overview tab) - "Parameter | Result", built only
## from values that are genuinely known this session; unloaded sub-modules
## show "Not loaded", never a fabricated value.
## ---------------------------------------------------------------------------

multi_analysis_summary_table <- function(multi_dataset, multi_results) {
  r <- multi_results %||% list()
  overview <- r$overview
  row <- function(parameter, result) data.frame(Parameter = parameter, Result = result, stringsAsFactors = FALSE)
  rbind(
    row("Omics layers", "Transcriptomics + Methylomics"),
    row("Samples analyzed (matched)", if (!is.null(overview$harmonization)) format(overview$harmonization$n_matched, big.mark = ",") else "Not loaded"),
    row("Active dataset table", if (!is.null(multi_dataset$table_label)) multi_dataset$table_label else "None loaded"),
    row("Integration cell", if (!is.null(r$integration)) r$integration$cell$label else "Not loaded"),
    row("Integration method(s)", if (!is.null(r$integration)) paste(c("DIABLO", if (!is.null(r$integration$snf_perf)) "SNF"), collapse = " + ") else "Not loaded"),
    row("Biomarker confidence thresholds", if (!is.null(r$biomarker)) "See Biomarker Discovery tab (user-adjustable)" else "Not loaded"),
    row("Pathway terms shown", if (!is.null(r$pathway)) format(nrow(r$pathway$df), big.mark = ",") else "Not loaded"),
    row("SNF cohort loaded", if (!is.null(r$stratification)) r$stratification$drug else "Not loaded")
  )
}

## ---------------------------------------------------------------------------
## Concordance table: optional BH-FDR recompute over the raw nominal
## expr_p/meth_p columns already retained in Table42/Table45 - same
## relabel-from-retained-raw-values idiom cx_bc_relabel()'s
## mqtl_sig_basis = "fdr" option already uses; never a fabricated value.
## ---------------------------------------------------------------------------

multi_concordance_add_fdr <- function(df) {
  if (is.null(df)) return(df)
  out <- df
  if ("expr_p" %in% colnames(out)) out$expr_fdr <- stats::p.adjust(out$expr_p, method = "BH")
  if ("meth_p" %in% colnames(out)) out$meth_fdr <- stats::p.adjust(out$meth_p, method = "BH")
  out
}

## ---------------------------------------------------------------------------
## Reproducibility: real installed package versions, never a fabricated
## version string - packages not installed show "not installed", not a guess.
## ---------------------------------------------------------------------------

multi_package_versions <- function() {
  pkgs <- c("mixOmics", "SNFtool", "MOFA2", "reticulate", "limma", "sva", "clusterProfiler")
  data.frame(
    Package = pkgs,
    Version = vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"), character(1)),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## Results Summary & Reproducibility - what's genuinely NOT implemented in
## this module (spec: never render a fake placeholder as if it were a
## result - state it plainly instead), and a pointer to the real pipeline
## scripts that produced whatever is currently loaded, for reproducibility.
## ---------------------------------------------------------------------------

MULTI_KNOWN_LIMITATIONS <- c(
  "MOFA2 was never run on the precomputed cohort by the source pipeline - DIABLO/SNF/concordance/pathway tabs stay precomputed-only for that cohort. A real MOFA2 fit is available for YOUR OWN uploaded data on the Live Analysis tab (matched samples only, run asynchronously so it doesn't freeze the app).",
  "The Live Analysis tab's DIABLO-style supervised classification is not offered for uploaded data - only unsupervised MOFA2 - since a properly leakage-safe nested-CV supervised fit (like the precomputed pipeline's own) is a materially larger undertaking than this delivery scopes; do not read MOFA2 factors as a diagnostic classifier.",
  "No dedicated held-out validation sub-module: reported AUROCs (precomputed tabs) are the pipeline's own leave-one-out / nested cross-validation performance, not performance on an independent replication cohort.",
  "No dedicated network/circos visualization beyond the gene<->CpG concordance scatter, DIABLO panel plots, and the Live Analysis correlation heatmap.",
  "Live Analysis supports up to 4 uploaded omics layers and matched-sample integration only (no partial-sample MOFA2 in this delivery).",
  "Most precomputed analysis cells' AUROC 95% CIs include chance performance (see the Overview tab) - treat those outputs as exploratory hypothesis generation, not a validated diagnostic panel."
)

MULTI_REPRODUCIBILITY_SCRIPTS <- c(
  "Shared upstream (sample matching, QC, normalization, annotation, genome-wide discovery, leakage-safe nested-CV benchmark): analyses/data_preparation/scripts/01-06_*.R",
  "SNF integration (classifier + unsupervised joint-biomarker discovery): analyses/01_female_male_adalimumab/scripts/07*.R, analyses/02_female_male_etanercept/scripts/07*.R",
  "DIABLO integration (drug x sex, response, drug-type): analyses/01_female_male_adalimumab/scripts/08_*.R, analyses/07_cross_analysis_summary/scripts/11_*.R, 14_*.R",
  "Gene<->CpG concordance: analyses/07_cross_analysis_summary/scripts/15_*.R, 18_*.R",
  "Pathway enrichment: analyses/07_cross_analysis_summary/scripts/16_*.R",
  "Cross-cell summary assembly: analyses/07_cross_analysis_summary/scripts/09_*.R, 12_*.R",
  "Independent audit (methodology, leakage findings/fixes, honest AUROC verdicts per cell): AUDIT.md"
)

## Plain-text report for the "Download everything loaded so far" bundle -
## lists which sub-modules had results loaded this session and where their
## numbers actually came from, not a regenerated analysis.
multi_build_report <- function(multi_results) {
  ids <- c("overview", "integration", "stratification", "biomarker", "concordance", "pathway", "live_qc", "live_mofa")
  lines <- c(
    "# ArthOMix Multi-Omics module - session report",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "This module browses Research_05_multiomics_sexstratified's own precomputed DIABLO/SNF/concordance/pathway results - nothing below was recomputed by this app.",
    ""
  )
  for (i in ids) {
    res <- multi_results[[i]]
    lines <- c(lines, sprintf("## %s", i), if (is.null(res)) "(not loaded this session)" else "Loaded - see the accompanying CSV(s) in this bundle for the exact rows shown.", "")
  }
  lines <- c(lines, "## Known limitations", paste0("- ", MULTI_KNOWN_LIMITATIONS), "",
             "## Reproducibility - source scripts", paste0("- ", MULTI_REPRODUCIBILITY_SCRIPTS))
  lines
}
