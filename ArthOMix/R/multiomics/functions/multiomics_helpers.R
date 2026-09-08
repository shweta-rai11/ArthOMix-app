## R/multiomics/functions/multiomics_helpers.R
## Shared low-level utilities for the Multi-Omics module: cell/cohort metadata, CSV loaders, sample-harmonization.

multi_read_table <- function(path) {
  if (is.null(path) || !MULTI_DATA_AVAILABLE || !file.exists(path)) {
    return(list(ok = FALSE, df = NULL, error = sprintf("Not available in this deployment (%s).", if (is.null(path)) "no path" else basename(path))))
  }
  df <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read table:", conditionMessage(df))))
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "This table has no rows."))
  list(ok = TRUE, df = df, error = NULL)
}

multi_read_registry_table <- function(label) {
  path <- MULTI_TABLE_REGISTRY[[label]]
  multi_read_table(path)
}

MULTI_CELLS <- list(
  list(key = "female_Adalimumab", label = "Female - Adalimumab (response)", sex = "female", drug = "Adalimumab", question = "female-Adalimumab", has_snf = FALSE),
  list(key = "male_Adalimumab",   label = "Male - Adalimumab (response)",   sex = "male",   drug = "Adalimumab", question = "male-Adalimumab",   has_snf = FALSE),
  list(key = "female_Etanercept", label = "Female - Etanercept (response)", sex = "female", drug = "Etanercept", question = "female-Etanercept", has_snf = TRUE),
  list(key = "male_Etanercept",   label = "Male - Etanercept (response)",   sex = "male",   drug = "Etanercept", question = "male-Etanercept",   has_snf = TRUE),
  list(key = "female_response",   label = "Female - drug-pooled (response)", sex = "female", drug = NA_character_, question = "female-response", has_snf = FALSE),
  list(key = "male_response",     label = "Male - drug-pooled (response)",   sex = "male",   drug = NA_character_, question = "male-response",   has_snf = FALSE)
)
MULTI_CELL_CHOICES <- setNames(vapply(MULTI_CELLS, function(c) c$key, character(1)), vapply(MULTI_CELLS, function(c) c$label, character(1)))
multi_cell_by_key <- function(key) Find(function(c) identical(c$key, key), MULTI_CELLS)

multi_norm_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("f", "female"), "female", ifelse(x %in% c("m", "male"), "male", x))
}

multi_filter_cell <- function(df, sex = NULL, drug = NULL) {
  if (is.null(df)) return(df)
  out <- df
  if (!is.null(sex) && "sex" %in% colnames(out)) out <- out[multi_norm_sex(out$sex) %in% multi_norm_sex(sex), , drop = FALSE]
  if (!is.null(drug) && !is.na(drug) && "drug" %in% colnames(out)) out <- out[out$drug %in% drug, , drop = FALSE]
  out
}

multi_sex_candidates <- function(sample_meta) {
  if (is.null(sample_meta) || ncol(sample_meta) == 0) return(character(0))
  colnames(sample_meta)[grepl("sex|gender", colnames(sample_meta), ignore.case = TRUE)]
}

multi_sex_groups <- function(sample_meta, sex_col, sample_ids) {
  if (is.null(sample_meta) || is.null(sex_col) || !sex_col %in% colnames(sample_meta)) return(NULL)
  vals <- stats::setNames(as.character(sample_meta[[sex_col]]), rownames(sample_meta))
  vals <- vals[intersect(sample_ids, names(vals))]
  vals <- vals[!is.na(vals) & nzchar(vals)]
  split(names(vals), vals)
}

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

multi_diablo_variance_df <- function(fit) {
  pv <- fit$prop_expl_var
  if (is.null(pv)) return(NULL)
  blocks <- setdiff(names(pv), "Y")
  if (length(blocks) == 0) return(NULL)
  do.call(rbind, lapply(blocks, function(b) {
    v <- pv[[b]]
    data.frame(block = b, component = factor(names(v), levels = names(v)), variance_explained = as.numeric(v))
  }))
}

## Direction-aware classification of an AUROC + its 95% CI. "excludes_chance"
## (CI does not contain 0.5) is a statement about *statistical significance
## only* - it says nothing about whether performance is genuine (above chance)
## or actually WORSE than chance (below 0.5), and treating both directions as
## an undifferentiated "pass" mislabels a classifier that performs
## significantly worse than random guessing as a validation success. This is
## the single source of truth for that distinction; every scorecard/UI
## consumer must use it rather than re-deriving pass/fail from excludes_chance.
##   "above_chance"  - CI entirely above 0.5: genuine, statistically supported discrimination.
##   "below_chance"  - CI entirely below 0.5: statistically supported, but WORSE than random -
##                      never a validation success; typically indicates label/orientation
##                      error or severe overfitting to a spuriously anti-correlated small sample.
##   "chance"        - CI spans 0.5: not distinguishable from random guessing.
##   "unsupported"   - AUROC/CI missing or non-finite: cannot be classified.
multi_auroc_call <- function(auroc, ci_lo, ci_hi) {
  if (is.null(auroc) || length(auroc) == 0 || is.na(auroc) ||
      is.null(ci_lo) || length(ci_lo) == 0 || is.na(ci_lo) ||
      is.null(ci_hi) || length(ci_hi) == 0 || is.na(ci_hi)) {
    return("unsupported")
  }
  if (ci_lo > 0.5) return("above_chance")
  if (ci_hi < 0.5) return("below_chance")
  "chance"
}

## Vectorised convenience for a data.frame of auroc/ci_lo/ci_hi columns (e.g. the
## precomputed Table36 registry table, which is read as-is and must never be
## trusted to already carry a direction-aware flag).
multi_auroc_call_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(character(0))
  vapply(seq_len(nrow(df)), function(i) multi_auroc_call(df$auroc[i], df$ci_lo[i], df$ci_hi[i]), character(1))
}

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
           calls <- multi_auroc_call_df(overview$summary36)
           n_below <- sum(calls == "below_chance")
           n_above <- sum(calls == "above_chance")
           ## A genuinely below-chance result is a red flag in its own right (label/orientation
           ## error or severe overfitting), never a "pass" merely because SOME result excludes
           ## chance - it must dominate the status even if an above-chance result also exists.
           if (n_below > 0) "fail" else if (n_above > 0) "pass" else "warn"
         } else "warn",
         if (!is.null(overview$summary36)) {
           calls <- multi_auroc_call_df(overview$summary36)
           n_below <- sum(calls == "below_chance"); n_above <- sum(calls == "above_chance")
           if (n_below > 0) {
             sprintf("%d of %d method x cell results are SIGNIFICANTLY BELOW CHANCE (AUROC CI entirely under 0.5) - not a validation success, a red flag (label/orientation error or overfitting). %d genuinely exclude chance in the expected direction.", n_below, nrow(overview$summary36), n_above)
           } else {
             sprintf("%d of %d method x cell results genuinely exclude chance performance (AUROC CI entirely above 0.5).", n_above, nrow(overview$summary36))
           }
         } else "Load cohort tables on the Overview tab to compute this."),
    item("Integration cell loaded", if (!is.null(r$integration)) "pass" else "warn",
         if (!is.null(r$integration)) sprintf("Loaded: %s", r$integration$cell$label) else "No cell loaded yet on the Integration tab."),
    item("Sex-stratified DIABLO loaded", if (!is.null(r$integration_stratified)) "pass" else "warn",
         if (!is.null(r$integration_stratified)) sprintf("Loaded: %s", r$integration_stratified$cell$label) else "No sex-stratified comparison loaded yet on the Integration tab."),
    item("Patient stratification loaded", if (!is.null(r$stratification)) "pass" else "warn",
         if (!is.null(r$stratification)) sprintf("Loaded: %s SNF clusters", r$stratification$drug) else "No SNF clusters loaded yet on the Stratification tab."),
    item("Biomarker Discovery signature loaded", if (!is.null(r$biomarker)) "pass" else "warn",
         if (!is.null(r$biomarker)) sprintf("%d selected features loaded.", length(unique(r$biomarker$df$feature))) else "No signature loaded yet on the Biomarker Discovery tab."),
    item("Pathway enrichment loaded", if (!is.null(r$pathway)) "pass" else "warn",
         if (!is.null(r$pathway)) sprintf("%d enriched terms loaded.", nrow(r$pathway$df)) else "No pathway table loaded yet on the Pathway tab.")
  )
}

multi_analysis_summary_table <- function(multi_dataset, multi_results) {
  r <- multi_results %||% list()
  overview <- r$overview
  row <- function(parameter, result) data.frame(Parameter = parameter, Result = result, stringsAsFactors = FALSE)
  active_layers <- if (!is.null(multi_dataset$layers) && length(multi_dataset$layers) > 0) paste(names(multi_dataset$layers), collapse = " + ") else "Transcriptomics + Methylomics (precomputed cohort)"
  rbind(
    row("Omics layers", active_layers),
    row("Samples analyzed (matched)", if (!is.null(overview$harmonization)) format(overview$harmonization$n_matched, big.mark = ",") else "Not loaded"),
    row("Active dataset table", if (!is.null(multi_dataset$table_label)) multi_dataset$table_label else "None loaded"),
    row("Active Multi-Omics Dataset source", if (isTRUE(multi_dataset$active)) switch(multi_dataset$source %||% "", preloaded = "Preloaded Dataset", upload = "User Upload", geo = "NCBI GEO", "Unknown") else "None selected yet"),
    row("Integration cell", if (!is.null(r$integration)) r$integration$cell$label else "Not loaded"),
    row("Integration method(s)", if (!is.null(r$integration)) paste(c("DIABLO", if (!is.null(r$integration$snf_perf)) "SNF"), collapse = " + ") else "Not loaded"),
    row("Sex-stratified DIABLO comparison", if (!is.null(r$integration_stratified)) sprintf("Loaded: %s", r$integration_stratified$cell$label) else "Not loaded"),
    row("Biomarker Discovery signature", if (!is.null(r$biomarker)) sprintf("%d selected features (DIABLO) - see Biomarker Discovery tab", length(unique(r$biomarker$df$feature))) else "Not loaded"),
    row("Pathway terms shown", if (!is.null(r$pathway)) format(nrow(r$pathway$df), big.mark = ",") else "Not loaded"),
    row("SNF cohort loaded", if (!is.null(r$stratification)) r$stratification$drug else "Not loaded")
  )
}

multi_mapping_add_fdr <- function(df) {
  if (is.null(df)) return(df)
  out <- df
  if ("expr_p" %in% colnames(out)) out$expr_fdr <- stats::p.adjust(out$expr_p, method = "BH")
  if ("meth_p" %in% colnames(out)) out$meth_fdr <- stats::p.adjust(out$meth_p, method = "BH")
  out
}

multi_active_dataset_banner <- function(multi_dataset, multi_results = NULL) {
  md <- multi_dataset %||% list()
  source <- md$source
  progress <- if (!is.null(multi_results)) multi_pipeline_progress_ui(multi_dataset, multi_results) else NULL
  if (is.null(source) || !isTRUE(md$active %||% FALSE)) {
    return(div(class = "empty-note", icon("circle-info"),
               "No active Multi-Omics dataset yet. Pick one on the Dataset tab, or choose Preloaded/Reference here to use the bundled RA anti-TNF cohort.",
               progress))
  }
  if (identical(source, "preloaded")) {
    return(div(class = "empty-note", icon("circle-check"),
               tags$strong("Data source: Preloaded Dataset."), " Existing results are available and shown below.", progress))
  }
  div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
      tags$strong(sprintf("Data source: %s.", if (identical(source, "geo")) "NCBI GEO" else "User Upload")),
      " No stored results for this dataset yet. Select \"Active Multi-Omics Dataset\" above to run this module directly.",
      progress)
}

multi_package_versions <- function() {
  pkgs <- c("mixOmics", "SNFtool", "MOFA2", "reticulate", "limma", "sva", "clusterProfiler",
            "ReactomePA", "fgsea", "org.Hs.eg.db", "KEGGREST", "msigdbr", "pathview")
  data.frame(
    Package = pkgs,
    Version = vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"), character(1)),
    stringsAsFactors = FALSE
  )
}


## Pipeline hand-offs: Cohort Harmonisation -> Integration (fused network) -> Biomarker Discovery/SNF Clustering -> Gene-CpG Mapping -> Pathways/Biomarker Card.

multi_harmonisation_state <- function(multi_results) {
  if (is.null(multi_results)) return(NULL)
  h <- tryCatch(multi_results$overview$harmonization, error = function(e) NULL)
  if (is.null(h) || !isTRUE(h$ok)) return(NULL)
  h
}

## One-line status of the Cohort Harmonisation hand-off for a downstream sub-module.
multi_harmonisation_note <- function(multi_results, n_here = NULL) {
  h <- multi_harmonisation_state(multi_results)
  if (is.null(h)) {
    return(div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
               " Cohort Harmonisation hasn't been run on this dataset, so readiness, identifier audit and phenotype classification aren't available here. Run it first (Sub-modules > Cohort Harmonisation)."))
  }
  rd <- h$readiness %||% list(label = "Not assessed", reason = "")
  color <- switch(rd$level %||% "unknown", ready = ARTHOMIX_COLORS$aqua, limited = ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)
  mismatch <- !is.null(n_here) && !is.na(h$n_matched) && n_here != h$n_matched
  div(class = "empty-note", icon("link"),
      tags$strong(" Cohort Harmonisation: "), span(style = sprintf("color:%s; font-weight:600;", color), rd$label),
      sprintf(" - %s", rd$reason %||% ""),
      sprintf(" %d of %d identifiers matched across %s; %d unmatched/duplicate/ambiguous.", h$n_matched, h$n_total, paste(h$modalities, collapse = " + "), h$n_unmatched_ids %||% 0L),
      if (length(h$phenotype_candidates) > 0) sprintf(" Phenotype column(s) classified: %s.", paste(h$phenotype_candidates, collapse = ", ")),
      if (mismatch) tags$span(style = "color: var(--color-danger, #d9534f);", sprintf(" This analysis currently matches %d samples, which differs from the harmonised %d - check the modality selection.", n_here, h$n_matched)))
}

## Default outcome column: first offered phenotype candidate, else first offered column.
multi_harmonisation_outcome_default <- function(multi_results, candidates) {
  if (length(candidates) == 0) return(NULL)
  h <- multi_harmonisation_state(multi_results)
  hit <- if (!is.null(h)) intersect(h$phenotype_candidates %||% character(0), candidates) else character(0)
  if (length(hit) > 0) hit[1] else candidates[1]
}

## Stage-by-stage status of the Multi-Omics pipeline for the shared banner.
MULTI_PIPELINE_STAGES <- c("Dataset Workspace", "Cohort Harmonisation", "Integration (DIABLO / SNF)", "SNF Clustering",
                           "Biomarker Discovery", "Gene-CpG Mapping", "Pathways", "Biomarker Card")

multi_pipeline_progress <- function(multi_dataset, multi_results) {
  md <- multi_dataset %||% list(); mr <- multi_results
  g <- function(expr) tryCatch(expr, error = function(e) NULL)
  done <- c(
    isTRUE(md$active) && length(md$layers %||% list()) >= 2,
    !is.null(multi_harmonisation_state(mr)),
    !is.null(g(mr$integration$perf)) || !is.null(g(mr$integration$snf)),
    !is.null(g(mr$stratification$clusters)),
    !is.null(g(mr$biomarker$df)),
    !is.null(g(mr$mapping)),
    !is.null(g(mr$pathway)),
    NA
  )
  stats::setNames(done, MULTI_PIPELINE_STAGES)
}

multi_pipeline_progress_ui <- function(multi_dataset, multi_results) {
  pr <- multi_pipeline_progress(multi_dataset, multi_results)
  items <- lapply(seq_along(pr), function(i) {
    st <- pr[[i]]
    icn <- if (isTRUE(st)) icon("circle-check") else if (is.na(st)) icon("circle") else icon("circle-xmark")
    col <- if (isTRUE(st)) ARTHOMIX_COLORS$aqua else "var(--color-ink-muted, #898781)"
    tagList(span(style = sprintf("color:%s; white-space:nowrap;", col), icn, " ", names(pr)[i]), if (i < length(pr)) span(style = "color:var(--color-ink-muted, #898781);", " → "))
  })
  div(style = "font-size:0.82em; margin-top:6px; display:flex; flex-wrap:wrap; gap:2px 4px; align-items:center;",
      tags$strong("Pipeline: "), items)
}
