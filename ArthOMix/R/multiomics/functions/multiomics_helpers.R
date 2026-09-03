## R/multiomics/functions/multiomics_helpers.R
## Shared low-level utilities for the Multi-Omics module: cell/cohort
## metadata (MULTI_CELLS), generic CSV table loaders, sample-harmonization

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
  colnames(sample_meta)[grepl("^sex$|^gender$", colnames(sample_meta), ignore.case = TRUE)]
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

multi_concordance_add_fdr <- function(df) {
  if (is.null(df)) return(df)
  out <- df
  if ("expr_p" %in% colnames(out)) out$expr_fdr <- stats::p.adjust(out$expr_p, method = "BH")
  if ("meth_p" %in% colnames(out)) out$meth_fdr <- stats::p.adjust(out$meth_p, method = "BH")
  out
}

multi_active_dataset_banner <- function(multi_dataset) {
  md <- multi_dataset %||% list()
  source <- md$source
  if (is.null(source) || !isTRUE(md$active %||% FALSE)) {
    return(div(class = "empty-note", icon("circle-info"),
               "No active Multi-Omics dataset yet - pick one on the Dataset tab, or choose the Preloaded/Reference option in this tab's own data-source selector to use the bundled RA anti-TNF cohort."))
  }
  if (identical(source, "preloaded")) {
    return(div(class = "empty-note", icon("circle-check"),
               tags$strong("Data source: Preloaded Dataset."), " Existing results are available and shown below."))
  }
  div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
      tags$strong(sprintf("Data source: %s.", if (identical(source, "geo")) "NCBI GEO" else "User Upload")),
      " No stored results for this dataset yet - select \"Active Multi-Omics Dataset\" above to run this module's analysis on it directly.")
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

MULTI_KNOWN_LIMITATIONS <- c(
  "No cell-type composition sub-module: methylation (EpiDISH) and expression (CIBERSORTx) deconvolution, composition-vs-disease testing, and composition-adjusted matrices are not implemented.",
  "Probe QC covers duplicate IDs, zero-variance and missingness filtering, and variance/MAD-ranked feature filtering - it does not include cross-reactive probe masking, SNP-overlapping probe filtering, sex-chromosome probe handling, or an expression/methylation sex check (XIST/RPS4Y1 vs. reported sex).",
  "No dedicated Covariates & Clinical Metadata sub-module: there is no clinical-variable coverage table, treatment-strata harmonization, sample x metadata missingness map, or covariate collinearity/VIF check.",
  "Validation & Stability covers leakage-safe nested cross-validation only - bootstrap feature-selection stability, a permutation null, external cohort scoring (e.g. GSE17755, GSE15573), calibration curves, and decision-curve analysis are not implemented.",
  "No unified Model Benchmarking panel comparing clinical-only, transcriptomics-only, methylomics-only, and combined models under identical cross-validation folds.",
  "The Machine Learning panel offers Elastic Net and Random Forest only - XGBoost, Boruta, SHAP explanations, and learning curves are not implemented.",
  "Sex differences are assessed by running the same analysis separately per sex stratum, not by a formal sex x disease interaction model, FDR-ranked interaction table, or sex-shared/sex-specific/interaction classification.",
  "No sex-chromosome biology annotation (X/Y/autosomal feature classification, X-inactivation escapees, hormone-responsive gene sets).",
  "Gene<->CpG mapping uses each CpG's Illumina manifest gene annotation, not a configurable genomic-distance window - and mediation analysis (methylation -> expression -> outcome) is not implemented.",
  "No Druggable Target Linker (Open Targets tractability, known drugs, or clinical-phase lookups for candidate genes).",
  "DIABLO exposes one cross-block design-weight control, not a design-matrix sweep grid, and does not render correlation-circle or arrow plots.",
  "No parameter-manifest JSON export or decision log - the session bundle below includes the loaded result tables, a plain-text report, and installed package versions only.",
  "Reported AUROCs (precomputed cohort tabs) are the source pipeline's own nested cross-validation performance, not performance on an independent replication cohort - see AUDIT.md."
)

MULTI_REPRODUCIBILITY_SCRIPTS <- c(
  "Shared upstream (sample matching, QC, normalization, annotation, genome-wide discovery, leakage-safe nested-CV benchmark): analyses/data_preparation/scripts/01-06_*.R",
  "SNF integration (classifier + unsupervised joint-biomarker discovery): analyses/01_female_male_adalimumab/scripts/07*.R, analyses/02_female_male_etanercept/scripts/07*.R",
  "DIABLO integration (drug x sex, response, drug-type): analyses/01_female_male_adalimumab/scripts/08_*.R, analyses/07_cross_analysis_summary/scripts/11_*.R, 14_*.R",
  "Gene<->CpG concordance: analyses/07_cross_analysis_summary/scripts/15_*.R, 18_*.R",
  "Pathway enrichment: analyses/07_cross_analysis_summary/scripts/16_*.R",
  "Cross-cell summary assembly: analyses/07_cross_analysis_summary/scripts/09_*.R, 12_*.R",
  "Independent audit (methodology, leakage findings/fixes, honest AUROC verdicts per cell): AUDIT.md",
  "Live-computation code (this app, not the source pipeline): R/multiomics/01_Data_Workspace/multiomics_dataset_helpers.R, R/multiomics/functions/multiomics_integration_helpers.R, R/multiomics/02_Cohort_Harmonization/cohort_harmonization_helpers.R, R/multiomics/06_Gene_CpG_Concordance/multiomics_concordance_helpers.R, R/multiomics/functions/multiomics_sexstratified_engine.R"
)

multi_build_report <- function(multi_results) {
  ids <- c("overview", "integration", "integration_stratified", "stratification", "biomarker", "concordance", "pathway", "live_qc", "live_mofa")
  lines <- c(
    "# ArthOMix Multi-Omics module - session report",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "This report lists which sub-modules had results loaded this session. Precomputed-cohort tabs browse Research_05_multiomics_sexstratified's own saved DIABLO/SNF/concordance/pathway output; Dataset Workspace/Live Analysis results are computed in this app from the currently active dataset.",
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

