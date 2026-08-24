## R/submodules_registry.R
## Assembles TX_MODULES from every mod_<id>.R file's config/ui/server trio.
## Sourced after all of them: Shiny sources R/*.R alphabetically
## (shiny:::loadSupport, sort_c()), and "submodules_registry.R" sorts after
## every "mod_*.R" filename, so every mod_<id>_config/_ui/_server referenced
## below already exists by the time this file runs.

TX_MODULES <- list(
  list(config = mod_overview_config,        ui = mod_overview_ui,        server = mod_overview_server),
  list(config = mod_preprocessing_config,   ui = mod_preprocessing_ui,   server = mod_preprocessing_server),
  list(config = mod_dge_config,             ui = mod_dge_ui,             server = mod_dge_server),
  list(config = mod_wgcna_config,           ui = mod_wgcna_ui,           server = mod_wgcna_server),
  list(config = mod_candidates_config,      ui = mod_candidates_ui,      server = mod_candidates_server),
  list(config = mod_mr_config,              ui = mod_mr_ui,              server = mod_mr_server),
  list(config = mod_coloc_config,           ui = mod_coloc_ui,           server = mod_coloc_server),
  list(config = mod_featureselection_config, ui = mod_featureselection_ui, server = mod_featureselection_server),
  list(config = mod_diagnostic_config,      ui = mod_diagnostic_ui,      server = mod_diagnostic_server),
  list(config = mod_interaction_config,     ui = mod_interaction_ui,     server = mod_interaction_server),
  list(config = mod_crosstissue_config,     ui = mod_crosstissue_ui,     server = mod_crosstissue_server),
  list(config = mod_crossancestry_config,   ui = mod_crossancestry_ui,   server = mod_crossancestry_server),
  list(config = mod_enrichment_config,      ui = mod_enrichment_ui,      server = mod_enrichment_server),
  list(config = mod_deconvolution_config,   ui = mod_deconvolution_ui,   server = mod_deconvolution_server),
  list(config = mod_nomogram_config,        ui = mod_nomogram_ui,        server = mod_nomogram_server),
  list(config = mod_biomarkercard_config,   ui = mod_biomarkercard_ui,   server = mod_biomarkercard_server)
)

TX_MODULES_BY_ID <- setNames(TX_MODULES, vapply(TX_MODULES, function(m) m$config$id, character(1)))

## Methylomics sub-modules - same config/ui/server trio shape as TX_MODULES
## above, added to as the Methylomics module grows past Dataset + Quality
## Control (see ui.R's methylomicsUI(), which reuses build_submodule_grid()
## against this list instead of TX_MODULES).
##
## Order follows the sex-stratified methylomics pipeline
## (Research_Q3_METHYLOMICS_sexstratified_COPY/methylomics/script0N_*):
## qc/normalization/dmp are built; celltype through diagnostic are
## registry-only placeholder scaffolds (mod_methyl_<id>_ui() renders a
## "not built yet" box) queued up to be built out one at a time.
MX_MODULES <- list(
  list(config = mod_methyl_qc_config,             ui = mod_methyl_qc_ui,             server = mod_methyl_qc_server),
  list(config = mod_methyl_normalization_config,  ui = mod_methyl_normalization_ui,  server = mod_methyl_normalization_server),
  list(config = mod_methyl_celltype_config,        ui = mod_methyl_celltype_ui,        server = mod_methyl_celltype_server),
  list(config = mod_methyl_dmp_config,             ui = mod_methyl_dmp_ui,             server = mod_methyl_dmp_server),
  list(config = mod_methyl_dmr_config,             ui = mod_methyl_dmr_ui,             server = mod_methyl_dmr_server),
  list(config = mod_methyl_wgcna_config,           ui = mod_methyl_wgcna_ui,           server = mod_methyl_wgcna_server),
  list(config = mod_methyl_candidates_config,      ui = mod_methyl_candidates_ui,      server = mod_methyl_candidates_server),
  list(config = mod_methyl_featureselection_config, ui = mod_methyl_featureselection_ui, server = mod_methyl_featureselection_server),
  list(config = mod_methyl_mr_config,              ui = mod_methyl_mr_ui,              server = mod_methyl_mr_server),
  list(config = mod_methyl_coloc_config,           ui = mod_methyl_coloc_ui,           server = mod_methyl_coloc_server),
  list(config = mod_methyl_diagnostic_config,      ui = mod_methyl_diagnostic_ui,      server = mod_methyl_diagnostic_server),
  list(config = mod_methyl_biomarkercard_config,   ui = mod_methyl_biomarkercard_ui,   server = mod_methyl_biomarkercard_server)
)

MX_MODULES_BY_ID <- setNames(MX_MODULES, vapply(MX_MODULES, function(m) m$config$id, character(1)))

## Cross-Omics sub-modules - same config/ui/server trio shape as TX_MODULES/
## MX_MODULES above (see ui.R's crossomicsUI(), which reuses
## build_submodule_grid() against this list instead of TX_MODULES/MX_MODULES).
##
## Order follows the cross-omics pipeline's own two scripts
## (Research_Q4_cross_Omics_sexstratified_COPY/cross_Omics_Sexstratified_COPY/
## scripts/01_*, 02_*): both are registry-only placeholder scaffolds
## (mod_cross_<id>_ui() renders a "not built yet" box) queued up to be built
## out one at a time, same convention MX_MODULES used while Methylomics was
## still partway built.
CX_MODULES <- list(
  list(config = mod_cross_integration_config,    ui = mod_cross_integration_ui,    server = mod_cross_integration_server),
  list(config = mod_cross_biomarker_conv_config, ui = mod_cross_biomarker_conv_ui, server = mod_cross_biomarker_conv_server),
  list(config = mod_cross_mr_stage_config,       ui = mod_cross_mr_stage_ui,       server = mod_cross_mr_stage_server)
)

CX_MODULES_BY_ID <- setNames(CX_MODULES, vapply(CX_MODULES, function(m) m$config$id, character(1)))

## Multi-Omics sub-modules - same config/ui/server trio shape as TX_MODULES/
## MX_MODULES/CX_MODULES above (see ui.R's multiomicsUI(), which reuses
## build_submodule_grid() against this list). Order follows
## Research_05_multiomics_sexstratified's own pipeline stages: cohort/sample
## harmonization first, then DIABLO+SNF integration, patient stratification,
## joint biomarker discovery, gene<->CpG concordance, pathway enrichment,
## and finally the session summary/reproducibility rollup.
MULTI_MODULES <- list(
  list(config = mod_multi_overview_config,       ui = mod_multi_overview_ui,       server = mod_multi_overview_server),
  list(config = mod_multi_integration_config,    ui = mod_multi_integration_ui,    server = mod_multi_integration_server),
  list(config = mod_multi_stratification_config, ui = mod_multi_stratification_ui, server = mod_multi_stratification_server),
  list(config = mod_multi_biomarker_config,      ui = mod_multi_biomarker_ui,      server = mod_multi_biomarker_server),
  list(config = mod_multi_concordance_config,    ui = mod_multi_concordance_ui,    server = mod_multi_concordance_server),
  list(config = mod_multi_pathway_config,        ui = mod_multi_pathway_ui,        server = mod_multi_pathway_server),
  list(config = mod_multi_live_config,           ui = mod_multi_live_ui,           server = mod_multi_live_server),
  list(config = mod_multi_summary_config,        ui = mod_multi_summary_ui,        server = mod_multi_summary_server)
)

MULTI_MODULES_BY_ID <- setNames(MULTI_MODULES, vapply(MULTI_MODULES, function(m) m$config$id, character(1)))

## ---------------------------------------------------------------------------
## AI assistant context
## ---------------------------------------------------------------------------
## Turns the shared `dataset` reactiveValues + `results` reactiveValues into a
## flat text block for ArthOChat's system prompt. Modules that haven't
## published a results$<id> yet (never run this session) render as
## "not yet run" rather than being silently omitted, so the assistant can
## tell the user what to go compute instead of guessing.
##
## Lives here, not in global.R, because it reads TX_MODULES/TX_MODULES_BY_ID
## above: shiny:::loadSupport() sources global.R into .GlobalEnv but every
## R/*.R file (this one included) into a separate child environment, so a
## global.R-defined function can never see a binding made in R/*.R - only
## the other way around. Defining it here instead keeps it in the same
## environment as the TX_MODULES it depends on.

## `cross_results` is optional (defaults to NULL, so every existing call site
## keeps producing byte-identical output) - only mod_arthochat_server()'s
## Cross-Omics-aware call site passes it, appending a section grounded in
## cross_results$integration (set by mod_cross_integration_server() once an
## "Expression x Methylation" integration has actually been run).
build_assistant_context <- function(dataset, results, cross_results = NULL) {
  meta <- dataset$meta
  grp_tbl <- table(meta$group)
  lines <- c(
    "## Currently loaded dataset",
    sprintf("- Source: %s", dataset$source),
    sprintf("- %s genes x %s samples", format(nrow(dataset$expr), big.mark = ","), ncol(dataset$expr)),
    sprintf("- Groups: %s", paste(sprintf("%s (n=%d)", names(grp_tbl), grp_tbl), collapse = ", ")),
    if ("sex" %in% names(meta) && any(!is.na(meta$sex))) {
      sex_tbl <- table(meta$sex)
      sprintf("- Sex: %s", paste(sprintf("%s (n=%d)", names(sex_tbl), sex_tbl), collapse = ", "))
    },
    "",
    "## Computed analysis results (this session)"
  )

  module_ids <- vapply(TX_MODULES, function(m) m$config$id, character(1))
  for (mid in module_ids) {
    title <- TX_MODULES_BY_ID[[mid]]$config$title
    res <- results[[mid]]
    if (is.null(res)) {
      lines <- c(lines, sprintf("### %s", title), "(not yet run in this session)")
    } else {
      kv <- vapply(names(res), function(nm) {
        v <- res[[nm]]
        sprintf("- %s: %s", nm, paste(utils::head(as.character(v), 20), collapse = ", "))
      }, character(1))
      lines <- c(lines, sprintf("### %s", title), kv)
    }
  }

  cx <- cross_results$integration
  lines <- c(lines, "", "## Cross-Omics: Expression x Methylation integration")
  if (is.null(cx)) {
    lines <- c(lines, "(not yet run in this session)")
  } else {
    s <- cx$summary
    counts <- s$counts
    lines <- c(
      lines,
      sprintf("- Dataset: %s", cx$params$input_mode %||% "(unknown)"),
      sprintf("- Genes analyzed: %s", format(s$n_genes, big.mark = ",")),
      sprintf("- Significant DEGs: %s, Significant DMGs: %s, Significant in both: %s",
              format(s$n_deg, big.mark = ","), format(s$n_dmg, big.mark = ","), format(s$n_integrated, big.mark = ",")),
      sprintf("- Hyper + Down (potential methylation-associated repression): %s", counts[["Hyper + Down"]] %||% 0),
      sprintf("- Hypo + Up (potential methylation-associated activation): %s", counts[["Hypo + Up"]] %||% 0),
      sprintf("- Hyper + Up (concordant-direction/noncanonical): %s", counts[["Hyper + Up"]] %||% 0),
      sprintf("- Hypo + Down (concordant-direction/noncanonical): %s", counts[["Hypo + Down"]] %||% 0),
      sprintf("- Sample matching: %s", cx$params$sample_matching %||% "Not available"),
      "- These are statistical associations, not established causal relationships - never state that methylation \"causes\" an expression change."
    )
  }
  paste(lines, collapse = "\n")
}
