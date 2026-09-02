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
  list(config = mod_methyl_validation_config,      ui = mod_methyl_validation_ui,      server = mod_methyl_validation_server),
  list(config = mod_methyl_biomarkercard_config,   ui = mod_methyl_biomarkercard_ui,   server = mod_methyl_biomarkercard_server)
)

MX_MODULES_BY_ID <- setNames(MX_MODULES, vapply(MX_MODULES, function(m) m$config$id, character(1)))

## Cross-Omics sub-modules - same config/ui/server trio shape as TX_MODULES/
## MX_MODULES above (see ui.R's crossomicsUI(), which reuses
## build_submodule_grid() against this list instead of TX_MODULES/MX_MODULES).
##
## All three are fully built (see mod_cross_integration.R, mod_cross_biomarker_conv.R,
## mod_cross_mr_stage.R for their real UI/server implementations) - this list
## no longer holds registry-only placeholder scaffolds, unlike when it followed
## the same convention MX_MODULES used while Methylomics was still partway built.
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
## biomarker discovery, gene<->CpG concordance, and pathway enrichment.
## The MOFA2 module (mod_multi_mofa.R / mod_multi_mofa_engine.R)
## is no longer a standalone sub-module here; it's mounted directly inside
## the Dataset Workspace tab (mod_multi_dataset.R) as "Integrated Analysis
## (MOFA2)", since it operates on the same Active Multi-Omics Dataset that
## tab builds.
MULTI_MODULES <- list(
  list(config = mod_multi_overview_config,       ui = mod_multi_overview_ui,       server = mod_multi_overview_server),
  list(config = mod_multi_integration_config,    ui = mod_multi_integration_ui,    server = mod_multi_integration_server),
  list(config = mod_multi_stratification_config, ui = mod_multi_stratification_ui, server = mod_multi_stratification_server),
  list(config = mod_multi_biomarker_config,      ui = mod_multi_biomarker_ui,      server = mod_multi_biomarker_server),
  list(config = mod_multi_concordance_config,    ui = mod_multi_concordance_ui,    server = mod_multi_concordance_server),
  list(config = mod_multi_pathway_config,        ui = mod_multi_pathway_ui,        server = mod_multi_pathway_server),
  list(config = mod_multi_biomarkercard_config,  ui = mod_multi_biomarkercard_ui,  server = mod_multi_biomarkercard_server),
  list(config = mod_multi_summary_config,        ui = mod_multi_summary_ui,        server = mod_multi_summary_server)
)

MULTI_MODULES_BY_ID <- setNames(MULTI_MODULES, vapply(MULTI_MODULES, function(m) m$config$id, character(1)))

## ---------------------------------------------------------------------------
## AI assistant context - module-scoped
## ---------------------------------------------------------------------------
## Turns each vertical's shared dataset/results reactiveValues into a flat
## text block for ArthOChat's system prompt. Modules that haven't published a
## results$<id> yet (never run this session) render as "not yet run" rather
## than being silently omitted, so the assistant can tell the user what to go
## compute instead of guessing.
##
## Every build_*_context() below takes an optional `focus_id`: a single
## TX_MODULES/MX_MODULES/CX_MODULES/MULTI_MODULES config$id narrows the loop
## to just that one sub-module - this is what lets ArthOChat's default
## per-message context be scoped to whatever sub-module the user is actually
## looking at (see mod_arthochat.R's `current_context`), instead of always
## dumping every sub-module of every vertical regardless of what's on screen.
## `focus_id = NULL` (the fallback used for the Dataset/Sub-modules-picker
## tabs, the Home/Modules landing pages, and the other_module_context() tool
## an explicit cross-module question can call) dumps every sub-module of that
## vertical, same as this file's original single build_assistant_context().
##
## Lives here, not in global.R, because it reads TX_MODULES/TX_MODULES_BY_ID
## etc above: shiny:::loadSupport() sources global.R into .GlobalEnv but every
## R/*.R file (this one included) into a separate child environment, so a
## global.R-defined function can never see a binding made in R/*.R - only
## the other way around. Defining it here instead keeps it in the same
## environment as the *_MODULES lists it depends on.

## Renders one sub-module's results block - shared by every build_*_context()
## below so "not yet run" / key-value formatting stays identical everywhere.
.format_results_block <- function(title, res) {
  if (is.null(res)) {
    return(c(
      sprintf("### %s", title),
      "NOT YET RUN IN THIS SESSION - no numbers exist for this sub-module yet.",
      "If asked for a result from it, say it hasn't been run in this session -",
      "do not answer with a number from methodology/literature instead."
    ))
  }
  kv <- vapply(names(res), function(nm) {
    v <- res[[nm]]
    sprintf("- %s: %s", nm, paste(utils::head(as.character(v), 20), collapse = ", "))
  }, character(1))
  c(sprintf("### %s", title), kv)
}

build_tx_context <- function(dataset, results, focus_id = NULL) {
  meta <- dataset$meta
  grp_tbl <- table(meta$group)
  lines <- c(
    "## Transcriptomics: currently loaded dataset",
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
  ids <- focus_id %||% vapply(TX_MODULES, function(m) m$config$id, character(1))
  for (mid in ids) {
    lines <- c(lines, .format_results_block(TX_MODULES_BY_ID[[mid]]$config$title, results[[mid]]))
  }
  paste(lines, collapse = "\n")
}

build_mx_context <- function(methyl_dataset, methyl_results, focus_id = NULL) {
  if (is.null(methyl_dataset$source)) {
    return("## Methylomics: no dataset loaded yet in this session (visit the Methylomics Dataset tab).")
  }
  lines <- c(
    "## Methylomics: currently loaded dataset",
    sprintf("- Source: %s", methyl_dataset$source),
    sprintf("- Array type: %s", methyl_dataset$array_type %||% "(unknown)"),
    if (!is.null(methyl_dataset$beta)) {
      sprintf("- %s probes x %s samples", format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta))
    },
    "",
    "## Computed analysis results (this session)"
  )
  ids <- focus_id %||% vapply(MX_MODULES, function(m) m$config$id, character(1))
  for (mid in ids) {
    lines <- c(lines, .format_results_block(MX_MODULES_BY_ID[[mid]]$config$title, methyl_results[[mid]]))
  }
  paste(lines, collapse = "\n")
}

## The "integration" sub-module gets the bespoke summary-stats formatting it
## always had (cx$summary's counts aren't a flat named vector, so the generic
## .format_results_block() key-value loop can't render it usefully); every
## other Cross-Omics sub-module falls back to that generic loop.
.format_cx_integration <- function(cx) {
  if (is.null(cx)) return("(not yet run in this session)")
  s <- cx$summary
  counts <- s$counts
  c(
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

build_cx_context <- function(cross_results, focus_id = NULL) {
  ids <- focus_id %||% vapply(CX_MODULES, function(m) m$config$id, character(1))
  lines <- character(0)
  for (mid in ids) {
    title <- CX_MODULES_BY_ID[[mid]]$config$title
    if (identical(mid, "integration")) {
      lines <- c(lines, sprintf("## Cross-Omics: %s", title), .format_cx_integration(cross_results$integration))
    } else {
      lines <- c(lines, .format_results_block(paste("Cross-Omics:", title), cross_results[[mid]]))
    }
  }
  paste(lines, collapse = "\n")
}

build_mo_context <- function(multi_dataset, multi_results, focus_id = NULL) {
  if (!isTRUE(multi_dataset$active)) {
    return("## Multi-Omics: no active dataset loaded yet in this session (visit the Multi-Omics Dataset Workspace tab).")
  }
  lines <- c(
    "## Multi-Omics: currently loaded dataset",
    sprintf("- Source: %s", multi_dataset$source %||% "(unknown)"),
    sprintf("- Layers: %s", if (length(multi_dataset$layers)) paste(names(multi_dataset$layers), collapse = ", ") else "(none)"),
    sprintf("- Samples: %s", if (!is.null(multi_dataset$sample_meta)) nrow(multi_dataset$sample_meta) else "(unknown)"),
    "",
    "## Computed analysis results (this session)"
  )
  ids <- focus_id %||% vapply(MULTI_MODULES, function(m) m$config$id, character(1))
  for (mid in ids) {
    lines <- c(lines, .format_results_block(MULTI_MODULES_BY_ID[[mid]]$config$title, multi_results[[mid]]))
  }
  paste(lines, collapse = "\n")
}

## Whole-app fallback context - every Transcriptomics sub-module plus the
## Cross-Omics integration summary. Used when no single omics module is in
## view (Home/Modules landing pages) as ArthOChat's context builder.
build_assistant_context <- function(dataset, results, cross_results = NULL) {
  paste(build_tx_context(dataset, results), "", build_cx_context(cross_results %||% list(), focus_id = "integration"), sep = "\n")
}

## Dispatches to the one build_*_context() matching whichever top-level
## module (server.R's `input$sidebar_tabs` value) is currently open, scoped
## to `submodule_id` (that vertical's own tx_menu/mx_menu/cx_menu/mo_menu
## selection, resolved back to a config$id - NULL when the user is on that
## vertical's Dataset/Sub-modules-picker tab rather than inside one specific
## sub-module). Anything outside the four omics verticals (Home, Modules,
## the ArthOChat drawer's own tab id) falls back to the whole-app view.
build_scoped_assistant_context <- function(module, submodule_id,
                                            dataset, results,
                                            methyl_dataset, methyl_results,
                                            cross_dataset, cross_results,
                                            multi_dataset, multi_results) {
  switch(module,
    transcriptomics = build_tx_context(dataset, results, focus_id = submodule_id),
    methylomics = build_mx_context(methyl_dataset, methyl_results, focus_id = submodule_id),
    crossomics = build_cx_context(cross_results, focus_id = submodule_id),
    multiomics = build_mo_context(multi_dataset, multi_results, focus_id = submodule_id),
    build_assistant_context(dataset, results, cross_results)
  )
}
