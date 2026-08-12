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
  list(config = mod_nomogram_config,        ui = mod_nomogram_ui,        server = mod_nomogram_server)
)

TX_MODULES_BY_ID <- setNames(TX_MODULES, vapply(TX_MODULES, function(m) m$config$id, character(1)))

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

build_assistant_context <- function(dataset, results) {
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
  paste(lines, collapse = "\n")
}
