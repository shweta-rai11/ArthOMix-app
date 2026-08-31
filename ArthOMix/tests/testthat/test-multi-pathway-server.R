## Module 3 (Multiomics) - Pathways sub-module, via testServer(): the
## upload -> detect -> confirm-mapping -> run pipeline is SYNCHRONOUS (no
## ExtendedTask), so this drives a REAL, full click-through ORA run: a real
## CSV upload of immune genes, real column-role auto-detection, real
## identifier harmonization, and a real offline clusterProfiler::enrichGO()
## call - verifying a genuine immune-related GO term surfaces, not a
## fabricated result.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_concordance_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_pathway_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_pathway_plots.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_pathway.R"))

test_that("clicking 'Run' on a real uploaded immune-gene table (ORA, GO_BP, entire-database background) surfaces a real immune-related GO term", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  ## Same real, well-known immune/inflammation gene panel used in
  ## test-multi-pathway-helpers-functions.R's direct mp_run_ora_go() check -
  ## here it goes through the full upload/detect/confirm/run pipeline instead.
  genes <- c("IL6", "TNF", "CXCL8", "IL13", "CCL2", "IL15")
  df <- data.frame(gene_symbol = genes, log2FC = c(2.1, 1.8, -1.5, 1.2, 2.4, 1.1), pvalue = c(0.001, 0.002, 0.01, 0.02, 0.001, 0.03))
  path <- tempfile(fileext = ".csv")
  write.csv(df, path, row.names = FALSE)

  multi_dataset <- shiny::reactiveValues(active = FALSE, layers = list())
  multi_results <- shiny::reactiveValues()

  shiny::testServer(mod_multi_pathway_server, args = list(id = "mp", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(data_source = "upload")
    session$setInputs(upload_file = fx_mkfile(path))
    det <- upload_detected()
    expect_equal(det$detected$id_col, "gene_symbol")

    session$setInputs(map_id_col = "gene_symbol", map_effect_col = "log2FC", map_pvalue_col = "pvalue",
                        map_fdr_col = "(none)", map_direction_col = "(none)", map_omics_col = "(none)", map_sex_col = "(none)")
    session$setInputs(confirm_mapping = 1)
    expect_true(isTRUE(upload_confirmed()$ok))

    session$setInputs(database = "GO_BP", method = "ORA", background = "entire_database",
                        pval_cut = 0.2, min_size = 3, max_size = 500, fdr_cut = 0.25)
    session$setInputs(run_btn = 1)

    r <- result()
    expect_true(isTRUE(r$ok))
    expect_true(nrow(r$table) > 0)
    expect_true(any(grepl("immune|inflamm|cytokine", r$table$Description, ignore.case = TRUE)))
  })
})

test_that("the Run handler reports a clear error (never a crash) when no database is selected", {
  genes <- c("IL6", "TNF", "CXCL8")
  df <- data.frame(gene_symbol = genes, log2FC = c(2.1, 1.8, -1.5), pvalue = c(0.001, 0.002, 0.01))
  path <- tempfile(fileext = ".csv")
  write.csv(df, path, row.names = FALSE)

  multi_dataset <- shiny::reactiveValues(active = FALSE, layers = list())
  multi_results <- shiny::reactiveValues()

  shiny::testServer(mod_multi_pathway_server, args = list(id = "mp", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(data_source = "upload")
    session$setInputs(upload_file = fx_mkfile(path))
    session$setInputs(map_id_col = "gene_symbol", map_effect_col = "log2FC", map_pvalue_col = "pvalue",
                        map_fdr_col = "(none)", map_direction_col = "(none)", map_omics_col = "(none)", map_sex_col = "(none)")
    session$setInputs(confirm_mapping = 1)
    session$setInputs(database = character(0), method = "ORA")
    session$setInputs(run_btn = 1)
    expect_null(result())  ## returns early with a notification, never sets a fabricated result
  })
})
