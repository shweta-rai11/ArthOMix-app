## Module 1 (Transcriptomics) - Functional Enrichment: gene-panel source
## gating (live Feature Selection/Cross-Ancestry results take precedence
## over the bundled fallback, and the bundled fallback is only ever offered
## for the exact bundled reference dataset - a prior audit fix, regression-
## guarded here) plus a real, offline (no network - GO_BP uses the bundled
## org.Hs.eg.db annotation, not a live KEGG/Reactome API call) enrichGO()
## run verifying the scientific/data-contract: pathway identifiers,
## enrichment statistics, p/adjusted-p-values, non-empty valid output.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_enrichment.R"))

test_that("same_tissue_panel prefers this session's live Feature Selection consensus panel over any bundled fallback", {
  dataset <- shiny::reactiveValues(expr = load_default_dataset()$expr, meta = load_default_dataset()$meta,
                                     source = "test", source_type = "preloaded", is_bundled_reference = TRUE)
  results <- shiny::reactiveValues(featureselection = list(female = list(consensus_genes = c("TP53", "IL6"))))
  shiny::testServer(mod_enrichment_server, args = list(id = "enr", dataset = dataset, results = results), {
    p <- same_tissue_panel("female")
    expect_true(p$is_live)
    expect_setequal(p$genes, c("TP53", "IL6"))
  })
})

test_that("the bundled same-tissue/cross-tissue/cross-ancestry fallback panels are NEVER offered for a non-bundled (uploaded/GEO) dataset", {
  ## Regression guard: this fallback used to apply regardless of dataset
  ## provenance; it must only ever apply to the exact bundled reference
  ## cohort (dataset$is_bundled_reference), matching the same isolation
  ## fix already made in mod_wgcna.R/mod_diagnostic.R/mod_featureselection.R.
  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = "uploaded test",
                                     source_type = "uploaded", is_bundled_reference = FALSE)
  results <- shiny::reactiveValues()  ## no live results this session
  shiny::testServer(mod_enrichment_server, args = list(id = "enr", dataset = dataset, results = results), {
    st <- same_tissue_panel("female")
    ct <- cross_tissue_panel("female")
    ca <- cross_ancestry_panel("female")
    expect_length(st$genes, 0)
    expect_length(ct$genes, 0)
    expect_length(ca$genes, 0)
    expect_true(grepl("No live female panel yet", st$note))
  })
})

## NOTE (not written as a test): load_panel's observeEvent() calls
## updateTextAreaInput() on success, which is unobservable through
## testServer (shiny:::MockShinySession$sendInputMessage() is an explicit
## no-op), AND its validate(need(...)) failure for an empty panel has no
## externally observable effect either - unlike a reactive()/eventReactive()
## whose thrown condition propagates to an explicit caller, a validate()
## failure inside a bare observeEvent is contained by Shiny's normal per-
## observer error isolation and never reaches session$setInputs()'s caller.
## The real logic this observer depends on - which panel resolves to which
## genes, and the bundled/live/empty gating - is already fully covered by
## the two tests above (same_tissue_panel() etc., called directly).

## ---- Real, offline enrichGO() scientific-contract check --------------------

test_that("enrichGO (GO_BP, offline via bundled org.Hs.eg.db) on a real immune gene panel returns a valid, non-empty enrichment table", {
  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = "test",
                                     source_type = "preloaded", is_bundled_reference = TRUE)
  results <- shiny::reactiveValues()
  ## A well-known, strongly immune-associated gene set - real GO BP terms
  ## (e.g. cytokine/inflammatory response) are near-certain to come back.
  immune_genes <- c("IL6", "TNF", "IL1B", "CXCL8", "IL10", "STAT3", "NFKB1", "TLR4")

  shiny::testServer(mod_enrichment_server, args = list(id = "enr", dataset = dataset, results = results), {
    session$setInputs(gene_source = "own", gene_list = paste(immune_genes, collapse = "\n"), ontology = "GO_BP")
    ## result()'s eventReactive is ignoreInit=TRUE on a bare input$run_btn -
    ## needs priming (see feedback_shiny_testserver_ignoreinit_actionbutton_priming).
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)

    res <- result()
    df <- res$table
    expect_true(all(c("ID", "Description", "pvalue", "p.adjust", "qvalue", "geneID") %in% colnames(df)))
    expect_gt(nrow(df), 0)
    expect_true(all(df$pvalue >= 0 & df$pvalue <= 1))
    expect_true(all(df$p.adjust >= df$pvalue - 1e-9))  ## BH-adjusted p is never smaller than the raw p
  })
})
