## Module 1 (Transcriptomics) - Functional Enrichment: gene-panel source
## gating (live Feature Selection/Cross-Ancestry results take precedence
## over the bundled fallback, and the bundled fallback is only ever offered

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "transcriptomics", "14_Functional_Enrichment", "mod_enrichment.R"))

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
  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = "uploaded test",
                                     source_type = "uploaded", is_bundled_reference = FALSE)
  results <- shiny::reactiveValues()
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

test_that("enrichGO (GO_BP, offline via bundled org.Hs.eg.db) on a real immune gene panel returns a valid, non-empty enrichment table", {
  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = "test",
                                     source_type = "preloaded", is_bundled_reference = TRUE)
  results <- shiny::reactiveValues()
  immune_genes <- c("IL6", "TNF", "IL1B", "CXCL8", "IL10", "STAT3", "NFKB1", "TLR4")

  shiny::testServer(mod_enrichment_server, args = list(id = "enr", dataset = dataset, results = results), {
    session$setInputs(gene_source = "own", gene_list = paste(immune_genes, collapse = "\n"), ontology = "GO_BP")
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)

    res <- result()
    df <- res$table
    expect_true(all(c("ID", "Description", "pvalue", "p.adjust", "qvalue", "geneID") %in% colnames(df)))
    expect_gt(nrow(df), 0)
    expect_true(all(df$pvalue >= 0 & df$pvalue <= 1))
    expect_true(all(df$p.adjust >= df$pvalue - 1e-9))
  })
})
