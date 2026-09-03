## Module 3 (Multiomics) - Biomarker Card sub-module, via testServer(): a
## fully synchronous, read-only interpretation layer over
## multi_results$concordance$df. Verifies the real evidence-tier

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "08_Biomarker_Card", "mod_multi_biomarkercard.R"))

mbc_conc_fixture <- function() {
  data.frame(
    gene_symbol = c("TP53", "BRCA1", "EGFR"),
    cpg = c("cg001", "cg002", "cg003"),
    sig_expression = c(TRUE, TRUE, FALSE),
    sig_methylation = c(TRUE, TRUE, FALSE),
    log2fc = c(2, 1, 0.1), dbeta = c(-0.3, -0.2, 0.01),
    correlation_r = c(-0.8, NA, NA), correlation_fdr = c(0.001, NA, NA),
    direction_classification = c("Up expression + Hypomethylation", "Up expression + Hypomethylation", "Weak/uncertain"),
    diablo = c(TRUE, FALSE, FALSE), snf = c(FALSE, TRUE, FALSE),
    priority_score = c(85, 40, 10),
    stringsAsFactors = FALSE
  )
}

test_that("evidence_df() classifies a real inverse-direction, negatively-correlated pair as 'Strong candidate' via the real cx_classify_evidence()", {
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload")
  multi_results <- shiny::reactiveValues(concordance = list(df = mbc_conc_fixture()))
  shiny::testServer(mod_multi_biomarkercard_server, args = list(id = "bc", multi_dataset = multi_dataset, multi_results = multi_results), {
    df <- evidence_df()
    expect_equal(df$evidence_tier[df$gene_symbol == "TP53"], "Strong candidate")
    expect_equal(df$evidence_tier[df$gene_symbol == "BRCA1"], "Moderate candidate")
    expect_equal(df$evidence_tier[df$gene_symbol == "EGFR"], "Insufficient evidence")
  })
})

test_that("selecting a row in the browse table sets mbc_selected(), and clicking Generate flips has_card() to TRUE", {
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload")
  multi_results <- shiny::reactiveValues(concordance = list(df = mbc_conc_fixture()))
  shiny::testServer(mod_multi_biomarkercard_server, args = list(id = "bc", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(mbc_browse_table_rows_selected = 1)
    sel <- mbc_selected()
    expect_equal(sel$gene_symbol, "TP53")
    expect_false(has_card())

    session$setInputs(mbc_generate_btn = 0)
    session$setInputs(mbc_generate_btn = 1)
    expect_true(has_card())
  })
})

test_that("text search finds a real case-insensitive match by gene symbol or CpG ID, and reports no match honestly otherwise", {
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload")
  multi_results <- shiny::reactiveValues(concordance = list(df = mbc_conc_fixture()))
  shiny::testServer(mod_multi_biomarkercard_server, args = list(id = "bc", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(mbc_search_mode = "type", mbc_text_input = "tp53")
    hits <- mbc_text_hits()
    expect_equal(nrow(hits), 1L)
    expect_equal(hits$gene_symbol, "TP53")

    session$setInputs(mbc_text_input = "cg002")
    expect_equal(mbc_text_hits()$gene_symbol, "BRCA1")

    session$setInputs(mbc_text_input = "NOT_A_REAL_GENE")
    expect_null(mbc_text_hits())
  })
})

test_that("the stale-selection guard clears mbc_selected()/has_card() when Concordance is re-run or the dataset source changes", {
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload")
  multi_results <- shiny::reactiveValues(concordance = list(df = mbc_conc_fixture()))
  shiny::testServer(mod_multi_biomarkercard_server, args = list(id = "bc", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(mbc_browse_table_rows_selected = 1)
    session$setInputs(mbc_generate_btn = 0)
    session$setInputs(mbc_generate_btn = 1)
    expect_true(has_card())

    nudge <- mbc_conc_fixture(); nudge$priority_score <- nudge$priority_score + 1
    multi_results$concordance <- list(df = nudge)
    session$flushReact()
    multi_results$concordance <- list(df = mbc_conc_fixture())
    session$flushReact()
    expect_null(mbc_selected())
    expect_false(has_card())
  })
})

test_that("evidence_df()/conc_df() return NULL honestly when Gene-CpG Concordance hasn't been run this session", {
  multi_dataset <- shiny::reactiveValues(active = FALSE, source = NULL)
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_biomarkercard_server, args = list(id = "bc", multi_dataset = multi_dataset, multi_results = multi_results), {
    expect_null(conc_df())
    expect_null(evidence_df())
  })
})
