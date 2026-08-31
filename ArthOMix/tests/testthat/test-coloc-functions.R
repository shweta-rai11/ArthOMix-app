## Module 1 (Transcriptomics) - Colocalization: runs coloc.abf on the real
## bundled eQTL cis-window instrument against the bundled RA GWAS (the
## "project" data source), via testServer(), verifying the scientific
## contract (5 posterior-probability hypotheses summing to ~1, valid ranges,
## SNP-level table structure) and the <10-shared-SNP validation gate.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_coloc.R"))

test_that("running colocalisation on a real bundled gene produces 5 posterior probabilities summing to ~1", {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  gene <- sort(names(coloc_regions))[1]

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(data_source = "project", gene = gene, case_frac = 0.33)
    session$setInputs(run_btn = 1)

    res <- coloc_result()
    pp <- as.numeric(res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
    expect_length(pp, 5)
    expect_true(all(pp >= 0 & pp <= 1))
    expect_equal(sum(pp), 1, tolerance = 1e-6)

    expect_true(all(c("snp", "eqtl_beta", "gwas_beta", "snp_pp_h4") %in% colnames(res$snp_df)))
    expect_gte(res$n_snp, 10)

    expect_false(is.null(results$coloc$genes_tested[[gene]]))
    expect_equal(results$coloc$genes_tested[[gene]]$pp_h4, round(pp[5], 3))
  })
})

test_that("switching genes and re-running updates results$coloc$genes_tested for the new gene without dropping the previous one", {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  genes <- sort(names(coloc_regions))
  skip_if(length(genes) < 2, "Needs at least 2 bundled coloc regions to test multi-gene accumulation.")

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(data_source = "project", gene = genes[1], case_frac = 0.33)
    session$setInputs(run_btn = 1)
    session$setInputs(gene = genes[2])
    session$setInputs(run_btn = 2)

    expect_true(all(genes[1:2] %in% names(results$coloc$genes_tested)))
  })
})
