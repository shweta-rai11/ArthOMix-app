## Module 1 (Transcriptomics) - Colocalization: runs coloc.abf on the real
## bundled eQTL cis-window instrument against the bundled RA GWAS (the
## "project" data source), via testServer(), verifying the scientific
## contract (5 posterior-probability hypotheses summing to ~1, valid ranges,
## SNP-level table structure) and the <10-shared-SNP validation gate.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "08_Colocalization", "mod_coloc.R"))

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

## Builds a synthetic uploaded-GWAS file from the bundled eQTL region itself
## (guarantees rsID/allele overlap without depending on network access), for
## exercising the "upload your own GWAS" arm end to end.
make_synthetic_gwas_upload <- function(eqtl, beta_scale = 1.05, extra_dup_rows = 0) {
  df <- data.frame(
    SNP = eqtl$rsid, beta = eqtl$beta * beta_scale, se = eqtl$se, pval = eqtl$p,
    ea = eqtl$ea, oa = eqtl$nea, eaf = eqtl$eaf, N = eqtl$n
  )
  if (extra_dup_rows > 0) df <- rbind(df, df[seq_len(extra_dup_rows), , drop = FALSE])
  path <- tempfile(fileext = ".csv")
  write.csv(df, path, row.names = FALSE)
  path
}

test_that("uploading a case/control GWAS (default trait type) runs coloc.abf and returns 5 valid posterior probabilities", {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  gene <- sort(names(coloc_regions))[1]
  eqtl <- as.data.frame(coloc_regions[[gene]]$eqtl)
  gwas_path <- make_synthetic_gwas_upload(eqtl)

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload", gene = gene, case_frac = 0.33, gwas_type = "cc",
                       gwas_label = "Synthetic cc trait", gwas_file = list(datapath = gwas_path, name = "gwas.csv"))
    session$setInputs(gwas_snp = "SNP", gwas_beta = "beta", gwas_se = "se", gwas_pval = "pval",
                       gwas_ea = "ea", gwas_oa = "oa", gwas_eaf = "eaf", gwas_n = "N")
    session$setInputs(run_btn = 1)

    res <- coloc_result()
    expect_true(res$uploaded)
    expect_equal(res$gwas_type, "cc")
    pp <- as.numeric(res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
    expect_length(pp, 5)
    expect_true(all(pp >= 0 & pp <= 1))
    expect_equal(sum(pp), 1, tolerance = 1e-6)
  })
})

test_that("uploading a quantitative-trait GWAS is analysed as type='quant' (not silently forced to case/control)", {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  gene <- sort(names(coloc_regions))[1]
  eqtl <- as.data.frame(coloc_regions[[gene]]$eqtl)
  gwas_path <- make_synthetic_gwas_upload(eqtl)

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload", gene = gene, case_frac = 0.33, gwas_type = "quant",
                       gwas_label = "Synthetic quant trait", gwas_file = list(datapath = gwas_path, name = "gwas.csv"))
    session$setInputs(gwas_snp = "SNP", gwas_beta = "beta", gwas_se = "se", gwas_pval = "pval",
                       gwas_ea = "ea", gwas_oa = "oa", gwas_eaf = "eaf", gwas_n = "N")
    session$setInputs(run_btn = 1)

    res <- coloc_result()
    expect_equal(res$gwas_type, "quant")
    pp <- as.numeric(res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
    expect_equal(sum(pp), 1, tolerance = 1e-6)
  })
})

test_that("choosing a quantitative trait without mapping an effect-allele-frequency column fails with an actionable message, not a raw coloc.abf error", {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  gene <- sort(names(coloc_regions))[1]
  eqtl <- as.data.frame(coloc_regions[[gene]]$eqtl)
  df <- data.frame(SNP = eqtl$rsid, beta = eqtl$beta, se = eqtl$se, pval = eqtl$p, ea = eqtl$ea, oa = eqtl$nea, N = eqtl$n)
  gwas_path <- tempfile(fileext = ".csv")
  write.csv(df, gwas_path, row.names = FALSE)

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload", gene = gene, case_frac = 0.33, gwas_type = "quant",
                       gwas_label = "No EAF trait", gwas_file = list(datapath = gwas_path, name = "gwas.csv"))
    session$setInputs(gwas_snp = "SNP", gwas_beta = "beta", gwas_se = "se", gwas_pval = "pval",
                       gwas_ea = "ea", gwas_oa = "oa", gwas_eaf = "", gwas_n = "N")
    session$setInputs(run_btn = 1)

    err <- tryCatch(coloc_result(), error = function(e) e)
    expect_true(inherits(err, "shiny.silent.error"))
    expect_match(conditionMessage(err), "effect-allele-frequency", fixed = TRUE)
  })
})

test_that("duplicated SNP rows in the uploaded GWAS file are removed before coloc.abf, not double-counted", {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  gene <- sort(names(coloc_regions))[1]
  eqtl <- as.data.frame(coloc_regions[[gene]]$eqtl)
  gwas_path <- make_synthetic_gwas_upload(eqtl, extra_dup_rows = 5)

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload", gene = gene, case_frac = 0.33, gwas_type = "cc",
                       gwas_label = "Dup-row trait", gwas_file = list(datapath = gwas_path, name = "gwas.csv"))
    session$setInputs(gwas_snp = "SNP", gwas_beta = "beta", gwas_se = "se", gwas_pval = "pval",
                       gwas_ea = "ea", gwas_oa = "oa", gwas_eaf = "eaf", gwas_n = "N")
    session$setInputs(run_btn = 1)

    res <- coloc_result()
    expect_equal(res$n_snp, length(unique(res$snp_df$snp)))
    expect_false(any(duplicated(res$snp_df$snp)))
  })
})
