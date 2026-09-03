## Module 2 (Methylomics) - Colocalization server tests: exercises the Upload
## Data route (validate -> filters/priors -> run) via testServer(), on
## synthetic mQTL/GWAS summary statistics built from the real bundled eQTL
## rsIDs (so TwoSampleMR::format_data()/harmonise_data() see well-formed
## variant IDs and alleles). Focused on the p1/p2/p12 prior inputs: confirms
## they flow into coloc.abf(), are validated, and leave default behaviour
## (the pre-existing coloc conventions of 1e-4/1e-4/1e-5) unchanged.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "12_Colocalization", "mod_methyl_coloc.R"))

make_synthetic_meth_coloc_files <- function() {
  coloc_regions <- readRDS(COLOC_REGIONS_RDS)
  gene <- sort(names(coloc_regions))[1]
  eqtl <- as.data.frame(coloc_regions[[gene]]$eqtl)

  meth_df <- data.frame(
    cpg = "cg00000029", snp = eqtl$rsid, beta = eqtl$beta, se = eqtl$se, pval = eqtl$p,
    ea = eqtl$ea, oa = eqtl$nea, eaf = eqtl$eaf, n = eqtl$n,
    chr = eqtl$chr, pos = eqtl$position
  )
  ## Deliberately weakened, noisy GWAS signal (not a near-perfect echo of the
  ## eQTL effect): a slam-dunk concordant signal saturates PP.H4 near 1.0
  ## regardless of the p1/p2/p12 priors, which would make the prior-sensitivity
  ## test below vacuous. A noisier, only-moderately-supported signal keeps
  ## PP.H4 in a mid-range where changing the priors is actually detectable.
  set.seed(20260904)
  noise <- stats::rnorm(nrow(eqtl), mean = 0, sd = eqtl$se * 4)
  gwas_df <- data.frame(
    snp = eqtl$rsid, beta = eqtl$beta * 0.4 + noise, se = eqtl$se, pval = eqtl$p,
    ea = eqtl$ea, oa = eqtl$nea, n = eqtl$n
  )
  meth_path <- tempfile(fileext = ".csv"); write.csv(meth_df, meth_path, row.names = FALSE)
  gwas_path <- tempfile(fileext = ".csv"); write.csv(gwas_df, gwas_path, row.names = FALSE)
  list(meth_path = meth_path, gwas_path = gwas_path)
}

run_methyl_coloc_upload <- function(p1 = NULL, p2 = NULL, p12 = NULL, expect_run_error = FALSE) {
  files <- make_synthetic_meth_coloc_files()
  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  out <- NULL
  shiny::testServer(mod_methyl_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    session$setInputs(
      data_source = "upload",
      meth_file = list(datapath = files$meth_path, name = "meth.csv"),
      gwas_file = list(datapath = files$gwas_path, name = "gwas.csv"),
      gwas_label = "Synthetic GWAS", gwas_type = "cc", case_frac = 0.33
    )
    session$setInputs(
      meth_cpg = "cpg", meth_snp = "snp", meth_beta = "beta", meth_se = "se", meth_pval = "pval",
      meth_ea = "ea", meth_oa = "oa", meth_eaf = "eaf", meth_n = "n",
      meth_snp_chr = "chr", meth_snp_pos = "pos"
    )
    session$setInputs(
      gwas_snp = "snp", gwas_beta = "beta", gwas_se = "se", gwas_pval = "pval",
      gwas_ea = "ea", gwas_oa = "oa", gwas_n = "n"
    )
    session$setInputs(validate_btn = 1)

    prior_inputs <- list()
    if (!is.null(p1)) prior_inputs$p1 <- p1
    if (!is.null(p2)) prior_inputs$p2 <- p2
    if (!is.null(p12)) prior_inputs$p12 <- p12
    if (length(prior_inputs) > 0) do.call(session$setInputs, prior_inputs)

    if (isTRUE(expect_run_error)) {
      ## The run_btn observeEvent() wraps build_run_state_upload() in its own
      ## tryCatch() and only shows a notification on a validate() failure - it
      ## never re-throws, and run_state() is simply left unset. So to actually
      ## observe the validate() condition from a test, call the (module-local)
      ## builder function directly rather than going through run_btn/run_state.
      out <<- tryCatch(build_run_state_upload(), error = function(e) e)
    } else {
      session$setInputs(run_btn = 1)
      out <<- run_state()
    }
  })
  out
}

test_that("uploading synthetic mQTL/GWAS data runs coloc.abf with coloc's default priors when none are set", {
  rs <- run_methyl_coloc_upload()
  expect_equal(rs$mode, "upload")
  expect_equal(rs$priors, list(p1 = MCOL_DEFAULT_P1, p2 = MCOL_DEFAULT_P2, p12 = MCOL_DEFAULT_P12))
  pp <- rs$abf_res$summary
  expect_true(all(c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf") %in% names(pp)))
  expect_equal(sum(as.numeric(pp[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])), 1, tolerance = 1e-6)
})

test_that("changing p1/p2/p12 inputs away from coloc's defaults actually changes the reported posterior probabilities", {
  res_default <- run_methyl_coloc_upload(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_diff <- run_methyl_coloc_upload(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4)

  pp4_default <- unname(res_default$abf_res$summary["PP.H4.abf"])
  pp4_diff <- unname(res_diff$abf_res$summary["PP.H4.abf"])
  expect_false(isTRUE(all.equal(pp4_default, pp4_diff)))
  expect_equal(res_default$priors, list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  expect_equal(res_diff$priors, list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4))
})

test_that("omitting the prior inputs entirely (pre-existing behaviour) reproduces coloc's own conventional defaults exactly", {
  res_no_input <- run_methyl_coloc_upload()
  res_explicit_default <- run_methyl_coloc_upload(p1 = MCOL_DEFAULT_P1, p2 = MCOL_DEFAULT_P2, p12 = MCOL_DEFAULT_P12)

  expect_equal(res_no_input$priors, list(p1 = MCOL_DEFAULT_P1, p2 = MCOL_DEFAULT_P2, p12 = MCOL_DEFAULT_P12))
  expect_equal(unname(res_no_input$abf_res$summary), unname(res_explicit_default$abf_res$summary), tolerance = 1e-12)
  expect_equal(res_no_input$snp_df$snp_pp_h4, res_explicit_default$snp_df$snp_pp_h4, tolerance = 1e-12)
})

test_that("an out-of-range prior (p12 = 0.5) fails with a clear validate() message, not a raw coloc.abf crash", {
  err <- run_methyl_coloc_upload(p1 = 1e-4, p2 = 1e-4, p12 = 0.5, expect_run_error = TRUE)
  expect_true(inherits(err, "shiny.silent.error"))
  expect_match(conditionMessage(err), "p12", fixed = TRUE)
})

test_that("a prior outside (0,1) (e.g. p1 = 0) fails with a clear validate() message", {
  err <- run_methyl_coloc_upload(p1 = 0, p2 = 1e-4, p12 = 1e-5, expect_run_error = TRUE)
  expect_true(inherits(err, "shiny.silent.error"))
  expect_match(conditionMessage(err), "p1", fixed = TRUE)
})
