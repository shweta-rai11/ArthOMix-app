## Colocalization upload-route server tests on synthetic stats from real eQTL rsIDs, focused on the p1/p2/p12 priors.

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
  ## Use a deliberately noisy GWAS signal: a near-perfect echo saturates PP.H4 and makes the prior test vacuous.
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

make_identity_ld_file <- function(snp_ids) {
  n <- length(snp_ids)
  mat <- diag(1, n, n)
  df <- data.frame(SNP = snp_ids, mat)
  colnames(df) <- c("SNP", snp_ids)
  path <- tempfile(fileext = ".csv")
  write.csv(df, path, row.names = FALSE)
  path
}

run_methyl_coloc_upload <- function(p1 = NULL, p2 = NULL, p12 = NULL, expect_run_error = FALSE,
                                     use_susie = FALSE) {
  files <- make_synthetic_meth_coloc_files()
  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  out <- NULL
  susie_calls <- list()
  shiny::testServer(mod_methyl_coloc_server, args = list(id = "coloc", dataset = dataset, results = results), {
    if (isTRUE(use_susie)) {
      snp_ids <- read.csv(files$meth_path, stringsAsFactors = FALSE)$snp
      ld_path <- make_identity_ld_file(snp_ids)
      testthat::local_mocked_bindings(
        coloc.susie = function(...) {
          susie_calls[[length(susie_calls) + 1]] <<- list(...)
          list(summary = data.frame(nsnps = 3, PP.H4.abf = 0.5))
        },
        .package = "coloc"
      )
    }
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

    if (isTRUE(use_susie)) {
      session$setInputs(
        use_susie = TRUE,
        ld1_file = list(datapath = ld_path, name = "ld1.csv"),
        ld2_file = list(datapath = ld_path, name = "ld2.csv")
      )
    }

    if (isTRUE(expect_run_error)) {
      ## Call the module-local builder directly: run_btn wraps it in tryCatch, so validate() is unobservable via run_state.
      out <<- tryCatch(build_run_state_upload(), error = function(e) e)
    } else {
      session$setInputs(run_btn = 1)
      out <<- run_state()
    }
  })
  attr(out, "susie_calls") <- susie_calls
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

test_that("coloc.susie() is called with the user's chosen p12, not a hardcoded value diverging from coloc.abf's", {
  ## Guard: coloc.susie() must use the user's p12 (shown in the UI, stored in rs$priors), not a hardcoded 5e-6.
  rs <- run_methyl_coloc_upload(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4, use_susie = TRUE)
  calls <- attr(rs, "susie_calls")
  expect_length(calls, 1)
  expect_equal(calls[[1]]$p12, 1e-4)
  expect_equal(calls[[1]]$p12, rs$priors$p12)
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
