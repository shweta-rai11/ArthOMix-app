## Regression coverage for a gap found in the transcriptomics audit
## (2026-09-03): tests/testthat/test-txn-mr-functions.R exercises only the
## shared pure-math helper estimate_mr_set() (in global.R); the actual
## mod_mr_server() upload-path logic in mod_mr.R - TwoSampleMR harmonisation,
## F-statistic weak-instrument filtering, and pleiotropy-sensitivity output
## gated on instrument count - had zero test coverage. These tests drive the
## real upload data_source through testServer(), on synthetic but
## well-formed two-sample-MR-format exposure/outcome files.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "07_Mendelian_Randomization", "mod_mr.R"))

## 5 SNPs: snp1/snp2/snp3/snp5 are strong instruments (F-stat well above the
## default cutoff of 10); snp4 is a deliberately weak instrument (F-stat well
## below 10). Outcome effects are a fixed multiple (0.3) of exposure effects
## plus a small amount of noise, i.e. a real, roughly-linear (if noisy)
## causal relationship - not just F-stat-passing noise.
mr_write_upload_files <- function(dir, weak_beta = 0.02, seed = 55) {
  set.seed(seed)
  snps <- paste0("snp", 1:5)
  exp_beta <- c(0.50, 0.60, 0.55, weak_beta, 0.45)
  exp_se <- rep(0.05, 5)
  exp_z <- exp_beta / exp_se
  exp_pval <- 2 * stats::pnorm(-abs(exp_z))
  out_beta <- exp_beta * 0.3 + stats::rnorm(5, 0, 0.01)
  out_se <- rep(0.05, 5)
  out_pval <- 2 * stats::pnorm(-abs(out_beta / out_se))

  exp_df <- data.frame(snp = snps, beta = exp_beta, se = exp_se, pval = exp_pval,
                        ea = "A", oa = "G", eaf = 0.3, stringsAsFactors = FALSE)
  out_df <- data.frame(snp = snps, beta = out_beta, se = out_se, pval = out_pval,
                        ea = "A", oa = "G", eaf = 0.3, stringsAsFactors = FALSE)

  exp_path <- file.path(dir, "exposure.csv"); write.csv(exp_df, exp_path, row.names = FALSE)
  out_path <- file.path(dir, "outcome.csv"); write.csv(out_df, out_path, row.names = FALSE)
  list(exp_path = exp_path, out_path = out_path, snps = snps, exp_beta = exp_beta,
       f_stat = (exp_beta / exp_se)^2)
}

mr_mkfile <- function(path) {
  data.frame(name = basename(path), size = file.info(path)$size, type = "text/csv",
             datapath = path, stringsAsFactors = FALSE)
}

run_mr_upload <- function(fx, fstat_cut = "10", pval_cut = "1", include_mode = FALSE, run_presso = FALSE) {
  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  out <- NULL
  shiny::testServer(mod_mr_server, args = list(id = "mr", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload")
    session$setInputs(exp_file = mr_mkfile(fx$exp_path), out_file = mr_mkfile(fx$out_path))
    session$setInputs(exp_snp = "snp", exp_beta = "beta", exp_se = "se", exp_pval = "pval",
                       exp_ea = "ea", exp_oa = "oa", exp_eaf = "eaf")
    session$setInputs(out_snp = "snp", out_beta = "beta", out_se = "se", out_pval = "pval",
                       out_ea = "ea", out_oa = "oa", out_eaf = "eaf")
    session$setInputs(upload_label = "Test exposure", ci_level = "0.95",
                       pval_cut = pval_cut, fstat_cut = fstat_cut,
                       do_clump = FALSE, include_mode = include_mode, run_presso = run_presso)
    session$setInputs(run_btn = 1)
    out <<- mr_result()
  })
  out
}

test_that("F-statistic filtering excludes a deliberately weak instrument from the final instrument set", {
  dir <- withr::local_tempdir()
  fx <- mr_write_upload_files(dir)
  expect_lt(fx$f_stat[4], 10)              # snp4 is the weak one, by construction
  expect_true(all(fx$f_stat[-4] > 10))     # every other SNP clears the default cutoff

  res <- run_mr_upload(fx, fstat_cut = "10", pval_cut = "1")
  expect_equal(res$n_before, 5L)           # all 5 SNPs harmonise successfully
  expect_equal(res$n_after, 4L)            # exactly the weak one is filtered out
  expect_false("snp4" %in% res$d$SNP)
  expect_true(all(c("snp1", "snp2", "snp3", "snp5") %in% res$d$SNP))
  expect_true(all(res$d$Fstat >= 10))
})

test_that("raising the F-statistic cutoff further excludes more instruments, in the expected direction", {
  dir <- withr::local_tempdir()
  fx <- mr_write_upload_files(dir)
  ## F-stats: snp1=100, snp2=144, snp3=121, snp4=0.16, snp5=81. A cutoff of
  ## 110 retains only snp2 and snp3 - both the weak snp4 AND the merely
  ## moderate snp1/snp5 are now excluded too, unlike the default cutoff of 10
  ## which retains everything except snp4.
  res <- run_mr_upload(fx, fstat_cut = "110", pval_cut = "1")
  expect_true(all(res$d$Fstat >= 110))
  expect_setequal(res$d$SNP, c("snp2", "snp3"))
})

test_that("with >=3 retained instruments, heterogeneity and pleiotropy sensitivity statistics are populated and numeric", {
  dir <- withr::local_tempdir()
  fx <- mr_write_upload_files(dir)
  res <- run_mr_upload(fx, fstat_cut = "10", pval_cut = "1")
  expect_gte(res$n_after, 3L)
  expect_false(is.null(res$est$heterogeneity))
  expect_false(is.null(res$est$pleiotropy))
  expect_true(is.numeric(res$est$heterogeneity$Q) && !is.na(res$est$heterogeneity$Q))
  expect_true(is.numeric(res$est$heterogeneity$Q_pval) && !is.na(res$est$heterogeneity$Q_pval))
  expect_true(is.numeric(res$est$pleiotropy$intercept) && !is.na(res$est$pleiotropy$intercept))
  expect_true(is.numeric(res$est$pleiotropy$p) && !is.na(res$est$pleiotropy$p))
  expect_true("Weighted median" %in% res$est$res_table$method)
  expect_true("MR-Egger" %in% res$est$res_table$method)
})

test_that("with <3 retained instruments, heterogeneity/pleiotropy are correctly reported as unavailable rather than omitted or erroring", {
  dir <- withr::local_tempdir()
  ## F-stat cutoff of 110 retains only snp2 (F ~= 144) and snp3 (F ~= 121) -
  ## exactly 2 instruments.
  fx <- mr_write_upload_files(dir)
  res <- run_mr_upload(fx, fstat_cut = "110", pval_cut = "1")
  expect_equal(res$n_after, 2L)
  expect_null(res$est$heterogeneity)
  expect_null(res$est$pleiotropy)
  expect_true("IVW" %in% res$est$res_table$method)
  expect_false("MR-Egger" %in% res$est$res_table$method)
})

test_that("a single retained instrument produces a Wald-ratio estimate, not IVW/heterogeneity, without erroring", {
  dir <- withr::local_tempdir()
  fx <- mr_write_upload_files(dir)
  ## Cutoff high enough that only the single strongest instrument (snp2, F ~= 144) survives.
  res <- run_mr_upload(fx, fstat_cut = "143", pval_cut = "1")
  expect_equal(res$n_after, 1L)
  expect_true("Wald ratio" %in% res$est$res_table$method)
  expect_null(res$est$heterogeneity)
  expect_null(res$est$pleiotropy)
})

test_that("a palindromic SNP is harmonised (kept or dropped via mr_keep) without crashing the upload path", {
  dir <- withr::local_tempdir()
  fx <- mr_write_upload_files(dir)
  ## Append a 6th, palindromic (A/T) SNP with an ambiguous (near-0.5) EAF on
  ## the exposure side - TwoSampleMR::harmonise_data(action = 2) should not
  ## error on this; it will either resolve it via EAF or drop it as
  ## unresolvable (mr_keep = FALSE), and mod_mr.R already filters to
  ## mr_keep == TRUE before returning - either outcome is acceptable here,
  ## the point is that it must not crash the pipeline.
  exp_extra <- data.frame(snp = "snp6", beta = 0.5, se = 0.05,
                           pval = 2 * stats::pnorm(-abs(0.5 / 0.05)),
                           ea = "A", oa = "T", eaf = 0.5, stringsAsFactors = FALSE)
  out_extra <- data.frame(snp = "snp6", beta = 0.15, se = 0.05,
                           pval = 2 * stats::pnorm(-abs(0.15 / 0.05)),
                           ea = "A", oa = "T", eaf = 0.5, stringsAsFactors = FALSE)
  exp_df <- rbind(utils::read.csv(fx$exp_path, stringsAsFactors = FALSE), exp_extra)
  out_df <- rbind(utils::read.csv(fx$out_path, stringsAsFactors = FALSE), out_extra)
  write.csv(exp_df, fx$exp_path, row.names = FALSE)
  write.csv(out_df, fx$out_path, row.names = FALSE)

  res <- tryCatch(run_mr_upload(fx, fstat_cut = "10", pval_cut = "1"), error = function(e) e)
  expect_false(inherits(res, "error"))
  expect_true(is.list(res) && !is.null(res$d))
  ## snp6 either survived harmonisation (kept, with F-stat filtering applied
  ## like any other SNP) or was dropped as an unresolved palindrome - either
  ## way n_before/n_after must be internally consistent (no crash, no silent
  ## corruption of the other 4-5 instruments).
  expect_true(res$n_before %in% c(5L, 6L))
  expect_true(all(c("snp1", "snp2", "snp3") %in% res$d$SNP))
})
