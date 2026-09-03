## Regression coverage for a gap found in the forensic audit (2026-09-03):
## the Methylomics Mendelian Randomization module (mod_methyl_mr.R) had NO
## test file at all, despite being the most methodologically dense module in
## Methylomics (real TwoSampleMR-based harmonisation, F-statistic
## weak-instrument flagging, and heterogeneity/pleiotropy/leave-one-out
## sensitivity analysis gated on per-CpG instrument count). These tests drive
## the real upload data_source through testServer(), walking the module's own
## stage sequence (Data -> Instruments -> Clump -> Harmonise -> MR -> Sensitivity)
## exactly as the UI would, on synthetic but well-formed two-sample-MR-format
## mQTL/GWAS files for a single CpG.
##
## LD clumping (ieugwasr::ld_clump()) requires network access to the IEU
## OpenGWAS API, unavailable in this environment; the module already handles
## that failure gracefully (falls back to the unclumped instrument set with a
## status note) - clicking through the Clump stage below exercises exactly
## that graceful-fallback path, which is itself worth confirming works.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "11_Mendelian_Randomization", "mod_methyl_mr.R"))

## One CpG, 5 SNP instruments: snp1/snp2/snp3/snp5 are strong (F well above
## the default cutoff of 10); snp4 is a deliberately weak instrument (F well
## below 10). Outcome effects are exposure effects x a fixed causal multiplier
## (0.3) plus small noise - a real, roughly-linear relationship.
mmr_write_upload_files <- function(dir, cpg = "cg00000001", weak_beta = 0.02, seed = 77) {
  set.seed(seed)
  snps <- paste0("snp", 1:5)
  exp_beta <- c(0.50, 0.60, 0.55, weak_beta, 0.45)
  exp_se <- rep(0.05, 5)
  exp_pval <- 2 * stats::pnorm(-abs(exp_beta / exp_se))
  out_beta <- exp_beta * 0.3 + stats::rnorm(5, 0, 0.01)
  out_se <- rep(0.05, 5)
  out_pval <- 2 * stats::pnorm(-abs(out_beta / out_se))

  exp_df <- data.frame(cpg = cpg, snp = snps, beta = exp_beta, se = exp_se, pval = exp_pval,
                        ea = "A", oa = "G", eaf = 0.3, stringsAsFactors = FALSE)
  out_df <- data.frame(snp = snps, beta = out_beta, se = out_se, pval = out_pval,
                        ea = "A", oa = "G", eaf = 0.3, stringsAsFactors = FALSE)

  exp_path <- file.path(dir, "mqtl_exposure.csv"); write.csv(exp_df, exp_path, row.names = FALSE)
  out_path <- file.path(dir, "gwas_outcome.csv"); write.csv(out_df, out_path, row.names = FALSE)
  list(exp_path = exp_path, out_path = out_path, snps = snps, cpg = cpg,
       f_stat = (exp_beta / exp_se)^2)
}

mmr_mkfile <- function(path) {
  data.frame(name = basename(path), size = file.info(path)$size, type = "text/csv",
             datapath = path, stringsAsFactors = FALSE)
}

## Walks the module's real stage sequence (mirroring the UI's own tab order)
## via testServer(), stopping after whichever stage(s) the caller needs and
## returning every intermediate reactiveVal so tests can inspect each stage.
run_mmr_upload_stages <- function(fx, f_min_f = 10, f_pval = 1,
                                   through = c("data", "instruments", "clump", "harmonise", "mr", "sensitivity")) {
  through <- match.arg(through, several.ok = FALSE)
  stage_order <- c("data", "instruments", "clump", "harmonise", "mr", "sensitivity")
  stop_at <- which(stage_order == through)

  methyl_dataset <- shiny::reactiveValues()
  methyl_results <- shiny::reactiveValues()
  out <- list()
  shiny::testServer(mod_methyl_mr_server, args = list(id = "mmr", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(data_source = "upload")
    session$setInputs(exp_file = mmr_mkfile(fx$exp_path), out_file = mmr_mkfile(fx$out_path))
    session$setInputs(exp_cpg = "cpg", exp_snp = "snp", exp_beta = "beta", exp_se = "se", exp_pval = "pval",
                       exp_ea = "ea", exp_oa = "oa", exp_eaf = "eaf")
    session$setInputs(out_snp = "snp", out_beta = "beta", out_se = "se", out_pval = "pval",
                       out_ea = "ea", out_oa = "oa", out_eaf = "eaf")
    session$setInputs(outcome_label = "Test outcome", up_binary_outcome = FALSE)
    session$setInputs(validate_btn = 1)
    out$data_state <<- data_state()

    if (stop_at >= 2) {
      ## f_exclude_weak = TRUE: by default weak instruments are only FLAGGED
      ## (weak_instrument = TRUE) but still counted as "retained" - a
      ## deliberate "flag, don't silently drop" design. Excluding them from
      ## "retained" (what actually reaches MR) requires this explicit opt-in.
      session$setInputs(f_pval = f_pval, f_min_f = f_min_f, f_maf = 0, f_require_eaf = FALSE,
                         f_region_mode = "cis", f_min_instruments = 1, f_exclude_weak = TRUE)
      session$setInputs(instruments_btn = 1)
      out$instruments_state <<- instruments_state()
    }

    if (stop_at >= 3) {
      session$setInputs(clump_r2 = 0.001, clump_kb = 10000, clump_p = 1, clump_pop = "EUR")
      session$setInputs(clump_btn = 1)
      out$clump_state <<- clump_state()
    }

    if (stop_at >= 4) {
      session$setInputs(harmonise_action = "2")
      session$setInputs(harmonise_btn = 1)
      out$harmonise_state <<- harmonise_state()
    }

    if (stop_at >= 5) {
      session$setInputs(mr_methods = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                                        "mr_weighted_mode", "mr_simple_mode"),
                         mr_ci_level = 0.95, mr_run_presso = FALSE)
      session$setInputs(mr_run_btn = 1)
      out$mr_state <<- mr_state()
    }

    if (stop_at >= 6) {
      session$setInputs(sensitivity_run_btn = 1)
      out$sensitivity_state <<- sensitivity_state()
    }
  })
  out
}

test_that("uploading synthetic mQTL/GWAS data validates and reaches the Data stage", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  out <- run_mmr_upload_stages(fx, through = "data")
  expect_equal(out$data_state$mode, "upload")
  expect_equal(nrow(out$data_state$exp_raw), 5L)
})

test_that("F-statistic filtering flags the deliberately weak instrument as excluded, others as retained", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  expect_lt(fx$f_stat[4], 10); expect_true(all(fx$f_stat[-4] > 10))

  out <- run_mmr_upload_stages(fx, f_min_f = 10, f_pval = 1, through = "instruments")
  d <- out$instruments_state$d
  weak_row <- d[d$SNP == "snp4", ]
  expect_true(weak_row$weak_instrument)
  expect_false(weak_row$retained)
  expect_true(all(d$retained[d$SNP != "snp4"]))
})

test_that("the weak instrument is genuinely absent from the harmonised instrument set reaching MR (not merely flagged and forgotten)", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  out <- run_mmr_upload_stages(fx, f_min_f = 10, f_pval = 1, through = "harmonise")
  h <- out$harmonise_state$harmonised
  expect_false("snp4" %in% h$SNP)
  expect_true(all(c("snp1", "snp2", "snp3", "snp5") %in% h$SNP))
})

test_that("raising the F-statistic cutoff further excludes more instruments, in the expected direction", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  ## F-stats: snp1=100, snp2=144, snp3=121, snp4=0.16, snp5=81. A cutoff of
  ## 110 retains only snp2 and snp3 (both the weak snp4 AND the merely
  ## moderate snp1/snp5 are now excluded too), unlike the default cutoff of
  ## 10 which retains everything except snp4.
  out <- run_mmr_upload_stages(fx, f_min_f = 110, f_pval = 1, through = "harmonise")
  h <- out$harmonise_state$harmonised
  expect_setequal(h$SNP, c("snp2", "snp3"))
})

test_that("MR analysis runs via TwoSampleMR::mr() and produces a real IVW estimate for the retained instrument set", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  out <- run_mmr_upload_stages(fx, f_min_f = 10, f_pval = 1, through = "mr")
  res <- out$mr_state$results
  expect_true(nrow(res) > 0)
  expect_true("Inverse variance weighted" %in% res$method | "mr_ivw" %in% res$method)
  ivw_row <- res[grepl("nverse variance", res$method) | res$method == "mr_ivw", ][1, ]
  ## The true causal multiplier injected into the fixture is 0.3 - the IVW
  ## estimate on 4 clean instruments should land in its ballpark, not be wildly
  ## off (a broken harmonisation - e.g. an unflipped allele - would typically
  ## produce a wrong-signed or wildly-scaled estimate instead).
  expect_gt(ivw_row$b, 0.15); expect_lt(ivw_row$b, 0.45)
})

test_that("with >=3 retained instruments for a CpG, heterogeneity/pleiotropy/leave-one-out are populated", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  out <- run_mmr_upload_stages(fx, f_min_f = 10, f_pval = 1, through = "sensitivity")
  ss <- out$sensitivity_state
  expect_true(fx$cpg %in% ss$cpgs_ge3)
  expect_false(is.null(ss$het)); expect_true(nrow(ss$het) > 0)
  expect_false(is.null(ss$pleio)); expect_true(nrow(ss$pleio) > 0)
  expect_false(is.null(ss$loo)); expect_true(nrow(ss$loo) > 0)
  expect_false(is.null(ss$single_snp))
})

test_that("with <3 retained instruments for the only CpG, heterogeneity/pleiotropy/leave-one-out are correctly unavailable, not omitted or erroring", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  ## Cutoff high enough that only snp2 (F ~= 144) and snp1 (F = 100) clear it -
  ## exactly 2 retained instruments for the CpG's only tested locus.
  out <- run_mmr_upload_stages(fx, f_min_f = 110, f_pval = 1, through = "sensitivity")
  ss <- out$sensitivity_state
  expect_false(fx$cpg %in% ss$cpgs_ge3)
  expect_null(ss$het)
  expect_null(ss$pleio)
  expect_null(ss$loo)
  ## Single-SNP estimates are documented as "always shown", regardless of tier.
  expect_false(is.null(ss$single_snp))
})

test_that("LD clumping gracefully falls back to the unclumped instrument set when the OpenGWAS API is unreachable (no network in this environment)", {
  dir <- withr::local_tempdir()
  fx <- mmr_write_upload_files(dir)
  out <- run_mmr_upload_stages(fx, f_min_f = 10, f_pval = 1, through = "clump")
  cs <- out$clump_state
  expect_equal(cs$mode, "upload")
  expect_equal(cs$status, "api_error")
  ## Falls back to the pre-clump retained set (4 instruments after F-stat filtering), not an empty/broken one.
  expect_equal(nrow(cs$d), 4L)
})
