## Module 1 (Transcriptomics) - Mendelian Randomization's core estimator,
## estimate_mr_set() (global.R, shared with mod_crossancestry.R's live
## replication/transfer arm) - Wald ratio for a single instrument SNP, IVW +
## median + Egger for >=3, with heterogeneity/pleiotropy statistics.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))

mr_snp_fixture <- function(n_snp, seed = 150, beta_causal = 0.3) {
  set.seed(seed)
  beta.exposure <- runif(n_snp, 0.1, 0.5)
  se.exposure <- runif(n_snp, 0.01, 0.03)
  ## Each SNP's outcome effect follows the same causal beta, plus small noise.
  beta.outcome <- beta_causal * beta.exposure + rnorm(n_snp, sd = 0.005)
  se.outcome <- runif(n_snp, 0.01, 0.03)
  data.frame(SNP = paste0("rs", 1:n_snp), gene = "TESTGENE",
              beta.exposure = beta.exposure, se.exposure = se.exposure,
              beta.outcome = beta.outcome, se.outcome = se.outcome, stringsAsFactors = FALSE)
}

test_that("estimate_mr_set() with a single instrument computes the exact Wald ratio formula", {
  d <- mr_snp_fixture(1)
  out <- estimate_mr_set(d)
  expect_equal(out$primary_method, "Wald ratio")
  expect_equal(out$n_snp, 1L)
  b_expected <- d$beta.outcome[1] / d$beta.exposure[1]
  se_expected <- abs(d$se.outcome[1] / d$beta.exposure[1])
  row <- out$res_table[out$res_table$method == "Wald ratio", ]
  expect_equal(row$estimate, b_expected)
  expect_equal(row$se, se_expected)
  expect_equal(row$ci_low, b_expected - stats::qnorm(0.975) * se_expected)
})

test_that("estimate_mr_set() with 2 instruments runs IVW only (no Egger/median, which need >=3)", {
  d <- mr_snp_fixture(2)
  out <- estimate_mr_set(d)
  expect_equal(out$primary_method, "IVW")
  expect_setequal(out$res_table$method, "IVW")
  expect_null(out$heterogeneity)
})

test_that("estimate_mr_set() with >=3 instruments runs IVW + Weighted median + MR-Egger, with heterogeneity/pleiotropy stats", {
  d <- mr_snp_fixture(6)
  out <- estimate_mr_set(d)
  expect_setequal(out$res_table$method, c("IVW", "Weighted median", "MR-Egger"))
  expect_true(sum(out$res_table$primary) == 1)
  expect_equal(out$res_table$method[out$res_table$primary], "IVW")
  expect_true(all(c("Q", "Q_df", "Q_pval") %in% names(out$heterogeneity)))
  expect_true(all(c("intercept", "se", "p", "I2") %in% names(out$pleiotropy)))
  ## Every method's CI should bracket its own point estimate.
  expect_true(all(out$res_table$ci_low <= out$res_table$estimate))
  expect_true(all(out$res_table$estimate <= out$res_table$ci_high))
})

test_that("estimate_mr_set() recovers a strong, real causal signal with a significant IVW p-value", {
  d <- mr_snp_fixture(8, beta_causal = 0.5, seed = 151)
  out <- estimate_mr_set(d)
  ivw_row <- out$res_table[out$res_table$method == "IVW", ]
  expect_true(ivw_row$p < 0.05)
  expect_equal(sign(ivw_row$estimate), sign(0.5))
})

test_that("estimate_mr_set()'s confidence interval widens as ci_level increases (e.g. 99% wider than 95%)", {
  d <- mr_snp_fixture(1, seed = 152)
  ci95 <- estimate_mr_set(d, ci_level = 0.95)$res_table
  ci99 <- estimate_mr_set(d, ci_level = 0.99)$res_table
  width95 <- ci95$ci_high - ci95$ci_low
  width99 <- ci99$ci_high - ci99$ci_low
  expect_gt(width99, width95)
})

test_that("estimate_mr_set()'s Weighted mode is included only when include_mode=TRUE and n>=3", {
  d <- mr_snp_fixture(6, seed = 153)
  without_mode <- estimate_mr_set(d, include_mode = FALSE)
  with_mode <- estimate_mr_set(d, include_mode = TRUE)
  expect_false("Weighted mode" %in% without_mode$res_table$method)
  expect_true("Weighted mode" %in% with_mode$res_table$method)
})
