## Regression guards for two issues found in the transcriptomics audit
## (2026-08-26) in mod_nomogram.R's nom_fit_core():
##
## 1. The ridge penalty used to be hardcoded 0 (female) / 5 (male) with no
##    citation or cross-validation, confounding any male-vs-female comparison
##    of the resulting nomograms with an arbitrary sex-keyed regularization
##    choice. It is now selected from the fit's own data via rms::pentrace(),
##    the same procedure regardless of sex.
## 2. rms::datadist was registered under one fixed .GlobalEnv name, which
##    races under concurrent multi-session use (one R process, standard
##    Shiny deployment). It is now a per-call-unique name.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "16_Nomogram", "mod_nomogram.R"))

nomogram_test_df <- function(seed, n = 50) {
  set.seed(seed)
  X <- matrix(rnorm(n * 4), n, 4, dimnames = list(paste0("S", 1:n), paste0("gene", 1:4)))
  y01 <- rbinom(n, 1, plogis(0.9 * X[, 1] - 0.6 * X[, 2]))
  df <- as.data.frame(X); names(df) <- make.names(names(df)); df$y <- y01
  df
}

test_that("penalty = 'auto' resolves to a single numeric value chosen from the data", {
  df <- nomogram_test_df(3)
  params <- list(seed = 1234, calibrate_B = 40, dca_step = 0.05, impact_N = 100, impact_B = 20)
  res <- nom_fit_core(df, setdiff(names(df), "y"), "auto", "RA", params)
  expect_true(is.numeric(res$penalty))
  expect_length(res$penalty, 1)
  expect_false(identical(res$penalty, "auto"))
})

test_that("two sequential nom_fit_core() calls don't leak datadist state into each other", {
  df_a <- nomogram_test_df(3, n = 50)
  df_b <- nomogram_test_df(3, n = 40)
  params <- list(seed = 1234, calibrate_B = 40, dca_step = 0.05, impact_N = 100, impact_B = 20)

  res_a <- nom_fit_core(df_a, setdiff(names(df_a), "y"), "auto", "RA", params)
  res_b <- nom_fit_core(df_b, setdiff(names(df_b), "y"), "auto", "RA", params)

  expect_equal(res_a$n_samples, 50)
  expect_equal(res_b$n_samples, 40)
  expect_false(identical(res_a$c_stat, res_b$c_stat))
  ## Neither call should leave a stray datadist object behind in .GlobalEnv.
  expect_length(grep("^\\.arthomix_nomogram_dd", ls(envir = .GlobalEnv)), 0)
})

## Module 1 additions (2026-08-30): structural/scientific-contract checks on
## nom_dca()'s net-benefit formula (Vickers & Elkin 2006) and
## nom_clinical_impact()'s bootstrap bands, extending this file's existing
## nom_fit_core() regression coverage rather than duplicating a new file.

test_that("nom_dca() computes net benefit matching the Vickers & Elkin formula on a hand-worked example", {
  y <- c(1, 0, 1, 0)
  p <- c(0.9, 0.2, 0.8, 0.1)
  out <- nom_dca(y, p, th = 0.5)
  ## At threshold 0.5: positives flagged are indices 1,3 (both true cases) ->
  ## model NB = 2/4 - 0*(0.5/0.5) = 0.5; "treat all" NB = ev - (1-ev)*1 = 0.
  expect_equal(out$model, 0.5)
  expect_equal(out$all, 0)
})

test_that("nom_dca() 'treat none' is implicitly zero and 'treat all' degrades as the threshold rises", {
  set.seed(130)
  y <- rbinom(200, 1, 0.3)
  p <- runif(200)
  th <- c(0.1, 0.5, 0.9)
  out <- nom_dca(y, p, th)
  expect_length(out$model, 3)
  expect_length(out$all, 3)
  ## "Treat all"'s net benefit is monotonically decreasing in threshold (the
  ## false-positive penalty term pt/(1-pt) grows while the true-positive
  ## term - prevalence - is fixed).
  expect_true(all(diff(out$all) < 0))
})

test_that("nom_clinical_impact() computes deterministic nhigh/nevent counts exactly, with valid (lo<=hi) bootstrap bands", {
  set.seed(131)
  n <- 60
  X <- rnorm(n)
  y <- rbinom(n, 1, plogis(X))
  df <- data.frame(x = X, y = y)
  p <- as.numeric(plogis(X))
  th <- c(0.3, 0.5, 0.7)
  N <- 1000

  impact <- nom_clinical_impact(df, y ~ x, penalty = 0, p = p, y = y, th = th, N = N, B = 30, seed = 1234)

  ## nhigh/nevent are computed directly from p/y (not bootstrapped) - exactly reproducible.
  expect_equal(impact$nhigh, vapply(th, function(pt) mean(p >= pt) * N, numeric(1)))
  expect_equal(impact$nevent, vapply(th, function(pt) mean(p >= pt & y == 1) * N, numeric(1)))
  ## Bootstrap CI bands must bracket sensibly (2.5th <= 97.5th percentile).
  expect_true(all(impact$nhigh_lo <= impact$nhigh_hi))
  expect_true(all(impact$nevent_lo <= impact$nevent_hi))
  expect_equal(impact$N, N)
  expect_equal(impact$B, 30)
})
