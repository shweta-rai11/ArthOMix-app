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
source_from_app_root(file.path("R", "transcriptomics", "mod_nomogram.R"))

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
