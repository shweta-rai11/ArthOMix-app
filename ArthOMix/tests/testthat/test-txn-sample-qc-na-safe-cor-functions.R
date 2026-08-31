## Regression guard for a gap found in the transcriptomics audit (2026-08-26):
## compute_sample_qc()'s cohort-correlation check called cor() with the
## default use="everything", so a single NA anywhere in the top-variance
## subset made every pairwise correlation - and so the whole flag_cor QC
## check for the entire cohort - silently NA, with no warning to the user.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))

test_that("a single NA in the expression matrix does not NA out the whole cohort-correlation QC check", {
  d0 <- load_default_dataset()
  expr <- d0$expr[1:500, 1:20]
  expr[5, 3] <- NA
  qc <- compute_sample_qc(expr, mad_k = 3)
  expect_false(any(is.na(qc$mean_cor)))
})
