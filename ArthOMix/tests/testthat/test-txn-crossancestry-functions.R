## Module 1 (Transcriptomics) - Cross-Ancestry Validation's single source-of-
## truth classifier, ca_classify(), covering every branch of its 2x2
## (replicated_EUR x transferable_EAS) -> ancestry_class mapping plus its

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "13_Cross_Ancestral_Validation", "mod_crossancestry.R"))

ca_fixture <- function() {
  data.frame(
    gene = c("SHARED", "UNTESTABLE_EAS", "NOT_TRANSFERABLE", "NOT_REPLICATED", "NA_STAHL"),
    dir_okada_eq_stahl = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    dir_okada_eq_bbj   = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    p_stahl = c(0.001, 0.001, 0.001, 0.001, NA_real_),
    p_bbj   = c(0.001, 0.001, 0.5, 0.001, 0.001),
    testable_EAS = c(TRUE, FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("ca_classify() assigns 'shared EUR+EAS' when both replicated and transferable", {
  out <- ca_classify(ca_fixture(), p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  row <- out[out$gene == "SHARED", ]
  expect_true(row$replicated_EUR)
  expect_true(row$transferable_EAS)
  expect_true(row$biomarker)
  expect_equal(row$ancestry_class, "shared EUR+EAS")
})

test_that("ca_classify() assigns 'EUR-replicated, untestable in EAS' when testable_EAS is FALSE", {
  out <- ca_classify(ca_fixture(), p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  row <- out[out$gene == "UNTESTABLE_EAS", ]
  expect_true(row$replicated_EUR)
  expect_false(row$transferable_EAS)
  expect_equal(row$ancestry_class, "EUR-replicated, untestable in EAS")
})

test_that("ca_classify() assigns 'EUR-replicated, not EAS' when testable but p_bbj fails the cutoff", {
  out <- ca_classify(ca_fixture(), p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  row <- out[out$gene == "NOT_TRANSFERABLE", ]
  expect_true(row$replicated_EUR)
  expect_false(row$transferable_EAS)
  expect_equal(row$ancestry_class, "EUR-replicated, not EAS")
})

test_that("ca_classify() assigns 'not EUR-replicated' when direction is required and doesn't match, regardless of EAS", {
  out <- ca_classify(ca_fixture(), p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  row <- out[out$gene == "NOT_REPLICATED", ]
  expect_false(row$replicated_EUR)
  expect_equal(row$ancestry_class, "not EUR-replicated")
})

test_that("ca_classify() coerces an NA p-value (e.g. zero surviving instrument SNPs) to FALSE/not-replicated, not NA", {
  out <- ca_classify(ca_fixture(), p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  row <- out[out$gene == "NA_STAHL", ]
  expect_false(is.na(row$replicated_EUR))
  expect_false(row$replicated_EUR)
  expect_equal(row$ancestry_class, "not EUR-replicated")
})

test_that("require_dir = FALSE ignores direction concordance entirely", {
  df <- ca_fixture()
  out <- ca_classify(df, p_eur = 0.05, p_eas = 0.05, require_dir = FALSE)
  row <- out[out$gene == "NOT_REPLICATED", ]
  expect_true(row$replicated_EUR)
})

test_that("threshold sensitivity: a gene just above p_eur is not replicated, just below is", {
  df <- ca_fixture()[1, ]
  above <- ca_classify(df, p_eur = 0.0005, p_eas = 0.05, require_dir = TRUE)
  expect_false(above$replicated_EUR)
  below <- ca_classify(df, p_eur = 0.01, p_eas = 0.05, require_dir = TRUE)
  expect_true(below$replicated_EUR)
})

test_that("ca_classify() defaults testable_EUR to TRUE when the column is absent (bundled panel)", {
  out <- ca_classify(ca_fixture(), p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  expect_true(all(out$testable_EUR))
})

test_that("ca_classify() assigns 'untestable in EUR' rather than 'not EUR-replicated' when testable_EUR is FALSE", {
  df <- ca_fixture()[1, ]
  df$testable_EUR <- FALSE
  out <- ca_classify(df, p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  expect_false(out$replicated_EUR)
  expect_false(out$biomarker)
  expect_equal(out$ancestry_class, "untestable in EUR")
})

test_that("ca_classify() never labels a gene 'untestable in EUR' when testable_EUR is TRUE and it fails on p-value", {
  row <- ca_fixture()[ca_fixture()$gene == "NOT_REPLICATED", ]
  out <- ca_classify(row, p_eur = 0.05, p_eas = 0.05, require_dir = TRUE)
  expect_equal(out$ancestry_class, "not EUR-replicated")
})
