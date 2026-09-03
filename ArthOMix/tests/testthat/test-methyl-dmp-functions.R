## Module 2 (Methylomics) - DMP tab's pure functions: significance filter,
## genomic-inflation factor, sex-column/covariate detection, top-CpG
## ranking, and the chunked limma::lmFit() reused by mod_methyl_dmr.R - its

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "05_Differential_Methylation_Position", "mod_methyl_dmp.R"))

dmp_fixture_df <- function() {
  data.frame(
    cpg = paste0("cg", 1:6),
    fdr = c(0.001, 0.001, 0.2, NA, 0.01, 0.01),
    dbeta = c(0.1, -0.1, 0.3, 0.3, 0.005, -0.02),
    stringsAsFactors = FALSE
  )
}

test_that("mod_methyl_dmp_filter() applies FDR + effect-size thresholds and drops NA FDR rows", {
  df <- dmp_fixture_df()
  out <- mod_methyl_dmp_filter(df, "fdr", "dbeta", fdr_max = 0.05, effect_min = 0.05, direction = "any")
  expect_setequal(out$cpg, c("cg1", "cg2"))
})

test_that("mod_methyl_dmp_filter() direction='hyper'/'hypo' restrict to positive/negative dbeta", {
  df <- dmp_fixture_df()
  hyper <- mod_methyl_dmp_filter(df, "fdr", "dbeta", fdr_max = 0.05, effect_min = 0.05, direction = "hyper")
  hypo  <- mod_methyl_dmp_filter(df, "fdr", "dbeta", fdr_max = 0.05, effect_min = 0.05, direction = "hypo")
  expect_equal(hyper$cpg, "cg1")
  expect_equal(hypo$cpg, "cg2")
})

test_that("mod_methyl_lambda_gc() is ~1 for a null (uniform p-value) distribution", {
  set.seed(250)
  p <- runif(5000)
  lambda <- mod_methyl_lambda_gc(p)
  expect_true(abs(lambda - 1) < 0.1)
})

test_that("mod_methyl_lambda_gc() is well above 1 for an inflated (systematically small) p-value distribution", {
  set.seed(251)
  p <- rbeta(5000, 0.5, 5)
  lambda <- mod_methyl_lambda_gc(p)
  expect_gt(lambda, 1.2)
})

test_that("mod_methyl_lambda_gc() returns NA for an empty/all-invalid input rather than erroring", {
  expect_true(is.na(mod_methyl_lambda_gc(c(NA, NA, -1, 2))))
})

test_that("mod_methyl_dmp_sex_col() finds the first recognized sex-column name, or NULL", {
  expect_equal(mod_methyl_dmp_sex_col(data.frame(id = 1, Sex = "F")), "Sex")
  expect_null(mod_methyl_dmp_sex_col(data.frame(id = 1, notsex = "F")))
  expect_null(mod_methyl_dmp_sex_col(NULL))
})

test_that("mod_methyl_dmp_sex_choices() maps F/M-like values to Female/Male labels", {
  sheet <- data.frame(sex = c("F", "M", "F"))
  choices <- mod_methyl_dmp_sex_choices(sheet, "sex")
  expect_equal(unname(choices[["Female"]]), "F")
  expect_equal(unname(choices[["Male"]]), "M")
  expect_true("All samples" %in% names(choices))
})

test_that("mod_methyl_dmp_sex_choices() falls back to raw values when they don't cleanly map to F/M", {
  sheet <- data.frame(sex = c("XX", "XY", "XX"))
  choices <- mod_methyl_dmp_sex_choices(sheet, "sex")
  expect_setequal(unname(choices), c("__all__", "XX", "XY"))
})

test_that("mod_methyl_dmp_covariate_cols() excludes all-unique (ID-like) character columns but keeps numeric ones", {
  sheet <- data.frame(
    sample_id = paste0("S", 1:10),
    age = 1:10,
    batch = rep(c("A", "B"), 5),
    constant = rep("X", 10),
    stringsAsFactors = FALSE
  )
  cols <- mod_methyl_dmp_covariate_cols(sheet, exclude = character(0))
  expect_true("age" %in% cols)
  expect_true("batch" %in% cols)
  expect_false("sample_id" %in% cols)
  expect_false("constant" %in% cols)
})

test_that("mod_methyl_dmp_topplot() ranks by FDR (then |dbeta|) and returns exactly n CpGs", {
  set.seed(252)
  df <- data.frame(cpg = paste0("cg", 1:50), gene = NA_character_,
                     fdr = runif(50), dbeta = rnorm(50), stringsAsFactors = FALSE)
  out <- mod_methyl_dmp_topplot(df, rank_by = "fdr", n = 10)
  expect_length(out$cpgs, 10)
  expect_equal(out$cpgs, df$cpg[order(df$fdr, -abs(df$dbeta))][1:10])
})

test_that("mod_methyl_dmp_topplot() errors clearly when no CpG has both fdr and dbeta", {
  df <- data.frame(cpg = "cg1", gene = NA_character_, fdr = NA_real_, dbeta = NA_real_)
  expect_error(mod_methyl_dmp_topplot(df), class = "validation")
})

test_that("methyl_chunked_lmfit() produces bit-for-bit identical topTable() output to a whole-matrix fit", {
  set.seed(253)
  n_probes <- 250; n_samples <- 20
  m <- matrix(rnorm(n_probes * n_samples, mean = 0, sd = 1), n_probes, n_samples,
               dimnames = list(paste0("cg", 1:n_probes), paste0("S", 1:n_samples)))
  grp <- factor(rep(c("HC", "RA"), each = n_samples / 2))
  design <- model.matrix(~grp)

  whole_fit <- limma::eBayes(limma::lmFit(m, design))
  chunked_fit <- limma::eBayes(methyl_chunked_lmfit(m, design, chunk_size = 40))

  tt_whole <- limma::topTable(whole_fit, coef = 2, number = Inf, sort.by = "none")
  tt_chunked <- limma::topTable(chunked_fit, coef = 2, number = Inf, sort.by = "none")
  expect_equal(tt_whole$logFC, tt_chunked$logFC)
  expect_equal(tt_whole$P.Value, tt_chunked$P.Value)
  expect_equal(tt_whole$adj.P.Val, tt_chunked$adj.P.Val)
})

test_that("methyl_chunked_lmfit() takes the direct (non-chunked) path when the matrix already fits in one chunk", {
  m <- matrix(rnorm(100), 10, 10)
  design <- model.matrix(~factor(rep(c("A", "B"), 5)))
  out <- methyl_chunked_lmfit(m, design, chunk_size = 20000)
  expect_true(is(out, "MArrayLM"))
})
