## Module 1 (Transcriptomics) - Cross-Tissue Validation's pure functions
## specific to this submodule (ct_cv_eval/ct_fit_sex reuse the same CV/model-
## fitting pattern already thoroughly tested in test-diagnostic-functions.R
## for mod_diagnostic.R's diag_cv_auc/diag_fit_sex, so are not re-verified
## here in the same depth - covering ct_biomarker_flag/ct_gene_auc/
## ct_discovery_table/ct_voom_de_table/ct_build_uploaded_val instead, the
## logic genuinely unique to cross-tissue validation).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "12_Cross_Tissue_Validation", "mod_crosstissue.R"))

## ---- ct_biomarker_flag() ---------------------------------------------------

test_that("ct_biomarker_flag() requires concordant direction, significant synovium DE, and AUC >= 0.70 all at once", {
  d <- data.frame(
    concordant = c(TRUE, FALSE, TRUE, TRUE, TRUE),
    syn_adjP   = c(0.01, 0.01, 0.20, 0.01, 0.01),
    auc_bestdir = c(0.75, 0.75, 0.75, 0.75, 0.65)
  )
  flags <- ct_biomarker_flag(d, sig_cut = 0.05)
  ## row1: all three pass -> TRUE. row2: discordant -> FALSE. row3: syn_adjP
  ## not significant -> FALSE. row4: identical to row1 (all pass) -> TRUE.
  ## row5: AUC below 0.70 -> FALSE.
  expect_equal(flags, c(TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that("ct_biomarker_flag() treats any NA input as not-a-biomarker rather than erroring", {
  d <- data.frame(concordant = NA, syn_adjP = 0.01, auc_bestdir = 0.9)
  expect_false(ct_biomarker_flag(d, sig_cut = 0.05))
  d2 <- data.frame(concordant = TRUE, syn_adjP = NA_real_, auc_bestdir = 0.9)
  expect_false(ct_biomarker_flag(d2, sig_cut = 0.05))
})

## ---- ct_gene_auc() ----------------------------------------------------------

test_that("ct_gene_auc() returns the best-direction AUC (always >= 0.5) regardless of the sign of separation", {
  set.seed(140)
  y <- factor(rep(c("HC", "RA"), each = 20), levels = c("HC", "RA"))
  up <- rnorm(40) + ifelse(y == "RA", 3, 0)     ## RA higher
  down <- rnorm(40) - ifelse(y == "RA", 3, 0)   ## RA lower
  auc_up <- ct_gene_auc(up, y)
  auc_down <- ct_gene_auc(down, y)
  expect_true(auc_up >= 0.5 && auc_up <= 1)
  expect_true(auc_down >= 0.5 && auc_down <= 1)
  ## Both should recover the same strong separation, just oriented differently.
  expect_equal(auc_up, auc_down, tolerance = 0.15)
})

## ---- ct_discovery_table() ---------------------------------------------------

test_that("ct_discovery_table() flags a requested gene absent from the synovium data as present=FALSE with all-NA stats", {
  val <- list(
    sex = factor(rep(c("F", "M"), each = 10)),
    grp = factor(rep(c("HC", "RA"), 10), levels = c("HC", "RA")),
    tt = data.frame(gene = c("GENE1", "GENE2"), logFC = c(1.2, -0.8), adj.P.Val = c(0.01, 0.2)),
    logcpm = matrix(rnorm(2 * 20), 2, 20, dimnames = list(c("GENE1", "GENE2"), NULL))
  )
  blood_dir <- list(logfc = c(GENE1 = 1.0, GENE2 = -0.5))
  out <- ct_discovery_table(c("GENE1", "GENE3"), sex_code = NULL, val = val, blood_dir = blood_dir)
  expect_equal(out$present, c(TRUE, FALSE))
  expect_true(is.na(out$syn_log2FC[2]))
  expect_true(is.na(out$concordant[2]))
})

test_that("ct_discovery_table() correctly flags direction concordance between blood and synovium logFC sign", {
  val <- list(
    sex = factor(rep(c("F", "M"), each = 10)),
    grp = factor(rep(c("HC", "RA"), 10), levels = c("HC", "RA")),
    tt = data.frame(gene = c("GENE1", "GENE2"), logFC = c(1.2, -0.8), adj.P.Val = c(0.01, 0.01)),
    logcpm = matrix(rnorm(2 * 20, mean = 8), 2, 20, dimnames = list(c("GENE1", "GENE2"), NULL))
  )
  ## GENE1: blood up, synovium up -> concordant. GENE2: blood up, synovium down -> discordant.
  blood_dir <- list(logfc = c(GENE1 = 1.0, GENE2 = 0.5))
  out <- ct_discovery_table(c("GENE1", "GENE2"), sex_code = NULL, val = val, blood_dir = blood_dir)
  expect_equal(out$concordant, c(TRUE, FALSE))
})

## ---- ct_voom_de_table() ------------------------------------------------------

test_that("ct_voom_de_table() rejects non-count (negative or fractional) input before attempting any fit", {
  counts <- matrix(rnorm(100 * 10), 100, 10)  ## negative + fractional values
  grp <- factor(rep(c("HC", "RA"), each = 5), levels = c("HC", "RA"))
  sex <- rep(c("F", "M"), 5)
  expect_error(ct_voom_de_table(counts, grp, sex), class = "validation")
})

test_that("ct_voom_de_table() fits successfully on real raw counts and returns a topTable with expected columns", {
  set.seed(141)
  n_genes <- 200; n_samples <- 16
  grp <- factor(rep(c("HC", "RA"), each = n_samples / 2), levels = c("HC", "RA"))
  sex <- rep(c("F", "M"), length.out = n_samples)
  counts <- matrix(rpois(n_genes * n_samples, lambda = 200), n_genes, n_samples,
                    dimnames = list(paste0("G", 1:n_genes), paste0("S", 1:n_samples)))
  out <- ct_voom_de_table(counts, grp, sex)
  expect_true(all(c("gene", "logFC", "P.Value", "adj.P.Val") %in% colnames(out$tt)))
  expect_equal(nrow(out$logcpm), nrow(out$tt))
})

## ---- ct_build_uploaded_val() --------------------------------------------------

ct_val_fixture <- function(n_per_cell = 4, seed = 142) {
  set.seed(seed)
  n <- n_per_cell * 4  ## 2 sexes x 2 groups
  genes <- paste0("G", 1:200)
  samples <- paste0("S", 1:n)
  grp <- rep(rep(c("HC", "RA"), each = n_per_cell), 2)
  sex <- rep(c("F", "M"), each = n_per_cell * 2)
  counts <- matrix(rpois(200 * n, lambda = 200), 200, n, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, sex = sex, group = grp, stringsAsFactors = FALSE)
  list(expr = counts, meta = meta)
}

test_that("ct_build_uploaded_val() assembles a val_synovium.rds-shaped object from a valid uploaded cohort", {
  fx <- ct_val_fixture()
  out <- ct_build_uploaded_val(fx$expr, fx$meta, id_col = "sample", sex_col = "sex", group_col = "group",
                                 ref_group = "HC", comp_group = "RA")
  expect_true(all(c("logcpm", "grp", "sex", "tt") %in% names(out)))
  expect_equal(levels(out$grp), c("HC", "RA"))
  expect_setequal(unique(out$sex), c("F", "M"))
})

test_that("ct_build_uploaded_val() rejects a cohort where one sex has fewer than 4 samples after group filtering", {
  fx <- ct_val_fixture(n_per_cell = 4)
  ## Recode all but one male sample to female, leaving only 1 male.
  male_idx <- which(fx$meta$sex == "M")
  fx$meta$sex[male_idx[-1]] <- "F"
  err <- tryCatch(
    ct_build_uploaded_val(fx$expr, fx$meta, id_col = "sample", sex_col = "sex", group_col = "group",
                            ref_group = "HC", comp_group = "RA"),
    error = function(e) e
  )
  expect_s3_class(err, "validation")
  expect_true(grepl("Each sex needs at least 4 samples", conditionMessage(err)))
})

test_that("ct_build_uploaded_val() rejects a sex column whose values don't start with F/M", {
  fx <- ct_val_fixture()
  fx$meta$sex <- "unknown"
  err <- tryCatch(
    ct_build_uploaded_val(fx$expr, fx$meta, id_col = "sample", sex_col = "sex", group_col = "group",
                            ref_group = "HC", comp_group = "RA"),
    error = function(e) e
  )
  expect_s3_class(err, "validation")
  expect_true(grepl('must contain values starting with', conditionMessage(err)))
})
