## Module 2 (Methylomics) - DMP tab's pure functions: significance filter,
## genomic-inflation factor, sex-column/covariate detection, top-CpG
## ranking, and the chunked limma::lmFit() reused by mod_methyl_dmr.R - its

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
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

## --- Cell-type-fraction covariates in mod_methyl_dmp_prepare_subset() -----------------------
## Wires methyl_results$celltype$fractions (Cell-Type Deconvolution's per-sample output) into
## the same covariate-adjustment path already used for sample-sheet columns.

mod_methyl_dmp_celltype_fixture <- function(n_per_group = 15, seed = 400) {
  set.seed(seed)
  n <- n_per_group * 2
  samples <- paste0("S", 1:n)
  group <- rep(c("HC", "RA"), each = n_per_group)

  # CD4T fraction is confounded with disease group (higher in RA); Mono varies independently of
  # group/CD4T; Bcell fills to 1 - a realistic compositional cell-fraction table (rows sum to ~1).
  # Mono must carry its own independent variation (not a constant): a 3-part composition has 2
  # true degrees of freedom, and dropping exactly one column (the standard fix for the sum-to-1
  # collinearity) only recovers a full-rank design if the two *remaining* columns are not, in
  # turn, an exact affine pair of each other - which they would be if the third part never varied.
  cd4t <- pmin(pmax(c(rnorm(n_per_group, 0.30, 0.03), rnorm(n_per_group, 0.55, 0.03)), 0.05), 0.9)
  mono <- pmin(pmax(rnorm(n, 0.2, 0.03), 0.05), 0.35)
  bcell <- 1 - cd4t - mono
  fractions <- cbind(CD4T = cd4t, Mono = mono, Bcell = bcell)
  rownames(fractions) <- samples

  # 49 pure-noise probes + one probe ("cg1") driven by CD4T fraction with essentially no direct
  # group term - its apparent group association is entirely mediated through cell composition,
  # the textbook EWAS confound that cell-type adjustment exists to correct.
  n_probes <- 50
  m <- matrix(rnorm(n_probes * n, 0, 0.3), n_probes, n, dimnames = list(paste0("cg", 1:n_probes), samples))
  m["cg1", ] <- 4 * cd4t + rnorm(n, 0, 0.05)

  sheet <- data.frame(sample = samples, group = group, stringsAsFactors = FALSE)
  list(methyl_dataset = list(beta = m, sample_sheet = sheet, input_scale = "m"), fractions = fractions)
}

mod_methyl_dmp_fit_group_effect <- function(sub) {
  grp <- sub$grp
  cov_df <- sub$cov_df
  design_grp <- stats::model.matrix(~0 + grp); colnames(design_grp) <- levels(grp)
  if (!is.null(cov_df)) {
    safe <- sprintf("`%s`", colnames(cov_df))
    design_cov <- stats::model.matrix(stats::as.formula(paste("~", paste(safe, collapse = " + "))), data = cov_df)
    design_cov <- design_cov[, setdiff(colnames(design_cov), "(Intercept)"), drop = FALSE]
    design <- cbind(design_grp, design_cov)
  } else {
    design <- design_grp
  }
  fit <- limma::lmFit(sub$m, design)
  cm <- limma::makeContrasts(contrasts = "RA-HC", levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  list(design = design, tt = tt)
}

test_that("mod_methyl_dmp_prepare_subset() folds a selected cell-type-fraction covariate into cov_df, and it measurably changes the fitted group effect at a cell-type-confounded CpG", {
  fx <- mod_methyl_dmp_celltype_fixture()
  anno_none <- list(ok = FALSE, reason = "n/a")

  sub_unadj <- mod_methyl_dmp_prepare_subset(
    fx$methyl_dataset, "__all__", NULL, "group", "HC", "RA", character(0),
    min_valid_pct = 80, min_variance = 0, snp_filter = FALSE, anno = anno_none)
  sub_adj <- mod_methyl_dmp_prepare_subset(
    fx$methyl_dataset, "__all__", NULL, "group", "HC", "RA", "celltype__CD4T",
    min_valid_pct = 80, min_variance = 0, snp_filter = FALSE, anno = anno_none,
    celltype_fractions = fx$fractions)

  # (a) the design/covariate matrix actually contains the cell-type covariate column - not
  # silently dropped.
  expect_null(sub_unadj$cov_df)
  expect_true("celltype_CD4T" %in% colnames(sub_adj$cov_df))

  res_unadj <- mod_methyl_dmp_fit_group_effect(sub_unadj)
  res_adj <- mod_methyl_dmp_fit_group_effect(sub_adj)
  expect_true("celltype_CD4T" %in% colnames(res_adj$design))

  # (b) the fit succeeds without rank deficiency for this reasonable covariate subset.
  expect_equal(qr(res_adj$design)$rank, ncol(res_adj$design))

  # (c) adjusting for the confounded cell-type covariate measurably changes the fitted group
  # effect at the confounded CpG (cg1, driven by CD4T fraction with no direct group term) -
  # proving the covariate is genuinely absorbed into the model, not a cosmetic no-op.
  t_unadj <- res_unadj$tt["cg1", "t"]
  t_adj   <- res_adj$tt["cg1", "t"]
  expect_gt(abs(t_unadj), abs(t_adj) + 1)
})

test_that("mod_methyl_dmp_prepare_subset() aligns cell-type fractions to the DMP sample subset by sample ID, not row position", {
  fx <- mod_methyl_dmp_celltype_fixture()
  anno_none <- list(ok = FALSE, reason = "n/a")

  # Shuffle the fractions matrix's row order relative to the beta matrix/sample sheet - if
  # alignment were done positionally instead of match()-ing sample IDs, this would silently
  # scramble which CD4T value gets attached to which sample.
  shuffled <- fx$fractions[sample(nrow(fx$fractions)), , drop = FALSE]

  sub_inorder <- mod_methyl_dmp_prepare_subset(
    fx$methyl_dataset, "__all__", NULL, "group", "HC", "RA", "celltype__CD4T",
    min_valid_pct = 80, min_variance = 0, snp_filter = FALSE, anno = anno_none,
    celltype_fractions = fx$fractions)
  sub_shuffled <- mod_methyl_dmp_prepare_subset(
    fx$methyl_dataset, "__all__", NULL, "group", "HC", "RA", "celltype__CD4T",
    min_valid_pct = 80, min_variance = 0, snp_filter = FALSE, anno = anno_none,
    celltype_fractions = shuffled)

  expect_equal(sub_inorder$cov_df$celltype_CD4T, sub_shuffled$cov_df$celltype_CD4T)
})

test_that("mod_methyl_dmp_prepare_subset() auto-drops one cell type (with an explanatory note) when every estimated cell-type fraction is selected, instead of crashing or leaving a rank-deficient design", {
  fx <- mod_methyl_dmp_celltype_fixture()
  anno_none <- list(ok = FALSE, reason = "n/a")

  all_ct_covs <- paste0("celltype__", colnames(fx$fractions))
  sub <- mod_methyl_dmp_prepare_subset(
    fx$methyl_dataset, "__all__", NULL, "group", "HC", "RA", all_ct_covs,
    min_valid_pct = 80, min_variance = 0, snp_filter = FALSE, anno = anno_none,
    celltype_fractions = fx$fractions)

  # One fewer column than the number of cell types selected - one was dropped as the implicit
  # reference, avoiding perfect collinearity from the compositional (sum-to-1) constraint.
  expect_equal(ncol(sub$cov_df), length(all_ct_covs) - 1)
  expect_false(is.null(sub$celltype_note))
  expect_match(sub$celltype_note, "automatically dropped", fixed = TRUE)

  res <- mod_methyl_dmp_fit_group_effect(sub)
  expect_equal(qr(res$design)$rank, ncol(res$design))
})

test_that("mod_methyl_dmp_prepare_subset() errors clearly (not a crash) when a cell-type covariate is selected but no fractions are available", {
  fx <- mod_methyl_dmp_celltype_fixture()
  anno_none <- list(ok = FALSE, reason = "n/a")
  err <- tryCatch(
    mod_methyl_dmp_prepare_subset(
      fx$methyl_dataset, "__all__", NULL, "group", "HC", "RA", "celltype__CD4T",
      min_valid_pct = 80, min_variance = 0, snp_filter = FALSE, anno = anno_none,
      celltype_fractions = NULL),
    error = function(e) e)
  expect_s3_class(err, "validation")
  expect_true(grepl("run Cell-Type Deconvolution first", conditionMessage(err)))
})

test_that("mod_methyl_dmp_split_covariates()/mod_methyl_dmp_celltype_choices()/mod_methyl_dmp_covariate_display() round-trip cell-type covariate identifiers", {
  fractions <- matrix(0, 2, 2, dimnames = list(c("S1", "S2"), c("CD4T", "Bcell")))
  choices <- mod_methyl_dmp_celltype_choices(fractions)
  expect_setequal(unname(choices), c("celltype__CD4T", "celltype__Bcell"))
  expect_equal(unname(choices[["Estimated cell-type: CD4T"]]), "celltype__CD4T")

  split <- mod_methyl_dmp_split_covariates(c("age", "celltype__CD4T", "batch"))
  expect_setequal(split$sheet_cols, c("age", "batch"))
  expect_equal(split$celltype_names, "CD4T")

  expect_equal(mod_methyl_dmp_covariate_display(c("age", "celltype_CD4T")),
               c("age", "Estimated cell-type: CD4T"))
})
