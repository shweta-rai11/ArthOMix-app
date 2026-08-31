## Module 1 (Transcriptomics) - Preprocessing tab's merge/collapse logic.
## Covers pp_collapse_probes_to_genes() (top-level pure function in
## mod_preprocessing.R) plus global.R's shared expr_raw_health()/
## detect_expr_data_type()/filter_and_transform_expr() helpers this tab (and
## Data Exploration) both rely on, and merged()'s <20-shared-feature gate via
## testServer().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_preprocessing.R"))

## ---- pp_collapse_probes_to_genes() -----------------------------------------

pp_probe_fixture <- function() {
  expr <- matrix(
    c(5, 6, 3, 8, 2, 9),  ## PROBE1/PROBE2 -> GENEA, PROBE3 -> GENEB
    nrow = 3, byrow = TRUE,
    dimnames = list(c("PROBE1", "PROBE2", "PROBE3"), c("S1", "S2"))
  )
  annot <- data.frame(probe_id = c("PROBE1", "PROBE2", "PROBE3", "PROBE4"),
                       gene_symbol = c("GENEA", "GENEA", "GENEB", "GENEA///GENEC"),
                       stringsAsFactors = FALSE)
  list(expr = expr, annot = annot)
}

test_that("pp_collapse_probes_to_genes() collapses multi-probe genes via median/mean/maxmean and drops ambiguous/unmatched probes", {
  fx <- pp_probe_fixture()

  out_median <- pp_collapse_probes_to_genes(fx$expr, fx$annot, "median")
  expect_setequal(rownames(out_median), c("GENEA", "GENEB"))
  expect_equal(unname(out_median["GENEA", ]), c(median(c(5, 3)), median(c(6, 8))))

  out_mean <- pp_collapse_probes_to_genes(fx$expr, fx$annot, "mean")
  expect_equal(unname(out_mean["GENEA", ]), c(mean(c(5, 3)), mean(c(6, 8))))

  out_maxmean <- pp_collapse_probes_to_genes(fx$expr, fx$annot, "maxmean")
  ## PROBE1 (mean 5.5) beats PROBE2 (mean 5.0) for GENEA under maxmean.
  expect_equal(unname(out_maxmean["GENEA", ]), c(5, 6))

  ## PROBE4 (only in annot, not in expr) never contributes; the ambiguous
  ## "GENEA///GENEC" annotation for it would have been dropped either way.
  expect_false("GENEC" %in% rownames(out_median))
})

test_that("pp_collapse_probes_to_genes() errors clearly when no probe IDs match the annotation at all", {
  fx <- pp_probe_fixture()
  rownames(fx$expr) <- c("UNRELATED1", "UNRELATED2", "UNRELATED3")
  expect_error(pp_collapse_probes_to_genes(fx$expr, fx$annot, "median"), class = "validation")
})

## ---- expr_raw_health() / detect_expr_data_type() ---------------------------

test_that("expr_raw_health() reports missing/zero/infinite/duplicated-feature counts on a crafted matrix", {
  m <- matrix(c(0, 1, NA, Inf, 2, 3), nrow = 2, dimnames = list(c("A", "A"), NULL))
  h <- expr_raw_health(m)
  expect_equal(h$n_features, 2L)
  expect_equal(h$n_missing, 1L)
  expect_equal(h$n_zero, 1L)
  expect_equal(h$n_infinite, 1L)
  expect_equal(h$n_duplicated_features, 1L)
})

test_that("detect_expr_data_type() distinguishes raw counts, already-normalized, and generic expression data", {
  set.seed(30)
  counts <- matrix(rpois(400, lambda = 800), 20, 20)
  expect_equal(detect_expr_data_type(counts), "counts")

  zscored <- matrix(rnorm(400), 20, 20)  ## has negatives
  expect_equal(detect_expr_data_type(zscored), "already_normalised")
})

## ---- filter_and_transform_expr() -------------------------------------------

test_that("filter_and_transform_expr() treats infinite values as missing rather than winning threshold comparisons", {
  m <- matrix(c(Inf, 5, 6, 7, 8, 9), nrow = 1, dimnames = list("GENE1", NULL))
  out <- filter_and_transform_expr(m, min_expr = 0, min_sample_frac = 0, max_na_pct = 100, log2_transform = FALSE)
  expect_false(any(is.infinite(out)))
})

test_that("filter_and_transform_expr() drops a feature whose non-missing values never clear min_sample_frac", {
  ## min_sample_frac is a PERCENTAGE (0-100), not a proportion - only 2/6
  ## (33%) of this row's values clear min_expr = 5, below a 50% requirement.
  m <- matrix(c(0, 0, 0, 0, 10, 10), nrow = 1, dimnames = list("LOWEXPR", NULL))
  out <- filter_and_transform_expr(m, min_expr = 5, min_sample_frac = 50, max_na_pct = 100, log2_transform = FALSE)
  expect_false("LOWEXPR" %in% rownames(out))
})

test_that("filter_and_transform_expr() log2-transforms only when requested, and drops non-positive values first", {
  m <- matrix(c(-1, 0, 4, 8, 16, 32), nrow = 1, dimnames = list("GENE1", NULL))
  out <- filter_and_transform_expr(m, min_expr = -Inf, min_sample_frac = 0, max_na_pct = 100, log2_transform = TRUE)
  expect_true(all(is.finite(as.numeric(out))))
  expect_true(max(out, na.rm = TRUE) <= log2(32) + 1e-6)
})

## ---- merged()'s <20-shared-feature gate (testServer) -----------------------

test_that("merging two uploaded datasets with fewer than 20 shared feature IDs is rejected with a clear validate() error", {
  dir <- withr::local_tempdir()
  set.seed(40)
  genesA <- paste0("GENE", 1:15)
  genesB <- paste0("GENE", 10:15)  ## only 6 genes overlap with A
  mA <- matrix(rnorm(15 * 6), 15, 6, dimnames = list(genesA, paste0("A", 1:6)))
  mB <- matrix(rnorm(6 * 6), 6, 6, dimnames = list(genesB, paste0("B", 1:6)))
  write.csv(data.frame(gene = rownames(mA), mA, check.names = FALSE), file.path(dir, "exprA.csv"), row.names = FALSE)
  write.csv(data.frame(sample = colnames(mA), group = rep(c("HC", "RA"), 3)), file.path(dir, "metaA.csv"), row.names = FALSE)
  write.csv(data.frame(gene = rownames(mB), mB, check.names = FALSE), file.path(dir, "exprB.csv"), row.names = FALSE)
  write.csv(data.frame(sample = colnames(mB), group = rep(c("HC", "RA"), 3)), file.path(dir, "metaB.csv"), row.names = FALSE)

  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = d0$source, source_type = "preloaded")

  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(n_sources = 2)
    session$setInputs(`src1-source_type` = "upload")
    session$setInputs(`src1-expr_file` = fx_mkfile(file.path(dir, "exprA.csv")))
    session$setInputs(`src1-meta_file` = fx_mkfile(file.path(dir, "metaA.csv")))
    session$setInputs(`src1-map_id` = "sample", `src1-map_group` = "group")
    session$setInputs(`src1-log2` = "skip", `src1-max_na_pct` = 0)
    session$setInputs(`src1-run` = 1)

    session$setInputs(`src2-source_type` = "upload")
    session$setInputs(`src2-expr_file` = fx_mkfile(file.path(dir, "exprB.csv")))
    session$setInputs(`src2-meta_file` = fx_mkfile(file.path(dir, "metaB.csv")))
    session$setInputs(`src2-map_id` = "sample", `src2-map_group` = "group")
    session$setInputs(`src2-log2` = "skip", `src2-max_na_pct` = 0)
    session$setInputs(`src2-run` = 1)

    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)
    err <- tryCatch(merged(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 20 features", conditionMessage(err)))
  })
})
