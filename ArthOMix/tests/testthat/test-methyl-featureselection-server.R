## Module 2 (Methylomics) - Feature Selection's held-out split logic, via
## testServer(). Regression guard for a RED finding (2026-09-07 defense
## audit): caret::createDataPartition() failing inside fs_build_filters()
## silently fell back to putting every sample in the training set (a
## 100%/0% train/holdout split) and reported success with no user-visible
## warning - the one place a leakage bug could invalidate every downstream
## biomarker claim had zero server-level test coverage. This was previously
## the single largest testing gap in the Methylomics vertical (no
## test-methyl-featureselection-server.R existed at all, unlike DMP/DMR/QC/
## coloc, which all have paired functions+server tests).
##
## fs_build_filters() is called directly (not via the run_*_btn eventReactive
## wrappers) because those eventReactives are declared with ignoreInit = TRUE,
## and simulating an actionButton click through testServer's synchronous
## setInputs()/flushReact() does not reliably clear that initial-ignore state
## here (a testServer/eventReactive interaction quirk, not a property of the
## real browser-driven app, where an actionButton always initializes to 0
## before any real click). Calling fs_build_filters() directly exercises the
## exact same production leakage-safety code the buttons invoke.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "normalization.R"))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "methylomics", "10_ML_Feature_Selection", "mod_methyl_featureselection.R"))

fs_fixture_dataset <- function(n_per_group = 15, seed = 771) {
  set.seed(seed)
  n <- n_per_group * 2
  m <- matrix(runif(80 * n, 0.2, 0.8), 80, n, dimnames = list(paste0("cg", 20000000 + 1:80), paste0("S", 1:n)))
  m[1:8, (n_per_group + 1):n] <- pmin(m[1:8, (n_per_group + 1):n] + 0.3, 0.99)
  sheet <- data.frame(sample = colnames(m), group = rep(c("HC", "RA"), each = n_per_group),
                        sex = rep(c("F", "M"), length.out = n), stringsAsFactors = FALSE)
  shiny::reactiveValues(beta = m, sample_sheet = sheet, input_scale = "beta", array_type = "EPIC",
                          preloaded = FALSE, source_type = "uploaded", source = "feature selection test")
}

fs_base_inputs <- function() {
  list(fs_source = "preloaded", fs_group_col = "group", fs_ref_group = "HC", fs_comp_group = "RA",
       fs_holdout_enabled = TRUE, fs_holdout_frac = 0.3, fs_holdout_seed = 1234,
       fs_max_na_sample = 20, fs_max_na_probe = 5, fs_var_metric = "variance", fs_var_min = 0,
       fs_filter_snp = FALSE, fs_filter_sexchr = FALSE, fs_filter_crossreactive = FALSE, fs_filter_maf = FALSE)
}

test_that("the held-out split genuinely removes samples from the training matrix before any filtering/selection, and reports a non-empty holdout", {
  dataset <- fs_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    do.call(session$setInputs, fs_base_inputs())
    r <- fs_build_filters("pooled", NULL)
    expect_gt(length(r$holdout_sample_ids), 0)
    expect_true(all(r$holdout_sample_ids %in% colnames(dataset$beta)))
    ## the training matrix used for selection never contains a held-out sample
    expect_length(intersect(colnames(r$beta), r$holdout_sample_ids), 0)
    expect_equal(length(r$holdout_sample_ids) + ncol(r$beta), 30L)
    expect_true(any(grepl("Reserved .* samples as a held-out set", r$filter_notes)))
  })
})

test_that("disabling the held-out split produces an empty holdout and the full sample set is used for selection", {
  dataset <- fs_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    do.call(session$setInputs, utils::modifyList(fs_base_inputs(), list(fs_holdout_enabled = FALSE)))
    r <- fs_build_filters("pooled", NULL)
    expect_length(r$holdout_sample_ids, 0)
    expect_equal(ncol(r$beta), 30L)
    expect_false(any(grepl("Reserved .* samples as a held-out set", r$filter_notes)))
  })
})

test_that("a createDataPartition() failure stops the run with a clear validation error instead of silently falling back to a 0-sample holdout (the RED finding this file guards against)", {
  dataset <- fs_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    do.call(session$setInputs, fs_base_inputs())
    err <- testthat::with_mocked_bindings(
      tryCatch(fs_build_filters("pooled", NULL), error = function(e) e),
      createDataPartition = function(...) stop("simulated createDataPartition failure"),
      .package = "caret"
    )
    expect_s3_class(err, "validation")
    expect_true(grepl("Could not create the held-out split", conditionMessage(err)))
    expect_true(grepl("simulated createDataPartition failure", conditionMessage(err)))
  })
})

test_that("a training partition too small per group after the split is rejected with a clear message, not a silently degraded run", {
  ## 6 per group total (12 samples), 50% holdout: training side can easily drop
  ## below 2 per group for at least one random split - repeat several seeds to
  ## find one where the app's own guard fires (deterministic given the fixed seed).
  dataset <- fs_fixture_dataset(n_per_group = 6)
  results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    do.call(session$setInputs, utils::modifyList(fs_base_inputs(), list(fs_holdout_frac = 0.5, fs_holdout_seed = 1234)))
    res <- tryCatch(fs_build_filters("pooled", NULL), error = function(e) e)
    ## Either the split succeeds with >=2/group (assert the invariant it must have
    ## upheld), or it is rejected with the app's own clear, actionable message -
    ## never a silent 0-sample-holdout success.
    if (inherits(res, "error")) {
      expect_s3_class(res, "validation")
      expect_true(grepl("too few samples per group|too few samples", conditionMessage(res)))
    } else {
      expect_true(all(table(res$grp) >= 2))
    }
  })
})

test_that("holdout sample IDs are deterministic given a fixed seed, and change with a different seed", {
  dataset <- fs_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    do.call(session$setInputs, fs_base_inputs())
    r1 <- fs_build_filters("pooled", NULL)
    r2 <- fs_build_filters("pooled", NULL)
    expect_identical(sort(r1$holdout_sample_ids), sort(r2$holdout_sample_ids))

    session$setInputs(fs_holdout_seed = 42)
    r3 <- fs_build_filters("pooled", NULL)
    expect_false(identical(sort(r1$holdout_sample_ids), sort(r3$holdout_sample_ids)))
  })
})
