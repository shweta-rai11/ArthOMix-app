## Module 1 (Transcriptomics) - Dataset tab "Reserve validation samples" (sealed hold-out taken BEFORE any
## analysis), via a pure helper plus testServer(): stratified split, removal from the active dataset,
## bundled shortcuts disabled, release restores, and re-sealing when another module replaces dataset$expr.
## Added 2026-09-05 (transcriptomics audit, finding 1 long-term fix).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))

test_that("tx_reserve_validation_ids() splits by group and sex, is deterministic, and refuses tiny or one-group data", {
  meta <- data.frame(sample = paste0("S", 1:80),
                     group = rep(c("HC", "RA"), each = 40),
                     sex = rep(c("F", "M", "F", "M"), each = 20), stringsAsFactors = FALSE)
  ids <- tx_reserve_validation_ids(meta, 0.3, seed = 7, stratify_sex = TRUE)
  expect_true(all(ids %in% meta$sample))
  expect_equal(length(ids), 24)
  sub <- meta[meta$sample %in% ids, ]
  expect_equal(as.integer(table(sub$group)), c(12, 12))
  expect_equal(as.integer(table(sub$sex)), c(12, 12))
  expect_identical(ids, tx_reserve_validation_ids(meta, 0.3, seed = 7, stratify_sex = TRUE))
  expect_false(identical(ids, tx_reserve_validation_ids(meta, 0.3, seed = 8, stratify_sex = TRUE)))

  expect_error(tx_reserve_validation_ids(meta[1:10, ], 0.3), "At least 20 samples")
  expect_error(tx_reserve_validation_ids(transform(meta, group = "RA"), 0.3), "two distinct values")
  small <- data.frame(sample = paste0("S", 1:24), group = c(rep("HC", 4), rep("RA", 20)), stringsAsFactors = FALSE)
  expect_error(tx_reserve_validation_ids(small, 0.5, stratify_sex = FALSE), "fewer than 4 samples")
})

test_that("Reserve removes a sealed set from the active bundled dataset, disables bundled shortcuts, tags the source, and Release restores everything", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    arthomix_provenance_clear(session)
    session$setInputs(preloaded_choice = "__default_merged__")
    session$setInputs(load_preloaded_btn = 1)
    n_all <- ncol(dataset$expr)
    all_ids <- colnames(dataset$expr)
    expect_true(dataset$is_bundled_reference)
    expect_equal(dataset$reserved_ids, character(0))

    session$setInputs(reserve_frac_pct = 30, reserve_seed = 1234, reserve_stratify_sex = TRUE, reserve_btn = 0)
    session$setInputs(reserve_btn = 1)
    ids <- dataset$reserved_ids
    expect_gt(length(ids), 0)
    expect_equal(length(ids) + ncol(dataset$expr), n_all)
    expect_length(intersect(ids, colnames(dataset$expr)), 0)
    expect_setequal(c(ids, colnames(dataset$expr)), all_ids)
    expect_identical(colnames(dataset$reserved_expr), ids)
    expect_identical(dataset$reserved_meta$sample, ids)
    expect_identical(dataset$meta$sample, colnames(dataset$expr))
    expect_true(all(table(dataset$reserved_meta$group) >= 4))
    expect_true(all(table(dataset$meta$group) >= 4))
    expect_false(dataset$is_bundled_reference)
    expect_match(dataset$source, sprintf("\\[%d validation samples reserved\\]$", length(ids)))
    expect_true(isTRUE(dataset$reserved_info$was_bundled_reference))
    recs <- arthomix_provenance_records(session)
    expect_true("mod_dataset_reserve_validation" %in% vapply(recs, `[[`, character(1), "module"))

    ## a second Reserve click is refused, nothing changes
    session$setInputs(reserve_btn = 2)
    expect_identical(dataset$reserved_ids, ids)

    session$setInputs(release_btn = 0)
    session$setInputs(release_btn = 1)
    expect_equal(dataset$reserved_ids, character(0))
    expect_null(dataset$reserved_expr)
    expect_equal(ncol(dataset$expr), n_all)
    expect_setequal(colnames(dataset$expr), all_ids)
    expect_identical(dataset$meta$sample, colnames(dataset$expr))
    expect_true(dataset$is_bundled_reference)
    expect_false(grepl("validation samples reserved", dataset$source))
  })
})

test_that("reserved samples are sealed again when another module replaces the active matrix with one that contains them, and a new dataset load clears the reservation", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(preloaded_choice = "__default_merged__")
    session$setInputs(load_preloaded_btn = 1)
    full_expr <- dataset$expr; full_meta <- dataset$meta
    session$setInputs(reserve_frac_pct = 30, reserve_seed = 1234, reserve_stratify_sex = TRUE, reserve_btn = 0)
    session$setInputs(reserve_btn = 1)
    ids <- dataset$reserved_ids
    n_disc <- ncol(dataset$expr)

    ## simulate Preprocessing "activate" rebuilding the whole cohort (shifted values so we can see the re-seal used them)
    dataset$expr <- full_expr + 1
    dataset$meta <- full_meta
    dataset$source <- "Preprocessed dataset (rebuilt from raw)"
    session$flushReact()
    expect_equal(ncol(dataset$expr), n_disc)
    expect_length(intersect(ids, colnames(dataset$expr)), 0)
    expect_identical(colnames(dataset$reserved_expr), ids)
    expect_equal(unname(dataset$reserved_expr[1, ids[1]]), unname(full_expr[1, ids[1]]) + 1)
    expect_equal(dataset$reserved_info$resealed, 1L)
    expect_match(dataset$source, "validation samples reserved\\]$")

    ## loading another dataset drops the reservation entirely
    session$setInputs(preloaded_choice = "GSE15573")
    session$setInputs(load_preloaded_btn = 2)
    expect_equal(dataset$reserved_ids, character(0))
    expect_null(dataset$reserved_expr)
    expect_false(grepl("validation samples reserved", dataset$source))
  })
})
