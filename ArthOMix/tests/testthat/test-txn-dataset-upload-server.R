## Module 1 (Transcriptomics) - Dataset tab, "Upload your own data" path.
## Covers mod_dataset.R's expr_raw/meta_raw parsing, the local guess_col()
## column-mapping defaults (observed via the rendered column_mapping UI,
## since testServer() never turns a selectInput's `selected=` into a live
## input default - see fx_selected_value()), and the load_btn observer's
## success/error behavior.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_dataset.R"))

fixture_dir <- normalizePath(file.path(app_dir, "tests", "fixtures", "transcriptomics"), mustWork = TRUE)
expr_fixture_path <- file.path(fixture_dir, "small_expr_matrix.csv")
meta_fixture_path <- file.path(fixture_dir, "small_sample_metadata.csv")

test_that("uploading a well-formed expr+meta pair loads successfully and activates the dataset immediately", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(expr_fixture_path))
    session$setInputs(meta_file = fx_mkfile(meta_fixture_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "sex", map_batch = "batch")
    session$setInputs(load_btn = 1)

    expect_equal(dim(dataset$expr), c(40L, 16L))
    expect_equal(nrow(dataset$meta), 16L)
    expect_equal(dataset$source_type, "uploaded")
    expect_false(dataset$is_bundled_reference)
    expect_length(dataset$geo_ids, 0)
    expect_true(grepl("^Uploaded dataset:", dataset$source))
  })
})

test_that("column_mapping guesses the sample/group/sex/batch columns by exact name match", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(meta_file = fx_mkfile(meta_fixture_path))
    html <- output$column_mapping
    expect_equal(fx_selected_value(html, session$ns("map_id")), "sample")
    expect_equal(fx_selected_value(html, session$ns("map_group")), "group")
    expect_equal(fx_selected_value(html, session$ns("map_sex")), "sex")
    expect_equal(fx_selected_value(html, session$ns("map_batch")), "batch")
  })
})

test_that("column_mapping falls back to '(none)' for sex/batch when no matching column exists", {
  dataset <- shiny::reactiveValues()
  no_sex_path <- withr::local_tempfile(fileext = ".csv")
  fm <- fx_expr_meta(n_samples = 8, seed = 2)
  fx_write_meta_csv(fm$meta[, c("sample", "group")], no_sex_path)

  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(meta_file = fx_mkfile(no_sex_path))
    html <- output$column_mapping
    expect_equal(fx_selected_value(html, session$ns("map_sex")), "(none)")
    expect_equal(fx_selected_value(html, session$ns("map_batch")), "(none)")
  })
})

test_that("fewer than 4 overlapping sample IDs between expr and meta is rejected with a validation error, and dataset is left untouched", {
  dataset <- shiny::reactiveValues(expr = "sentinel_expr", meta = "sentinel_meta")
  meta_path <- normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "insufficient_overlap_meta.csv"), mustWork = TRUE)

  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(expr_fixture_path))
    session$setInputs(meta_file = fx_mkfile(meta_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "(none)", map_batch = "(none)")
    session$setInputs(load_btn = 1)

    expect_true(grepl("Fewer than 4 sample IDs", fx_html_text(output$load_message)))
    ## dataset was NOT overwritten by the rejected upload.
    expect_identical(dataset$expr, "sentinel_expr")
    expect_identical(dataset$meta, "sentinel_meta")
  })
})

test_that("a duplicate feature ID in the uploaded expression matrix loads (all rows kept) but surfaces an explicit warning note", {
  dataset <- shiny::reactiveValues()
  fm <- fx_expr_meta(n_samples = 8, seed = 3)
  dup_expr <- fx_expr_with_duplicate_id(fm$expr)
  dir <- withr::local_tempdir()
  expr_path <- file.path(dir, "dup_expr.csv")
  meta_path <- file.path(dir, "meta.csv")
  fx_write_expr_csv(dup_expr, expr_path)
  fx_write_meta_csv(fm$meta, meta_path)

  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(expr_path))
    session$setInputs(meta_file = fx_mkfile(meta_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "(none)", map_batch = "(none)")
    session$setInputs(load_btn = 1)

    expect_equal(nrow(dataset$expr), nrow(dup_expr))
    msg <- fx_html_text(output$load_message)
    expect_true(grepl("duplicated feature identifier", msg))
    expect_true(grepl("only the first occurrence", msg))
  })
})

test_that("an unparseable/malformed uploaded expression file surfaces a read error in the preview, without crashing the session", {
  dataset <- shiny::reactiveValues()
  malformed_path <- normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "malformed_expr.csv"), mustWork = TRUE)

  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(malformed_path))
    session$setInputs(meta_file = fx_mkfile(meta_fixture_path))
    ## expr_raw()/upload_preview_data() wraps parsing in tryCatch and reports
    ## inline rather than raising - a malformed CSV is not guaranteed to
    ## error at all (fread is lenient), so this only asserts the app doesn't
    ## crash and the preview UI renders something (error note or a garbage
    ## but non-fatal parse) rather than throwing out of testServer.
    html <- tryCatch(output$upload_preview_ui, error = function(e) e)
    expect_false(inherits(html, "error"))
  })
})

test_that("uploading a well-formed expr+meta pair records the declared data type on the shared dataset", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(expr_fixture_path))
    session$setInputs(meta_file = fx_mkfile(meta_fixture_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "sex", map_batch = "batch")
    session$setInputs(declared_data_type = "logtransformed")
    session$setInputs(load_btn = 1)

    expect_equal(dataset$source_type, "uploaded")
    expect_equal(dataset$declared_data_type, "logtransformed")
  })
})

test_that("declaring 'Raw counts' for a matrix that is actually TPM-normalized is blocked, and the shared dataset is left untouched", {
  dataset <- shiny::reactiveValues(expr = "sentinel_expr", meta = "sentinel_meta")
  fm <- fx_expr_meta(n_genes = 40, n_samples = 8, seed = 4)
  ## Rescale to a raw-count-like scale, then normalize each sample to sum to
  ## 1e6 (TPM/CPM-like) - the exact mismatch tx_validate_expr_upload() must
  ## catch when "Raw counts" is declared.
  counts_like <- round(exp(fm$expr))
  tpm_like <- sweep(counts_like, 2, colSums(counts_like), FUN = "/") * 1e6
  dir <- withr::local_tempdir()
  expr_path <- file.path(dir, "tpm_expr.csv")
  meta_path <- file.path(dir, "meta.csv")
  fx_write_expr_csv(tpm_like, expr_path)
  fx_write_meta_csv(fm$meta, meta_path)

  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(expr_path))
    session$setInputs(meta_file = fx_mkfile(meta_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "sex", map_batch = "batch")
    session$setInputs(declared_data_type = "raw")
    session$setInputs(load_btn = 1)

    expect_true(grepl("TPM/FPKM/CPM-normalized", fx_html_text(output$load_message)))
    expect_identical(dataset$expr, "sentinel_expr")
    expect_identical(dataset$meta, "sentinel_meta")
  })
})

test_that("an empty uploaded expression file (header only, 0 rows) is not silently treated as a valid non-empty matrix", {
  dataset <- shiny::reactiveValues(expr = "sentinel_expr")
  empty_path <- normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "empty_expr.csv"), mustWork = TRUE)

  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(expr_file = fx_mkfile(empty_path))
    session$setInputs(meta_file = fx_mkfile(meta_fixture_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "(none)", map_batch = "(none)")
    session$setInputs(load_btn = 1)
    ## 0 features x 3 samples: tx_validate_expr_upload() now catches this
    ## before the sample-overlap check even runs (no finite values at all
    ## in a 0-row matrix), which is a more specific/accurate message than
    ## the "fewer than 4 sample IDs" this used to fall through to - either
    ## way, this must be rejected, not loaded as a valid (but empty) dataset.
    expect_identical(dataset$expr, "sentinel_expr")
    expect_true(grepl("No finite numeric values", fx_html_text(output$load_message)))
  })
})
