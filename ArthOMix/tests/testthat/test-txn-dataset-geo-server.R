## Module 1 (Transcriptomics) - Dataset tab, "Fetch from NCBI GEO" path.
##
## mod_dataset.R's geo_fetch_result() calls GEOquery::getGEO(acc, GSEMatrix =
## TRUE), which - even when the target series_matrix.txt.gz is already
## cached on disk - always makes one network call first
## (GEOquery:::getAndParseGSEMatrices() calls getDirListing() over HTTPS to
## enumerate available matrix files, unconditionally, before ever checking
## its local cache; confirmed by reading that internal function's source).
## There is therefore no way to exercise the *actual* getGEO() call fully
## offline without either network access or modifying the production
## fetch path in mod_dataset.R - which the testing brief explicitly
## disallows. Per the brief's two-tier GEO testing requirement, the offline
## tier below instead mocks GEOquery::getGEO() itself at the package
## boundary (testthat::local_mocked_bindings(), restored automatically at
## the end of each test - never touches R/transcriptomics/01_Data/mod_dataset.R) so
## every line of ArthOMix's OWN post-fetch logic (platform-picker,
## probe-to-gene collapsing, guess_col()'s exact/contains/fallback tiers,
## sample-overlap validation, dataset activation) still runs for real
## against a small deterministic ExpressionSet fixture. The live tier
## actually calls the real network path against a real, small, stable GEO
## series and is skipped by default.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))

geo_fixture_eset <- fx_geo_eset()

## ---- Offline, deterministic (mocked GEOquery::getGEO) --------------------

test_that("fetching a valid accession (mocked, single platform) populates the platform mapping UI with the right guessed columns", {
  testthat::local_mocked_bindings(getGEO = function(...) list(GPL_FIXTURE = geo_fixture_eset), .package = "GEOquery")

  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)

    em <- geo_expr_meta()
    expect_false(inherits(em, "error"))
    expect_equal(dim(em$expr), c(30L, 10L))

    ## Single platform -> no platform-choice selector needed.
    expect_null(output$geo_platform_ui)

    ## guess_col()'s tiering (see mod_dataset.R comment): "status" and
    ## "characteristics_ch1" must NOT win the group guess over the real
    ## "disease state:ch1" field, even though "status" is a literal exact
    ## column name here - this fixture is deliberately built to include
    ## both decoys, matching the GSE93272 case the code comment describes.
    html <- output$geo_column_mapping
    expect_equal(fx_selected_value(html, session$ns("geo_map_group")), "disease state:ch1")
    expect_equal(fx_selected_value(html, session$ns("geo_map_sex")), "Sex:ch1")
  })
})

test_that("loading a mocked GEO fetch activates the dataset with source_type = 'geo' and the fetched accession as geo_ids", {
  testthat::local_mocked_bindings(getGEO = function(...) list(GPL_FIXTURE = geo_fixture_eset), .package = "GEOquery")

  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "gse99999")  ## lowercase - acc is upper-cased internally
    session$setInputs(geo_fetch_btn = 1)
    session$setInputs(geo_map_group = "disease state:ch1", geo_map_sex = "Sex:ch1", geo_map_batch = "(none)")
    session$setInputs(geo_load_btn = 1)

    expect_equal(ncol(dataset$expr), 10L)
    expect_equal(dataset$source_type, "geo")
    expect_false(dataset$is_bundled_reference)
    expect_equal(dataset$geo_ids, "GSE99999")
    expect_true(grepl("^NCBI GEO: GSE99999", dataset$source))
  })
})

test_that("loading a GEO dataset clears any stale declared_data_type left over from a previous upload", {
  testthat::local_mocked_bindings(getGEO = function(...) list(GPL_FIXTURE = geo_fixture_eset), .package = "GEOquery")

  dataset <- shiny::reactiveValues(declared_data_type = "raw")
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)
    session$setInputs(geo_map_group = "disease state:ch1", geo_map_sex = "Sex:ch1", geo_map_batch = "(none)")
    session$setInputs(geo_load_btn = 1)
    expect_true(is.na(dataset$declared_data_type))
  })
})

test_that("a multi-platform series (mocked) shows a platform picker, and switching platforms switches the parsed expr/meta", {
  eset_b <- geo_fixture_eset
  Biobase::exprs(eset_b) <- Biobase::exprs(eset_b) + 100  ## distinguishable from platform A
  Biobase::annotation(eset_b) <- "GPL_FIXTURE_B"
  testthat::local_mocked_bindings(
    getGEO = function(...) list(GPL_FIXTURE = geo_fixture_eset, GPL_FIXTURE_B = eset_b),
    .package = "GEOquery"
  )

  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)
    expect_false(is.null(output$geo_platform_ui))

    session$setInputs(geo_platform_choice = "GPL_FIXTURE")
    em_a <- geo_expr_meta()
    session$setInputs(geo_platform_choice = "GPL_FIXTURE_B")
    em_b <- geo_expr_meta()
    expect_false(isTRUE(all.equal(as.numeric(em_a$expr[1, 1]), as.numeric(em_b$expr[1, 1]))))
  })
})

test_that("geo_fetch_status shows no error banner for a multi-platform series before a platform is picked", {
  ## Regression test: geo_eset()'s req(input$geo_platform_choice) throws a
  ## shiny.silent.error, which DOES inherit from "error" - a tryCatch that
  ## isn't careful to treat that as "still waiting" instead of a real
  ## failure re-surfaces it through geo_fetch_status's generic
  ## inherits(em, "error") branch as a visible (near-blank) warning banner,
  ## right after every multi-platform fetch, before the user has had a
  ## chance to pick a platform.
  eset_b <- geo_fixture_eset
  Biobase::annotation(eset_b) <- "GPL_FIXTURE_B"
  testthat::local_mocked_bindings(
    getGEO = function(...) list(GPL_FIXTURE = geo_fixture_eset, GPL_FIXTURE_B = eset_b),
    .package = "GEOquery"
  )

  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)
    expect_false(is.null(output$geo_platform_ui))
    ## No geo_platform_choice set yet. Correct behavior is either (a) Shiny's
    ## own silent req()-halt (accessing the output re-throws its stored
    ## shiny.silent.error - normal testServer behavior for an unrendered,
    ## still-waiting output), or (b) real rendered content with no error
    ## banner in it. Before the fix, this used to render REAL content - a
    ## visible triangle-exclamation warning div - so branch (b) below is
    ## what would have caught the regression.
    status_result <- tryCatch(output$geo_fetch_status, error = function(e) e)
    if (inherits(status_result, "error")) {
      expect_true(inherits(status_result, "shiny.silent.error"))
    } else {
      expect_false(grepl("triangle-exclamation|Could not fetch", fx_html_text(status_result)))
    }
  })
})

test_that("fewer than 4 overlapping samples between a (mocked) fetched GEO series and its own pData is rejected", {
  ## Not reachable through the normal fetch path (pData rownames always equal
  ## exprs() colnames for a real ExpressionSet), but geo_load_btn's own
  ## common-sample validate() is still exercised directly here for
  ## completeness/regression safety, the same way mod_dataset.R's upload
  ## path's identical check is covered in test-dataset-upload.R.
  degenerate_eset <- geo_fixture_eset[, 1:3]
  testthat::local_mocked_bindings(getGEO = function(...) list(GPL_FIXTURE = degenerate_eset), .package = "GEOquery")

  dataset <- shiny::reactiveValues(expr = "sentinel_expr")
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)
    ## Only 3 samples total in this fixture - always below the 4-ID minimum
    ## regardless of mapping, so this exercises the same "Fewer than 4"
    ## validate() guard as the upload path.
    session$setInputs(geo_map_group = "disease state:ch1", geo_map_sex = "(none)", geo_map_batch = "(none)")
    session$setInputs(geo_load_btn = 1)
    expect_identical(dataset$expr, "sentinel_expr")
    expect_true(grepl("Fewer than 4 sample IDs", fx_html_text(output$load_message)))
  })
})

test_that("an invalid GEO accession format is rejected before any network call is attempted", {
  called <- FALSE
  testthat::local_mocked_bindings(
    getGEO = function(...) { called <<- TRUE; stop("should not be called") },
    .package = "GEOquery"
  )
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "not-a-real-accession")
    session$setInputs(geo_fetch_btn = 1)
    expect_true(grepl("Enter a valid GEO Series accession", fx_html_text(output$geo_fetch_status)))
  })
  expect_false(called)
})

test_that("a series with an empty expression matrix (RNA-seq raw-counts-only deposit) is rejected with a guiding error, not loaded as empty", {
  empty_eset <- geo_fixture_eset[0, ]
  testthat::local_mocked_bindings(getGEO = function(...) list(GPL_FIXTURE = empty_eset), .package = "GEOquery")

  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)
    expect_true(grepl("no expression matrix", fx_html_text(output$geo_fetch_status)))
  })
})

## ---- Live network integration (explicitly opt-in) ------------------------

test_that("live: fetching a real, small, stable GEO series succeeds end to end", {
  skip_if_not(
    identical(Sys.getenv("ARTHOMIX_RUN_LIVE_GEO_TESTS"), "1"),
    "Set ARTHOMIX_RUN_LIVE_GEO_TESTS=1 to run this network-dependent test."
  )
  skip_if_offline <- function() {
    ok <- tryCatch({
      con <- url("https://ftp.ncbi.nlm.nih.gov", open = "rb")
      close(con)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) skip("No network access to NCBI right now.")
  }
  skip_if_offline()

  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    ## GSE15573 is one of this app's own bundled validation cohorts (small,
    ## stable, already vetted by this project) - reused here instead of an
    ## arbitrary external accession.
    session$setInputs(geo_accession = "GSE15573")
    session$setInputs(geo_fetch_btn = 1)
    em <- geo_expr_meta()
    expect_false(inherits(em, "error"))
    expect_gt(nrow(em$expr), 0)
    expect_gt(ncol(em$expr), 0)
  })
})
