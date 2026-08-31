## Module 2 (Methylomics) - Dataset tab: preloaded/upload/GEO paths, via
## testServer(). The preloaded path normally loads its ~2.1GB live matrix
## asynchronously (ExtendedTask/future_promise) - pre-warming
## .arthomix_cache[["meth_default_matrix"]] by calling
## load_default_meth_matrix() once before the click routes it through the
## same function's own synchronous "already cached" branch instead, so the
## real finish_preloaded_load() path is exercised without needing to drive
## an actual async promise resolution through testServer().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "parse_upload.R"))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_dataset.R"))

test_that("loading the preloaded whole-blood dataset (pre-warmed cache) populates methyl_dataset with a real beta matrix", {
  invisible(load_default_meth_matrix())  ## pre-warm .arthomix_cache
  methyl_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mds", methyl_dataset = methyl_dataset), {
    session$setInputs(preloaded_choice = "gse42861_wholeblood")
    session$setInputs(load_preloaded_btn = 1)

    expect_false(is.null(methyl_dataset$beta))
    expect_gt(nrow(methyl_dataset$beta), 0)
    expect_equal(methyl_dataset$input_scale, "beta")
    expect_equal(methyl_dataset$array_type, "450K")
    expect_equal(methyl_dataset$source_type, "preloaded")
    expect_true(methyl_dataset$preloaded)
  })
})

test_that("uploading a matrix + sample sheet loads methyl_dataset with correctly mapped sample/group/sex columns", {
  set.seed(230)
  mat <- matrix(runif(200, 0, 1), 20, 10, dimnames = list(paste0("cg", 10000000 + 1:20), paste0("S", 1:10)))
  sheet <- data.frame(sample = colnames(mat), group = rep(c("HC", "RA"), 5), sex = rep(c("F", "M"), 5), stringsAsFactors = FALSE)
  dir <- withr::local_tempdir()
  mat_path <- file.path(dir, "beta.csv")
  sheet_path <- file.path(dir, "sheet.csv")
  data.table::fwrite(data.frame(probe = rownames(mat), mat, check.names = FALSE), mat_path)
  write.csv(sheet, sheet_path, row.names = FALSE)

  methyl_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mds", methyl_dataset = methyl_dataset), {
    session$setInputs(upload_format = "matrix", array_type = "EPIC", input_scale = "beta")
    session$setInputs(matrix_file = fx_mkfile(mat_path))
    session$setInputs(sheet_file = fx_mkfile(sheet_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "sex", map_batch = "(none)")
    session$setInputs(load_btn = 1)

    expect_equal(dim(methyl_dataset$beta), c(20L, 10L))
    expect_equal(methyl_dataset$source_type, "upload")
    expect_false(methyl_dataset$preloaded)
    expect_setequal(methyl_dataset$sample_sheet$sex, c("F", "M"))
  })
})

test_that("uploading a matrix with NO sample sheet still loads successfully (sheet is optional for methylomics)", {
  set.seed(231)
  mat <- matrix(runif(80, 0, 1), 20, 4, dimnames = list(paste0("cg", 10000000 + 1:20), paste0("S", 1:4)))
  path <- withr::local_tempfile(fileext = ".csv")
  data.table::fwrite(data.frame(probe = rownames(mat), mat, check.names = FALSE), path)

  methyl_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mds", methyl_dataset = methyl_dataset), {
    session$setInputs(upload_format = "matrix", array_type = "EPIC", input_scale = "beta")
    session$setInputs(matrix_file = fx_mkfile(path))
    session$setInputs(load_btn = 1)

    expect_equal(dim(methyl_dataset$beta), c(20L, 4L))
    expect_null(methyl_dataset$sample_sheet)
  })
})

test_that("an upload declared as beta values but mostly outside [0,1] is rejected", {
  mat <- matrix(rnorm(80, mean = 0, sd = 3), 20, 4, dimnames = list(paste0("cg", 10000000 + 1:20), paste0("S", 1:4)))
  path <- withr::local_tempfile(fileext = ".csv")
  data.table::fwrite(data.frame(probe = rownames(mat), mat, check.names = FALSE), path)

  methyl_dataset <- shiny::reactiveValues(beta = "sentinel")
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mds", methyl_dataset = methyl_dataset), {
    session$setInputs(upload_format = "matrix", array_type = "EPIC", input_scale = "beta")
    session$setInputs(matrix_file = fx_mkfile(path))
    session$setInputs(load_btn = 1)

    expect_identical(methyl_dataset$beta, "sentinel")
    expect_true(grepl("Beta values", fx_html_text(output$load_message)))
  })
})

test_that("fewer than 4 matched sample IDs between the matrix and sample sheet is rejected", {
  mat <- matrix(runif(80, 0, 1), 20, 4, dimnames = list(paste0("cg", 10000000 + 1:20), paste0("S", 1:4)))
  sheet <- data.frame(sample = c("ZZ1", "ZZ2", "ZZ3"), group = c("HC", "RA", "HC"), stringsAsFactors = FALSE)
  dir <- withr::local_tempdir()
  mat_path <- file.path(dir, "beta.csv"); sheet_path <- file.path(dir, "sheet.csv")
  data.table::fwrite(data.frame(probe = rownames(mat), mat, check.names = FALSE), mat_path)
  write.csv(sheet, sheet_path, row.names = FALSE)

  methyl_dataset <- shiny::reactiveValues(beta = "sentinel")
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mds", methyl_dataset = methyl_dataset), {
    session$setInputs(upload_format = "matrix", array_type = "EPIC", input_scale = "beta")
    session$setInputs(matrix_file = fx_mkfile(mat_path))
    session$setInputs(sheet_file = fx_mkfile(sheet_path))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "(none)", map_batch = "(none)")
    session$setInputs(load_btn = 1)

    expect_identical(methyl_dataset$beta, "sentinel")
    expect_true(grepl("Fewer than 4 sample IDs", fx_html_text(output$load_message)))
  })
})

## ---- GEO fetch (offline, mocked GEOquery::getGEO) --------------------------

test_that("fetching a GEO series (mocked, recognized 450K platform) loads with the correct array_type via MX_METHYLATION_GPL", {
  set.seed(232)
  n_probes <- 20; n_samples <- 6
  mat <- matrix(runif(n_probes * n_samples, 0, 1), n_probes, n_samples,
                 dimnames = list(paste0("cg", 10000000 + 1:n_probes), paste0("GSM", 1:n_samples)))
  pdat <- data.frame(geo_accession = colnames(mat), `disease state:ch1` = rep(c("RA", "HC"), 3),
                       check.names = FALSE, row.names = colnames(mat))
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL13534"  ## maps to "450K" per MX_METHYLATION_GPL

  testthat::local_mocked_bindings(getGEO = function(...) list(GPL13534 = eset), .package = "GEOquery")

  methyl_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mds", methyl_dataset = methyl_dataset), {
    session$setInputs(geo_accession = "GSE99999")
    session$setInputs(geo_fetch_btn = 1)
    em <- geo_expr_meta()
    expect_false(inherits(em, "error"))
    expect_equal(em$array_type, "450K")
  })
})
