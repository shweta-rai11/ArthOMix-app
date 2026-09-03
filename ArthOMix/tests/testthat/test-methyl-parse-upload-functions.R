## Module 2 (Methylomics) - Dataset tab's upload-parsing helpers
## (parse_upload.R). Every parser here promises "never throws, returns
## list(ok=FALSE, error=...)" - both the happy path and that fail-soft

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "parse_upload.R"))

meth_write_matrix_csv <- function(mat, path) {
  df <- data.frame(probe = rownames(mat), mat, check.names = FALSE)
  data.table::fwrite(df, path)
  invisible(path)
}

test_that("methyl_parse_matrix() parses a well-formed probe x sample beta matrix", {
  set.seed(200)
  mat <- matrix(runif(30, 0, 1), 10, 3, dimnames = list(paste0("cg", 10000000 + 1:10), paste0("S", 1:3)))
  path <- withr::local_tempfile(fileext = ".csv")
  meth_write_matrix_csv(mat, path)
  res <- methyl_parse_matrix(path, "test.csv")
  expect_true(res$ok)
  expect_equal(dim(res$mat), c(10L, 3L))
})

test_that("methyl_parse_matrix() rejects duplicate probe IDs", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("probe,S1,S2", "cg001,0.1,0.2", "cg001,0.3,0.4", "cg002,0.5,0.6"), path)
  res <- methyl_parse_matrix(path, "dup.csv")
  expect_false(res$ok)
  expect_true(grepl("duplicated probe ID", res$error))
})

test_that("methyl_parse_matrix() rejects a file with fewer than 2 columns or 0 rows", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines("probe", path)
  res <- methyl_parse_matrix(path, "onecol.csv")
  expect_false(res$ok)
})

test_that("methyl_parse_matrix() never throws on a genuinely unparseable/malformed file", {
  path <- normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "malformed_expr.csv"), mustWork = TRUE)
  res <- tryCatch(methyl_parse_matrix(path, "malformed_expr.csv"), error = function(e) e)
  expect_false(inherits(res, "error"))
  if (isTRUE(res$ok)) expect_true(is.matrix(res$mat)) else expect_true(nzchar(res$error))
})

test_that("methyl_parse_sample_sheet() parses a well-formed sheet and rejects an empty one", {
  path <- withr::local_tempfile(fileext = ".csv")
  write.csv(data.frame(sample = c("S1", "S2"), sex = c("F", "M")), path, row.names = FALSE)
  res <- methyl_parse_sample_sheet(path, "sheet.csv")
  expect_true(res$ok)
  expect_equal(nrow(res$df), 2L)

  empty_path <- normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "empty_expr.csv"), mustWork = TRUE)
  res2 <- methyl_parse_sample_sheet(empty_path, "empty.csv")
  expect_false(res2$ok)
})

test_that("methyl_detect_orientation() correctly identifies a transposed (samples x probes) matrix", {
  mat_correct <- matrix(runif(20), 10, 2, dimnames = list(paste0("cg", 1:10), paste0("S", 1:2)))
  expect_false(methyl_detect_orientation(mat_correct)$transposed)

  mat_transposed <- t(mat_correct)
  expect_true(methyl_detect_orientation(mat_transposed)$transposed)
})

test_that("methyl_validate_matrix_upload() auto-transposes a backwards matrix and preserves values", {
  mat <- matrix(c(0.1, 0.5, 0.9, 0.2), 2, 2, dimnames = list(c("cg001", "cg002"), c("S1", "S2")))
  backwards <- t(mat)
  res <- methyl_validate_matrix_upload(backwards, declared_scale = "beta")
  expect_true(res$ok)
  expect_equal(dim(res$mat), dim(mat))
  expect_true(grepl("transposed automatically", res$note))
})

test_that("methyl_validate_matrix_upload() rejects a matrix declared as beta values that is mostly out of [0,1]", {
  mat <- matrix(rnorm(40, mean = 0, sd = 2), 10, 4, dimnames = list(paste0("cg", 1:10), paste0("S", 1:4)))
  res <- methyl_validate_matrix_upload(mat, declared_scale = "beta")
  expect_false(res$ok)
  expect_true(grepl("Beta values", res$error))
})

test_that("methyl_validate_matrix_upload() accepts a genuine beta matrix within [0,1]", {
  set.seed(201)
  mat <- matrix(runif(40, 0, 1), 10, 4, dimnames = list(paste0("cg", 1:10), paste0("S", 1:4)))
  res <- methyl_validate_matrix_upload(mat, declared_scale = "beta")
  expect_true(res$ok)
})

test_that("methyl_validate_matrix_upload() warns (does not block) when M-values look suspiciously beta-like", {
  set.seed(202)
  mat <- matrix(runif(40, 0.3, 0.7), 10, 4, dimnames = list(paste0("cg", 1:10), paste0("S", 1:4)))
  res <- methyl_validate_matrix_upload(mat, declared_scale = "m")
  expect_true(res$ok)
  expect_true(grepl("could already be beta values", res$note))
})

test_that("methyl_validate_matrix_upload() rejects an all-non-finite matrix", {
  mat <- matrix(NA_real_, 5, 3, dimnames = list(paste0("cg", 1:5), paste0("S", 1:3)))
  res <- methyl_validate_matrix_upload(mat, declared_scale = "beta")
  expect_false(res$ok)
  expect_true(grepl("No finite numeric values", res$error))
})

test_that("methyl_parse_probe_list() reads one probe ID per line and dedupes", {
  path <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("cg001", "cg002", "cg001", "cg003"), path)
  res <- methyl_parse_probe_list(path, "probes.txt")
  expect_true(res$ok)
  expect_setequal(res$ids, c("cg001", "cg002", "cg003"))
})

test_that("methyl_parse_probe_list() also accepts the first column of a CSV", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("cg001,extra1", "cg002,extra2"), path)
  res <- methyl_parse_probe_list(path, "probes.csv")
  expect_true(res$ok)
  expect_setequal(res$ids, c("cg001", "cg002"))
})

test_that("methyl_parse_probe_list() reports failure for an empty file rather than an empty success", {
  path <- withr::local_tempfile(fileext = ".txt")
  writeLines(character(0), path)
  res <- methyl_parse_probe_list(path, "empty.txt")
  expect_false(res$ok)
})

test_that("methyl_read_idat() rejects an upload with no .idat files, and one missing a Grn/Red pair", {
  files_no_idat <- data.frame(name = c("a.csv", "b.txt"), datapath = c("/tmp/a.csv", "/tmp/b.txt"), stringsAsFactors = FALSE)
  res1 <- methyl_read_idat(files_no_idat)
  expect_false(res1$ok)
  expect_true(grepl("No .idat files", res1$error))

  files_grn_only <- data.frame(name = c("S1_Grn.idat"), datapath = c("/tmp/S1_Grn.idat"), stringsAsFactors = FALSE)
  res2 <- methyl_read_idat(files_grn_only)
  expect_false(res2$ok)
  expect_true(grepl("_Grn.idat and a _Red.idat", res2$error))
})
