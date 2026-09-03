## Provenance manifest helper (R/provenance.R): arthomix_provenance_record()
## and its JSON-serialization/download-handler support. This is the shared
## helper wired into mod_dge.R, mod_methyl_dmp.R, and (as an extension of

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "provenance.R"))

test_that("arthomix_provenance_record()'s checksum is deterministic for identical input and changes with different input", {
  input_a <- list(mat_dim = c(100, 10), group = c("HC", "RA", "HC", "RA"))
  r1 <- arthomix_provenance_record(module = "test_mod", checksum_input = input_a)
  r2 <- arthomix_provenance_record(module = "test_mod", checksum_input = input_a)
  expect_identical(r1$checksum, r2$checksum)

  input_b <- list(mat_dim = c(100, 11), group = c("HC", "RA", "HC", "RA"))
  r3 <- arthomix_provenance_record(module = "test_mod", checksum_input = input_b)
  expect_false(isTRUE(all.equal(r1$checksum, r3$checksum)))

  expect_type(r1$checksum, "character")
  expect_true(nzchar(r1$checksum))
})

test_that("arthomix_provenance_record()'s checksum also changes when only params differ but checksum_input is identical", {
  input_a <- list(mat_dim = c(50, 5))
  r1 <- arthomix_provenance_record(module = "m", checksum_input = input_a, params = list(x = 1))
  r2 <- arthomix_provenance_record(module = "m", checksum_input = input_a, params = list(x = 2))
  expect_identical(r1$checksum, r2$checksum)
  expect_false(identical(r1$params, r2$params))
})

test_that("arthomix_provenance_record() always returns every documented top-level field", {
  r <- arthomix_provenance_record(
    module = "mod_dge",
    checksum_input = list(a = 1),
    params = list(method = "limma"),
    seed = 42,
    packages = c("digest"),
    extra = list(note = "x")
  )
  expect_named(r, c("schema_version", "module", "run_at", "checksum", "params", "seed", "software", "extra"),
               ignore.order = FALSE)
  expect_identical(r$schema_version, "1.0")
  expect_identical(r$module, "mod_dge")
  expect_s3_class(r$run_at, "POSIXct")
  expect_identical(r$seed, 42)
  expect_named(r$software, c("r_version", "packages"))
  expect_identical(r$software$r_version, as.character(getRversion()))
  expect_true("digest" %in% names(r$software$packages))
})

test_that("arthomix_provenance_record() defaults produce an empty-but-present params/extra and NULL seed", {
  r <- arthomix_provenance_record(module = "m", checksum_input = list(a = 1))
  expect_identical(r$params, list())
  expect_identical(r$extra, list())
  expect_null(r$seed)
  expect_identical(r$software$packages, stats::setNames(list(), character(0)))
})

test_that("an uninstalled/nonexistent package name degrades to a note instead of erroring the whole record", {
  r <- expect_no_error(arthomix_provenance_record(
    module = "m", checksum_input = list(a = 1),
    packages = c("digest", "this_package_definitely_does_not_exist_12345")
  ))
  expect_identical(r$software$packages[["this_package_definitely_does_not_exist_12345"]], "not installed")
  expect_true(grepl("^[0-9]+\\.[0-9]+", r$software$packages[["digest"]]))
})

test_that("duplicate package names are only resolved/listed once", {
  r <- arthomix_provenance_record(module = "m", checksum_input = list(a = 1), packages = c("digest", "digest"))
  expect_identical(names(r$software$packages), "digest")
})

test_that("arthomix_provenance_json_safe() formats POSIXct/Date/factor fields instead of leaving them opaque for jsonlite", {
  x <- list(run_at = as.POSIXct("2026-01-01 12:00:00", tz = "UTC"), grp = factor(c("a", "b")), n = 5L, note = NULL)
  safe <- arthomix_provenance_json_safe(x)
  expect_type(safe$run_at, "character")
  expect_identical(safe$grp, c("a", "b"))
  expect_identical(safe$n, 5L)
})

test_that("a full provenance record round-trips through jsonlite::toJSON()/fromJSON() without erroring", {
  r <- arthomix_provenance_record(
    module = "mod_dge", checksum_input = list(dims = c(20, 8), grp = c("A", "B")),
    params = list(method = "limma", padj_cut = 0.05, lfc_cut = 0.1),
    seed = NULL, packages = c("limma", "digest"),
    extra = list(warnings = c("low sample size"))
  )
  safe <- arthomix_provenance_json_safe(r)
  json <- expect_no_error(jsonlite::toJSON(safe, pretty = TRUE, auto_unbox = TRUE, force = TRUE))
  parsed <- expect_no_error(jsonlite::fromJSON(json))
  expect_identical(parsed$schema_version, "1.0")
  expect_identical(parsed$module, "mod_dge")
  expect_identical(parsed$checksum, r$checksum)
})

provenance_handler_content_fn <- function(handler) {
  environment(environment(handler)$renderFunc)$content
}

test_that("arthomix_provenance_download_handler() returns a shiny download handler wrapping record_fn()", {
  record_fn <- function() arthomix_provenance_record(module = "m", checksum_input = list(a = 1), seed = 7)
  handler <- arthomix_provenance_download_handler(record_fn, "test_record")
  expect_true(inherits(handler, "shiny.render.function"))

  content_fn <- provenance_handler_content_fn(handler)
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  content_fn(tmp)
  expect_true(file.exists(tmp))
  parsed <- jsonlite::fromJSON(tmp)
  expect_identical(parsed$module, "m")
  expect_equal(parsed$seed, 7)
})

test_that("arthomix_provenance_download_handler() degrades gracefully instead of erroring when record_fn() itself throws", {
  record_fn <- function() stop("boom - simulated reactive failure")
  handler <- arthomix_provenance_download_handler(record_fn, "test_record")
  content_fn <- provenance_handler_content_fn(handler)
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  expect_no_error(content_fn(tmp))
  parsed <- jsonlite::fromJSON(tmp)
  expect_true(grepl("Could not build", parsed$error))
})
