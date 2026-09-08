## Regression guard for a Git LFS reproducibility finding (2026-09-07 defense
## audit): if a clone/checkout never ran `git lfs pull` (or git-lfs isn't
## installed), bundled precomputed files under data/preloaded/ are left as
## tiny LFS pointer-stub text files instead of the real data. Before this
## check, the app had no way to detect that up front - it would fail deep
## inside whichever module first tried to read the affected file, with a
## cryptic parse error that gave no hint the real cause was a missing
## `git lfs pull`. arthomix_check_lfs_pointers() (global.R) scans once at
## startup and raises one clear, actionable error instead.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))

test_that("a directory of real (non-pointer) files passes silently", {
  dir <- withr::local_tempdir()
  writeLines("not a pointer file, just some real content", file.path(dir, "real.rds"))
  expect_null(arthomix_check_lfs_pointers(dir))
})

test_that("a Git LFS pointer stub is detected and raises one clear, actionable error naming the file", {
  dir <- withr::local_tempdir()
  sub <- file.path(dir, "methylomics", "matrix")
  dir.create(sub, recursive = TRUE)
  pointer_path <- file.path(sub, "beta_raw.rds")
  writeLines(c(
    "version https://git-lfs.github.com/spec/v1",
    "oid sha256:0000000000000000000000000000000000000000000000000000000000000",
    "size 2141154843"
  ), pointer_path)

  err <- tryCatch(arthomix_check_lfs_pointers(dir), error = function(e) e)
  expect_s3_class(err, "error")
  expect_true(grepl("Git LFS pointer stub", conditionMessage(err)))
  expect_true(grepl("git lfs pull", conditionMessage(err)))
  expect_true(grepl("beta_raw.rds", conditionMessage(err)))
})

test_that("a small real file that happens to be under 200 bytes but isn't an LFS pointer is not flagged", {
  dir <- withr::local_tempdir()
  writeLines("gene,value\nTNF,1.23\n", file.path(dir, "tiny_real.csv"))
  expect_null(arthomix_check_lfs_pointers(dir))
})

test_that("a missing/non-existent root directory is a silent no-op, not an error", {
  expect_null(arthomix_check_lfs_pointers(file.path(withr::local_tempdir(), "does_not_exist")))
})

test_that("multiple pointer stubs are all listed in the single raised error", {
  dir <- withr::local_tempdir()
  make_pointer <- function(name) {
    writeLines(c("version https://git-lfs.github.com/spec/v1", "oid sha256:abc", "size 100"),
               file.path(dir, name))
  }
  make_pointer("a.rds")
  make_pointer("b.csv")
  err <- tryCatch(arthomix_check_lfs_pointers(dir), error = function(e) e)
  expect_true(grepl("2 bundled data file", conditionMessage(err)))
  expect_true(grepl("a.rds", conditionMessage(err)))
  expect_true(grepl("b.csv", conditionMessage(err)))
})
