## Authentication module (R/auth/) - pure, non-network functions only.
## R/auth/mod_auth_server.R and R/auth/auth_api.R are mostly thin wrappers
## around live Supabase GoTrue REST calls (supabase_sign_up()/sign_in()/

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "auth", "auth_api.R"))
source_from_app_root(file.path("R", "auth", "mod_auth_server.R"))

test_that("is_valid_email() accepts ordinary well-formed addresses", {
  expect_true(is_valid_email("user@example.com"))
  expect_true(is_valid_email("first.last+tag@sub.example.co.uk"))
})

test_that("is_valid_email() rejects strings with no '@', no domain dot, or embedded whitespace", {
  expect_false(is_valid_email("not-an-email"))
  expect_false(is_valid_email("user@example"))
  expect_false(is_valid_email("user name@example.com"))
  expect_false(is_valid_email("user@exam ple.com"))
})

test_that("is_valid_email() rejects NULL/NA/empty input rather than erroring", {
  expect_false(is_valid_email(NULL))
  expect_false(is_valid_email(""))
  expect_false(is_valid_email(NA_character_))
})

test_that("is_valid_password() enforces an 8-character minimum, matching Supabase's own default", {
  expect_false(is_valid_password("short1"))
  expect_false(is_valid_password("seven77"))
  expect_true(is_valid_password("exactly8"))
  expect_true(is_valid_password("wellabovethelimit123"))
})

test_that("is_valid_password() rejects NULL/empty input rather than erroring", {
  expect_false(is_valid_password(NULL))
  expect_false(is_valid_password(""))
})

test_that("is_valid_password() imposes no arbitrary uppercase/number/symbol requirement", {
  expect_true(is_valid_password("alllowercase"))
})

test_that("auth_status_ui() renders NULL for a NULL status (no message to show)", {
  expect_null(auth_status_ui(NULL))
})

test_that("auth_status_ui() renders the message text and an error/success CSS class", {
  err_html <- as.character(auth_status_ui(list(msg = "Invalid credentials", type = "error")))
  expect_match(err_html, "Invalid credentials", fixed = TRUE)
  expect_match(err_html, "auth-alert-error", fixed = TRUE)

  ok_html <- as.character(auth_status_ui(list(msg = "Signed in", type = "success")))
  expect_match(ok_html, "Signed in", fixed = TRUE)
  expect_match(ok_html, "auth-alert-success", fixed = TRUE)
})

test_that(".supabase_extract_error() prefers msg, then error_description, then error", {
  expect_equal(.supabase_extract_error(list(msg = "m", error_description = "d", error = "e")), "m")
  expect_equal(.supabase_extract_error(list(error_description = "d", error = "e")), "d")
  expect_equal(.supabase_extract_error(list(error = "e")), "e")
})

test_that(".supabase_extract_error() falls back to a generic message when nothing usable is present", {
  expect_equal(.supabase_extract_error(list()), "Something went wrong. Please try again.")
  expect_equal(.supabase_extract_error(list(msg = "")), "Something went wrong. Please try again.")
  expect_equal(.supabase_extract_error(NULL), "Something went wrong. Please try again.")
})

test_that("supabase_configured() reflects nzchar() of both SUPABASE_URL and SUPABASE_ANON_KEY", {
  expect_identical(supabase_configured(), nzchar(SUPABASE_URL) && nzchar(SUPABASE_ANON_KEY))
})
