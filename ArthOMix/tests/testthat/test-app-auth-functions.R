## Authentication module (R/auth/) - pure, non-network functions only.
## R/auth/mod_auth_server.R and R/auth/auth_api.R are mostly thin wrappers
## around live Supabase GoTrue REST calls (supabase_sign_up()/sign_in()/
## verify_otp()/recover()/resend()/update_password()/sign_out(), and the
## internal .supabase_req() they all funnel through) - none of those are
## exercised here, deliberately: they require a live, configured Supabase
## project (SUPABASE_URL/SUPABASE_ANON_KEY) and a real network call, which
## is out of scope for this file. E2E coverage of the actual login/signup
## flow already exists via login_test_user() (helper-setup.R) in the
## shinytest2-based app tests, gated on ARTHOMIX_TEST_EMAIL/PASSWORD being
## set.
##
## What IS covered here is every pure, input -> output validation/formatting
## function in R/auth/ that needs no network and no Supabase configuration:
## is_valid_email(), is_valid_password(), auth_status_ui() (mod_auth_server.R)
## and .supabase_extract_error() (auth_api.R).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "auth", "auth_api.R"))
source_from_app_root(file.path("R", "auth", "mod_auth_server.R"))

## ---- is_valid_email() ------------------------------------------------------

test_that("is_valid_email() accepts ordinary well-formed addresses", {
  expect_true(is_valid_email("user@example.com"))
  expect_true(is_valid_email("first.last+tag@sub.example.co.uk"))
})

test_that("is_valid_email() rejects strings with no '@', no domain dot, or embedded whitespace", {
  expect_false(is_valid_email("not-an-email"))
  expect_false(is_valid_email("user@example"))          ## no '.' after the '@'
  expect_false(is_valid_email("user name@example.com")) ## whitespace before '@'
  expect_false(is_valid_email("user@exam ple.com"))     ## whitespace after '@'
})

test_that("is_valid_email() rejects NULL/NA/empty input rather than erroring", {
  expect_false(is_valid_email(NULL))
  expect_false(is_valid_email(""))
  expect_false(is_valid_email(NA_character_))
})

## ---- is_valid_password() ---------------------------------------------------

test_that("is_valid_password() enforces an 8-character minimum, matching Supabase's own default", {
  expect_false(is_valid_password("short1"))    ## 6 chars
  expect_false(is_valid_password("seven77"))   ## 7 chars
  expect_true(is_valid_password("exactly8"))   ## exactly 8 chars - boundary
  expect_true(is_valid_password("wellabovethelimit123"))
})

test_that("is_valid_password() rejects NULL/empty input rather than erroring", {
  expect_false(is_valid_password(NULL))
  expect_false(is_valid_password(""))
})

test_that("is_valid_password() imposes no arbitrary uppercase/number/symbol requirement", {
  ## Deliberately matches the comment above is_valid_password(): only length
  ## is checked - an all-lowercase, all-letter 8+ char string must pass.
  expect_true(is_valid_password("alllowercase"))
})

## ---- auth_status_ui() -------------------------------------------------------

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

## ---- .supabase_extract_error() ---------------------------------------------

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

## ---- supabase_configured() --------------------------------------------------
## Pure logic (nzchar() on two already-read env vars), but reads module-level
## globals set once at source time from Sys.getenv() - not appropriate to
## flip via Sys.setenv() here (auth_api.R reads them only at source time, not
## per-call, and other test files/E2E tests in this same session rely on
## whatever real SUPABASE_URL/SUPABASE_ANON_KEY this machine has configured
## staying intact). Instead, just check the function is consistent with its
## own documented inputs, whatever they currently are - a real config-vs-
## no-config behavioral check belongs in the E2E suite (login_test_user()).
test_that("supabase_configured() reflects nzchar() of both SUPABASE_URL and SUPABASE_ANON_KEY", {
  expect_identical(supabase_configured(), nzchar(SUPABASE_URL) && nzchar(SUPABASE_ANON_KEY))
})
