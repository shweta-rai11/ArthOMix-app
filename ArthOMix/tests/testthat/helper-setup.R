## Sourced once, automatically, before every test file in this directory
## (testthat's helper-*.R convention). testthat::test_dir()/test_file() set
## the working directory to this file's own directory (tests/testthat/)
## while tests run, regardless of where Rscript/R was originally invoked
## from - so the app root is computed relative to THIS file's location, not
## assumed from the caller's getwd(). data_paths.R/global.R themselves still
## resolve their data root off getwd() exactly as the real app does under
## shiny::runApp() - it's just that "getwd()" during a test run is this
## directory two levels down, not the app root, so callers below always
## explicitly source with app_dir prepended rather than relying on cwd.
app_dir <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)

if (!file.exists(file.path(app_dir, "global.R")) || !file.exists(file.path(app_dir, "data_paths.R"))) {
  stop(
    "Could not locate the ArthOMix app root two levels above tests/testthat/. ",
    "Expected global.R and data_paths.R at: ", app_dir
  )
}

## data_paths.R computes DATA_DIR off getwd() (matching how the real app
## resolves it under shiny::runApp(), which sets cwd to the app dir before
## sourcing global.R) - so it must actually be sourced WITH cwd = app_dir,
## not just referenced by an absolute path, or DATA_DIR would resolve under
## tests/testthat/data instead. Restore the original (test) cwd afterward
## since the rest of testthat's own machinery expects it.
source_from_app_root <- function(relative_path, local = FALSE) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(app_dir)
  source(relative_path, local = local)
}

## data_paths.R itself only uses base R (no package loads) - safe and cheap
## to source once here so every test file can reference its constants
## without individually re-sourcing it or paying for global.R's ~25 package
## loads just to check a path exists.
source_from_app_root("data_paths.R")

## shinytest2::AppDriver$new() calls pkgload::pkg_path(), which walks UP
## from the calling session's cwd looking for any DESCRIPTION file - it
## finds ArthOMix/DESCRIPTION (present for dependency documentation, not an
## installable package) and, since app_dir is that same directory, flips on
## a strict "dev package" mode that sets options(warn = 2) for the entire
## app-launch subprocess. Every package this app loads that happens to be
## built for a slightly newer R patch version than is currently running
## emits an ordinary "was built under R version X" warning - harmless under
## warn=0 (the default, confirmed by directly sourcing global.R without
## shinytest2), but fatal under warn=2. Instantiating AppDriver with cwd
## pointed somewhere with no DESCRIPTION above it (R's own tempdir(), well
## outside this whole project tree) stops pkg_path() from finding one,
## which skips that branch entirely - app_dir itself is still passed as an
## absolute path, so the app launches from the right place regardless.
## .Rprofile pins the app to a fixed host:port (127.0.0.1:7788) so a human
## developer always finds it at the same URL - convenient for manual dev
## work, but that same fixed port is a real collision hazard for E2E tests:
## this project's dev server on 7788 is routinely already running (shared
## across peer Claude sessions - see feedback_arthomix_shared_dev_server
## memory), and a launched-for-test app process still sources that same
## .Rprofile. Confirmed live that without an explicit port override,
## app$get_values() (which needs the app to be running with
## shiny.testmode = TRUE) can 404 against Shiny's test-mode introspection
## endpoint - not because the freshly-launched test app itself lacks test
## mode, but because something in that startup path ends up serving through
## the pre-existing shared instance on 7788 instead, which was never
## started with test.mode = TRUE. Passing an explicit free port here
## isolates every E2E test from that shared instance entirely. Callers may
## still override shiny_args themselves (e.g. to add other runApp args) -
## in that case respect whatever port they asked for instead of inserting
## our own.
new_app_driver <- function(..., shiny_args = list()) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tempdir())
  if (is.null(shiny_args$port)) shiny_args$port <- httpuv::randomPort()
  shinytest2::AppDriver$new(app_dir = app_dir, ..., shiny_args = shiny_args)
}

## The auth gate (R/auth/mod_auth_*.R) means every AppDriver-based E2E test
## now lands on the login screen first, not the app itself - call this
## right after new_app_driver(), before any set_inputs(sidebar_tabs = ...),
## to get past it with one real Supabase account. Requires
## ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD to already be a confirmed
## account in whichever Supabase project SUPABASE_URL/SUPABASE_ANON_KEY (see
## .Renviron) point at - callers should skip_if() when those aren't set
## (see test-app-smoke.R/test-upload-transcriptomics.R) rather than fail.
login_test_user <- function(app) {
  app$set_inputs(`auth-login_email` = Sys.getenv("ARTHOMIX_TEST_EMAIL"))
  app$set_inputs(`auth-login_password` = Sys.getenv("ARTHOMIX_TEST_PASSWORD"))
  app$click("auth-login_btn")
  app$wait_for_idle(timeout = 20 * 1000)
}
