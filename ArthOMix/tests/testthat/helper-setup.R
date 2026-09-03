## Sourced once, automatically, before every test file in this directory
## (testthat's helper-*.R convention). testthat::test_dir()/test_file() set
## the working directory to this file's own directory (tests/testthat/)
app_dir <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)

if (!file.exists(file.path(app_dir, "global.R")) || !file.exists(file.path(app_dir, "data_paths.R"))) {
  stop(
    "Could not locate the ArthOMix app root two levels above tests/testthat/. ",
    "Expected global.R and data_paths.R at: ", app_dir
  )
}

source_from_app_root <- function(relative_path, local = FALSE) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(app_dir)
  source(relative_path, local = local)
}

source_from_app_root("data_paths.R")

new_app_driver <- function(..., shiny_args = list()) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tempdir())
  if (is.null(shiny_args$port)) shiny_args$port <- httpuv::randomPort()
  shinytest2::AppDriver$new(app_dir = app_dir, ..., shiny_args = shiny_args)
}

login_test_user <- function(app) {
  app$set_inputs(`auth-login_email` = Sys.getenv("ARTHOMIX_TEST_EMAIL"))
  app$set_inputs(`auth-login_password` = Sys.getenv("ARTHOMIX_TEST_PASSWORD"))
  app$click("auth-login_btn")
  app$wait_for_idle(timeout = 20 * 1000)
}

wait_for_html_containing <- function(app, selector, pattern, timeout = 30, interval = 0.5) {
  deadline <- Sys.time() + timeout
  repeat {
    html <- app$get_html(selector)
    if (grepl(pattern, html, fixed = TRUE)) return(html)
    if (Sys.time() >= deadline) return(html)
    Sys.sleep(interval)
  }
}

wait_for_input_value <- function(app, id, timeout = 30, interval = 0.5) {
  deadline <- Sys.time() + timeout
  repeat {
    val <- app$get_values(input = id)$input[[id]]
    if (!is.null(val) && nzchar(val)) return(val)
    if (Sys.time() >= deadline) return(val)
    Sys.sleep(interval)
  }
}

retry_click <- function(app, id, timeout = 30, interval = 0.5) {
  deadline <- Sys.time() + timeout
  repeat {
    result <- tryCatch({ app$click(id); TRUE }, error = function(e) e)
    if (isTRUE(result)) return(invisible(TRUE))
    if (Sys.time() >= deadline) stop(result)
    Sys.sleep(interval)
  }
}
