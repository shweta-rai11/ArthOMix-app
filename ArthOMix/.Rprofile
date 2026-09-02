## .Rprofile
## Sourced by R before any app code runs (unlike global.R, which Shiny
## sources only after it has already resolved the port to bind). Pins
## the app to a stable local URL: http://127.0.0.1:7788/
options(shiny.host = "127.0.0.1", shiny.port = 7788)

## Guarantee a CRAN mirror is always set. Without this, CI's bare
## `install.packages("renv")` bootstrap step (run before renv itself is
## available to honor RENV_CONFIG_REPOS_OVERRIDE) fails with "trying to
## use CRAN without setting a mirror", since the RSPM env var from
## setup-r's use-public-rspm doesn't reach that call.
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = Sys.getenv("RSPM", "https://cloud.r-project.org")))
}
