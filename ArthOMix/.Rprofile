## .Rprofile
## Sourced by R before any app code runs (unlike global.R, which Shiny
## sources only after it has already resolved the port to bind). Pins
options(shiny.host = "127.0.0.1", shiny.port = 7788)

repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = Sys.getenv("RSPM", "https://cloud.r-project.org")))
}

options(renv.install.timeout = 9000L)
