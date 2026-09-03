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

## renv::restore()'s own parallel install graph has an internal deadline
## separate from (and shorter than) the CI job's own timeout-minutes -
## confirmed live: renv/R/graph.R hardcodes
## getOption("renv.install.timeout", default = 3600L), a flat 1-hour cap
## on the ENTIRE restore, not per-package. A from-scratch CI restore of
## this project's ~95 direct (many more transitive) Bioconductor/CRAN
## packages, several minutes each to build from source with no cache yet,
## genuinely needs longer than that once - confirmed live, it was still
## installing cleanly with zero errors when the 1-hour mark killed an
## in-progress package. No env var equivalent exists for this option, so
## it can only be set here, before restore() runs.
options(renv.install.timeout = 9000L)
