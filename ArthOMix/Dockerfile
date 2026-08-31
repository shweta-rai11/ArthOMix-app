# ArthOMix Explorer - Dockerfile
#
# Build context is THIS directory (ArthOMix/ArthOMix/, the app root
# containing global.R/server.R/ui.R/renv.lock/DESCRIPTION) - build from
# here, e.g.:
#   docker build -t arthomix .
#   docker run --rm -p 3838:3838 arthomix
#
# This is a real, working starting point, not an exhaustive solution to
# every possible Bioconductor system-library edge case: this app's
# DESCRIPTION pulls in a large, Bioconductor-heavy dependency set (limma,
# minfi/IlluminaHumanMethylation*, org.Hs.eg.db, mixOmics, WGCNA, SNFtool,
# EpiDISH, ...), and different subsets of those can need additional -dev
# system libraries beyond the common ones installed below. Watch the
# `renv::restore()` step's output on the first real build and add whatever
# it reports missing.

# Pinned to match this project's renv.lock R version (renv.lock -> R ->
# Version: "4.4.2").
FROM rocker/r-ver:4.4.2

# Common system libraries this app's dependency stack needs to compile from
# source (curl/openssl/xml for httr2/xml2-family packages talking to
# Supabase/OpenGWAS/PubMed, pandoc for any rmarkdown-based rendering,
# fontconfig/freetype/harfbuzz for ragg/systemfonts-based plotting, BLAS/
# LAPACK + gfortran for the numerical/statistical packages (WGCNA, mixOmics,
# SNFtool, MRPRESSO, ...)). NOT an exhaustive Bioconductor system-dependency
# list - see the header comment above.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    pandoc \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install renv itself first and restore the exact locked package set before
# copying the rest of the app - so `docker build` only re-runs the (slow)
# package restore step when renv.lock actually changes, not on every source
# edit.
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
COPY renv.lock renv.lock
RUN R -e "renv::restore(lockfile = 'renv.lock')"

# Now bring in the rest of the app - source, data_paths.R, global.R/
# server.R/ui.R, R/, www/, and data/ (the bundled preloaded datasets this
# app is deliberately fully self-contained on - see data_paths.R's own
# header comment; this is several GB, so the resulting image is large by
# design, not by accident).
COPY . .

# Shiny's default port; overridden at runtime via -e PORT=... if needed.
EXPOSE 3838

# Explicit host="0.0.0.0"/port here deliberately override this project's
# own .Rprofile (options(shiny.host = "127.0.0.1", shiny.port = 7788) - a
# convenience pin for a human developer's local browser, see .Rprofile's own
# comment) - inside a container, binding to loopback-only would make the
# app unreachable through Docker's port mapping no matter what -p flag is
# used to run it.
CMD ["R", "-e", "shiny::runApp(appDir = '/app', host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', 3838)))"]
