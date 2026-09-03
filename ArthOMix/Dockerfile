# ArthOMix Explorer - Dockerfile
FROM rocker/r-ver:4.4.2


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

RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
COPY renv.lock renv.lock
RUN R -e "renv::restore(lockfile = 'renv.lock')"

COPY . .

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp(appDir = '/app', host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', 3838)))"]
