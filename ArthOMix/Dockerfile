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
    libmagick++-dev \
    libglpk-dev \
    libglu1-mesa-dev \
    libgl1-mesa-dev \
    freeglut3-dev \
    libhdf5-dev \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

## MOFA2 runs its factor model through the Python package mofapy2 (the app calls
## MOFA2::run_mofa(use_basilisk = FALSE) and resolves the interpreter from
## RETICULATE_PYTHON). Install it into the image so the MOFA sub-module works
## in the container.
RUN pip3 install --no-cache-dir --break-system-packages "mofapy2==0.7.5" || pip3 install --no-cache-dir "mofapy2==0.7.5"
ENV RETICULATE_PYTHON=/usr/bin/python3

WORKDIR /app

RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
COPY renv.lock renv.lock
RUN R -e "renv::restore(lockfile = 'renv.lock')"

COPY . .

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp(appDir = '/app', host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', 3838)))"]
