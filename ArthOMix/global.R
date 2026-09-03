## global.R
## ArthOMix Explorer
## Shared setup: packages, data paths, and the module/submodule registry.

source("data_paths.R")

if (file.exists(".Renviron")) readRenviron(".Renviron")

options(shiny.maxRequestSize = 3072 * 1024^2)

ARTHOMIX_TX_ML_SEED <- 1234

safe_read_rds <- function(path, max_size_mb = 1024,
                           allowed_classes = c("matrix", "array", "data.frame", "numeric",
                                                "integer", "character", "list", "factor")) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(list(ok = FALSE, value = NULL, error = "File not found."))
  }
  sz_mb <- as.numeric(file.info(path)$size) / (1024^2)
  if (!is.finite(sz_mb) || sz_mb > max_size_mb) {
    return(list(ok = FALSE, value = NULL,
                error = sprintf("This .rds file is larger than the %sMB limit for uploads.", max_size_mb)))
  }
  if (!requireNamespace("callr", quietly = TRUE)) {
    return(list(ok = FALSE, value = NULL, error = "The .rds upload path requires the 'callr' package, which is not installed."))
  }
  obj <- tryCatch(
    callr::r(function(p) readRDS(p), args = list(path), timeout = 60),
    error = function(e) e
  )
  if (inherits(obj, "condition")) {
    return(list(ok = FALSE, value = NULL, error = "Could not read this .rds file (it may be corrupt, oversized, or in an unsupported format)."))
  }
  cls <- class(obj)
  if (!any(cls %in% allowed_classes)) {
    return(list(ok = FALSE, value = NULL, error = sprintf(
      "This .rds file contains a %s object; only plain matrices, data frames, or vectors are accepted here.",
      paste(cls, collapse = "/"))))
  }
  list(ok = TRUE, value = obj, error = NULL)
}

ARTHOMIX_ASYNC_AVAILABLE <- requireNamespace("future", quietly = TRUE) && requireNamespace("promises", quietly = TRUE)
if (ARTHOMIX_ASYNC_AVAILABLE) {
  future::plan(future::multisession, workers = 2)
  invisible(future::value(lapply(seq_len(2), function(i) {
    future::future({
      requireNamespace("mixOmics", quietly = TRUE)
      requireNamespace("SNFtool", quietly = TRUE)
      requireNamespace("MOFA2", quietly = TRUE)
      TRUE
    }, seed = FALSE)
  })))
}

library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(ggrepel)
library(plotly)
library(dplyr)
library(tidyr)
library(data.table)

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(WGCNA)
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(caret)
  library(pROC)
  library(sva)
  library(VennDiagram)
  library(ggVennDiagram)
  library(MCPcounter)
  library(IOBR)
  library(rms)
  library(MendelianRandomization)
  library(coloc)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
})

## parallel::detectCores() can return NA in restricted/containerized
## environments (confirmed on GitHub Actions' runner), which propagates
## through max() and fails WGCNA's own "must be numeric and at least 2"
## check - so NA needs an explicit fallback. Separately, WGCNA requires
## nThreads >= 2, so the floor here must be max(2, ...), not max(1, ...):
## on any machine with 3 or fewer cores, detectCores() - 2 is 0 or 1,
## which would fail that check even with a valid, non-NA core count.
n_wgcna_cores <- parallel::detectCores()
if (is.na(n_wgcna_cores)) n_wgcna_cores <- 2L
suppressMessages(
  WGCNA::enableWGCNAThreads(nThreads = max(2, n_wgcna_cores - 2))
)

library(ellmer)
library(shinychat)

validate <- shiny::validate
cor <- WGCNA::cor

list_gene_panels <- function() {
  files <- list.files(GENE_PANELS_DIR, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(character(0))
  labels <- tools::toTitleCase(gsub("_", " ", tools::file_path_sans_ext(basename(files))))
  stats::setNames(files, labels)
}

load_gene_panel <- function(path) {
  genes <- trimws(readLines(path, warn = FALSE))
  unique(genes[nzchar(genes)])
}

addResourcePath("figures", FIGURES_DIR)

load_default_dataset <- function() {
  expr <- readRDS(DEFAULT_EXPR_RDS)
  meta <- as.data.frame(data.table::fread(DEFAULT_META_CSV, showProgress = FALSE))
  common <- intersect(colnames(expr), meta$sample)
  expr <- expr[, common, drop = FALSE]
  meta <- meta[match(common, meta$sample), , drop = FALSE]
  list(expr = expr, meta = meta, source = "Example dataset: sex-stratified RA blood cohort (GSE93272 + GSE110169)")
}

.arthomix_cache <- new.env(parent = emptyenv())

get_raw_eset <- function(gse_id) {
  key <- paste0("eset_", gse_id)
  if (is.null(.arthomix_cache[[key]])) {
    path <- file.path(RAW_DIR, paste0(gse_id, "_raw.rds"))
    if (!file.exists(path)) return(NULL)
    e <- readRDS(path)
    if (is(e, "list")) e <- e[[1]]
    .arthomix_cache[[key]] <- e
  }
  .arthomix_cache[[key]]
}

get_or_compute_wgcna_blocks <- function(key_parts, compute_fn) {
  cache_key <- paste0("wgcna_", digest::digest(key_parts, algo = "xxhash64"))
  if (!is.null(.arthomix_cache[[cache_key]])) return(.arthomix_cache[[cache_key]])
  disk_path <- file.path(WGCNA_CACHE_DIR, paste0(cache_key, ".rds"))
  result <- if (file.exists(disk_path)) {
    readRDS(disk_path)
  } else {
    result <- compute_fn()
    tryCatch(saveRDS(result, disk_path), error = function(e) NULL)
    result
  }
  .arthomix_cache[[cache_key]] <- result
  result
}

get_collapsed_genes <- function(gse_id) {
  key <- paste0("collapsed_", gse_id)
  if (is.null(.arthomix_cache[[key]])) {
    disk_path <- file.path(COLLAPSED_CACHE_DIR, paste0(gse_id, "_collapsed.rds"))
    if (file.exists(disk_path)) {
      .arthomix_cache[[key]] <- readRDS(disk_path)
    } else {
      eset <- get_raw_eset(gse_id)
      if (is.null(eset)) return(NULL)
      collapsed <- collapse_probes_to_genes(eset)
      .arthomix_cache[[key]] <- collapsed
      tryCatch(saveRDS(collapsed, disk_path), error = function(e) NULL)
    }
  }
  .arthomix_cache[[key]]
}

METH_QC_PROBE_CASCADE <- data.frame(
  step = c("Raw (post-parse)", "cg-prefix only", "MASK_general removal", "Multi-hit removal", "Sex-chromosome removal", "Missingness > 5% removal"),
  retained = c(485577L, 482421L, 422520L, 422520L, 412492L, 412492L),
  removed  = c(NA_integer_, 3156L, 59901L, 0L, 10028L, 0L)
)

load_default_dmp <- function(stage = c("plain", "sva"), sex = c("female", "male", "all")) {
  stage <- match.arg(stage); sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  dir <- if (stage == "plain") METH_DMP_PLAIN_DIR else METH_DMP_SVA_DIR
  path <- file.path(dir, sprintf("dmp_%s_full.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_dmr <- function(sex = c("female", "male", "all")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DMR_DIR, sprintf("dmr_%s_full.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_meth_pheno <- function() {
  if (!METH_DATA_AVAILABLE || !file.exists(METH_QC_PHENO_CSV)) return(NULL)
  as.data.frame(data.table::fread(METH_QC_PHENO_CSV, showProgress = FALSE))
}

load_default_meth_qc_sexcheck <- function() {
  if (!METH_DATA_AVAILABLE || !file.exists(METH_QC_PCA_SEXCHECK_CSV)) return(NULL)
  as.data.frame(data.table::fread(METH_QC_PCA_SEXCHECK_CSV, showProgress = FALSE))
}

load_default_meth_matrix <- function() {
  key <- "meth_default_matrix"
  if (!is.null(.arthomix_cache[[key]])) return(.arthomix_cache[[key]])
  if (!METH_RAW_DATA_AVAILABLE) return(NULL)

  beta  <- readRDS(METH_BETA_RAW_RDS)
  pheno <- readRDS(METH_PHENO_RDS)
  pheno$sample <- as.character(pheno$gsm)

  common <- intersect(colnames(beta), pheno$sample)
  beta  <- beta[, common, drop = FALSE]
  pheno <- pheno[match(common, pheno$sample), , drop = FALSE]

  result <- list(beta = beta, pheno = pheno)
  .arthomix_cache[[key]] <- result
  result
}

load_default_cx_table <- function(label) {
  path <- CX_TABLE_REGISTRY[[label]]
  if (is.null(path) || !CX_DATA_AVAILABLE || !file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

MULTI_MOFA_AVAILABLE <- requireNamespace("MOFA2", quietly = TRUE)

MULTI_DIABLO_LIVE_AVAILABLE <- requireNamespace("mixOmics", quietly = TRUE)
MULTI_SNF_LIVE_AVAILABLE <- requireNamespace("SNFtool", quietly = TRUE)

get_or_compute_meth_wgcna_blocks <- function(key_parts, compute_fn) {
  cache_key <- paste0("methwgcna_", digest::digest(key_parts, algo = "xxhash64"))
  if (!is.null(.arthomix_cache[[cache_key]])) return(.arthomix_cache[[cache_key]])
  disk_path <- file.path(METH_WGCNA_CACHE_DIR, paste0(cache_key, ".rds"))
  result <- if (file.exists(disk_path)) {
    readRDS(disk_path)
  } else {
    result <- compute_fn()
    tryCatch(saveRDS(result, disk_path), error = function(e) NULL)
    result
  }
  .arthomix_cache[[cache_key]] <- result
  result
}

load_default_wgcna_module_trait <- function(sex = c("female", "male"), merged = FALSE) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_WGCNA_DIR, sprintf("module_trait_%s%s.csv", sex, if (merged) "_merged10" else ""))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_wgcna_module_assignment <- function(sex = c("female", "male"), merged = FALSE) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_WGCNA_DIR, sprintf("module_assignment_%s%s.csv", sex, if (merged) "_merged10" else ""))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_dmr_biomarker_panel <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DMR_DIR, sprintf("biomarker_panel_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_mr_estimates <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, sprintf("mr_estimates_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_mr_harmonised <- function() {
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, "mr_harmonised_all_cpgs.csv")
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_mr_steiger <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, sprintf("mr_steiger_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_mr_instrument_counts <- function() {
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, "instrument_counts.csv")
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_meth_coloc_results <- function() {
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, "coloc_results.csv")
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_diagnostic_ensemble_votes <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DIAGNOSTIC_VOTES_DIR, sprintf("ensemble_votes_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_diagnostic_panel_auc <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DIAGNOSTIC_DIR, sprintf("diagnostic_panel_auc_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_diagnostic_perprobe_auc <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DIAGNOSTIC_DIR, sprintf("diagnostic_perprobe_auc_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_diagnostic_train_test <- function() {
  key <- "diag_train_test"
  if (!is.null(.arthomix_cache[[key]])) return(.arthomix_cache[[key]])
  if (!METH_DIAG_DATA_AVAILABLE) return(NULL)

  internal <- readRDS(METH_DIAG_INTERNAL_RDS)
  external <- readRDS(METH_DIAG_EXTERNAL_RDS)
  internal$pheno <- as.data.frame(internal$pheno)
  external$pheno <- as.data.frame(external$pheno)

  result <- list(internal = internal, external = external)
  .arthomix_cache[[key]] <- result
  result
}

GEO_SOURCES <- list(
  list(gse = "GSE93272",  role = "Training (whole blood)",   used_in = "Merged into the example cohort"),
  list(gse = "GSE110169", role = "Training (whole blood)",   used_in = "Merged into the example cohort"),
  list(gse = "GSE15573",  role = "Validation (blood, PBMC)", used_in = "Used later, for cross-ancestry validation"),
  list(gse = "GSE89408",  role = "Validation (synovium)",    used_in = "Used later, for cross-tissue validation")
)
geo_link <- function(gse) paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse)

collapse_probes_to_genes <- function(eset) {
  fd <- Biobase::fData(eset)
  col <- grep("^gene[ ._]?symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  ex <- Biobase::exprs(eset)
  if (is.na(col)) return(structure(ex, collapsed = FALSE))
  sym <- as.character(fd[[col]])
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
  ex <- ex[keep, , drop = FALSE]; sym <- sym[keep]
  ok <- rowSums(is.na(ex)) < ncol(ex)
  ex <- ex[ok, , drop = FALSE]; sym <- sym[ok]
  if (nrow(ex) == 0) return(structure(ex, collapsed = FALSE))
  out <- WGCNA::collapseRows(
    ex, rowGroup = sym, rowID = rownames(ex), method = "MaxMean",
    connectivityBasedCollapsing = FALSE, connectivityPower = 1,
    selectFewestMissing = TRUE, thresholdCombine = NA
  )$datETcollapsed
  structure(out, collapsed = TRUE)
}

eset_harmonize_meta <- function(eset, name) {
  pd <- Biobase::pData(eset)
  clin <- grep(":ch1$", colnames(pd), value = TRUE)
  discol <- grep("disease|status", clin, value = TRUE, ignore.case = TRUE)[1]
  sexcol <- grep("gender|sex", clin, value = TRUE, ignore.case = TRUE)[1]
  dis <- if (!is.na(discol)) as.character(pd[[discol]]) else rep(NA_character_, nrow(pd))
  grp <- ifelse(grepl("sle|lupus", dis, ignore.case = TRUE), "SLE",
          ifelse(grepl("normal|healthy|control", dis, ignore.case = TRUE), "HC",
          ifelse(grepl("rheumatoid|(^|[^a-z])ra([^a-z]|$)", dis, ignore.case = TRUE), "RA", "other")))
  sx <- if (!is.na(sexcol)) toupper(substr(gsub("[^A-Za-z]", "", as.character(pd[[sexcol]])), 1, 1)) else rep(NA_character_, nrow(pd))
  sx[!sx %in% c("F", "M")] <- NA
  data.frame(sample = rownames(pd), dataset = name, group = grp, sex = sx, stringsAsFactors = FALSE)
}

merged_training_subset <- function(gse_id) {
  d <- load_default_dataset()
  keep <- d$meta$dataset == gse_id
  if (!any(keep)) return(NULL)
  samples <- d$meta$sample[keep]
  list(
    expr = d$expr[, samples, drop = FALSE],
    meta = d$meta[keep, , drop = FALSE],
    label = paste0(gse_id, " (from the merged, batch-corrected training cohort - source-level raw probe data is not available in this deployment)")
  )
}

load_individual_dataset <- function(gse_id) {
  key <- paste0("indiv_", gse_id)
  if (!is.null(.arthomix_cache[[key]])) return(.arthomix_cache[[key]])

  result <- if (identical(gse_id, "GSE89408")) {
    path <- file.path(RAW_DIR, "GSE89408_counts.txt.gz")
    eset <- get_raw_eset("GSE89408")
    if (!file.exists(path) || is.null(eset)) {
      NULL
    } else {
      cnt <- as.matrix(utils::read.delim(gzfile(path), row.names = 1, check.names = FALSE))
      pref <- gsub("_[0-9]+$", "", colnames(cnt))
      dis_map <- c(normal_tissue = "HC", RA_tissue = "RA", OA_tissue = "other",
                    AG_tissue = "other", undiff_tissue = "other")
      grp <- unname(dis_map[pref])
      sx <- Biobase::pData(eset)[["Sex:ch1"]]
      meta <- data.frame(sample = colnames(cnt), dataset = "GSE89408",
                          group = grp, sex = sx, stringsAsFactors = FALSE)
      keep <- !is.na(meta$group)
      list(expr = cnt[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE],
           label = "GSE89408 (synovium, RNA-seq raw counts)")
    }
  } else {
    eset <- get_raw_eset(gse_id)
    if (is.null(eset) || nrow(Biobase::exprs(eset)) == 0) {
      merged_training_subset(gse_id)
    } else {
      expr <- Biobase::exprs(eset)
      meta <- eset_harmonize_meta(eset, gse_id)
      keep <- !is.na(meta$group)
      list(expr = expr[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE],
           label = paste0(gse_id, " (", Biobase::annotation(eset), ", microarray, raw probes)"))
    }
  }
  .arthomix_cache[[key]] <- result
  result
}

`%||%` <- function(a, b) if (is.null(a)) b else a

guess_column_by_name <- function(cols, exact, contains = exact, fallback = cols[1]) {
  hit <- cols[tolower(cols) %in% tolower(exact)]
  if (length(hit) > 0) return(hit[1])
  hit <- cols[grepl(paste(contains, collapse = "|"), cols, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[1])
  fallback
}

ARTHOMIX_OLLAMA_MODEL <- "qwen3:8b"
ollama_base_url <- function() Sys.getenv("OLLAMA_BASE_URL", "http://localhost:11434")

ollama_available <- function() {
  tryCatch({
    con <- url(paste0(ollama_base_url(), "/api/tags"), open = "rb")
    on.exit(close(con), add = TRUE)
    txt <- paste(readLines(con, warn = FALSE), collapse = "")
    grepl(ARTHOMIX_OLLAMA_MODEL, txt, fixed = TRUE)
  }, error = function(e) FALSE)
}

ARTHOCHAT_ANTHROPIC_MODEL <- "claude-sonnet-5"
anthropic_available <- function() nzchar(Sys.getenv("ANTHROPIC_API_KEY", ""))

## ArthOChat prefers a hosted Anthropic model when ANTHROPIC_API_KEY is set -
## this is what a deployed app (e.g. shinyapps.io) uses, since it has no
## local Ollama server reachable and needs to work for any visitor with
## nothing left running on a developer's own machine. It falls back to local
## Ollama for offline development where no API key is configured.
arthochat_backend <- function() {
  if (anthropic_available()) "anthropic"
  else if (ollama_available()) "ollama"
  else "none"
}

pubmed_search <- function(query, max_results = 5) {
  max_results <- min(max(as.integer(max_results %||% 5), 1L), 10L)
  esearch <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") %>%
    httr2::req_url_query(db = "pubmed", term = query, retmode = "json", retmax = max_results, sort = "relevance") %>%
    httr2::req_url_query(tool = "arthomix-explorer") %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()
  ids <- unlist(esearch$esearchresult$idlist)
  if (length(ids) == 0) return("No PubMed results found for this query.")

  esummary <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
    httr2::req_url_query(db = "pubmed", id = paste(ids, collapse = ","), retmode = "json") %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()
  recs <- esummary$result

  refs <- lapply(ids, function(pmid) {
    r <- recs[[pmid]]
    if (is.null(r)) return(NULL)
    authors <- vapply(r$authors, function(a) a$name %||% "?", character(1))
    author_str <- if (length(authors) > 1) paste0(authors[1], " et al.") else if (length(authors) == 1) authors[1] else "Unknown author"
    year <- substr(r$pubdate %||% "", 1, 4)
    sprintf("%s (%s). %s. %s. PMID: %s (https://pubmed.ncbi.nlm.nih.gov/%s)",
            author_str, year, r$title, r$fulljournalname %||% r$source %||% "", pmid, pmid)
  })
  refs <- Filter(Negate(is.null), refs)
  if (length(refs) == 0) return("No PubMed results found for this query.")
  paste(refs, collapse = "\n")
}

opengwas_token_configured <- function() nzchar(Sys.getenv("OPENGWAS_JWT", ""))

gwas_catalog_search <- function(query, max_results = 10) {
  if (!opengwas_token_configured()) {
    return(paste(
      "OpenGWAS dataset search isn't configured in this deployment - it needs a free",
      "personal access token. Sign in at https://api.opengwas.io/profile/, copy the JWT",
      "shown there, and set it as the OPENGWAS_JWT environment variable before launching",
      "this app (Sys.setenv(OPENGWAS_JWT = \"...\") or in .Renviron), then restart the app.",
      "In the meantime, the Mendelian Randomization tab's \"Upload your own GWAS summary",
      "statistics\" option lets you run MR directly on your own files instead of an",
      "OpenGWAS ID."
    ))
  }
  max_results <- min(max(as.integer(max_results %||% 10), 1L), 25L)

  cat_df <- .arthomix_cache[["opengwas_catalog"]]
  if (is.null(cat_df)) {
    cat_df <- tryCatch(as.data.frame(ieugwasr::gwasinfo()), error = function(e) NULL)
    if (is.null(cat_df) || !nrow(cat_df)) {
      return("Could not reach the OpenGWAS catalogue (network issue, or an invalid/expired token). Try again shortly, or use the Mendelian Randomization tab's upload option instead.")
    }
    .arthomix_cache[["opengwas_catalog"]] <- cat_df
  }

  hit <- rep(FALSE, nrow(cat_df))
  if ("trait" %in% names(cat_df)) hit <- hit | grepl(query, cat_df$trait, ignore.case = TRUE)
  if ("consortium" %in% names(cat_df)) hit <- hit | grepl(query, cat_df$consortium, ignore.case = TRUE)
  if ("author" %in% names(cat_df)) hit <- hit | grepl(query, cat_df$author, ignore.case = TRUE)
  hits <- cat_df[hit, , drop = FALSE]
  if (!nrow(hits)) return(sprintf("No OpenGWAS datasets found matching \"%s\". Try a broader term (e.g. a disease name without qualifiers).", query))
  if ("sample_size" %in% names(hits)) hits <- hits[order(-hits$sample_size), , drop = FALSE]
  hits <- utils::head(hits, max_results)

  col_or_na <- function(df, col) if (col %in% names(df)) as.character(df[[col]]) else rep(NA_character_, nrow(df))
  lines <- sprintf(
    "%s | %s | population: %s | N: %s | year: %s",
    col_or_na(hits, "id"), col_or_na(hits, "trait"), col_or_na(hits, "population"),
    col_or_na(hits, "sample_size"), col_or_na(hits, "year")
  )
  paste(c(sprintf("Top %d OpenGWAS datasets matching \"%s\" (id | trait | population | N | year):", length(lines), query), lines), collapse = "\n")
}

ARTHOMIX_METHODS_TOPICS <- list(
  list(id = "overview",         title = "Overview and Datasets",       section = "2.1",  satellite = NULL,
       aliases = c("dataset", "datasets", "geo", "cohort")),
  list(id = "preprocessing",    title = "Preprocessing and Batch Correction", section = "2.2", satellite = NULL,
       aliases = c("preprocessing", "pre-processing", "normalisation", "normalization", "batch correction", "batch effect", "combat", "quantile normalisation")),
  list(id = "dge",              title = "Differential Expression",     section = "2.3",  satellite = NULL,
       aliases = c("dge", "differential gene expression", "limma", "deg")),
  list(id = "wgcna",            title = "WGCNA Co-expression Network", section = "2.4",  satellite = "METHODS_2.4_WGCNA_expanded.md",
       aliases = c("wgcna", "co-expression", "coexpression", "gene network", "module detection")),
  list(id = "candidates",       title = "Candidate Gene Identification", section = "2.5", satellite = NULL,
       aliases = c("candidate gene", "candidate genes")),
  list(id = "mr",               title = "Mendelian Randomization",     section = "2.6",  satellite = "METHODS_2.6_mendelian_randomisation.md",
       aliases = c("mr", "mendelian randomisation", "mendelian randomization", "instrumental variable", "causal inference", "eqtl")),
  list(id = "coloc",            title = "Colocalization",              section = "2.7",  satellite = "METHODS_2.7_colocalisation.md",
       aliases = c("coloc", "colocalisation", "colocalization", "coloc.abf", "coloc.susie")),
  list(id = "featureselection", title = "Feature Selection",           section = "2.8",  satellite = "METHODS_2.8_feature_selection.md",
       aliases = c("feature selection", "lasso", "random forest", "svm-rfe", "svm rfe", "biomarker panel")),
  list(id = "diagnostic",       title = "Diagnostic Model",            section = "2.9",  satellite = "METHODS_2.9_diagnostic_model.md",
       aliases = c("diagnostic model", "classifier", "auc", "roc", "cross-validation")),
  list(id = "interaction",      title = "Sex Interaction Analysis",    section = "2.10", satellite = NULL,
       aliases = c("interaction", "sex interaction", "sex-differential", "sex differential")),
  list(id = "crosstissue",      title = "Cross-Tissue Validation",     section = "2.11", satellite = "METHODS_2.11_crosstissue.md",
       aliases = c("cross-tissue", "cross tissue", "synovium", "synovial")),
  list(id = "crossancestry",    title = "Cross-Ancestry Validation",   section = "2.12", satellite = "METHODS_2.12_crossancestry.md",
       aliases = c("cross-ancestry", "cross ancestry", "ancestry")),
  list(id = "enrichment",       title = "Functional Enrichment",       section = "2.13", satellite = "METHODS_2.13_functional_enrichment.md",
       aliases = c("enrichment", "go term", "kegg", "pathway analysis", "gene ontology")),
  list(id = "deconvolution",    title = "Immune Deconvolution",        section = "2.14", satellite = "METHODS_2.14_deconvolution.md",
       aliases = c("deconvolution", "immune cell", "cell composition", "mcp-counter", "mcpcounter")),
  list(id = "nomogram",         title = "Clinical Utility Nomogram",   section = "2.15", satellite = "METHODS_2.15_clinical_utility_nomogram.md",
       aliases = c("nomogram", "clinical utility", "decision curve"))
)

extract_chapter_section <- function(section) {
  if (!file.exists(PROJECT_CHAPTER_MD)) return(NULL)
  lines <- readLines(PROJECT_CHAPTER_MD, warn = FALSE)
  marker <- paste0("^\\*\\*", gsub("\\.", "\\\\.", section), " ")
  start <- grep(marker, lines)
  if (length(start) == 0) return(NULL)
  start <- start[1]
  all_markers <- grep("^\\*\\*2\\.", lines)
  later <- all_markers[all_markers > start]
  end <- if (length(later) > 0) min(later) - 1 else length(lines)
  paste(lines[start:end], collapse = "\n")
}

extract_references_section <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  start <- grep("^## References cited", lines, ignore.case = TRUE)
  if (length(start) == 0) return(NULL)
  paste(lines[start[1]:length(lines)], collapse = "\n")
}

.truncate_with_note <- function(txt, max_chars) {
  if (nchar(txt) <= max_chars) return(txt)
  paste0(substr(txt, 1, max_chars), "\n[...truncated - ask a more specific follow-up for more detail...]")
}

project_methods <- function(module) {
  q <- tolower(trimws(module %||% ""))
  hit <- Find(function(t) {
    q == tolower(t$id) || q == t$section ||
      grepl(q, tolower(t$title), fixed = TRUE) || grepl(tolower(t$title), q, fixed = TRUE) ||
      any(vapply(t$aliases, function(a) grepl(a, q, fixed = TRUE) || grepl(q, a, fixed = TRUE), logical(1)))
  }, ARTHOMIX_METHODS_TOPICS)

  if (is.null(hit)) {
    return(paste0(
      "No module matched \"", module, "\". Available topics: ",
      paste(vapply(ARTHOMIX_METHODS_TOPICS, function(t) sprintf("%s (Section %s)", t$title, t$section), character(1)), collapse = "; ")
    ))
  }

  narrative <- extract_chapter_section(hit$section)
  refs <- if (!is.null(hit$satellite)) extract_references_section(file.path(DATA_ROOT, "results", hit$satellite)) else NULL

  parts <- c(
    sprintf("## %s (Section %s) - this project's own methodology", hit$title, hit$section),
    if (!is.null(narrative)) .truncate_with_note(narrative, 3000) else "(No chapter write-up found for this section.)",
    "",
    if (!is.null(refs)) .truncate_with_note(refs, 3000) else "(No dedicated reference list for this section - use pubmed_search for supporting literature.)"
  )
  paste(parts, collapse = "\n")
}

METH_METHODS_TOPICS <- list(
  list(id = "qc",             title = "Quality Control",
       files = "script01_dataload_QC/METHODS_load_qc.md",
       aliases = c("qc", "quality control", "detection p-value", "bead count", "sex check", "probe filtering")),
  list(id = "normalization",  title = "Normalization",
       files = "script01_dataload_QC/METHODS_load_qc.md",
       aliases = c("normalization", "normalisation", "dasen", "bmiq", "noob", "funnorm", "swan", "pbc", "quantile normalisation")),
  list(id = "celltype",       title = "Cell-Type Deconvolution",
       files = "script02_celltype/METHODS_celltype.md",
       aliases = c("celltype", "cell type", "cell composition", "houseman", "deconvolution")),
  list(id = "dmp",            title = "Differentially Methylated Positions (DMPs)",
       files = c("script03_dmp_sexstratified/METHODS_dmp_sexstratified.md", "script03_dmp_sva_sexstratified/METHODS_dmp_sva_sexstratified.md"),
       aliases = c("dmp", "differentially methylated position", "differential methylation", "sva", "delta beta")),
  list(id = "dmr",            title = "Differentially Methylated Regions (DMRs)",
       files = "script04_dmr_sexstratified/METHODS_dmr_sexstratified.md",
       aliases = c("dmr", "differentially methylated region", "dmrcate", "comb-p", "bumphunter")),
  list(id = "wgcna",          title = "WGCNA (Co-Methylation Network)",
       files = "script05_wgcna_sexstratified/METHODS_wgcna_sexstratified.md",
       aliases = c("wgcna", "co-methylation", "comethylation", "module eigengene", "soft threshold", "network")),
  list(id = "candidates",     title = "Candidate CpGs (Module-DMR Overlap)",
       files = "script06_module_dmr_venn/METHODS_module_dmr_venn.md",
       aliases = c("candidate cpg", "candidate cpgs", "module-dmr", "venn", "overlap")),
  list(id = "featureselection", title = "ML Feature Selection",
       files = "script07_ml_feature_selection/METHODS_ml_feature_selection.md",
       aliases = c("feature selection", "lasso", "random forest", "svm-rfe", "svm rfe", "elastic net", "boruta")),
  list(id = "mr",             title = "Mendelian Randomization",
       files = "script08_mendelian_randomization/METHODS_mendelian_randomization.md",
       aliases = c("mr", "mendelian randomisation", "mendelian randomization", "mqtl", "instrumental variable", "causal inference")),
  list(id = "diagnostic",     title = "Diagnostic Classifier",
       files = "script09_diagnostic_classifier/METHODS_diagnostic_classifier.md",
       aliases = c("diagnostic model", "diagnostic classifier", "classifier", "auc", "roc", "panel"))
)

project_methods_methylomics <- function(module) {
  q <- tolower(trimws(module %||% ""))
  hit <- Find(function(t) {
    q == tolower(t$id) ||
      grepl(q, tolower(t$title), fixed = TRUE) || grepl(tolower(t$title), q, fixed = TRUE) ||
      any(vapply(t$aliases, function(a) grepl(a, q, fixed = TRUE) || grepl(q, a, fixed = TRUE), logical(1)))
  }, METH_METHODS_TOPICS)

  if (is.null(hit)) {
    return(paste0(
      "No Methylomics module matched \"", module, "\". Available topics: ",
      paste(vapply(METH_METHODS_TOPICS, function(t) t$title, character(1)), collapse = "; ")
    ))
  }

  paths <- file.path(METH_DATA_ROOT, hit$files)
  bodies <- Filter(Negate(is.null), lapply(paths, function(p) {
    if (!file.exists(p)) return(NULL)
    paste(readLines(p, warn = FALSE), collapse = "\n")
  }))

  header <- sprintf("## %s - this project's own methodology", hit$title)
  if (length(bodies) == 0) {
    return(paste(header, "(No write-up available for this stage in this deployment.)", sep = "\n"))
  }
  paste(header, "", .truncate_with_note(paste(bodies, collapse = "\n\n---\n\n"), 4000), sep = "\n")
}

read_table_safe <- function(filename, dir = TABLES_DIR) {
  path <- file.path(dir, filename)
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

figure_exists <- function(filename) {
  file.exists(file.path(FIGURES_DIR, filename))
}

GWAS_COL_PATTERNS <- list(
  snp  = c("^snp$", "^rsid$", "^rs_?id$", "variant"),
  ea   = c("^effect_allele$", "^ea$", "^a1$", "^alt$"),
  oa   = c("^other_allele$", "^oa$", "^a2$", "^ref$", "^nea$"),
  beta = c("^beta$", "^b$", "^effect$"),
  se   = c("^se$", "^standard_error$", "^stderr$"),
  pval = c("^p$", "^pval$", "^p_value$", "^pvalue$"),
  eaf  = c("^eaf$", "^freq$", "^maf$", "^effect_allele_freq"),
  n    = c("^n$", "^samplesize$", "^sample_size$", "^total_n$")
)

guess_gwas_col <- function(cols, patterns) {
  for (p in patterns) { hit <- grep(p, cols, ignore.case = TRUE, value = TRUE); if (length(hit)) return(hit[1]) }
  NA_character_
}

read_uploaded_table <- function(path) tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) NULL)

tx_parse_expr_matrix_rds <- function(datapath) {
  loaded <- safe_read_rds(datapath)
  if (!isTRUE(loaded$ok)) {
    return(list(ok = FALSE, mat = NULL, error = loaded$error))
  }
  x <- loaded$value
  if (is.matrix(x) && is.numeric(x)) return(list(ok = TRUE, mat = x, error = NULL))
  if (is.data.frame(x) && ncol(x) >= 2) {
    ids <- as.character(x[[1]])
    if (any(duplicated(ids))) {
      return(list(ok = FALSE, mat = NULL, error = sprintf("%d duplicated feature ID(s) in the first column.", sum(duplicated(ids)))))
    }
    rest <- x[, -1, drop = FALSE]
    all_numeric <- all(vapply(rest, is.numeric, logical(1)))
    if (!all_numeric) return(list(ok = FALSE, mat = NULL, error = "Every column after the first must be numeric."))
    m <- tryCatch({
      mm <- as.matrix(rest); storage.mode(mm) <- "double"; rownames(mm) <- ids; mm
    }, error = function(e) NULL)
    if (is.null(m)) return(list(ok = FALSE, mat = NULL, error = "Every column after the first must be numeric."))
    return(list(ok = TRUE, mat = m, error = NULL))
  }
  list(ok = FALSE, mat = NULL, error = "The .rds file must contain either a numeric feature x sample matrix, or a data.frame with feature IDs in the first column.")
}

gwas_col_map_ui <- function(ns, file_input, df_reactive, prefix, label, extra_fields = character(0)) {
  extra_labels <- c(n = "Sample size (N)")
  renderUI({
    fi <- file_input()
    req(fi)
    df <- df_reactive()
    validate(need(!is.null(df) && ncol(df) >= 4, sprintf("Could not read the %s file as a delimited table with at least 4 columns.", label)))
    cols <- colnames(df)
    guess <- function(key) { g <- guess_gwas_col(cols, GWAS_COL_PATTERNS[[key]]); if (is.na(g)) cols[1] else g }
    tagList(
      div(class = "empty-note", style = "font-size:12.5px; margin: 4px 0; padding: 8px 12px;",
        icon("circle-check"), strong(sprintf(" Uploaded: %s", fi$name)),
        sprintf(" - %s rows, %d columns. Columns auto-matched below - check them.", format(nrow(df), big.mark = ","), ncol(df))),
      fluidRow(
        column(6, selectInput(ns(paste0(prefix, "_snp")), "SNP ID", cols, selected = guess("snp"))),
        column(6, selectInput(ns(paste0(prefix, "_beta")), "Beta (effect size)", cols, selected = guess("beta")))
      ),
      fluidRow(
        column(6, selectInput(ns(paste0(prefix, "_se")), "Standard error", cols, selected = guess("se"))),
        column(6, selectInput(ns(paste0(prefix, "_pval")), "p-value", cols, selected = guess("pval")))
      ),
      fluidRow(
        column(6, selectInput(ns(paste0(prefix, "_ea")), "Effect allele", cols, selected = guess("ea"))),
        column(6, selectInput(ns(paste0(prefix, "_oa")), "Other allele", cols, selected = guess("oa")))
      ),
      selectInput(ns(paste0(prefix, "_eaf")), "Effect allele frequency (optional - improves palindromic SNP resolution)",
                  c("(none)" = "", cols), selected = { g <- guess_gwas_col(cols, GWAS_COL_PATTERNS$eaf); if (is.na(g)) "" else g }),
      lapply(extra_fields, function(key) selectInput(ns(paste0(prefix, "_", key)), extra_labels[[key]] %||% key, cols, selected = guess(key)))
    )
  })
}

estimate_mr_set <- function(d, include_mode = FALSE, full = TRUE, ci_level = 0.95, run_presso = FALSE) {
  n_snp <- nrow(d)
  methods <- list()
  heterogeneity <- NULL
  pleiotropy <- NULL
  presso <- NULL
  alpha <- 1 - ci_level
  z <- stats::qnorm(1 - alpha / 2)

  if (n_snp == 1) {
    b <- d$beta.outcome[1] / d$beta.exposure[1]
    se <- abs(d$se.outcome[1] / d$beta.exposure[1])
    methods[["Wald ratio"]] <- c(estimate = b, se = se, ci_low = b - z * se, ci_high = b + z * se,
                                  p = 2 * stats::pnorm(-abs(b / se)))
    primary_method <- "Wald ratio"
  } else {
    mrobj <- MendelianRandomization::mr_input(
      bx = d$beta.exposure, bxse = d$se.exposure,
      by = d$beta.outcome, byse = d$se.outcome,
      snps = d$SNP, exposure = unique(d$gene)[1], outcome = "RA risk (Okada 2014 GWAS)"
    )
    ivw <- MendelianRandomization::mr_ivw(mrobj, model = "random", alpha = alpha)
    methods[["IVW"]] <- c(estimate = ivw@Estimate, se = ivw@StdError, ci_low = ivw@CILower, ci_high = ivw@CIUpper, p = ivw@Pvalue)
    primary_method <- "IVW"
    if (n_snp >= 3 && full) {
      med <- MendelianRandomization::mr_median(mrobj, alpha = alpha)
      methods[["Weighted median"]] <- c(estimate = med@Estimate, se = med@StdError, ci_low = med@CILower, ci_high = med@CIUpper, p = med@Pvalue)
      egg <- MendelianRandomization::mr_egger(mrobj, alpha = alpha)
      egger_df <- n_snp - 2
      egger_slope_p <- 2 * stats::pt(-abs(egg@Estimate / egg@StdError.Est), df = egger_df)
      egger_int_p <- 2 * stats::pt(-abs(egg@Intercept / egg@StdError.Int), df = egger_df)
      methods[["MR-Egger"]] <- c(estimate = egg@Estimate, se = egg@StdError.Est, ci_low = egg@CILower.Est, ci_high = egg@CIUpper.Est, p = egger_slope_p)
      heterogeneity <- list(Q = unname(ivw@Heter.Stat[1]), Q_df = n_snp - 1, Q_pval = unname(ivw@Heter.Stat[2]))
      pleiotropy <- list(intercept = egg@Intercept, se = egg@StdError.Int, p = egger_int_p, I2 = egg@I.sq)
      if (isTRUE(include_mode)) {
        mode_est <- tryCatch(MendelianRandomization::mr_mbe(mrobj, alpha = alpha), error = function(e) NULL)
        if (!is.null(mode_est)) {
          methods[["Weighted mode"]] <- c(estimate = mode_est@Estimate, se = mode_est@StdError,
                                           ci_low = mode_est@CILower, ci_high = mode_est@CIUpper, p = mode_est@Pvalue)
        }
      }
    }
    if (n_snp >= 4 && full && isTRUE(run_presso)) {
      presso <- tryCatch({
        pr <- MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure", data = d,
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
          NbDistribution = 1000, SignifThreshold = alpha, seed = 2024
        )
        mrr <- pr$`Main MR results`
        raw <- mrr[mrr$`MR Analysis` == "Raw", ]
        corr <- mrr[mrr$`MR Analysis` == "Outlier-corrected", ]
        list(
          global_p = pr$`MR-PRESSO results`$`Global Test`$Pvalue,
          raw_estimate = raw$`Causal Estimate`[1], raw_p = raw$`P-value`[1],
          corrected_estimate = if (!is.na(corr$`Causal Estimate`[1])) corr$`Causal Estimate`[1] else NA_real_,
          corrected_p = if (!is.na(corr$`P-value`[1])) corr$`P-value`[1] else NA_real_,
          n_outliers = length(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue[
            !is.na(suppressWarnings(as.numeric(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue))) &
            suppressWarnings(as.numeric(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue)) < alpha
          ])
        )
      }, error = function(e) list(error = conditionMessage(e)))
    }
  }

  res_table <- as.data.frame(do.call(rbind, methods))
  res_table <- data.frame(method = names(methods), res_table, row.names = NULL, check.names = FALSE)
  res_table$primary <- res_table$method == primary_method

  list(res_table = res_table, n_snp = n_snp, heterogeneity = heterogeneity, pleiotropy = pleiotropy,
       primary_method = primary_method, presso = presso, ci_level = ci_level)
}

load_mr_instrument_table <- function() {
  mr_obj <- readRDS(MR_PRIMARY_OBJECTS_RDS)

  dat <- as.data.frame(mr_obj$dat)
  dat$gene <- mr_obj$inst$gene[match(
    paste(dat$SNP, dat$id.exposure),
    paste(mr_obj$inst$SNP, mr_obj$inst$id.exposure)
  )]
  dat <- dat[dat$mr_keep, , drop = FALSE]

  inst_dt <- as.data.table(mr_obj$inst)
  dat_dt  <- as.data.table(mr_obj$dat)
  chk <- merge(dat_dt[, .(SNP, id.exposure, gene_cached = gene)],
               inst_dt[, .(SNP, id.exposure, gene_correct = gene)],
               by = c("SNP", "id.exposure"))
  mism <- chk[gene_cached != gene_correct]
  relabel_check <- list(n_snp = nrow(mism), genes = unique(c(mism$gene_cached, mism$gene_correct)))

  list(dat = dat, primary = as.data.frame(mr_obj$primary), relabel_check = relabel_check)
}

ARTHOMIX_COLORS <- list(
  blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
  magenta = "#e87ba4", violet = "#4a3aa7", red = "#e34948",
  ink = "#0b0b0b", ink_secondary = "#52514e", ink_muted = "#898781",
  grid = "#e1e0d9", axis = "#c3c2b7"
)

ARTHOMIX_STATUS <- list(good = "#0ca30c", warning = "#fab219", critical = "#d03b3b")

arthomix_pair <- function(levels) {
  pal <- c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$red, ARTHOMIX_COLORS$aqua,
            ARTHOMIX_COLORS$orange, ARTHOMIX_COLORS$violet, ARTHOMIX_COLORS$magenta,
            ARTHOMIX_COLORS$yellow)
  levels <- as.character(sort(unique(levels)))
  setNames(pal[seq_along(levels)], levels)
}

theme_arthomix <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = ARTHOMIX_COLORS$grid, linewidth = 0.3),
      axis.line = element_line(color = ARTHOMIX_COLORS$axis, linewidth = 0.3),
      axis.ticks = element_line(color = ARTHOMIX_COLORS$axis, linewidth = 0.3),
      axis.text = element_text(color = ARTHOMIX_COLORS$ink_secondary),
      axis.title = element_text(color = ARTHOMIX_COLORS$ink_secondary),
      plot.title = element_text(color = ARTHOMIX_COLORS$ink, size = base_size, hjust = 0),
      plot.subtitle = element_text(color = ARTHOMIX_COLORS$ink_muted, size = base_size - 2, hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(color = ARTHOMIX_COLORS$ink_secondary, size = base_size - 2),
      legend.text = element_text(color = ARTHOMIX_COLORS$ink_secondary, size = base_size - 2),
      strip.background = element_blank(),
      strip.text = element_text(color = ARTHOMIX_COLORS$ink, size = base_size - 1)
    )
}

compute_sample_qc <- function(expr, mad_k = 3, top_n_cor = 2000) {
  detect_cutoff <- stats::quantile(expr, 0.25, na.rm = TRUE)
  signal   <- colSums(expr, na.rm = TRUE)
  detected <- colSums(expr > detect_cutoff, na.rm = TRUE)

  gene_var <- apply(expr, 1, stats::var, na.rm = TRUE)
  top_idx  <- order(gene_var, decreasing = TRUE)[seq_len(min(top_n_cor, nrow(expr)))]
  sample_cor <- cor(expr[top_idx, , drop = FALSE], use = "pairwise.complete.obs")
  mean_cor <- (colSums(sample_cor, na.rm = TRUE) - 1) / (ncol(sample_cor) - 1)

  is_outlier <- function(x, low_only = FALSE) {
    m <- stats::median(x); s <- stats::mad(x)
    if (s == 0) return(rep(FALSE, length(x)))
    if (low_only) (m - x) > mad_k * s else abs(x - m) > mad_k * s
  }

  data.frame(
    sample = colnames(expr),
    signal = signal,
    detected = detected,
    mean_cor = mean_cor,
    flag_signal = is_outlier(signal),
    flag_detected = is_outlier(detected),
    flag_cor = is_outlier(mean_cor, low_only = TRUE),
    row.names = NULL
  )
}

summarize_norm_diagnostics <- function(m) {
  sample_medians <- apply(m, 2, stats::median, na.rm = TRUE)
  sample_iqr <- apply(m, 2, stats::IQR, na.rm = TRUE)
  data.frame(
    n_samples = ncol(m), n_genes = nrow(m),
    max_value = max(m, na.rm = TRUE), min_value = min(m, na.rm = TRUE),
    median_sd = stats::sd(sample_medians, na.rm = TRUE), iqr_sd = stats::sd(sample_iqr, na.rm = TRUE),
    median_range = diff(range(sample_medians, na.rm = TRUE)),
    iqr_range = diff(range(sample_iqr, na.rm = TRUE))
  )
}

needs_quantile_norm <- function(diag) {
  diag$max_value > 100 || diag$median_sd > 0.5 || diag$iqr_sd > 0.5
}

expr_raw_health <- function(expr) {
  m <- as.matrix(expr)
  vals <- as.numeric(m)
  list(
    n_features = nrow(m), n_samples = ncol(m),
    n_missing = sum(is.na(vals)),
    n_zero = sum(vals == 0, na.rm = TRUE),
    n_infinite = sum(is.infinite(vals)),
    n_duplicated_features = sum(duplicated(rownames(m)))
  )
}

detect_expr_data_type <- function(expr) {
  m <- as.matrix(expr)
  finite_vals <- m[is.finite(m)]
  if (length(finite_vals) == 0) return("expression")
  has_negative <- any(finite_vals < 0)
  looks_like_counts <- !has_negative &&
    mean(abs(finite_vals - round(finite_vals)) < 1e-6) > 0.6 &&
    max(finite_vals) > 100
  if (looks_like_counts) return("counts")
  diag <- summarize_norm_diagnostics(m)
  if (has_negative || !needs_quantile_norm(diag)) return("already_normalised")
  "expression"
}

filter_and_transform_expr <- function(expr, min_expr = 0, min_sample_frac = 0,
                                        max_na_pct = 100, log2_transform = FALSE) {
  m <- as.matrix(expr)
  m[is.infinite(m)] <- NA
  keep_na <- (rowMeans(is.na(m)) * 100) <= max_na_pct
  above <- m >= min_expr
  above[is.na(above)] <- FALSE
  keep_expr <- (rowMeans(above) * 100) >= min_sample_frac
  m <- m[keep_na & keep_expr, , drop = FALSE]
  if (anyNA(m)) {
    row_med <- apply(m, 1, stats::median, na.rm = TRUE)
    na_idx <- which(is.na(m), arr.ind = TRUE)
    m[na_idx] <- row_med[na_idx[, 1]]
  }
  if (isTRUE(log2_transform)) {
    m[m <= 0] <- NA
    m <- log2(m)
    m <- m[stats::complete.cases(m), , drop = FALSE]
  }
  m
}

apply_chosen_norm <- function(expr, method = c("none", "quantile", "tmm")) {
  method <- match.arg(method)
  m <- as.matrix(expr)
  if (identical(method, "none")) {
    return(list(expr = m, label = "None - used as filtered/transformed, no normalisation applied"))
  }
  if (identical(method, "tmm")) {
    validate(need(all(m >= 0, na.rm = TRUE),
                  "TMM normalisation expects raw, non-negative counts, but this data has negative values."))
    counts <- round(m)
    storage.mode(counts) <- "integer"
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")
    return(list(expr = edgeR::cpm(dge, log = TRUE, prior.count = 1),
                label = "TMM (edgeR::calcNormFactors) plus log2-CPM"))
  }
  mtx <- limma::normalizeBetweenArrays(m, method = "quantile")
  rownames(mtx) <- rownames(m); colnames(mtx) <- colnames(m)
  list(expr = mtx, label = "Quantile normalisation (limma::normalizeBetweenArrays)")
}

pca_of <- function(m, n_pc = 5) {
  m <- m[apply(m, 1, stats::sd) > 0, , drop = FALSE]
  p <- prcomp(t(m), scale. = TRUE)
  var_exp <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  n_pc <- min(n_pc, ncol(p$x))
  df <- as.data.frame(p$x[, seq_len(n_pc), drop = FALSE])
  colnames(df) <- paste0("PC", seq_len(n_pc))
  df$sample <- rownames(df)
  list(df = df, var_exp = var_exp, n_pc = n_pc)
}

scree_plot <- function(var_exp, max_pc = 10) {
  n <- min(max_pc, length(var_exp))
  df <- data.frame(pc = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))),
                    pct = var_exp[seq_len(n)])
  ggplot(df, aes(x = pc, y = pct)) +
    geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.6) +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.5, size = 3, color = ARTHOMIX_COLORS$ink_secondary) +
    labs(x = NULL, y = "% variance explained") +
    theme_arthomix() +
    theme(panel.grid.major.x = element_blank())
}

plot_pca_advanced <- function(pca_obj, meta, color_by, pc_x = 1, pc_y = 2,
                                title_suffix = NULL, show_ellipse = TRUE, show_labels = FALSE) {
  xcol <- paste0("PC", pc_x); ycol <- paste0("PC", pc_y)
  df <- merge(pca_obj$df[, c("sample", xcol, ycol)], meta[, c("sample", color_by)], by = "sample")
  xlab <- paste0(xcol, " (", pca_obj$var_exp[pc_x], "%)")
  ylab <- paste0(ycol, " (", pca_obj$var_exp[pc_y], "%)")
  is_num <- is.numeric(df[[color_by]])

  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = .data[[color_by]])) +
    geom_point(size = 2.2, alpha = 0.85) +
    labs(x = xlab, y = ylab, color = color_by, title = title_suffix) +
    theme_arthomix()

  if (!is_num && show_ellipse) {
    ok_groups <- names(which(table(df[[color_by]]) >= 4))
    if (length(ok_groups) > 0) {
      p <- p + stat_ellipse(data = df[df[[color_by]] %in% ok_groups, , drop = FALSE],
                              type = "norm", level = 0.68, linewidth = 0.4, show.legend = FALSE)
    }
  }
  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 2.6, show.legend = FALSE, max.overlaps = 15)
  }

  if (is_num) p + scale_color_gradient(low = "#cde2fb", high = ARTHOMIX_COLORS$blue)
  else p + scale_color_manual(values = arthomix_pair(df[[color_by]]))
}

overlap_region_sizes <- function(sets, max_regions = NULL) {
  nms <- names(sets)
  combos <- unlist(lapply(seq_along(nms), function(k) utils::combn(nms, k, simplify = FALSE)), recursive = FALSE)
  sizes <- vapply(combos, function(cmb) {
    inc <- Reduce(intersect, sets[cmb])
    excl_from <- setdiff(nms, cmb)
    if (length(excl_from) > 0) inc <- setdiff(inc, Reduce(union, sets[excl_from]))
    length(inc)
  }, numeric(1))
  df <- data.frame(combination = vapply(combos, paste, character(1), collapse = " & "), n_features = sizes)
  df <- df[df$n_features > 0, , drop = FALSE]
  df <- df[order(-df$n_features, df$combination), , drop = FALSE]
  if (!is.null(max_regions) && nrow(df) > max_regions) df <- df[seq_len(max_regions), ]
  rownames(df) <- NULL
  df
}

draw_overlap_venn <- function(sets, title = NULL, max_bars = 30, fill_low = "#EAF3FB", fill_high = ARTHOMIX_COLORS$blue) {
  sets <- sets[lengths(sets) > 0]
  validate(need(length(sets) >= 2, "Need at least two non-empty datasets to compare."))

  common_n <- length(Reduce(intersect, sets))
  auto_title <- title %||% sprintf("Common to all %d datasets: %s features", length(sets), format(common_n, big.mark = ","))
  pal <- c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$orange, ARTHOMIX_COLORS$aqua,
            ARTHOMIX_COLORS$violet, ARTHOMIX_COLORS$magenta, ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)

  if (length(sets) <= 7) {
    p <- ggVennDiagram::ggVennDiagram(
      sets, label = "both", label_alpha = 0, label_size = 3.3, label_percent_digit = 1,
      edge_size = 0.9, set_size = 4.2, set_color = pal[seq_along(sets)]
    )
    two_set_asymmetric <- length(sets) == 2
    if (two_set_asymmetric) {
      for (i in seq_along(p$layers)) {
        l <- p$layers[[i]]
        if (inherits(l$geom, "GeomText") && !is.null(l$mapping$label) &&
            grepl("name", rlang::as_label(l$mapping$label), fixed = TRUE)) {
          p$layers[[i]]$aes_params$hjust <- 1
        }
      }
    }
    p <- p +
      scale_fill_gradient(low = fill_low, high = fill_high, name = "Features") +
      scale_x_continuous(expand = expansion(mult = 0.12)) +
      scale_y_continuous(expand = expansion(mult = 0.12)) +
      labs(title = auto_title) +
      theme(legend.position = "right", legend.text = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 13, hjust = 0))
    p$coordinates <- coord_equal(clip = "off")
    if (two_set_asymmetric) {
      p <- p + theme(plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5, l = 130, unit = "pt"))
    }
    p
  } else {
    df <- overlap_region_sizes(sets, max_regions = max_bars)
    ggplot(df, aes(x = reorder(combination, n_features), y = n_features)) +
      geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.7) +
      geom_text(aes(label = format(n_features, big.mark = ",")), hjust = -0.15, size = 3, color = ARTHOMIX_COLORS$ink_secondary) +
      coord_flip(clip = "off") +
      labs(x = NULL, y = "Features in exactly this combination of datasets", title = auto_title) +
      theme_arthomix() +
      theme(panel.grid.major.y = element_blank())
  }
}

venn_region_at_point <- function(sets, x, y) {
  sets <- sets[lengths(sets) > 0]
  if (length(sets) < 2 || length(sets) > 7 || is.null(x) || is.null(y)) return(NULL)
  d <- tryCatch(ggVennDiagram::process_data(ggVennDiagram::Venn(sets)), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  edge <- d$regionEdge
  label <- d$regionLabel

  point_in_polygon <- function(px, py, poly_x, poly_y) {
    n <- length(poly_x)
    inside <- FALSE
    j <- n
    for (i in seq_len(n)) {
      if (((poly_y[i] > py) != (poly_y[j] > py)) &&
          (px < (poly_x[j] - poly_x[i]) * (py - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])) {
        inside <- !inside
      }
      j <- i
    }
    inside
  }

  ids <- unique(edge$id)
  ids <- ids[order(-lengths(strsplit(ids, "/", fixed = TRUE)))]
  for (rid in ids) {
    seg <- edge[edge$id == rid, , drop = FALSE]
    if (nrow(seg) < 3 || !point_in_polygon(x, y, seg$X, seg$Y)) next
    lab <- label[label$id == rid, , drop = FALSE]
    if (nrow(lab) == 0) next
    return(list(name = lab$name[1], item = lab$item[[1]], count = lab$count[1]))
  }
  NULL
}

qc_bar_plot <- function(df, y_col, flag_col, y_label, subtitle = NULL) {
  df$.flag <- ifelse(df[[flag_col]], "Flagged", "Within range")
  ggplot(df, aes(x = reorder(sample, .data[[y_col]]), y = .data[[y_col]], fill = .flag)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("Within range" = ARTHOMIX_COLORS$blue, "Flagged" = ARTHOMIX_STATUS$critical)) +
    labs(x = NULL, y = y_label, fill = NULL, subtitle = subtitle) +
    theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
}

MODULE_REGISTRY <- list(
  list(
    id = "transcriptomics", tab = "transcriptomics",
    title = "Transcriptomics",
    tagline = "Differential expression, WGCNA, Mendelian randomisation and colocalisation, feature selection,machine learning modelling and cross-tissue validation, cross-ancestry validation and biomarker card",
    icon = "dna", status = "available", kind = "Single-omics"
  ),
  list(
    id = "methylomics", tab = "methylomics",
    title = "Methylomics",
    tagline = "Differential methylation, WGCNA, Mendelian randomisation and colocalisation, feature selection,machine learning modelling and external validation and biomarker card",
    icon = "circle-nodes", status = "available", kind = "Single-omics"
  ),
  list(
    id = "crossomics", tab = "crossomics",
    title = "Cross-Omics",
    tagline = "Perform the cross-omics integration.",
    icon = "arrows-left-right", status = "available", kind = "Multi-omics"
  ),
  list(
    id = "multiomics", tab = "multiomics",
    title = "Multi-Omics",
    tagline = "Integrate the two omics data.",
    icon = "layer-group", status = "available", kind = "Multi-omics"
  ),
  list(
    id = "arthochat", tab = "arthochat",
    title = "ArthOChat",
    tagline = "Ask anything about your dataset, your results so far, or the underlying methodology, grounded in a live PubMed search with citations. Automatically scoped to whichever module/sub-module you have open, with other modules available on request.",
    icon = "comments", status = "available", kind = "Assistant"
  )
)

ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT <- "document.getElementById('arthochat_drawer').classList.add('open')"
ARTHOCHAT_DRAWER_OPEN_JS <- paste0(ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT, "; return false;")

arthochat_shortcut_ui <- function(hint = NULL, compact = FALSE) {
  jump_link <- tags$a(
    "Open ArthOChat", href = "#", class = paste("btn btn-primary btn-sm", if (compact) "btn-block"),
    onclick = ARTHOCHAT_DRAWER_OPEN_JS
  )
  if (compact) {
    return(div(
      class = "arthochat-shortcut arthochat-shortcut-compact",
      div(icon("comments"), strong(" Ask ArthOChat")),
      p(hint %||% "Have a methodology question? Ask ArthOChat."),
      jump_link
    ))
  }
  div(
    class = "arthochat-shortcut",
    icon("comments"),
    span(
      strong("Have a methodology question? "),
      hint %||% "Ask ArthOChat. It can see this dataset and cites PubMed for any scientific claim.",
      " "
    ),
    jump_link
  )
}

