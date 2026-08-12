## R/mod_preprocessing.R
## Submodule: Preprocessing and Batch Correction (Section 2.2)
## Everything on this page runs live, on whatever data is loaded: each
## source dataset is cleaned individually first (log2 scale check, sample
## filters), then datasets are merged onto their shared genes/probes
## (Venn/overlap diagram), then the merged cohort is normalised and
## batch-corrected (ComBat, limma::removeBatchEffect, or ComBat-seq). Three
## tabs, in that order:
##   "Preprocessing"    - per-dataset upload/preloaded-GEO-dataset/currently
##                         loaded dataset, plus filters.
##   "Merge datasets"    - overlap of features across every configured
##                         dataset, and the merge itself.
##   "Batch correction"  - normalisation + batch correction on the merged
##                         (or single) dataset, with a full set of tunable
##                         filters.
## The analysis adapts to whatever columns/scale/platform(s) you give it
## rather than assuming any particular dataset's layout. Picking "A
## preloaded GEO dataset" for two sources (e.g. GSE93272 and GSE110169)
## reproduces the example cohort from scratch, computed live rather than
## read from a precomputed file.

MAX_PP_SOURCES <- 6

mod_preprocessing_config <- list(
  id = "preprocessing", group = "Data",
  title = "Preprocessing and Batch Correction",
  description = "Preprocess and merge one or more datasets on their shared genes and probes, then normalise and batch-correct - live, on whatever you load.",
  icon = "filter"
)

## ---------------------------------------------------------------------------
## One data source: upload + column mapping + per-dataset filters. Instantiated
## MAX_PP_SOURCES times up front (server side) so state survives changing how
## many are shown; the UI only reveals the first N via conditionalPanel, which
## keeps every already-configured block's inputs intact when N changes.
## ---------------------------------------------------------------------------

## Small SaaS-style hover tooltip: a circular info icon that reveals a dark
## floating card of `text` on hover/focus, pure CSS (see .field-hint* in
## www/custom.css) - no JS tooltip library needed.
mod_pp_field_hint <- function(text) {
  tags$span(class = "field-hint", tabindex = "0",
            icon("circle-info"),
            tags$span(class = "field-hint-box", text))
}

## Same name-based guess mod_pp_source_server's own colmap uses internally
## (see its comment) - pulled out to module scope so the batch-upload path
## below can reuse it without a UI-bound copy per dataset.
pp_guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
  hit <- cols[tolower(cols) %in% tolower(exact)]
  if (length(hit) > 0) return(hit[1])
  hit <- cols[grepl(paste(contains, collapse = "|"), cols, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[1])
  fallback
}

## Collapses a probe-level expression matrix to one row per gene, given an
## external probe -> gene symbol annotation table (two columns, auto-
## guessed by name like everything else in this file). Three interchangeable
## methods, since published methods vary: "median" (per gene, per sample,
## the median across every probe mapping to that gene - e.g. Zhu et al.
## 2021's stated method), "maxmean" (one representative probe per gene -
## the one with the highest mean expression - ArthOMix's own preloaded-GEO-
## dataset convention, global.R::collapse_probes_to_genes()), or "mean"
## (per gene, per sample, the mean across probes). Probes with no
## annotation match, or mapped to more than one gene symbol (e.g.
## "MIR4640///DDR1"), are dropped either way - never split or duplicated.
pp_collapse_probes_to_genes <- function(expr, annot, method = c("median", "maxmean", "mean")) {
  method <- match.arg(method)
  cols <- colnames(annot)
  probe_col <- pp_guess_col(cols, c("probe", "probe_id", "probeid", "probeset", "probesetid", "id"))
  gene_col  <- pp_guess_col(cols, c("gene_symbol", "genesymbol", "symbol", "gene"),
                             fallback = if (length(cols) >= 2) cols[2] else cols[1])
  map <- stats::setNames(as.character(annot[[gene_col]]), as.character(annot[[probe_col]]))
  sym <- unname(map[rownames(expr)])
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym, fixed = TRUE)
  validate(need(any(keep), "None of this expression matrix's row IDs matched the annotation file's probe-ID column, or every match was to more than one gene. Check the annotation file's columns."))
  ex <- expr[keep, , drop = FALSE]; sym <- sym[keep]

  if (method == "maxmean") {
    row_mean <- rowMeans(ex, na.rm = TRUE)
    keep_idx <- tapply(seq_along(sym), sym, function(idx) idx[which.max(row_mean[idx])])
    out <- ex[unlist(keep_idx), , drop = FALSE]
    rownames(out) <- names(keep_idx)
    return(out[order(rownames(out)), , drop = FALSE])
  }
  agg_fn <- if (method == "median") stats::median else base::mean
  apply(ex, 2, function(col_vals) tapply(col_vals, sym, agg_fn, na.rm = TRUE))
}

## Display names for this tab only. mod_dataset.R's own dropdown and
## Overview and Datasets keep showing the raw GEO accession (full
## traceability lives there); this tab is a working cohort picker for
## someone running the pipeline, so it reads as tissue + role rather than
## an accession lookup. Add an entry here if a new bundled source is ever
## added to GEO_SOURCES (global.R) - anything missing just falls back to
## its raw ID (a plain character vector's `[[` errors on an unknown name
## rather than returning NULL, so the fallback is a `%in%` check, not `%||%`).
PP_COHORT_LABELS <- c(
  "GSE93272"  = "Whole Blood Training Cohort A",
  "GSE110169" = "Whole Blood Training Cohort B",
  "GSE15573"  = "PBMC Validation Cohort",
  "GSE89408"  = "Synovial Tissue Validation Cohort"
)
PP_MERGED_COHORT_LABEL <- "Whole Blood Training Cohort (Merged)"

pp_cohort_label <- function(id) {
  if (identical(id, default_dataset_entry$id)) return(PP_MERGED_COHORT_LABEL)
  if (id %in% names(PP_COHORT_LABELS)) unname(PP_COHORT_LABELS[[id]]) else id
}

## Same id -> value pairing as mod_dataset.R's own preloaded_choices(), just
## with this tab's professional display names instead of the raw-accession
## ones.
pp_cohort_choices <- function() {
  ids <- unname(preloaded_choices())
  stats::setNames(ids, vapply(ids, pp_cohort_label, character(1)))
}

## Reads one dataset for the Preprocessing tab's "Preloaded Data" box -
## either one of the app's bundled raw GEO sources (same read path as
## mod_dataset.R's PRELOADED_DATASETS / mod_pp_source_server's own
## "preloaded" source_type: GSE89408 is already gene-level RNA-seq counts,
## every other source is microarray and gets collapsed to gene symbol via
## get_collapsed_genes()), or whatever is currently loaded via the Dataset
## tab (choice_id == "__current__", read from the shared `dataset`
## reactiveValues). Returns the same list(label=, expr=, meta=, ...) shape
## every dataset entry feeding merge_inputs() uses.
pp_preloaded_read <- function(choice_id, log2_choice, dataset = NULL) {
  if (identical(choice_id, "__current__")) {
    validate(need(!is.null(dataset) && !is.null(dataset$expr),
                  "No dataset is currently loaded. Pick one on the Dataset tab first."))
    expr <- dataset$expr
    meta <- dataset$meta
    label <- dataset$source %||% "Currently Loaded Dataset"
  } else {
    gse <- choice_id
    if (identical(gse, "GSE89408")) {
      d <- load_individual_dataset(gse)
      validate(need(!is.null(d), paste("Raw data for", gse, "was not found on disk.")))
      expr <- d$expr; meta <- d$meta
    } else {
      eset <- get_raw_eset(gse)
      validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
      expr <- get_collapsed_genes(gse)
      meta <- eset_harmonize_meta(eset, gse)
      keep <- !is.na(meta$group)
      meta <- meta[keep, , drop = FALSE]
      expr <- expr[, meta$sample, drop = FALSE]
    }
    label <- pp_cohort_label(gse)
  }
  if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_

  n_samples_before <- ncol(expr); n_genes_before <- nrow(expr)
  q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
  needs_log <- if (identical(log2_choice, "force")) TRUE
               else if (identical(log2_choice, "skip")) FALSE
               else isTRUE(!is.na(q99) && q99 > 100)
  if (needs_log) {
    expr[expr <= 0] <- NA
    expr <- log2(expr)
    expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  }

  list(label = label, expr = as.matrix(expr), meta = meta,
       n_samples_before = n_samples_before, n_samples_after = ncol(expr),
       n_genes_before = n_genes_before, n_genes_after = nrow(expr),
       log2_applied = needs_log)
}

mod_pp_source_ui <- function(id, default_gse = NULL, n_sources_id = NULL) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("source_type_ui")),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'current'", ns("source_type")),
      uiOutput(ns("current_note"))
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'preloaded'", ns("source_type")),
      if (!is.null(default_gse)) {
        ## This slot has a fixed training-cohort default (Dataset 1 ->
        ## GSE93272, Dataset 2 -> GSE110169) - no dropdown to pick a
        ## different one, just load it. `preloaded_choice` is still a real
        ## selectInput (shinyjs::hidden(), not absent) so raw_pair()'s
        ## req(input$preloaded_choice) keeps working unchanged.
        tagList(
          div(class = "empty-note", icon("check"),
              sprintf("Using %s.", pp_cohort_label(default_gse))),
          shinyjs::hidden(selectInput(ns("preloaded_choice"), "Choose a cohort", choices = pp_cohort_choices(),
                                       selected = default_gse, selectize = FALSE))
        )
      } else {
        selectInput(ns("preloaded_choice"), "Choose a cohort", choices = pp_cohort_choices(),
                    selected = character(0), selectize = FALSE)
      },
      p(class = "empty-note", icon("triangle-exclamation"),
        "Raw, single-platform data - not yet merged or normalised. Sample, group, and sex are mapped automatically.")
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("source_type")),
      if (!is.null(n_sources_id)) {
        conditionalPanel(
          condition = sprintf("input['%s'] > 1", n_sources_id),
          p(class = "empty-note", icon("circle-info"),
            "Combining more than one dataset: each one needs its own expression matrix + sample metadata pair, uploaded separately in its own box below. Every dataset's feature IDs must use the same identifier type across the board (e.g. all gene symbols, or all the same probe IDs) - the Merge step keeps only features common to every dataset, and fails below 20 shared features. Each dataset's Sample ID column should also be unique within that dataset.")
        )
      },
      textInput(ns("label"), "Label for this dataset", value = "", placeholder = "e.g. GSE12345"),
      div(class = "field-label-with-hint", span("Expression matrix"),
          mod_pp_field_hint("CSV or RDS file: genes or probes in rows, samples in columns. For CSV files, the first column must contain the feature ID.")),
      fileInput(ns("expr_file"), NULL, accept = c(".csv", ".rds", ".Rds")),
      div(class = "field-label-with-hint", span("Sample metadata"),
          mod_pp_field_hint("CSV or RDS file: a data frame with one row per sample.")),
      fileInput(ns("meta_file"), NULL, accept = c(".csv", ".rds", ".Rds")),
      uiOutput(ns("upload_preview_ui")),
      uiOutput(ns("colmap"))
    ),
    tags$hr(),
    div(class = "filter-section-header", icon("users"), "Sample filters"),
    uiOutput(ns("group_filter_ui")),
    uiOutput(ns("numeric_filter_ui")),
    checkboxInput(ns("dedup"), "Deduplicate samples by an ID column (keep first occurrence)", value = FALSE),
    conditionalPanel(condition = sprintf("input['%s']", ns("dedup")), uiOutput(ns("dedup_col_ui"))),
    tags$hr(),
    div(class = "filter-section-header", icon("dna"), "Feature filters"),
    uiOutput(ns("log2_ui")),
    sliderInput(ns("max_na_pct"), "Drop features with more than this % of samples missing (remaining gaps are median-imputed)",
                min = 0, max = 80, value = 0, step = 5),
    textInput(ns("exclude_pattern"), "Exclude features whose ID matches this pattern (regex, optional)",
              value = "", placeholder = "e.g. ^AFFX for Affymetrix control probes"),
    actionButton(ns("run"), "Preprocess this dataset", icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "margin-top:8px;", uiOutput(ns("status_ui")))
  )
}

mod_pp_source_server <- function(id, default_label = "Dataset", default_gse = NULL, dataset = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$source_type_ui <- renderUI({
      type_choices <- c("Upload files" = "upload", "A bundled cohort" = "preloaded",
                         "Currently loaded dataset" = "current")
      current <- isolate(input$source_type)
      selected <- if (!is.null(current) && current %in% type_choices) {
        current
      } else if (!is.null(default_gse)) {
        "preloaded"
      } else {
        "upload"
      }
      radioButtons(ns("source_type"), "Data source", choices = type_choices, selected = selected)
    })

    source_type <- reactive(input$source_type %||% "upload")
    use_preloaded <- reactive(identical(source_type(), "preloaded"))
    use_upload    <- reactive(identical(source_type(), "upload"))
    use_current   <- reactive(identical(source_type(), "current"))

    ## "Currently loaded dataset" reuses whatever the Dataset tab (or another
    ## sub-module that writes to the shared `dataset` reactiveValues) already
    ## has loaded - the same object the sidebar's "Current dataset" panel and
    ## every other sub-module read from. No file to parse, no GSE to fetch.
    output$current_note <- renderUI({
      if (is.null(dataset) || is.null(dataset$expr)) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    "Nothing is currently loaded. Pick a dataset on the Dataset tab first, then come back here."))
      }
      div(class = "empty-note", icon("check"),
          sprintf("Using %s: %s genes x %s samples.", dataset$source,
                   format(nrow(dataset$expr), big.mark = ","), ncol(dataset$expr)))
    })

    meta_raw <- reactive({
      req(use_upload(), input$meta_file)
      path <- input$meta_file$datapath
      if (grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
        validate(need(is.data.frame(d), "The uploaded metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    ## Fires the moment both files are selected and readable, before the
    ## column-mapping dropdowns even appear - otherwise a large upload can
    ## look like it silently did nothing while it's actually just parsing.
    expr_raw_preview <- reactive({
      req(use_upload(), input$expr_file)
      if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
        readRDS(input$expr_file$datapath)
      } else {
        m <- as.data.frame(data.table::fread(input$expr_file$datapath, showProgress = FALSE))
        as.matrix(m[, -1, drop = FALSE])
      }
    })

    output$upload_preview_ui <- renderUI({
      req(use_upload(), input$expr_file, input$meta_file)
      preview <- tryCatch(list(expr = expr_raw_preview(), meta = meta_raw()), error = function(e) e)
      if (inherits(preview, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not read the uploaded file(s):", conditionMessage(preview))))
      }
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s features x %s samples. Read %s: %s rows. Map the columns below.",
                  input$expr_file$name, format(nrow(preview$expr), big.mark = ","), ncol(preview$expr),
                  input$meta_file$name, nrow(preview$meta)))
    })

    ## Same name-based guess as mod_dataset.R's upload form - see its
    ## comment for why a plain selectInput() default (first column in the
    ## file) is actively wrong here, not just imprecise.
    guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
      hit <- cols[tolower(cols) %in% tolower(exact)]
      if (length(hit) > 0) return(hit[1])
      hit <- cols[grepl(paste(contains, collapse = "|"), cols, ignore.case = TRUE)]
      if (length(hit) > 0) return(hit[1])
      fallback
    }

    output$colmap <- renderUI({
      req(use_upload(), input$meta_file)
      cols <- colnames(meta_raw())
      tagList(
        selectInput(ns("map_id"), "Sample ID column", choices = cols,
                    selected = guess_col(cols, c("sample", "sample_id", "id", "geo_accession", "accession")),
                    selectize = FALSE),
        selectInput(ns("map_group"), "Group / diagnosis column", choices = cols,
                    selected = guess_col(cols, c("group", "diagnosis", "disease", "condition", "status", "phenotype")),
                    selectize = FALSE),
        selectInput(ns("map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "cohort", "platform", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })

    ## Loaded + column-mapped, before any filtering.
    raw_pair <- reactive({
      if (use_current()) {
        req(dataset, dataset$expr, dataset$meta)
        expr <- dataset$expr; meta <- dataset$meta
        if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_
        list(expr = expr, meta = meta, label = dataset$source %||% default_label)
      } else if (use_preloaded()) {
        req(input$preloaded_choice)
        gse <- input$preloaded_choice
        if (identical(gse, "GSE89408")) {
          ## RNA-seq counts, already gene-level - no probe/platform to collapse.
          d <- load_individual_dataset(gse)
          validate(need(!is.null(d), paste("Raw data for", gse, "was not found on disk.")))
          expr <- d$expr; meta <- d$meta
        } else {
          ## Microarray: collapse probes to gene symbol first, so this
          ## dataset can actually share features with one on a different
          ## platform when merged (see get_collapsed_genes()/
          ## collapse_probes_to_genes() in global.R - cached per GSE ID,
          ## since this is a genuinely slow step).
          eset <- get_raw_eset(gse)
          validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
          expr <- get_collapsed_genes(gse)
          meta <- eset_harmonize_meta(eset, gse)
          keep <- !is.na(meta$group)
          meta <- meta[keep, , drop = FALSE]
          expr <- expr[, meta$sample, drop = FALSE]
        }
        if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_
        list(expr = expr, meta = meta, label = pp_cohort_label(gse))
      } else {
        req(input$expr_file, input$meta_file, input$map_id, input$map_group)
        path <- input$expr_file$datapath
        expr <- if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
          readRDS(path)
        } else {
          m <- as.data.frame(data.table::fread(path, showProgress = FALSE))
          rn <- as.character(m[[1]])
          m <- as.matrix(m[, -1, drop = FALSE])
          rownames(m) <- rn
          m
        }
        meta <- meta_raw()
        meta$sample <- as.character(meta[[input$map_id]])
        meta$group  <- as.character(meta[[input$map_group]])
        meta$sex    <- if (!identical(input$map_sex %||% "(none)", "(none)")) as.character(meta[[input$map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$map_batch %||% "(none)", "(none)")) as.character(meta[[input$map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(length(common) >= 3, "Fewer than 3 sample IDs in the expression matrix match the metadata sample-ID column. Check the column mapping."))
        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        label <- if (nzchar(trimws(input$label %||% ""))) input$label else default_label
        list(expr = expr, meta = meta, label = label)
      }
    })

    output$group_filter_ui <- renderUI({
      pair <- raw_pair()
      cols <- colnames(pair$meta)
      tagList(
        selectInput(ns("filter_col"), "Keep only samples where", choices = c("(no filter)", cols), selectize = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] != '(no filter)'", ns("filter_col")),
          uiOutput(ns("filter_val_ui"))
        )
      )
    })

    output$filter_val_ui <- renderUI({
      req(input$filter_col, !identical(input$filter_col, "(no filter)"))
      pair <- raw_pair()
      vals <- sort(unique(stats::na.omit(as.character(pair$meta[[input$filter_col]]))))
      checkboxGroupInput(ns("filter_vals"), "...equals one of", choices = vals, selected = vals, inline = TRUE)
    })

    output$numeric_filter_ui <- renderUI({
      pair <- raw_pair()
      num_cols <- names(pair$meta)[vapply(pair$meta, is.numeric, logical(1))]
      req(length(num_cols) > 0)
      tagList(
        selectInput(ns("num_filter_col"), "...and within this numeric range (optional, for example age or RIN)",
                    choices = c("(no filter)", num_cols), selectize = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] != '(no filter)'", ns("num_filter_col")),
          uiOutput(ns("num_filter_range_ui"))
        )
      )
    })

    output$num_filter_range_ui <- renderUI({
      req(input$num_filter_col, !identical(input$num_filter_col, "(no filter)"))
      pair <- raw_pair()
      v <- suppressWarnings(as.numeric(pair$meta[[input$num_filter_col]]))
      rng <- suppressWarnings(range(v, na.rm = TRUE))
      validate(need(all(is.finite(rng)), "This column has no numeric values to filter on."))
      sliderInput(ns("num_filter_range"), "...between", min = floor(rng[1]), max = ceiling(rng[2]), value = rng)
    })

    output$dedup_col_ui <- renderUI({
      pair <- raw_pair()
      selectInput(ns("dedup_col"), "Deduplicate by", choices = colnames(pair$meta), selectize = FALSE)
    })

    ## Two pipelines, one control: same three choices either way, but the
    ## PRELOADED path's default stays "auto" exactly as already built -
    ## untouched. The UPLOAD path defaults to "skip" instead, because
    ## "auto" only looks at value magnitude (q99 > 100), which can't tell a
    ## large-valued RAW RNA-seq count matrix (correct DESeq2 input, should
    ## stay unlogged) apart from large-valued already-normalised data that
    ## genuinely needs logging - a real upload of raw counts was getting
    ## auto-log2'd here and then rejected by DESeq2 downstream for having
    ## negative/non-integer values. A user can still override either
    ## default by picking a different radio option themselves.
    output$log2_ui <- renderUI({
      default <- if (use_upload()) "skip" else "auto"
      radioButtons(ns("log2"), "Log2 transform",
                   choiceNames = list("Auto-detect (recommended for preloaded/normalised data)", "Force log2",
                                       "Skip (recommended for raw RNA-seq counts, e.g. for DESeq2)"),
                   choiceValues = list("auto", "force", "skip"), selected = default)
    })

    result <- eventReactive(input$run, {
      pair <- raw_pair()
      expr <- pair$expr; meta <- pair$meta; label <- pair$label

      n_samples_before <- ncol(expr); n_genes_before <- nrow(expr)

      if (!is.null(input$filter_col) && !identical(input$filter_col, "(no filter)") && length(input$filter_vals) > 0) {
        keep <- as.character(meta[[input$filter_col]]) %in% input$filter_vals
        meta <- meta[keep, , drop = FALSE]
        expr <- expr[, meta$sample, drop = FALSE]
      }
      if (!is.null(input$num_filter_col) && !identical(input$num_filter_col, "(no filter)") && !is.null(input$num_filter_range)) {
        v <- suppressWarnings(as.numeric(meta[[input$num_filter_col]]))
        keep <- !is.na(v) & v >= input$num_filter_range[1] & v <= input$num_filter_range[2]
        meta <- meta[keep, , drop = FALSE]
        expr <- expr[, meta$sample, drop = FALSE]
      }
      validate(need(ncol(expr) >= 3, "Fewer than 3 samples remain after the sample filters."))

      if (isTRUE(input$dedup) && !is.null(input$dedup_col) && input$dedup_col %in% colnames(meta)) {
        keep <- !duplicated(meta[[input$dedup_col]])
        meta <- meta[keep, , drop = FALSE]
        expr <- expr[, meta$sample, drop = FALSE]
      }

      pattern <- trimws(input$exclude_pattern %||% "")
      if (nzchar(pattern)) {
        matched <- tryCatch(grepl(pattern, rownames(expr), perl = TRUE), error = function(e) NULL)
        validate(need(!is.null(matched), paste("Invalid regex pattern:", pattern)))
        expr <- expr[!matched, , drop = FALSE]
        validate(need(nrow(expr) > 0, "The exclusion pattern matched every feature. Check the regular expression and try again."))
      }

      ## missing-data tolerance: drop features missing in more than the
      ## chosen % of samples, median-impute whatever gaps remain in the rest
      ## - a real (if simple) missing-data technique, and strictly more
      ## permissive than the old hard complete.cases() when the slider is
      ## above 0, identical to it at the default of 0.
      na_pct <- rowMeans(is.na(expr)) * 100
      expr <- expr[na_pct <= (input$max_na_pct %||% 0), , drop = FALSE]
      validate(need(nrow(expr) > 0, "No features remain within the missing-data tolerance. Raise the missing-data slider and try again."))
      if (anyNA(expr)) {
        row_med <- apply(expr, 1, stats::median, na.rm = TRUE)
        na_idx <- which(is.na(expr), arr.ind = TRUE)
        expr[na_idx] <- row_med[na_idx[, 1]]
      }

      q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
      needs_log <- if (identical(input$log2, "force")) TRUE
                   else if (identical(input$log2, "skip")) FALSE
                   else isTRUE(!is.na(q99) && q99 > 100)
      if (needs_log) {
        expr[expr <= 0] <- NA
        expr <- log2(expr)
        expr <- expr[stats::complete.cases(expr), , drop = FALSE]
      }

      list(
        label = label, expr = as.matrix(expr), meta = meta,
        n_samples_before = n_samples_before, n_samples_after = ncol(expr),
        n_genes_before = n_genes_before, n_genes_after = nrow(expr),
        log2_applied = needs_log
      )
    })

    output$status_ui <- renderUI({
      if (is.null(input$run) || input$run == 0) {
        return(p(class = "empty-note", icon("circle-info"),
                  "Not preprocessed yet. Set the filters above, then click \"Preprocess this dataset\"."))
      }
      res <- tryCatch(result(), error = function(e) e)
      if (inherits(res, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not preprocess this dataset:", conditionMessage(res))))
      }
      div(class = "empty-note", icon("check"),
          sprintf("%s: %s of %s samples kept, %s of %s features kept%s.",
                  res$label, res$n_samples_after, res$n_samples_before,
                  format(res$n_genes_after, big.mark = ","), format(res$n_genes_before, big.mark = ","),
                  if (res$log2_applied) ", log2-transformed" else ""))
    })

    result
  })
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

## Step title for the workflow-stepper nav (number + icon + label + a short
## sub-label), matching the visual reference. `value` is set explicitly on
## each tabPanel() below so this richer title markup doesn't change the tab's
## selectable value - anything elsewhere that does
## updateTabsetPanel(session, ns("tabs"), selected = "Merge datasets") still
## works exactly as before.
pp_step_title <- function(number, ic, label, sublabel) {
  tagList(
    span(class = "step-number", as.character(number)),
    icon(ic),
    span(class = "step-text",
         span(class = "step-label", label),
         span(class = "step-sublabel", sublabel))
  )
}

mod_preprocessing_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      8,
      div(
        class = "workflow-stepper-wrap",
        tabsetPanel(
          id = ns("tabs"), type = "tabs",
          header = tagList(
            tags$hr()
          ),
          tabPanel(
            value = "Preprocessing", title = pp_step_title(1, "broom", "Preprocessing", "Clean & filter data"),
            br(), uiOutput(ns("preprocessing_tab_ui"))
          ),
          tabPanel(
            value = "Merge datasets", title = pp_step_title(2, "code-merge", "Merge Datasets", "Combine datasets"),
            br(), uiOutput(ns("merge_tab_ui"))
          ),
          tabPanel(
            value = "Batch correction", title = pp_step_title(3, "wand-magic-sparkles", "Batch Correction", "Correct batch effects"),
            br(), uiOutput(ns("batch_tab_ui"))
          )
        )
      )
    ),
    column(
      4,
      div(
        class = "pipeline-rail",
        arthochat_shortcut_ui(
          "Not sure whether to quantile-normalise, which ComBat prior to use, or how to read a PCA batch-effect plot?",
          compact = TRUE
        )
      )
    )
  )
}

mod_preprocessing_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Shared progress signal for the right-column "Pipeline summary"
    ## timeline. `pp_sources`, `merged` and `result` are all defined further
    ## down this function; referencing them here is fine since R resolves
    ## the closure lazily, only when this reactive actually fires.
    ## Nothing here is ever marked ready/done automatically - every step
    ## needs its own explicit button click: "Preprocess this dataset" per
    ## source for step 1, merge_use_example_btn/merge_btn for step 2,
    ## run_btn for step 3.
    pp_progress <- reactive({
      merged_ok <- !is.null(tryCatch(merged(), error = function(e) NULL))
      batch_ok  <- !is.null(tryCatch(result(), error = function(e) NULL))
      pl_res <- preloaded_results()
      n_pl <- length(input$preloaded_selected %||% character(0))
      n <- max(1, n_pl)
      n_ready <- sum(vapply(pl_res, function(r) isTRUE(r$ok), logical(1)))
      list(n = n, n_ready = n_ready, merged_ok = merged_ok, batch_ok = batch_ok)
    })

    ## Right-column "Pipeline summary" timeline (see R/ui_shell.R), scoped to
    ## exactly the 3 steps this page covers - reuses pp_progress() above.
    output$pipeline_summary <- renderUI({
      pr <- pp_progress()
      step_state <- function(done, current) if (done) "done" else if (current) "current" else "future"

      steps <- list(
        list(number = 1, label = "Preprocessing", sublabel = "Filter low-quality samples & genes",
             state = step_state(pr$n_ready == pr$n, pr$n_ready < pr$n)),
        list(number = 2, label = "Merge Datasets", sublabel = "Harmonize and merge cohorts",
             state = step_state(pr$merged_ok, pr$n_ready == pr$n && !pr$merged_ok)),
        list(number = 3, label = "Batch Correction", sublabel = "Remove batch effects (ComBat)",
             state = step_state(pr$batch_ok, pr$merged_ok && !pr$batch_ok))
      )
      pipeline_summary_ui(steps)
    })

    ## The two GEO series the example cohort is built from, raw - matching
    ## the "Datasets" tab of Overview and Datasets. Read once via the same
    ## cached load_individual_dataset() every other view of these datasets
    ## already uses.
    PP_TRAINING_GEO_IDS <- c("GSE93272", "GSE110169")
    ## Joined display name for the two training cohorts above - "Whole Blood
    ## Training Cohort A and Whole Blood Training Cohort B" would repeat
    ## itself, so this is its own short-hand rather than built from
    ## pp_cohort_label() twice and pasted together.
    PP_TRAINING_COHORT_LABEL <- "Whole Blood Training Cohorts A and B"

    ## Every non-missing diagnosis group actually present across the two
    ## training GEO series, meta-only (cheap - no expression matrix read).
    ## GSE110169 carries SLE samples alongside RA/HC (GSE93272 has only
    ## RA/HC); this project's own training cohort is RA-vs-control only
    ## (Chapter_2_subchapter2_sexstratified.md: "183 samples: 145 female
    ## (86 RA, 59 healthy control) and 38 male (17 RA, 21 healthy
    ## control)" - no SLE), so example_live_merge() below defaults to
    ## HC+RA only. Every group actually present is still individually
    ## selectable for anyone who wants a different comparison.
    available_example_groups <- reactive({
      grps <- unlist(lapply(PP_TRAINING_GEO_IDS, function(gse) {
        eset <- get_raw_eset(gse)
        if (is.null(eset)) return(character(0))
        unique(stats::na.omit(eset_harmonize_meta(eset, gse)$group))
      }))
      sort(unique(grps))
    })

    ## =====================================================================
    ## Preprocessing tab: N per-dataset upload/filter blocks
    ## =====================================================================

    pp_source_default_gse <- function(i) if (i <= length(PP_TRAINING_GEO_IDS)) PP_TRAINING_GEO_IDS[i] else NULL

    pp_sources <- lapply(seq_len(MAX_PP_SOURCES), function(i) {
      mod_pp_source_server(paste0("src", i), default_label = paste("Dataset", i),
                            default_gse = pp_source_default_gse(i), dataset = dataset)
    })

    ## observeEvent()+reactiveVal() for the "Use a preloaded dataset" box:
    ## each selected choice (a GSE ID, or "__current__" for whatever's loaded
    ## via the Dataset tab) is read and preprocessed independently via
    ## pp_preloaded_read(), so one bad selection doesn't block the others.
    preloaded_results_val <- reactiveVal(NULL)
    observeEvent(input$preloaded_run, {
      req(length(input$preloaded_selected) > 0)
      choices <- input$preloaded_selected
      out <- lapply(choices, function(choice_id) {
        tryCatch(
          list(ok = TRUE, value = pp_preloaded_read(choice_id, input$preloaded_log2, dataset)),
          error = function(e) list(ok = FALSE,
                                    label = if (identical(choice_id, "__current__")) "Currently Loaded Dataset" else pp_cohort_label(choice_id),
                                    error = conditionMessage(e))
        )
      })
      names(out) <- vapply(seq_along(out), function(i) {
        if (isTRUE(out[[i]]$ok)) out[[i]]$value$label else out[[i]]$label
      }, character(1))
      preloaded_results_val(out)
    }, ignoreInit = TRUE)
    preloaded_results <- function() preloaded_results_val()

    output$preloaded_status_ui <- renderUI({
      if (is.null(input$preloaded_run) || input$preloaded_run == 0) {
        return(p(class = "empty-note", icon("circle-info"),
                  "Select one or more cohorts above, then click \"Load and Preprocess Selected Cohorts\"."))
      }
      res <- preloaded_results()
      if (is.null(res)) {
        return(div(class = "empty-note", icon("triangle-exclamation"), "Select at least one cohort above first."))
      }
      tagList(lapply(res, function(r) {
        if (isTRUE(r$ok)) {
          v <- r$value
          div(class = "empty-note", icon("check"),
              sprintf("%s: %s of %s samples kept, %s of %s features kept%s.",
                      v$label, v$n_samples_after, v$n_samples_before,
                      format(v$n_genes_after, big.mark = ","), format(v$n_genes_before, big.mark = ","),
                      if (v$log2_applied) ", log2-transformed" else ""))
        } else {
          div(class = "empty-note", icon("triangle-exclamation"),
              sprintf("%s: %s", r$label, r$error))
        }
      }))
    })

    ## Feeds merge_inputs() below: a preloaded box (any of the app's bundled
    ## raw GEO sources, or whatever's currently loaded via the Dataset tab -
    ## including your own upload, made there) - no upload here. Whatever's
    ## loaded from it becomes selectable in the Merge tab.
    output$preprocessing_tab_ui <- renderUI({
      tagList(
        box(
          width = 12, title = tagList(icon("database"), " Preloaded Data"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc",
            "Select one or more bundled blood or synovium cohorts, or use whatever is currently loaded on the Dataset tab. No upload required - each selection is preprocessed independently."),
          checkboxGroupInput(ns("preloaded_selected"), "Cohorts",
                              choices = c(pp_cohort_choices(), "Currently Loaded Dataset" = "__current__"),
                              inline = TRUE),
          p(class = "empty-note", icon("triangle-exclamation"),
            "These cohorts are raw and platform-specific. Normalization and merging are handled in the next two steps."),
          radioButtons(ns("preloaded_log2"), "Log2 Transform (applied to every selected cohort)",
                       choiceNames = list("Auto-detect per cohort (recommended)", "Force log2", "Skip (raw RNA-seq counts)"),
                       choiceValues = list("auto", "force", "skip"), selected = "auto", inline = TRUE),
          actionButton(ns("preloaded_run"), "Load and Preprocess Selected Cohorts", icon = icon("play"), class = "btn-primary btn-sm"),
          div(style = "margin-top:8px;", uiOutput(ns("preloaded_status_ui")))
        )
      )
    })

    ## =====================================================================
    ## Merge datasets tab
    ## =====================================================================

    ## The Preprocessing tab's preloaded selections, in the list(label=,
    ## expr=, meta=, ...) shape produced by pp_preloaded_read(); everything
    ## successfully loaded becomes selectable in the Merge tab below.
    merge_inputs <- reactive({
      res <- preloaded_results()
      validate(need(length(res) > 0, "Load at least one dataset in the Preprocessing tab first."))
      failed <- Filter(Negate(function(r) isTRUE(r$ok)), res)
      validate(need(length(failed) == 0,
                    sprintf("%s failed to load: %s. Fix and re-run before merging.",
                            paste(vapply(failed, `[[`, character(1), "label"), collapse = ", "),
                            paste(vapply(failed, `[[`, character(1), "error"), collapse = "; "))))
      lapply(res, `[[`, "value")
    })

    output$merge_tab_ui <- renderUI({
      tagList(
        div(
          class = "card",
          div(class = "card-title", icon("code-merge"), "Merge data for this step"),
          radioButtons(
            ns("merge_mode"), NULL,
            choiceNames = list(
              tagList(
                div(class = "data-source-option-title", "Merge the example pipeline's training datasets"),
                div(class = "data-source-option-desc", "Merges the two raw training datasets on shared genes, before batch correction.")
              ),
              tagList(
                div(class = "data-source-option-title", "Merge your own data"),
                div(class = "data-source-option-desc", "Merges the datasets from the Preprocessing tab.")
              )
            ),
            choiceValues = list("example", "own"), selected = "example"
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'example'", ns("merge_mode")),
          uiOutput(ns("merge_example_ui"))
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'own'", ns("merge_mode")),
          box(
            width = 12, title = tagList(icon("dna"), " Probe-to-gene collapsing (optional)"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc",
              "Turn this on if the datasets you loaded are still at probe level (e.g. raw Affymetrix probe IDs) rather than one row per gene already - useful for merging same-platform sources the way a published method describes, without leaving the app. Applies to every dataset selected below, using the same annotation file for all of them - for different platforms, collapse each dataset separately before uploading instead."),
            checkboxInput(ns("collapse_probes"), "My selected data is at probe level - collapse to one row per gene before merging", value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s']", ns("collapse_probes")),
              div(class = "field-label-with-hint", span("Probe → gene annotation file"),
                  mod_pp_field_hint("CSV or RDS with at least two columns: the probe/feature ID (matching the expression matrix's row IDs) and the gene symbol. Column names are auto-guessed.")),
              fileInput(ns("collapse_annot_file"), NULL, accept = c(".csv", ".rds", ".Rds")),
              radioButtons(
                ns("collapse_method"), "Collapsing method",
                choiceNames = list(
                  "Median across probes per gene (recommended - matches most published methods, e.g. Zhu et al. 2021)",
                  "Highest-mean probe represents the gene (this app's own preloaded-dataset convention)",
                  "Mean across probes per gene"
                ),
                choiceValues = list("median", "maxmean", "mean"), selected = "median"
              ),
              p(class = "empty-note", icon("circle-info"),
                "Probes with no match in the annotation file, or matched to more than one gene, are dropped - never split or duplicated.")
            )
          ),
          uiOutput(ns("merge_select_ui")),
          uiOutput(ns("merge_venn_ui")),
          box(width = 12, title = tagList(icon("code-merge"), " Merge"), status = "primary", solidHeader = FALSE,
              actionButton(ns("merge_btn"), "Merge datasets", icon = icon("code-merge"), class = "btn-primary btn-sm"),
              withSpinner(uiOutput(ns("merge_summary_ui")), color = "#2563EB", type = 6))
        )
      )
    })

    ## Live merge of the two raw training datasets (example_live_merge()
    ## above) - deliberately NOT the app's already-batch-corrected `dataset`
    ## reactiveValues, so Batch Correction downstream has a genuine,
    ## uncorrected batch effect to remove instead of a near-identical
    ## before/after. Still shows the same feature-overlap Venn/region-
    ## breakdown/download the manual "own data" merge path shows, computed
    ## from the same two training GEO sources, so "how many probes in one
    ## dataset and the overlap" is visible here too.
    output$merge_example_ui <- renderUI({
      groups <- available_example_groups()
      default_groups <- if (length(intersect(c("HC", "RA"), groups)) > 0) intersect(c("HC", "RA"), groups) else groups
      ## Named generically (whatever's excluded by default), not hardcoded
      ## to "SLE" - GSE110169 happens to be the source of that extra group
      ## today, but the exclusion itself is by group identity (anything
      ## that isn't HC/RA), not by that specific label, so this stays
      ## accurate if the training sources or their diagnosis categories
      ## ever change.
      excluded_groups <- setdiff(groups, default_groups)
      excluded_note <- if (length(excluded_groups) > 0) {
        sprintf(" %s also present (%s) - tick above to include for a different comparison.",
                if (length(excluded_groups) == 1) "One other group is" else "Other groups are",
                paste(excluded_groups, collapse = ", "))
      } else ""
      tagList(
        div(
          class = "card",
          div(class = "empty-note", icon("circle-info"),
              sprintf(
                "Merges %s on shared genes, before batch correction.",
                PP_TRAINING_COHORT_LABEL
              )),
          if (length(groups) > 0) tagList(
            div(style = "margin-top:12px;",
                strong("Diagnosis groups to include"),
                p(class = "submodule-desc",
                  sprintf("Defaults to control (HC) + RA only, matching this project's own training cohort.%s", excluded_note)),
                checkboxGroupInput(ns("example_groups"), NULL, choices = groups, selected = default_groups, inline = TRUE))
          ),
          div(style = "margin-top:12px;",
              actionButton(ns("merge_use_example_btn"), "Merge these datasets", icon = icon("code-merge"), class = "btn-primary btn-sm"))
        ),
        uiOutput(ns("merge_venn_example_ui"))
      )
    })

    output$merge_venn_example_ui <- renderUI({
      if (!isTRUE((input$merge_use_example_btn %||% 0) > 0)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Click “Merge these datasets” above to merge and see the feature overlap between the two training datasets."))
      }
      tagList(
        box(width = 12, title = tagList(icon("diagram-project"), " Feature overlap across datasets"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", sprintf(
              "%s, collapsed to gene symbol.",
              PP_TRAINING_COHORT_LABEL
            )),
            fluidRow(
              column(7,
                withSpinner(plotOutput(ns("venn_plot_example"), height = 440), color = "#2563EB", type = 6),
                div(class = "table-toolbar", downloadButton(ns("download_venn_example_png"), "Download diagram (PNG)", class = "btn-primary btn-sm"))
              ),
              column(5, DT::dataTableOutput(ns("venn_table_example")))
            )),
        box(width = 12, title = tagList(icon("table-list"), " Region breakdown"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Every exact combination of datasets a feature belongs to, largest first. These are the same regions the diagram above draws, shown here as a full table."),
            div(class = "table-toolbar", downloadButton(ns("download_venn_example"), "Download region counts (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("venn_region_table_example")))
      )
    })

    ## Gene-symbol-collapsed feature sets for the two training GEO sources -
    ## the same collapse_probes_to_genes() step mod_pp_source_server's own
    ## "preloaded GEO dataset" option already applies (raw probe IDs almost
    ## never overlap across the two platforms; gene symbols do).
    example_overlap_sets <- reactive({
      sets <- lapply(PP_TRAINING_GEO_IDS, function(gse) {
        collapsed <- get_collapsed_genes(gse)
        if (is.null(collapsed)) return(NULL)
        rownames(collapsed)
      })
      names(sets) <- PP_TRAINING_GEO_IDS
      sets <- Filter(Negate(is.null), sets)
      validate(need(length(sets) >= 2, "Raw source data for the example cohort's training datasets was not found on disk."))
      sets
    })

    ## Live merge of the two RAW training datasets - NOT the app's default
    ## `dataset` reactiveValues, which is already batch-corrected
    ## (combined_expr_batchcorrected.rds). Adopting that directly here would
    ## feed already-corrected data back into the Batch Correction step,
    ## where a second ComBat pass has almost no residual batch signal left
    ## to remove - "before" and "after" PCA would look nearly identical.
    ## Rebuilding from raw (same collapse_probes_to_genes()/
    ## eset_harmonize_meta() steps mod_pp_source_server's own "preloaded GEO
    ## dataset" option already applies per-source, and the same common-gene
    ## merge mod_pp_source_server's own "own data" path below already does)
    ## gives Batch Correction a genuine, uncorrected batch effect to remove,
    ## so the before/after contrast is real - matching how the example
    ## cohort was actually built in the first place.
    example_live_merge <- reactive({
      parts <- lapply(PP_TRAINING_GEO_IDS, function(gse) {
        eset <- get_raw_eset(gse)
        validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
        expr <- get_collapsed_genes(gse)
        meta <- eset_harmonize_meta(eset, gse)
        ## input$example_groups (checkboxGroupInput above, defaults to
        ## HC+RA) - not just !is.na(meta$group) - so GSE110169's SLE
        ## samples are excluded by default instead of silently merging in
        ## alongside RA/HC, matching this project's own RA-vs-control
        ## training cohort. Falls back to every non-NA group if the input
        ## hasn't rendered yet (e.g. very first reactive pass).
        wanted_groups <- input$example_groups %||% available_example_groups()
        keep <- !is.na(meta$group) & meta$group %in% wanted_groups
        meta <- meta[keep, , drop = FALSE]
        expr <- expr[, meta$sample, drop = FALSE]

        ## Same auto-detect log2 rule every per-source preprocessing path in
        ## this app uses (see mod_pp_source_server's result()) - a no-op for
        ## these two datasets (already log2, per the Scale check panel
        ## above), but applied for correctness on whatever's actually on disk.
        q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
        if (isTRUE(!is.na(q99) && q99 > 100)) {
          expr[expr <= 0] <- NA
          expr <- log2(expr)
          expr <- expr[stats::complete.cases(expr), , drop = FALSE]
        }
        if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_
        list(expr = expr, meta = meta, label = pp_cohort_label(gse))
      })

      common <- Reduce(intersect, lapply(parts, function(x) rownames(x$expr)))
      validate(need(length(common) >= 20, "Fewer than 20 common genes between the two training datasets."))
      merged_expr <- do.call(cbind, lapply(parts, function(x) x$expr[common, , drop = FALSE]))
      metas <- lapply(parts, function(x) { m <- x$meta; m$dataset <- x$label; m })
      all_cols <- unique(unlist(lapply(metas, colnames)))
      metas <- lapply(metas, function(m) {
        missing_cols <- setdiff(all_cols, colnames(m))
        for (cl in missing_cols) m[[cl]] <- NA
        m[, all_cols, drop = FALSE]
      })
      merged_meta <- do.call(rbind, metas)
      if (!"batch" %in% colnames(merged_meta) || all(is.na(merged_meta$batch))) merged_meta$batch <- merged_meta$dataset
      stopifnot(identical(colnames(merged_expr), merged_meta$sample))
      list(expr = merged_expr, meta = merged_meta, sources = PP_TRAINING_COHORT_LABEL)
    })

    venn_plot_example_obj <- reactive({
      draw_overlap_venn(example_overlap_sets())
    })

    output$venn_plot_example <- renderPlot({
      venn_plot_example_obj()
    })

    output$download_venn_example_png <- downloadHandler(
      filename = function() "example_feature_overlap_venn.png",
      content = function(file) ggsave(file, plot = venn_plot_example_obj(), width = 7.5, height = 6.5, dpi = 300, bg = "white")
    )

    output$venn_table_example <- DT::renderDataTable({
      sets <- example_overlap_sets()
      common <- Reduce(intersect, sets)
      df <- data.frame(dataset = c(names(sets), "Common to all"), n_features = c(lengths(sets), length(common)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    venn_regions_example <- reactive(overlap_region_sizes(example_overlap_sets()))

    output$venn_region_table_example <- DT::renderDataTable({
      DT::datatable(venn_regions_example(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_venn_example <- downloadHandler(
      filename = function() "example_feature_overlap_regions.csv",
      content = function(file) write.csv(venn_regions_example(), file, row.names = FALSE)
    )

    ## Which configured, preprocessed datasets actually go into this merge is
    ## the user's call, not automatic - unchecking one here drops it from
    ## both the overlap diagram below and the merge itself, so what you see
    ## is exactly what you'll get.
    output$merge_select_ui <- renderUI({
      ## Silent (not validate()'d) when sources aren't ready yet - merge_venn_ui
      ## below is what surfaces that "preprocess dataset X first" message, so
      ## it isn't shown twice.
      lst <- tryCatch(merge_inputs(), error = function(e) NULL)
      req(length(lst) >= 2)
      labels <- vapply(lst, `[[`, character(1), "label")
      box(width = 12, title = tagList(icon("check-double"), " Choose which datasets to merge"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "All preprocessed datasets are selected by default. Uncheck any dataset you do not want included in this merge."),
          checkboxGroupInput(ns("merge_selected"), NULL, choices = labels, selected = labels, inline = TRUE))
    })

    ## Shared annotation table for the optional probe-collapsing step below -
    ## read once, reused for every selected dataset (same-platform merges
    ## only; see the Merge tab's own note on this).
    collapse_annot <- reactive({
      validate(need(!is.null(input$collapse_annot_file),
                    "Probe-to-gene collapsing is turned on, but no annotation file has been uploaded yet - add one above, or turn the checkbox off if your data is already at gene level."))
      path <- input$collapse_annot_file$datapath
      if (grepl("\\.rds$", input$collapse_annot_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
        validate(need(is.data.frame(d), "The annotation RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    selected_lst <- reactive({
      lst <- merge_inputs()
      labels <- vapply(lst, `[[`, character(1), "label")
      sel <- if (length(lst) < 2) labels else (input$merge_selected %||% labels)
      validate(need(length(sel) >= 1, "Select at least one dataset to merge."))
      lst <- lst[labels %in% sel]
      if (isTRUE(input$collapse_probes)) {
        annot <- collapse_annot()
        lst <- lapply(lst, function(x) {
          x$expr <- pp_collapse_probes_to_genes(x$expr, annot, input$collapse_method %||% "median")
          x
        })
      }
      lst
    })

    output$merge_venn_ui <- renderUI({
      ## tryCatch, not a bare selected_lst() call: merge_select_ui above
      ## already surfaces "load a dataset first" nicely (via its own
      ## tryCatch + req()), but this output has no such guard, so an
      ## unready merge_inputs()/selected_lst() (e.g. no dataset loaded, or
      ## no probe annotation uploaded yet) was bubbling up as Shiny's bare,
      ## unstyled default validation message instead of this app's normal
      ## empty-note styling - easy to miss/read as "nothing rendered".
      lst <- tryCatch(selected_lst(), error = function(e) e)
      if (inherits(lst, "error")) {
        return(div(class = "empty-note", icon("circle-info"), conditionMessage(lst)))
      }
      if (length(lst) < 2) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Fewer than two datasets are selected, so there is nothing to compare. Click \"Merge datasets\" below to continue with just this one, or select more datasets above."))
      }
      tagList(
        box(width = 12, title = tagList(icon("diagram-project"), " Feature overlap across datasets"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Only features (genes/probes) present in every selected dataset are kept in the merge."),
            fluidRow(
              column(7,
                withSpinner(plotOutput(ns("venn_plot_custom"), height = 440), color = "#2563EB", type = 6),
                div(class = "table-toolbar", downloadButton(ns("download_venn_custom_png"), "Download diagram (PNG)", class = "btn-primary btn-sm"))
              ),
              column(5, DT::dataTableOutput(ns("venn_table_custom")))
            )),
        box(width = 12, title = tagList(icon("table-list"), " Region breakdown"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Every exact combination of datasets a feature belongs to, largest first. These are the same regions the diagram above draws, shown here as a full table."),
            div(class = "table-toolbar", downloadButton(ns("download_venn_custom"), "Download region counts (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("venn_region_table_custom")))
      )
    })

    overlap_sets <- reactive({
      lst <- selected_lst()
      validate(need(length(lst) >= 2, "Fewer than two datasets are selected, so there is nothing to compare."))
      setNames(lapply(lst, function(x) rownames(x$expr)), vapply(lst, `[[`, character(1), "label"))
    })

    venn_plot_custom_obj <- reactive({
      draw_overlap_venn(overlap_sets())
    })

    output$venn_plot_custom <- renderPlot({
      venn_plot_custom_obj()
    })

    output$download_venn_custom_png <- downloadHandler(
      filename = function() "feature_overlap_venn.png",
      content = function(file) ggsave(file, plot = venn_plot_custom_obj(), width = 7.5, height = 6.5, dpi = 300, bg = "white")
    )

    output$venn_table_custom <- DT::renderDataTable({
      sets <- overlap_sets()
      common <- Reduce(intersect, sets)
      df <- data.frame(dataset = c(names(sets), "Common to all"), n_features = c(lengths(sets), length(common)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    venn_regions_custom <- reactive(overlap_region_sizes(overlap_sets()))

    output$venn_region_table_custom <- DT::renderDataTable({
      DT::datatable(venn_regions_custom(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_venn_custom <- downloadHandler(
      filename = function() "feature_overlap_regions.csv",
      content = function(file) write.csv(venn_regions_custom(), file, row.names = FALSE)
    )

    ## Triggered by either button: "Use this dataset" (example path, added
    ## below) or "Merge datasets" (the original, unchanged own-data path).
    ## ignoreInit = TRUE keeps this from firing at session start, the same
    ## guarantee a single eventReactive(input$merge_btn, ...) already had.
    merged <- eventReactive(list(input$merge_btn, input$merge_use_example_btn), {
      if (identical(input$merge_mode, "example")) {
        return(example_live_merge())
      }

      lst <- selected_lst()
      if (length(lst) == 1) {
        x <- lst[[1]]
        meta <- x$meta
        if (!"dataset" %in% colnames(meta)) meta$dataset <- x$label
        return(list(expr = x$expr, meta = meta, sources = x$label))
      }
      sets <- lapply(lst, function(x) rownames(x$expr))
      common <- Reduce(intersect, sets)
      validate(need(length(common) >= 20,
                    "Fewer than 20 features are in common across the selected datasets. Check that every uploaded dataset uses the same type of row name, for example all gene symbols or all the same probe IDs."))
      merged_expr <- do.call(cbind, lapply(lst, function(x) x$expr[common, , drop = FALSE]))
      metas <- lapply(lst, function(x) { m <- x$meta; m$dataset <- x$label; m })
      all_cols <- unique(unlist(lapply(metas, colnames)))
      metas <- lapply(metas, function(m) {
        missing_cols <- setdiff(all_cols, colnames(m))
        for (cl in missing_cols) m[[cl]] <- NA
        m[, all_cols, drop = FALSE]
      })
      merged_meta <- tryCatch(do.call(rbind, metas), error = function(e) {
        validate("Could not combine metadata across datasets. A column with the same name has a different type in different datasets, for example numeric in one and text in another. Rename or fix that column, then preprocess again.")
      })
      if (!"batch" %in% colnames(merged_meta) || all(is.na(merged_meta$batch))) merged_meta$batch <- merged_meta$dataset
      stopifnot(identical(colnames(merged_expr), merged_meta$sample))
      list(expr = merged_expr, meta = merged_meta, sources = paste(vapply(lst, `[[`, character(1), "label"), collapse = " + "))
    }, ignoreInit = TRUE)

    output$merge_summary_ui <- renderUI({
      if (!isTRUE((input$merge_btn %||% 0) > 0) && !isTRUE((input$merge_use_example_btn %||% 0) > 0)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Select datasets above, then click \"Merge datasets\"."))
      }
      m <- tryCatch(merged(), error = function(e) e)
      if (inherits(m, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), conditionMessage(m)))
      }
      tagList(
        div(class = "empty-note", icon("check"),
            sprintf("Merged: %s samples x %s common features, from %s.", ncol(m$expr), format(nrow(m$expr), big.mark = ","), m$sources)),
        div(class = "table-toolbar",
            downloadButton(ns("download_merged_expr"), "Download expression matrix (CSV)", class = "btn-sm"),
            downloadButton(ns("download_merged_meta"), "Download sample metadata (CSV)", class = "btn-sm"),
            downloadButton(ns("download_merged_rds"), "Download both (RDS)", class = "btn-sm")),
        p(class = "submodule-desc", "Composition by dataset (and group, if mapped):"),
        DT::dataTableOutput(ns("merge_composition_table"))
      )
    })

    output$download_merged_expr <- downloadHandler(
      filename = function() "merged_expression.csv",
      content = function(file) {
        m <- merged()
        write.csv(data.frame(feature = rownames(m$expr), m$expr, check.names = FALSE), file, row.names = FALSE)
      }
    )
    output$download_merged_meta <- downloadHandler(
      filename = function() "merged_metadata.csv",
      content = function(file) write.csv(merged()$meta, file, row.names = FALSE)
    )
    output$download_merged_rds <- downloadHandler(
      filename = function() "merged_dataset.rds",
      content = function(file) {
        m <- merged()
        saveRDS(list(expr = m$expr, meta = m$meta, sources = m$sources), file)
      }
    )

    output$merge_composition_table <- DT::renderDataTable({
      m <- merged()
      grp <- if ("group" %in% colnames(m$meta)) m$meta$group else rep("(no group column)", nrow(m$meta))
      df <- as.data.frame(table(dataset = m$meta$dataset, group = grp))
      colnames(df)[3] <- "n_samples"
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    ## =====================================================================
    ## Batch correction tab: normalisation + ComBat, live, on the merged (or
    ## single, un-merged) dataset from the Merge datasets tab.
    ## =====================================================================

    active_meta_df <- reactive({
      m <- merged()
      m$meta
    })

    output$settings_ui <- renderUI({
      m <- tryCatch(merged(), error = function(e) NULL)
      if (is.null(m)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Finish the Merge datasets tab first (or preprocess and merge just one dataset there) before configuring batch correction."))
      }
      cols <- colnames(m$meta)
      batch_default <- intersect(c("batch", "batch_full", "dataset"), cols)
      protect_default <- intersect(c("group", "sex"), cols)
      tagList(
        selectInput(ns("color_by"), "Color PCA by", choices = cols,
                    selected = if ("group" %in% cols) "group" else cols[1], selectize = FALSE),
        fluidRow(
          column(6, selectInput(ns("pc_x"), "X axis", choices = setNames(1:5, paste0("PC", 1:5)), selected = 1, selectize = FALSE)),
          column(6, selectInput(ns("pc_y"), "Y axis", choices = setNames(1:5, paste0("PC", 1:5)), selected = 2, selectize = FALSE))
        ),
        checkboxInput(ns("show_ellipse"), "Show group confidence ellipses", value = TRUE),
        checkboxInput(ns("show_labels"), "Label points with sample ID", value = FALSE),
        selectInput(ns("batch_col"), "Batch column to correct for",
                    choices = cols,
                    selected = if (length(batch_default) > 0) batch_default[1] else cols[1], selectize = FALSE),
        {
          tagList(
            radioButtons(
              ns("norm_method"), "Normalisation",
              choiceNames = list(
                "Auto-detect (recommended)",
                "Already normalised, skip",
                "Quantile normalisation (microarray or log-scale data)",
                "TMM plus log2-CPM (RNA-seq raw counts)"
              ),
              choiceValues = list("auto", "skip", "quantile", "tmm"),
              selected = "auto"
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] != 'tmm'", ns("norm_method")),
              sliderInput(ns("min_pct"), "Drop genes below this percentile of mean expression", min = 0, max = 90, value = 0, step = 5)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'tmm'", ns("norm_method")),
              p(class = "empty-note", icon("circle-info"),
                "Low-count genes are filtered with edgeR::filterByExpr() by group instead of the percentile slider. Set log2 to \"Skip\" for this dataset in the Preprocessing tab first."),
              radioButtons(
                ns("tmm_correction_stage"), "Batch-correct raw counts, or the normalised log2-CPM?",
                choiceNames = list(
                  "After TMM: ComBat or limma on log2-CPM below (typical batch effects)",
                  "Before TMM: ComBat-seq on raw counts (Zhang, Parmigiani and Johnson 2020, better for large, count-driven batch effects)"
                ),
                choiceValues = list("post", "pre"), selected = "post"
              ),
              p(class = "empty-note", icon("circle-info"),
                "ComBat-seq (before TMM) ignores the correction method, prior, reference batch and exclude-outliers options below. It always protects the group column directly, and uses the batch column and biological covariates chosen here.")
            ),
            selectInput(ns("protect_cols"), "Biological covariates to protect (won't be treated as batch)",
                        choices = cols, selected = protect_default, multiple = TRUE),
            checkboxInput(ns("skip_combat"), "Skip batch correction (normalise only)", value = FALSE),
            tags$hr(),
            checkboxInput(ns("show_advanced"), strong("Show advanced batch-correction options"), value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s']", ns("show_advanced")),
              selectInput(ns("correction_method"), "Correction method",
                          choices = c("ComBat (empirical Bayes, recommended)" = "combat",
                                      "limma::removeBatchEffect (simple linear adjustment)" = "limma",
                                      "Surrogate Variable Analysis - SVA (unknown/hidden sources, no batch column needed)" = "sva"),
                          selected = "combat", selectize = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'combat'", ns("correction_method")),
                radioButtons(ns("combat_prior"), "ComBat empirical Bayes prior",
                             choiceNames = list("Parametric (faster, default)", "Non-parametric (slower, more robust for small or uneven batches)"),
                             choiceValues = list("param", "nonparam"), selected = "param"),
                checkboxInput(ns("combat_mean_only"), "Adjust batch mean only (skip variance adjustment)", value = FALSE),
                uiOutput(ns("ref_batch_ui"))
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'sva'", ns("correction_method")),
                p(class = "empty-note", icon("circle-info"),
                  "SVA estimates unwanted variation directly from the data instead of using the batch column above - useful when the real source of batch effects is unknown or only partly captured by a column you have. It still protects the biological covariates chosen below. Leek JT et al., Bioinformatics 2012;28(6):882-883."),
                numericInput(ns("sva_n_sv"), "Number of surrogate variables (0 = auto-estimate)", value = 0, min = 0, max = 20, step = 1)
              ),
              selectInput(ns("batch_col2"), "Combine with a second column into an interaction batch (optional, for example dataset by scan batch)",
                          choices = c("(none)", cols), selectize = FALSE),
              sliderInput(ns("variance_pct"), "Also drop genes below this percentile of variance", min = 0, max = 90, value = 0, step = 5),
              checkboxInput(ns("exclude_outliers"), "Exclude samples flagged as QC outliers before correcting (uses the outlier sensitivity below)", value = FALSE)
            )
          )
        },
        sliderInput(ns("mad_k"), "Outlier sensitivity (MADs from the cohort median)", min = 2, max = 6, value = 3, step = 0.5),
        actionButton(ns("run_btn"), "Run normalisation and batch correction",
                      icon = icon("play"), class = "btn-primary btn-sm")
      )
    })

    output$ref_batch_ui <- renderUI({
      req(input$batch_col)
      meta <- tryCatch(active_meta_df(), error = function(e) NULL)
      req(meta)
      lvls <- sort(unique(stats::na.omit(as.character(meta[[input$batch_col]]))))
      selectInput(ns("ref_batch"), "Reference batch (optional). Other batches are shifted to match this one instead of a pooled average.",
                  choices = c("(none)", lvls), selected = "(none)", selectize = FALSE)
    })

    ## A column can't be both the thing corrected for and a thing protected
    ## from correction - keep it out of reach in the UI instead of letting
    ## ComBat/limma resolve the resulting confound unpredictably (see the
    ## defensive setdiff() a bit further down, which still applies if this
    ## observer hasn't fired yet).
    observeEvent(list(input$batch_col, input$batch_col2), {
      meta <- tryCatch(active_meta_df(), error = function(e) NULL)
      req(meta)
      cols <- colnames(meta)
      exclude <- c(input$batch_col, if (!identical(input$batch_col2 %||% "(none)", "(none)")) input$batch_col2 else NULL)
      choices <- setdiff(cols, exclude)
      current <- intersect(input$protect_cols %||% character(0), choices)
      updateSelectInput(session, "protect_cols", choices = choices, selected = current)
    }, ignoreInit = TRUE)

    result <- eventReactive(input$run_btn, {
      req(input$batch_col, input$color_by)

      {
        m <- merged()
        expr <- m$expr
        meta <- m$meta
        sources <- m$sources

        norm_method <- input$norm_method %||% "auto"
        skip_combat <- isTRUE(input$skip_combat)
        already_corrected <- FALSE

        if (identical(norm_method, "tmm")) {
          ## TMM + log2-CPM, matching scripts/goal2_sex_stratified/
          ## 20_testing_synovium_external.R: edgeR::filterByExpr() by group,
          ## calcNormFactors(method = "TMM"), then log2-CPM. For raw RNA-seq
          ## counts, not log-scale microarray/already-normalised data.
          validate(need(all(expr >= 0, na.rm = TRUE),
                        "TMM normalisation expects raw, non-negative counts, but this data has negative values, which suggests it is already log-transformed. Preprocess this dataset again with log2 set to \"Skip\"."))
          validate(need("group" %in% colnames(meta), "TMM normalisation needs a group column to filter low-count genes by."))
          grp <- factor(meta$group)
          validate(need(length(unique(na.omit(grp))) >= 2, "TMM normalisation needs at least two group levels."))
          counts <- round(as.matrix(expr))
          storage.mode(counts) <- "integer"
          dge0 <- edgeR::DGEList(counts = counts)
          keepg <- edgeR::filterByExpr(dge0, group = grp)
          n_before <- nrow(expr)
          validate(need(sum(keepg) >= 50, "Fewer than 50 genes pass edgeR's expression filter for this group split."))
          counts_f <- counts[keepg, , drop = FALSE]

          tmm_stage <- input$tmm_correction_stage %||% "post"
          if (!skip_combat && identical(tmm_stage, "pre")) {
            ## ComBat-seq: batch-correct the raw counts directly with a
            ## negative-binomial model, THEN TMM-normalise the corrected
            ## counts - better suited to batch effects that are themselves
            ## count-scale (library-size/depth-driven) rather than effects
            ## that only show up after log-CPM.
            ## Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect
            ## adjustment for RNA-seq count data. Brief Bioinform
            ## 2020;21(6):1954-1964.
            cs_batch_primary <- as.character(meta[[input$batch_col]])
            cs_use_batch2 <- !identical(input$batch_col2 %||% "(none)", "(none)") && (input$batch_col2 %in% colnames(meta))
            cs_batch <- if (cs_use_batch2) paste(cs_batch_primary, as.character(meta[[input$batch_col2]]), sep = "_") else cs_batch_primary
            validate(need(length(unique(na.omit(cs_batch))) >= 2, "The chosen batch column (or combination) needs at least two levels for ComBat-seq."))
            validate(need(all(table(cs_batch) >= 2), "Every level of the chosen batch column (or combination) needs at least 2 samples for ComBat-seq."))

            cs_batch_cols_used <- c(input$batch_col, if (cs_use_batch2) input$batch_col2 else NULL)
            cs_protect <- intersect(input$protect_cols %||% character(0), colnames(meta))
            protect_dropped_for_batch <- intersect(setdiff(cs_protect, "group"), cs_batch_cols_used)
            cs_covar_cols <- setdiff(cs_protect, c("group", cs_batch_cols_used))
            cs_covar_cols <- cs_covar_cols[vapply(cs_covar_cols, function(cl) length(unique(na.omit(meta[[cl]]))) >= 2, logical(1))]
            cs_covar_mod <- if (length(cs_covar_cols) > 0) {
              meta_mod <- meta
              for (cl in cs_covar_cols) meta_mod[[cl]] <- ifelse(is.na(meta_mod[[cl]]), "Unknown", meta_mod[[cl]])
              stats::model.matrix(stats::as.formula(paste("~", paste(cs_covar_cols, collapse = " + "))), data = meta_mod)
            } else {
              NULL
            }

            counts_adj <- tryCatch(
              sva::ComBat_seq(counts = counts_f, batch = cs_batch, group = as.character(grp), covar_mod = cs_covar_mod),
              error = function(e) sva::ComBat_seq(counts = counts_f, batch = cs_batch_primary, group = as.character(grp))
            )
            dge_before <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_f), method = "TMM")
            dge_after  <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_adj), method = "TMM")
            expr_prenorm <- counts_f
            expr_qnorm  <- edgeR::cpm(dge_before, log = TRUE, prior.count = 1)
            expr_combat <- edgeR::cpm(dge_after,  log = TRUE, prior.count = 1)
            already_corrected <- TRUE
            correction_method <- "combat_seq"
            protect <- intersect(c("group", cs_covar_cols), c("group", cs_protect))
            use_batch2 <- cs_use_batch2; ref_batch <- NULL
            combat_prior <- NA_character_; combat_mean_only <- FALSE
            norm_label <- "TMM plus log2-CPM, batch-corrected on raw counts with ComBat-seq before normalisation"
          } else {
            dge <- edgeR::calcNormFactors(dge0[keepg, , keep.lib.sizes = FALSE], method = "TMM")
            expr_prenorm <- counts_f
            expr_qnorm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
            norm_label <- "TMM (edgeR::calcNormFactors) plus log2-CPM"
          }
          needs_log <- FALSE; q99 <- NA_real_
          apply_qnorm <- TRUE
        } else {
          ## low-expression + low-variance gene filter, on the merged matrix
          n_before <- nrow(expr)
          gene_mean <- rowMeans(expr, na.rm = TRUE)
          gene_var  <- apply(expr, 1, stats::var)
          mean_cutoff <- stats::quantile(gene_mean, input$min_pct / 100)
          var_cutoff  <- stats::quantile(gene_var, (input$variance_pct %||% 0) / 100)
          keep <- gene_mean >= mean_cutoff & gene_var > 0 & gene_var >= var_cutoff
          expr <- expr[keep, , drop = FALSE]
          validate(need(nrow(expr) >= 50, "Fewer than 50 genes remain after filtering. Lower the expression/variance percentile cutoffs."))
          expr_prenorm <- expr
          needs_log <- FALSE; q99 <- NA_real_

          diag_before <- summarize_norm_diagnostics(expr_prenorm)
          apply_qnorm <- switch(norm_method,
            skip = FALSE,
            quantile = TRUE,
            needs_quantile_norm(diag_before) # "auto"
          )
          expr_qnorm <- if (apply_qnorm) {
            mtx <- limma::normalizeBetweenArrays(as.matrix(expr_prenorm), method = "quantile")
            rownames(mtx) <- rownames(expr_prenorm); colnames(mtx) <- colnames(expr_prenorm)
            mtx
          } else {
            as.matrix(expr_prenorm)
          }
          norm_label <- switch(norm_method,
            skip = "None, used as loaded",
            quantile = "Quantile normalisation (forced)",
            if (apply_qnorm) "Quantile normalisation (auto-detected as needed)" else "None, auto-detected as already normalised"
          )
        }

        ## Everything below only applies to the standard, post-normalisation
        ## correction path - ComBat-seq (already_corrected == TRUE) already
        ## produced expr_combat directly, on the raw counts, above.
        n_excluded_outliers <- 0L
        if (!already_corrected) {
          ## optional: drop samples flagged as QC outliers before correcting,
          ## instead of only flagging them afterwards
          if (isTRUE(input$exclude_outliers)) {
            qc_pre <- compute_sample_qc(expr_qnorm, mad_k = input$mad_k)
            flagged <- qc_pre$sample[qc_pre$flag_signal | qc_pre$flag_detected | qc_pre$flag_cor]
            if (length(flagged) > 0) {
              validate(need(ncol(expr_qnorm) - length(flagged) >= 6,
                            "Excluding QC-flagged samples would leave fewer than 6 samples. Lower the outlier sensitivity, or turn this option off."))
              keep_samples <- setdiff(colnames(expr_qnorm), flagged)
              expr_prenorm <- expr_prenorm[, keep_samples, drop = FALSE]
              expr_qnorm   <- expr_qnorm[, keep_samples, drop = FALSE]
              meta <- meta[match(keep_samples, meta$sample), , drop = FALSE]
              n_excluded_outliers <- length(flagged)
            }
          }

          ## batch correction, protecting whichever covariates were selected -
          ## the batch column only needs to be a real, 2+ level split when a
          ## correction is actually going to run against it.
          batch_primary <- as.character(meta[[input$batch_col]])
          use_batch2 <- !identical(input$batch_col2 %||% "(none)", "(none)") && (input$batch_col2 %in% colnames(meta))
          batch <- if (use_batch2) paste(batch_primary, as.character(meta[[input$batch_col2]]), sep = "_") else batch_primary
          if (!skip_combat) {
            validate(need(length(unique(na.omit(batch))) >= 2, "The chosen batch column (or combination) needs at least two levels. If you don't need batch correction, tick \"Skip batch correction\" above."))
            validate(need(all(table(batch) >= 2), "Every level of the chosen batch column (or combination) needs at least 2 samples for correction."))
          }

          protect <- intersect(input$protect_cols %||% character(0), colnames(meta))
          protect <- protect[vapply(protect, function(cl) length(unique(na.omit(meta[[cl]]))) >= 2, logical(1))]
          ## Protecting the exact column being corrected for is a degenerate
          ## design - batch and mod become collinear, and depending on which
          ## fallback ComBat/limma lands on, that can either silently strip
          ## the "protected" signal (mod gets absorbed into batch) or refuse
          ## to remove any batch effect at all. Drop the overlap up front
          ## rather than let ComBat resolve it unpredictably.
          batch_cols_used <- c(input$batch_col, if (use_batch2) input$batch_col2 else NULL)
          protect_dropped_for_batch <- intersect(protect, batch_cols_used)
          protect <- setdiff(protect, batch_cols_used)
          mod <- if (length(protect) > 0) {
            meta_mod <- meta
            for (cl in protect) meta_mod[[cl]] <- ifelse(is.na(meta_mod[[cl]]), "Unknown", meta_mod[[cl]])
            stats::model.matrix(stats::as.formula(paste("~", paste(protect, collapse = " + "))), data = meta_mod)
          } else {
            NULL
          }

          combat_prior <- input$combat_prior %||% "param"
          combat_mean_only <- isTRUE(input$combat_mean_only)
          correction_method <- input$correction_method %||% "combat"
          ref_batch <- if (!use_batch2 && !identical(input$ref_batch %||% "(none)", "(none)")) input$ref_batch else NULL

          run_combat <- function(b, use_mod = TRUE, use_ref = TRUE) {
            sva::ComBat(dat = expr_qnorm, batch = b, mod = if (use_mod) mod else NULL,
                        par.prior = identical(combat_prior, "param"), mean.only = combat_mean_only,
                        ref.batch = if (use_ref) ref_batch else NULL)
          }
          run_limma <- function(b) {
            design <- if (!is.null(mod)) mod else matrix(1, ncol(expr_qnorm), 1)
            limma::removeBatchEffect(expr_qnorm, batch = b, design = design)
          }
          ## SVA doesn't take a batch label at all - it estimates unknown/
          ## residual sources of unwanted variation straight from the data
          ## (protecting `mod`'s biological covariates from being absorbed
          ## into those surrogate variables), then the estimated surrogate
          ## variables are regressed out the same way limma::removeBatchEffect
          ## regresses out a known batch - this is the sva package's own
          ## documented recipe for producing SVA-adjusted data for
          ## visualisation/downstream use, not a novel technique.
          ## Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. The sva
          ## package for removing batch effects. Bioinformatics 2012;28(6):882-883.
          run_sva <- function() {
            mod_full <- if (!is.null(mod)) mod else matrix(1, ncol(expr_qnorm), 1)
            mod0 <- matrix(1, ncol(expr_qnorm), 1)
            ## num.sv's permutation-based "be" method and sva()'s own eigen
            ## decomposition are O(genes) per iteration; vfilter restricts
            ## both to the most variable genes so estimation stays fast even
            ## on tens of thousands of features, as recommended by the sva
            ## package vignette for large expression matrices.
            n_genes <- nrow(expr_qnorm)
            vfilt <- if (n_genes > 2000) 2000L else NULL
            n_sv <- as.integer(input$sva_n_sv %||% 0)
            if (n_sv <= 0) {
              n_sv <- tryCatch(sva::num.sv(as.matrix(expr_qnorm), mod_full, method = "be", vfilter = vfilt),
                                error = function(e) NA_integer_)
            }
            n_sv <- if (is.na(n_sv)) 1L else max(1L, min(n_sv, ncol(expr_qnorm) - ncol(mod_full) - 1L, 20L))
            sv_obj <- sva::sva(as.matrix(expr_qnorm), mod_full, mod0, n.sv = n_sv, vfilter = vfilt)
            validate(need(sv_obj$n.sv >= 1, "SVA did not find any significant hidden sources of variation to correct for - try ComBat or limma instead, or set the number of surrogate variables manually."))
            limma::removeBatchEffect(expr_qnorm, covariates = sv_obj$sv, design = mod_full)
          }

          expr_combat <- if (skip_combat) {
            expr_qnorm
          } else if (identical(correction_method, "limma")) {
            tryCatch(run_limma(batch), error = function(e) run_limma(batch_primary))
          } else if (identical(correction_method, "sva")) {
            run_sva()
          } else {
            ## ladder of fallbacks, same idea as scripts/00_shared/03_normalize_batch.R's
            ## run_combat(): try the requested batch/mod/ref.batch combination first,
            ## then progressively drop ref.batch, the interaction, and finally mod.
            tryCatch(run_combat(batch, use_mod = TRUE, use_ref = TRUE),
              error = function(e) tryCatch(run_combat(batch_primary, use_mod = TRUE, use_ref = FALSE),
                error = function(e2) run_combat(batch_primary, use_mod = FALSE, use_ref = FALSE)))
          }
        }
        n_after <- nrow(expr_prenorm)
      }

      diag_before <- summarize_norm_diagnostics(expr_prenorm)
      diag_after  <- summarize_norm_diagnostics(expr_qnorm)
      norm_diag <- rbind(cbind(stage = "before normalisation", diag_before),
                          cbind(stage = "after normalisation", diag_after))

      qc <- compute_sample_qc(expr_combat, mad_k = input$mad_k)
      qc <- merge(qc, meta[, intersect(c("sample", "group"), colnames(meta))], by = "sample", all.x = TRUE)

      before <- pca_of(expr_qnorm)
      after  <- pca_of(expr_combat)

      list(
        expr_prenorm = expr_prenorm, expr_qnorm = expr_qnorm, expr_combat = expr_combat,
        before = before, after = after, meta = meta, qc = qc,
        n_before = n_before, n_after = n_after,
        needs_log = needs_log, q99 = q99, apply_qnorm = apply_qnorm, norm_label = norm_label, norm_diag = norm_diag,
        protect = protect, protect_dropped_for_batch = protect_dropped_for_batch,
        batch_col = input$batch_col, color_by = input$color_by,
        skip_combat = skip_combat, combat_prior = combat_prior,
        combat_mean_only = combat_mean_only, sources = sources, correction_method = correction_method,
        n_excluded_outliers = n_excluded_outliers, use_batch2 = use_batch2, ref_batch = ref_batch
      )
    })

    ## ---- Value boxes ------------------------------------------------------

    output$vb_samples <- renderValueBox({
      res <- result()
      valueBox(nrow(res$meta), "Samples", icon = icon("users"), color = "light-blue")
    })
    output$vb_genes_kept <- renderValueBox({
      res <- result()
      valueBox(format(res$n_after, big.mark = ","), "Genes retained", icon = icon("check"), color = "green")
    })
    output$vb_genes_dropped <- renderValueBox({
      res <- result()
      valueBox(format(res$n_before - res$n_after, big.mark = ","), "Genes filtered out", icon = icon("filter"), color = "purple")
    })
    output$vb_flagged <- renderValueBox({
      res <- result()
      n_flagged <- sum(res$qc$flag_signal | res$qc$flag_detected | res$qc$flag_cor)
      valueBox(n_flagged, "Samples flagged", icon = icon("triangle-exclamation"),
                color = if (n_flagged > 0) "red" else "green")
    })

    output$decisions_ui <- renderUI({
      res <- result()
      tagList(
        p(icon("circle-info"), " Normalisation: ", strong(res$norm_label), "."),
        p(icon("circle-info"), " Batch correction: ",
          if (isTRUE(res$skip_combat)) strong("skipped, normalisation only") else tagList(
            strong(switch(res$correction_method %||% "combat",
              combat_seq = "ComBat-seq, on raw counts before TMM normalisation",
              limma = "limma::removeBatchEffect",
              sva = "Surrogate Variable Analysis (unknown/hidden sources), regressed out via limma::removeBatchEffect",
              if (identical(res$combat_prior, "param")) "ComBat, parametric prior" else "ComBat, non-parametric prior"
            )),
            if (identical(res$correction_method, "combat") && isTRUE(res$combat_mean_only)) ", mean only" else "",
            if (isTRUE(res$use_batch2)) " on an interaction of the two chosen batch columns" else "",
            if (!is.null(res$ref_batch)) paste0(", referenced to batch \"", res$ref_batch, "\"") else ""
          ), "."),
        if (!isTRUE(res$skip_combat)) p(icon("circle-info"), " Protected ",
          if (length(res$protect) > 0) paste(res$protect, collapse = " and ") else "no biological covariates, since none were selected, or none had two or more levels in the current metadata",
          ", so that signal is not removed as if it were batch."),
        if (length(res$protect_dropped_for_batch) > 0) p(icon("triangle-exclamation"), " ",
          strong(paste(res$protect_dropped_for_batch, collapse = " and ")),
          if (length(res$protect_dropped_for_batch) > 1) " were" else " was",
          " also selected as the batch column (or part of it), so ", if (length(res$protect_dropped_for_batch) > 1) "they were" else "it was",
          " dropped from the protected covariates instead - a column can't be both corrected for and protected from correction at the same time."),
        if (isTRUE(res$n_excluded_outliers > 0)) p(icon("circle-info"), " Excluded ", strong(res$n_excluded_outliers),
          " sample(s) flagged as QC outliers before correcting.")
      )
    })

    ## ---- QC plots -----------------------------------------------------------

    output$signal_plot <- renderPlot({
      res <- result()
      qc_bar_plot(res$qc, "signal", "flag_signal", "Total signal")
    })
    output$detected_plot <- renderPlot({
      res <- result()
      qc_bar_plot(res$qc, "detected", "flag_detected", "Detected features")
    })
    output$cor_plot <- renderPlot({
      res <- result()
      qc_bar_plot(res$qc, "mean_cor", "flag_cor", "Mean correlation")
    })

    dist_summary <- function(expr, meta) {
      qs <- apply(expr, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      df <- data.frame(sample = colnames(expr), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ],
                        upper = qs[4, ], ymax = qs[5, ])
      merge(df, meta[, intersect(c("sample", "group"), colnames(meta))], by = "sample", all.x = TRUE)
    }

    output$dist_plot <- renderPlot({
      res <- result()
      before_df <- dist_summary(res$expr_prenorm, res$meta); before_df$stage <- "Before normalisation"
      after_df  <- dist_summary(res$expr_qnorm, res$meta);   after_df$stage  <- "After normalisation"
      box_df <- rbind(before_df, after_df)
      box_df$stage <- factor(box_df$stage, levels = c("Before normalisation", "After normalisation"))
      if (!"group" %in% colnames(box_df)) box_df$group <- "all samples"

      ggplot(box_df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle,
                           upper = upper, ymax = ymax, fill = group)) +
        geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15) +
        scale_fill_manual(values = arthomix_pair(box_df$group)) +
        facet_wrap(~stage, ncol = 1, scales = "free_x") +
        labs(x = NULL, y = "Expression", fill = NULL) +
        theme_arthomix() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
    })

    output$pca_before <- renderPlot({
      res <- result()
      plot_pca_advanced(res$before, res$meta, res$color_by,
                          pc_x = as.integer(input$pc_x), pc_y = as.integer(input$pc_y),
                          title_suffix = "Before batch correction",
                          show_ellipse = input$show_ellipse, show_labels = input$show_labels)
    })
    output$pca_after <- renderPlot({
      res <- result()
      suffix <- if (isTRUE(res$skip_combat)) "After batch correction (skipped - identical to \"before\")"
                else paste0("After batch correction (vs. ", res$batch_col, ")")
      plot_pca_advanced(res$after, res$meta, res$color_by,
                          pc_x = as.integer(input$pc_x), pc_y = as.integer(input$pc_y),
                          title_suffix = suffix,
                          show_ellipse = input$show_ellipse, show_labels = input$show_labels)
    })
    output$scree_plot <- renderPlot({
      res <- result()
      scree_plot(res$after$var_exp)
    })

    assoc_pvalue <- function(pc1, batch) {
      df <- data.frame(pc1 = pc1, b = batch)
      df <- df[!is.na(df$b), ]
      tryCatch(summary(aov(pc1 ~ b, data = df))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    }

    output$summary_ui <- renderUI({
      res <- result()
      before_p <- assoc_pvalue(res$before$df$PC1, res$meta[[res$batch_col]])
      after_p  <- assoc_pvalue(res$after$df$PC1,  res$meta[[res$batch_col]])
      if (isTRUE(res$skip_combat)) {
        return(tagList(
          p("PC1 versus ", strong(res$batch_col), " association (one-way ANOVA p-value): ",
            sprintf("%.2g", before_p), " before, and ", sprintf("%.2g", after_p), " after."),
          p(class = "empty-note", icon("triangle-exclamation"),
            "These two numbers are identical because \"Skip batch correction\" is ticked in Settings, so \"after\" is the same uncorrected data as \"before\". Untick it to see whether correction actually reduces this column's effect on PC1.")
        ))
      }
      tagList(
        p("PC1 versus ", strong(res$batch_col), " association (one-way ANOVA p-value): ",
          sprintf("%.2g", before_p), " before correction, and ", sprintf("%.2g", after_p), " after."),
        p(class = "empty-note", icon("circle-info"),
          "A larger p-value after correction means PC1 is less explained by that column: the correction reduced its effect.")
      )
    })

    ## ---- Normalisation diagnostics table --------------------------------

    output$norm_table <- DT::renderDataTable({
      res <- result()
      DT::datatable(res$norm_diag, rownames = FALSE,
                     options = list(pageLength = 4, dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatRound(columns = c("max_value", "min_value", "median_sd", "iqr_sd", "median_range", "iqr_range"), digits = 3)
    })

    output$download_norm <- downloadHandler(
      filename = function() "normalization_diagnostics.csv",
      content = function(file) write.csv(result()$norm_diag, file, row.names = FALSE)
    )

    ## ---- QC table -----------------------------------------------------------

    qc_table_display <- reactive({
      res <- result()
      df <- res$qc
      df$reason <- apply(df[, c("flag_signal", "flag_detected", "flag_cor")], 1, function(r) {
        reasons <- c("low or high signal", "low or high detected features", "low cohort correlation")[r]
        if (length(reasons) == 0) "" else paste(reasons, collapse = "; ")
      })
      df$signal <- round(df$signal, 1)
      df$mean_cor <- round(df$mean_cor, 3)
      df[, c("sample", intersect("group", colnames(df)), "signal", "detected", "mean_cor", "reason")]
    })

    output$qc_table <- DT::renderDataTable({
      df <- qc_table_display()
      flagged_only <- df[df$reason != "", , drop = FALSE]
      DT::datatable(
        if (nrow(flagged_only) > 0) flagged_only else df[0, ],
        rownames = FALSE, filter = "top",
        options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact"
      )
    })

    output$download_qc <- downloadHandler(
      filename = function() "qc_metrics.csv",
      content = function(file) write.csv(qc_table_display(), file, row.names = FALSE)
    )

    ## ---- PCA table -----------------------------------------------------------

    pca_table <- reactive({
      res <- result()
      rename_pc <- function(df, suffix) {
        pc_cols <- setdiff(colnames(df), "sample")
        colnames(df)[match(pc_cols, colnames(df))] <- paste0(pc_cols, "_", suffix)
        df
      }
      merge(rename_pc(res$before$df, "before"), rename_pc(res$after$df, "after"), by = "sample")
    })

    output$pca_table <- DT::renderDataTable({
      DT::datatable(pca_table(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_pca <- downloadHandler(
      filename = function() "pca_before_after.csv",
      content = function(file) write.csv(pca_table(), file, row.names = FALSE)
    )

    ## ---- Activate result app-wide -------------------------------------------

    output$activate_ui <- renderUI({
      box(width = 12, title = "Use this dataset app-wide", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Once you're happy with the result above, make it the active dataset for every other sub-module (WGCNA, differential expression, feature selection, etc.)."),
          actionButton(ns("activate_btn"), "Use this as the active dataset", icon = icon("check"), class = "btn-primary btn-sm"),
          uiOutput(ns("activate_status_ui")))
    })

    observeEvent(input$activate_btn, {
      res <- result()
      dataset$expr <- res$expr_combat
      dataset$meta <- res$meta
      dataset$source <- paste0("Uploaded dataset (preprocessed + ",
                                if (isTRUE(res$skip_combat)) "normalised" else "batch-corrected", "): ",
                                res$sources %||% "your data")
      output$activate_status_ui <- renderUI(
        div(class = "empty-note", icon("check"), "This is now the active dataset. Every other sub-module will use it.")
      )
    }, ignoreInit = TRUE)

    ## =====================================================================
    ## Assemble the page
    ## =====================================================================

    ## Nothing below the settings panel is meaningful until result() has
    ## actually run once, so none of it renders until then - no empty boxes,
    ## no spinners stuck with nothing to spin for.
    ## Plain, unboxed section - a small heading rule plus its content -
    ## used throughout Batch Correction's results INSTEAD of shinydashboard
    ## box() cards, per explicit feedback that stacked bordered cards read
    ## as an unprofessional dashboard look here. Batch Correction only;
    ## every other tab/module still uses box() as before.
    bc_section <- function(icon_name, title, ..., desc = NULL) {
      tagList(
        tags$h4(
          style = "margin: 20px 0 4px 0; padding-top: 14px; border-top: 1px solid var(--color-border); font-size: 15px; font-weight: 600; color: var(--color-ink); display: flex; align-items: center; gap: 8px;",
          icon(icon_name), title
        ),
        if (!is.null(desc)) p(class = "submodule-desc", desc),
        ...
      )
    }

    output$results_top_ui <- renderUI({
      res <- tryCatch(result(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"),
            "Set the options on the left, then click \"Run normalisation and batch correction\" to see results here."))
      }
      tagList(
        ## 2-per-row, not 4: this row lives inside a column(6) beside
        ## Settings, so 4-across would squeeze each value box to an
        ## eighth of the page.
        fluidRow(
          valueBoxOutput(ns("vb_samples"), width = 6),
          valueBoxOutput(ns("vb_genes_kept"), width = 6)
        ),
        fluidRow(
          valueBoxOutput(ns("vb_genes_dropped"), width = 6),
          valueBoxOutput(ns("vb_flagged"), width = 6)
        ),
        bc_section("clipboard-list", "Pipeline decisions",
          withSpinner(uiOutput(ns("decisions_ui")), color = "#2563EB", type = 6)
        )
      )
    })

    output$results_rest_ui <- renderUI({
      res <- tryCatch(result(), error = function(e) NULL)
      req(res)
      tagList(
        ## Stacked full-width, not nested into 3- and 2-column sub-rows:
        ## this whole tagList already lives inside a column(6) beside
        ## Settings (see batch_content below), and subdividing that further
        ## squeezed each plot down to roughly a sixth of the page - visibly
        ## clipped legends and unreadable per-sample bar plots. Full width
        ## of the results column is what these need to render legibly.
        bc_section("chart-simple", "Per-sample signal",
          withSpinner(plotOutput(ns("signal_plot"), height = 210), color = "#2563EB", type = 6)
        ),
        bc_section("chart-simple", "Detected features per sample",
          withSpinner(plotOutput(ns("detected_plot"), height = 210), color = "#2563EB", type = 6)
        ),
        bc_section("chart-simple", "Mean correlation to cohort",
          withSpinner(plotOutput(ns("cor_plot"), height = 210), color = "#2563EB", type = 6)
        ),
        bc_section("chart-area", "Expression distribution: before and after normalisation",
          withSpinner(plotOutput(ns("dist_plot"), height = 300), color = "#2563EB", type = 6)
        ),
        bc_section("chart-line", "Scree plot",
          desc = "Percentage of variance explained per component, after batch correction. A sharp drop after PC2 or PC3 means most of the structure is captured by the axes plotted below.",
          withSpinner(plotOutput(ns("scree_plot"), height = 220), color = "#2563EB", type = 6)
        ),
        if (isTRUE(res$skip_combat)) div(class = "empty-note", icon("triangle-exclamation"),
          strong(" \"Skip batch correction\" is ticked in Settings: "),
          "no correction was run, so the panels below are the same normalised data twice, and every before/after comparison is necessarily identical. Untick it and re-run to see an actual before/after."),
        bc_section("braille", "PCA before batch correction",
          withSpinner(plotOutput(ns("pca_before"), height = 340), color = "#2563EB", type = 6)
        ),
        bc_section("braille", "PCA after batch correction",
          withSpinner(plotOutput(ns("pca_after"), height = 340), color = "#2563EB", type = 6)
        ),
        bc_section("scale-balanced", "Batch effect, before and after",
          withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6)
        ),
        bc_section("table-list", "Normalisation diagnostics",
          desc = "Per-sample median and IQR agreement across the cohort. Large 'after' values mean normalisation didn't fully align the samples.",
          div(class = "table-toolbar", downloadButton(ns("download_norm"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("norm_table"))
        ),
        bc_section("triangle-exclamation", "Flagged samples",
          desc = "Samples flagged by signal, detected-feature count, or cohort correlation, with the reason.",
          div(class = "table-toolbar", downloadButton(ns("download_qc"), "Download full QC table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("qc_table"))
        ),
        bc_section("table-cells", "PCA coordinates",
          div(class = "table-toolbar", downloadButton(ns("download_pca"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("pca_table"))
        ),
        uiOutput(ns("activate_ui"))
      )
    })

    ## Settings (column(6), unboxed - a plain section header, matching the
    ## unboxed style bc_section() gives every results section below, so
    ## nothing in this tab looks visually inconsistent with the rest) beside
    ## just the compact top-of-results value boxes + decisions
    ## (results_top_ui) in the other column(6) - the same "settings beside
    ## the first result card" row mod_dge.R uses (Contrast beside Result).
    ## results_rest_ui (the plots and tables - the bulk of the page) breaks
    ## out of that row to full width below, exactly like mod_dge.R's
    ## Heatmap/Result table/References boxes do beneath its own settings+
    ## result row - previously it was squeezed into the same column(6) as
    ## Settings, clipping every plot and table to a quarter of the page.
    ## Pipeline summary is removed from this tab entirely so its space goes
    ## to Results - it is not shown elsewhere either (it was already pulled
    ## out of the shared rail in mod_preprocessing_ui in an earlier change);
    ## ArthOChat in that rail is untouched throughout.
    batch_content <- tagList(
      fluidRow(
        column(
          6,
          bc_section("sliders", "Settings", uiOutput(ns("settings_ui")))
        ),
        column(
          6,
          withSpinner(uiOutput(ns("results_top_ui")), color = "#2563EB", type = 6)
        )
      ),
      withSpinner(uiOutput(ns("results_rest_ui")), color = "#2563EB", type = 6)
    )

    output$batch_tab_ui <- renderUI({
      tagList(
        div(class = "empty-note", icon("circle-info"),
            "Runs on whatever you merged in the Merge datasets tab (or your single preprocessed dataset, if you only configured one)."),
        batch_content
      )
    })
  })
}
