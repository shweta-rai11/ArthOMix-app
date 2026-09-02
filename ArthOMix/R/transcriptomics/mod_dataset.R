## R/mod_dataset.R
## Dataset tab: pick a preloaded dataset, upload your own, or fetch from NCBI
## GEO - each of the three is an independent pipeline, and whichever one you
## load becomes dataset$expr/meta/source/source_type immediately (every
## sub-module below reads it), not just dataset$staged_* (kept in parallel as
## the preview shown on this tab and as the input Preprocessing's "Currently
## loaded dataset" option reads before any merge/batch-correction step runs).

mod_dataset_config <- list(
  id = "dataset", 
  title = "Dataset", 
  icon = "database",
  description = "Pick a preloaded dataset or upload your own - either way it's what every sub-module below reads from."
)

## Tissue-first display names for the preloaded dropdown - no GEO accession shown;
## falls back to the raw GSE ID if a source has no entry here.
INDIVIDUAL_DATASET_LABELS <- c(
  "GSE93272"  = "Whole Blood Training Cohort A",
  "GSE110169" = "Whole Blood Training Cohort B",
  "GSE15573"  = "PBMC Validation Cohort",
  "GSE89408"  = "Synovial Tissue Validation Cohort"
)

individual_dataset_entry <- function(gse_id) {
  list(
    id = gse_id,
    label = if (gse_id %in% names(INDIVIDUAL_DATASET_LABELS)) INDIVIDUAL_DATASET_LABELS[[gse_id]] else gse_id,
    load = function() {
      d <- load_individual_dataset(gse_id)
      validate(need(!is.null(d), paste("Raw data for", gse_id, "was not found on disk.")))
      list(expr = d$expr, meta = d$meta, source = paste0("Individual dataset: ", d$label))
    }
  )
}

## Reloads the same merged, batch-corrected cohort loaded at startup, so switching
## to an individual raw dataset below isn't a one-way trip.
default_dataset_entry <- list(
  id = "__default_merged__",
  label = "Merged Data",
  load = function() {
    d <- load_default_dataset()
    list(expr = d$expr, meta = d$meta, source = d$source)
  }
)

## Catalog: the merged/batch-corrected default cohort plus each individual raw GEO source.
PRELOADED_DATASETS <- c(
  list(default_dataset_entry),
  lapply(vapply(GEO_SOURCES, `[[`, character(1), "gse"), individual_dataset_entry)
)

preloaded_choices <- function() {
  setNames(vapply(PRELOADED_DATASETS, `[[`, character(1), "id"),
           vapply(PRELOADED_DATASETS, `[[`, character(1), "label"))
}

## The two GEO series load_default_dataset()'s merged cohort is built from training data

MERGED_DEFAULT_GEO_IDS <- c("GSE93272", "GSE110169")

entry_geo_ids <- function(entry) {
  if (identical(entry$id, default_dataset_entry$id)) MERGED_DEFAULT_GEO_IDS else entry$id
}

mod_dataset_ui <- function(id) {
  ns <- NS(id) # NS means namespace
  tagList(  # returns multiple UI components together
    fluidRow(  # for bootstrap grid layout
      column(
        6,
        box( # Box is a UI container in shiny dashboard
          width = NULL, title = "Preloaded Datasets", status = "primary", solidHeader = FALSE,
          selectInput(ns("preloaded_choice"), "Individual dataset",
                      choices = preloaded_choices(), selected = character(0), width = "100%"),
          uiOutput(ns("preloaded_note")),
          uiOutput(ns("preloaded_geo_card_ui")),
          div(style = "display: flex; align-items: center; gap: 10px; flex-wrap: wrap;",
              actionButton(ns("load_preloaded_btn"), "Load this dataset", icon = icon("rotate-left"), class = "btn-primary btn-sm"),
              uiOutput(ns("preloaded_load_message"), inline = TRUE))
        ),
        box(
          width = NULL, title = "Fetch from NCBI GEO", status = "primary", solidHeader = FALSE,
          p(strong("GEO Series accession"), " (e.g. ", code("GSE12345"), "). Fetches the expression data and metadata from NCBI. For your own data, use \"Upload your own data\" instead."),
          fluidRow(
            column(4, textInput(ns("geo_accession"), NULL, placeholder = "GSE12345", width = "100%")),
            column(3, actionButton(ns("geo_fetch_btn"), "Fetch", icon = icon("download"), class = "btn-primary btn-sm", width = "100%"))
          ),

          tags$script(HTML(sprintf(
            "document.addEventListener('click', function(e) {
               var btn = e.target.closest('#%s');
               if (!btn) return;
               var statusDiv = document.getElementById('%s');
               if (!statusDiv) return;
               btn.disabled = true;
               statusDiv.innerHTML = '<div class=\"empty-note\"><i class=\"fas fa-spinner fa-spin\" role=\"presentation\"></i> Fetching from NCBI GEO - this can take anywhere from a few seconds to a couple of minutes depending on the series size and NCBI\\'s current load. This app is unresponsive to other actions while it fetches - please don\\'t navigate away or click Fetch again.</div>';
               var observer = new MutationObserver(function() {
                 var liveBtn = document.getElementById('%s');
                 if (liveBtn) liveBtn.disabled = false;
                 observer.disconnect();
               });
               observer.observe(statusDiv, { childList: true, subtree: true });
             });",
            ns("geo_fetch_btn"), ns("geo_fetch_status"), ns("geo_fetch_btn")
          ))),
          uiOutput(ns("geo_fetch_status")),
          uiOutput(ns("geo_platform_ui")),
          uiOutput(ns("geo_column_mapping")),
          actionButton(ns("geo_load_btn"), "Load this dataset", icon = icon("upload"), class = "btn-primary btn-sm")
        )
      ),
      column(
        6,
        box(
          width = NULL, title = "Upload your own data", status = "primary", solidHeader = FALSE,
          div(class = "upload-step-label", "STEP 1 · Choose your files"),
          p(strong("Expression matrix"), " - CSV or RDS. Genes in rows, samples in columns; for CSV, the first column is the gene ID."),
          radioButtons(ns("declared_data_type"), "Data type", inline = TRUE,
                       choices = c("Raw counts" = "raw", "Normalized (TPM/FPKM/CPM)" = "normalized", "Already log-transformed" = "logtransformed"),
                       selected = "normalized"),
          fileInput(ns("expr_file"), "Expression matrix", accept = c(".csv", ".rds", ".Rds")),
          p(strong("Sample metadata"), " - CSV or RDS data frame, one row per sample."),
          fileInput(ns("meta_file"), "Sample metadata", accept = c(".csv", ".rds", ".Rds")),
          uiOutput(ns("upload_preview_ui")),
          uiOutput(ns("upload_preview_tables_ui")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 2 · Map the columns"),
          uiOutput(ns("column_mapping")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 3 · Confirm"),
          actionButton(ns("load_btn"), "Upload Data", icon = icon("upload"), class = "btn-primary btn-sm"),
          div(class = "empty-note", style = "margin-top:6px;", icon("circle-info"),
              "Loads exactly what you provide, as-is - no merging, normalising, or batch correction.")
        )
      )
    ),
    uiOutput(ns("load_message"))
  )
}

mod_dataset_server <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$preloaded_note <- renderUI({
      req(input$preloaded_choice)
      if (identical(input$preloaded_choice, default_dataset_entry$id)) {
        return(p(class = "empty-note", icon("circle-info"),
          "The same merged, batch-corrected training cohort the app loads by default on startup - pick this to switch back to it after loading something else."))
      }
      ## GSE93272/GSE110169 only fall back to samples filtered out of the
      ## merged cohort (via load_individual_dataset()'s own
      ## merged_training_subset() fallback in global.R) when their raw
      ## ExpressionSet isn't actually on disk in this deployment - checked
      ## live here (matching preloaded_geo_card_ui's own get_raw_eset() check
      ## just below, and mod_overview.R's identical card) rather than assumed,
      ## so this note can't drift out of sync with which raw files this
      ## deployment actually ships.
      raw_unavailable <- input$preloaded_choice %in% c("GSE93272", "GSE110169") &&
        {
          eset <- get_raw_eset(input$preloaded_choice)
          is.null(eset) || nrow(Biobase::exprs(eset)) == 0
        }
      if (raw_unavailable) {
        p(class = "empty-note", icon("circle-info"),
          "This source's raw probe-level file isn't available in this deployment, so this loads its samples only, filtered out of the merged, batch-corrected training cohort - not raw, single-platform data. To see it merged with the other training source instead, pick \"Merged Data\" above.")
      } else {
        p(class = "empty-note", icon("triangle-exclamation"),
          "Raw, single-platform data - probe-level, not merged or normalised. You can run any sub-module directly against it, but most were built assuming the merged cohort, so results may look different. To just look at it without changing what every sub-module runs on, use the QC tab on Overview and Datasets instead.")
      }
    })

    ## GEO info-card for the 4 individual sources (not "Merged Data"), matching the
    ## card the Overview and Datasets tab uses for the same sources.
    output$preloaded_geo_card_ui <- renderUI({
      req(input$preloaded_choice)
      if (identical(input$preloaded_choice, default_dataset_entry$id)) return(NULL)
      gse_id <- input$preloaded_choice
      src <- Find(function(s) identical(s$gse, gse_id), GEO_SOURCES)
      req(src)
      eset <- get_raw_eset(gse_id)
      div(
        class = "info-card",
        div(
          class = "module-card-title-row",
          h4(gse_id),
          tags$a(href = geo_link(gse_id), target = "_blank", rel = "noopener",
                  icon("up-right-from-square"), " NCBI GEO")
        ),
        if (!is.null(eset)) {
          tagList(
            p(class = "module-card-tagline",
              tryCatch(Biobase::experimentData(eset)@title, error = function(e) NULL)),
            p(strong("Role: "), src$role, br(), strong("Used for: "), src$used_in),
            p(strong("Platform: "), Biobase::annotation(eset), br(),
              strong("Samples: "), ncol(eset), ", ", strong("Probes: "), format(nrow(eset), big.mark = ","))
          )
        } else {
          tagList(
            p(strong("Role: "), src$role, br(), strong("Used for: "), src$used_in),
            div(class = "empty-note", icon("triangle-exclamation"), "Raw file not found on disk.")
          )
        }
      )
    })

    meta_raw <- reactive({
      req(input$meta_file)
      path <- input$meta_file$datapath
      if (grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        d <- loaded$value
        validate(need(is.data.frame(d), "The uploaded metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    ## Parses the expression file once; preview and load_btn both read from this
    ## instead of re-parsing a potentially large CSV.
    expr_raw <- reactive({
      req(input$expr_file)
      path <- input$expr_file$datapath
      if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        loaded$value
      } else {
        m <- as.data.frame(data.table::fread(path, showProgress = FALSE))
        rn <- as.character(m[[1]])
        m <- as.matrix(m[, -1, drop = FALSE])
        ## as.matrix() on a data.frame with any non-numeric data cell (stray
        ## text, a thousands-separated "1,234", a footnote) silently coerces
        ## the whole matrix to character (or, for an all-blank/0-row file,
        ## logical) instead of erroring - that would otherwise surface later
        ## as an opaque low-level error out of tx_validate_expr_upload(). Only
        ## block on genuine non-numeric TEXT here; an all-NA/empty matrix is
        ## intentionally let through so tx_validate_expr_upload's own "No
        ## finite numeric values" check - a more specific message - catches
        ## it downstream instead.
        if (!is.numeric(m)) {
          storage.mode(m) <- "character"
          m_num <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m)))
          validate(need(!any(!is.na(m) & m != "" & is.na(m_num)),
            "This file has non-numeric values outside the first (gene ID) column. Expression matrices must be purely numeric after the ID column - check for stray text, footnotes, or thousands separators in the data cells."))
          m <- m_num
        }
        rownames(m) <- rn
        m
      }
    })

    ## Fires as soon as both files are readable, before Load is clicked, so a large
    ## upload doesn't look stuck; shared by the summary line and preview tables below.
    upload_preview_data <- reactive({
      req(input$expr_file, input$meta_file)
      tryCatch(list(expr = expr_raw(), meta = meta_raw()), error = function(e) e)
    })

    output$upload_preview_ui <- renderUI({
      preview <- upload_preview_data()
      if (inherits(preview, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not read the uploaded file(s):", conditionMessage(preview))))
      }
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s features x %s samples. Read %s: %s rows. Map the columns below, then click \"Load dataset\".",
                  input$expr_file$name, format(nrow(preview$expr), big.mark = ","), ncol(preview$expr),
                  input$meta_file$name, nrow(preview$meta)))
    })

    ## Shows the first rows of the uploaded metadata/expression matrix as read,
    ## before any column mapping.
    output$upload_preview_tables_ui <- renderUI({
      preview <- upload_preview_data()
      req(!inherits(preview, "error"))
      tagList(
        div(class = "upload-step-label", "Preview of what you uploaded"),
        p(class = "submodule-desc", "First 5 rows/columns, exactly as read - before any column mapping below."),
        strong("Sample metadata (first 5 rows)"),
        DT::dataTableOutput(ns("upload_preview_meta_table")),
        br(),
        strong("Expression matrix (first 5 features x 5 samples)"),
        DT::dataTableOutput(ns("upload_preview_expr_table"))
      )
    })

    output$upload_preview_meta_table <- DT::renderDataTable({
      preview <- upload_preview_data()
      req(!inherits(preview, "error"))
      DT::datatable(head(preview$meta, 5), rownames = FALSE,
                     options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$upload_preview_expr_table <- DT::renderDataTable({
      preview <- upload_preview_data()
      req(!inherits(preview, "error"))
      m <- preview$expr
      cols <- seq_len(min(5, ncol(m)))
      df <- data.frame(feature = rownames(m), round(m[, cols, drop = FALSE], 3), check.names = FALSE)
      DT::datatable(head(df, 5), rownames = FALSE,
                     options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    ## Guesses a mapping dropdown's default column by name, not position: each
    ## exact-match term tried in order, then each substring term tried in
    ## order, then falls back to the first column. Each tier must go
    ## term-by-term (not one combined OR regex checked column-by-column) - a
    ## combined regex returns whichever candidate COLUMN happens to come
    ## first in the data, not whichever TERM the caller ranked highest, and
    ## for a real GEO series that's often wrong (e.g. an administrative
    ## "status" column - not the caller's 5th-priority term but simply the
    ## column GEO happened to place earliest - beating the actual "disease
    ## state:ch1" column for the group/diagnosis guess).
    ## Return value carries a "matched" attribute (TRUE = a real heuristic hit,
    ## FALSE = the positional fallback was used) so required-field callers
    ## below can warn the user instead of silently pre-selecting an arbitrary
    ## column - a native (selectize = FALSE) <select> always shows some
    ## option selected, so "no default" isn't renderable here; flagging the
    ## guess as unconfident is the achievable alternative.
    guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
      for (term in exact) {
        hit <- cols[tolower(cols) == tolower(term)]
        if (length(hit) > 0) return(structure(hit[1], matched = TRUE))
      }
      for (term in contains) {
        hit <- cols[grepl(term, cols, ignore.case = TRUE)]
        if (length(hit) > 0) return(structure(hit[1], matched = TRUE))
      }
      structure(fallback, matched = FALSE)
    }

    ## Small "double check this" note shown under a required-field dropdown
    ## whose selection is only the positional fallback, not a real name-based
    ## guess - see guess_col()'s header comment.
    unconfident_guess_note <- function(guess, label) {
      if (isTRUE(attr(guess, "matched"))) return(NULL)
      div(class = "empty-note", icon("triangle-exclamation"),
          sprintf("Couldn't confidently guess the %s from its column name - double-check the selection above.", label))
    }

    output$column_mapping <- renderUI({
      req(input$meta_file)
      cols <- colnames(meta_raw())
      id_guess <- guess_col(cols, c("sample", "sample_id", "id", "geo_accession", "accession"))
      group_guess <- guess_col(cols, c("group", "diagnosis", "disease", "condition", "status", "phenotype"))
      tagList(
        selectInput(ns("map_id"), "Sample ID column", choices = cols, selected = id_guess, selectize = FALSE),
        unconfident_guess_note(id_guess, "sample-ID column"),
        selectInput(ns("map_group"), "Group / diagnosis column", choices = cols, selected = group_guess, selectize = FALSE),
        unconfident_guess_note(group_guess, "group/diagnosis column"),
        selectInput(ns("map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "cohort", "platform", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })

    ## Keeps "Upload Data" disabled until files and required column mappings are set.
    observe({
      ready <- !is.null(input$expr_file) && !is.null(input$meta_file) &&
        !is.null(input$map_id) && !is.null(input$map_group)
      if (isTRUE(ready)) shinyjs::enable("load_btn") else shinyjs::disable("load_btn")
    })

    ## Mirrors the enable/disable guards on load_btn/geo_load_btn below, so
    ## all three "Load this dataset" buttons behave the same way instead of
    ## this one being always-clickable while unset (a blank preloaded_choice
    ## used to be a silent no-op with no feedback at all).
    observe({
      shinyjs::toggleState("load_preloaded_btn", condition = !is.null(input$preloaded_choice))
    })

    ## Stages + immediately activates a freshly loaded dataset - the seven
    ## fields below used to be hand-copied identically into each of the three
    ## load handlers (preloaded/upload/GEO), risking the three copies
    ## silently drifting apart. Also runs the duplicate-feature-ID check that
    ## used to run only on the upload path, so all three paths warn the same
    ## way about a condition the Preprocessing merge tab silently papers over
    ## (keeps only the first occurrence of a duplicate row name).
    activate_dataset <- function(expr, meta, source, source_type, is_bundled_reference,
                                  geo_ids, declared_data_type = NA_character_) {
      dataset$staged_expr <- expr
      dataset$staged_meta <- meta
      dataset$staged_source <- source
      dataset$expr <- expr
      dataset$meta <- meta
      dataset$source <- source
      dataset$source_type <- source_type
      dataset$is_bundled_reference <- is_bundled_reference
      dataset$geo_ids <- geo_ids
      dataset$declared_data_type <- declared_data_type
      sum(duplicated(rownames(expr)))
    }

    duplicate_feature_note <- function(n_dup) {
      if (n_dup == 0) return(NULL)
      div(class = "empty-note", icon("triangle-exclamation"),
          sprintf("%d duplicated feature identifier(s) were detected in this dataset. All rows are kept here, but downstream row-name-keyed steps (e.g. the Preprocessing merge tab) will keep only the first occurrence of each - rename duplicates in your source file if this is unintended.", n_dup))
    }

    observeEvent(input$load_preloaded_btn, {
      req(input$preloaded_choice)
      entry <- Find(function(d) d$id == input$preloaded_choice, PRELOADED_DATASETS)
      req(entry)
      d <- entry$load()
      ## is_bundled_reference is TRUE only for the exact default merged
      ## cohort - the one the app's bundled/precomputed reference-panel CSVs
      ## (FS_input_*.csv etc.) were actually computed from - not for the 4
      ## individual raw GSE picks. No declaration UI for preloaded data, so
      ## declared_data_type is left at its NA default - clears any stale
      ## value left over from a previous upload (see the upload handler
      ## below), falling back to heuristic-only inference downstream.
      n_dup <- activate_dataset(
        expr = d$expr, meta = d$meta, source = d$source,
        source_type = "preloaded",
        is_bundled_reference = identical(entry$id, default_dataset_entry$id),
        geo_ids = entry_geo_ids(entry)
      )
      output$preloaded_load_message <- renderUI(tagList(
        span(style = "color: var(--color-success); font-size: 13px; font-weight: 600;", icon("check"), " ",
             sprintf("Loaded %s: %s genes across %s samples. Every sub-module now runs on this dataset - optionally go to Preprocessing and pick \"Currently loaded dataset\" to merge, normalise, or batch-correct it first.",
                      entry$label, format(nrow(d$expr), big.mark = ","), ncol(d$expr))),
        duplicate_feature_note(n_dup)
      ))
    })

    observeEvent(input$load_btn, {
      req(input$expr_file, input$meta_file, input$map_id, input$map_group)

      result <- tryCatch({
        expr <- expr_raw()
        ## Declare-then-verify (mirrors R/methylomics/parse_upload.R's
        ## methyl_validate_matrix_upload()): checks the "Data type" radio
        ## button above against the matrix's actual value range/shape, hard-
        ## blocking on an unambiguous mismatch (negative/TPM-pinned values
        ## declared "raw", raw-count-shaped values declared "normalized"/
        ## "logtransformed", or a differential-results table shaped upload)
        ## rather than only discovering the mismatch later, when
        ## Differential Expression's method gate runs at analysis time.
        checked <- tx_validate_expr_upload(expr, input$declared_data_type)
        validate(need(isTRUE(checked$ok), checked$error))
        expr <- checked$mat
        meta <- meta_raw()
        meta$sample <- as.character(meta[[input$map_id]])
        meta$group  <- as.character(meta[[input$map_group]])
        meta$sex    <- if (!identical(input$map_sex, "(none)")) as.character(meta[[input$map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$map_batch, "(none)")) as.character(meta[[input$map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(
          length(common) >= 4,
          "Fewer than 4 sample IDs in the expression matrix match the metadata sample-ID column. Check the column mapping."
        ))

        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        list(expr = expr, meta = meta, type_note = checked$note)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this dataset:", conditionMessage(result)))
        )
      } else {
        source_label <- paste0("Uploaded dataset: ", input$expr_file$name, " + ", input$meta_file$name)
        ## This is Pipeline 1's isolation moment: every downstream module now
        ## reads only this uploaded data. No known GEO provenance for an
        ## arbitrary upload, so geo_ids is empty - the Overview and Datasets
        ## tab's "Datasets" tab shows "No GEO ID found" for this.
        ## declared_data_type (raw/normalized/logtransformed) is read by
        ## Differential Expression's/Immune Deconvolution's run-time method
        ## gates in preference to live heuristic inference.
        n_dup <- activate_dataset(
          expr = result$expr, meta = result$meta, source = source_label,
          source_type = "uploaded", is_bundled_reference = FALSE,
          geo_ids = character(0), declared_data_type = input$declared_data_type
        )
        output$load_message <- renderUI(
          tagList(
            div(class = "empty-note", icon("check"),
                sprintf("Loaded %s genes across %s samples, as-is. Every sub-module now runs on this dataset.",
                        format(nrow(result$expr), big.mark = ","), ncol(result$expr))),
            duplicate_feature_note(n_dup),
            if (!is.null(result$type_note) && nzchar(result$type_note)) div(class = "empty-note", icon("triangle-exclamation"), result$type_note)
          )
        )
      }
    })

    ## Fetches one GEO Series as an ExpressionSet via GEOquery on Fetch click; returns
    ## an error condition object instead of raising, so it can be shown inline.
    geo_fetch_result <- eventReactive(input$geo_fetch_btn, {
      if (!requireNamespace("GEOquery", quietly = TRUE)) {
        return(simpleError("The GEOquery package is not installed in this deployment. Install it with BiocManager::install(\"GEOquery\") to enable fetching by GEO accession, or use \"Upload your own data\" instead."))
      }
      acc <- toupper(trimws(input$geo_accession %||% ""))
      if (!grepl("^GSE[0-9]+$", acc)) {
        return(simpleError("Enter a valid GEO Series accession, e.g. GSE12345."))
      }
      tryCatch({
        gpl_skipped <- FALSE
        ## GEOquery's default getGPL=TRUE also downloads the platform annotation file;
        ## when NCBI rate-limits/captchas that request, it hands back an HTML challenge
        ## page that GEOquery's SOFT parser chokes on with this same generic message even
        ## though the series matrix itself downloaded fine - retry once without it.
        ## "series_data_table_begin" is GEOquery's own SOFT-parser sentinel string
        ## (verified present in installed GEOquery 2.74.0) - a GEOquery upgrade that
        ## changes this wording would silently disable this retry; re-verify against
        ## GEOquery::getGEO's parser on any GEOquery version bump.
        gset <- tryCatch(
          suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE)),
          error = function(e) {
            if (!grepl("series_data_table_begin", conditionMessage(e), fixed = TRUE)) stop(e)
            gpl_skipped <<- TRUE
            suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE))
          }
        )
        if (!is.list(gset) || length(gset) == 0) {
          stop(paste(acc, "returned no series matrix from GEO - check the accession is a Series (GSExxxxx), not a Sample (GSM) or Platform (GPL) ID."))
        }
        list(acc = acc, platforms = gset, gpl_skipped = gpl_skipped)
      }, error = function(e) e)
    })

    output$geo_platform_ui <- renderUI({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) <= 1) return(NULL)
      selectInput(ns("geo_platform_choice"), "This series spans multiple platforms - pick one",
                  choices = names(res$platforms), width = "100%")
    })

    geo_eset <- reactive({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) > 1) {
        req(input$geo_platform_choice)
        res$platforms[[input$geo_platform_choice]]
      } else {
        res$platforms[[1]]
      }
    })

    ## Extracts expr/meta from the fetched ExpressionSet, collapsing probes to genes
    ## via collapse_probes_to_genes() when the platform annotation allows it.
    ##
    ## "Still waiting on a platform pick" for a multi-platform series is
    ## deliberately checked BEFORE calling geo_eset() (whose own
    ## req(input$geo_platform_choice) would otherwise throw a Shiny
    ## shiny.silent.error - which DOES inherit from "error", so the tryCatch
    ## below would catch it and re-surface it through the generic
    ## inherits(em, "error") branch in geo_fetch_status as a visible,
    ## near-blank warning banner instead of the silent "not ready yet" Shiny
    ## intends req() failures to produce). Returning NULL here instead lets
    ## req(em) downstream do the silent wait correctly.
    geo_expr_meta <- reactive({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) > 1 && is.null(input$geo_platform_choice)) return(NULL)
      tryCatch({
        eset <- geo_eset()
        ex <- Biobase::exprs(eset)
        if (nrow(ex) == 0 || ncol(ex) == 0) {
          stop("This GEO series has no expression matrix in its series matrix file - common for RNA-seq series that only deposit raw counts as supplementary files. Download that file from the GEO page and use \"Upload your own data\" instead.")
        }
        collapsed <- tryCatch(collapse_probes_to_genes(eset), error = function(e) NULL)
        used_collapse <- !is.null(collapsed) && nrow(collapsed) > 0 && isTRUE(attr(collapsed, "collapsed"))
        expr <- if (used_collapse) collapsed else ex
        list(expr = expr, meta = as.data.frame(Biobase::pData(eset)),
             collapsed = used_collapse, platform = Biobase::annotation(eset))
      }, error = function(e) e)
    })

    ## Compact fetch-confirmation line, not a full preview - once loaded, the
    ## detailed info card (title, platform, GEO link) lives on the Overview
    ## and Datasets sub-module's "Datasets" tab (dataset$geo_ids-driven, see
    ## mod_overview.R), and the full sample metadata table on its "Metadata"
    ## tab, rather than duplicating both here before the user has even
    ## clicked "Load this dataset".
    output$geo_fetch_status <- renderUI({
      req(input$geo_fetch_btn)
      res <- geo_fetch_result()
      if (inherits(res, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), paste("Could not fetch from GEO:", conditionMessage(res))))
      }
      em <- geo_expr_meta()
      req(em) # NULL = still waiting on a platform pick above; render nothing yet, not an error
      if (inherits(em, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), conditionMessage(em)))
      }
      tagList(
        div(class = "empty-note", icon("circle-info"),
            sprintf("Fetched %s (%s): %s samples x %s %s. Map the columns below, then click \"Load this dataset\".",
                    res$acc, em$platform, ncol(em$expr), format(nrow(em$expr), big.mark = ","),
                    if (em$collapsed) "genes" else "probes/features")),
        if (!em$collapsed) {
          div(class = "empty-note", icon("triangle-exclamation"),
              if (isTRUE(res$gpl_skipped))
                "GEO's platform-annotation file couldn't be fetched just now (NCBI may be rate-limiting automated requests) - loaded at probe/feature-ID level without gene symbols. Try fetching again later, or use Preprocessing's own probe-collapse step afterward."
              else
                "No gene-symbol annotation found for this platform - left at probe/feature-ID level. You can still load it as-is, or use Preprocessing's own probe-collapse step afterward.")
        }
      )
    })

    ## No sample-ID mapping needed here - pData() rownames are already the same
    ## GSM accessions exprs() is indexed by.
    output$geo_column_mapping <- renderUI({
      em <- geo_expr_meta()
      req(em); req(!inherits(em, "error"))
      cols <- colnames(em$meta)
      geo_group_guess <- guess_col(
        cols,
        exact = c("group", "diagnosis", "disease", "condition", "phenotype"),
        contains = c("group", "diagnosis", "disease", "condition", "phenotype", "characteristics_ch1")
      )
      tagList(
        ## "status" excluded, and "characteristics_ch1" moved to contains-only
        ## (never exact) - unlike the upload mapping's own guess_col() call
        ## below, GEO's own pData() always has a literal administrative
        ## "status" column ("Public on <date>") AND, for any series GEOquery
        ## didn't auto-expand into named "<key>:ch1" columns, a literal
        ## column named exactly "characteristics_ch1" (the raw, unparsed
        ## first characteristics field) - either would win the exact-match
        ## tier over the real disease/diagnosis column (almost always one of
        ## the *:ch1-suffixed fields, e.g. "disease state:ch1") for
        ## essentially every real GEO series fetched this way. Confirmed live
        ## against GSE93272 (has both a literal "status" and a literal
        ## "characteristics_ch1" column, alongside the correct
        ## "disease state:ch1").
        selectInput(ns("geo_map_group"), "Group / diagnosis column", choices = cols,
                    selected = geo_group_guess, selectize = FALSE),
        unconfident_guess_note(geo_group_guess, "group/diagnosis column"),
        selectInput(ns("geo_map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("geo_map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "cohort", "platform", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })

    observe({
      em <- tryCatch(geo_expr_meta(), error = function(e) e)
      ready <- !is.null(em) && !inherits(em, "error") && !is.null(input$geo_map_group)
      if (isTRUE(ready)) shinyjs::enable("geo_load_btn") else shinyjs::disable("geo_load_btn")
    })

    observeEvent(input$geo_load_btn, {
      req(input$geo_map_group)

      result <- tryCatch({
        em <- geo_expr_meta()
        validate(need(!inherits(em, "error"), "No GEO data fetched yet."))
        expr <- em$expr
        meta <- em$meta
        meta$sample <- rownames(meta)
        meta$group  <- as.character(meta[[input$geo_map_group]])
        meta$sex    <- if (!identical(input$geo_map_sex, "(none)")) as.character(meta[[input$geo_map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$geo_map_batch, "(none)")) as.character(meta[[input$geo_map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(
          length(common) >= 4,
          "Fewer than 4 sample IDs matched between the fetched expression matrix and metadata."
        ))
        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        acc <- geo_fetch_result()$acc
        label <- sprintf("%s (%s, %s)", acc, em$platform,
                          if (em$collapsed) "collapsed to genes" else "probe-level, raw")
        list(expr = expr, meta = meta, label = label, acc = acc)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this GEO dataset:", conditionMessage(result)))
        )
      } else {
        source_label <- paste0("NCBI GEO: ", result$label)
        ## This is Pipeline 2's isolation moment. No declaration UI for
        ## GEO-fetched data either - declared_data_type is left at its NA
        ## default, same stale-value guard as the preloaded handler above.
        n_dup <- activate_dataset(
          expr = result$expr, meta = result$meta, source = source_label,
          source_type = "geo", is_bundled_reference = FALSE, geo_ids = result$acc
        )
        n_samples <- ncol(result$expr)
        wgcna_note <- if (n_samples < 15) {
          " Note: WGCNA needs at least 15 samples to detect modules reliably, so its tab will stay blank until a larger dataset (e.g. merged with another GEO series or the preloaded cohort) is loaded."
        } else {
          ""
        }
        ## Echoes the pre-Load "no gene-symbol annotation" banner (see geo_fetch_status
        ## above) here too - a user who skimmed past that banner and clicked Load
        ## should still see, in the one message that persists after loading, that
        ## rows are raw probe/feature IDs rather than gene symbols.
        probe_note <- if (!isTRUE(em$collapsed)) " Note: no gene-symbol annotation was found for this platform, so rows are raw probe/feature IDs, not gene symbols." else ""
        output$load_message <- renderUI(tagList(
          div(class = "empty-note", icon("check"),
              sprintf("Loaded %s genes across %s samples. Every sub-module now runs on this GEO dataset.%s%s",
                      format(nrow(result$expr), big.mark = ","), n_samples, wgcna_note, probe_note)),
          duplicate_feature_note(n_dup)
        ))
      }
    })

  })
}
