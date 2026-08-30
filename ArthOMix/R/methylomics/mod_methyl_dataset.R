## R/methylomics/mod_methyl_dataset.R
## Methylomics Dataset tab: loads the shared `methyl_dataset` reactiveValues every
## Methylomics sub-module reads from. Mirrors the Transcriptomics Dataset tab's
## layout/interaction pattern (R/transcriptomics/mod_dataset.R) exactly - preloaded
## card + GEO-fetch card on the left, upload card (Step 1/2/3) on the right - but
## every source and validation is methylation-specific: beta/M-value matrices,
## probe (CpG) IDs, Illumina array platforms, IDAT raw files.

mod_methyl_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Pick a preloaded methylation dataset, upload your own, or fetch from NCBI GEO - either way it's what every sub-module below reads from."
)

## ---- Preloaded catalog -----------------------------------------------------
## Only one full, sub-module-ready preloaded matrix ships today (the live
## GSE42861 beta matrix load_default_meth_matrix() reads) - the GSE111942
## "external panel" bundled elsewhere is a curated 21-CpG diagnostic-classifier
## panel, not a full raw matrix, so it isn't offered here as a swappable pick.
## Structured as a list (like transcriptomics' PRELOADED_DATASETS) so a second
## full preloaded methylation cohort can be added later without a UI change.
MX_PRELOADED_DATASETS <- list(
  list(id = "gse42861_wholeblood", label = "Whole-blood sex-stratified RA cohort", gse = "GSE42861")
)

mx_preloaded_choices <- function() {
  setNames(vapply(MX_PRELOADED_DATASETS, `[[`, character(1), "id"),
           vapply(MX_PRELOADED_DATASETS, `[[`, character(1), "label"))
}

## ---- GEO platform recognition ----------------------------------------------
## GEO Platform (GPL) accession -> this app's array_type label, for every
## Illumina methylation array with any real GEO usage. A fetched series on an
## unrecognized platform falls back to a value-range heuristic (see
## mx_geo_expr_meta below) rather than being rejected outright, since custom/
## unlisted array platforms still submit ordinary beta-value matrices.
MX_METHYLATION_GPL <- c(
  "GPL8490"  = "27K",
  "GPL13534" = "450K",
  "GPL21145" = "EPIC",
  "GPL33022" = "EPICv2"
)

mod_methyl_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        6,
        box(
          width = NULL, title = "Switch to preloaded data", status = "primary", solidHeader = FALSE,
          selectInput(ns("preloaded_choice"), "Individual dataset",
                      choices = mx_preloaded_choices(), selected = character(0), width = "100%"),
          uiOutput(ns("preloaded_note")),
          uiOutput(ns("preloaded_geo_card_ui")),
          div(style = "display: flex; align-items: center; gap: 10px; flex-wrap: wrap;",
              actionButton(ns("load_preloaded_btn"), "Load this dataset", icon = icon("rotate-left"), class = "btn-primary btn-sm"),
              uiOutput(ns("preloaded_load_message"), inline = TRUE))
        ),
        box(
          width = NULL, title = "Fetch from NCBI GEO", status = "primary", solidHeader = FALSE,
          p(strong("GEO Series accession"), " (e.g. ", code("GSE12345"), "). Fetches a methylation beta-value matrix and sample metadata from NCBI. For your own data, use \"Upload your own data\" instead."),
          fluidRow(
            column(4, textInput(ns("geo_accession"), NULL, placeholder = "GSE12345", width = "100%")),
            column(3, actionButton(ns("geo_fetch_btn"), "Fetch", icon = icon("download"), class = "btn-primary btn-sm", width = "100%"))
          ),
          ## Same synchronous-fetch spinner workaround as the Transcriptomics
          ## Dataset tab's GEO card (see its own comment) - GEOquery::getGEO()
          ## blocks this app's single R process for the length of the fetch,
          ## so only a client-side script can show a spinner before that block
          ## starts.
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
          radioButtons(ns("upload_format"), "File format", inline = TRUE,
                       choices = c("Beta/M-value matrix (CSV/TSV/RDS)" = "matrix", "Raw IDAT files (.idat)" = "idat"),
                       selected = "matrix"),
          selectInput(ns("array_type"), "Dataset type", choices = METHYL_ARRAY_TYPES, selected = "EPIC", width = "100%"),
          uiOutput(ns("array_type_note")),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'matrix'", ns("upload_format")),
            p(strong("Methylation matrix"), " - CSV/TSV/RDS. Probes (CpGs) in rows, samples in columns; for CSV/TSV, the first column is the probe ID."),
            radioButtons(ns("input_scale"), "Input scale", inline = TRUE,
                         choices = c("Beta values (0-1)" = "beta", "M-values" = "m"), selected = "beta"),
            fileInput(ns("matrix_file"), "Methylation matrix", accept = c(".csv", ".tsv", ".txt", ".rds", ".Rds"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'idat'", ns("upload_format")),
            p(strong("IDAT files"), ": select every ", code("_Grn.idat"), " and ", code("_Red.idat"), " file at once. Beta values, detection p-values, and bead counts are all derived from these, unnormalized."),
            fileInput(ns("idat_files"), "IDAT files", multiple = TRUE, accept = c(".idat", ".idat.gz"))
          ),
          p(strong("Sample sheet / phenotype metadata"), " - CSV/TSV/RDS, one row per sample (optional, but needed to map group/sex/batch below)."),
          fileInput(ns("sheet_file"), "Sample sheet", accept = c(".csv", ".tsv", ".txt", ".rds", ".Rds")),
          uiOutput(ns("upload_preview_ui")),
          uiOutput(ns("upload_preview_tables_ui")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 2 · Map the columns"),
          uiOutput(ns("column_mapping")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 3 · Confirm"),
          actionButton(ns("load_btn"), "Upload Data", icon = icon("upload"), class = "btn-primary btn-sm"),
          div(class = "empty-note", style = "margin-top:6px;", icon("circle-info"),
              "Loads exactly what you provide, as-is - no normalising or batch correction. For that, use this module's own Normalization and Batch Correction sub-modules instead.")
        )
      )
    ),
    uiOutput(ns("load_message"))
  )
}

mod_methyl_dataset_server <- function(id, methyl_dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Guesses a mapping dropdown's default column by name, not position - same
    ## term-tier logic as the Transcriptomics Dataset tab's own guess_col()
    ## (see its comment for why exact-match and contains-match must each be
    ## tried term-by-term, not as one combined regex).
    guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
      for (term in exact) {
        hit <- cols[tolower(cols) == tolower(term)]
        if (length(hit) > 0) return(hit[1])
      }
      for (term in contains) {
        hit <- cols[grepl(term, cols, ignore.case = TRUE)]
        if (length(hit) > 0) return(hit[1])
      }
      fallback
    }

    output$array_type_note <- renderUI({
      req(input$array_type)
      if (input$array_type %in% c("450K", "EPIC")) {
        p(class = "empty-note", icon("circle-check"), "Manifest annotation available for this array type.")
      } else {
        p(class = "empty-note", icon("circle-info"), "No manifest annotation for this array type - sample-level QC still works.")
      }
    })

    ## ---- Preloaded path: loads the ~2.1GB matrix in the background and caches it per session.

    ## Tracks which preloaded choice has actually finished loading (via "Load
    ## this dataset"), as opposed to merely being picked in the dropdown - the
    ## GSE42861 info-card/note below must only appear once loading completes,
    ## not the instant a choice is selected.
    preloaded_loaded_choice <- reactiveVal(NULL, label = "preloaded_loaded_choice")

    finish_preloaded_load <- function(live) {
      pheno <- load_default_meth_pheno()
      validate(need(!is.null(pheno), "Could not read the preloaded dataset's sample metadata."))

      methyl_dataset$beta <- if (!is.null(live)) live$beta else NULL
      methyl_dataset$input_scale <- "beta"
      methyl_dataset$array_type <- "450K"
      methyl_dataset$sample_sheet <- if (!is.null(live)) live$pheno else pheno
      methyl_dataset$rg_set <- NULL
      methyl_dataset$mset <- NULL
      methyl_dataset$detp <- NULL
      methyl_dataset$beadcount <- NULL
      methyl_dataset$preloaded <- TRUE
      methyl_dataset$source_type <- "preloaded"
      methyl_dataset$source <- "Preloaded whole-blood dataset: sex-stratified RA cohort (GSE42861)"
      preloaded_loaded_choice(input$preloaded_choice)

      msg <- if (!is.null(live)) {
        tagList(
          sprintf("Loaded the preloaded whole-blood dataset - %s probes x %s samples, live.",
                  format(nrow(live$beta), big.mark = ","), ncol(live$beta)),
          br(),
          "Raw IDAT-only metrics (detection p-value, bead count, bisulfite conversion, sex-check) stay unavailable, since only the derived beta matrix is bundled."
        )
      } else {
        sprintf("Loaded the preloaded whole-blood dataset's metadata (%d samples). The live beta matrix isn't available in this deployment, so Quality Control, Normalization, Differential Methylation, and Differentially Methylated Regions show their default, already-completed pipeline analysis rather than a recomputable live tool.", nrow(pheno))
      }
      output$preloaded_load_message <- renderUI(
        span(style = "color: var(--color-success); font-size: 13px; font-weight: 600;", icon("check"), " ", msg)
      )
      showNotification("Preloaded whole-blood dataset loaded.", type = "message", duration = 5)
    }

    output$preloaded_note <- renderUI({
      req(input$preloaded_choice)
      req(identical(preloaded_loaded_choice(), input$preloaded_choice))
      p(class = "empty-note", icon("circle-info"),
        "689-sample whole-blood cohort (Liu et al. 2013), sex-stratified in every downstream sub-module's default analysis. The live beta matrix (2.1GB) loads in the background so the rest of the app stays usable.")
    })

    ## GEO info-card, matching the one the Transcriptomics Dataset tab shows for
    ## its own individual raw sources - static, verified facts (no live
    ## ExpressionSet is cached for this preloaded source, unlike transcriptomics'
    ## get_raw_eset()), except sample count, which reads the real bundled pheno
    ## table rather than being hardcoded. Only shown once "Load this dataset"
    ## has actually finished for the currently-selected choice, not merely once
    ## it's picked in the dropdown.
    output$preloaded_geo_card_ui <- renderUI({
      req(input$preloaded_choice)
      req(identical(preloaded_loaded_choice(), input$preloaded_choice))
      entry <- Find(function(d) d$id == input$preloaded_choice, MX_PRELOADED_DATASETS)
      req(entry)
      pheno <- load_default_meth_pheno()
      div(
        class = "info-card",
        div(
          class = "module-card-title-row",
          h4(entry$gse),
          tags$a(href = geo_link(entry$gse), target = "_blank", rel = "noopener",
                  icon("up-right-from-square"), " NCBI GEO")
        ),
        p(class = "module-card-tagline", "Differential DNA methylation in Rheumatoid arthritis"),
        p(strong("Role: "), "Internal training/discovery cohort", br(),
          strong("Used for: "), "Quality Control, Normalization, Differential Methylation, DMR calling, WGCNA, Diagnostic Classifier"),
        p(strong("Platform: "), "GPL13534 (Illumina HumanMethylation450 BeadChip)", br(),
          strong("Samples: "), if (!is.null(pheno)) nrow(pheno) else "unknown")
      )
    })

    preloaded_matrix_cache <- reactiveVal(NULL, label = "preloaded_matrix_cache")
    preloaded_matrix_fetched <- reactiveVal(FALSE, label = "preloaded_matrix_fetched")

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      preloaded_task <- ExtendedTask$new(function() {
        promises::future_promise(load_default_meth_matrix())
      })

      observeEvent(input$load_preloaded_btn, {
        req(METH_DATA_AVAILABLE)
        if (isTRUE(preloaded_matrix_fetched()) || !METH_RAW_DATA_AVAILABLE) {
          finish_preloaded_load(preloaded_matrix_cache())
        } else if (!is.null(.arthomix_cache[["meth_default_matrix"]])) {
          live <- load_default_meth_matrix()
          preloaded_matrix_cache(live)
          preloaded_matrix_fetched(TRUE)
          finish_preloaded_load(live)
        } else {
          output$preloaded_load_message <- renderUI(span(
            style = "font-size: 13px;", icon("spinner", class = "fa-spin"),
            " Loading the preloaded dataset's ~2.1GB live methylation matrix in the background - the rest of the app stays usable while this runs."))
          preloaded_task$invoke()
        }
      })

      observe({
        live <- tryCatch(preloaded_task$result(), error = function(e) e)
        if (inherits(live, "shiny.silent.error")) return()
        if (inherits(live, "error")) {
          output$preloaded_load_message <- renderUI(span(
            style = "color: var(--color-danger); font-size: 13px;", icon("triangle-exclamation"),
            paste(" Could not load the preloaded dataset's live matrix:", conditionMessage(live))))
          return()
        }
        .arthomix_cache[["meth_default_matrix"]] <- live
        preloaded_matrix_cache(live)
        preloaded_matrix_fetched(TRUE)
        finish_preloaded_load(live)
      })
    } else {
      observeEvent(input$load_preloaded_btn, {
        req(METH_DATA_AVAILABLE)
        live <- if (isTRUE(preloaded_matrix_fetched())) {
          preloaded_matrix_cache()
        } else {
          result <- withProgress(message = "Loading preloaded dataset", value = 0.2, {
            load_default_meth_matrix()
          })
          preloaded_matrix_cache(result)
          preloaded_matrix_fetched(TRUE)
          result
        }
        finish_preloaded_load(live)
      })
    }

    ## ---- Upload path: shared Step 1 (files) / Step 2 (column mapping) / Step 3 (confirm),
    ## branching on input$upload_format ("matrix" vs "idat") only where the two formats
    ## genuinely differ (parsing) - sample-sheet mapping and sample-ID matching are identical
    ## either way, matching the Transcriptomics Dataset tab's own upload pipeline.

    sheet_raw <- reactive({
      req(input$sheet_file)
      path <- input$sheet_file$datapath
      if (grepl("\\.rds$", input$sheet_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
        validate(need(is.data.frame(d), "The uploaded sample sheet RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        s <- methyl_parse_sample_sheet(path, input$sheet_file$name)
        validate(need(isTRUE(s$ok), s$error))
        s$df
      }
    })

    ## Cheap preview only - the CSV/TSV/RDS matrix parse itself is fast even for a
    ## full 450K/EPIC matrix, but IDAT reading (minfi::read.metharray.exp) is not,
    ## so that stays deferred to the "Upload Data" click below, exactly like the
    ## previous implementation's own comment on this ("Reading is deferred until
    ## \"Load dataset\" is clicked").
    matrix_preview <- reactive({
      req(input$upload_format == "matrix", input$matrix_file)
      path <- input$matrix_file$datapath
      if (grepl("\\.rds$", input$matrix_file$name, ignore.case = TRUE)) {
        m <- tryCatch(readRDS(path), error = function(e) NULL)
        if (is.null(m) || !is.matrix(m)) return(list(ok = FALSE, error = "The uploaded matrix RDS file must contain a numeric matrix."))
        list(ok = TRUE, mat = m)
      } else {
        methyl_parse_matrix(path, input$matrix_file$name)
      }
    })

    output$upload_preview_ui <- renderUI({
      if (identical(input$upload_format, "idat")) {
        req(input$idat_files)
        n_grn <- sum(grepl("_Grn\\.idat", input$idat_files$name, ignore.case = TRUE))
        n_red <- sum(grepl("_Red\\.idat", input$idat_files$name, ignore.case = TRUE))
        return(div(class = "empty-note", icon("circle-info"),
                   sprintf("%d file(s) selected (%d Grn, %d Red). Reading is deferred until \"Upload Data\" is clicked.", nrow(input$idat_files), n_grn, n_red)))
      }
      req(input$matrix_file)
      p <- matrix_preview()
      if (!isTRUE(p$ok)) {
        return(div(class = "empty-note", icon("triangle-exclamation"), p$error))
      }
      orient <- methyl_detect_orientation(p$mat)
      orient_note <- if (isTRUE(orient$transposed)) " This looks like it's oriented samples x probes - it will be transposed automatically on upload." else ""
      sheet_note <- if (!is.null(input$sheet_file)) {
        s <- tryCatch(list(ok = TRUE, df = sheet_raw()), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
        if (isTRUE(s$ok)) sprintf(" Sample sheet: %d row(s).", nrow(s$df)) else paste(" Sample sheet could not be read:", s$error)
      } else ""
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s probes x %s samples.%s%s", input$matrix_file$name,
                  format(nrow(p$mat), big.mark = ","), ncol(p$mat), sheet_note, orient_note))
    })

    output$upload_preview_tables_ui <- renderUI({
      req(!is.null(input$sheet_file))
      s <- tryCatch(sheet_raw(), error = function(e) NULL)
      req(!is.null(s))
      tagList(
        div(class = "upload-step-label", "Preview of what you uploaded"),
        p(class = "submodule-desc", "First 5 rows of the sample sheet, exactly as read - before any column mapping below."),
        DT::dataTableOutput(ns("upload_preview_meta_table"))
      )
    })

    output$upload_preview_meta_table <- DT::renderDataTable({
      s <- sheet_raw()
      DT::datatable(head(s, 5), rownames = FALSE,
                     options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$column_mapping <- renderUI({
      req(input$sheet_file)
      cols <- colnames(sheet_raw())
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
                    selected = guess_col(cols, c("batch", "chip", "plate", "slide", "sentrix", "cohort", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })

    ## Keeps "Upload Data" disabled until files are chosen and, if a sample sheet
    ## was provided, its required column mappings are set - a sheet is optional
    ## for methylomics (unlike the Transcriptomics Dataset tab, which requires
    ## metadata), so mapping is only required when there's a sheet to map.
    observe({
      files_ready <- if (identical(input$upload_format, "idat")) !is.null(input$idat_files) else !is.null(input$matrix_file)
      mapping_ready <- is.null(input$sheet_file) || (!is.null(input$map_id) && !is.null(input$map_group))
      if (isTRUE(files_ready) && isTRUE(mapping_ready)) shinyjs::enable("load_btn") else shinyjs::disable("load_btn")
    })

    ## Shared by both upload_format branches below: maps the sample sheet's
    ## chosen columns onto canonical sample/group/sex/batch columns (added
    ## alongside the original columns, never replacing them - every other
    ## sub-module's own column pickers, e.g. mod_methyl_qc.R's live_group_col,
    ## read the raw sheet as-is and would otherwise lose access to it), then
    ## intersects sample IDs against the matrix's column names via
    ## methyl_sheet_sample_ids() - the same ID-resolution convention already
    ## used throughout qc.R - so every downstream tool sees a consistently
    ## matched matrix + sheet.
    map_and_match_sheet <- function(mat) {
      if (is.null(input$sheet_file)) return(list(mat = mat, sheet = NULL))
      sheet <- sheet_raw()
      sheet$sample <- as.character(sheet[[input$map_id]])
      sheet$group  <- as.character(sheet[[input$map_group]])
      sheet$sex    <- if (!identical(input$map_sex, "(none)")) as.character(sheet[[input$map_sex]]) else NA_character_
      sheet$batch  <- if (!identical(input$map_batch, "(none)")) as.character(sheet[[input$map_batch]]) else NA_character_

      all_ids <- colnames(mat)
      sheet_ids <- methyl_sheet_sample_ids(sheet, all_ids)
      common <- intersect(all_ids, sheet_ids)
      validate(need(
        length(common) >= 4,
        "Fewer than 4 sample IDs in the methylation matrix match the sample sheet's sample-ID column. Check the column mapping."
      ))
      list(mat = mat[, common, drop = FALSE], sheet = sheet[match(common, sheet_ids), , drop = FALSE])
    }

    observeEvent(input$load_btn, {
      result <- tryCatch({
        if (identical(input$upload_format, "idat")) {
          validate(need(!is.null(input$idat_files), "No IDAT files selected."))
          res <- withProgress(message = "Reading IDAT files", value = 0.3, { methyl_read_idat(input$idat_files) })
          validate(need(isTRUE(res$ok), res$error))
          derived <- withProgress(message = "Deriving beta values, detection p-values, bead counts", value = 0.7, { methyl_idat_derive(res$rg) })
          validate(need(isTRUE(derived$ok), derived$reason))

          detected_array <- tryCatch({
            arr <- minfi::annotation(res$rg)[["array"]]
            if (grepl("450k", arr, ignore.case = TRUE)) "450K"
            else if (grepl("EPICv2", arr, ignore.case = TRUE)) "EPICv2"
            else if (grepl("EPIC", arr, ignore.case = TRUE)) "EPIC"
            else NA_character_
          }, error = function(e) NA_character_)
          array_type_used <- if (!is.na(detected_array)) detected_array else input$array_type
          mismatch_note <- if (!is.na(detected_array) && !identical(detected_array, input$array_type)) sprintf(
            " Note: the IDAT files themselves are %s, not the %s selected above - using %s (the file-detected type) for annotation instead.",
            detected_array, input$array_type, detected_array) else ""

          matched <- map_and_match_sheet(derived$beta)
          rg_matched <- tryCatch(res$rg[, colnames(matched$mat)], error = function(e) NULL)
          list(mat = matched$mat, sheet = matched$sheet, rg_set = rg_matched,
               mset = tryCatch(derived$mset[, colnames(matched$mat)], error = function(e) NULL),
               detp = if (!is.null(derived$detp)) derived$detp[, colnames(matched$mat), drop = FALSE] else NULL,
               beadcount = if (!is.null(derived$beadcount)) derived$beadcount[, colnames(matched$mat), drop = FALSE] else NULL,
               input_scale = "beta", array_type = array_type_used,
               source = sprintf("Uploaded IDAT files: %s (%s), raw/unnormalized beta values", array_type_used, ncol(matched$mat)),
               note = mismatch_note,
               n_probes_msg = sprintf("Loaded %s probes across %s samples from raw IDAT (unnormalized). Detection p-values%s and bead counts%s are available for Quality Control.",
                                       format(nrow(matched$mat), big.mark = ","), ncol(matched$mat),
                                       if (is.null(derived$detp)) " (not available)" else "",
                                       if (is.null(derived$beadcount)) " (not available)" else ""))
        } else {
          validate(need(!is.null(input$matrix_file), "No methylation matrix file selected."))
          p <- matrix_preview()
          validate(need(isTRUE(p$ok), p$error))
          checked <- methyl_validate_matrix_upload(p$mat, input$input_scale)
          validate(need(isTRUE(checked$ok), checked$error))

          matched <- map_and_match_sheet(checked$mat)
          list(mat = matched$mat, sheet = matched$sheet, rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
               input_scale = input$input_scale, array_type = input$array_type,
               source = sprintf("Uploaded %s matrix: %s (%s)",
                                 if (input$input_scale == "beta") "beta-value" else "M-value",
                                 input$matrix_file$name, input$array_type),
               note = checked$note %||% "",
               n_probes_msg = sprintf("Loaded %s probes across %s samples.", format(nrow(matched$mat), big.mark = ","), ncol(matched$mat)))
        }
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this dataset:", conditionMessage(result)))
        )
        showNotification(paste("Could not load this dataset:", conditionMessage(result)), type = "error", duration = NULL)
        return()
      }

      methyl_dataset$beta <- result$mat
      methyl_dataset$input_scale <- result$input_scale
      methyl_dataset$array_type <- result$array_type
      methyl_dataset$sample_sheet <- result$sheet
      methyl_dataset$rg_set <- result$rg_set
      methyl_dataset$mset <- result$mset
      methyl_dataset$detp <- result$detp
      methyl_dataset$beadcount <- result$beadcount
      methyl_dataset$preloaded <- FALSE
      methyl_dataset$source_type <- "upload"
      methyl_dataset$source <- result$source

      n_dup <- sum(duplicated(rownames(result$mat)))
      output$load_message <- renderUI(
        tagList(
          div(class = "empty-note", icon("check"), paste(result$n_probes_msg, result$note)),
          if (n_dup > 0) div(class = "empty-note", icon("triangle-exclamation"),
              sprintf("%d duplicated probe ID(s) were detected - only the first occurrence of each is kept.", n_dup))
        )
      )
      showNotification("Uploaded dataset loaded - Quality Control, Normalization, Differential Methylation and every other sub-module now run against it.", type = "message", duration = 5)
      ## A corner toast is easy to miss on a long page - a modal forces an
      ## explicit acknowledgement that the upload actually succeeded.
      showModal(modalDialog(
        title = tagList(icon("circle-check", style = "color: var(--color-success);"), " Data loaded"),
        p(result$n_probes_msg),
        if (nzchar(result$note %||% "")) p(class = "empty-note", result$note),
        if (n_dup > 0) p(class = "empty-note", icon("triangle-exclamation"),
            sprintf("%d duplicated probe ID(s) were detected - only the first occurrence of each is kept.", n_dup)),
        p("Quality Control, Normalization, Differential Methylation, and every other sub-module now run against this dataset."),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })

    ## ---- GEO fetch path ---------------------------------------------------

    geo_fetch_result <- eventReactive(input$geo_fetch_btn, {
      if (!requireNamespace("GEOquery", quietly = TRUE)) {
        return(simpleError("The GEOquery package is not installed in this deployment. Install it with BiocManager::install(\"GEOquery\") to enable fetching by GEO accession, or use \"Upload your own data\" instead."))
      }
      acc <- toupper(trimws(input$geo_accession %||% ""))
      if (!grepl("^GSE[0-9]+$", acc)) {
        return(simpleError("Enter a valid GEO Series accession, e.g. GSE12345."))
      }
      tryCatch({
        gset <- suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE))
        if (!is.list(gset) || length(gset) == 0) {
          stop(paste(acc, "returned no series matrix from GEO - check the accession is a Series (GSExxxxx), not a Sample (GSM) or Platform (GPL) ID."))
        }
        list(acc = acc, platforms = gset)
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

    ## Extracts beta/meta from the fetched series and confirms it's actually
    ## methylation data - either a recognized Illumina methylation GPL, or
    ## (for an unlisted/custom platform) values that mostly fall in the [0,1]
    ## beta-value range. Neither test passing means this is very likely a
    ## gene-expression or other non-methylation series, so it's rejected with
    ## a pointer to the Transcriptomics module instead.
    geo_expr_meta <- reactive({
      tryCatch({
        eset <- geo_eset()
        ex <- Biobase::exprs(eset)
        if (nrow(ex) == 0 || ncol(ex) == 0) {
          stop("This GEO series has no data matrix in its series matrix file - common for series that only deposit raw IDAT/supplementary files. Download those from the GEO page and use \"Upload your own data\" instead.")
        }
        gpl <- Biobase::annotation(eset)
        array_type <- unname(MX_METHYLATION_GPL[gpl])
        platform_note <- NULL
        if (is.na(array_type) || is.null(array_type)) {
          vals <- ex[is.finite(ex)]
          frac_in_unit <- if (length(vals) > 0) mean(vals >= -0.05 & vals <= 1.05) else 0
          if (frac_in_unit < 0.95) {
            stop(sprintf("Platform %s is not a recognized Illumina methylation array (450K/EPIC/EPICv2/27K), and its values don't look like methylation beta values (0-1 range) either - this looks like a gene-expression or other non-methylation dataset. Use the Transcriptomics module's GEO fetch instead.", gpl))
          }
          array_type <- "Custom array"
          platform_note <- sprintf("Platform %s was not recognized as a standard Illumina methylation array, but its values fall within the expected 0-1 beta-value range, so this was accepted as methylation data - verify the platform manually before relying on downstream results.", gpl)
        }
        list(expr = ex, meta = as.data.frame(Biobase::pData(eset)), array_type = array_type,
             platform = gpl, platform_note = platform_note)
      }, error = function(e) e)
    })

    output$geo_fetch_status <- renderUI({
      req(input$geo_fetch_btn)
      res <- geo_fetch_result()
      if (inherits(res, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), paste("Could not fetch from GEO:", conditionMessage(res))))
      }
      em <- geo_expr_meta()
      if (inherits(em, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), conditionMessage(em)))
      }
      tagList(
        div(class = "empty-note", icon("circle-info"),
            sprintf("Fetched %s (%s, %s): %s samples x %s probes. Beta values are assumed (the near-universal convention for GEO methylation series matrices). Map the columns below, then click \"Load this dataset\".",
                    res$acc, em$platform, em$array_type, ncol(em$expr), format(nrow(em$expr), big.mark = ","))),
        if (!is.null(em$platform_note)) div(class = "empty-note", icon("triangle-exclamation"), em$platform_note)
      )
    })

    ## No sample-ID mapping needed here - pData() rownames are already the same
    ## GSM accessions exprs() is indexed by (same as the Transcriptomics Dataset
    ## tab's own GEO path).
    output$geo_column_mapping <- renderUI({
      em <- geo_expr_meta()
      req(em); req(!inherits(em, "error"))
      cols <- colnames(em$meta)
      tagList(
        selectInput(ns("geo_map_group"), "Group / diagnosis column", choices = cols,
                    selected = guess_col(
                      cols,
                      exact = c("group", "diagnosis", "disease", "condition", "phenotype"),
                      contains = c("group", "diagnosis", "disease", "condition", "phenotype", "characteristics_ch1")
                    ),
                    selectize = FALSE),
        selectInput(ns("geo_map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("geo_map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "chip", "plate", "slide", "sentrix", "cohort", "dataset"), fallback = "(none)"),
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
          "Fewer than 4 sample IDs matched between the fetched methylation matrix and metadata."
        ))
        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        acc <- geo_fetch_result()$acc
        list(expr = expr, meta = meta, acc = acc, array_type = em$array_type, platform = em$platform)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this GEO dataset:", conditionMessage(result)))
        )
        showNotification(paste("Could not load this GEO dataset:", conditionMessage(result)), type = "error", duration = NULL)
      } else {
        methyl_dataset$beta <- result$expr
        methyl_dataset$input_scale <- "beta"
        methyl_dataset$array_type <- result$array_type
        methyl_dataset$sample_sheet <- result$meta
        methyl_dataset$rg_set <- NULL
        methyl_dataset$mset <- NULL
        methyl_dataset$detp <- NULL
        methyl_dataset$beadcount <- NULL
        methyl_dataset$preloaded <- FALSE
        methyl_dataset$source_type <- "geo"
        methyl_dataset$source <- sprintf("NCBI GEO: %s (%s, %s)", result$acc, result$platform, result$array_type)
        n_samples <- ncol(result$expr)
        output$load_message <- renderUI(
          div(class = "empty-note", icon("check"),
              sprintf("Loaded %s probes across %s samples. Every sub-module now runs on this GEO dataset.",
                      format(nrow(result$expr), big.mark = ","), n_samples))
        )
        showNotification("GEO dataset loaded - Quality Control, Normalization, Differential Methylation and every other sub-module now run against it.", type = "message", duration = 5)
      }
    })

  })
}
