## R/mod_methyl_dataset.R
## Methylomics "Dataset" tab: loads the shared `methyl_dataset` reactiveValues
## every Methylomics sub-module (starting with Quality Control, see
## mod_methyl_qc.R) reads from - the methylation-side equivalent of
## mod_dataset.R's `dataset`, but its own separate reactiveValues (never the
## transcriptomics `dataset`), since the two omics layers are independent
## datasets until the Cross-Omics module links them.
##
## Two upload paths, picked with a radio button rather than offered as two
## always-visible boxes (unlike mod_dataset.R's Dataset tab) because they
## produce genuinely different downstream state - a matrix upload has no
## RGChannelSet, so IDAT-only QC metrics (detection p-value, bead count,
## bisulfite conversion, median intensity, minfi::getSex()) are unavailable
## until this module is re-run with an IDAT upload instead.
##
## Per this app's "no guessing" methylomics requirement: manifest/array
## annotation is looked up from the installed Bioconductor annotation
## packages for the array type picked below (see annotation.R, sourced
## automatically alongside every other file in this folder), not from a
## user-uploaded manifest file - a raw manifest CSV is not something this
## module parses itself.

mod_methyl_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Upload a methylation dataset - either a beta/M-value matrix or raw IDAT files - for every Methylomics sub-module below to read from."
)

mod_methyl_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "empty-note", icon("circle-info"),
        "Independent of the Transcriptomics dataset - nothing loaded here changes what any Transcriptomics sub-module reads, and nothing you did there changes this."),
    fluidRow(
      column(
        4,
        box(
          width = NULL, title = "1. Data source", status = "primary", solidHeader = FALSE,
          radioButtons(ns("upload_mode"), NULL, width = "100%",
                       choices = c("Preloaded whole-blood dataset (sex-stratified RA)" = "preloaded",
                                   "Upload: beta or M-value matrix" = "matrix",
                                   "Upload: raw IDAT files" = "idat"),
                       selected = if (METH_DATA_AVAILABLE) "preloaded" else "matrix"),
          conditionalPanel(
            condition = sprintf("input['%s'] != 'preloaded'", ns("upload_mode")),
            selectInput(ns("array_type"), "Dataset type", choices = METHYL_ARRAY_TYPES, selected = "EPIC", width = "100%"),
            uiOutput(ns("array_type_note"))
          )
        )
      ),
      column(
        8,
        conditionalPanel(
          condition = sprintf("input['%s'] == 'preloaded'", ns("upload_mode")),
          box(
            width = NULL, title = "2. Preloaded whole-blood dataset", status = "primary", solidHeader = FALSE,
            if (METH_DATA_AVAILABLE) tagList(
              p(strong("Whole blood, 689 samples"), " (354 rheumatoid arthritis / 335 control; 492 female / 197 male), Illumina 450K array. Each sub-module below (QC, Normalization, Differential Methylation) opens to this cohort's own completed, sex-stratified analysis by default."),
              p(class = "empty-note", icon("circle-info"),
                "These results are reproduced from a completed pipeline run rather than recomputed live, since the underlying QC'd beta matrix (~2.1GB) isn't bundled with this deployment. Upload your own dataset for a fully live, rerunnable analysis."),
              actionButton(ns("load_preloaded_btn"), "Load preloaded dataset", icon = icon("database"), class = "btn-primary btn-sm")
            ) else p(class = "empty-note", icon("triangle-exclamation"),
                     "The preloaded whole-blood dataset isn't available in this deployment - upload your own dataset instead.")
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'matrix'", ns("upload_mode")),
          box(
            width = NULL, title = "2. Upload a beta/M-value matrix", status = "primary", solidHeader = FALSE,
            p(strong("Methylation matrix"), " — CSV/TSV. Probes (CpGs) in rows, samples in columns; first column the probe ID."),
            radioButtons(ns("input_scale"), "Input scale", inline = TRUE,
                         choices = c("Beta values (0-1)" = "beta", "M-values" = "m"), selected = "beta"),
            fileInput(ns("matrix_file"), "Methylation matrix", accept = c(".csv", ".tsv", ".txt")),
            p(strong("Sample sheet / phenotype metadata"), " — CSV/TSV, one row per sample (optional but needed for group-aware analyses later)."),
            fileInput(ns("sheet_file"), "Sample sheet / phenotype metadata (optional)", accept = c(".csv", ".tsv", ".txt")),
            uiOutput(ns("matrix_preview_ui")),
            actionButton(ns("load_matrix_btn"), "Load dataset", icon = icon("upload"), class = "btn-primary btn-sm")
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'idat'", ns("upload_mode")),
          box(
            width = NULL, title = "2. Upload raw IDAT files", status = "primary", solidHeader = FALSE,
            p(strong("IDAT files"), " — select every ", code("_Grn.idat"), " and ", code("_Red.idat"),
              " file for the samples to load, all at once. Enables detection p-value, bead count, bisulfite conversion, and raw-intensity sex-check QC that a matrix upload can't provide."),
            fileInput(ns("idat_files"), "IDAT files", multiple = TRUE, accept = c(".idat", ".idat.gz")),
            p(strong("Sample sheet"), " — CSV/TSV, one row per sample (optional; used for group-aware analyses later)."),
            fileInput(ns("idat_sheet_file"), "Sample sheet (optional)", accept = c(".csv", ".tsv", ".txt")),
            uiOutput(ns("idat_preview_ui")),
            actionButton(ns("load_idat_btn"), "Load dataset", icon = icon("upload"), class = "btn-primary btn-sm")
          )
        ),
        uiOutput(ns("load_message"))
      )
    )
  )
}

mod_methyl_dataset_server <- function(id, methyl_dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$array_type_note <- renderUI({
      req(input$array_type)
      if (input$array_type %in% c("450K", "EPIC")) {
        p(class = "empty-note", icon("circle-check"),
          "Manifest annotation (chromosome position, SNP overlap) is available for this array type.")
      } else {
        p(class = "empty-note", icon("triangle-exclamation"),
          "No Bioconductor manifest annotation is available for this array type in this deployment - SNP/sex-chromosome probe filters and the manifest-based sex check will be unavailable; missingness, variance, mean-range, and sample-level QC still work.")
      }
    })

    ## ---- Preloaded path ------------------------------------------------------
    ## The ~2.1GB live beta matrix read (METH_RAW_DATA_ROOT) is the one
    ## genuinely slow step here - previously a plain synchronous
    ## load_default_meth_matrix() call, which blocked this entire app (every
    ## tab, every open session) for however long readRDS() took. Pushed to a
    ## background process via shiny::ExtendedTask (see global.R's
    ## ARTHOMIX_ASYNC_AVAILABLE/future::plan()) so the app stays fully
    ## usable while it loads; falls back to the previous synchronous
    ## behavior if `future`/`promises` aren't installed. Cached in
    ## preloaded_matrix_cache() for the life of this session so a second
    ## "Load preloaded dataset" click doesn't re-read the file.

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
      methyl_dataset$source <- "Preloaded whole-blood dataset: sex-stratified RA cohort"

      msg <- if (!is.null(live)) {
        sprintf("Loaded the preloaded whole-blood dataset - %s probes x %s samples, live. Every filter/control on Quality Control, Normalization, Differential Methylation, and Differentially Methylated Regions now recomputes against this real matrix (raw IDAT-only metrics - detection p-value, bead count, bisulfite conversion, minfi::getSex() - stay unavailable, since only the derived beta matrix is bundled, not the raw intensities).",
                format(nrow(live$beta), big.mark = ","), ncol(live$beta))
      } else {
        sprintf("Loaded the preloaded whole-blood dataset's metadata (%d samples). The live beta matrix isn't available in this deployment, so Quality Control, Normalization, Differential Methylation, and Differentially Methylated Regions show their default, already-completed pipeline analysis rather than a recomputable live tool.", nrow(pheno))
      }
      output$load_message <- renderUI(div(class = "empty-note", icon("check"), msg))
    }

    preloaded_matrix_cache <- reactiveVal(NULL, label = "preloaded_matrix_cache")
    preloaded_matrix_fetched <- reactiveVal(FALSE, label = "preloaded_matrix_fetched")

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      preloaded_task <- ExtendedTask$new(function() {
        promises::future_promise(load_default_meth_matrix())
      })

      observeEvent(input$load_preloaded_btn, {
        req(METH_DATA_AVAILABLE)
        if (isTRUE(preloaded_matrix_fetched()) || !METH_RAW_DATA_AVAILABLE) {
          ## Already fetched once this session, or there's no big matrix to
          ## fetch at all in this deployment - both are already fast, so no
          ## need to go through the background task.
          finish_preloaded_load(preloaded_matrix_cache())
        } else {
          output$load_message <- renderUI(div(class = "empty-note", icon("spinner", class = "fa-spin"),
            "Loading the preloaded dataset's ~2.1GB live methylation matrix in the background - the rest of the app (including this tab) stays usable while this runs. This can take a minute or more; you'll see a confirmation here when it's done."))
          preloaded_task$invoke()
        }
      })

      ## A plain observe() (not observeEvent/isolate - the ExtendedTask docs
      ## are explicit that result()'s reactive invalidation is ignored by
      ## those) that reacts to the task's status. result() throws a special
      ## silent "not ready yet" condition while "initial"/"running" - caught
      ## and ignored below - and either the real value ("success") or the
      ## real error ("error").
      observe({
        live <- tryCatch(preloaded_task$result(), error = function(e) e)
        if (inherits(live, "shiny.silent.error")) return()
        if (inherits(live, "error")) {
          output$load_message <- renderUI(div(class = "empty-note", icon("triangle-exclamation"),
            paste("Could not load the preloaded dataset's live matrix:", conditionMessage(live))))
          return()
        }
        preloaded_matrix_cache(live)
        preloaded_matrix_fetched(TRUE)
        finish_preloaded_load(live)
      })
    } else {
      ## `future`/`promises` not installed - same behavior as before
      ## (synchronous, blocks the app while the file loads).
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

    ## ---- Matrix upload path ------------------------------------------------

    matrix_parsed <- reactive({
      req(input$matrix_file)
      methyl_parse_matrix(input$matrix_file$datapath, input$matrix_file$name)
    })
    sheet_parsed <- reactive({
      req(input$sheet_file)
      methyl_parse_sample_sheet(input$sheet_file$datapath, input$sheet_file$name)
    })

    output$matrix_preview_ui <- renderUI({
      req(input$matrix_file)
      p <- matrix_parsed()
      if (!isTRUE(p$ok)) {
        return(div(class = "empty-note", icon("triangle-exclamation"), p$error))
      }
      sheet_note <- if (!is.null(input$sheet_file)) {
        s <- sheet_parsed()
        if (isTRUE(s$ok)) sprintf(" Sample sheet: %d row(s).", nrow(s$df)) else paste(" Sample sheet could not be read:", s$error)
      } else ""
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s probes x %s samples.%s", input$matrix_file$name,
                  format(nrow(p$mat), big.mark = ","), ncol(p$mat), sheet_note))
    })

    observeEvent(input$load_matrix_btn, {
      p <- matrix_parsed()
      if (!isTRUE(p$ok)) {
        output$load_message <- renderUI(div(class = "empty-note", icon("triangle-exclamation"), p$error))
        return()
      }
      sheet_df <- NULL
      if (!is.null(input$sheet_file)) {
        s <- sheet_parsed()
        if (isTRUE(s$ok)) sheet_df <- s$df
      }
      methyl_dataset$beta <- p$mat
      methyl_dataset$input_scale <- input$input_scale
      methyl_dataset$array_type <- input$array_type
      methyl_dataset$sample_sheet <- sheet_df
      methyl_dataset$rg_set <- NULL
      methyl_dataset$mset <- NULL
      methyl_dataset$detp <- NULL
      methyl_dataset$beadcount <- NULL
      methyl_dataset$preloaded <- FALSE
      methyl_dataset$source <- sprintf("Uploaded %s matrix: %s (%s)",
                                        if (input$input_scale == "beta") "beta-value" else "M-value",
                                        input$matrix_file$name, input$array_type)
      output$load_message <- renderUI(
        div(class = "empty-note", icon("check"),
            sprintf("Loaded %s probes across %s samples.", format(nrow(p$mat), big.mark = ","), ncol(p$mat)))
      )
    })

    ## ---- IDAT upload path ---------------------------------------------------

    idat_parsed <- reactive({
      req(input$idat_files)
      methyl_read_idat(input$idat_files)
    })
    idat_sheet_parsed <- reactive({
      req(input$idat_sheet_file)
      methyl_parse_sample_sheet(input$idat_sheet_file$datapath, input$idat_sheet_file$name)
    })

    output$idat_preview_ui <- renderUI({
      req(input$idat_files)
      n_grn <- sum(grepl("_Grn\\.idat", input$idat_files$name, ignore.case = TRUE))
      n_red <- sum(grepl("_Red\\.idat", input$idat_files$name, ignore.case = TRUE))
      div(class = "empty-note", icon("circle-info"),
          sprintf("%d file(s) selected (%d Grn, %d Red). Reading is deferred until \"Load dataset\" is clicked.", nrow(input$idat_files), n_grn, n_red))
    })

    observeEvent(input$load_idat_btn, {
      res <- withProgress(message = "Reading IDAT files", value = 0.3, {
        methyl_read_idat(input$idat_files)
      })
      if (!isTRUE(res$ok)) {
        output$load_message <- renderUI(div(class = "empty-note", icon("triangle-exclamation"), res$error))
        return()
      }
      derived <- withProgress(message = "Deriving beta values, detection p-values, bead counts", value = 0.7, {
        methyl_idat_derive(res$rg)
      })
      if (!isTRUE(derived$ok)) {
        output$load_message <- renderUI(div(class = "empty-note", icon("triangle-exclamation"), derived$reason))
        return()
      }
      sheet_df <- NULL
      if (!is.null(input$idat_sheet_file)) {
        s <- idat_sheet_parsed()
        if (isTRUE(s$ok)) sheet_df <- s$df
      }
      ## The IDAT files themselves record their true array type (read off
      ## the physical chip barcode, not a guess) - prefer that over
      ## input$array_type, which for this upload path is just a dropdown
      ## the user may not have changed from its default. Using the wrong
      ## type here silently breaks every downstream manifest lookup
      ## (SNP/sex-chromosome filters, the sex check, BMIQ/PBC's probe-type
      ## vector) since probes that don't resolve against the wrong
      ## package's manifest are just dropped, not flagged.
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

      methyl_dataset$beta <- derived$beta
      methyl_dataset$input_scale <- "beta"
      methyl_dataset$array_type <- array_type_used
      methyl_dataset$sample_sheet <- sheet_df
      methyl_dataset$rg_set <- res$rg
      methyl_dataset$mset <- derived$mset
      methyl_dataset$detp <- derived$detp
      methyl_dataset$beadcount <- derived$beadcount
      methyl_dataset$preloaded <- FALSE
      methyl_dataset$source <- sprintf("Uploaded IDAT files: %s (%s), raw/unnormalized beta values", array_type_used, ncol(derived$beta))
      output$load_message <- renderUI(
        div(class = "empty-note", icon("check"),
            sprintf("Loaded %s probes across %s samples from raw IDAT (unnormalized). Detection p-values%s and bead counts%s are available for Quality Control.%s",
                    format(nrow(derived$beta), big.mark = ","), ncol(derived$beta),
                    if (is.null(derived$detp)) " (not available)" else "",
                    if (is.null(derived$beadcount)) " (not available)" else "",
                    mismatch_note))
      )
    })
  })
}
