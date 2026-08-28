## R/mod_dataset.R
## Dataset tab: preview a preloaded dataset, upload, or GEO fetch - written only to
## dataset$staged_* (a candidate); Preprocessing promotes a staged pick to the active
## dataset$expr/meta/source that the rest of the app runs on.

mod_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
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

mod_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        6,
        box(
          width = NULL, title = "Switch to preloaded data", status = "primary", solidHeader = FALSE,
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
              "Loads exactly what you provide, as-is - no merging, normalising, or batch correction. For that, use Preprocessing and Batch Correction instead.")
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
        p(class = "empty-note", icon("circle-info"),
          "The same merged, batch-corrected training cohort the app loads by default on startup - pick this to switch back to it after loading something else.")
      } else if (input$preloaded_choice %in% c("GSE93272", "GSE110169")) {
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
        d <- readRDS(path)
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
        readRDS(path)
      } else {
        m <- as.data.frame(data.table::fread(path, showProgress = FALSE))
        rn <- as.character(m[[1]])
        m <- as.matrix(m[, -1, drop = FALSE])
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

    ## Guesses a mapping dropdown's default column by name, not position: exact match
    ## first, then substring, then falls back to the first column.
    guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
      hit <- cols[tolower(cols) %in% tolower(exact)]
      if (length(hit) > 0) return(hit[1])
      hit <- cols[grepl(paste(contains, collapse = "|"), cols, ignore.case = TRUE)]
      if (length(hit) > 0) return(hit[1])
      fallback
    }

    output$column_mapping <- renderUI({
      req(input$meta_file)
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

    ## Keeps "Upload Data" disabled until files and required column mappings are set.
    observe({
      ready <- !is.null(input$expr_file) && !is.null(input$meta_file) &&
        !is.null(input$map_id) && !is.null(input$map_group)
      if (isTRUE(ready)) shinyjs::enable("load_btn") else shinyjs::disable("load_btn")
    })

    observeEvent(input$load_preloaded_btn, {
      req(input$preloaded_choice)
      entry <- Find(function(d) d$id == input$preloaded_choice, PRELOADED_DATASETS)
      req(entry)
      d <- entry$load()
      dataset$staged_expr <- d$expr
      dataset$staged_meta <- d$meta
      dataset$staged_source <- d$source
      output$preloaded_load_message <- renderUI(
        span(style = "color: var(--color-success); font-size: 13px; font-weight: 600;", icon("check"), " ",
             sprintf("Previewed %s: %s genes across %s samples. This doesn't change what any sub-module runs on - go to Preprocessing and pick \"Currently loaded dataset\" to analyze it.",
                      entry$label, format(nrow(d$expr), big.mark = ","), ncol(d$expr)))
      )
    })

    observeEvent(input$load_btn, {
      req(input$expr_file, input$meta_file, input$map_id, input$map_group)

      result <- tryCatch({
        expr <- expr_raw()
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
        list(expr = expr, meta = meta)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this dataset:", conditionMessage(result)))
        )
      } else {
        dataset$staged_expr <- result$expr
        dataset$staged_meta <- result$meta
        dataset$staged_source <- paste0("Uploaded dataset: ", input$expr_file$name, " + ", input$meta_file$name)
        n_dup <- sum(duplicated(rownames(result$expr)))
        output$load_message <- renderUI(
          tagList(
            div(class = "empty-note", icon("check"),
                sprintf("Previewed %s genes across %s samples. This doesn't change what any sub-module runs on - go to Preprocessing and pick \"Currently loaded dataset\" to analyze it.",
                        format(nrow(result$expr), big.mark = ","), ncol(result$expr))),
            if (n_dup > 0) div(class = "empty-note", icon("triangle-exclamation"),
                sprintf("%d duplicated feature identifier(s) were detected in this dataset. All rows are kept here, but downstream row-name-keyed steps (e.g. the Preprocessing merge tab) will keep only the first occurrence of each - rename duplicates in your source file if this is unintended.", n_dup))
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

    ## Extracts expr/meta from the fetched ExpressionSet, collapsing probes to genes
    ## via collapse_probes_to_genes() when the platform annotation allows it.
    geo_expr_meta <- reactive({
      tryCatch({
        eset <- geo_eset()
        ex <- Biobase::exprs(eset)
        if (nrow(ex) == 0 || ncol(ex) == 0) {
          stop("This GEO series has no expression matrix in its series matrix file - common for RNA-seq series that only deposit raw counts as supplementary files. Download that file from the GEO page and use \"Upload your own data\" instead.")
        }
        collapsed <- tryCatch(collapse_probes_to_genes(eset), error = function(e) NULL)
        used_collapse <- !is.null(collapsed) && nrow(collapsed) > 0 && nrow(collapsed) < nrow(ex)
        expr <- if (used_collapse) collapsed else ex
        list(expr = expr, meta = as.data.frame(Biobase::pData(eset)),
             collapsed = used_collapse, platform = Biobase::annotation(eset))
      }, error = function(e) e)
    })

    ## GEO fetch info-card plus metadata preview, matching the Overview and Datasets
    ## tab's card format for the app's fixed GEO sources.
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
      eset <- geo_eset()
      tagList(
        div(
          class = "info-card",
          div(
            class = "module-card-title-row",
            h4(res$acc),
            tags$a(href = geo_link(res$acc), target = "_blank", rel = "noopener",
                    icon("up-right-from-square"), " NCBI GEO")
          ),
          p(class = "module-card-tagline",
            tryCatch(Biobase::experimentData(eset)@title, error = function(e) NULL)),
          p(strong("Platform: "), em$platform, br(),
            strong("Samples: "), ncol(em$expr), ", ",
            strong(if (em$collapsed) "Genes: " else "Probes/features: "), format(nrow(em$expr), big.mark = ",")),
          if (!em$collapsed) {
            div(class = "empty-note", icon("triangle-exclamation"),
                "No gene-symbol annotation found for this platform - left at probe/feature-ID level. You can still load it as-is, or use Preprocessing's own probe-collapse step afterward.")
          }
        ),
        strong("Sample metadata preview (first 5 rows)"),
        DT::dataTableOutput(ns("geo_preview_meta_table"))
      )
    })

    output$geo_preview_meta_table <- DT::renderDataTable({
      em <- geo_expr_meta()
      req(!inherits(em, "error"))
      DT::datatable(head(em$meta, 5), rownames = FALSE,
                     options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    ## No sample-ID mapping needed here - pData() rownames are already the same
    ## GSM accessions exprs() is indexed by.
    output$geo_column_mapping <- renderUI({
      em <- geo_expr_meta()
      req(em); req(!inherits(em, "error"))
      cols <- colnames(em$meta)
      tagList(
        selectInput(ns("geo_map_group"), "Group / diagnosis column", choices = cols,
                    selected = guess_col(cols, c("group", "diagnosis", "disease", "condition", "status", "phenotype", "characteristics_ch1")),
                    selectize = FALSE),
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
        label <- sprintf("%s (%s, %s)", geo_fetch_result()$acc, em$platform,
                          if (em$collapsed) "collapsed to genes" else "probe-level, raw")
        list(expr = expr, meta = meta, label = label)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this GEO dataset:", conditionMessage(result)))
        )
      } else {
        dataset$staged_expr <- result$expr
        dataset$staged_meta <- result$meta
        dataset$staged_source <- paste0("NCBI GEO: ", result$label)
        output$load_message <- renderUI(
          div(class = "empty-note", icon("check"),
              sprintf("Previewed %s genes across %s samples. This doesn't change what any sub-module runs on - go to Preprocessing and pick \"Currently loaded dataset\" to analyze it.",
                      format(nrow(result$expr), big.mark = ","), ncol(result$expr)))
        )
      }
    })

  })
}
