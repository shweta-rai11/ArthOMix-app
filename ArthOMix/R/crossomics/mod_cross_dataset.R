## R/crossomics/mod_cross_dataset.R
## Cross-Omics "Dataset" tab: the module's single data-entry point for
## "Expression and Methylation". Two ways to arrive at the exact same shape -
## a standardized Transcriptomics table (gene, log2fc, pvalue, fdr) and/or
## Methylomics table (cpg, gene, dbeta, pvalue, fdr) - through one shared
## preview and one shared hand-off:
##
## - "Example data" - this app's own real, sex-stratified Transcriptomics
##   (DEG) and Methylomics (DMP) results (cx_load_default_deg()/
##   cx_load_default_methylation(), crossomics_integration_helpers.R) - the
##   SAME files "Expression and Methylation"'s own "Preloaded data" mode uses,
##   standardized to the identical shape "Upload your own data" produces.
##   Shown through the identical preview, so it doubles as a worked example
##   of what an uploaded file should look like - not a different kind of
##   object from what you'd upload.
## - "Upload your own data" - your own files, auto-detected and standardized
##   with cx_read_and_detect()/cx_standardize_expression()/
##   cx_standardize_methylation() (crossomics_integration_upload.R,
##   crossomics_integration_helpers.R) - the exact same standardized shape
##   the example path's own loaders already produce.
##
## Either way, "Use this data" publishes the same standardized
## user_expr_df/user_meth_df pair into the shared `cross_dataset` store,
## which "Expression and Methylation"'s "From Dataset tab" input mode reads
## directly - loading the example or uploading your own is one real,
## working path, not two incompatible ones.
##
## Earlier version of this tab browsed the pipeline's own already-joined
## eQTL-MR x mQTL-MR x DEG x DMP x DMR biomarker-convergence tables
## (CX_TABLE_REGISTRY, still defined in data_paths.R and still exercised by
## tests/test-data-loaders.R) - removed from here because that shape has
## nothing to do with what Expression and Methylation needs, so it could never
## usefully serve as a worked example for the Upload path. Those tables are
## unrelated to this workflow: Biomarker Convergence and Cross-Omics MR
## below still load them directly and independently via their own Load
## buttons, unaffected by this change.
##
## If a file's required columns can't be auto-detected, this tab says so
## explicitly and points to Expression and Methylation's own Upload option,
## which supports manual column mapping; it never guesses at an ambiguous
## column.

mod_cross_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Load example Transcriptomics/Methylomics data, or upload your own in the same standardized format, for the Expression and Methylation sub-module to run on."
)

mod_cross_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "empty-note", icon("circle-info"),
        "Loads the data for Expression and Methylation, separate from the Transcriptomics and Methylomics tabs."),
    fluidRow(
      column(
        4,
        box(
          width = NULL, title = "1. Data Source", status = "primary", solidHeader = FALSE,
          radioButtons(ns("source_mode"), "Data source",
                       choices = c("Example data" = "example", "Upload your own data" = "upload"),
                       selected = "example"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'example'", ns("source_mode")),
            radioButtons(ns("sex_stratum"), "Analysis group",
                         choices = c("ALL" = "all", "FEMALE" = "female", "MALE" = "male"),
                         selected = "female", inline = TRUE),
            radioButtons(ns("meth_level"), "Methylation data",
                         choices = c("CpG-level (DMP)" = "dmp", "Region-level (DMR)" = "dmr"),
                         selected = "dmp", inline = TRUE),
            p(class = "submodule-desc", "Female/Male Transcriptomics (DEG) and Methylomics (DMP or DMR) example data, in the same format as \"Upload your own data.\""),
            actionButton(ns("load_example_btn"), "Load example data", icon = icon("database"), class = "btn-primary btn-sm")
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload'", ns("source_mode")),
            p(class = "submodule-desc", "Columns are auto-detected (gene symbol, log2FC/Δβ, P-value, FDR)."),
            fileInput(ns("expr_file"), "Transcriptomics file (optional)", accept = c(".csv", ".tsv", ".txt", ".xlsx"),
                      placeholder = "CSV / TSV / TXT / XLSX"),
            p(class = "empty-note", icon("circle-info"), "Differentially Expressed Genes (DEG) format - one row per gene, with a gene symbol/ID and a log2 fold-change column."),
            fileInput(ns("meth_file"), "Methylomics file (optional)", accept = c(".csv", ".tsv", ".txt", ".xlsx"),
                      placeholder = "CSV / TSV / TXT / XLSX"),
            p(class = "empty-note", icon("circle-info"), "Differentially Methylated Position/Region (DMP/DMR) format - one row per CpG or region, with a gene symbol/ID and a Δβ (methylation change) column.")
          ),
          tags$hr(),
          fluidRow(
            column(6, actionButton(ns("use_data_btn"), "Use this data", icon = icon("check"), class = "btn-primary btn-sm", width = "100%")),
            column(6, actionButton(ns("clear_btn"), "Clear", icon = icon("xmark"), class = "btn-sm", width = "100%"))
          )
        )
      ),
      column(
        8,
        box(width = NULL, title = "Preview", status = "primary", solidHeader = FALSE,
            uiOutput(ns("preview_ui")))
      )
    )
  )
}

mod_cross_dataset_server <- function(id, cross_dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Single shared representation regardless of source - list(df =
    ## standardized, source = display label, raw = original wide df or NULL,
    ## mapping = column mapping or NULL). `raw`/`mapping` stay NULL for
    ## example data (it's already gene/CpG-level, no per-sample columns to
    ## detect), matching Expression and Methylation's own Preloaded mode.
    expr_data <- reactiveVal(NULL)
    meth_data <- reactiveVal(NULL)

    ## ---- Example data -------------------------------------------------------

    observeEvent(input$load_example_btn, {
      sex <- input$sex_stratum
      deg <- cx_load_default_deg(sex = sex)
      if (is.null(deg)) {
        showNotification("Could not read the example Transcriptomics (DEG) table for this sex stratum.", type = "error")
        expr_data(NULL)
      } else {
        std <- cx_standardize_expression(deg, mapping = c(gene = "gene", log2fc = "logFC", pvalue = "P.Value", fdr = "adj.P.Val"))
        if (!std$ok) { showNotification(std$error, type = "error"); expr_data(NULL) }
        else expr_data(list(df = std$df, source = sprintf("Example data (%s, %s DEG)", toupper(sex), if (identical(sex, "all")) "pooled" else "sex-stratified"), raw = NULL, mapping = NULL))
      }

      meth <- if (identical(input$meth_level, "dmr")) cx_load_default_dmr(sex = sex) else cx_load_default_methylation(sex = sex)
      if (!meth$ok) {
        showNotification(meth$error, type = "warning", duration = 10)
        meth_data(NULL)
      } else {
        strat_word <- if (identical(sex, "all")) "pooled" else "sex-stratified"
        meth_label <- if (identical(input$meth_level, "dmr")) sprintf("%s DMR", strat_word) else sprintf("%s DMP, SVA/bacon-adjusted", strat_word)
        meth_data(list(df = meth$df, source = sprintf("Example data (%s, %s)", toupper(sex), meth_label), raw = NULL, mapping = NULL))
      }
    })

    ## ---- Upload your own data ---------------------------------------------
    ## Reuses the exact same read/auto-detect/standardize helpers Expression
    ## x Methylation's own Upload mode uses (crossomics_integration_upload.R,
    ## crossomics_integration_helpers.R) - no new parsing/detection logic.

    observeEvent(input$expr_file, {
      res <- cx_read_and_detect(input$expr_file$datapath, input$expr_file$name, kind = "expression")
      if (!res$ok) { showNotification(res$error, type = "error"); expr_data(NULL); return() }
      std <- cx_standardize_expression(res$df, res$mapping)
      if (!std$ok) {
        showNotification(
          sprintf("Transcriptomics file: %s Required columns could not be auto-detected here - use Expression and Methylation's own Upload option instead, which supports manual column mapping.", std$error),
          type = "warning", duration = 15
        )
        expr_data(NULL)
        return()
      }
      expr_data(list(df = std$df, source = sprintf("Uploaded: %s", input$expr_file$name), raw = res$df, mapping = res$mapping))
    })

    observeEvent(input$meth_file, {
      res <- cx_read_and_detect(input$meth_file$datapath, input$meth_file$name, kind = "methylation")
      if (!res$ok) { showNotification(res$error, type = "error"); meth_data(NULL); return() }
      std <- cx_standardize_methylation(res$df, res$mapping)
      if (!std$ok) {
        showNotification(
          sprintf("Methylomics file: %s Required columns could not be auto-detected here - use Expression and Methylation's own Upload option instead, which supports manual column mapping.", std$error),
          type = "warning", duration = 15
        )
        meth_data(NULL)
        return()
      }
      meth_data(list(df = std$df, source = sprintf("Uploaded: %s", input$meth_file$name), raw = res$df, mapping = res$mapping))
    })

    ## Switching data source clears whatever the other mode had loaded, so
    ## stale example/uploaded data can't be silently carried into "Use this
    ## data" for the mode you're no longer looking at.
    observeEvent(input$source_mode, { expr_data(NULL); meth_data(NULL) }, ignoreInit = TRUE)

    ## ---- Hand-off into the shared cross_dataset store (shared by both modes) --

    observeEvent(input$use_data_btn, {
      validate(need(!is.null(expr_data()) || !is.null(meth_data()), "Load example data or upload a file first."))
      if (!is.null(expr_data())) {
        cross_dataset$user_expr_df <- expr_data()$df
        cross_dataset$user_expr_source <- expr_data()$source
        cross_dataset$user_expr_wide <- expr_data()$raw
        cross_dataset$user_expr_mapping <- expr_data()$mapping
        cross_dataset$user_expr_sample_cols <- if (!is.null(expr_data()$raw)) cx_detect_sample_columns(expr_data()$raw, expr_data()$mapping) else character(0)
      }
      if (!is.null(meth_data())) {
        cross_dataset$user_meth_df <- meth_data()$df
        cross_dataset$user_meth_source <- meth_data()$source
        cross_dataset$user_meth_wide <- meth_data()$raw
        cross_dataset$user_meth_mapping <- meth_data()$mapping
        cross_dataset$user_meth_sample_cols <- if (!is.null(meth_data()$raw)) cx_detect_sample_columns(meth_data()$raw, meth_data()$mapping) else character(0)
      }
      showNotification("Ready for Expression and Methylation.", type = "message")
    })

    observeEvent(input$clear_btn, {
      expr_data(NULL); meth_data(NULL)
      cross_dataset$user_expr_df <- NULL; cross_dataset$user_expr_source <- NULL
      cross_dataset$user_expr_wide <- NULL; cross_dataset$user_expr_mapping <- NULL; cross_dataset$user_expr_sample_cols <- character(0)
      cross_dataset$user_meth_df <- NULL; cross_dataset$user_meth_source <- NULL
      cross_dataset$user_meth_wide <- NULL; cross_dataset$user_meth_mapping <- NULL; cross_dataset$user_meth_sample_cols <- character(0)
    }, ignoreInit = TRUE)

    ## ---- Preview (identical for both modes) ----------------------------------

    output$preview_ui <- renderUI({
      if (is.null(expr_data()) && is.null(meth_data())) {
        msg <- if (identical(input$source_mode, "upload")) "Upload a Transcriptomics and/or Methylomics file to preview it here."
               else "Click \"Load example data\" to preview it here."
        return(div(class = "empty-note", icon("circle-info"), msg))
      }
      tagList(
        if (!is.null(expr_data())) tagList(
          h5("Transcriptomics"),
          p(class = "submodule-desc", sprintf("%s - %s genes standardized (gene, log2FC, P-value, FDR).",
                                                expr_data()$source, format(nrow(expr_data()$df), big.mark = ","))),
          fluidRow(
            column(7, DT::dataTableOutput(ns("expr_table"))),
            column(5, plotOutput(ns("expr_plot"), height = "260px"))
          ),
          tags$hr()
        ),
        if (!is.null(meth_data())) tagList(
          h5("Methylomics"),
          p(class = "submodule-desc", sprintf("%s - %s CpG/gene records standardized (gene, Δβ, P-value, FDR).",
                                                meth_data()$source, format(nrow(meth_data()$df), big.mark = ","))),
          fluidRow(
            column(7, DT::dataTableOutput(ns("meth_table"))),
            column(5, plotOutput(ns("meth_plot"), height = "260px"))
          )
        )
      )
    })

    output$expr_table <- DT::renderDataTable({
      req(expr_data())
      DT::datatable(expr_data()$df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5), class = "stripe hover compact")
    })
    output$expr_plot <- renderPlot({
      req(expr_data())
      df <- expr_data()$df
      validate(need(any(!is.na(df$log2fc)), "No numeric log2FC values to plot."))
      ggplot2::ggplot(df, ggplot2::aes(x = log2fc)) +
        ggplot2::geom_histogram(bins = 40, fill = ARTHOMIX_COLORS$blue, na.rm = TRUE) +
        ggplot2::labs(x = "log2 Fold Change", y = "Genes", title = "Distribution of log2FC") +
        theme_arthomix()
    })

    output$meth_table <- DT::renderDataTable({
      req(meth_data())
      DT::datatable(meth_data()$df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5), class = "stripe hover compact")
    })
    output$meth_plot <- renderPlot({
      req(meth_data())
      df <- meth_data()$df
      validate(need(any(!is.na(df$dbeta)), "No numeric Δβ values to plot."))
      ggplot2::ggplot(df, ggplot2::aes(x = dbeta)) +
        ggplot2::geom_histogram(bins = 40, fill = ARTHOMIX_COLORS$aqua, na.rm = TRUE) +
        ggplot2::labs(x = "Δβ (methylation change)", y = "CpG/gene records", title = "Distribution of Δβ") +
        theme_arthomix()
    })
  })
}
