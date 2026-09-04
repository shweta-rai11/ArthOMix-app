## R/crossomics/01_Data/mod_cross_dataset.R
## Cross-Omics "Dataset" tab: the module's single data-entry point for
## "Expression and Methylation". Two ways to arrive at the exact same shape -

mod_cross_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Load your own Transcriptomics/Methylomics analysis results from this session, or upload data in the same standardized format, for the Expression and Methylation sub-module to run on."
)

mod_cross_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "empty-note", icon("circle-info"),
        "Loads the data for Expression and Methylation. Pull in the Transcriptomics/Methylomics results you've already run in this session, or upload files below."),
    fluidRow(
      column(
        4,
        box(
          width = NULL, title = "1. Data Source", status = "primary", solidHeader = FALSE,
          radioButtons(ns("source_mode"), "Data source",
                       choices = c("My analysis results" = "example", "Upload your own data" = "upload"),
                       selected = "example"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'example'", ns("source_mode")),
            uiOutput(ns("live_source_ui"))
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

mod_cross_dataset_server <- function(id, cross_dataset, results = NULL, methyl_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    expr_data <- reactiveVal(NULL)
    meth_data <- reactiveVal(NULL)

    ## ---- "My analysis results" data source: live DGE/DMP runs via the same adapters as mod_cross_integration.R.
    live_dge_choices <- reactive({
      runs <- (results %||% list())$dge_runs %||% list()
      if (length(runs) == 0) return(NULL)
      stats::setNames(names(runs), vapply(runs, function(r) r$contrast %||% "(unnamed run)", character(1)))
    })

    live_dmp_run <- reactive({
      tbl <- (methyl_results %||% list())$dmp_table
      if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(NULL)
      list(comparison = ((methyl_results %||% list())$dmp %||% list())$comparison %||% "Live Methylomics DMP run", table = tbl)
    })

    output$live_source_ui <- renderUI({
      ch <- live_dge_choices()
      dmp <- live_dmp_run()
      if (is.null(ch) && is.null(dmp)) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    "Run Differential Expression in Transcriptomics and DMP Analysis in Methylomics first, then come back here to load your results."))
      }
      tagList(
        tags$div(style = "font-weight:600; margin-bottom:2px;", "Transcriptomics"),
        if (is.null(ch)) p(class = "empty-note", icon("triangle-exclamation"), "Run Differential Expression in Transcriptomics first.")
        else selectInput(ns("live_dge_run"), "DGE run", choices = ch, selected = unname(utils::tail(ch, 1))),
        tags$div(style = "font-weight:600; margin-bottom:2px;", "Methylomics"),
        if (is.null(dmp)) p(class = "empty-note", icon("triangle-exclamation"), "Run DMP Analysis in Methylomics first.")
        else p(class = "empty-note", icon("circle-info"), sprintf("Using \"%s\" (latest live DMP run this session).", dmp$comparison)),
        if (!is.null(ch) && !is.null(dmp))
          actionButton(ns("load_live_btn"), "Load my analysis results (Transcriptomics + Methylomics)", icon = icon("bolt"), class = "btn-primary btn-sm", width = "100%")
      )
    })

    observeEvent(input$load_live_btn, {
      req(input$live_dge_run)
      runs <- (results %||% list())$dge_runs %||% list()
      dge_run <- runs[[input$live_dge_run]]
      dmp_run <- live_dmp_run()

      expr_res <- cx_build_live_expr_df(dge_run)
      if (!expr_res$ok) { showNotification(expr_res$error, type = "error"); expr_data(NULL) }
      else expr_data(list(df = expr_res$df, source = sprintf("My analysis: Transcriptomics DGE run \"%s\"", dge_run$contrast %||% input$live_dge_run), raw = NULL, mapping = NULL))

      meth_res <- cx_build_live_meth_df(dmp_run)
      if (!meth_res$ok) { showNotification(meth_res$error, type = "error"); meth_data(NULL) }
      else meth_data(list(df = meth_res$df, source = sprintf("My analysis: Methylomics DMP run (%s)", dmp_run$comparison), raw = NULL, mapping = NULL))
    }, ignoreInit = TRUE)

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

    observeEvent(input$source_mode, { expr_data(NULL); meth_data(NULL) }, ignoreInit = TRUE)

    observeEvent(input$use_data_btn, {
      validate(need(!is.null(expr_data()) || !is.null(meth_data()), "Load your analysis results or upload a file first."))
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
    }, ignoreInit = TRUE)

    observeEvent(input$clear_btn, {
      expr_data(NULL); meth_data(NULL)
      cross_dataset$user_expr_df <- NULL; cross_dataset$user_expr_source <- NULL
      cross_dataset$user_expr_wide <- NULL; cross_dataset$user_expr_mapping <- NULL; cross_dataset$user_expr_sample_cols <- character(0)
      cross_dataset$user_meth_df <- NULL; cross_dataset$user_meth_source <- NULL
      cross_dataset$user_meth_wide <- NULL; cross_dataset$user_meth_mapping <- NULL; cross_dataset$user_meth_sample_cols <- character(0)
    }, ignoreInit = TRUE)

    output$preview_ui <- renderUI({
      if (is.null(expr_data()) && is.null(meth_data())) {
        msg <- if (identical(input$source_mode, "upload")) "Upload a Transcriptomics and/or Methylomics file to preview it here."
               else "Click \"Load my analysis results\" to preview both the Transcriptomics (DEG) and Methylomics (DMP) tables here."
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
