## R/multiomics/mod_multi_concordance.R
## Submodule: Gene <-> CpG Concordance - the pipeline's own genomic-region-
## aware correlation analysis (every CpG annotated to a gene, region -
## TSS200/TSS1500/5'UTR/1stExon/Body/3'UTR/ExonBnd - and CpG-island context,
## then classified into the canonical promoter-hypomethylation-> higher-
## expression / gene-body-hypermethylation-> higher-expression patterns, or
## flagged non-canonical where the data doesn't fit that textbook direction
## - spec: "do not assume every methylation site should negatively correlate
## with expression"). Descriptive, not independent statistical evidence
## (panel genes were selected on the same patients this recomputes logFC/
## delta-M for) - the pipeline's own circularity caveat is shown verbatim.

mod_multi_concordance_config <- list(
  id = "concordance", title = "Gene <-> CpG Concordance", icon = "arrows-left-right", group = "Biomarker modeling",
  description = "Expression log2FC vs. methylation delta-M for every CpG mapped to each candidate panel gene, colored by genomic region and by the canonical vs. non-canonical direction of association."
)

MULTI_CONCORDANCE_COHORTS <- c(
  "Drug x sex (Etanercept panel)" = "Gene <-> CpG concordance — drug x sex (Etanercept panel)",
  "Response (drug-pooled)" = "Gene <-> CpG concordance — response (drug-pooled)"
)

mod_multi_concordance_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Cohort", status = "primary", solidHeader = FALSE,
        selectInput(ns("cohort"), NULL, choices = MULTI_CONCORDANCE_COHORTS, width = "100%"),
        actionButton(ns("load_btn"), "Load table", icon = icon("database"), class = "btn-primary btn-sm", width = "100%")
      ),
      box(
        width = NULL, title = "2. Filters", status = "primary", solidHeader = FALSE, collapsible = TRUE,
        selectizeInput(ns("region"), "Genomic region", choices = NULL, multiple = TRUE, options = list(placeholder = "All regions")),
        selectizeInput(ns("sex"), "Sex", choices = c("Female" = "female", "Male" = "male"), multiple = TRUE, options = list(placeholder = "Both sexes")),
        checkboxInput(ns("add_fdr"), "Recompute BH-FDR over this table's tested pairs (from its retained raw p-values)", value = FALSE)
      ),
      div(class = "empty-note", icon("circle-info"),
          "Descriptive, not independent statistical evidence: panel genes were selected on these same patients, so this does not re-test significance - it characterizes the direction and genomic context of methylation-expression association for genes already flagged as candidates."),
      div(class = "empty-note", icon("circle-info"),
          "One row per (gene, CpG) pair, not per CpG - a CpG overlapping more than one gene's annotated region appears once per gene it maps to (one-to-many mapping is shown explicitly, never silently collapsed).")
    ),
    column(
      8,
      uiOutput(ns("summary_ui")),
      tabsetPanel(
        id = ns("tabs"), type = "tabs",
        tabPanel("Scatter", br(), uiOutput(ns("scatter_ui"))),
        tabPanel("Pattern counts", br(), uiOutput(ns("pattern_ui"))),
        tabPanel("Table", br(), uiOutput(ns("table_ui")))
      )
    )
  )
}

mod_multi_concordance_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    raw <- reactiveValues(df = NULL, cohort = NULL)

    observeEvent(input$load_btn, {
      res <- multi_read_registry_table(input$cohort)
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      raw$df <- res$df
      raw$cohort <- names(MULTI_CONCORDANCE_COHORTS)[MULTI_CONCORDANCE_COHORTS == input$cohort]
      if ("region" %in% colnames(res$df)) updateSelectizeInput(session, "region", choices = sort(unique(res$df$region)), server = FALSE)
      showNotification(sprintf("Loaded %s gene-CpG pairs.", format(nrow(res$df), big.mark = ",")), type = "message")
    })

    filtered <- reactive({
      req(raw$df)
      df <- raw$df
      if (length(input$sex) > 0) df <- df[tolower(df$sex) %in% tolower(input$sex), , drop = FALSE]
      if (isTRUE(input$add_fdr)) df <- multi_concordance_add_fdr(df)
      df
    })

    output$summary_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state("Load a cohort to see gene<->CpG concordance here."))
      n_genes <- length(unique(df$SYMBOL %||% df$gene))
      n_cpg <- length(unique(df$CpG))
      n_canon <- sum(grepl("^concordant", df$biological_pattern))
      div(class = "empty-note", icon("chart-simple"),
          sprintf(" %s genes, %s CpGs, %d/%d (%.0f%%) pairs in the canonical direction.",
                  format(n_genes, big.mark = ","), format(n_cpg, big.mark = ","), n_canon, nrow(df), 100 * n_canon / nrow(df)))
    })

    output$scatter_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state())
      tagList(
        multi_plot_or_empty(scatter_plot_fn, ns("scatter_plot"), "No gene-CpG pairs match the current region/sex filters.", height = "440px"),
        div(class = "table-toolbar", downloadButton(ns("dl_scatter_png"), "Download plot (PNG)", class = "btn-sm"))
      )
    })
    scatter_plot_fn <- reactive(multi_concordance_scatter(filtered(), region_filter = input$region))
    output$scatter_plot <- renderPlot(scatter_plot_fn())
    output$dl_scatter_png <- multi_png_download(scatter_plot_fn, function() sprintf("multiomics_concordance_scatter_%s.png", gsub("[^A-Za-z0-9]+", "_", raw$cohort %||% "cohort")))

    output$pattern_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state())
      DT::dataTableOutput(ns("pattern_table"))
    })
    output$pattern_table <- DT::renderDataTable({
      tab <- multi_concordance_pattern_tally(req(filtered()))
      req(tab)
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    output$table_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state())
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_table_csv"), "Download table (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("conc_table"))
      )
    })
    output$conc_table <- DT::renderDataTable({
      df <- req(filtered())
      if (length(input$region) > 0 && "region" %in% colnames(df)) df <- df[df$region %in% input$region, , drop = FALSE]
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_table_csv <- downloadHandler(function() sprintf("multiomics_concordance_%s.csv", gsub("[^A-Za-z0-9]+", "_", raw$cohort %||% "cohort")),
                                             function(file) utils::write.csv(filtered(), file, row.names = FALSE))

    observe({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df) || is.null(multi_results)) return()
      multi_results$concordance <- list(df = df, cohort = raw$cohort)
    })
  })
}
