## R/multiomics/mod_multi_pathway.R
## Submodule: Pathway-Level Integration - the pipeline's own GO/KEGG/
## Reactome over-representation analysis (clusterProfiler/ReactomePA,
## org.Hs.eg.db) on the DIABLO Etanercept-panel genes (script 16), against a
## PBMC-expressed-gene background universe. Turns a gene list into
## interpretable mechanism; no enrichment test is re-run here.

mod_multi_pathway_config <- list(
  id = "pathway", title = "Pathway-Level Integration", icon = "sitemap", group = "Interpretation",
  description = "GO/KEGG/Reactome enrichment on the multi-omics candidate biomarker panel - which biological processes and pathways the panel points to."
)

mod_multi_pathway_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Load", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Pathway enrichment is available for the drug x sex Etanercept panel (the only cohort the pipeline ran enrichment on)."),
        actionButton(ns("load_btn"), "Load pathway table", icon = icon("sitemap"), class = "btn-primary btn-sm", width = "100%")
      ),
      box(
        width = NULL, title = "2. Filters", status = "primary", solidHeader = FALSE, collapsible = TRUE,
        selectizeInput(ns("sex"), "Sex", choices = c("Female" = "female", "Male" = "male"), multiple = TRUE, options = list(placeholder = "Both sexes")),
        selectizeInput(ns("source"), "Source", choices = NULL, multiple = TRUE, options = list(placeholder = "All sources (GO/KEGG/Reactome)")),
        numericInput(ns("padj"), "Adjusted P <", value = 0.25, min = 0, max = 1, step = 0.01),
        selectInput(ns("top_n"), "Show top", choices = c("10" = 10, "20" = 20, "50" = 50), selected = 20),
        p(class = "submodule-desc", "Small candidate panels rarely survive strict multiple-testing correction - a looser default is used so the pathway list isn't empty by construction; every row still shows its own p.adjust for the reader to judge.")
      )
    ),
    column(
      8,
      uiOutput(ns("summary_ui")),
      tabsetPanel(
        id = ns("tabs"), type = "tabs",
        tabPanel("Dotplot", br(), uiOutput(ns("dot_ui"))),
        tabPanel("Table", br(), uiOutput(ns("table_ui")))
      )
    )
  )
}

mod_multi_pathway_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    raw <- reactiveVal(NULL)

    observeEvent(input$load_btn, {
      res <- multi_read_registry_table("Pathway enrichment — drug x sex (Etanercept panel)")
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      raw(res$df)
      if ("source" %in% colnames(res$df)) updateSelectizeInput(session, "source", choices = sort(unique(res$df$source)), server = FALSE)
      showNotification(sprintf("Loaded %s enriched terms.", format(nrow(res$df), big.mark = ",")), type = "message")
    })

    filtered <- reactive({
      req(raw())
      df <- raw()
      if (length(input$sex) > 0) df <- df[tolower(df$sex) %in% tolower(input$sex), , drop = FALSE]
      if (length(input$source) > 0) df <- df[df$source %in% input$source, , drop = FALSE]
      df <- df[!is.na(df$p.adjust) & df$p.adjust < (input$padj %||% 1), , drop = FALSE]
      df
    })

    output$summary_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state("Load the pathway table to see results here."))
      div(class = "empty-note", icon("chart-simple"), sprintf(" %s terms pass the current filters (of %s total).", format(nrow(df), big.mark = ","), format(nrow(raw()), big.mark = ",")))
    })

    output$dot_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) return(multi_empty_state("No terms pass the current filters."))
      tagList(
        multi_plot_or_empty(dot_plot_fn, ns("dotplot"), height = "520px"),
        div(class = "table-toolbar", downloadButton(ns("dl_dot_png"), "Download plot (PNG)", class = "btn-sm")),
        p(class = "submodule-desc", tags$em("How to read this: ", "dot size is the number of panel genes in that term; color is the adjusted P-value (bluer = smaller). A significant term is a statistical association with the gene set, not proof that the pathway drives the phenotype."))
      )
    })
    dot_plot_fn <- reactive(multi_pathway_dotplot(filtered(), top_n = as.integer(input$top_n %||% 20)))
    output$dotplot <- renderPlot(dot_plot_fn())
    output$dl_dot_png <- multi_png_download(dot_plot_fn, function() "multiomics_pathway_dotplot.png")

    output$table_ui <- renderUI({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state())
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_table_csv"), "Download table (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("pathway_table"))
      )
    })
    output$pathway_table <- DT::renderDataTable({
      df <- req(filtered())
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_table_csv <- downloadHandler(function() "multiomics_pathway_enrichment.csv", function(file) utils::write.csv(filtered(), file, row.names = FALSE))

    observe({
      df <- tryCatch(filtered(), error = function(e) NULL)
      if (is.null(df) || is.null(multi_results)) return()
      multi_results$pathway <- list(df = df)
    })
  })
}
