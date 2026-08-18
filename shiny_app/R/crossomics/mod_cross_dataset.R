## R/crossomics/mod_cross_dataset.R
## Cross-Omics "Dataset" tab: browses the already-computed biomarker-
## convergence and cross-omics MR tables produced by the cross-omics
## pipeline (Research_Q4_cross_Omics_sexstratified_COPY/cross_Omics_
## Sexstratified_COPY/scripts/01_cross_omics_eqtl_mqtl_biomarkers.R and
## 02_mr_stage_cross_omics.R - see global.R's CX_TABLE_REGISTRY), and loads
## the selected one into the shared `cross_dataset` reactiveValues every
## Cross-Omics sub-module below reads from - the cross-omics equivalent of
## mod_dataset.R's `dataset` / mod_methyl_dataset.R's `methyl_dataset`.
##
## Unlike those two, there is no upload path here: Cross-Omics doesn't run
## its own eQTL-MR/mQTL-MR/DEG/DMP joins live, it browses the join the
## pipeline already computed once both single-omics pipelines finished.
## Live re-computation (picking a different FDR/effect-size threshold and
## re-joining) is future work for the Biomarker Convergence and Cross-Omics
## MR sub-modules below (see mod_cross_biomarker_conv.R / mod_cross_mr_stage.R),
## not this tab.

mod_cross_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Browse the precomputed eQTL/mQTL/DEG/DMP biomarker-convergence and cross-omics MR tables for every Cross-Omics sub-module below to read from."
)

mod_cross_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "empty-note", icon("circle-info"),
        "Independent of the Transcriptomics and Methylomics datasets - this reads the pipeline's own already-joined output tables, not whatever is currently loaded on either single-omics Dataset tab."),
    if (!CX_DATA_AVAILABLE) {
      div(class = "empty-note", icon("triangle-exclamation"),
          "The cross-omics pipeline's output folder isn't available in this deployment - nothing to browse here.")
    } else {
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = "1. Table", status = "primary", solidHeader = FALSE,
            selectInput(ns("table_pick"), NULL, choices = names(CX_TABLE_REGISTRY), width = "100%"),
            p(class = "submodule-desc", "The preview on the right updates automatically as you change the table above - reading it is cheap, so there's nothing to wait on. Click \"Load table\" only when you want THIS table to become the one every Cross-Omics sub-module below reads from."),
            actionButton(ns("load_btn"), "Load table", icon = icon("upload"), class = "btn-primary btn-sm")
          ),
          box(
            width = NULL, title = "About this pipeline", status = "primary", solidHeader = FALSE,
            p("Stage 1 joins the eQTL-MR gene panel and mQTL-MR CpG panel (each already causally screened against its own single-omics GWAS) with the raw DEG and DMP tables into one biomarker/annotation table per sex."),
            p("Stage 2 runs a follow-up, single-instrument mQTL-MR analysis for the eQTL-MR-significant genes, since the two original panels share zero genes by design.")
          )
        ),
        column(
          8,
          box(width = NULL, title = "2. Preview", status = "primary", solidHeader = FALSE,
              uiOutput(ns("preview_ui"))),
          uiOutput(ns("load_message"))
        )
      )
    }
  )
}

mod_cross_dataset_server <- function(id, cross_dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    picked <- reactive({
      req(input$table_pick)
      load_default_cx_table(input$table_pick)
    })

    output$preview_ui <- renderUI({
      df <- tryCatch(picked(), error = function(e) NULL)
      if (is.null(df)) return(div(class = "empty-note", icon("circle-info"), "Select a table to preview it here."))
      tagList(
        p(class = "empty-note", icon("table"),
          sprintf("%s: %s rows x %s columns.", input$table_pick, format(nrow(df), big.mark = ","), ncol(df))),
        DT::dataTableOutput(ns("preview_table"))
      )
    })
    output$preview_table <- DT::renderDataTable({
      req(picked())
      DT::datatable(picked(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    outputOptions(output, "preview_table", suspendWhenHidden = FALSE)

    observeEvent(input$load_btn, {
      df <- picked()
      validate(need(!is.null(df), "Could not read this table."))
      cross_dataset$table_label <- input$table_pick
      cross_dataset$df <- df
      cross_dataset$source <- sprintf("Cross-omics pipeline: %s", input$table_pick)
      output$load_message <- renderUI(
        div(class = "empty-note", icon("check"),
            sprintf("Loaded \"%s\" (%s rows) for every Cross-Omics sub-module to read.", input$table_pick, format(nrow(df), big.mark = ",")))
      )
    })
  })
}
