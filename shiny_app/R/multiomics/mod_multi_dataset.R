## R/multiomics/mod_multi_dataset.R
## Multi-Omics "Dataset" tab: browses the already-computed DIABLO/SNF
## integration, joint biomarker, gene<->CpG concordance, pathway enrichment,
## and patient/sample matching tables produced by
## Research_05_multiomics_sexstratified's own pipeline (see global.R's
## MULTI_TABLE_REGISTRY) and loads the selected one into the shared
## `multi_dataset` reactiveValues every Multi-Omics sub-module below reads
## from - the multi-omics equivalent of mod_cross_dataset.R.
##
## Like Cross-Omics, there is no upload path here: this module browses a
## pipeline's own precomputed, already-integrated output, it does not run
## DIABLO/SNF/nested-CV live over arbitrary uploaded matrices (that
## leakage-safe re-implementation is exactly the kind of work the source
## pipeline's own AUDIT.md shows required real bug-fixing effort even for
## its original authors - out of scope for this module, not silently
## skipped).

mod_multi_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Browse the precomputed DIABLO, SNF, joint-biomarker, concordance, and pathway tables every Multi-Omics sub-module below reads from."
)

mod_multi_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "empty-note", icon("circle-info"),
        "Independent of the Transcriptomics/Methylomics/Cross-Omics datasets - this reads the multi-omics pipeline's own already-integrated output tables (Research_05_multiomics_sexstratified), not whatever is currently loaded elsewhere."),
    if (!MULTI_DATA_AVAILABLE) {
      div(class = "empty-note", icon("triangle-exclamation"),
          "The multi-omics pipeline's output folder isn't available in this deployment - nothing to browse here.")
    } else {
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = "1. Table", status = "primary", solidHeader = FALSE,
            selectInput(ns("table_pick"), NULL, choices = names(MULTI_TABLE_REGISTRY), width = "100%"),
            p(class = "submodule-desc", "The preview on the right updates automatically - reading a table is cheap. Click \"Load table\" only when you want THIS table to become the one every Multi-Omics sub-module below reads from."),
            actionButton(ns("load_btn"), "Load table", icon = icon("upload"), class = "btn-primary btn-sm")
          ),
          box(
            width = NULL, title = "About this pipeline", status = "primary", solidHeader = FALSE,
            p("Independent RA anti-TNF cohort (Tao et al. 2021; GSE138653 methylation EPIC array + GSE138746 RNA-seq; 80 patients, baseline, pre-treatment)."),
            p("Both DIABLO (mixOmics, supervised) and SNF (SNFtool, unsupervised patient-similarity fusion) were run per sex x drug cell, with strictly train-fold-only feature selection under nested/leave-one-out cross-validation."),
            p(tags$strong("Read before trusting these numbers:"), " the pipeline's own audit found that most cells' AUROC 95% CIs include chance performance - see the Overview tab for the honest per-cell breakdown.")
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

mod_multi_dataset_server <- function(id, multi_dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    picked <- reactive({
      req(input$table_pick)
      res <- multi_read_registry_table(input$table_pick)
      if (!res$ok) return(NULL)
      res$df
    })

    output$preview_ui <- renderUI({
      df <- tryCatch(picked(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state("Select a table to preview it here."))
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
      multi_dataset$table_label <- input$table_pick
      multi_dataset$df <- df
      multi_dataset$source <- sprintf("Multi-omics pipeline: %s", input$table_pick)
      output$load_message <- renderUI(
        div(class = "empty-note", icon("check"),
            sprintf("Loaded \"%s\" (%s rows) for every Multi-Omics sub-module to read.", input$table_pick, format(nrow(df), big.mark = ",")))
      )
    })
  })
}
