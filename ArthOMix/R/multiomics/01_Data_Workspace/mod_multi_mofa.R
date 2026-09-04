## R/multiomics/01_Data_Workspace/mod_multi_mofa.R
## "Integrated Analysis (MOFA2)" - the one part of the Multi-Omics module
## that computes on data the user supplies, rather than browsing the

mod_multi_mofa_config <- list(
  id = "mofa", title = "MOFA", icon = "chart-line", group = "Data",
  description = "Run a real MOFA2 factor analysis on the Active Multi-Omics Dataset built on the Dataset Workspace tab (matched samples only)."
)

mod_multi_mofa_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(width = NULL, title = "Dataset Summary", status = "primary", solidHeader = FALSE,
        DT::dataTableOutput(ns("summary_table"))),
    uiOutput(ns("live_body_ui"))
  )
}

mod_multi_mofa_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    live_state <- reactiveValues(mats = NULL, meta = NULL)

    output$summary_table <- DT::renderDataTable({
      df <- mo_summary_table(multi_dataset$layer_meta %||% list())
      if (is.null(df)) df <- data.frame(Dataset = character(0), `Omics Type` = character(0), Samples = integer(0),
                                          Features = integer(0), Processing = character(0), Status = character(0), check.names = FALSE)
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    observe({
      if (!is.null(multi_dataset) && isTRUE(multi_dataset$active) && length(multi_dataset$layers %||% list()) >= 2) {
        live_state$mats <- multi_dataset$layers
        live_state$meta <- multi_dataset$sample_meta
      } else {
        live_state$mats <- NULL
        live_state$meta <- NULL
      }
    })

    output$live_body_ui <- renderUI({
      if (is.null(live_state$mats)) {
        return(tagList(
          div(class = "empty-note", icon("circle-info"),
              "No Active Multi-Omics Dataset yet. Build one above, then click \"Use Selected Datasets for Multi-Omics Analysis\"."),
          if (!is.null(multi_dataset) && identical(multi_dataset$source %||% "", "preloaded"))
            div(class = "empty-note", icon("triangle-exclamation"),
                "MOFA2 was not run on the preloaded dataset. See Multi-omics Integration for the preloaded cohort's DIABLO/SNF results.")
        ))
      }
      mod_multi_mofa_engine_ui(ns("mofa"))
    })

    mod_multi_mofa_engine_server("mofa", live_state, multi_results)

    observe({
      if (is.null(multi_results)) return()
      multi_results$live_qc <- list(n_layers = length(live_state$mats %||% list()), active_source = multi_dataset$source %||% NA_character_)
    })
  })
}
