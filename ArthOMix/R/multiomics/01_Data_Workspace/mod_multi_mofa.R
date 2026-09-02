## R/multiomics/01_Data_Workspace/mod_multi_mofa.R
## "Integrated Analysis (MOFA2)" - the one part of the Multi-Omics module
## that computes on data the user supplies, rather than browsing the
## precomputed cohort's results. Mounted directly inside the Dataset
## Workspace tab (mod_multi_dataset.R), not as its own top-level sub-module -
## upload/validation/sample matching/preprocessing/batch diagnostics all
## happen there too, and publish one validated Active Multi-Omics Dataset
## into the shared `multi_dataset` reactiveValues; this file just adapts
## that shared dataset into the shape the nested MOFA2 sub-module
## (mod_multi_mofa_engine.R, unchanged) already expects and mounts it.
## (This module and its data structures still use the historical "live"
## prefix internally - e.g. multiomics_dataset_helpers.R, live_state,
## multi_results$live_qc/live_mofa - only the file/exported-function names
## and user-facing text were renamed.)
##
## Per spec section 34's "preloaded vs. uploaded must never contaminate each
## other": this only trains on `multi_dataset$layers` when
## `multi_dataset$active` is TRUE, whatever its source - it never reads the
## precomputed-cohort tables directly.

mod_multi_mofa_config <- list(
  id = "mofa", title = "Integrated Analysis (MOFA2)", icon = "chart-line",
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

    ## Dataset Summary table - relocated here from the Dataset Workspace tab's
    ## own box (mod_multi_dataset.R); still built from the same
    ## `multi_dataset$layer_meta` set by that tab's "Use Selected Datasets for
    ## Multi-Omics Analysis" (activate_btn) handler, via the shared
    ## `mo_summary_table()` helper (mod_multi_dataset.R).
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
