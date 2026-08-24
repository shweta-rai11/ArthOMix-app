## R/multiomics/mod_multi_biomarker.R
## Submodule: Joint Biomarker Discovery - the pipeline's own candidate
## multi-omics biomarker tables (DIABLO panel features, each row already
## carrying panel-level nested-CV performance and the pipeline's own
## biomarker_status call). This tab lets the confidence tier shown be
## relabeled live against adjustable thresholds, applied to the ALREADY-
## COMPUTED performance_auroc/ci_lo/ci_hi columns - no combined statistic is
## invented, no p-value is created, and DIABLO is never re-fit.
##
## Language throughout mirrors the source pipeline's own honest labels
## ("statistically_supported" vs. "exploratory_not_significant") rather than
## a stronger word like "confirmed" or "validated".

mod_multi_biomarker_config <- list(
  id = "biomarker", title = "Joint Biomarker Discovery", icon = "star", group = "Biomarker modeling",
  description = "Candidate multi-omics biomarkers (DIABLO panel features spanning transcriptomics + methylomics), with live-adjustable confidence relabeling against the panel's own nested-CV performance."
)

MULTI_BIOMARKER_COHORTS <- c(
  "Drug x sex (Etanercept panel)" = "Candidate multi-omics biomarkers — drug x sex (Etanercept panel)",
  "Response (drug-pooled)" = "Candidate multi-omics biomarkers — response (drug-pooled)"
)

mod_multi_biomarker_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Cohort", status = "primary", solidHeader = FALSE,
        selectInput(ns("cohort"), NULL, choices = MULTI_BIOMARKER_COHORTS, width = "100%"),
        actionButton(ns("load_btn"), "Load table", icon = icon("database"), class = "btn-primary btn-sm", width = "100%")
      ),
      box(
        width = NULL, title = "2. Confidence thresholds", status = "primary", solidHeader = FALSE, collapsible = TRUE,
        numericInput(ns("min_auroc"), "Minimum panel AUROC", value = 0.6, min = 0.5, max = 1, step = 0.01),
        checkboxInput(ns("require_excludes_chance"), "Require the panel's 95% CI to exclude chance (0.5)", value = TRUE),
        p(class = "submodule-desc", "Relabels which rows count as \"high confidence\" from the panel's own already-computed AUROC/CI - no new statistic is calculated here.")
      )
    ),
    column(
      8,
      uiOutput(ns("summary_cards")),
      tabsetPanel(
        id = ns("tabs"), type = "tabs",
        tabPanel("Candidates", br(), uiOutput(ns("table_ui"))),
        tabPanel("By omics layer", br(), uiOutput(ns("omics_ui")))
      )
    )
  )
}

mod_multi_biomarker_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    raw <- reactiveValues(df = NULL, cohort = NULL)

    observeEvent(input$load_btn, {
      res <- multi_read_registry_table(input$cohort)
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      raw$df <- res$df
      raw$cohort <- names(MULTI_BIOMARKER_COHORTS)[MULTI_BIOMARKER_COHORTS == input$cohort]
      showNotification(sprintf("Loaded %s candidate features.", format(nrow(res$df), big.mark = ",")), type = "message")
    })

    labeled <- reactive({
      req(raw$df)
      multi_biomarker_relabel(raw$df, list(min_auroc = input$min_auroc, require_excludes_chance = input$require_excludes_chance))
    })

    output$summary_cards <- renderUI({
      df <- tryCatch(labeled(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state("Load a cohort to see candidate biomarkers here."))
      n_features <- length(unique(df$feature))
      n_high <- length(unique(df$feature[df$display_confidence == "High confidence (panel-level)"]))
      by_omics <- table(df$omics)
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          lapply(list(
            list(label = "Candidate features", value = n_features, color = "blue"),
            list(label = "High confidence (panel-level)", value = n_high, color = "aqua"),
            list(label = "Transcriptomic features", value = sum(by_omics[grepl("transcript", names(by_omics))]), color = "violet"),
            list(label = "Methylomic features", value = sum(by_omics[grepl("methyl", names(by_omics))]), color = "orange")
          ), function(c) div(class = "card", style = "flex:1 1 150px; text-align:center; padding:10px;",
                              div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[c$color]]), format(c$value, big.mark = ",")),
                              div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", c$label)))
      )
    })

    output$table_ui <- renderUI({
      df <- tryCatch(labeled(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state())
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_csv"), "Download candidates (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("cand_table"))
      )
    })
    output$cand_table <- DT::renderDataTable({
      df <- req(labeled())
      cols <- intersect(c("sex", "drug", "omics", "feature", "loading", "rank_within_block", "panel_source_model",
                           "performance_model", "performance_auroc", "performance_ci_lo", "performance_ci_hi",
                           "biomarker_status", "display_confidence"), colnames(df))
      DT::datatable(df[, cols, drop = FALSE], rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_csv <- downloadHandler(function() sprintf("multiomics_biomarker_candidates_%s.csv", gsub("[^A-Za-z0-9]+", "_", raw$cohort %||% "cohort")),
                                       function(file) utils::write.csv(labeled(), file, row.names = FALSE))

    output$omics_ui <- renderUI({
      df <- tryCatch(labeled(), error = function(e) NULL)
      if (is.null(df)) return(multi_empty_state())
      tagList(
        p(class = "submodule-desc", "How many candidate features come from each omics layer, per sex - a compact multi-omics panel should draw evidence from more than one layer where possible, not just whichever has more available features."),
        DT::dataTableOutput(ns("omics_table"))
      )
    })
    output$omics_table <- DT::renderDataTable({
      df <- req(labeled())
      tab <- as.data.frame(table(sex = df$sex, omics = df$omics))
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    observe({
      df <- tryCatch(labeled(), error = function(e) NULL)
      if (is.null(df) || is.null(multi_results)) return()
      multi_results$biomarker <- list(df = df, cohort = raw$cohort)
    })
  })
}
