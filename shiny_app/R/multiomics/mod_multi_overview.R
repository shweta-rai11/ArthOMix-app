## R/multiomics/mod_multi_overview.R
## Submodule: Cohort Overview & Sample Harmonization - reads the pipeline's
## own patient/sample matching table and QC summaries directly (independent
## of whichever table is "loaded" on the Dataset tab, since this is
## foundational cohort information every other sub-module assumes), and
## surfaces the pipeline's own excludes_chance verdicts across all six
## analysis cells up front - "read this before trusting these numbers"
## rather than only mentioning it in a caption once. Spec: sample
## harmonization must never be silently merged - shared vs. omics-specific
## sample counts are shown explicitly, per cell.

mod_multi_overview_config <- list(
  id = "overview", title = "Cohort Overview & Sample Harmonization", icon = "users", group = "Data",
  description = "Which patients have both RNA-seq and methylation, per analysis cell - and an honest summary of which cells' integrated models actually beat chance."
)

mod_multi_overview_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "About sample harmonization", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Multi-omics integration below only ever uses patients with BOTH RNA-seq and methylation available (\"matched-sample integration\") - no dataset is silently merged on mismatched samples."),
        actionButton(ns("load_btn"), "Load cohort tables", icon = icon("users"), class = "btn-primary btn-sm", width = "100%")
      )
    ),
    column(
      8,
      box(width = NULL, title = "Quality-control scorecard", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Computed from what's actually loaded this session - not all green by default."),
          uiOutput(ns("qc_scorecard"))),
      uiOutput(ns("honesty_banner")),
      uiOutput(ns("harmonization_ui")),
      uiOutput(ns("qc_ui")),
      box(width = NULL, title = "Analysis summary", status = "primary", solidHeader = FALSE,
          uiOutput(ns("analysis_summary_ui")))
    )
  )
}

mod_multi_overview_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    matching <- reactiveVal(NULL)
    qc <- reactiveValues(rna = NULL, meth = NULL)
    summary36 <- reactiveVal(NULL)

    observeEvent(input$load_btn, {
      m <- multi_read_registry_table("Patient sample matching (all 80 patients)")
      r <- multi_read_registry_table("RNA-seq QC summary")
      me <- multi_read_registry_table("Methylation QC summary")
      s36 <- multi_read_registry_table("Master six-part summary (integrated vs single-omics)")
      if (!m$ok) { showNotification(m$error, type = "error"); return() }
      matching(m$df)
      qc$rna <- if (r$ok) r$df else NULL
      qc$meth <- if (me$ok) me$df else NULL
      summary36(if (s36$ok) s36$df else NULL)
      showNotification("Loaded cohort and QC tables.", type = "message")
    }, ignoreInit = TRUE)

    harmonization <- reactive({
      req(matching())
      multi_sample_harmonization(matching())
    })

    output$honesty_banner <- renderUI({
      s <- summary36()
      if (is.null(s) || !"excludes_chance" %in% colnames(s)) return(NULL)
      n_total <- nrow(s)
      n_ex <- sum(s$excludes_chance %in% TRUE)
      ex_rows <- s[s$excludes_chance %in% TRUE, , drop = FALSE]
      div(
        class = "empty-note", style = "border-color: var(--color-warning, #eda100);",
        icon("triangle-exclamation"),
        tags$strong(sprintf(" %d of %d method x cell results below exclude chance performance (95%% CI entirely above 0.5).", n_ex, n_total)),
        if (nrow(ex_rows) > 0) tags$span(
          sprintf(" The rest should be treated as exploratory, not validated diagnostic models. Cell(s) that do exclude chance: %s.",
                  paste(sprintf("%s-%s (%s, AUROC %.2f)", ex_rows$sex, ex_rows$drug, ex_rows$method, ex_rows$auroc), collapse = "; "))
        ) else tags$span(" None of the six analysis cells' integrated models currently exclude chance performance at this sample size.")
      )
    })

    output$harmonization_ui <- renderUI({
      h <- tryCatch(harmonization(), error = function(e) NULL)
      if (is.null(h) || !isTRUE(h$ok)) return(multi_empty_state("Click \"Load cohort tables\" to see sample harmonization here."))
      tagList(
        box(
          width = NULL, title = "Sample harmonization (whole cohort)", status = "primary", solidHeader = FALSE,
          div(style = "display:flex; gap:10px; flex-wrap:wrap;",
              lapply(list(
                list(label = "Total patients", value = h$n_total, color = "blue"),
                list(label = "RNA-seq available", value = h$n_rna, color = "aqua"),
                list(label = "Methylation available", value = h$n_meth, color = "aqua"),
                list(label = "Matched (both omics)", value = h$n_matched, color = "violet"),
                list(label = "RNA-only", value = h$n_rna_only, color = "orange"),
                list(label = "Methylation-only", value = h$n_meth_only, color = "orange")
              ), function(c) {
                div(style = "flex: 1 1 140px; text-align:center; padding:10px;", class = "card",
                    div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[c$color]]), format(c$value, big.mark = ",")),
                    div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", c$label))
              })
          ),
          p(class = "submodule-desc", style = "margin-top:10px;",
            sprintf("Integration below uses matched-sample integration: %s of %s patients have both omics layers.", format(h$n_matched, big.mark = ","), format(h$n_total, big.mark = ",")))
        ),
        box(
          width = NULL, title = "Sample overlap", status = "primary", solidHeader = FALSE,
          multi_plot_or_empty(venn_plot_fn, ns("overlap_venn"), "Not enough data to draw a sample-overlap diagram.", height = "300px"),
          div(class = "table-toolbar", downloadButton(ns("dl_venn_png"), "Download plot (PNG)", class = "btn-sm"))
        ),
        box(
          width = NULL, title = "Sample counts per analysis cell", status = "primary", solidHeader = FALSE,
          DT::dataTableOutput(ns("cell_table"))
        )
      )
    })
    venn_plot_fn <- reactive(multi_sample_overlap_venn(matching()))
    output$overlap_venn <- renderPlot(venn_plot_fn())
    output$dl_venn_png <- multi_png_download(venn_plot_fn, function() "multiomics_sample_overlap_venn.png")
    output$cell_table <- DT::renderDataTable({
      h <- req(harmonization())
      DT::datatable(h$by_cell, rownames = FALSE, options = list(dom = "t", pageLength = 10), class = "stripe hover compact")
    })

    output$qc_ui <- renderUI({
      if (is.null(qc$rna) && is.null(qc$meth)) return(NULL)
      tagList(
        if (!is.null(qc$rna)) box(width = NULL, title = "RNA-seq QC summary (per cell type)", status = "primary", solidHeader = FALSE,
                                    DT::dataTableOutput(ns("rna_qc_table"))),
        if (!is.null(qc$meth)) box(width = NULL, title = "Methylation QC summary", status = "primary", solidHeader = FALSE,
                                     DT::dataTableOutput(ns("meth_qc_table")))
      )
    })
    output$rna_qc_table <- DT::renderDataTable({
      req(qc$rna)
      DT::datatable(qc$rna, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })
    output$meth_qc_table <- DT::renderDataTable({
      req(qc$meth)
      DT::datatable(qc$meth, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    observe({
      h <- tryCatch(harmonization(), error = function(e) NULL)
      if (is.null(h) || is.null(multi_results)) return()
      multi_results$overview <- list(harmonization = h, summary36 = summary36())
    })

    output$qc_scorecard <- renderUI({
      items <- multi_qc_scorecard(multi_results)
      color_for <- c(pass = ARTHOMIX_COLORS$aqua, warn = ARTHOMIX_COLORS$yellow, fail = ARTHOMIX_COLORS$red)
      icon_for <- c(pass = "circle-check", warn = "circle-exclamation", fail = "circle-xmark")
      tags$ul(
        style = "list-style:none; padding-left:0; margin:0;",
        lapply(items, function(it) tags$li(
          style = "display:flex; align-items:baseline; gap:8px; padding:4px 0; border-bottom:1px solid var(--color-grid, #eee);",
          span(style = sprintf("color:%s; width:16px; flex-shrink:0;", color_for[[it$status]]), icon(icon_for[[it$status]])),
          span(style = "font-weight:600; min-width:220px;", it$label),
          span(style = "color:var(--color-ink-muted, #898781); font-size:0.88em;", it$detail)
        ))
      )
    })

    output$analysis_summary_ui <- renderUI({
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_summary_csv"), "Download summary (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("analysis_summary_table"))
      )
    })
    analysis_summary_df <- reactive(multi_analysis_summary_table(multi_dataset %||% list(), multi_results))
    output$analysis_summary_table <- DT::renderDataTable({
      DT::datatable(analysis_summary_df(), rownames = FALSE, options = list(dom = "t", pageLength = 20), class = "stripe hover compact")
    })
    output$dl_summary_csv <- downloadHandler(function() "multiomics_analysis_summary.csv", function(file) utils::write.csv(analysis_summary_df(), file, row.names = FALSE))
  })
}
