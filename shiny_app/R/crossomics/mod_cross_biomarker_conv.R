## R/crossomics/mod_cross_biomarker_conv.R
## Submodule: Biomarker Convergence - a live, reconfigurable reimplementation
## of cross_Omics_Sexstratified_COPY/scripts/01_cross_omics_eqtl_mqtl_biomarkers.R
## (see crossomics_biomarkerconv_helpers.R for the join logic itself). Lets
## the join that script already ran once (eQTL-MR gene panel x mQTL-MR CpG
## panel x raw DEG/DMP/DMR tables, per sex) be reconfigured live - different
## FDR/effect-size thresholds, a different vote cutoff - rather than only
## browsing its already-computed output, which the Dataset tab already does.
##
## No new statistics are computed here - only joining and threshold-based
## labeling of already-computed per-layer results (spec/report language
## throughout: "evidence present in N of 3 layers", never "confirmed
## biomarker").

mod_cross_biomarker_conv_config <- list(
  id = "biomarkerconv", title = "Biomarker Convergence", icon = "diagram-project", group = "Data",
  description = "Joins the eQTL-MR gene panel, mQTL-MR CpG panel, and raw DEG/DMP/DMR tables into one sex-stratified biomarker/annotation table, live and reconfigurable."
)

mod_cross_biomarker_conv_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Cohort", status = "primary", solidHeader = FALSE,
        radioButtons(ns("sex"), "Sex", choices = c("Female" = "female", "Male" = "male"), selected = "female", inline = TRUE),
        p(class = "empty-note", icon("circle-info"),
          "Reads this project's own precomputed eQTL-MR gene panel, mQTL-MR CpG panel, and raw DEG/DMP/DMR tables directly - independent of whatever is currently loaded on the Transcriptomics/Methylomics tabs or the Expression x Methylation sub-module."),
        actionButton(ns("run_join"), "Run Join", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")
      ),
      box(
        width = NULL, title = "2. Reconfigure Thresholds", status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        numericInput(ns("vote_threshold"), "mQTL candidate CpG: min ensemble votes", value = 2, min = 1, max = 3, step = 1),
        numericInput(ns("dmp_genomewide_fdr"), "DMP genome-wide-significant: FDR <", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("dmp_suggestive_p"), "DMP suggestive: raw P <", value = 1e-4, min = 0, max = 1, step = 1e-5),
        radioButtons(ns("mqtl_sig_basis"), "mQTL-MR significance basis", choices = c("Nominal P (matches original script)" = "nominal_p", "BH-FDR (recomputed over this run)" = "fdr"), selected = "nominal_p"),
        numericInput(ns("mqtl_sig_cutoff"), "mQTL-MR significance cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("dmr_fdr"), "DMR significant: FDR <", value = 0.05, min = 0, max = 1, step = 0.01),
        p(class = "submodule-desc", "Defaults match the original pipeline script's own hardcoded values. Changing these only relabels which rows count as significant - it never re-runs the upstream eQTL-MR/mQTL-MR/DEG/DMP/DMR analyses themselves.")
      ),
      box(
        width = NULL, title = "Analysis Settings / Reproducibility", status = "primary", solidHeader = FALSE,
        uiOutput(ns("provenance_ui"))
      )
    ),
    column(
      8,
      uiOutput(ns("summary_cards")),
      tabsetPanel(
        id = ns("result_tabs"), type = "tabs",
        tabPanel("Results Table", br(), uiOutput(ns("results_table_ui"))),
        tabPanel("Evidence Overlap", br(), uiOutput(ns("overlap_ui"))),
        tabPanel("Downloads", br(), uiOutput(ns("downloads_ui")))
      )
    )
  )
}

mod_cross_biomarker_conv_server <- function(id, cross_dataset, cross_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    bc <- reactiveValues(df = NULL, params = NULL, sex = NULL, run_at = NULL, provenance = NULL)
    active_filter <- reactiveVal("All")

    observeEvent(input$run_join, {
      params <- list(
        vote_threshold = input$vote_threshold, dmp_genomewide_fdr = input$dmp_genomewide_fdr,
        dmp_suggestive_p = input$dmp_suggestive_p, mqtl_sig_basis = input$mqtl_sig_basis,
        mqtl_sig_cutoff = input$mqtl_sig_cutoff, dmr_fdr = input$dmr_fdr
      )
      res <- cx_bc_build_join(input$sex, params)
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      bc$df <- res$df
      bc$params <- res$params
      bc$sex <- input$sex
      bc$run_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      bc$provenance <- cx_bc_build_provenance(input$sex, res$params, bc$run_at)
      active_filter("All")
      if (!is.null(cross_results)) {
        cross_results$biomarkerconv <- list(df = res$df, sex = input$sex, params = res$params, run_at = bc$run_at)
      }
      showNotification(sprintf("Join complete - %s genes (%s).", format(nrow(res$df), big.mark = ","), toupper(input$sex)), type = "message")
    })

    output$provenance_ui <- renderUI({
      if (is.null(bc$provenance)) return(div(class = "empty-note", icon("circle-info"), "Run the join to see its parameters here."))
      tags$ul(style = "padding-left: 18px; font-size: 0.85em;", lapply(bc$provenance, tags$li))
    })

    output$summary_cards <- renderUI({
      if (is.null(bc$df)) return(cx_empty_state())
      df <- bc$df
      card <- function(label, value, key, color = "blue") {
        div(class = "card", style = "flex: 1 1 150px; cursor:pointer; text-align:center; padding:10px;",
            onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("card_click"), key),
            div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink), format(value, big.mark = ",")),
            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label))
      }
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          card("Genes (union)", nrow(df), "All"),
          card("In eQTL-MR panel", sum(df$in_eQTL_MR_panel), "eqtl", "aqua"),
          card("In mQTL-MR panel", sum(df$in_mQTL_MR_panel), "mqtl", "aqua"),
          card("DEG significant", sum(df$DEG_significant, na.rm = TRUE), "deg", "violet"),
          card("Methylation significant (DMP or DMR)", sum(df$methylation_significant, na.rm = TRUE), "meth", "violet"),
          card("Evidence in all applicable layers", sum(df$n_evidence_layers >= 3), "layers3", "red")
      )
    })
    observeEvent(input$card_click, {
      key <- input$card_click
      active_filter(switch(key,
        "eqtl" = "eqtl", "mqtl" = "mqtl", "deg" = "deg", "meth" = "meth", "layers3" = "layers3", "All"))
    })

    filtered_df <- reactive({
      req(bc$df)
      df <- bc$df
      switch(active_filter(),
        "eqtl" = df[df$in_eQTL_MR_panel %in% TRUE, , drop = FALSE],
        "mqtl" = df[df$in_mQTL_MR_panel %in% TRUE, , drop = FALSE],
        "deg" = df[df$DEG_significant %in% TRUE, , drop = FALSE],
        "meth" = df[df$methylation_significant %in% TRUE, , drop = FALSE],
        "layers3" = df[df$n_evidence_layers >= 3, , drop = FALSE],
        df
      )
    })

    output$results_table_ui <- renderUI({
      if (is.null(bc$df)) return(cx_empty_state())
      display_cols <- c("gene", "in_eQTL_MR_panel", "eQTL_MR_OR", "eQTL_MR_FDR", "eQTL_MHC_region",
                         "in_mQTL_MR_panel", "mQTL_candidate_cpg", "mQTL_MR_pval", "mQTL_MR_significant",
                         "DEG_logFC", "DEG_adjP", "DEG_direction", "DEG_significant",
                         "DMP_top_cpg", "DMP_dbeta", "DMP_fdr_bacon", "DMP_genomewide_significant", "DMP_suggestive",
                         "dmr_id", "DMR_meandiff", "DMR_fdr", "DMR_significant", "n_evidence_layers")
      display_cols <- intersect(display_cols, colnames(bc$df))
      tagList(
        fluidRow(
          column(6, selectInput(ns("filter_pick"), "Quick filter", choices = c("All genes" = "All", "In eQTL-MR panel" = "eqtl", "In mQTL-MR panel" = "mqtl", "DEG significant" = "deg", "Methylation significant" = "meth", "Evidence in all applicable layers" = "layers3"), selected = "All")),
          column(6, selectizeInput(ns("visible_cols"), "Visible columns", choices = display_cols, selected = display_cols, multiple = TRUE))
        ),
        div(class = "table-toolbar", downloadButton(ns("dl_table_csv"), "CSV", class = "btn-sm"), downloadButton(ns("dl_table_xlsx"), "XLSX", class = "btn-sm")),
        DT::dataTableOutput(ns("results_table"))
      )
    })
    observeEvent(input$filter_pick, active_filter(input$filter_pick), ignoreInit = TRUE)

    table_df <- reactive({
      df <- filtered_df()
      cols <- intersect(input$visible_cols %||% colnames(df), colnames(df))
      if (length(cols) == 0) cols <- colnames(df)
      df[, cols, drop = FALSE]
    })
    output$results_table <- DT::renderDataTable({
      req(table_df())
      DT::datatable(table_df(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    output$overlap_ui <- renderUI({
      if (is.null(bc$df)) return(cx_empty_state())
      tagList(
        p(class = "submodule-desc", "How many genes have evidence in each combination of layers - eQTL-MR panel membership, DEG significance, and methylation significance (DMP genome-wide or DMR)."),
        plotOutput(ns("overlap_plot"), height = "420px")
      )
    })
    output$overlap_plot <- renderPlot({
      df <- bc$df
      validate(need(!is.null(df), "Run the join first."))
      sets <- list(`eQTL-MR panel` = df$gene[df$in_eQTL_MR_panel %in% TRUE],
                   `DEG significant` = df$gene[df$DEG_significant %in% TRUE],
                   `Methylation significant` = df$gene[df$methylation_significant %in% TRUE])
      ggVennDiagram::ggVennDiagram(sets, label = "count") +
        ggplot2::scale_fill_gradient(low = "white", high = ARTHOMIX_COLORS$blue) + theme_arthomix() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
    })

    output$downloads_ui <- renderUI({
      if (is.null(bc$df)) return(cx_empty_state())
      tagList(
        downloadButton(ns("dl_table_csv"), "Results (CSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_xlsx"), "Results (XLSX)", class = "btn-sm")
      )
    })
    output$dl_table_csv <- downloadHandler(function() paste0("biomarker_convergence_", bc$sex, "_", Sys.Date(), ".csv"), function(file) utils::write.csv(bc$df, file, row.names = FALSE))
    output$dl_table_xlsx <- downloadHandler(function() paste0("biomarker_convergence_", bc$sex, "_", Sys.Date(), ".xlsx"), function(file) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) { utils::write.csv(bc$df, file, row.names = FALSE); return() }
      openxlsx::write.xlsx(bc$df, file)
    })
  })
}
