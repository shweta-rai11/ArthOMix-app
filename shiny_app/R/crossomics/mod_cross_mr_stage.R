## R/crossomics/mod_cross_mr_stage.R
## Submodule: Cross-Omics MR - loads the pipeline's own already-run
## single-instrument mQTL-MR results (Wald ratio, GoDMC exposure -> Ishigaki
## 2022 RA outcome; cross_Omics_Sexstratified_COPY/results/mr_stage_eqtl_
## significant_genes_mqtl_mr.csv - see crossomics_mrstage_helpers.R) and
## INTEGRATES them with the Biomarker Convergence sub-module's DEG/DMP/DMR
## evidence to compute the Tier 1/2/3 classification described in the
## pipeline's own report but never implemented in code by either analysis
## script - see crossomics_mrstage_helpers.R's cx_mr_tier_rule_text() for the
## exact rule, shown verbatim in this module's UI. No MR statistics are
## re-run here, and no file outside cross_Omics_Sexstratified_COPY is read -
## every "Integration Filter" below relabels already-computed values, so
## changing one updates the Tier table instantly.

mod_cross_mr_stage_config <- list(
  id = "mrstage", title = "Cross-Omics MR", icon = "arrow-right-arrow-left", group = "Genetics",
  description = "Integrates the pipeline's own precomputed single-instrument mQTL-MR results with Biomarker Convergence's DEG/DMP/DMR evidence, computing the Tier 1/2/3 classification with instantly reconfigurable filters."
)

mod_cross_mr_stage_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. MR Data", status = "primary", solidHeader = FALSE,
        p(class = "empty-note", icon("circle-info"),
          "Loads this project's own precomputed single-instrument (Wald-ratio) mQTL-MR results directly - GoDMC exposure, Ishigaki et al. 2022 RA outcome. The MR analysis itself already ran; nothing is re-computed here."),
        actionButton(ns("load_mr"), "Load MR Results", icon = icon("database"), class = "btn-primary btn-sm", width = "100%"),
        uiOutput(ns("panel_info_ui"))
      ),
      box(
        width = NULL, title = "2. Integration Filters", status = "primary", solidHeader = FALSE,
        radioButtons(ns("sex"), "DEG/DMP/DMR evidence from", choices = c("Female" = "female", "Male" = "male"), selected = "female", inline = TRUE),
        numericInput(ns("fdr_cutoff"), "\"Credible\" mQTL-MR: FDR <", value = 0.05, min = 0, max = 1, step = 0.01),
        checkboxInput(ns("exclude_mhc"), "Exclude MHC-flagged genes from \"credible\"", value = TRUE),
        checkboxInput(ns("require_steiger"), "Require Steiger-consistent direction for \"credible\"", value = TRUE),
        p(class = "submodule-desc", "Defaults match the pipeline's own report. These relabel already-computed evidence - changing them updates the Tier table instantly, no re-run needed.")
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
        tabPanel("Tier Classification", br(), uiOutput(ns("tier_ui"))),
        tabPanel("Results Table", br(), uiOutput(ns("results_table_ui"))),
        tabPanel("Volcano", br(), uiOutput(ns("volcano_ui"))),
        tabPanel("Downloads", br(), uiOutput(ns("downloads_ui")))
      )
    )
  )
}

mod_cross_mr_stage_server <- function(id, cross_dataset, cross_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    mrs <- reactiveValues(df = NULL, loaded_at = NULL)

    observeEvent(input$load_mr, {
      res <- cx_mr_load_precomputed()
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      mrs$df <- res$df
      mrs$loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      showNotification(sprintf("Loaded %s MR instrument results.", format(nrow(res$df), big.mark = ",")), type = "message")
    })

    output$panel_info_ui <- renderUI({
      if (!CX_MR_DATA_AVAILABLE) return(div(class = "empty-note", icon("triangle-exclamation"), "Cross-Omics MR source data is not available in this deployment."))
      if (is.null(mrs$df)) return(NULL)
      div(class = "empty-note", icon("dna"), sprintf(
        "%s CpG-instrument results covering %s genes.", format(nrow(mrs$df), big.mark = ","), format(length(unique(mrs$df$gene)), big.mark = ",")))
    })

    ## Biomarker Convergence's own table for the selected sex - reused
    ## directly if that tab has already been loaded in this session (fast
    ## path), otherwise loaded fresh here via the same shared helpers (no
    ## dependency on that tab having been opened first).
    join_df <- reactive({
      if (!is.null(cross_results) && !is.null(cross_results$biomarkerconv) && identical(cross_results$biomarkerconv$sex, input$sex)) {
        return(cross_results$biomarkerconv$df)
      }
      j <- cx_bc_load_precomputed(input$sex)
      if (!isTRUE(j$ok)) return(NULL)
      cx_bc_relabel(j$df)
    })

    tier_filters <- reactive(list(fdr_cutoff = input$fdr_cutoff %||% 0.05, exclude_mhc = isTRUE(input$exclude_mhc), require_steiger = isTRUE(input$require_steiger)))

    tiers <- reactive({
      req(mrs$df)
      cx_mr_classify_tier(mrs$df, join_df(), tier_filters())
    })

    output$provenance_ui <- renderUI({
      if (is.null(mrs$df)) return(div(class = "empty-note", icon("circle-info"), "Load MR results to see provenance here."))
      tags$ul(style = "padding-left: 18px; font-size: 0.85em;",
              lapply(cx_mr_build_provenance(tier_filters(), nrow(mrs$df), if (!is.null(join_df())) input$sex else NULL, mrs$loaded_at), tags$li))
    })

    output$summary_cards <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      df <- mrs$df; t <- tiers()
      card <- function(label, value, color = "blue") {
        div(class = "card", style = "flex: 1 1 150px; text-align:center; padding:10px;",
            div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink), format(value, big.mark = ",")),
            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label))
      }
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          card("Instruments", nrow(df)), card("Genes", length(unique(df$gene)), "aqua"),
          card("Significant (current FDR)", sum(!is.na(df$FDR) & df$FDR < (input$fdr_cutoff %||% 0.05), na.rm = TRUE), "red"),
          card("Tier 1", sum(t$tier == "Tier 1", na.rm = TRUE), "violet"),
          card("Tier 2", sum(t$tier == "Tier 2", na.rm = TRUE), "orange"),
          card("Tier 3", sum(t$tier == "Tier 3", na.rm = TRUE), "ink_muted")
      )
    })

    output$tier_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      t <- tiers()
      tagList(
        h5("How Tier is computed"),
        tags$ul(style = "padding-left: 18px; font-size: 0.88em;", lapply(cx_mr_tier_rule_text(tier_filters()), tags$li)),
        if (!isTRUE(t$tier_evidence_available[1])) div(class = "empty-note", icon("triangle-exclamation"), sprintf("DEG/methylation evidence for %s was not available - every gene defaults to Tier 3 (no gene can be shown as Tier 1/2 without it).", toupper(input$sex))),
        DT::dataTableOutput(ns("tier_table"))
      )
    })
    output$tier_table <- DT::renderDataTable({
      t <- tiers()
      cols <- intersect(c("gene", "DEG_significant", "methylation_significant", "credible_mQTL_MR", "tier"), colnames(t))
      DT::datatable(t[order(t$tier), cols, drop = FALSE], rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    output$results_table_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      tagList(
        p(class = "submodule-desc", "Mendelian randomization estimates - valid under the standard instrumental-variable assumptions. A single instrument per exposure cannot be tested for validity via heterogeneity; treat individual estimates as an association consistent with a causal effect, not as proof of one."),
        div(class = "table-toolbar", downloadButton(ns("dl_table_csv"), "CSV", class = "btn-sm"), downloadButton(ns("dl_table_xlsx"), "XLSX", class = "btn-sm")),
        DT::dataTableOutput(ns("results_table"))
      )
    })
    output$results_table <- DT::renderDataTable({
      req(mrs$df)
      cols <- intersect(c("gene", "cpg", "SNP", "nsnp", "b", "se", "pval", "OR", "OR_lo", "OR_hi", "FDR", "mr_significant", "steiger_dir", "steiger_pval"), colnames(mrs$df))
      DT::datatable(mrs$df[, cols, drop = FALSE], rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    output$volcano_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      plotOutput(ns("volcano_plot"), height = "460px")
    })
    output$volcano_plot <- renderPlot({
      df <- mrs$df
      validate(need(!is.null(df) && nrow(df) > 0, "Load MR results first."))
      cutoff <- input$fdr_cutoff %||% 0.05
      d <- df[!is.na(df$b) & !is.na(df$FDR), , drop = FALSE]
      d$neglog10fdr <- -log10(pmax(d$FDR, 1e-300))
      d$sig <- d$FDR < cutoff
      ggplot2::ggplot(d, ggplot2::aes(x = b, y = neglog10fdr, color = sig)) +
        ggplot2::geom_point(alpha = 0.75, size = 1.8) +
        ggplot2::geom_hline(yintercept = -log10(cutoff), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
        ggrepel::geom_text_repel(data = d[d$sig %in% TRUE, , drop = FALSE], ggplot2::aes(label = gene), color = ARTHOMIX_COLORS$ink, size = 3, max.overlaps = 30) +
        ggplot2::scale_color_manual(values = c(`TRUE` = ARTHOMIX_COLORS$red, `FALSE` = ARTHOMIX_COLORS$ink_muted), labels = c(`TRUE` = sprintf("FDR < %s", cutoff), `FALSE` = "Not significant"), name = NULL) +
        ggplot2::labs(x = "MR estimate (b, log-odds scale)", y = "-log10(FDR)", title = "Cross-Omics MR (Wald ratio) volcano") +
        theme_arthomix()
    })

    output$downloads_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      tagList(
        downloadButton(ns("dl_table_csv"), "MR results (CSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_xlsx"), "MR results (XLSX)", class = "btn-sm"),
        downloadButton(ns("dl_tier_csv"), "Tier classification (CSV)", class = "btn-sm")
      )
    })
    output$dl_table_csv <- downloadHandler(function() paste0("cross_omics_mr_", Sys.Date(), ".csv"), function(file) utils::write.csv(mrs$df, file, row.names = FALSE))
    output$dl_table_xlsx <- downloadHandler(function() paste0("cross_omics_mr_", Sys.Date(), ".xlsx"), function(file) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) { utils::write.csv(mrs$df, file, row.names = FALSE); return() }
      openxlsx::write.xlsx(mrs$df, file)
    })
    output$dl_tier_csv <- downloadHandler(function() paste0("cross_omics_mr_tiers_", Sys.Date(), ".csv"), function(file) utils::write.csv(tiers(), file, row.names = FALSE))

    observe({
      if (is.null(mrs$df) || is.null(cross_results)) return()
      cross_results$mrstage <- list(df = mrs$df, tiers = tiers(), filters = tier_filters(), run_at = mrs$loaded_at)
    })
  })
}
