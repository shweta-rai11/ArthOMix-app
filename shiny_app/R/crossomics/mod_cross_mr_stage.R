## R/crossomics/mod_cross_mr_stage.R
## Submodule: Cross-Omics MR - a live, reconfigurable reimplementation of
## cross_Omics_Sexstratified_COPY/scripts/02_mr_stage_cross_omics.R:
## single-instrument mQTL-MR (Wald ratio, lead cis-SNP) for the eQTL-MR-
## significant genes, bridging the eQTL-MR gene panel and mQTL-MR CpG panel
## that share zero genes by design (see crossomics_mrstage_helpers.R for the
## MR logic itself). Also computes the Tier 1/2/3 classification described in
## the pipeline's own report but never implemented in code by either script -
## see crossomics_mrstage_helpers.R's CX_MR_TIER_RULE_TEXT for the exact rule,
## shown verbatim in this module's UI.
##
## Scope boundary (also shown in the UI): the instrument panel is the
## pipeline's own pre-extracted GoDMC lead-SNP set - this module never reads
## the raw ~6.3GB GoDMC association file live, so it cannot add arbitrary
## genes outside that existing panel.

mod_cross_mr_stage_config <- list(
  id = "mrstage", title = "Cross-Omics MR", icon = "arrow-right-arrow-left", group = "Genetics",
  description = "Single-instrument mQTL-MR for the eQTL-MR-significant genes (GoDMC exposure, Ishigaki 2022 RA outcome), bridging the two panels that share zero genes by design."
)

mod_cross_mr_stage_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Instrument Panel", status = "primary", solidHeader = FALSE,
        p(class = "empty-note", icon("circle-info"),
          "Fixed panel: one lead cis-mQTL SNP per CpG, pre-extracted from GoDMC for the genes already significant in the eQTL-MR panel. This module cannot add arbitrary genes - doing so would require re-extracting instruments live from the raw ~6.3GB genome-wide GoDMC association file, which it never reads."),
        uiOutput(ns("panel_info_ui"))
      ),
      box(
        width = NULL, title = "2. MR Setup", status = "primary", solidHeader = FALSE,
        numericInput(ns("min_f_stat"), "Minimum instrument F-statistic", value = 10, min = 0, step = 1),
        radioButtons(ns("harmonise_action"), "Harmonisation action (TwoSampleMR)",
                     choices = c("2 - drop palindromic SNPs (matches original script; outcome file has no EAF)" = 2, "1 - assume all alleles on the forward strand" = 1, "3 - infer strand from allele frequency" = 3),
                     selected = 2),
        radioButtons(ns("steiger_mode"), "Steiger-filtering handling",
                     choices = c("Flag only (matches original script)" = "flag", "Drop instruments failing the direction test" = "drop"),
                     selected = "flag"),
        tags$hr(),
        actionButton(ns("run_mr"), "Run MR", icon = icon("play"), class = "btn-primary btn-sm", width = "100%"),
        p(class = "submodule-desc", style = "margin-top:8px;", "Reads a filtered subset of the outcome GWAS file on every run - typically takes about a minute.")
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
        tabPanel("Tier Classification", br(), uiOutput(ns("tier_ui"))),
        tabPanel("Volcano", br(), uiOutput(ns("volcano_ui"))),
        tabPanel("Downloads", br(), uiOutput(ns("downloads_ui")))
      )
    )
  )
}

mod_cross_mr_stage_server <- function(id, cross_dataset, cross_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    mrs <- reactiveValues(df = NULL, tiers = NULL, params = NULL, run_at = NULL, provenance = NULL,
                           n_tested = NULL, n_harmonised = NULL, join_used_sex = NULL)

    output$panel_info_ui <- renderUI({
      if (!CX_MR_DATA_AVAILABLE) return(div(class = "empty-note", icon("triangle-exclamation"), "Cross-Omics MR source data is not available in this deployment."))
      inst <- cx_mr_load_instruments(min_f_stat = input$min_f_stat %||% 10)
      if (!inst$ok) return(div(class = "empty-note", icon("triangle-exclamation"), inst$error))
      div(class = "empty-note", icon("dna"), sprintf(
        "%s instruments with a resolvable rsID; %s pass the current F-statistic threshold, covering %s genes.",
        format(inst$n_before_fstat, big.mark = ","), format(inst$n_after_fstat, big.mark = ","), format(length(unique(inst$df$gene)), big.mark = ",")))
    })

    observeEvent(input$run_mr, {
      if (!CX_MR_DATA_AVAILABLE) { showNotification("Cross-Omics MR source data is not available in this deployment.", type = "error"); return() }
      params <- list(min_f_stat = input$min_f_stat, harmonise_action = as.integer(input$harmonise_action), steiger_mode = input$steiger_mode)
      withProgress(message = "Running Cross-Omics MR...", detail = "Loading instruments...", value = 0.1, {
        inst <- cx_mr_load_instruments(params$min_f_stat)
        if (!inst$ok) { showNotification(inst$error, type = "error"); return() }
        incProgress(0.2, detail = "Reading outcome GWAS (this is the slow step, ~1 min)...")
        res <- cx_mr_run_wald_ratio(inst$df, params)
        if (!res$ok) { showNotification(res$error, type = "error"); return() }
        incProgress(0.6, detail = "Classifying evidence tiers...")

        ## Tier classification needs DEG/methylation evidence per gene - reuse
        ## a completed Biomarker Convergence run for the majority sex among
        ## the instrument genes if one is already available (fast path);
        ## otherwise compute it directly via the shared helpers (no
        ## dependency on that tab having been opened).
        join_sex <- if (!is.null(cross_results) && !is.null(cross_results$biomarkerconv)) cross_results$biomarkerconv$sex else "female"
        join_df <- NULL
        if (!is.null(cross_results) && !is.null(cross_results$biomarkerconv) && identical(cross_results$biomarkerconv$sex, join_sex)) {
          join_df <- cross_results$biomarkerconv$df
        } else {
          j <- cx_bc_build_join(join_sex)
          if (isTRUE(j$ok)) join_df <- j$df
        }
        mrs$join_used_sex <- if (!is.null(join_df)) join_sex else NULL
        mrs$tiers <- cx_mr_classify_tier(res$df, join_df)

        mrs$df <- res$df
        mrs$params <- params
        mrs$n_tested <- res$n_tested
        mrs$n_harmonised <- res$n_harmonised
        mrs$run_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        mrs$provenance <- cx_mr_build_provenance(params, res$n_tested, res$n_harmonised, mrs$run_at)
        if (!is.null(cross_results)) {
          cross_results$mrstage <- list(df = res$df, tiers = mrs$tiers, params = params, run_at = mrs$run_at)
        }
        incProgress(1)
      })
      if (!is.null(mrs$df)) showNotification(sprintf("MR complete - %s instruments tested, %s significant at FDR<0.05.", format(mrs$n_tested, big.mark = ","), sum(mrs$df$mr_significant, na.rm = TRUE)), type = "message")
    })

    output$provenance_ui <- renderUI({
      if (is.null(mrs$provenance)) return(div(class = "empty-note", icon("circle-info"), "Run MR to see its parameters here."))
      tagList(
        tags$ul(style = "padding-left: 18px; font-size: 0.85em;", lapply(mrs$provenance, tags$li)),
        if (!is.null(mrs$join_used_sex)) p(class = "submodule-desc", sprintf("Tier classification used the %s Biomarker Convergence join for DEG/methylation evidence.", toupper(mrs$join_used_sex)))
      )
    })

    output$summary_cards <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      df <- mrs$df; t <- mrs$tiers
      card <- function(label, value, color = "blue") {
        div(class = "card", style = "flex: 1 1 150px; text-align:center; padding:10px;",
            div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink), format(value, big.mark = ",")),
            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label))
      }
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          card("Instruments tested", mrs$n_tested), card("Harmonised", mrs$n_harmonised, "aqua"),
          card("Significant (FDR<0.05)", sum(df$mr_significant, na.rm = TRUE), "red"),
          card("Tier 1", if (!is.null(t)) sum(t$tier == "Tier 1", na.rm = TRUE) else 0, "violet"),
          card("Tier 2", if (!is.null(t)) sum(t$tier == "Tier 2", na.rm = TRUE) else 0, "orange"),
          card("Tier 3", if (!is.null(t)) sum(t$tier == "Tier 3", na.rm = TRUE) else 0, "ink_muted")
      )
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

    output$tier_ui <- renderUI({
      if (is.null(mrs$tiers)) return(cx_empty_state())
      tagList(
        h5("How Tier is computed"),
        tags$ul(style = "padding-left: 18px; font-size: 0.88em;", lapply(CX_MR_TIER_RULE_TEXT, tags$li)),
        if (!isTRUE(mrs$tiers$tier_evidence_available[1])) div(class = "empty-note", icon("triangle-exclamation"), "DEG/methylation evidence was not available - every gene defaults to Tier 3 (no gene can be shown as Tier 1/2 without it)."),
        DT::dataTableOutput(ns("tier_table"))
      )
    })
    output$tier_table <- DT::renderDataTable({
      req(mrs$tiers)
      cols <- intersect(c("gene", "DEG_significant", "methylation_significant", "credible_mQTL_MR", "tier"), colnames(mrs$tiers))
      DT::datatable(mrs$tiers[order(mrs$tiers$tier), cols, drop = FALSE], rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    output$volcano_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state())
      plotOutput(ns("volcano_plot"), height = "460px")
    })
    output$volcano_plot <- renderPlot({
      df <- mrs$df
      validate(need(!is.null(df) && nrow(df) > 0, "Run MR first."))
      d <- df[!is.na(df$b) & !is.na(df$FDR), , drop = FALSE]
      d$neglog10fdr <- -log10(pmax(d$FDR, 1e-300))
      ggplot2::ggplot(d, ggplot2::aes(x = b, y = neglog10fdr, color = mr_significant)) +
        ggplot2::geom_point(alpha = 0.75, size = 1.8) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
        ggrepel::geom_text_repel(data = d[d$mr_significant %in% TRUE, , drop = FALSE], ggplot2::aes(label = gene), color = ARTHOMIX_COLORS$ink, size = 3, max.overlaps = 30) +
        ggplot2::scale_color_manual(values = c(`TRUE` = ARTHOMIX_COLORS$red, `FALSE` = ARTHOMIX_COLORS$ink_muted), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "Not significant"), name = NULL) +
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
    output$dl_tier_csv <- downloadHandler(function() paste0("cross_omics_mr_tiers_", Sys.Date(), ".csv"), function(file) utils::write.csv(mrs$tiers, file, row.names = FALSE))
  })
}
