## R/crossomics/03_Biomarker_Convergence/mod_cross_biomarker_conv.R
## Submodule: Biomarker Convergence - loads the pipeline's own already-joined
## eQTL-MR x mQTL-MR x DEG x DMP x DMR table (cross_Omics_Sexstratified_COPY/
## results/cross_omics_eQTL_mQTL_{sex}.csv - see crossomics_biomarkerconv_
## helpers.R) and lets the significance thresholds that decide "which rows
## count as evidence" be reconfigured instantly, by relabeling the table's
## own retained raw p/FDR/effect-size values. No join is re-run, no
## statistics are recomputed, and no file outside cross_Omics_Sexstratified_
## COPY is ever read - this module integrates already-computed
## Transcriptomics/Methylomics potential-biomarker evidence, it does not
## reproduce the analyses that produced it.
##
## Language throughout: "evidence present in N of 3 layers", never
## "confirmed biomarker" - joining/relabeling never upgrades association
## evidence into a causal or diagnostic claim.

mod_cross_biomarker_conv_config <- list(
  id = "biomarkerconv", title = "Biomarker Convergence", icon = "diagram-project", group = "Data",
  description = "Loads and relabels precomputed eQTL-MR and mQTL-MR evidence at configurable significance thresholds."
)

mod_cross_biomarker_conv_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Cohort", status = "primary", solidHeader = FALSE,
        radioButtons(ns("data_source"), "Data source",
                     choices = c("Preloaded (this project's pipeline)" = "preloaded", "Upload your own data" = "upload"),
                     selected = "preloaded"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
          radioButtons(ns("sex"), "Sex", choices = c("Female" = "female", "Male" = "male", "Both (female and male)" = "combined"), selected = "female", inline = TRUE),
          p(class = "empty-note", icon("circle-info"),
            "Already-joined eQTL-MR, mQTL-MR, DEG, DMP, and DMR results, loaded as one table."),
          actionButton(ns("load_table"), "Load Table", icon = icon("database"), class = "btn-primary btn-sm", width = "100%")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
          fileInput(ns("upload_eqtl_file"), "eQTL-MR results (optional)",
                    accept = c(".csv", ".tsv", ".txt", ".xlsx"), placeholder = "CSV / TSV / TXT / XLSX"),
          p(class = "empty-note", icon("circle-info"),
            "One row per gene. Required: gene. Optional: eQTL_MR_OR, eQTL_MR_pval, eQTL_MR_FDR, eQTL_MHC_region."),
          fileInput(ns("upload_mqtl_file"), "mQTL-MR results (optional)",
                    accept = c(".csv", ".tsv", ".txt", ".xlsx"), placeholder = "CSV / TSV / TXT / XLSX"),
          p(class = "empty-note", icon("circle-info"),
            "One row per gene. Required: gene, mQTL_MR_pval. Optional: mQTL_candidate_cpg, mQTL_cpg_chr, mQTL_cpg_pos_hg19, mQTL_MR_beta."),
          p(class = "empty-note", icon("triangle-exclamation"),
            "Upload at least one of the two files - the eQTL-MR and mQTL-MR tabs each work from their own file alone. But the eQTL-mQTL tab (genes with evidence in BOTH) needs both files uploaded - with only one, it will always be empty, since there's nothing to intersect it against. The two files are merged live by gene here (a plain outer join, not a re-run of either MR analysis). DEG/DMP/DMR significance isn't available from this path - this module never re-derives those; use the Cross-Omics Dataset tab / Expression and Methylation module for that."),
          actionButton(ns("load_table_upload"), "Merge & Load", icon = icon("upload"), class = "btn-primary btn-sm", width = "100%")
        )
      )
    ),
    column(
      8,
      uiOutput(ns("summary_cards")),
      tabsetPanel(
        id = ns("result_tabs"), type = "tabs",
        tabPanel("eQTL-MR", br(), uiOutput(ns("eqtl_tab_ui"))),
        tabPanel("mQTL-MR", br(), uiOutput(ns("mqtl_tab_ui"))),
        tabPanel("eQTL-mQTL", br(), uiOutput(ns("eqtl_mqtl_tab_ui"))),
        tabPanel("Downloads", br(), uiOutput(ns("downloads_ui")))
      )
    )
  )
}

mod_cross_biomarker_conv_server <- function(id, cross_dataset, cross_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw <- reactiveValues(df = NULL, sex = NULL, loaded_at = NULL, missing_layer = NULL)

    ## Preloaded and Upload are independent pipelines - switching the radio
    ## clears whatever the other mode had loaded, so stale preloaded data
    ## can't keep showing under "Upload your own data" (or vice versa) just
    ## because its own Load button hasn't been clicked yet in this mode.
    observeEvent(input$data_source, {
      raw$df <- NULL; raw$sex <- NULL; raw$loaded_at <- NULL; raw$missing_layer <- NULL
    }, ignoreInit = TRUE)

    observeEvent(input$load_table, {
      res <- cx_bc_load_precomputed(input$sex)
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      raw$df <- res$df
      raw$sex <- input$sex
      raw$loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      raw$missing_layer <- NULL
      showNotification(sprintf("Loaded %s genes (%s).", format(nrow(res$df), big.mark = ","), toupper(input$sex)), type = "message")
    })

    observeEvent(input$load_table_upload, {
      validate(need(!is.null(input$upload_eqtl_file) || !is.null(input$upload_mqtl_file),
                    "Upload at least one file (eQTL-MR and/or mQTL-MR)."))
      eqtl_res <- if (!is.null(input$upload_eqtl_file)) cx_bc_load_eqtl_upload(input$upload_eqtl_file$datapath, input$upload_eqtl_file$name) else NULL
      if (!is.null(eqtl_res) && !eqtl_res$ok) { showNotification(eqtl_res$error, type = "error", duration = 12); return() }
      mqtl_res <- if (!is.null(input$upload_mqtl_file)) cx_bc_load_mqtl_upload(input$upload_mqtl_file$datapath, input$upload_mqtl_file$name) else NULL
      if (!is.null(mqtl_res) && !mqtl_res$ok) { showNotification(mqtl_res$error, type = "error", duration = 12); return() }

      merged <- cx_bc_merge_eqtl_mqtl(if (!is.null(eqtl_res)) eqtl_res$df else NULL, if (!is.null(mqtl_res)) mqtl_res$df else NULL)
      if (!merged$ok) { showNotification(merged$error, type = "error"); return() }

      raw$df <- merged$df
      raw$sex <- "uploaded"
      raw$loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      files_desc <- paste(c(
        if (!is.null(input$upload_eqtl_file)) sprintf("eQTL-MR: \"%s\"", input$upload_eqtl_file$name),
        if (!is.null(input$upload_mqtl_file)) sprintf("mQTL-MR: \"%s\"", input$upload_mqtl_file$name)
      ), collapse = "; ")
      ## Only one of the two files given - flagged so the eQTL-mQTL tab can
      ## explain WHY it's empty (nothing to intersect against), rather than
      ## looking like the same "these panels just don't overlap" case the
      ## preloaded data can genuinely show.
      raw$missing_layer <- if (is.null(input$upload_eqtl_file)) "eQTL-MR" else if (is.null(input$upload_mqtl_file)) "mQTL-MR" else NULL
      showNotification(sprintf("Loaded %s genes from %s.", format(nrow(merged$df), big.mark = ","), files_desc), type = "message")
      if (!is.null(raw$missing_layer)) {
        showNotification(sprintf("You only uploaded %s data - the eQTL-mQTL tab needs both files to find genes with evidence in both, so it will show 0. Upload %s results too if you want that tab populated.",
                                  if (raw$missing_layer == "eQTL-MR") "mQTL-MR" else "eQTL-MR", raw$missing_layer),
                          type = "warning", duration = 15)
      }
    })

    ## No filter UI left in this module - every significance threshold
    ## (DEG/DMP/DMR FDR, mQTL-MR basis/cutoff) is fixed at the pipeline's
    ## own default (CX_BC_DEFAULT_PARAMS). Relabeling itself still runs
    ## through cx_bc_relabel() so the *_significant/n_evidence_layers
    ## columns published to other Cross-Omics sub-modules are unchanged.
    bc_df <- reactive({
      req(raw$df)
      cx_bc_relabel(raw$df, CX_BC_DEFAULT_PARAMS)
    })

    output$summary_cards <- renderUI({
      df <- tryCatch(bc_df(), error = function(e) NULL)
      if (is.null(df)) return(NULL)
      card <- function(label, value, color = "blue") {
        div(class = "card", style = "flex: 1 1 150px; text-align:center; padding:10px;",
            div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink), format(value, big.mark = ",")),
            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label))
      }
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          card("Genes (union)", nrow(df)),
          card("In eQTL-MR panel", sum(df$in_eQTL_MR_panel), "aqua"),
          card("In mQTL-MR panel", sum(df$in_mQTL_MR_panel), "aqua"),
          card("In both (eQTL-mQTL)", sum(df$in_eQTL_MR_panel %in% TRUE & df$in_mQTL_MR_panel %in% TRUE), "red"),
          card("DEG significant", sum(df$DEG_significant, na.rm = TRUE), "violet"),
          card("Methylation significant (DMP or DMR)", sum(df$methylation_significant, na.rm = TRUE), "violet")
      )
    })

    ## ---- eQTL-MR tab: genes in the eQTL-MR panel, that layer's own columns ----
    eqtl_df <- reactive({
      df <- bc_df()
      d <- df[df$in_eQTL_MR_panel %in% TRUE, , drop = FALSE]
      cols <- intersect(c("gene", "eQTL_MR_OR", "eQTL_MR_pval", "eQTL_MR_FDR", "eQTL_MHC_region", "eQTL_MR_significant"), colnames(d))
      d[, cols, drop = FALSE]
    })
    ## Shared by the eQTL-MR and eQTL-mQTL tabs below: how many rows in a table that
    ## carries eQTL_MHC_region sit in the MHC region (chr6, ~25-34Mb) - the single most
    ## notorious horizontal-pleiotropy hotspot in autoimmune-disease genetics, where an
    ## MR estimate is more likely to reflect regional LD than the gene's own causal
    ## effect. NULL (renders nothing) when there's nothing to warn about.
    cx_bc_mhc_warning <- function(d) {
      if (is.null(d) || !"eQTL_MHC_region" %in% colnames(d)) return(NULL)
      n_mhc <- sum(d$eQTL_MHC_region %in% TRUE)
      if (n_mhc == 0) return(NULL)
      p(class = "empty-note", icon("triangle-exclamation"), style = "border-color: var(--color-warning, #e0a800);",
        sprintf("%d of %d row(s) here fall in the MHC region - treat these as considerably less reliable than non-MHC hits, regardless of how significant they look.", n_mhc, nrow(d)))
    }

    output$eqtl_tab_ui <- renderUI({
      df <- tryCatch(bc_df(), error = function(e) NULL)
      if (is.null(df)) return(cx_empty_state("Load a table on the \"1. Cohort\" panel to see results here."))
      tagList(
        p(class = "submodule-desc", sprintf("%s of %s genes have an eQTL-MR causal-expression instrument in this panel.", format(sum(df$in_eQTL_MR_panel %in% TRUE), big.mark = ","), format(nrow(df), big.mark = ","))),
        cx_bc_mhc_warning(eqtl_df()),
        div(class = "table-toolbar", downloadButton(ns("dl_eqtl_csv"), "CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("eqtl_table"))
      )
    })
    output$eqtl_table <- DT::renderDataTable({
      req(eqtl_df())
      DT::datatable(eqtl_df(), rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })
    output$dl_eqtl_csv <- downloadHandler(function() paste0("biomarker_convergence_eQTL-MR_", raw$sex, "_", Sys.Date(), ".csv"), function(file) utils::write.csv(eqtl_df(), file, row.names = FALSE))

    ## ---- mQTL-MR tab: genes in the mQTL-MR panel, that layer's own columns ----
    mqtl_df <- reactive({
      df <- bc_df()
      d <- df[df$in_mQTL_MR_panel %in% TRUE, , drop = FALSE]
      cols <- intersect(c("gene", "mQTL_candidate_cpg", "mQTL_cpg_chr", "mQTL_cpg_pos_hg19",
                           "mQTL_instrument_available", "mQTL_MR_beta", "mQTL_MR_pval", "mQTL_MR_significant"), colnames(d))
      d[, cols, drop = FALSE]
    })
    output$mqtl_tab_ui <- renderUI({
      df <- tryCatch(bc_df(), error = function(e) NULL)
      if (is.null(df)) return(cx_empty_state("Load a table on the \"1. Cohort\" panel to see results here."))
      tagList(
        p(class = "submodule-desc", sprintf("%s of %s genes have an mQTL-MR causal-methylation instrument in this panel.", format(sum(df$in_mQTL_MR_panel %in% TRUE), big.mark = ","), format(nrow(df), big.mark = ","))),
        div(class = "table-toolbar", downloadButton(ns("dl_mqtl_csv"), "CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("mqtl_table"))
      )
    })
    output$mqtl_table <- DT::renderDataTable({
      req(mqtl_df())
      DT::datatable(mqtl_df(), rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })
    output$dl_mqtl_csv <- downloadHandler(function() paste0("biomarker_convergence_mQTL-MR_", raw$sex, "_", Sys.Date(), ".csv"), function(file) utils::write.csv(mqtl_df(), file, row.names = FALSE))

    ## ---- eQTL-mQTL tab: genes with BOTH a causal-expression AND a causal- ----
    ## methylation instrument - the genetic convergence set, a stronger
    ## candidate list than either MR panel alone.
    eqtl_mqtl_df <- reactive({
      df <- bc_df()
      d <- df[df$in_eQTL_MR_panel %in% TRUE & df$in_mQTL_MR_panel %in% TRUE, , drop = FALSE]
      cols <- intersect(c("gene", "eQTL_MR_OR", "eQTL_MR_pval", "eQTL_MR_FDR", "eQTL_MHC_region", "eQTL_MR_significant",
                           "mQTL_candidate_cpg", "mQTL_MR_beta", "mQTL_MR_pval", "mQTL_MR_significant"), colnames(d))
      d[, cols, drop = FALSE]
    })
    output$eqtl_mqtl_tab_ui <- renderUI({
      df <- tryCatch(bc_df(), error = function(e) NULL)
      if (is.null(df)) return(cx_empty_state("Load a table on the \"1. Cohort\" panel to see results here."))
      n_both <- sum(df$in_eQTL_MR_panel %in% TRUE & df$in_mQTL_MR_panel %in% TRUE)
      tagList(
        if (!is.null(raw$missing_layer)) p(class = "empty-note", icon("triangle-exclamation"), style = "border-color: var(--color-warning, #e0a800);",
          sprintf("This tab is empty because you only uploaded %s data - there's nothing to intersect it against. Upload %s results too (on the \"1. Cohort\" panel) and click \"Merge & Load\" again to see genes with evidence in both.",
                  if (identical(raw$missing_layer, "eQTL-MR")) "mQTL-MR" else "eQTL-MR", raw$missing_layer))
        else if (n_both == 0) p(class = "empty-note", icon("circle-info"),
          "0 here doesn't necessarily mean something is wrong: both files/panels are loaded, but this loaded data's eQTL-MR and mQTL-MR gene sets simply don't share any genes (this is expected for independently-uploaded files with unrelated candidate gene lists; the preloaded pipeline data always has some overlap, once genes with real mQTL-MR evidence that this table's own join had dropped are backfilled)."),
        p(class = "submodule-desc", sprintf("%s genes have BOTH an eQTL-MR and an mQTL-MR instrument - genetic evidence for a causal effect on both expression and methylation, independent of the DEG/DMP/DMR observational layers.", format(n_both, big.mark = ","))),
        cx_bc_mhc_warning(eqtl_mqtl_df()),
        div(class = "table-toolbar", downloadButton(ns("dl_eqtl_mqtl_csv"), "CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("eqtl_mqtl_table"))
      )
    })
    output$eqtl_mqtl_table <- DT::renderDataTable({
      req(eqtl_mqtl_df())
      DT::datatable(eqtl_mqtl_df(), rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })
    output$dl_eqtl_mqtl_csv <- downloadHandler(function() paste0("biomarker_convergence_eQTL-mQTL_", raw$sex, "_", Sys.Date(), ".csv"), function(file) utils::write.csv(eqtl_mqtl_df(), file, row.names = FALSE))

    output$downloads_ui <- renderUI({
      if (is.null(raw$df)) return(cx_empty_state("Load a table on the \"1. Cohort\" panel to see results here."))
      tagList(
        downloadButton(ns("dl_table_csv"), "Results (CSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_xlsx"), "Results (XLSX)", class = "btn-sm")
      )
    })
    output$dl_table_csv <- downloadHandler(function() paste0("biomarker_convergence_", raw$sex, "_", Sys.Date(), ".csv"), function(file) utils::write.csv(bc_df(), file, row.names = FALSE))
    output$dl_table_xlsx <- downloadHandler(function() paste0("biomarker_convergence_", raw$sex, "_", Sys.Date(), ".xlsx"), function(file) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) { utils::write.csv(bc_df(), file, row.names = FALSE); return() }
      openxlsx::write.xlsx(bc_df(), file)
    })

    ## Published for other Cross-Omics sub-modules (Cross-Omics MR's Tier
    ## classification) and ArthOChat to read - kept live/reactive to filter
    ## changes, not just the initial load.
    observe({
      df <- tryCatch(bc_df(), error = function(e) NULL)
      if (is.null(df) || is.null(cross_results)) return()
      cross_results$biomarkerconv <- list(df = df, sex = raw$sex, run_at = raw$loaded_at)
    })
  })
}
