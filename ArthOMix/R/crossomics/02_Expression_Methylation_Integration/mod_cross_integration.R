## R/crossomics/02_Expression_Methylation_Integration/mod_cross_integration.R
## Cross-Omics sub-module: "Expression and Methylation" - integrates the
## Transcriptomics DGE output (gene, log2FC, FDR) and the Methylomics DMP

mod_cross_integration_config <- list(
  id = "integration", title = "Expression and Methylation", icon = "dna", group = "Data",
  description = "Gene-level integration of Transcriptomics differential expression and Methylomics differential methylation."
)

.cx_detect_sex_stratum <- function(expr_source, meth_source) {
  s <- toupper(paste(expr_source %||% "", meth_source %||% ""))
  if (grepl("FEMALE", s)) "female" else if (grepl("MALE", s)) "male" else "all"
}

mod_cross_integration_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("live_source_ui")),
    uiOutput(ns("status_bar")),
    uiOutput(ns("summary_cards")),
    tabsetPanel(
      id = ns("result_tabs"), type = "tabs",
      tabPanel("Expression data", br(), uiOutput(ns("expr_data_ui"))),
      tabPanel("Methylation data", br(), uiOutput(ns("meth_data_ui"))),
      tabPanel("Integration", br(), fluidRow(
        column(
          7,
          h4("Integration Setup"),
          fluidRow(
            column(6, numericInput(ns("expr_thresh"), "Min |log2FC|", value = 0.585, min = 0, step = 0.1)),
            column(6, numericInput(ns("expr_fdr_thresh"), "Max expression FDR", value = 0.05, min = 0, max = 1, step = 0.01))
          ),
          fluidRow(
            column(6, numericInput(ns("meth_thresh"), "Min |Δβ|", value = 0.02, min = 0, max = 1, step = 0.01)),
            column(6, numericInput(ns("meth_fdr_thresh"), "Max methylation FDR", value = 0.05, min = 0, max = 1, step = 0.01))
          ),
          selectInput(ns("agg_method"), "Methylation aggregation (multiple CpGs/gene)",
                      choices = setNames(names(CX_AGGREGATION_METHODS), CX_AGGREGATION_METHODS), selected = "mean"),
          radioButtons(ns("cor_method"), "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"), inline = TRUE),
          radioButtons(ns("padj_method"), "Multiple-testing adjustment", choices = c("Benjamini-Hochberg" = "BH", "Bonferroni" = "bonferroni"), inline = TRUE),
          tags$hr(),
          checkboxInput(ns("show_sig_only"), "Show significant genes only (plots)", value = FALSE),
          checkboxInput(ns("show_labels"), "Show gene labels", value = FALSE),
          checkboxInput(ns("show_nonsig"), "Show non-significant genes", value = TRUE),
          checkboxInput(ns("show_quadrant_lines"), "Show quadrant boundaries", value = TRUE),
          tags$hr(),
          actionButton(ns("run_integration"), "Run Integration", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        column(
          5,
          h4("Advanced Filters"),
          selectizeInput(ns("filter_region"), "Genomic region", choices = CX_REGION_FINE_VOCAB, multiple = TRUE, options = list(placeholder = "Any region")),
          selectizeInput(ns("filter_island"), "CpG island status", choices = CX_ISLAND_VOCAB, multiple = TRUE, options = list(placeholder = "Any island status")),
          selectizeInput(ns("filter_evidence"), "Evidence level", choices = CX_EVIDENCE_LEVELS, multiple = TRUE, options = list(placeholder = "Any evidence level")),
          numericInput(ns("filter_min_cpg"), "Minimum CpG count", value = 0, min = 0, step = 1),
          radioButtons(ns("filter_cor_direction"), "Correlation direction (where computed)", choices = c("Any" = "any", "Positive" = "pos", "Negative" = "neg"), inline = TRUE)
        )
      )),
      tabPanel("Quadrant plot", br(), uiOutput(ns("quadrant_ui"))),
      tabPanel("Heatmap", br(), uiOutput(ns("heatmap_ui"))),
      tabPanel("Network analysis", br(), uiOutput(ns("network_ui"))),
      tabPanel("Export", br(), tagList(
        uiOutput(ns("downloads_ui"))
      ))
    )
  )
}

mod_cross_integration_server <- function(id, cross_dataset, cross_results,
                                          dataset = NULL, results = NULL,
                                          methyl_dataset = NULL, methyl_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw <- reactiveValues(
      expr_df = NULL, expr_source = NULL, expr_sample_cols = character(0), expr_wide = NULL, expr_mapping = NULL,
      meth_df = NULL, meth_source = NULL, meth_sample_cols = character(0), meth_wide = NULL, meth_mapping = NULL,
      meth_unavailable_reason = NULL
    )
    integ <- reactiveValues(df = NULL, params = NULL, pairing = NULL, provenance = NULL, run_at = NULL,
                             universe = NULL, cpg_level = NULL, id_harmonization = NULL, validation = NULL)
    integ_by_sex <- reactiveValues(all = NULL, female = NULL, male = NULL)
    active_filter <- reactiveVal("All")
    network_hint <- reactiveVal(NULL)

    observe({
      raw$expr_df <- cross_dataset$user_expr_df
      raw$expr_source <- cross_dataset$user_expr_source
      raw$expr_wide <- cross_dataset$user_expr_wide
      raw$expr_mapping <- cross_dataset$user_expr_mapping
      raw$expr_sample_cols <- cross_dataset$user_expr_sample_cols %||% character(0)
      raw$meth_df <- cross_dataset$user_meth_df
      raw$meth_source <- cross_dataset$user_meth_source
      raw$meth_wide <- cross_dataset$user_meth_wide
      raw$meth_mapping <- cross_dataset$user_meth_mapping
      raw$meth_sample_cols <- cross_dataset$user_meth_sample_cols %||% character(0)
      raw$meth_unavailable_reason <- NULL
    })

    ## ---- "Use live session results" data source ------------------------
    ## Writes adapted live results/methyl_results into cross_dataset, same as mod_cross_dataset_server's "Use this data" button.
    live_dge_choices <- reactive({
      runs <- (results %||% list())$dge_runs %||% list()
      if (length(runs) == 0) return(NULL)
      stats::setNames(names(runs), vapply(runs, function(r) r$contrast %||% "(unnamed run)", character(1)))
    })

    live_dmp_run <- reactive({
      tbl <- (methyl_results %||% list())$dmp_table
      if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(NULL)
      list(comparison = ((methyl_results %||% list())$dmp %||% list())$comparison %||% "Live Methylomics DMP run", table = tbl)
    })

    output$live_source_ui <- renderUI({
      ch <- live_dge_choices()
      dmp <- live_dmp_run()
      div(
        class = "card", style = "margin-bottom: 10px;",
        div(class = "card-title", icon("bolt"), "Data source"),
        radioButtons(ns("cx_source_mode"), NULL,
                     choices = c("From the Cross-Omics Dataset tab (example/uploaded)" = "dataset_tab",
                                 "Use live Transcriptomics/Methylomics session results" = "live"),
                     selected = "dataset_tab"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'live'", ns("cx_source_mode")),
          if (is.null(ch) && is.null(dmp)) {
            div(class = "empty-note", icon("triangle-exclamation"),
                "Run Differential Expression in Transcriptomics and DMP Analysis in Methylomics first.")
          } else if (is.null(ch)) {
            div(class = "empty-note", icon("triangle-exclamation"), "Run Differential Expression in Transcriptomics first.")
          } else if (is.null(dmp)) {
            div(class = "empty-note", icon("triangle-exclamation"), "Run DMP Analysis in Methylomics first.")
          } else {
            tagList(
              selectInput(ns("live_dge_run"), "Transcriptomics DGE run", choices = ch, selected = unname(utils::tail(ch, 1))),
              p(class = "submodule-desc", sprintf("Methylomics: \"%s\" (latest live DMP run this session).", dmp$comparison)),
              actionButton(ns("use_live_btn"), "Use live results", icon = icon("bolt"), class = "btn-primary btn-sm")
            )
          }
        )
      )
    })

    observeEvent(input$use_live_btn, {
      req(input$live_dge_run)
      runs <- (results %||% list())$dge_runs %||% list()
      dge_run <- runs[[input$live_dge_run]]
      dmp_run <- live_dmp_run()
      if (is.null(dge_run) || is.null(dmp_run)) {
        showNotification("Live Transcriptomics/Methylomics results are not both available.", type = "error")
        return()
      }
      expr_res <- cx_build_live_expr_df(dge_run)
      if (!expr_res$ok) { showNotification(expr_res$error, type = "error"); return() }
      meth_res <- cx_build_live_meth_df(dmp_run)
      if (!meth_res$ok) { showNotification(meth_res$error, type = "error"); return() }

      cross_dataset$user_expr_df <- expr_res$df
      cross_dataset$user_expr_source <- sprintf("Live Transcriptomics DGE run \"%s\"", dge_run$contrast %||% input$live_dge_run)
      cross_dataset$user_expr_wide <- NULL
      cross_dataset$user_expr_mapping <- NULL
      cross_dataset$user_expr_sample_cols <- character(0)

      cross_dataset$user_meth_df <- meth_res$df
      cross_dataset$user_meth_source <- sprintf("Live Methylomics DMP run (%s)", dmp_run$comparison)
      cross_dataset$user_meth_wide <- NULL
      cross_dataset$user_meth_mapping <- NULL
      cross_dataset$user_meth_sample_cols <- character(0)

      showNotification("Live session results loaded - ready to Run Integration.", type = "message")
    }, ignoreInit = TRUE)

    output$status_bar <- renderUI({
      tx_ok <- !is.null(raw$expr_df)
      mx_ok <- !is.null(raw$meth_df)
      map_ok <- tx_ok && mx_ok
      pairing <- integ$pairing
      sample_txt <- if (is.null(pairing)) "Not available" else if (isTRUE(pairing$paired)) sprintf("✓ (%d matched samples)", pairing$n_common) else "Not available (unpaired)"
      integ_ok <- !is.null(integ$df)
      item <- function(label, ok, extra = NULL) {
        span(class = if (isTRUE(ok)) "cx-status-ok" else "cx-status-pending",
             icon(if (isTRUE(ok)) "check" else "circle-notch"), " ", label, if (!is.null(extra)) paste0(" ", extra) else NULL)
      }
      tagList(
        div(class = "card", style = "padding: 10px 14px; margin-bottom: 10px; display:flex; gap:18px; flex-wrap:wrap; font-size: 0.92em;",
            item("Transcriptomics", tx_ok), item("Methylomics", mx_ok),
            item("Gene mapping", map_ok), item("Sample matching", isTRUE(pairing$paired), sample_txt),
            item("Integration", integ_ok)
        ),
        if (!is.null(raw$meth_unavailable_reason)) div(class = "empty-note", icon("triangle-exclamation"), raw$meth_unavailable_reason)
      )
    })

    observeEvent(input$run_integration, withProgress(message = "Running Integration - harmonizing gene identifiers across both panels can take up to a minute for a genome-wide (pooled/ALL) dataset...", value = 0.2, {
      if (is.null(raw$expr_df)) { showNotification("No Transcriptomics data loaded - go to the Cross-Omics \"Dataset\" tab first.", type = "error"); return() }
      if (is.null(raw$meth_df)) { showNotification("No Methylomics data loaded - go to the Cross-Omics \"Dataset\" tab first.", type = "error"); return() }

      id_harm <- cx_harmonize_gene_ids(c(raw$expr_df$gene, unique(raw$meth_df$gene)))
      integ$id_harmonization <- id_harm
      expr_h <- raw$expr_df
      meth_h <- raw$meth_df
      if (isTRUE(id_harm$ok)) {
        expr_h$gene <- cx_apply_harmonization(expr_h$gene, id_harm$df)
        expr_h <- cx_dedup_by_gene(expr_h)
        meth_h$gene <- cx_apply_harmonization(meth_h$gene, id_harm$df)
      }

      agg <- cx_aggregate_methylation(meth_h, method = input$agg_method, meth_thresh = input$meth_thresh, meth_fdr_thresh = input$meth_fdr_thresh)
      if (!agg$ok) { showNotification(agg$error, type = "error"); return() }
      integ$cpg_level <- agg$cpg_level

      expr_j <- data.frame(gene = expr_h$gene, log2fc = expr_h$log2fc,
                            expr_pvalue = expr_h$pvalue, expr_fdr = expr_h$fdr, stringsAsFactors = FALSE)
      meth_j <- agg$df
      colnames(meth_j)[colnames(meth_j) == "pvalue"] <- "meth_pvalue"
      colnames(meth_j)[colnames(meth_j) == "fdr"] <- "meth_fdr"
      joined <- merge(expr_j, meth_j, by = "gene", all = TRUE)

      if (all(is.na(joined$expr_fdr)) && any(!is.na(joined$expr_pvalue))) {
        joined$expr_fdr <- cx_adjust_p(joined$expr_pvalue, input$padj_method)
      }
      if (all(is.na(joined$meth_fdr)) && any(!is.na(joined$meth_pvalue))) {
        joined$meth_fdr <- cx_adjust_p(joined$meth_pvalue, input$padj_method)
      }

      classified <- cx_classify(joined, input$expr_thresh, input$expr_fdr_thresh, input$meth_thresh, input$meth_fdr_thresh)

      pairing <- cx_detect_sample_pairing(raw$expr_sample_cols, raw$meth_sample_cols)
      if (isTRUE(pairing$paired)) {
        expr_mat <- cx_build_gene_sample_matrix(raw$expr_wide, mapping_gene_col(raw$expr_mapping, raw$expr_wide), raw$expr_sample_cols)
        meth_mat <- cx_build_gene_sample_matrix(raw$meth_wide, mapping_gene_col(raw$meth_mapping, raw$meth_wide), raw$meth_sample_cols)
        cor_res <- cx_gene_correlation(expr_mat, meth_mat, pairing$common_samples, method = input$cor_method)
        if (isTRUE(cor_res$ok)) {
          cor_res$df$correlation_fdr <- cx_adjust_p(cor_res$df$p, input$padj_method)
          colnames(cor_res$df) <- c("gene", "correlation_r", "correlation_p", "correlation_n", "correlation_fdr")
          classified <- merge(classified, cor_res$df, by = "gene", all.x = TRUE)
        }
      }
      classified$evidence_level <- cx_classify_evidence(classified, has_correlation = isTRUE(pairing$paired))

      integ$df <- classified
      integ$pairing <- pairing
      integ$validation <- cx_validate_dataset(raw$expr_df, raw$meth_df, id_harm)
      sex_key <- .cx_detect_sex_stratum(raw$expr_source, raw$meth_source)
      integ_by_sex[[sex_key]] <- classified
      integ$universe <- classified$gene
      integ$run_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      integ$params <- list(
        sex_stratum = toupper(sex_key),
        input_mode = if (grepl("^Live Transcriptomics DGE run", raw$expr_source %||% "")) "Live Transcriptomics/Methylomics session results" else "From Dataset tab",
        expr_source = raw$expr_source, meth_source = raw$meth_source,
        expr_thresh = input$expr_thresh, expr_fdr_thresh = input$expr_fdr_thresh,
        meth_thresh = input$meth_thresh, meth_fdr_thresh = input$meth_fdr_thresh,
        agg_method = input$agg_method, cor_method = if (isTRUE(pairing$paired)) input$cor_method else NULL,
        padj_method = input$padj_method,
        sample_matching = if (isTRUE(pairing$paired)) sprintf("Paired (%d common samples)", pairing$n_common) else "Not available (unpaired datasets)",
        gene_annotation_source = if (isTRUE(id_harm$ok)) "org.Hs.eg.db (Bioconductor) - exact ID/alias lookup only, no fuzzy matching" else "Not available - matched on exact provided text only",
        methylation_platform = if (grepl("bacon-adjusted|DMR", raw$meth_source %||% "")) "Illumina 450K (IlluminaHumanMethylation450kanno.ilmn12.hg19)" else "Not specified by uploaded data",
        run_at = integ$run_at
      )
      integ$provenance <- cx_build_provenance(integ$params)
      active_filter("All")

      counts <- table(classified$category)
      cross_results$integration <- list(
        df = classified,
        summary = list(
          n_genes = nrow(classified),
          n_deg = sum(classified$sig_expression, na.rm = TRUE),
          n_dmg = sum(classified$sig_methylation, na.rm = TRUE),
          n_integrated = sum(classified$sig_expression & classified$sig_methylation, na.rm = TRUE),
          counts = as.list(counts)
        ),
        provenance = integ$provenance, params = integ$params, run_at = integ$run_at
      )
      showNotification("Integration complete.", type = "message")
    }))

    mapping_gene_col <- function(mapping, wide_df) {
      if (!is.null(mapping) && !is.na(mapping[["gene"]])) return(mapping[["gene"]])
      "gene"
    }

    output$summary_cards <- renderUI({
      req(integ$df)
      df <- integ$df
      counts <- table(df$category)
      card <- function(label, value, key, color = "blue") {
        div(
          class = "card", style = "flex: 1 1 140px; cursor:pointer; text-align:center; padding:10px;",
          onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("card_click"), key),
          div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink),
              format(value, big.mark = ",")),
          div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label)
        )
      }
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          card("Genes analyzed", nrow(df), "All"),
          card("Significant DEGs", sum(df$sig_expression, na.rm = TRUE), "sig_expr_only", "aqua"),
          card("Significant DMGs", sum(df$sig_methylation, na.rm = TRUE), "sig_meth_only", "aqua"),
          card("Integrated genes", sum(df$sig_expression & df$sig_methylation, na.rm = TRUE), "sig_both", "violet"),
          card("Hyper + Down", counts[["Hyper + Down"]] %||% 0L, "Hyper + Down", "red"),
          card("Hypo + Up", counts[["Hypo + Up"]] %||% 0L, "Hypo + Up", "blue"),
          card("Hyper + Up", counts[["Hyper + Up"]] %||% 0L, "Hyper + Up", "orange"),
          card("Hypo + Down", counts[["Hypo + Down"]] %||% 0L, "Hypo + Down", "violet")
      )
    })
    observeEvent(input$card_click, active_filter(input$card_click))

    filtered_df <- reactive({
      req(integ$df)
      df <- cx_filter_by_category(integ$df, active_filter())
      if (isTRUE(input$show_sig_only)) df <- df[df$sig_expression | df$sig_methylation, , drop = FALSE]
      df <- cx_filter_by_region(df, input$filter_region)
      df <- cx_filter_by_island(df, input$filter_island)
      df <- cx_filter_by_evidence(df, input$filter_evidence)
      df <- cx_filter_by_min_cpg(df, input$filter_min_cpg)
      df <- cx_filter_by_correlation_direction(df, input$filter_cor_direction)
      df
    })

    output$expr_data_ui <- renderUI({
      if (is.null(raw$expr_df) && is.null(integ$df)) {
        return(div(class = "empty-note", icon("circle-info"),
          "No Transcriptomics data loaded yet - go to the Cross-Omics \"Dataset\" tab and load or upload it there. It will appear here automatically."))
      }
      n <- if (!is.null(integ$df)) nrow(integ$df) else nrow(raw$expr_df)
      status <- if (is.null(integ$df)) "not yet integrated - click \"Run Integration\" in the Integration tab" else "genes analyzed"
      tagList(
        p(class = "submodule-desc", sprintf("%s - %s %s.", raw$expr_source %||% "(unknown source)", format(n, big.mark = ","), status)),
        DT::dataTableOutput(ns("expr_data_table"))
      )
    })

    expr_data_display <- reactive({
      if (!is.null(integ$df)) {
        df <- filtered_df()
        cols <- intersect(c("gene", "log2fc", "expr_pvalue", "expr_fdr", "expression_direction", "sig_expression"), colnames(df))
        return(df[, cols, drop = FALSE])
      }
      req(raw$expr_df)
      raw$expr_df
    })

    output$expr_data_table <- DT::renderDataTable({
      req(expr_data_display())
      DT::datatable(expr_data_display(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    output$meth_data_ui <- renderUI({
      if (is.null(raw$meth_df) && is.null(integ$df)) {
        return(div(class = "empty-note", icon("circle-info"),
          "No Methylomics data loaded yet - go to the Cross-Omics \"Dataset\" tab and load or upload it there. It will appear here automatically."))
      }
      n <- if (!is.null(integ$df)) nrow(integ$df) else nrow(raw$meth_df)
      status <- if (is.null(integ$df)) "not yet integrated - click \"Run Integration\" in the Integration tab" else "genes analyzed"
      tagList(
        p(class = "submodule-desc", sprintf("%s - %s %s.", raw$meth_source %||% "(unknown source)", format(n, big.mark = ","), status)),
        DT::dataTableOutput(ns("meth_data_table"))
      )
    })

    meth_data_display <- reactive({
      if (!is.null(integ$df)) {
        df <- filtered_df()
        cols <- intersect(c("gene", "dbeta_mean", "dbeta_median", "meth_pvalue", "meth_fdr", "methylation_direction",
                             "sig_methylation", "primary_region", "region", "n_island",
                             "n_cpg_total", "n_cpg_significant", "cpg"), colnames(df))
        return(df[, cols, drop = FALSE])
      }
      req(raw$meth_df)
      raw$meth_df
    })

    output$meth_data_table <- DT::renderDataTable({
      req(meth_data_display())
      DT::datatable(meth_data_display(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    gene_cpg_rows <- function(gene) {
      if (is.null(integ$cpg_level)) return(NULL)
      integ$cpg_level[integ$cpg_level$gene == gene, , drop = FALSE]
    }

    output$quadrant_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        fluidRow(
          column(6, textInput(ns("quadrant_search"), "Search / label gene(s)", placeholder = "Comma-separated symbols, e.g. TP53, BRCA1")),
          column(6, div(style = "padding-top: 24px;",
                        downloadButton(ns("dl_quadrant_png"), "PNG", class = "btn-sm"),
                        downloadButton(ns("dl_quadrant_pdf"), "PDF", class = "btn-sm"),
                        downloadButton(ns("dl_quadrant_svg"), "SVG", class = "btn-sm")))
        ),
        p(class = "submodule-desc", "Each point is a gene. This shows a statistical association between methylation and expression change - it does not establish that one causes the other."),
        plotly::plotlyOutput(ns("quadrant_plot"), height = "520px")
      )
    })

    quadrant_source_df <- reactive({
      req(integ$df)
      df <- integ$df
      if (isTRUE(input$show_nonsig) == FALSE) df <- df[df$category != "Not significant", , drop = FALSE]
      if (isTRUE(input$show_sig_only)) df <- df[df$sig_expression | df$sig_methylation, , drop = FALSE]
      df[!is.na(df$dbeta) & !is.na(df$log2fc), , drop = FALSE]
    })

    quadrant_search_genes <- reactive({
      g <- trimws(strsplit(input$quadrant_search %||% "", ",")[[1]])
      g[nzchar(g)]
    })

    output$quadrant_plot <- plotly::renderPlotly({
      df <- quadrant_source_df()
      validate(need(nrow(df) > 0, "No genes to plot at the current filters/thresholds."))
      cx_quadrant_plotly(df, input$expr_thresh, input$meth_thresh, input$show_quadrant_lines,
                         highlight = quadrant_search_genes(), show_labels = isTRUE(input$show_labels), source_id = ns("quadrant_src"))
    })

    observeEvent(plotly::event_data("plotly_click", source = ns("quadrant_src")), {
      ed <- plotly::event_data("plotly_click", source = ns("quadrant_src"))
      req(ed$key)
      row <- integ$df[integ$df$gene == ed$key[1], , drop = FALSE]
      validate(need(nrow(row) >= 1, "Gene not found."))
      showModal(cx_gene_detail_modal(row[1, , drop = FALSE], integ$pairing, gene_cpg_rows(row$gene[1])))
    })

    output$heatmap_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        fluidRow(
          column(4, selectInput(ns("heatmap_selection"), "Genes to show",
                                  choices = c("Top N by |log2FC|" = "topn", "All significant" = "sig", "By category" = "category"))),
          column(4, conditionalPanel(condition = sprintf("input['%s'] == 'topn'", ns("heatmap_selection")),
                                       numericInput(ns("heatmap_topn"), "N", value = 30, min = 5, max = 200))),
          column(4, conditionalPanel(condition = sprintf("input['%s'] == 'category'", ns("heatmap_selection")),
                                       selectInput(ns("heatmap_category"), "Category", choices = CX_CATEGORY_ORDER)))
        ),
        fluidRow(
          column(6, checkboxInput(ns("heatmap_cluster"), "Hierarchical clustering", value = TRUE)),
          column(6, radioButtons(ns("heatmap_scale"), "Scaling", choices = c("None" = "none", "Row" = "row", "Column" = "column"), inline = TRUE))
        ),
        downloadButton(ns("dl_heatmap_png"), "PNG", class = "btn-sm"),
        downloadButton(ns("dl_heatmap_pdf"), "PDF", class = "btn-sm"),
        plotOutput(ns("heatmap_plot"), height = "560px")
      )
    })

    heatmap_genes <- reactive({
      req(integ$df)
      df <- integ$df
      df <- switch(input$heatmap_selection %||% "topn",
        topn = df[order(-abs(df$log2fc %||% 0)), , drop = FALSE][seq_len(min(input$heatmap_topn %||% 30, nrow(df))), , drop = FALSE],
        sig = df[df$sig_expression & df$sig_methylation, , drop = FALSE],
        category = df[df$category == input$heatmap_category, , drop = FALSE]
      )
      df
    })

    cx_heatmap_matrix <- function() {
      df <- heatmap_genes()
      validate(need(nrow(df) >= 2, "Not enough genes selected to build a heatmap (need at least 2)."))
      m <- as.matrix(df[, c("log2fc", "dbeta")])
      rownames(m) <- df$gene
      colnames(m) <- c("Expression (log2FC)", "Methylation (Δβ)")
      m
    }

    output$heatmap_plot <- renderPlot({
      m <- cx_heatmap_matrix()
      pheatmap::pheatmap(m, cluster_rows = isTRUE(input$heatmap_cluster), cluster_cols = FALSE,
                          scale = input$heatmap_scale %||% "none",
                          color = grDevices::colorRampPalette(c(ARTHOMIX_COLORS$blue, "white", ARTHOMIX_COLORS$red))(100),
                          fontsize_row = if (nrow(m) > 60) 5 else 8)
    })

    output$network_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        p(class = "submodule-desc", "A static gene-CpG association network (this deployment doesn't have an interactive network graphing package installed)."),
        fluidRow(
          column(6, selectInput(ns("network_category"), "Gene set", choices = c("All significant (both)" = "sig_both", CX_CATEGORY_ORDER[CX_CATEGORY_ORDER != "Not significant"]))),
          column(6, numericInput(ns("network_topn"), "Max genes shown", value = 15, min = 3, max = 40))
        ),
        plotOutput(ns("network_plot"), height = "520px")
      )
    })
    output$network_plot <- renderPlot({
      df <- integ$df
      sub <- if (identical(input$network_category, "sig_both")) df[df$sig_expression & df$sig_methylation, , drop = FALSE] else df[df$category == input$network_category, , drop = FALSE]
      sub <- sub[order(-abs(sub$log2fc %||% 0)), , drop = FALSE]
      sub <- sub[seq_len(min(input$network_topn %||% 15, nrow(sub))), , drop = FALSE]
      validate(need(nrow(sub) >= 2, "Not enough genes in this set to draw a network."))
      cx_gene_cpg_network_plot(sub)
    })

    output$downloads_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      strata_available <- names(Filter(Negate(is.null), list(all = integ_by_sex$all, female = integ_by_sex$female, male = integ_by_sex$male)))
      tagList(
        h5("Gene-level results (current run)"),
        downloadButton(ns("dl_table_csv"), "Results (CSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_tsv"), "Results (TSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_xlsx"), "Results (XLSX)", class = "btn-sm"),
        tags$hr(),
        h5("CpG-level results (Level 2, current run)"),
        downloadButton(ns("dl_cpg_csv"), "CpG-level (CSV)", class = "btn-sm"),
        tags$hr(),
        h5("By analysis group (spec section 28)"),
        if (length(strata_available) == 0) div(class = "empty-note", icon("circle-info"), "Run at least one of ALL/FEMALE/MALE to enable per-group downloads.") else
          tagList(
            if ("all" %in% strata_available) downloadButton(ns("dl_all_stratum_csv"), "All-sample results (CSV)", class = "btn-sm"),
            if ("female" %in% strata_available) downloadButton(ns("dl_female_csv"), "Female results (CSV)", class = "btn-sm"),
            if ("male" %in% strata_available) downloadButton(ns("dl_male_csv"), "Male results (CSV)", class = "btn-sm")
          ),
        tags$hr(),
        h5("Report"),
        downloadButton(ns("dl_report"), "Analysis report (Markdown)", class = "btn-sm"),
        tags$hr(),
        h5("Everything"),
        downloadButton(ns("dl_all"), "Download All Results (.zip)", class = "btn-primary btn-sm")
      )
    })

    output$dl_table_csv <- downloadHandler(paste0("cross_omics_integration_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ$df, file, row.names = FALSE))
    output$dl_table_tsv <- downloadHandler(paste0("cross_omics_integration_", Sys.Date(), ".tsv"), function(file) utils::write.table(integ$df, file, row.names = FALSE, sep = "\t", quote = FALSE))
    output$dl_table_xlsx <- downloadHandler(paste0("cross_omics_integration_", Sys.Date(), ".xlsx"), function(file) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) { utils::write.csv(integ$df, file, row.names = FALSE); return() }
      openxlsx::write.xlsx(integ$df, file)
    })
    output$dl_cpg_csv <- downloadHandler(paste0("cross_omics_cpg_level_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ$cpg_level, file, row.names = FALSE))
    output$dl_all_stratum_csv <- downloadHandler(paste0("cross_omics_ALL_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ_by_sex$all, file, row.names = FALSE))
    output$dl_female_csv <- downloadHandler(paste0("cross_omics_FEMALE_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ_by_sex$female, file, row.names = FALSE))
    output$dl_male_csv <- downloadHandler(paste0("cross_omics_MALE_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ_by_sex$male, file, row.names = FALSE))

    output$dl_quadrant_png <- downloadHandler(paste0("quadrant_plot_", Sys.Date(), ".png"), function(file) ggplot2::ggsave(file, cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6, dpi = 300))
    output$dl_quadrant_pdf <- downloadHandler(paste0("quadrant_plot_", Sys.Date(), ".pdf"), function(file) ggplot2::ggsave(file, cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6))
    output$dl_quadrant_svg <- downloadHandler(paste0("quadrant_plot_", Sys.Date(), ".svg"), function(file) ggplot2::ggsave(file, cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6))

    output$dl_heatmap_png <- downloadHandler(paste0("heatmap_", Sys.Date(), ".png"), function(file) { grDevices::png(file, width = 7, height = 7, units = "in", res = 300); pheatmap::pheatmap(cx_heatmap_matrix(), cluster_rows = isTRUE(input$heatmap_cluster), cluster_cols = FALSE, scale = input$heatmap_scale %||% "none"); grDevices::dev.off() })
    output$dl_heatmap_pdf <- downloadHandler(paste0("heatmap_", Sys.Date(), ".pdf"), function(file) { grDevices::pdf(file, width = 7, height = 7); pheatmap::pheatmap(cx_heatmap_matrix(), cluster_rows = isTRUE(input$heatmap_cluster), cluster_cols = FALSE, scale = input$heatmap_scale %||% "none"); grDevices::dev.off() })

    output$dl_report <- downloadHandler(paste0("cross_omics_report_", Sys.Date(), ".md"), function(file) writeLines(cx_build_report(integ$df, integ$provenance), file))

    output$dl_all <- downloadHandler(
      filename = paste0("cross_omics_all_results_", Sys.Date(), ".zip"),
      content = function(file) {
        tmp <- tempfile(); dir.create(tmp)
        utils::write.csv(integ$df, file.path(tmp, "results.csv"), row.names = FALSE)
        writeLines(cx_build_report(integ$df, integ$provenance), file.path(tmp, "report.md"))
        tryCatch(ggplot2::ggsave(file.path(tmp, "quadrant_plot.png"), cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6, dpi = 300), error = function(e) NULL)
        old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
        setwd(tmp)
        utils::zip(file, files = list.files(tmp))
      }
    )
  })
}
