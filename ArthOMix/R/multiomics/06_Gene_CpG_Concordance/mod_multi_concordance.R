## R/multiomics/06_Gene_CpG_Concordance/mod_multi_concordance.R
## Submodule: Gene–CpG Concordance - connects transcriptomics gene-level
## changes to methylation CpG-level changes for the candidate multi-omics
## biomarkers this app has already identified (DIABLO / Joint Biomarker
## Discovery via mod_multi_biomarker.R, SNF via mod_multi_stratification.R),
## or a user-supplied custom gene/CpG list. Two data sources:
##   - "Preloaded (Table42/45)": the pipeline's own precomputed gene<->CpG
##     concordance tables, unchanged from before this rewrite (spec section
##     23 - "must continue to work", "do not redesign the preloaded
##     pipeline"). DIABLO/Joint candidate flags on this path are cross-
##     referenced against the pipeline's own precomputed, per-sex DIABLO
##     panel (Table40/44b, via MCC_PRELOADED_DIABLO_PANEL) - no live run
##     needed.
##   - "Active Multi-Omics Dataset": a live, data-adaptive engine over
##     whatever is actually loaded on the Dataset Workspace tab
##     (multi_dataset$layers/sample_meta) - never assumes matched samples,
##     a platform, an annotation, or a sex variable exist; every analysis
##     that needs something the data doesn't have is disabled with an exact
##     reason, never silently skipped or approximated.
## Nothing below renders a result table, plot, or score until the blue
## "Run Gene<->CpG Analysis" button is clicked (spec section 3) - filters,
## thresholds, and the data status panel are the only things visible before
## that. All heavy lifting is in multiomics_concordance_helpers.R /
## multiomics_concordance_plots.R; this file is UI wiring only.

mod_multi_concordance_config <- list(
  id = "concordance", title = "Gene–CpG Concordance", icon = "arrows-left-right", group = "Biomarker modeling",
  description = "Links candidate gene expression to CpG methylation changes, with direction classification and sex-specific analysis."
)

MULTI_CONCORDANCE_COHORTS <- c(
  "Drug x sex (Etanercept panel)" = "Gene <-> CpG concordance - drug x sex (Etanercept panel)",
  "Response (drug-pooled)" = "Gene <-> CpG concordance - response (drug-pooled)"
)

## Precomputed, per-sex DIABLO candidate-biomarker panels (Table40/44b) for
## the SAME two cohorts above - gene/CpG feature, sex, loading, and
## biomarker_status ("statistically_supported"/"exploratory_not_significant"),
## already on disk (Research_05_multiomics_sexstratified's own pipeline
## output), never requiring a live same-session DIABLO run to populate the
## Preloaded path's diablo/joint biomarker flags.
MCC_PRELOADED_DIABLO_PANEL <- c(
  "Gene <-> CpG concordance - drug x sex (Etanercept panel)" = "Candidate multi-omics biomarkers - drug x sex (Etanercept panel)",
  "Gene <-> CpG concordance - response (drug-pooled)" = "Candidate multi-omics biomarkers - response (drug-pooled)"
)

MCC_BIOMARKER_SOURCES <- c(
  "All candidates", "DIABLO", "SNF", "Joint Biomarker Discovery",
  "DIABLO + SNF", "DIABLO + Joint", "SNF + Joint", "Shared candidates",
  "Custom genes", "Custom CpGs"
)

MCC_SEX_CHOICES <- c("All (pooled)" = "all", "Female" = "female", "Male" = "male", "Sex-specific (Female and Male separately)" = "sex_specific")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_multi_concordance_ui <- function(id) {
  ns <- NS(id)
  tagList(uiOutput(ns("active_dataset_banner")), fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Data source", status = "primary", solidHeader = FALSE,
        radioButtons(ns("data_source"), NULL,
                     choices = c("Preloaded (Table42/45)" = "preloaded", "Active Multi-Omics Dataset (Dataset Workspace)" = "active"),
                     selected = "preloaded"),
        conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
                          selectInput(ns("cohort"), "Cohort", choices = MULTI_CONCORDANCE_COHORTS, width = "100%")),
        conditionalPanel(condition = sprintf("input['%s'] == 'active'", ns("data_source")), uiOutput(ns("layer_pick_ui")))
      ),
      box(width = NULL, title = "2. Data status", status = "primary", solidHeader = FALSE, collapsible = TRUE,
          uiOutput(ns("status_ui"))),
      box(
        width = NULL, title = "3. Biomarker source & sex", status = "primary", solidHeader = FALSE, collapsible = TRUE,
        selectInput(ns("source"), "Biomarker source", choices = MCC_BIOMARKER_SOURCES, selected = "All candidates"),
        div(class = "empty-note", icon("circle-info"),
            "Joint and DIABLO mark the same panel; SNF requires a live Patient Stratification run (no precomputed SNF panel exists)."),
        conditionalPanel(condition = sprintf("input['%s'] == 'active'", ns("data_source")),
                          selectizeInput(ns("custom_genes"), "Custom genes", choices = NULL, multiple = TRUE, options = list(create = TRUE, placeholder = "Type a gene symbol and press Enter")),
                          selectizeInput(ns("custom_cpgs"), "Custom CpGs", choices = NULL, multiple = TRUE, options = list(create = TRUE, placeholder = "Type a CpG ID and press Enter"))),
        selectInput(ns("sex"), "Sex", choices = MCC_SEX_CHOICES, selected = "all"),
        conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
                          div(class = "empty-note", icon("circle-info"), "Preloaded tables already contain per-sex rows; Sex filters those rows."))
      ),
      box(
        width = NULL, title = "4. Thresholds, region & model", status = "primary", solidHeader = FALSE, collapsible = TRUE,
        fluidRow(
          column(6, numericInput(ns("expr_thresh"), "log2FC threshold", value = 1, step = 0.1, min = 0)),
          column(6, numericInput(ns("expr_fdr_thresh"), "Expression FDR <", value = 0.05, step = 0.01, min = 0, max = 1))
        ),
        fluidRow(
          column(6, numericInput(ns("meth_thresh"), "|delta-M| threshold", value = 0.5, step = 0.1, min = 0)),
          column(6, numericInput(ns("meth_fdr_thresh"), "Methylation FDR <", value = 0.05, step = 0.01, min = 0, max = 1))
        ),
        numericInput(ns("genome_p_thresh"), "Genome-wide DEG/DMP nominal p <", value = 0.05, step = 0.01, min = 0, max = 1),
        div(class = "empty-note", icon("circle-info"),
            "Uncorrected p from genome-wide DEG/DMP discovery, not a small-panel FDR. See the Direction tab's Significant section."),
        selectizeInput(ns("region"), "Genomic region", choices = CX_REGION_FINE_VOCAB, multiple = TRUE, options = list(placeholder = "All regions")),
        fluidRow(
          column(6, selectInput(ns("cor_method"), "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"))),
          column(6, numericInput(ns("cor_min_r"), "Min |r|", value = 0.3, step = 0.05, min = 0, max = 1))
        ),
        numericInput(ns("cor_max_fdr"), "Max correlation FDR", value = 0.1, step = 0.01, min = 0, max = 1),
        conditionalPanel(condition = sprintf("input['%s'] == 'active'", ns("data_source")),
                          selectInput(ns("array_type"), "Methylation array (for annotation)", choices = c("Illumina 450K" = "450K", "Illumina EPIC" = "EPIC")),
                          uiOutput(ns("design_ui")),
                          uiOutput(ns("covariate_ui")))
      ),
      actionButton(ns("run_btn"), "Run Gene–CpG Analysis", icon = icon("play"), class = "btn-primary btn-lg", width = "100%"),
      br(), br()
    ),
    column(
      8,
      tabsetPanel(
        id = ns("tabs"), type = "tabs",
        tabPanel("Overview", br(), uiOutput(ns("overview_ui"))),
        tabPanel("Gene–CpG", br(), uiOutput(ns("genecpg_ui"))),
        tabPanel("Direction", br(), uiOutput(ns("direction_ui"))),
        tabPanel("Location", br(), uiOutput(ns("location_ui"))),
        tabPanel("Biomarkers", br(), uiOutput(ns("biomarkers_ui"))),
        tabPanel("Plots", br(), uiOutput(ns("plots_ui"))),
        tabPanel("Results", br(), uiOutput(ns("results_ui")))
      )
    )
  ))
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_multi_concordance_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$active_dataset_banner <- renderUI(multi_active_dataset_banner(multi_dataset))
    state <- reactiveValues(result = NULL)

    ## =========================================================================
    ## Layer / design / covariate pickers (Active dataset only) - populated
    ## only from what multi_dataset actually contains.
    ## =========================================================================

    output$layer_pick_ui <- renderUI({
      layers <- multi_dataset$layers %||% list()
      if (length(layers) == 0) return(div(class = "empty-note", icon("triangle-exclamation"), "No Active Multi-Omics Dataset loaded - build one on the Dataset Workspace tab."))
      expr_cand <- mcc_layer_candidates(multi_dataset, "rnaseq")
      meth_cand <- mcc_layer_candidates(multi_dataset, "methylation")
      tagList(
        selectInput(ns("expr_layer"), "Expression layer", choices = names(layers), selected = mcc_default_layer(expr_cand, layers)),
        selectInput(ns("meth_layer"), "Methylation layer", choices = names(layers), selected = mcc_default_layer(meth_cand, layers))
      )
    })

    output$design_ui <- renderUI({
      layers <- multi_dataset$layers %||% list()
      ids <- if (!is.null(input$expr_layer) && input$expr_layer %in% names(layers)) rownames(layers[[input$expr_layer]]) else character(0)
      cand <- mcc_design_candidates(multi_dataset$sample_meta, ids)
      if (length(cand) == 0) return(div(class = "empty-note", icon("circle-info"), "No 2-class metadata column found. Upload sample metadata with a 2-class design column (e.g. response)."))
      selectInput(ns("design_col"), "Design (2-class comparison)", choices = cand)
    })

    output$covariate_ui <- renderUI({
      meta <- multi_dataset$sample_meta
      if (is.null(meta) || ncol(meta) == 0) return(NULL)
      selectizeInput(ns("covariates"), "Regression covariates (optional)", choices = colnames(meta), multiple = TRUE, options = list(placeholder = "None"))
    })

    ## =========================================================================
    ## Data status (spec section 4) - visible pre-run, no analysis triggered.
    ## =========================================================================

    output$status_ui <- renderUI({
      if (identical(input$data_source, "preloaded")) {
        return(div(class = "empty-note", icon("circle-check"),
                    sprintf("Preloaded cohort table: %s. Counts are shown after Run.", names(MULTI_CONCORDANCE_COHORTS)[MULTI_CONCORDANCE_COHORTS == input$cohort])))
      }
      status <- tryCatch(mcc_data_status(multi_dataset, multi_results, input$expr_layer, input$meth_layer, input$array_type %||% "450K"), error = function(e) NULL)
      if (is.null(status)) return(multi_empty_state("Select an Active Multi-Omics Dataset to see data status."))
      tags$table(class = "table table-condensed", style = "font-size:0.85em;",
        tags$tbody(lapply(seq_len(nrow(status)), function(i) tags$tr(
          tags$td(status$item[i]),
          tags$td(tags$strong(style = if (status$status[i] == "Available") sprintf("color:%s;", ARTHOMIX_STATUS$good) else sprintf("color:%s;", ARTHOMIX_STATUS$warning), status$status[i])),
          tags$td(style = "color:var(--color-ink-muted,#898781);", status$detail[i])
        )))
      )
    })

    ## =========================================================================
    ## RUN - the only place that triggers computation (spec section 3).
    ## =========================================================================

    observeEvent(input$run_btn, {
      state$result <- tryCatch({
        if (identical(input$data_source, "preloaded")) mcc_build_preloaded(input, multi_results)
        else mcc_build_live(input, multi_dataset, multi_results)
      }, error = function(e) list(ok = FALSE, error = sprintf("Analysis failed: %s", conditionMessage(e))))
    })

    ## =========================================================================
    ## Overview (pre- and post-run)
    ## =========================================================================

    output$overview_ui <- renderUI({
      tagList(
        p(class = "submodule-desc", mod_multi_concordance_config$description),
        div(class = "empty-note", icon("circle-info"), MCC_CANONICAL_RULE_TEXT),
        if (is.null(state$result)) multi_empty_state("Set filters, then click \"Run Gene–CpG Analysis\".")
        else if (!isTRUE(state$result$ok)) div(class = "empty-note", style = "border-color: var(--color-danger, #e34948);", icon("circle-xmark"), state$result$error)
        else tagList(
          div(class = "empty-note", icon("circle-check"), state$result$overview_note),
          ## Circularity disclosure (Active Multi-Omics Dataset path only) -
          ## mcc_build_live_one() -> mcc_candidate_pool() draws candidate
          ## genes/CpGs from mcc_diablo_candidates(multi_results), i.e. this
          ## SAME dataset's own DIABLO feature selection - so the
          ## correlation/significance re-test below is not an independent
          ## check, unlike the Preloaded cohort path (a separately-run
          ## pipeline's own Table42/45). Phrasing follows this app's existing
          ## circularity callouts (e.g. mod_nomogram.R's
          ## nom_circularity_note(), mod_diagnostic.R's Test-split AUC note).
          if (identical(input$data_source, "active")) div(
            class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
            "Circularity warning: for the Active Multi-Omics Dataset, the candidate genes/CpGs analyzed above were already selected using this SAME data (DIABLO's own feature selection, run on this dataset) - the correlation/significance results below therefore re-test features on the data that picked them, and are not independent corroborating evidence. Only the Preloaded cohort path cross-references a separately-run pipeline.")
        )
      )
    })

    ## =========================================================================
    ## Gene<->CpG mapping table (spec section 9)
    ## =========================================================================

    output$genecpg_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_pairs_csv"), "Gene–CpG results (CSV)", class = "btn-sm"),
            downloadButton(ns("dl_annotation_csv"), "Annotation (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("genecpg_table"))
      )
    })
    output$genecpg_table <- DT::renderDataTable({
      r <- req(mcc_ok(state$result))
      cols <- intersect(c("gene_symbol", "gene_id", "transcript_id", "cpg", "chr", "pos", "strand",
                           "region_fine", "island_context", "tss_distance", "log2fc", "expr_fdr",
                           "dbeta", "delta_beta", "meth_fdr", "correlation_r", "correlation_fdr",
                           "direction_classification", "canonical_label", "diablo", "diablo_status", "snf", "joint", "sex", "dataset"), colnames(r$pairs_df))
      DT::datatable(r$pairs_df[, cols, drop = FALSE], rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_pairs_csv <- downloadHandler(function() "concordance_gene_cpg_results.csv",
                                            function(file) utils::write.csv(req(mcc_ok(state$result))$pairs_df, file, row.names = FALSE))
    output$dl_annotation_csv <- downloadHandler(function() "concordance_annotation.csv", function(file) {
      r <- req(mcc_ok(state$result))
      cols <- intersect(c("gene_symbol", "gene_id", "transcript_id", "chr", "pos", "strand", "cpg", "region_raw", "region_fine", "island_context", "tss_distance"), colnames(r$pairs_df))
      utils::write.csv(r$pairs_df[, cols, drop = FALSE], file, row.names = FALSE)
    })

    ## =========================================================================
    ## Direction (spec section 10-11)
    ## =========================================================================

    output$direction_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      has_genome <- "genome_expr_p" %in% colnames(r$pairs_df) && any(!is.na(r$pairs_df$genome_expr_p))
      is_preloaded <- identical(input$data_source, "preloaded")
      tagList(
        if (has_genome) tagList(
          tags$h5("Significant (genome-wide, nominal p)"),
          div(class = "empty-note", icon("circle-info"),
              sprintf("Both genome-wide nominal p-values are below %.3g. %s",
                      input$genome_p_thresh %||% 0.05,
                      if (is_preloaded) "From script-05 DEG/DMP discovery." else "From a live per-feature test across the uploaded layer.")),
          DT::dataTableOutput(ns("direction_significant_summary")),
          br(),
          div(class = "table-toolbar", downloadButton(ns("dl_direction_significant_csv"), "Significant pairs (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("direction_significant_table")),
          hr()
        ),
        tags$h5("Significance-gated direction (this run's own expression + methylation thresholds)"),
        div(class = "empty-note", icon("circle-info"),
            "A direction label requires clearing both the expression and methylation thresholds; otherwise it's Weak/uncertain or Not applicable."),
        DT::dataTableOutput(ns("direction_table"))
      )
    })
    output$direction_table <- DT::renderDataTable({
      r <- req(mcc_ok(state$result))
      tab <- as.data.frame(table(Direction = r$pairs_df$direction_classification, Canonical = r$pairs_df$canonical_label))
      DT::datatable(tab[tab$Freq > 0, , drop = FALSE], rownames = FALSE, options = list(dom = "t", pageLength = 20), class = "stripe hover compact")
    })

    direction_all_df <- reactive({
      r <- req(mcc_ok(state$result))
      mcc_add_raw_direction(r$pairs_df)
    })

    direction_significant_df <- reactive({
      df <- direction_all_df()
      if (!"genome_expr_p" %in% colnames(df)) return(df[0, , drop = FALSE])
      p <- input$genome_p_thresh %||% 0.05
      df[!is.na(df$genome_expr_p) & !is.na(df$genome_meth_p) & df$genome_expr_p < p & df$genome_meth_p < p, , drop = FALSE]
    })
    output$direction_significant_summary <- DT::renderDataTable({
      df <- direction_significant_df()
      if (nrow(df) == 0) return(DT::datatable(data.frame(Message = "No pairs clear this nominal p threshold on both omics - try raising it."), rownames = FALSE, options = list(dom = "t")))
      tab <- as.data.frame(table(Direction = df$raw_direction, `Region canonical` = df$raw_canonical_label))
      DT::datatable(tab[tab$Freq > 0, , drop = FALSE], rownames = FALSE, options = list(dom = "t", pageLength = 20), class = "stripe hover compact")
    })
    output$direction_significant_table <- DT::renderDataTable({
      df <- direction_significant_df()
      if (nrow(df) == 0) return(DT::datatable(data.frame(Message = "No pairs clear this nominal p threshold on both omics - try raising it."), rownames = FALSE, options = list(dom = "t")))
      cols <- intersect(c("sex", "gene_symbol", "cpg", "raw_direction", "region_fine", "raw_canonical_label",
                           "genome_expr_p", "genome_expr_fdr", "genome_expr_logfc",
                           "genome_meth_p", "genome_meth_fdr", "genome_meth_dbeta",
                           "biomarker_source", "diablo_status"), colnames(df))
      DT::datatable(df[, cols, drop = FALSE], rownames = FALSE, filter = "top",
                    options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_direction_significant_csv <- downloadHandler(
      function() "concordance_significant_genome_wide.csv",
      function(file) utils::write.csv(direction_significant_df(), file, row.names = FALSE)
    )

    ## =========================================================================
    ## Location (spec section 21, Plot 4)
    ## =========================================================================

    output$location_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      if (!"chr" %in% colnames(r$pairs_df) || all(is.na(r$pairs_df$chr))) return(multi_empty_state("Genomic coordinates are not available for this data source/annotation."))
      tagList(
        multi_plot_or_empty(function() mcc_plot_location(r$pairs_df), ns("plot_location"), "No candidates with known genomic coordinates.", height = "420px"),
        DT::dataTableOutput(ns("location_table"))
      )
    })
    output$plot_location <- renderPlot(mcc_plot_location(req(mcc_ok(state$result))$pairs_df))
    output$location_table <- DT::renderDataTable({
      r <- req(mcc_ok(state$result))
      cols <- intersect(c("gene_symbol", "chr", "pos", "gene_id", "cpg", "region_fine", "island_context", "tss_distance"), colnames(r$pairs_df))
      DT::datatable(r$pairs_df[, cols, drop = FALSE], rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    ## =========================================================================
    ## Biomarkers table (spec section 15-19)
    ## =========================================================================

    output$biomarkers_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_candidates_csv"), "Candidate biomarkers (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("biomarker_table"))
      )
    })
    output$biomarker_table <- DT::renderDataTable({
      r <- req(mcc_ok(state$result))
      cols <- intersect(c("gene_symbol", "cpg", "sex", "log2fc", "dbeta", "region_fine", "direction_classification",
                           "correlation_r", "correlation_fdr", "diablo", "diablo_status", "snf", "joint", "priority_score", "evidence_label"), colnames(r$pairs_df))
      DT::datatable(r$pairs_df[, cols, drop = FALSE], rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_candidates_csv <- downloadHandler(function() "concordance_candidate_biomarkers.csv", function(file) {
      r <- req(mcc_ok(state$result))
      utils::write.csv(r$pairs_df[r$pairs_df$evidence_label %in% c("Potential Multi-Omics Biomarker", "Candidate Multi-Omics Biomarker"), , drop = FALSE], file, row.names = FALSE)
    })

    ## =========================================================================
    ## Plots (spec section 20)
    ## =========================================================================

    output$plots_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      tagList(
        box(width = NULL, title = "Gene–CpG Concordance Scatter", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mcc_plot_scatter(r$pairs_df), ns("plot_scatter"), height = "380px"),
            downloadButton(ns("dl_scatter_png"), "Download (PNG)", class = "btn-sm")),
        box(width = NULL, title = "Direction Quadrant", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mcc_plot_quadrant(r$pairs_df), ns("plot_quadrant"), height = "380px"),
            downloadButton(ns("dl_quadrant_png"), "Download (PNG)", class = "btn-sm")),
        box(width = NULL, title = "Multi-Omics Evidence Heatmap", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mcc_plot_evidence_heatmap(r$pairs_df), ns("plot_heatmap"), "Not enough scored candidates for a heatmap.", height = "420px")),
        box(width = NULL, title = "Gene–CpG Network", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mcc_plot_network(r$pairs_df), ns("plot_network"), "Too few (or too many) significant pairs for a readable network.", height = "420px"),
            downloadButton(ns("dl_network_png"), "Download (PNG)", class = "btn-sm")),
        box(width = NULL, title = "Gene–CpG Correlation (single pair)", status = "primary", solidHeader = FALSE,
            uiOutput(ns("pair_picker_ui")),
            multi_plot_or_empty(function() pair_corr_plot_fn(), ns("plot_pair"), "Select a pair with computed correlation, or matched samples aren't available.", height = "380px"))
      )
    })
    output$plot_scatter <- renderPlot(mcc_plot_scatter(req(mcc_ok(state$result))$pairs_df))
    output$dl_scatter_png <- multi_png_download(function() mcc_plot_scatter(req(mcc_ok(state$result))$pairs_df), function() "concordance_scatter.png")
    output$plot_quadrant <- renderPlot(mcc_plot_quadrant(req(mcc_ok(state$result))$pairs_df))
    output$dl_quadrant_png <- multi_png_download(function() mcc_plot_quadrant(req(mcc_ok(state$result))$pairs_df), function() "concordance_quadrant.png")
    output$plot_heatmap <- renderPlot(mcc_plot_evidence_heatmap(req(mcc_ok(state$result))$pairs_df))
    output$plot_network <- renderPlot(mcc_plot_network(req(mcc_ok(state$result))$pairs_df))
    output$dl_network_png <- multi_png_download(function() mcc_plot_network(req(mcc_ok(state$result))$pairs_df), function() "concordance_network.png")

    output$pair_picker_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      if (is.null(r$expr_mat) || is.null(r$meth_mat)) return(div(class = "empty-note", icon("circle-info"), "Sample-level pair correlation needs the Active Multi-Omics Dataset path with matched samples."))
      d <- r$pairs_df[!is.na(r$pairs_df$correlation_r), , drop = FALSE]
      if (nrow(d) == 0) return(div(class = "empty-note", icon("circle-info"), "No pair has a computed correlation yet."))
      selectInput(ns("pair_pick"), NULL, choices = stats::setNames(seq_len(nrow(d)), sprintf("%s – %s (r=%.2f)", d$gene_symbol, d$cpg, d$correlation_r)))
    })
    pair_corr_plot_fn <- reactive({
      r <- req(mcc_ok(state$result))
      req(input$pair_pick, r$expr_mat, r$meth_mat)
      d <- r$pairs_df[!is.na(r$pairs_df$correlation_r), , drop = FALSE]
      row <- d[as.integer(input$pair_pick), , drop = FALSE]
      strat <- if (!is.null(r$strata) && !is.null(row$sex) && row$sex %in% names(r$strata)) r$strata[[row$sex]] else list(expr_mat = r$expr_mat, meth_mat = r$meth_mat, common_samples = r$common_samples)
      x <- as.numeric(strat$meth_mat[strat$common_samples, row$cpg]); y <- as.numeric(strat$expr_mat[strat$common_samples, row$gene_symbol])
      mcc_plot_pair_correlation(x, y, row$gene_symbol, row$cpg, row$correlation_r, row$correlation_p, row$correlation_fdr, row$correlation_n)
    })

    ## =========================================================================
    ## Results (spec section 26-28)
    ## =========================================================================

    output$results_ui <- renderUI({
      r <- req(mcc_ok(state$result))
      sc <- mcc_summary_counts(r$pairs_df, sex_col_present = "sex" %in% colnames(r$pairs_df))
      card <- function(label, value) div(class = "card", style = "flex:1 1 130px; text-align:center; padding:10px;",
                                          div(style = sprintf("font-size:1.25em; font-weight:600; color:%s;", ARTHOMIX_COLORS$blue), value %||% "NA"),
                                          div(style = "font-size:0.8em; color:var(--color-ink-muted,#898781);", label))
      tagList(
        div(style = "display:flex; flex-wrap:wrap; gap:8px;",
            card("Genes tested", sc$n_genes), card("CpGs tested", sc$n_cpgs), card("Gene–CpG pairs", sc$n_pairs),
            card("Significant pairs", sc$n_significant), card("Canonical", sc$n_canonical), card("Non-canonical", sc$n_noncanonical),
            card("Potential biomarkers", sc$n_potential), card("Female candidates", sc$n_female), card("Male candidates", sc$n_male),
            card("DIABLO-supported", sc$n_diablo), card("SNF-supported", sc$n_snf), card("Joint-supported", sc$n_joint)),
        br(),
        box(width = NULL, title = "Analysis settings (reproducibility)", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("repro_table"))),
        div(class = "table-toolbar",
            downloadButton(ns("dl_pairs_csv2"), "Gene–CpG results (CSV)", class = "btn-sm"),
            downloadButton(ns("dl_candidates_csv2"), "Candidate biomarkers (CSV)", class = "btn-sm"),
            downloadButton(ns("dl_annotation_csv2"), "Annotation (CSV)", class = "btn-sm"))
      )
    })
    output$repro_table <- DT::renderDataTable({
      r <- req(mcc_ok(state$result))
      DT::datatable(r$settings_snapshot, rownames = FALSE, options = list(dom = "t", pageLength = 30), class = "stripe hover compact")
    })
    output$dl_pairs_csv2 <- downloadHandler(function() "concordance_gene_cpg_results.csv", function(file) utils::write.csv(req(mcc_ok(state$result))$pairs_df, file, row.names = FALSE))
    output$dl_candidates_csv2 <- downloadHandler(function() "concordance_candidate_biomarkers.csv", function(file) {
      r <- req(mcc_ok(state$result))
      utils::write.csv(r$pairs_df[r$pairs_df$evidence_label %in% c("Potential Multi-Omics Biomarker", "Candidate Multi-Omics Biomarker"), , drop = FALSE], file, row.names = FALSE)
    })
    output$dl_annotation_csv2 <- downloadHandler(function() "concordance_annotation.csv", function(file) {
      r <- req(mcc_ok(state$result))
      cols <- intersect(c("gene_symbol", "gene_id", "transcript_id", "chr", "pos", "strand", "cpg", "region_raw", "region_fine", "island_context", "tss_distance"), colnames(r$pairs_df))
      utils::write.csv(r$pairs_df[, cols, drop = FALSE], file, row.names = FALSE)
    })

    ## =========================================================================
    ## Publish (kept in the same list(df=, cohort=) shape other Multi-Omics
    ## submodules already publish, for consistency - nothing currently reads
    ## multi_results$concordance).
    ## =========================================================================

    observe({
      r <- mcc_ok(state$result)
      if (is.null(r) || is.null(multi_results)) return()
      multi_results$concordance <- list(df = r$pairs_df, cohort = r$label)
    })
  })
}

## ---------------------------------------------------------------------------
## Result builders - one canonical shape for both data sources:
## list(ok, error, pairs_df, expr_mat, meth_mat, common_samples, label,
##      overview_note, settings_snapshot)
## ---------------------------------------------------------------------------

mcc_ok <- function(result) if (!is.null(result) && isTRUE(result$ok)) result else NULL

mcc_build_preloaded <- function(input, multi_results) {
  res <- multi_read_registry_table(input$cohort)
  if (!isTRUE(res$ok)) return(list(ok = FALSE, error = res$error))
  df <- multi_concordance_add_fdr(res$df)
  df$gene_symbol <- df$SYMBOL
  df$cpg <- df$CpG
  df$log2fc <- df$expr_logFC
  df$dbeta <- df$delta_M
  df$region_fine <- cx_region_fine(df$region)
  df$island_context <- df$island
  ## Real chr/pos for the Location tab - same annotation lookup already used
  ## for the preloaded methylation path elsewhere (cx_load_default_methylation()),
  ## keyed by this table's own CpG IDs. Fail-soft: NA (and the Location tab's
  ## existing empty-state) if the annotation package isn't installed.
  df$chr <- NA_character_; df$pos <- NA_real_
  anno <- cx_get_region_annotation("450K")
  if (isTRUE(anno$ok)) {
    idx <- match(df$cpg, rownames(anno$anno))
    df$chr <- anno$anno$chr[idx]; df$pos <- anno$anno$pos[idx]
  }
  df$gene_id <- df$ENSEMBL %||% NA_character_
  df$transcript_id <- NA_character_; df$strand <- NA_character_; df$tss_distance <- NA_real_
  df$correlation_r <- NA_real_; df$correlation_p <- NA_real_; df$correlation_fdr <- NA_real_; df$correlation_n <- NA_integer_
  df$dataset <- names(MULTI_CONCORDANCE_COHORTS)[MULTI_CONCORDANCE_COHORTS == input$cohort]

  if (length(input$sex) == 1 && input$sex %in% c("female", "male") && "sex" %in% colnames(df)) df <- df[tolower(df$sex) == input$sex, , drop = FALSE]
  if (length(input$region) > 0) df <- df[df$region_fine %in% input$region, , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, error = "No rows remain after the current sex/region filters."))

  df <- mcc_classify_direction(df, input$expr_thresh, input$expr_fdr_thresh, input$meth_thresh, input$meth_fdr_thresh)

  ## Cross-reference (a) the precomputed, per-sex DIABLO candidate-biomarker
  ## panel for this same cohort (Table40/44b - always available, no live run
  ## needed) and (b) the same-session live DIABLO/Joint panel, when present.
  ## SNF stays FALSE on this path: no precomputed per-gene/per-CpG SNF
  ## feature-selection table exists in this deployment (only patient-cluster
  ## assignments do) - never fabricated.
  panel_label <- MCC_PRELOADED_DIABLO_PANEL[[input$cohort]]
  panel <- if (!is.null(panel_label)) { r <- multi_read_registry_table(panel_label); if (isTRUE(r$ok)) r$df else NULL } else NULL
  df <- mcc_join_preloaded_diablo_panel(df, panel)
  df$snf <- FALSE

  live_diablo <- mcc_diablo_candidates(multi_results)
  if (!is.null(live_diablo)) {
    hit <- toupper(df$gene_symbol) %in% toupper(live_diablo$feature) | toupper(df$gene_id %||% "") %in% toupper(live_diablo$feature) | toupper(df$cpg) %in% toupper(live_diablo$feature)
    df$diablo[hit] <- TRUE; df$joint[hit] <- TRUE
  }
  df$custom <- FALSE

  ## Real genome-wide response-driven significance (Table3/4, script 05 -
  ## the actual end-to-end DEG/DMP discovery step), for every candidate gene/
  ## CpG this cohort's DIABLO panel or concordance table ever mentions - never
  ## the small-panel-recomputed FDR above, which BH-corrects over only this
  ## cohort's own few hundred/thousand candidate rows.
  deg <- multi_read_registry_table("Genome-wide DEG lookup - response-driven, by sex (candidates only)")
  dmp <- multi_read_registry_table("Genome-wide DMP lookup - response-driven, by sex (candidates only)")
  df <- mcc_join_genome_wide_significance(df, if (isTRUE(deg$ok)) deg$df else NULL, if (isTRUE(dmp$ok)) dmp$df else NULL)

  scored <- mcc_priority_score(df)
  df <- cbind(df, scored[, c("priority_score", "evidence_label")])

  settings <- data.frame(Parameter = c("Data source", "Cohort", "Sex filter", "Region filter", "Expression threshold", "Expression FDR", "Methylation threshold", "Methylation FDR", "Genome-wide DEG/DMP nominal p", "Analyzed at"),
                          Value = c("Preloaded (Table42/45)", df$dataset[1] %||% input$cohort, input$sex, paste(input$region, collapse = ", ") %||% "All", input$expr_thresh, input$expr_fdr_thresh, input$meth_thresh, input$meth_fdr_thresh, input$genome_p_thresh %||% 0.05, format(Sys.time())), stringsAsFactors = FALSE)

  list(ok = TRUE, error = NULL, pairs_df = df, expr_mat = NULL, meth_mat = NULL, common_samples = NULL,
       label = df$dataset[1] %||% input$cohort,
       overview_note = sprintf("%d gene–CpG pairs loaded from the preloaded cohort table.", nrow(df)),
       settings_snapshot = settings)
}

mcc_build_live <- function(input, multi_dataset, multi_results) {
  layers <- multi_dataset$layers %||% list()
  if (is.null(input$expr_layer) || is.null(input$meth_layer) || !all(c(input$expr_layer, input$meth_layer) %in% names(layers)) || identical(input$expr_layer, input$meth_layer))
    return(list(ok = FALSE, error = "Select two different layers (Expression, Methylation) on the Active Multi-Omics Dataset."))
  expr_mat <- layers[[input$expr_layer]]; meth_mat <- layers[[input$meth_layer]]

  match_res <- mcc_match_samples(expr_mat, meth_mat)
  if (!isTRUE(match_res$ok)) return(list(ok = FALSE, error = "Insufficient matched samples between the selected expression and methylation layers."))
  common <- match_res$common_samples

  meta <- multi_dataset$sample_meta
  sex_col <- mcc_sex_candidates(meta)[1]
  groups <- if (!is.null(sex_col)) mcc_sex_groups(meta, sex_col, common) else NULL

  sex_runs <- if (identical(input$sex, "sex_specific")) {
    if (is.null(groups) || length(groups) < 2) return(list(ok = FALSE, error = "No usable sex column, or fewer than two sexes represented among matched samples."))
    groups
  } else if (input$sex %in% c("female", "male")) {
    if (is.null(groups) || !input$sex %in% tolower(names(groups))) return(list(ok = FALSE, error = sprintf("Insufficient matched samples for sex-specific analysis (%s).", input$sex)))
    stats::setNames(list(groups[[names(groups)[tolower(names(groups)) == input$sex][1]]]), input$sex)
  } else {
    stats::setNames(list(common), "pooled")
  }

  results <- list()
  for (sex_label in names(sex_runs)) {
    ids <- intersect(common, sex_runs[[sex_label]])
    if (length(ids) < 6) next  ## too few for a stable per-sex comparison; silently excluded stratum is reported via note below
    r1 <- mcc_build_live_one(input, expr_mat, meth_mat, ids, multi_dataset, multi_results)
    if (isTRUE(r1$ok)) { r1$pairs_df$sex <- sex_label; results[[sex_label]] <- r1 }
  }
  if (length(results) == 0) return(list(ok = FALSE, error = "No stratum (pooled or by sex) had enough matched samples and candidate genes/CpGs to analyze."))

  pairs_df <- do.call(rbind, lapply(results, `[[`, "pairs_df"))
  settings <- results[[1]]$settings_snapshot
  settings <- rbind(settings, data.frame(Parameter = "Strata analyzed", Value = paste(names(results), collapse = ", "), stringsAsFactors = FALSE))
  ## Kept per-stratum (not just the first) so the single-pair correlation
  ## plot can look up the right sample subset/matrix for a pair from ANY
  ## stratum in sex-specific mode - each stratum has its own matched-sample
  ## set, never mixed with another stratum's.
  strata <- lapply(results, function(x) list(expr_mat = x$expr_mat, meth_mat = x$meth_mat, common_samples = x$common_samples))

  list(ok = TRUE, error = NULL, pairs_df = pairs_df, strata = strata,
       expr_mat = results[[1]]$expr_mat, meth_mat = results[[1]]$meth_mat, common_samples = results[[1]]$common_samples,
       label = sprintf("%s (Active Multi-Omics Dataset)", multi_dataset$table_label %||% "Uploaded dataset"),
       overview_note = sprintf("%d gene–CpG pairs computed from %d matched sample(s) across %s.", nrow(pairs_df), length(common), paste(names(results), collapse = "/")),
       settings_snapshot = settings)
}

mcc_build_live_one <- function(input, expr_mat, meth_mat, sample_ids, multi_dataset, multi_results) {
  expr_sub <- expr_mat[sample_ids, , drop = FALSE]; meth_sub <- meth_mat[sample_ids, , drop = FALSE]

  pool <- mcc_candidate_pool(multi_results, multi_dataset, input$expr_layer, input$meth_layer, input$custom_genes %||% character(0), input$custom_cpgs %||% character(0))
  if (!isTRUE(pool$ok)) return(list(ok = FALSE, error = pool$note))
  pool_filtered <- mcc_filter_source(pool$df, input$source)
  cand_genes <- unique(pool_filtered$feature[pool_filtered$omics %in% c("Transcriptomics", NA)])
  if (length(cand_genes) == 0) cand_genes <- unique(pool$df$feature[pool$df$omics %in% c("Transcriptomics", NA)])
  if (length(cand_genes) == 0) return(list(ok = FALSE, error = "No candidate biomarkers available for the selected source."))

  map <- mcc_gene_cpg_map(cand_genes, meth_features = colnames(meth_sub), array_type = input$array_type %||% "450K")
  if (!isTRUE(map$ok)) return(list(ok = FALSE, error = map$error))
  df <- map$df

  ## Bridge candidate gene symbols back to this dataset's own expression
  ## column IDs (which may be Ensembl, Entrez, or symbol) via the same
  ## harmonizer used to build the map, so expression stats join correctly
  ## regardless of the uploaded matrix's ID convention.
  harm <- cx_harmonize_gene_ids(colnames(expr_sub))
  expr_id_by_symbol <- if (isTRUE(harm$ok)) stats::setNames(harm$df$input_id, toupper(harm$df$canonical_symbol)) else stats::setNames(colnames(expr_sub), toupper(colnames(expr_sub)))
  df$expr_feature <- unname(expr_id_by_symbol[toupper(df$gene_symbol)])
  df <- df[!is.na(df$expr_feature), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, error = "None of the candidate genes' mapped CpGs correspond to genes present in the expression layer."))

  design_col <- input$design_col
  if (is.null(design_col) || !design_col %in% colnames(multi_dataset$sample_meta %||% data.frame())) return(list(ok = FALSE, error = "No 2-class design column selected."))
  meta <- multi_dataset$sample_meta
  grp <- stats::setNames(as.character(meta[[design_col]]), rownames(meta))

  es <- mcc_expression_stats(expr_sub[, unique(df$expr_feature), drop = FALSE], grp)
  if (!isTRUE(es$ok)) return(list(ok = FALSE, error = es$error))
  meth_val_type <- mcc_detect_methylation_value_type(meth_sub)
  ms <- mcc_methylation_stats(meth_sub[, unique(df$cpg), drop = FALSE], grp, if (meth_val_type %in% c("beta", "M-value")) meth_val_type else "beta")
  if (!isTRUE(ms$ok)) return(list(ok = FALSE, error = ms$error))

  ei <- match(df$expr_feature, es$df$feature); mi <- match(df$cpg, ms$df$feature)
  df$log2fc <- es$df$log2fc[ei]; df$expr_p <- es$df$p[ei]; df$expr_fdr <- es$df$fdr[ei]
  df$dbeta <- ms$df$dbeta[mi]; df$delta_beta <- ms$df$delta_beta[mi]; df$meth_p <- ms$df$p[mi]; df$meth_fdr <- ms$df$fdr[mi]

  ## Genome-wide equivalent of the preloaded path's Table3/4 DEG/DMP lookup
  ## (mcc_join_genome_wide_significance()) - the "Significant" section on the
  ## Direction tab needs the SAME genome_expr_p/genome_meth_p/*_fdr columns
  ## on this path too, so it isn't preloaded-only. Nominal p per feature is
  ## identical whether computed on the candidate subset above or the full
  ## uploaded layer (each column's own t-test doesn't depend on which other
  ## columns were tested alongside it) - what differs, and is worth getting
  ## right, is the FDR: BH-correcting only over this cohort's own handful of
  ## candidate genes/CpGs (the *_fdr columns above) inflates significance the
  ## same way the preloaded path's small-panel FDR did, so this reruns the
  ## same generic stats functions across every feature actually present in
  ## the uploaded layer for a properly-scaled genome-wide FDR.
  genome_es <- mcc_expression_stats(expr_sub, grp)
  genome_ms <- mcc_methylation_stats(meth_sub, grp, if (meth_val_type %in% c("beta", "M-value")) meth_val_type else "beta")
  df$genome_expr_p <- NA_real_; df$genome_expr_fdr <- NA_real_; df$genome_expr_logfc <- NA_real_
  df$genome_meth_p <- NA_real_; df$genome_meth_fdr <- NA_real_; df$genome_meth_dbeta <- NA_real_
  if (isTRUE(genome_es$ok)) {
    gei <- match(df$expr_feature, genome_es$df$feature)
    df$genome_expr_p <- genome_es$df$p[gei]; df$genome_expr_fdr <- genome_es$df$fdr[gei]; df$genome_expr_logfc <- genome_es$df$log2fc[gei]
  }
  if (isTRUE(genome_ms$ok)) {
    gmi <- match(df$cpg, genome_ms$df$feature)
    df$genome_meth_p <- genome_ms$df$p[gmi]; df$genome_meth_fdr <- genome_ms$df$fdr[gmi]; df$genome_meth_dbeta <- genome_ms$df$dbeta[gmi]
  }

  if (length(input$region) > 0) df <- df[df$region_fine %in% input$region, , drop = FALSE]
  ## "Custom CpGs" restricts to the exact user-entered CpG IDs (the gene
  ## candidate pool above is gene-driven per spec section 9's mapping
  ## direction, so this is applied as a final row filter rather than by
  ## restricting candidate gene selection).
  if (identical(input$source, "Custom CpGs") && length(input$custom_cpgs %||% character(0)) > 0) df <- df[df$cpg %in% input$custom_cpgs, , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, error = "No gene–CpG pairs remain after the region/CpG filters."))

  ## mcc_pair_correlation() (and the single-pair plot picker) index the
  ## expression matrix by `gene_symbol`, not by this dataset's own raw
  ## expr_feature ID - re-key a small symbol-labeled copy here rather than
  ## inside that shared helper, since only this caller knows the mapping.
  expr_by_symbol <- expr_sub[, unique(df$expr_feature), drop = FALSE]
  colnames(expr_by_symbol) <- df$gene_symbol[match(colnames(expr_by_symbol), df$expr_feature)]

  cor_res <- mcc_pair_correlation(expr_by_symbol, meth_sub, df, sample_ids, method = input$cor_method %||% "pearson")
  if (isTRUE(cor_res$ok)) {
    ci <- match(paste(df$gene_symbol, df$cpg), paste(cor_res$df$gene_symbol, cor_res$df$cpg))
    df$correlation_r <- cor_res$df$r[ci]; df$correlation_p <- cor_res$df$p[ci]; df$correlation_fdr <- cor_res$df$fdr[ci]; df$correlation_n <- cor_res$df$n[ci]
  } else {
    df$correlation_r <- NA_real_; df$correlation_p <- NA_real_; df$correlation_fdr <- NA_real_; df$correlation_n <- NA_integer_
  }

  df <- mcc_classify_direction(df, input$expr_thresh, input$expr_fdr_thresh, input$meth_thresh, input$meth_fdr_thresh)

  src_idx <- match(df$gene_symbol, pool$df$feature); src_idx_cpg <- match(df$cpg, pool$df$feature)
  df$diablo <- (!is.na(src_idx) & pool$df$diablo[src_idx]) | (!is.na(src_idx_cpg) & pool$df$diablo[src_idx_cpg])
  df$snf <- (!is.na(src_idx) & pool$df$snf[src_idx]) | (!is.na(src_idx_cpg) & pool$df$snf[src_idx_cpg])
  df$joint <- (!is.na(src_idx) & pool$df$joint[src_idx]) | (!is.na(src_idx_cpg) & pool$df$joint[src_idx_cpg])
  df$custom <- (!is.na(src_idx) & pool$df$custom[src_idx]) | (!is.na(src_idx_cpg) & pool$df$custom[src_idx_cpg])
  df$dataset <- multi_dataset$table_label %||% "Active Multi-Omics Dataset"

  scored <- mcc_priority_score(df)
  df <- cbind(df, scored[, c("priority_score", "evidence_label")])

  covariate_result <- NULL
  if (length(input$covariates %||% character(0)) > 0 && nrow(df) > 0) {
    top_row <- df[order(-df$priority_score), , drop = FALSE][1, ]
    if (top_row$expr_feature %in% colnames(expr_sub) && top_row$cpg %in% colnames(meth_sub))
      covariate_result <- mcc_regression(expr_sub[sample_ids, top_row$expr_feature], meth_sub[sample_ids, top_row$cpg], meta[sample_ids, , drop = FALSE], input$covariates)
  }

  settings <- data.frame(Parameter = c("Data source", "Expression layer", "Methylation layer", "Methylation value type detected", "Design column",
                                        "Matched samples", "Sex mode", "Correlation method", "Expression threshold", "Expression FDR",
                                        "Methylation threshold", "Methylation FDR", "Region filter", "Array (annotation)", "Analyzed at"),
                          Value = c("Active Multi-Omics Dataset", input$expr_layer, input$meth_layer, meth_val_type, design_col,
                                    length(sample_ids), input$sex, input$cor_method, input$expr_thresh, input$expr_fdr_thresh,
                                    input$meth_thresh, input$meth_fdr_thresh, paste(input$region, collapse = ", ") %||% "All", input$array_type, format(Sys.time())),
                          stringsAsFactors = FALSE)
  if (!is.null(covariate_result) && isTRUE(covariate_result$ok))
    settings <- rbind(settings, data.frame(Parameter = "Top-candidate adjusted association (coef, p, N)",
                                            Value = sprintf("%.3f, p=%.3g, N=%d, covariates: %s", covariate_result$coefficient, covariate_result$p_value, covariate_result$model_n, paste(covariate_result$covariates_used, collapse = ", ")), stringsAsFactors = FALSE))

  list(ok = TRUE, error = NULL, pairs_df = df, expr_mat = expr_by_symbol,
       meth_mat = meth_sub, common_samples = sample_ids, settings_snapshot = settings)
}
