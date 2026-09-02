## R/multiomics/02_Cohort_Harmonization/mod_multi_overview.R
## Submodule: Cohort Harmonization - data-adaptive report on the Active
## Multi-Omics Dataset (multi_dataset, built on the Dataset Workspace tab):
## which modalities are present, which samples are actually shared, and
## whether integration is feasible.
## Never assumes RNA-seq + methylation specifically - every check inspects
## whatever modalities the active dataset actually contains, and reports
## "Not detected"/"Insufficient information" rather than guessing when it
## can't tell (cohort_harmonization_helpers.R). Nothing renders until the
## relevant blue button is clicked (mirrors the Dataset Workspace's own
## button-gated pattern).

mod_multi_overview_config <- list(
  id = "overview", title = "Cohort Harmonization", icon = "users", group = "Data",
  description = "Modality availability, sample matching, and integration readiness for the active dataset."
)

mod_multi_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("active_dataset_banner")),
    uiOutput(ns("body_ui"))
  )
}

mod_multi_overview_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$active_dataset_banner <- renderUI(multi_active_dataset_banner(multi_dataset))

    descriptors <- reactive(ch_modality_descriptors(multi_dataset))
    id_sets_all <- reactive(stats::setNames(lapply(descriptors(), function(x) x$sample_ids), names(descriptors())))
    pheno_candidates <- reactive(ch_detect_candidate_columns(multi_dataset$sample_meta, "phenotype"))
    batch_candidates <- reactive(ch_detect_candidate_columns(multi_dataset$sample_meta, "batch"))

    ## Raw per-sample matrix source for PCA/correlation only -
    ## descriptors()/harmonization_result() above are untouched and keep
    ## reading multi_dataset directly (real, patient-ID-only for the
    ## preloaded cohort, per ch_modality_descriptors_preloaded()). The
    ## preloaded cohort's modality descriptors never carry a raw matrix
    ## (has_raw_matrix=FALSE for all ~80 patients - no bundled full-cohort
    ## matrix), so these two panels fall back to one analysis cell's own
    ## saved-fit matrices (mi_preloaded_cell_dataset(), the same live-recompute
    ## mechanism mod_multi_integration.R already uses) rather than staying
    ## permanently inert for the preloaded path.
    ov_raw_dataset <- reactive({
      if (identical(multi_dataset$source, "preloaded")) {
        req(input$preloaded_cell)
        mi_preloaded_cell_dataset(input$preloaded_cell)
      } else {
        list(ok = TRUE, layers = multi_dataset$layers %||% list(), sample_meta = multi_dataset$sample_meta, provenance = NULL)
      }
    })

    ## ---- Filters + buttons + tabset - only built once an active dataset
    ## with at least one detectable modality exists (spec section 15: don't
    ## show filters that don't apply to the actual data). ----
    output$body_ui <- renderUI({
      if (is.null(multi_dataset) || !isTRUE(multi_dataset$active)) {
        return(div(class = "empty-note", icon("circle-info"),
                   "No Active Multi-Omics Dataset yet - build one on the Dataset Workspace tab, then return here."))
      }
      d <- descriptors()
      if (length(d) == 0) {
        return(div(class = "empty-note", icon("triangle-exclamation"), "Could not determine any modalities for the active dataset."))
      }
      bc <- batch_candidates()
      fluidRow(
        column(
          4,
          box(width = NULL, title = "Filters", status = "primary", solidHeader = FALSE,
              checkboxGroupInput(ns("sel_modalities"), "Select modalities", choices = names(d), selected = names(d)),
              numericInput(ns("min_overlap"), "Minimum sample overlap", value = 3, min = 1),
              if (length(bc) > 0) selectInput(ns("batch_col"), "Batch/cohort column", choices = c("(none)" = "", bc)),
              if (identical(multi_dataset$source, "preloaded")) tagList(
                selectInput(ns("preloaded_cell"), "Analysis cell (used for PCA and correlation)", choices = MULTI_CELL_CHOICES),
                div(class = "empty-note", icon("circle-info"), "Preloaded cohort has no bundled raw matrix. PCA/correlation below use one analysis cell's matched-sample subset, recomputed from its saved DIABLO fit.")
              )
          ),
          actionButton(ns("analyze_btn"), "Analyze Cohort", icon = icon("magnifying-glass-chart"), class = "btn-primary btn-sm", width = "100%")
        ),
        column(
          8,
          tabsetPanel(
            id = ns("ch_tabs"), type = "tabs",
            tabPanel("Overview", br(), conditionalPanel(condition = sprintf("input['%s'] > 0", ns("analyze_btn")), uiOutput(ns("overview_ui")))),
            tabPanel("Sample Match", br(), conditionalPanel(condition = sprintf("input['%s'] > 0", ns("analyze_btn")), uiOutput(ns("sample_match_ui")))),
            tabPanel("Sample Explorer", br(), conditionalPanel(condition = sprintf("input['%s'] > 0", ns("analyze_btn")), uiOutput(ns("sample_explorer_ui"))))
          )
        )
      )
    })

    ## =========================================================================
    ## Analyze Cohort - Overview + Sample Match, computed once per click
    ## (spec section 14: never reactively on every input change).
    ## =========================================================================
    harmonization_result <- eventReactive(input$analyze_btn, {
      d <- descriptors()
      sel <- intersect(input$sel_modalities %||% names(d), names(d))
      validate(need(length(sel) >= 1, "Select at least one modality."))
      d_sel <- d[sel]
      ids <- stats::setNames(lapply(d_sel, function(x) x$sample_ids), names(d_sel))
      min_overlap <- input$min_overlap %||% 3

      cells_res <- ch_analysis_cells(ids, pheno_available = length(pheno_candidates()) > 0, min_integration = min_overlap, min_prediction = max(6, min_overlap))
      readiness <- lapply(cells_res$cells, ch_integration_readiness, min_limited = min_overlap)
      id_table <- ch_id_harmonization_table(ids)

      full_overlap <- if (length(ids) > 0) Reduce(intersect, ids) else character(0)
      union_all <- if (length(ids) > 0) Reduce(union, ids) else character(0)

      list(
        ok = TRUE, descriptors = d_sel, ids = ids,
        overlap_matrix = ch_pairwise_overlap_matrix(ids),
        cells = cells_res$cells, cells_omitted_note = cells_res$omitted_note, readiness = readiness,
        id_table = id_table, n_total = length(union_all), n_matched = length(full_overlap),
        matched_summary = ch_matched_sample_summary(ids)
      )
    }, ignoreInit = TRUE)

    output$overview_ui <- renderUI({
      h <- tryCatch(harmonization_result(), error = function(e) NULL)
      if (is.null(h) || !isTRUE(h$ok)) return(multi_empty_state("Click \"Analyze Cohort\" to see results here."))
      ms <- h$matched_summary
      badge_color <- switch(ms$status, Matched = ARTHOMIX_COLORS$aqua, `Partially matched` = ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)
      tagList(
        div(class = "card", style = sprintf("padding:12px 14px; border-left:4px solid %s; margin-bottom:14px;", badge_color),
            span(style = sprintf("color:%s; font-weight:700;", badge_color), ms$status), " - ", ms$sentence),
        div(style = "display:flex; gap:14px; flex-wrap:wrap;",
            lapply(names(h$descriptors), function(nm) {
              desc <- h$descriptors[[nm]]
              st <- desc$status %||% list(level = "ready", label = "Ready", reasons = character(0))
              mo_dataset_block_card(nm, desc$n_samples, desc$n_features %||% NA, st)
            })),
        box(width = NULL, title = "Integration Readiness", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("readiness_table"))),
        if (!is.null(h$cells_omitted_note)) div(class = "empty-note", icon("circle-info"), h$cells_omitted_note),
        box(width = NULL, title = "Analysis Cells", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("cells_table"))),
        box(width = NULL, title = "Batch and Cohort Summary", status = "primary", solidHeader = FALSE,
            uiOutput(ns("batch_summary_ui"))),
        box(width = NULL, title = "Data Completeness", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(completeness_plot_fn, ns("completeness_plot"), "Not enough data to draw a completeness plot.", height = "260px"),
            div(class = "table-toolbar", downloadButton(ns("dl_completeness_png"), "Download plot (PNG)", class = "btn-sm")))
      )
    })
    output$readiness_table <- DT::renderDataTable({
      h <- req(harmonization_result())
      req(length(h$cells) > 0)
      df <- do.call(rbind, Map(function(cl, rd) data.frame(Cell = cl$label, Modalities = length(cl$modalities), `Matched samples` = cl$n_matched, Status = rd$label, Reason = rd$reason, check.names = FALSE), h$cells, h$readiness))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 20), class = "stripe hover compact")
    })
    output$cells_table <- DT::renderDataTable({
      h <- req(harmonization_result())
      req(length(h$cells) > 0)
      df <- do.call(rbind, lapply(h$cells, function(cl) data.frame(Cell = cl$label, `Matched samples` = cl$n_matched, `Feasible methods` = paste(cl$methods, collapse = "; "), check.names = FALSE)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 20), class = "stripe hover compact")
    })
    output$batch_summary_ui <- renderUI({
      bc <- input$batch_col
      if (is.null(bc) || !nzchar(bc)) return(div(class = "empty-note", icon("circle-info"), "No batch/cohort column selected."))
      meta_df <- multi_dataset$sample_meta
      if (is.null(meta_df) || !bc %in% colnames(meta_df)) return(multi_empty_state("Not detected."))
      tab <- table(as.character(meta_df[[bc]]))
      tagList(
        tags$table(class = "table table-condensed", style = "font-size:0.88em;",
                    tags$tbody(lapply(names(tab), function(l) tags$tr(tags$td(l), tags$td(tab[[l]]))))),
        multi_plot_or_empty(batch_plot_fn, ns("batch_plot"), height = "220px"),
        div(class = "table-toolbar", downloadButton(ns("dl_batch_bar_png"), "Download plot (PNG)", class = "btn-sm"))
      )
    })
    batch_plot_fn <- reactive({
      bc <- req(input$batch_col); req(nzchar(bc))
      ch_category_bar_plot(req(multi_dataset$sample_meta), bc, "Batch distribution")
    })
    output$batch_plot <- renderPlot(batch_plot_fn())
    output$dl_batch_bar_png <- multi_png_download(batch_plot_fn, function() "cohort_harmonization_batch_distribution.png")

    completeness_plot_fn <- reactive({ h <- req(harmonization_result()); ch_completeness_heatmap_plot(h$ids) })
    output$completeness_plot <- renderPlot(completeness_plot_fn(), alt = "Modality by sample data-completeness heatmap")
    output$dl_completeness_png <- multi_png_download(completeness_plot_fn, function() "cohort_harmonization_data_completeness.png")

    output$sample_match_ui <- renderUI({
      h <- tryCatch(harmonization_result(), error = function(e) NULL)
      if (is.null(h) || !isTRUE(h$ok)) return(multi_empty_state("Click \"Analyze Cohort\" to see results here."))
      tagList(
        box(width = NULL, title = "Pairwise Sample Overlap", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("overlap_table")),
            multi_plot_or_empty(overlap_heatmap_fn, ns("overlap_heatmap"), height = "300px"),
            div(class = "table-toolbar", downloadButton(ns("dl_overlap_heatmap_png"), "Download plot (PNG)", class = "btn-sm"))),
        box(width = NULL, title = "Sample-ID Harmonization", status = "primary", solidHeader = FALSE, collapsible = TRUE,
            DT::dataTableOutput(ns("id_table"))),
        box(width = NULL, title = "Sample Structure (PCA)", status = "primary", solidHeader = FALSE, collapsible = TRUE,
            uiOutput(ns("pca_ui"))),
        box(width = NULL, title = "Cross-Modality Correlation", status = "primary", solidHeader = FALSE, collapsible = TRUE,
            uiOutput(ns("correlation_ui")))
      )
    })
    output$overlap_table <- DT::renderDataTable({
      h <- req(harmonization_result())
      m <- req(h$overlap_matrix)
      df <- as.data.frame(m)
      df <- cbind(Modality = rownames(df), df)
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })
    overlap_heatmap_fn <- reactive({ h <- req(harmonization_result()); ch_overlap_heatmap_plot(h$overlap_matrix) })
    output$overlap_heatmap <- renderPlot(overlap_heatmap_fn(), alt = "Pairwise sample overlap heatmap")
    output$dl_overlap_heatmap_png <- multi_png_download(overlap_heatmap_fn, function() "cohort_harmonization_sample_overlap.png")
    output$id_table <- DT::renderDataTable({
      h <- req(harmonization_result())
      df <- req(h$id_table)
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    output$pca_ui <- renderUI({
      h <- req(harmonization_result())
      rd <- ov_raw_dataset()
      live_names <- if (isTRUE(rd$ok)) names(rd$layers) else character(0)
      if (length(live_names) == 0) return(div(class = "empty-note", icon("circle-info"), if (!isTRUE(rd$ok)) rd$error else "Insufficient information - no per-sample matrix is available."))
      ## "Color PCA by" offers every classified metadata column for this raw
      ## dataset (not just batch-name-pattern matches, which was this panel's
      ## only color option before) - defaults to the same suggested variable
      ## the Filters box uses, fully overridable.
      meta_cls <- ch_classify_metadata_columns(rd$sample_meta)
      color_choices <- meta_cls$table$column %||% character(0)
      tagList(
        if (!is.null(rd$provenance)) div(class = "empty-note", icon("circle-info"), rd$provenance),
        selectInput(ns("pca_layer"), "Modality", choices = live_names),
        if (length(color_choices) > 0) selectInput(ns("pca_color_by"), "Color PCA by", choices = c("(none)" = "", color_choices),
                                                     selected = meta_cls$suggested_default %||% ""),
        multi_plot_or_empty(pca_plot_fn, ns("pca_plot"), height = "340px"),
        div(class = "table-toolbar", downloadButton(ns("dl_pca_png"), "Download plot (PNG)", class = "btn-sm"))
      )
    })
    outputOptions(output, "pca_ui", suspendWhenHidden = FALSE)
    pca_plot_fn <- reactive({
      req(input$pca_layer)
      rd <- req(ov_raw_dataset())
      mat <- rd$layers[[input$pca_layer]]
      req(mat)
      multi_live_pca_plot(multi_live_pca(mat), rd$sample_meta, if (nzchar(input$pca_color_by %||% "")) input$pca_color_by else NULL)
    })
    output$pca_plot <- renderPlot(pca_plot_fn(), alt = "PCA of the selected modality, colored by the selected metadata column when available")
    output$dl_pca_png <- multi_png_download(pca_plot_fn, function() sprintf("cohort_harmonization_pca_%s.png", make.names(input$pca_layer %||% "modality")))

    output$correlation_ui <- renderUI({
      h <- req(harmonization_result())
      rd <- ov_raw_dataset()
      live_names <- if (isTRUE(rd$ok)) names(rd$layers) else character(0)
      if (length(live_names) < 2) return(div(class = "empty-note", icon("circle-info"), if (!isTRUE(rd$ok)) rd$error else "Insufficient information - correlation requires at least two modalities with a raw matrix."))
      tagList(
        fluidRow(
          column(6, selectInput(ns("corr_a"), "Modality A", choices = live_names)),
          column(6, selectInput(ns("corr_b"), "Modality B", choices = live_names, selected = live_names[min(2, length(live_names))]))
        ),
        multi_plot_or_empty(corr_plot_fn, ns("corr_plot"), height = "340px"),
        div(class = "table-toolbar", downloadButton(ns("dl_corr_png"), "Download plot (PNG)", class = "btn-sm"))
      )
    })
    outputOptions(output, "correlation_ui", suspendWhenHidden = FALSE)
    corr_plot_fn <- reactive({
      req(input$corr_a, input$corr_b)
      validate(need(!identical(input$corr_a, input$corr_b), "Choose two different modalities."))
      rd <- req(ov_raw_dataset())
      mA <- rd$layers[[input$corr_a]]; mB <- rd$layers[[input$corr_b]]
      req(mA, mB)
      d <- multi_live_correlation_heatmap_data(mA, mB, top_n = 20)
      req(isTRUE(d$ok))
      multi_live_correlation_heatmap_plot(d$df)
    })
    output$corr_plot <- renderPlot(corr_plot_fn(), alt = "Cross-modality feature correlation heatmap")
    output$dl_corr_png <- multi_png_download(corr_plot_fn, function() sprintf("cohort_harmonization_correlation_%s_vs_%s.png", make.names(input$corr_a %||% "A"), make.names(input$corr_b %||% "B")))

    ## =========================================================================
    ## Sample Explorer - browse/search every sample individually (which
    ## modalities it's in, its metadata) and see where one selected sample
    ## sits on a PCA of any modality with a raw matrix.
    ## =========================================================================
    output$sample_explorer_ui <- renderUI({
      h <- tryCatch(harmonization_result(), error = function(e) NULL)
      if (is.null(h) || !isTRUE(h$ok)) return(multi_empty_state("Click \"Analyze Cohort\" to see results here."))
      all_ids <- sort(unique(unlist(h$ids)))
      rd <- ov_raw_dataset()
      live_names <- if (isTRUE(rd$ok)) names(rd$layers) else character(0)
      tagList(
        box(width = NULL, title = "All Samples", status = "primary", solidHeader = FALSE,
            div(class = "table-toolbar", downloadButton(ns("dl_sample_table_csv"), "Download (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("sample_master_table"))),
        box(width = NULL, title = "Sample Profile", status = "primary", solidHeader = FALSE,
            selectizeInput(ns("explore_sample"), "Select a sample", choices = all_ids, options = list(placeholder = "Type to search...")),
            uiOutput(ns("sample_profile_ui")),
            if (length(live_names) > 0) tagList(
              selectInput(ns("explore_pca_layer"), "Highlight in PCA (modality)", choices = live_names),
              if (identical(multi_dataset$source, "preloaded")) div(class = "empty-note", icon("circle-info"), "Uses the selected analysis cell's matched-sample subset - the sample must be in that cell to appear."),
              multi_plot_or_empty(explore_pca_fn, ns("explore_pca_plot"), "PCA needs at least 3 samples and 2 features.", height = "340px"),
              div(class = "table-toolbar", downloadButton(ns("dl_explore_pca_png"), "Download plot (PNG)", class = "btn-sm"))
            ) else div(class = "empty-note", icon("circle-info"), if (!isTRUE(rd$ok)) rd$error else "Insufficient information - no per-sample matrix is available."))
      )
    })
    output$sample_master_table <- DT::renderDataTable({
      h <- req(harmonization_result())
      df <- req(ch_sample_master_table(h$ids, multi_dataset$sample_meta))
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })
    output$dl_sample_table_csv <- downloadHandler(
      filename = function() "cohort_harmonization_sample_table.csv",
      content = function(file) {
        h <- harmonization_result()
        utils::write.csv(ch_sample_master_table(h$ids, multi_dataset$sample_meta), file, row.names = FALSE)
      }
    )
    output$sample_profile_ui <- renderUI({
      req(input$explore_sample)
      h <- req(harmonization_result())
      present_in <- names(h$ids)[vapply(h$ids, function(x) input$explore_sample %in% x, logical(1))]
      meta_df <- multi_dataset$sample_meta
      meta_row <- if (!is.null(meta_df) && input$explore_sample %in% rownames(meta_df)) meta_df[input$explore_sample, , drop = FALSE] else NULL
      tagList(
        p(class = "submodule-desc", tags$strong("Present in: "), if (length(present_in) > 0) paste(present_in, collapse = ", ") else "None"),
        if (!is.null(meta_row)) tags$table(class = "table table-condensed", style = "font-size:0.88em;",
            tags$tbody(lapply(colnames(meta_row), function(cn) tags$tr(tags$td(tags$strong(cn)), tags$td(as.character(meta_row[[cn]]))))))
        else div(class = "empty-note", icon("circle-info"), "No metadata available for this sample.")
      )
    })
    explore_pca_fn <- reactive({
      req(input$explore_pca_layer, input$explore_sample)
      rd <- req(ov_raw_dataset())
      mat <- rd$layers[[input$explore_pca_layer]]
      req(mat)
      pca <- multi_live_pca(mat)
      req(isTRUE(pca$ok))
      ch_sample_highlight_pca_plot(pca, input$explore_sample)
    })
    output$explore_pca_plot <- renderPlot(explore_pca_fn(), alt = "PCA with the selected sample highlighted")
    output$dl_explore_pca_png <- multi_png_download(explore_pca_fn, function() sprintf("cohort_harmonization_sample_highlight_%s.png", make.names(input$explore_sample %||% "sample")))

    ## Backward-compatible publish: multi_qc_scorecard()/
    ## multi_analysis_summary_table() (multiomics_helpers.R) read
    ## overview$harmonization$ok/n_total/n_matched and
    ## overview$summary36$excludes_chance - field names kept identical so
    ## those two cross-cutting functions need no changes.
    observe({
      h <- tryCatch(harmonization_result(), error = function(e) NULL)
      if (is.null(h) || !isTRUE(h$ok) || is.null(multi_results)) return()
      summary36 <- NULL
      if (identical(multi_dataset$source, "preloaded")) {
        s36 <- multi_read_registry_table("Master six-part summary (integrated vs single-omics)")
        if (s36$ok) summary36 <- s36$df
      }
      multi_results$overview <- list(
        harmonization = list(ok = TRUE, n_total = h$n_total, n_matched = h$n_matched),
        summary36 = summary36, cells = h$cells
      )
    })
  })
}
