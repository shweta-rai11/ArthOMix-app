## Overview and Datasets submodule: GEO source catalog, metadata/expression
## browsers, and QC (missing values, outlier detection, normalisation check,
## group filtering) for the selected raw dataset.

mod_overview_config <- list(
  id = "overview", group = "Data",
  title = "Overview and Datasets",
  description = "Quality control for the currently loaded dataset - preloaded or uploaded.",
  icon = "table-cells"
)

mod_overview_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    id = ns("tabs"), type = "tabs",
    tabPanel(
      "Datasets", br(),
      p(class = "submodule-desc", "The NCBI GEO series accession(s) behind whichever dataset is currently active on the Dataset tab - preloaded, GEO-fetched, or uploaded."),
      withSpinner(uiOutput(ns("sources_ui")), color = "#2c6fbb", type = 6)
    ),
    tabPanel(
      "Metadata", br(),
      box(
        width = 12, status = "primary", solidHeader = FALSE,
        fluidRow(
          column(5, uiOutput(ns("qc_source_ui_meta"))),
          column(7, uiOutput(ns("qc_source_info_ui_meta")))
        )
      ),
      p(class = "submodule-desc", "What's currently selected above, as a whole, before any filtering."),
      withSpinner(uiOutput(ns("understand_ui")), color = "#2c6fbb", type = 6),
      box(
        width = 12, title = "Sample metadata", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Every column is sortable and filterable (click a header, or use the box underneath it)."),
        div(class = "table-toolbar", downloadButton(ns("download_meta_full"), "Download CSV", class = "btn-sm")),
        withSpinner(DT::dataTableOutput(ns("meta_table_full")), color = "#2c6fbb", type = 6)
      )
    ),
    tabPanel(
      "Expression data", br(),
      box(
        width = 12, status = "primary", solidHeader = FALSE,
        fluidRow(
          column(5, uiOutput(ns("qc_source_ui_expr"))),
          column(7, uiOutput(ns("qc_source_info_ui_expr")))
        )
      ),
      p(class = "submodule-desc desc-aligned", "The actual expression matrix behind every check below - features as rows, samples as columns. Search for a gene/probe, sort by any sample, or download."),
      box(
        width = 12, title = "Expression matrix", status = "primary", solidHeader = FALSE,
        div(class = "table-toolbar", downloadButton(ns("download_expr"), "Download CSV", class = "btn-sm")),
        withSpinner(DT::dataTableOutput(ns("expr_table")), color = "#2c6fbb", type = 6)
      )
    ),
    tabPanel(
      "QC", br(),
      box(
        width = 12, status = "primary", solidHeader = FALSE,
        fluidRow(
          column(5, uiOutput(ns("qc_source_ui"))),
          column(7, uiOutput(ns("qc_source_info_ui")))
        )
      ),
      tabsetPanel(
        id = ns("qc_tabs"), type = "tabs",
        tabPanel(
          "Missing values", br(),
          p(class = "submodule-desc", "Percent missing per metadata field, across every sample in the selected dataset - the same audit run before anything else."),
          box(
            width = 12, title = "Percent missing by field", status = "primary", solidHeader = FALSE,
            withSpinner(plotOutput(ns("missing_plot"), height = 320), color = "#2c6fbb", type = 6)
          ),
          box(
            width = 12, title = "By field", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("missing_table"))
          )
        ),
        tabPanel(
          "Outliers", br(),
          p(class = "submodule-desc", "Signal, detected features and correlation to the rest of the cohort, per sample - the same checks used in Preprocessing."),
          fluidRow(
            column(
              4,
              sliderInput(ns("mad_k"), "Outlier sensitivity (MADs from the cohort median)", min = 2, max = 6, value = 3, step = 0.5),
              actionButton(ns("run_qc_btn"), "Run outlier detection", icon = icon("play"), class = "btn-primary btn-sm")
            ),
            column(8, withSpinner(uiOutput(ns("qc_summary_ui")), color = "#2c6fbb", type = 6))
          ),
          withSpinner(uiOutput(ns("qc_plots_ui")), color = "#2c6fbb", type = 6)
        ),
        tabPanel(
          "Normalised data", br(),
          p(class = "submodule-desc", "Checks whether samples are on a comparable scale: per-sample distribution, a spread diagnostic, a scree plot and PCA - the same checks used to decide on quantile normalisation in Preprocessing."),
          fluidRow(
            column(
              4,
              uiOutput(ns("norm_color_by_ui")),
              fluidRow(
                column(6, selectInput(ns("norm_pc_x"), "X axis", choices = setNames(1:5, paste0("PC", 1:5)), selected = 1, selectize = FALSE)),
                column(6, selectInput(ns("norm_pc_y"), "Y axis", choices = setNames(1:5, paste0("PC", 1:5)), selected = 2, selectize = FALSE))
              ),
              checkboxInput(ns("norm_show_ellipse"), "Show group confidence ellipses", value = TRUE),
              checkboxInput(ns("norm_show_labels"), "Label points with sample ID", value = FALSE),
              actionButton(ns("run_norm_btn"), "Run normalisation check", icon = icon("play"), class = "btn-primary btn-sm")
            ),
            column(8, withSpinner(uiOutput(ns("norm_summary_ui")), color = "#2c6fbb", type = 6))
          ),
          withSpinner(uiOutput(ns("norm_views_ui")), color = "#2c6fbb", type = 6),
          uiOutput(ns("norm_apply_ui"))
        ),
        tabPanel(
          "Group", br(),
          p(class = "submodule-desc", "Pick which samples to look at, then apply - composition and the sample table below reflect your selection."),
          fluidRow(
            column(4, uiOutput(ns("filters"))),
            column(8, uiOutput(ns("filter_summary_ui")))
          ),
          uiOutput(ns("filtered_views_ui"))
        )
      )
    )
  )
}

mod_overview_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ---- Datasets tab -----------------------------------------------------
    ## Always visible; content follows whichever pipeline is currently active
    ## (dataset$geo_ids, set alongside source_type wherever the Dataset tab's
    ## three pipelines write dataset$source - see mod_dataset.R). Shows the
    ## real NCBI GEO accession(s) behind the current dataset when there are
    ## any (every preloaded pick, any successful GEO fetch), or a clear
    ## "no GEO ID" message otherwise (the common case for an upload, which
    ## has no known GEO provenance) - never the old fixed 4-source catalog
    ## regardless of what's actually loaded.

    output$sources_ui <- renderUI({
      ids <- dataset$geo_ids %||% character(0)
      if (length(ids) == 0) {
        return(div(class = "empty-note", icon("circle-info"),
                    "No GEO ID found for the currently loaded dataset - it wasn't fetched from, or matched to, an NCBI GEO series."))
      }
      cards <- lapply(ids, function(gse_id) {
        src <- Find(function(s) identical(s$gse, gse_id), GEO_SOURCES)
        eset <- get_raw_eset(gse_id)
        div(
          class = "info-card",
          div(
            class = "module-card-title-row",
            h4(gse_id),
            tags$a(href = geo_link(gse_id), target = "_blank", rel = "noopener",
                    icon("up-right-from-square"), " NCBI GEO")
          ),
          if (!is.null(eset)) {
            tagList(
              p(class = "module-card-tagline",
                tryCatch(Biobase::experimentData(eset)@title, error = function(e) NULL)),
              if (!is.null(src)) p(strong("Role: "), src$role, br(), strong("Used for: "), src$used_in),
              p(strong("Platform: "), Biobase::annotation(eset), br(),
                strong("Samples: "), ncol(eset), ", ", strong("Probes: "), format(nrow(eset), big.mark = ","))
            )
          } else if (!is.null(src)) {
            tagList(
              p(strong("Role: "), src$role, br(), strong("Used for: "), src$used_in),
              div(class = "empty-note", icon("triangle-exclamation"), "Raw file not found on disk.")
            )
          } else {
            p(class = "empty-note", icon("circle-info"), "Fetched live from NCBI GEO for the currently loaded dataset.")
          }
        )
      })
      div(class = "module-grid", cards)
    })

    ## QC tab source picker: raw individual datasets only (no merge/ComBat),
    ## plus the user's own upload if present.

    ## Isolation: once an upload or GEO fetch is the active pipeline, these
    ## pickers collapse to that one active dataset only - no preloaded/GEO
    ## catalog choices at all. The preloaded pipeline keeps today's behavior
    ## (browse any of the 4 fixed reference sources) unchanged.
    qc_source_choices <- reactive({
      if (identical(dataset$source_type, "uploaded") || identical(dataset$source_type, "geo")) {
        label <- if (identical(dataset$source_type, "geo")) {
          paste0("Your GEO-fetched data (", dataset$source, ")")
        } else {
          paste0("Your uploaded data (", dataset$source, ")")
        }
        return(setNames("active", label))
      }
      setNames(vapply(GEO_SOURCES, `[[`, character(1), "gse"),
               vapply(GEO_SOURCES, function(s) sprintf("%s (%s, raw)", s$gse, s$role), character(1)))
    })

    qc_source_default <- reactive({
      if (identical(dataset$source_type, "uploaded") || identical(dataset$source_type, "geo")) "active" else GEO_SOURCES[[1]]$gse
    })

    resolve_qc_source <- function(source_id) {
      if (identical(source_id, "active")) {
        req(dataset$expr, dataset$meta)
        return(list(expr = dataset$expr, meta = dataset$meta, label = dataset$source %||% "Currently loaded dataset"))
      }
      d <- load_individual_dataset(source_id)
      validate(need(!is.null(d), paste("Raw data for", source_id, "was not found on disk.")))
      d
    }

    ## Each tab gets its own picker instead of sharing one, since a shared
    ## conditionalPanel picker doesn't reliably bind once inserted via insertTab.
    output$qc_source_ui_meta <- renderUI({
      selectInput(ns("qc_source_meta"), "Dataset to inspect", choices = qc_source_choices(),
                  selected = qc_source_default(), width = "100%")
    })
    qc_target_meta <- reactive({ req(input$qc_source_meta); resolve_qc_source(input$qc_source_meta) })
    output$qc_source_info_ui_meta <- renderUI({
      t <- qc_target_meta()
      div(class = "empty-note", icon("circle-info"), strong(t$label), " - ",
          format(nrow(t$meta), big.mark = ","), " samples, ",
          format(nrow(t$expr), big.mark = ","), " features.")
    })

    output$qc_source_ui_expr <- renderUI({
      selectInput(ns("qc_source_expr"), "Dataset to inspect", choices = qc_source_choices(),
                  selected = qc_source_default(), width = "100%")
    })
    qc_target_expr <- reactive({ req(input$qc_source_expr); resolve_qc_source(input$qc_source_expr) })
    output$qc_source_info_ui_expr <- renderUI({
      t <- qc_target_expr()
      div(class = "empty-note", icon("circle-info"), strong(t$label), " - ",
          format(nrow(t$meta), big.mark = ","), " samples, ",
          format(nrow(t$expr), big.mark = ","), " features.")
    })

    output$qc_source_ui <- renderUI({
      selectInput(ns("qc_source"), "Dataset to inspect", choices = qc_source_choices(),
                  selected = qc_source_default(), width = "100%")
    })
    qc_target <- reactive({ req(input$qc_source); resolve_qc_source(input$qc_source) })

    ## The three on-demand panels below cache an eventReactive keyed to a click
    ## count, which can't be reset server-side - so a stale flag per panel marks
    ## their cached results as belonging to a previous dataset until re-run.
    qc_stale <- reactiveVal(FALSE)
    norm_stale <- reactiveVal(FALSE)
    filter_stale <- reactiveVal(FALSE)
    observeEvent(qc_target(), {
      qc_stale(TRUE); norm_stale(TRUE); filter_stale(TRUE)
    }, ignoreInit = TRUE)
    observeEvent(input$run_qc_btn, qc_stale(FALSE), ignoreInit = TRUE)
    observeEvent(input$run_norm_btn, norm_stale(FALSE), ignoreInit = TRUE)
    observeEvent(input$apply_btn, filter_stale(FALSE), ignoreInit = TRUE)

    output$qc_source_info_ui <- renderUI({
      t <- qc_target()
      div(class = "empty-note", icon("circle-info"), strong(t$label), " - ",
          format(nrow(t$meta), big.mark = ","), " samples, ",
          format(nrow(t$expr), big.mark = ","), " features.")
    })

    ## ---- Metadata tab: overview of the selected dataset, unfiltered -------

    output$understand_ui <- renderUI({
      t <- qc_target_meta()
      meta <- t$meta
      tagList(
        p(strong("Source: "), t$label),
        fluidRow(
          valueBox(nrow(meta), "Samples", icon = icon("users"), color = "light-blue", width = 3),
          valueBox(length(unique(na.omit(meta$group))), "Groups", icon = icon("layer-group"), color = "purple", width = 3),
          valueBox(format(nrow(t$expr), big.mark = ","), "Features in matrix", icon = icon("dna"), color = "green", width = 3),
          valueBox(length(unique(na.omit(meta$sex))), "Sex categories", icon = icon("venus-mars"), color = "red", width = 3)
        ),
        p(class = "submodule-desc",
          "Next: browse the full metadata table below or the Expression data tab for the matrix itself, check Missing values and Outliers for data quality issues, Normalised data for whether samples are on a comparable scale, or filter by any metadata column to explore a subset.")
      )
    })

    output$meta_table_full <- DT::renderDataTable({
      DT::datatable(qc_target_meta()$meta, rownames = FALSE, filter = "top",
                     options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_meta_full <- downloadHandler(
      filename = function() paste0(qc_target_meta()$meta$dataset[1] %||% "metadata", "_metadata.csv"),
      content = function(file) write.csv(qc_target_meta()$meta, file, row.names = FALSE)
    )

    ## ---- QC tab: missing-value audit, whole dataset, always on -----------

    missing_audit <- reactive({
      df <- qc_target()$meta
      miss <- data.frame(
        field = colnames(df),
        n_missing = vapply(df, function(x) sum(is.na(x) | x %in% c("", "NA", "N/A", "unknown")), integer(1)),
        stringsAsFactors = FALSE
      )
      miss$pct_missing <- round(100 * miss$n_missing / nrow(df), 1)
      miss$status <- ifelse(miss$n_missing == 0, "Complete", "Has missing")
      miss[order(miss$pct_missing), ]
    })

    output$missing_plot <- renderPlot({
      miss <- missing_audit()
      ggplot(miss, aes(x = reorder(field, pct_missing), y = pct_missing, fill = status)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = c("Complete" = ARTHOMIX_STATUS$good, "Has missing" = ARTHOMIX_STATUS$critical)) +
        coord_flip() +
        labs(x = NULL, y = "% missing", fill = NULL) +
        theme_arthomix(base_size = 13) +
        theme(panel.grid.major.y = element_blank())
    })

    output$missing_table <- DT::renderDataTable({
      DT::datatable(missing_audit()[, c("field", "n_missing", "pct_missing", "status")],
                     rownames = FALSE, options = list(pageLength = 8, dom = "tp", scrollX = TRUE), class = "stripe hover compact")
    })

    ## ---- QC tab: sample-level outlier detection, run on demand -----------

    sample_qc <- eventReactive(input$run_qc_btn, {
      t <- qc_target()
      qc <- compute_sample_qc(t$expr, mad_k = input$mad_k)
      merge(qc, t$meta[, intersect(c("sample", "group"), colnames(t$meta))], by = "sample", all.x = TRUE)
    })

    output$qc_summary_ui <- renderUI({
      if (!isTruthy(input$run_qc_btn) || input$run_qc_btn == 0) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet - click Run outlier detection to check for technical outliers."))
      }
      if (qc_stale()) {
        return(div(class = "empty-note", icon("circle-info"), "Dataset changed - click Run outlier detection again to refresh."))
      }
      qc <- sample_qc()
      n_flagged <- sum(qc$flag_signal | qc$flag_detected | qc$flag_cor)
      valueBox(n_flagged, "Samples flagged", icon = icon("triangle-exclamation"),
                color = if (n_flagged > 0) "red" else "green", width = 12)
    })

    output$qc_plots_ui <- renderUI({
      if (qc_stale()) return(NULL)
      req(sample_qc())
      tagList(
        fluidRow(
          column(4, box(width = NULL, title = "Per-sample signal", status = "primary", solidHeader = FALSE,
                          plotOutput(ns("signal_plot"), height = 200))),
          column(4, box(width = NULL, title = "Detected features", status = "primary", solidHeader = FALSE,
                          plotOutput(ns("detected_plot"), height = 200))),
          column(4, box(width = NULL, title = "Mean correlation to cohort", status = "primary", solidHeader = FALSE,
                          plotOutput(ns("cor_plot"), height = 200)))
        ),
        box(
          width = 12, title = "Flagged samples", status = "primary", solidHeader = FALSE,
          div(class = "table-toolbar", downloadButton(ns("download_qc"), "Download full QC table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("qc_table"))
        )
      )
    })

    output$signal_plot   <- renderPlot(qc_bar_plot(sample_qc(), "signal", "flag_signal", "Total signal"))
    output$detected_plot <- renderPlot(qc_bar_plot(sample_qc(), "detected", "flag_detected", "Detected features"))
    output$cor_plot       <- renderPlot(qc_bar_plot(sample_qc(), "mean_cor", "flag_cor", "Mean correlation"))

    qc_table_display <- reactive({
      df <- sample_qc()
      df$reason <- apply(df[, c("flag_signal", "flag_detected", "flag_cor")], 1, function(r) {
        reasons <- c("low/high signal", "low/high detected features", "low cohort correlation")[r]
        if (length(reasons) == 0) "" else paste(reasons, collapse = "; ")
      })
      df$signal <- round(df$signal, 1); df$mean_cor <- round(df$mean_cor, 3)
      df[, c("sample", intersect("group", colnames(df)), "signal", "detected", "mean_cor", "reason")]
    })

    output$qc_table <- DT::renderDataTable({
      df <- qc_table_display()
      flagged_only <- df[df$reason != "", , drop = FALSE]
      DT::datatable(if (nrow(flagged_only) > 0) flagged_only else df[0, ], rownames = FALSE, filter = "top",
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_qc <- downloadHandler(
      filename = function() "qc_metrics.csv",
      content = function(file) write.csv(qc_table_display(), file, row.names = FALSE)
    )

    ## ---- QC tab: normalisation check, run on demand -----------------------

    ## Rendered (not updateSelectInput) so choices populate once the tab is
    ## actually visible, since the server starts before the tab exists in the DOM.
    output$norm_color_by_ui <- renderUI({
      cols <- setdiff(colnames(qc_target()$meta), "sample")
      selectInput(ns("norm_color_by"), "Color by", choices = cols,
                  selected = if ("group" %in% cols) "group" else cols[1])
    })

    norm_check <- eventReactive(input$run_norm_btn, {
      expr <- qc_target()$expr
      diag <- summarize_norm_diagnostics(expr)
      qs <- apply(expr, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      box_df <- data.frame(sample = colnames(expr), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ],
                            upper = qs[4, ], ymax = qs[5, ])
      persample <- data.frame(sample = colnames(expr), median = round(qs[3, ], 2),
                                q25 = round(qs[2, ], 2), q75 = round(qs[4, ], 2))
      list(diag = diag, box_df = box_df, persample = persample, pca = pca_of(expr))
    })

    output$norm_summary_ui <- renderUI({
      if (!isTruthy(input$run_norm_btn) || input$run_norm_btn == 0) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet - click Run normalisation check to see whether samples are on a comparable scale."))
      }
      if (norm_stale()) {
        return(div(class = "empty-note", icon("circle-info"), "Dataset changed - click Run normalisation check again to refresh."))
      }
      d <- norm_check()$diag
      needs <- needs_quantile_norm(d)
      tagList(
        p(icon("circle-info"), " Spread of per-sample medians (SD): ", strong(sprintf("%.3f", d$median_sd)),
          "; spread of per-sample IQRs (SD): ", strong(sprintf("%.3f", d$iqr_sd)),
          "; max value in matrix: ", strong(format(round(d$max_value, 1), big.mark = ",")), "."),
        p(class = "empty-note", icon(if (!needs) "check" else "triangle-exclamation"),
          if (!needs) "Samples look well aligned and log-scaled - quantile normalisation likely isn't needed."
          else "This looks like it needs quantile normalisation: either samples disagree by more than 0.5 on the log scale, or values are still on a linear (not log2) scale - the same rule Preprocessing uses to decide. See Normalise this dataset below.")
      )
    })

    output$norm_views_ui <- renderUI({
      if (norm_stale()) return(NULL)
      req(norm_check())
      tagList(
        box(
          width = 12, title = "Per-sample expression distribution", status = "primary", solidHeader = FALSE,
          plotOutput(ns("norm_dist_plot"), height = 280)
        ),
        box(
          width = 12, title = "Scree plot", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "% variance explained per component."),
          plotOutput(ns("norm_scree_plot"), height = 220)
        ),
        fluidRow(
          column(7, box(width = NULL, title = "PCA", status = "primary", solidHeader = FALSE,
                          plotOutput(ns("norm_pca_plot"), height = 340))),
          column(5, box(width = NULL, title = "Per-sample summary", status = "primary", solidHeader = FALSE,
                          DT::dataTableOutput(ns("norm_table"))))
        )
      )
    })

    output$norm_dist_plot <- renderPlot({
      req(norm_check(), input$norm_color_by)
      box_df <- norm_check()$box_df
      color_col <- input$norm_color_by
      meta <- qc_target()$meta
      join_cols <- intersect(c("sample", color_col), colnames(meta))
      box_df <- merge(box_df, meta[, join_cols, drop = FALSE], by = "sample", all.x = TRUE)
      if (!color_col %in% colnames(box_df)) { box_df[[color_col]] <- "all samples"; }
      ggplot(box_df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle,
                           upper = upper, ymax = ymax, fill = .data[[color_col]])) +
        geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15) +
        scale_fill_manual(values = arthomix_pair(box_df[[color_col]])) +
        labs(x = NULL, y = "Expression", fill = NULL) +
        theme_arthomix() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
    })

    output$norm_scree_plot <- renderPlot({
      scree_plot(norm_check()$pca$var_exp)
    })

    output$norm_pca_plot <- renderPlot({
      req(norm_check(), input$norm_color_by)
      plot_pca_advanced(norm_check()$pca, qc_target()$meta, input$norm_color_by,
                          pc_x = as.integer(input$norm_pc_x), pc_y = as.integer(input$norm_pc_y),
                          show_ellipse = input$norm_show_ellipse, show_labels = input$norm_show_labels)
    })

    output$norm_table <- DT::renderDataTable({
      DT::datatable(norm_check()$persample, rownames = FALSE,
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    ## Live quantile normalisation preview; only adoptable app-wide for the
    ## user's own upload, since the fixed GEO sources are read-only reference data.

    norm_apply_result <- eventReactive(input$apply_norm_btn, {
      expr <- as.matrix(qc_target()$expr)
      normalized <- limma::normalizeBetweenArrays(expr, method = "quantile")
      list(
        diag_before = summarize_norm_diagnostics(expr), diag_after = summarize_norm_diagnostics(normalized),
        box_before = norm_check()$box_df,
        box_after = {
          qs <- apply(normalized, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
          data.frame(sample = colnames(normalized), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ], upper = qs[4, ], ymax = qs[5, ])
        },
        expr_after = normalized
      )
    })

    output$norm_apply_ui <- renderUI({
      if (norm_stale()) return(NULL)
      req(norm_check())
      needs <- needs_quantile_norm(norm_check()$diag)
      tagList(
        box(
          width = 12, title = "Normalise this dataset", status = if (needs) "warning" else "primary", solidHeader = FALSE,
          p(class = "submodule-desc",
            "Runs the same quantile normalisation (limma::normalizeBetweenArrays) Preprocessing applies to the merged training cohort, live, on whatever's selected above - so if you uploaded your own, un-normalised data, you can check it and fix it right here."),
          actionButton(ns("apply_norm_btn"), "Apply quantile normalisation", icon = icon("wand-magic-sparkles"), class = "btn-primary btn-sm"),
          uiOutput(ns("norm_apply_result_ui"))
        )
      )
    })

    output$norm_apply_result_ui <- renderUI({
      req(norm_apply_result())
      r <- norm_apply_result()
      tagList(
        br(),
        fluidRow(
          valueBox(sprintf("%.3f -> %.3f", r$diag_before$median_sd, r$diag_after$median_sd), "Spread of medians (SD)",
                    icon = icon("chart-line"), color = if (r$diag_after$median_sd < 0.5) "green" else "red", width = 6),
          valueBox(sprintf("%.3f -> %.3f", r$diag_before$iqr_sd, r$diag_after$iqr_sd), "Spread of IQRs (SD)",
                    icon = icon("chart-line"), color = if (r$diag_after$iqr_sd < 0.5) "green" else "red", width = 6)
        ),
        fluidRow(
          column(6, box(width = NULL, title = "Before", status = "primary", solidHeader = FALSE,
                          plotOutput(ns("norm_before_plot"), height = 260))),
          column(6, box(width = NULL, title = "After", status = "primary", solidHeader = FALSE,
                          plotOutput(ns("norm_after_plot"), height = 260)))
        ),
        if (identical(input$qc_source, "active")) {
          div(
            actionButton(ns("adopt_norm_btn"), "Use this normalised version for every sub-module", icon = icon("check"), class = "btn-success btn-sm"),
            uiOutput(ns("adopt_norm_msg"))
          )
        } else {
          p(class = "empty-note", icon("circle-info"),
            "This is one of the app's fixed reference datasets, so it stays read-only here - only its diagnostics change, not the data every sub-module reads. Upload your own data, or fetch from NCBI GEO, on the Dataset tab to normalise it and use the result app-wide.")
        }
      )
    })

    norm_box_plot <- function(box_df) {
      ggplot(box_df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax)) +
        geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15, fill = ARTHOMIX_COLORS$blue) +
        labs(x = NULL, y = "Expression") +
        theme_arthomix() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
    }
    output$norm_before_plot <- renderPlot(norm_box_plot(norm_apply_result()$box_before))
    output$norm_after_plot  <- renderPlot(norm_box_plot(norm_apply_result()$box_after))

    observeEvent(input$adopt_norm_btn, {
      req(norm_apply_result(), identical(input$qc_source, "active"))
      dataset$expr <- norm_apply_result()$expr_after
      dataset$source <- paste0(dataset$source, " (quantile-normalised)")
      output$adopt_norm_msg <- renderUI(
        div(class = "empty-note", icon("check"), "Done - every sub-module now reads the quantile-normalised version of your data.")
      )
    })

    ## Server-side DT paging (server = TRUE below) since these matrices run to
    ## tens of thousands of rows.

    expr_table_data <- reactive({
      m <- qc_target_expr()$expr
      data.frame(feature = rownames(m), round(m, 3), check.names = FALSE, stringsAsFactors = FALSE)
    })

    output$expr_table <- DT::renderDataTable({
      df <- expr_table_data()
      DT::datatable(
        df, rownames = FALSE, extensions = "FixedColumns",
        options = list(pageLength = 15, dom = "lfrtip", scrollX = TRUE, fixedColumns = list(leftColumns = 1)),
        class = "stripe hover compact"
      )
    }, server = TRUE)

    output$download_expr <- downloadHandler(
      filename = function() "expression_matrix.csv",
      content = function(file) data.table::fwrite(expr_table_data(), file)
    )

    ## Builds filter widgets from whatever metadata columns exist: numeric ->
    ## range slider, low-cardinality categorical -> multi-select, high-cardinality skipped.

    filter_spec <- reactive({
      meta <- qc_target()$meta
      cols <- setdiff(colnames(meta), "sample")
      specs <- list()
      for (cl in cols) {
        x <- meta[[cl]]
        if (is.numeric(x)) {
          rng <- suppressWarnings(range(x, na.rm = TRUE))
          if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) {
            specs[[cl]] <- list(type = "numeric", min = floor(rng[1]), max = ceiling(rng[2]))
          }
        } else {
          u <- sort(unique(na.omit(as.character(x))))
          if (length(u) >= 2 && length(u) <= 30) {
            specs[[cl]] <- list(type = "categorical", choices = u)
          }
        }
      }
      specs
    })

    ## Filter widget input ids are built from each metadata column's position,
    ## never the raw column name - GEO's own raw pData column names routinely
    ## contain spaces and colons (e.g. "disease state:ch1", confirmed live off
    ## a real GSE93272 fetch), and Shiny's client dispatches messages keyed as
    ## "inputType:inputId" - a colon (or a space, which breaks the jQuery
    ## selector/DOM id) inside the id corrupts that routing and throws an
    ## uncaught "No handler registered for type ..." error that kills the
    ## whole session (confirmed live: the entire page goes grey/unresponsive -
    ## a Shiny disconnection, not a rendering glitch). names(specs) has a
    ## stable order across recomputation for a fixed active dataset, so the
    ## same column always maps back to the same positional id.
    filter_id <- function(specs, cl) paste0("f_", match(cl, names(specs)))

    output$filters <- renderUI({
      specs <- filter_spec()
      validate(need(length(specs) > 0, "No filterable columns in the current metadata."))
      tagList(
        lapply(names(specs), function(cl) {
          s <- specs[[cl]]
          fid <- filter_id(specs, cl)
          if (s$type == "categorical") {
            pickerInput(ns(fid), cl, choices = s$choices, selected = character(0),
                        multiple = TRUE, options = list(`actions-box` = TRUE, title = paste0("Choose ", cl, "... (optional)")))
          } else {
            sliderInput(ns(fid), cl, min = s$min, max = s$max, value = c(s$min, s$max))
          }
        }),
        fluidRow(
          column(6, actionButton(ns("apply_btn"), "Apply filters", icon = icon("filter"), class = "btn-primary btn-sm")),
          column(6, actionButton(ns("reset_btn"), "Reset", icon = icon("rotate-left"), class = "btn-default btn-sm"))
        )
      )
    })

    observeEvent(input$reset_btn, {
      specs <- filter_spec()
      for (cl in names(specs)) {
        s <- specs[[cl]]; fid <- filter_id(specs, cl)
        if (s$type == "categorical") updatePickerInput(session, fid, selected = character(0))
        else updateSliderInput(session, fid, value = c(s$min, s$max))
      }
    })

    filtered_meta <- eventReactive(input$apply_btn, {
      meta <- qc_target()$meta
      specs <- filter_spec()
      for (cl in names(specs)) {
        s <- specs[[cl]]; val <- input[[filter_id(specs, cl)]]
        if (s$type == "categorical") {
          if (length(val) > 0) meta <- meta[is.na(meta[[cl]]) | meta[[cl]] %in% val, , drop = FALSE]
        } else if (!is.null(val)) {
          meta <- meta[is.na(meta[[cl]]) | (meta[[cl]] >= val[1] & meta[[cl]] <= val[2]), , drop = FALSE]
        }
      }
      validate(need(nrow(meta) > 0, "No samples match the current filter combination - try Reset."))
      meta
    })

    output$filter_summary_ui <- renderUI({
      if (!isTruthy(input$apply_btn) || input$apply_btn == 0) {
        return(div(class = "empty-note", icon("circle-info"), "Set any filters on the left (all optional), then click Apply filters. Everything shows by default."))
      }
      if (filter_stale()) {
        return(div(class = "empty-note", icon("circle-info"), "Dataset changed - click Apply filters again to refresh."))
      }
      df <- filtered_meta()
      fluidRow(
        valueBox(nrow(df), "Samples selected", icon = icon("users"), color = "light-blue", width = 3),
        valueBox(length(unique(na.omit(df$group))), "Groups", icon = icon("layer-group"), color = "purple", width = 3),
        if ("dataset" %in% colnames(df)) valueBox(length(unique(na.omit(df$dataset))), "Datasets", icon = icon("database"), color = "green", width = 3),
        valueBox(length(unique(na.omit(df$sex))), "Sex categories", icon = icon("venus-mars"), color = "red", width = 3)
      )
    })

    output$filtered_views_ui <- renderUI({
      if (filter_stale()) return(NULL)
      req(filtered_meta())
      tagList(
        fluidRow(
          column(7, box(title = "Cohort composition", width = NULL, status = "primary", solidHeader = FALSE,
                          plotOutput(ns("group_plot"), height = 320))),
          if ("dataset" %in% colnames(filtered_meta())) {
            column(5, box(title = "By dataset", width = NULL, status = "primary", solidHeader = FALSE,
                            plotOutput(ns("dataset_plot"), height = 320)))
          }
        ),
        box(
          title = "Sample metadata", width = 12, status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Every column below is also individually sortable and filterable (click a header, or use the box underneath it)."),
          div(class = "table-toolbar", downloadButton(ns("download_meta"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("meta_table"))
        )
      )
    })

    output$group_plot <- renderPlot({
      df <- filtered_meta()
      counted <- df %>% count(group, sex)
      p <- ggplot(counted, aes(x = group, y = n, fill = group)) +
        geom_col(width = 0.6) +
        scale_fill_manual(values = arthomix_pair(counted$group)) +
        labs(x = NULL, y = "Samples", fill = "Group") +
        theme_arthomix(base_size = 13)
      if (length(unique(na.omit(df$sex))) > 0) p <- p + facet_wrap(~sex)
      p
    })

    output$dataset_plot <- renderPlot({
      df <- filtered_meta()
      req("dataset" %in% colnames(df))
      counted <- df %>% count(dataset)
      ggplot(counted, aes(x = reorder(dataset, n), y = n, fill = dataset)) +
        geom_col(width = 0.6) +
        scale_fill_manual(values = arthomix_pair(counted$dataset)) +
        coord_flip() +
        labs(x = NULL, y = "Samples", fill = NULL) +
        theme_arthomix(base_size = 13) +
        theme(legend.position = "none")
    })

    output$meta_table <- DT::renderDataTable({
      DT::datatable(filtered_meta(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_meta <- downloadHandler(
      filename = function() "sample_metadata_filtered.csv",
      content = function(file) write.csv(filtered_meta(), file, row.names = FALSE)
    )
  })
}
