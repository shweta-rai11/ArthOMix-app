## R/multiomics/mod_multi_live.R
## Submodule: Live Analysis (Upload & MOFA2) - the one part of the
## Multi-Omics module that computes on data the user supplies, rather than
## browsing Research_05_multiomics_sexstratified's precomputed results. This
## file covers sections 1-4 (Upload & Validate, Sample Matching & Missing
## Data, Normalization/Filtering/Scaling, Batch Diagnostics); section 5-8
## (MOFA2 integration, factor results, cross-omics correlation, export) are
## a nested sub-module (mod_multi_live_mofa.R) mounted at the bottom of this
## tab, reading the `live_state` reactiveValues this file publishes into.
##
## Up to 4 omics layers, each independently optional - supports the pipeline
## cohort's own transcriptomics+methylomics pair, and is extensible beyond
## it (proteomics/metabolomics/other), per spec section 34's "preloaded vs.
## uploaded must never contaminate each other": this tab never reads
## multi_dataset (the precomputed-cohort Dataset tab's own state), and
## nothing here writes back to it either - two fully independent pipelines.

mod_multi_live_config <- list(
  id = "live", title = "Live Analysis (Upload & MOFA2)", icon = "upload", group = "Live Analysis",
  description = "Upload your own omics matrices, validate and QC them, choose omics-appropriate normalization/filtering/scaling, check for batch effects, then run a real MOFA2 factor analysis."
)

MULTI_LIVE_LAYER_IDS <- c("layer1", "layer2", "layer3", "layer4")

mod_multi_live_layer_ui <- function(ns, i) {
  lid <- MULTI_LIVE_LAYER_IDS[i]
  box(
    width = NULL, title = sprintf("Layer %d%s", i, if (i <= 2) " (required)" else " (optional)"), status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = i > 2,
    selectInput(ns(paste0(lid, "_type")), "Omics type", choices = MULTI_LIVE_OMICS_TYPES, selected = if (i == 1) "rnaseq" else if (i == 2) "methylation" else "other"),
    textInput(ns(paste0(lid, "_label")), "Display label", value = sprintf("Layer %d", i)),
    fileInput(ns(paste0(lid, "_file")), "Matrix (CSV or RDS)", accept = c(".csv", ".rds", ".Rds")),
    radioButtons(ns(paste0(lid, "_orient")), "First column / matrix orientation",
                 choices = c("Rows = samples (first column = sample ID)" = "samples_rows", "Rows = features (first column = feature ID)" = "features_rows"),
                 selected = "samples_rows"),
    p(class = "submodule-desc", "Optional: rows must match sample IDs used in the metadata file below.")
  )
}

mod_multi_live_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "About Live Analysis", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "This tab is completely separate from the precomputed-cohort tabs above - nothing you upload here changes them, and nothing precomputed is read into this one."),
        p(class = "submodule-desc", "Matched-sample integration only: MOFA2 (below) trains on samples present in every uploaded layer. No live DIABLO-style supervised classifier is offered here - see the Summary tab's known limitations.")
      ),
      lapply(seq_along(MULTI_LIVE_LAYER_IDS), function(i) mod_multi_live_layer_ui(ns, i)),
      box(
        width = NULL, title = "Sample metadata (optional)", status = "primary", solidHeader = FALSE,
        fileInput(ns("meta_file"), "Metadata (CSV, first column = sample ID)", accept = c(".csv")),
        p(class = "submodule-desc", "Used to color PCA/factor plots and to run batch-effect diagnostics; not required to validate or normalize the matrices above.")
      ),
      actionButton(ns("validate_btn"), "1. Validate uploads", icon = icon("check-double"), class = "btn-primary btn-sm", width = "100%")
    ),
    column(
      8,
      tabsetPanel(
        id = ns("live_tabs"), type = "tabs",
        tabPanel("1. Upload & Validate", br(), uiOutput(ns("validate_ui"))),
        tabPanel("2. Sample Matching & Missing Data", br(), uiOutput(ns("matching_ui"))),
        tabPanel("3. Normalization, Filtering & Scaling", br(), uiOutput(ns("preprocess_ui"))),
        tabPanel("4. Batch Diagnostics", br(), uiOutput(ns("batch_ui"))),
        tabPanel("5-8. MOFA2 Integration & Results", br(), mod_multi_live_mofa_ui(ns("mofa")))
      )
    )
  )
}

mod_multi_live_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw <- reactiveValues(mats = list(), validations = list(), labels = list(), meta = NULL)
    proc <- reactiveValues(matched_mats = NULL, filtered_mats = NULL, scaled_mats = NULL, batch_corrected = NULL)
    live_state <- reactiveValues(mats = NULL, labels = NULL, meta = NULL, matched_ids = NULL)

    ## ---- 1. Upload & Validate ---------------------------------------------
    observeEvent(input$validate_btn, {
      mats <- list(); validations <- list(); labels <- list()
      for (i in seq_along(MULTI_LIVE_LAYER_IDS)) {
        lid <- MULTI_LIVE_LAYER_IDS[i]
        file_input <- input[[paste0(lid, "_file")]]
        if (is.null(file_input)) next
        orient <- input[[paste0(lid, "_orient")]] %||% "samples_rows"
        label <- input[[paste0(lid, "_label")]] %||% sprintf("Layer %d", i)
        res <- multi_live_read_matrix(file_input$datapath, orientation = orient)
        if (!res$ok) { showNotification(sprintf("%s: %s", label, res$error), type = "error"); next }
        v <- multi_live_validate_matrix(res$mat, layer_label = label)
        mats[[label]] <- res$mat
        validations[[label]] <- v
        labels[[lid]] <- label
      }
      if (length(mats) < 2) {
        showNotification("Upload at least two omics layers to proceed.", type = "warning")
      }
      raw$mats <- mats
      raw$validations <- validations
      raw$labels <- labels

      if (!is.null(input$meta_file)) {
        m <- tryCatch(as.data.frame(data.table::fread(input$meta_file$datapath, showProgress = FALSE)), error = function(e) NULL)
        if (!is.null(m) && ncol(m) >= 1) {
          rownames(m) <- as.character(m[[1]])
          raw$meta <- m
        }
      } else raw$meta <- NULL

      showNotification(sprintf("Validated %d omics layer(s).", length(mats)), type = "message")
    })

    output$validate_ui <- renderUI({
      if (length(raw$validations) == 0) return(multi_empty_state("Upload at least two omics matrices and click \"Validate uploads\"."))
      tbl <- multi_live_qc_summary_table(raw$validations)
      tagList(
        DT::dataTableOutput(ns("qc_table")),
        p(class = "submodule-desc", tags$em("Every value above is computed from your actual uploaded matrices - nothing here is a hardcoded example.")),
        lapply(raw$validations, function(v) {
          if (!isTRUE(v$ok)) return(NULL)
          issues <- c(
            if (v$n_duplicate_samples > 0) sprintf("%d duplicate sample ID(s)", v$n_duplicate_samples),
            if (v$n_duplicate_features > 0) sprintf("%d duplicate feature ID(s)", v$n_duplicate_features),
            if (v$n_zero_variance > 0) sprintf("%d zero-variance feature(s)", v$n_zero_variance),
            if (v$n_non_finite > 0) sprintf("%d non-finite (Inf/NaN) value(s)", v$n_non_finite)
          )
          if (length(issues) == 0) return(NULL)
          div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
              sprintf(" %s: %s.", v$layer, paste(issues, collapse = "; ")))
        })
      )
    })
    output$qc_table <- DT::renderDataTable({
      tbl <- req(multi_live_qc_summary_table(raw$validations))
      DT::datatable(tbl, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    ## ---- 2. Sample Matching & Missing Data ---------------------------------
    overlap <- reactive({
      req(length(raw$mats) >= 2)
      multi_live_sample_overlap(raw$mats)
    })

    output$matching_ui <- renderUI({
      ov <- tryCatch(overlap(), error = function(e) NULL)
      if (is.null(ov) || !isTRUE(ov$ok)) return(multi_empty_state("Validate at least two omics layers first."))
      tagList(
        div(style = "display:flex; gap:10px; flex-wrap:wrap;",
            lapply(names(ov$per_layer), function(nm) div(class = "card", style = "flex:1 1 140px; text-align:center; padding:10px;",
                                                            div(style = sprintf("font-size:1.3em; font-weight:600; color:%s;", ARTHOMIX_COLORS$blue), ov$per_layer[[nm]]),
                                                            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", nm))),
            div(class = "card", style = "flex:1 1 140px; text-align:center; padding:10px;",
                div(style = sprintf("font-size:1.3em; font-weight:600; color:%s;", ARTHOMIX_COLORS$aqua), ov$n_shared),
                div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", "Shared (matched)"))
        ),
        if (ov$n_shared < 3) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
                                   " Fewer than 3 samples are shared across all layers - integration below will not be possible until sample IDs are aligned.")
        else div(class = "empty-note", icon("circle-check"), sprintf(" Matched-sample integration below will use exactly these %d samples.", ov$n_shared)),
        hr(),
        h5("Missing-data QC"),
        selectInput(ns("miss_layer"), "Layer", choices = names(raw$mats)),
        multi_plot_or_empty(function() multi_live_missingness_by_omics_plot(raw$validations), ns("miss_by_omics"), height = "260px"),
        fluidRow(
          column(6, sliderInput(ns("max_sample_missing"), "Max sample missingness (%)", min = 0, max = 100, value = 50)),
          column(6, sliderInput(ns("max_feature_missing"), "Max feature missingness (%)", min = 0, max = 100, value = 50))
        ),
        multi_plot_or_empty(sample_miss_plot_fn, ns("sample_miss_plot"), height = "300px"),
        multi_plot_or_empty(feature_miss_plot_fn, ns("feature_miss_plot"), height = "260px"),
        selectInput(ns("impute_method"), "Missing-value handling (applied per layer before normalization)",
                    choices = c("Leave as-is" = "none", "Mean imputation" = "mean", "Median imputation" = "median",
                                "Remove samples/features exceeding thresholds" = "remove_rows"), selected = "none"),
        p(class = "submodule-desc", tags$strong("Missing-value handling can influence downstream factor structure and biological interpretation."), " Nothing is imputed until you pick a method above and it is applied in step 3.")
      )
    })
    output$miss_by_omics <- renderPlot(multi_live_missingness_by_omics_plot(raw$validations))
    sample_miss_plot_fn <- reactive(multi_live_sample_missingness_plot(multi_live_missingness(raw$mats[[req(input$miss_layer)]]), threshold = input$max_sample_missing))
    output$sample_miss_plot <- renderPlot(sample_miss_plot_fn())
    feature_miss_plot_fn <- reactive(multi_live_feature_missingness_plot(multi_live_missingness(raw$mats[[req(input$miss_layer)]])))
    output$feature_miss_plot <- renderPlot(feature_miss_plot_fn())

    ## ---- 3. Normalization, Filtering & Scaling -----------------------------
    output$preprocess_ui <- renderUI({
      if (length(raw$mats) < 2) return(multi_empty_state("Validate at least two omics layers first."))
      tagList(
        p(class = "submodule-desc", "Normalization choices are restricted to what's appropriate for each layer's declared omics type - never a one-size-fits-all list."),
        uiOutput(ns("norm_controls")),
        actionButton(ns("preprocess_btn"), "Apply normalization + filtering + scaling", icon = icon("gears"), class = "btn-primary btn-sm"),
        p(class = "submodule-desc", tags$em("Normalization should be selected according to assay type, measurement scale, experimental design, and downstream analysis - it is never applied silently.")),
        uiOutput(ns("preprocess_result_ui"))
      )
    })

    output$norm_controls <- renderUI({
      req(length(raw$mats) > 0)
      tagList(lapply(seq_along(MULTI_LIVE_LAYER_IDS), function(i) {
        lid <- MULTI_LIVE_LAYER_IDS[i]
        label <- input[[paste0(lid, "_label")]]
        if (is.null(label) || !label %in% names(raw$mats)) return(NULL)
        otype <- input[[paste0(lid, "_type")]] %||% "other"
        box(width = NULL, title = label, status = "primary", solidHeader = FALSE, collapsible = TRUE,
            selectInput(ns(paste0(lid, "_norm")), "Normalization", choices = MULTI_LIVE_NORM_CHOICES[[otype]]),
            numericInput(ns(paste0(lid, "_topvar")), "Feature filtering: keep top-N most variable features (blank = keep all)", value = NA, min = 10))
      }))
    })

    observeEvent(input$preprocess_btn, {
      ov <- tryCatch(overlap(), error = function(e) NULL)
      validate(need(!is.null(ov) && isTRUE(ov$ok) && ov$n_shared >= 3, "Need at least 3 matched samples across layers before preprocessing."))
      matched <- ov$shared_ids

      out_mats <- list()
      for (i in seq_along(MULTI_LIVE_LAYER_IDS)) {
        lid <- MULTI_LIVE_LAYER_IDS[i]
        label <- input[[paste0(lid, "_label")]]
        if (is.null(label) || !label %in% names(raw$mats)) next
        m <- raw$mats[[label]][matched, , drop = FALSE]

        imp <- multi_live_handle_missing(m, method = input$impute_method %||% "none",
                                          max_sample_missing_pct = input$max_sample_missing %||% 100,
                                          max_feature_missing_pct = input$max_feature_missing %||% 100)
        if (!imp$ok) { showNotification(sprintf("%s: %s", label, imp$error), type = "error"); next }
        m <- imp$mat

        norm <- multi_live_normalize(m, otype <- input[[paste0(lid, "_type")]] %||% "other", input[[paste0(lid, "_norm")]] %||% "none")
        m <- norm$mat

        top_n <- input[[paste0(lid, "_topvar")]]
        if (!is.null(top_n) && !is.na(top_n) && top_n > 0) {
          filt <- multi_live_filter_features(m, criterion = "variance", keep_top_n = top_n)
          m <- filt$mat
        }
        out_mats[[label]] <- m
      }
      req(length(out_mats) >= 2)
      proc$filtered_mats <- out_mats
      proc$scaled_mats <- lapply(out_mats, multi_live_scale)
      showNotification("Preprocessing applied.", type = "message")
    })

    output$preprocess_result_ui <- renderUI({
      if (is.null(proc$filtered_mats)) return(NULL)
      tagList(
        h5("Feature retention"),
        fluidRow(lapply(names(proc$filtered_mats), function(nm) column(6,
          multi_plot_or_empty(function() multi_live_retention_plot(ncol(raw$mats[[nm]]), ncol(proc$filtered_mats[[nm]])), ns(paste0("retention_", make.names(nm))), height = "220px")
        ))),
        h5("Distribution before / after normalization"),
        selectInput(ns("dist_layer"), "Layer", choices = names(proc$filtered_mats)),
        radioButtons(ns("dist_kind"), NULL, choices = c("Boxplot" = "box", "Density" = "density"), inline = TRUE),
        fluidRow(
          column(6, p(tags$strong("Before")), multi_plot_or_empty(before_dist_fn, ns("dist_before"), height = "280px")),
          column(6, p(tags$strong("After")), multi_plot_or_empty(after_dist_fn, ns("dist_after"), height = "280px"))
        ),
        h5("Cross-omics scale comparison (after z-score standardization)"),
        multi_plot_or_empty(scale_compare_fn, ns("scale_compare"), height = "300px")
      )
    })
    before_dist_fn <- reactive(multi_live_distribution_plot(raw$mats[[req(input$dist_layer)]], kind = input$dist_kind %||% "box"))
    output$dist_before <- renderPlot(before_dist_fn())
    after_dist_fn <- reactive(multi_live_distribution_plot(proc$filtered_mats[[req(input$dist_layer)]], kind = input$dist_kind %||% "box"))
    output$dist_after <- renderPlot(after_dist_fn())
    scale_compare_fn <- reactive({
      req(proc$scaled_mats)
      multi_live_scale_comparison_plot(proc$scaled_mats, names(proc$scaled_mats))
    })
    output$scale_compare <- renderPlot(scale_compare_fn())

    ## Dynamic per-layer retention plots (the number of layers is not fixed) -
    ## the UI above already references output ids "retention_<layer>"; this
    ## observe() is what actually assigns each one once the layer set is known.
    observe({
      req(proc$filtered_mats)
      for (nm in names(proc$filtered_mats)) {
        local({
          nm_local <- nm
          output[[paste0("retention_", make.names(nm_local))]] <- renderPlot(multi_live_retention_plot(ncol(raw$mats[[nm_local]]), ncol(proc$filtered_mats[[nm_local]])))
        })
      }
    })

    ## ---- 4. Batch Diagnostics ----------------------------------------------
    output$batch_ui <- renderUI({
      if (is.null(proc$scaled_mats)) return(multi_empty_state("Apply preprocessing (step 3) first."))
      meta_cols <- if (!is.null(raw$meta)) colnames(raw$meta) else character(0)
      tagList(
        selectInput(ns("batch_layer"), "Layer", choices = names(proc$scaled_mats)),
        if (length(meta_cols) == 0) div(class = "empty-note", icon("circle-info"), "Batch information is unavailable (no metadata uploaded) - batch-effect assessment cannot be performed. PCA below still shows structure without a color-by variable.")
        else tagList(
          selectInput(ns("color_by"), "Color PCA by", choices = c("(none)" = "", meta_cols)),
          selectInput(ns("batch_col"), "Batch column", choices = c("(none)" = "", meta_cols)),
          selectInput(ns("phenotype_col"), "Phenotype/group column", choices = c("(none)" = "", meta_cols))
        ),
        h5("PCA before correction"),
        multi_plot_or_empty(pca_before_fn, ns("pca_before"), height = "360px"),
        uiOutput(ns("confound_ui")),
        conditionalPanel(condition = sprintf("input['%s'] != '' && input['%s'] != ''", ns("batch_col"), ns("phenotype_col")), tagList(
          selectInput(ns("correct_method"), "Correction method", choices = c("ComBat (empirical Bayes)" = "combat", "limma::removeBatchEffect" = "limma")),
          actionButton(ns("correct_btn"), "Apply batch correction", icon = icon("play"), class = "btn-primary btn-sm"),
          h5("PCA after correction"),
          multi_plot_or_empty(pca_after_fn, ns("pca_after"), height = "360px"),
          uiOutput(ns("variance_diagnostic_ui"))
        ))
      )
    })
    pca_before_fn <- reactive(multi_live_pca_plot(multi_live_pca(proc$scaled_mats[[req(input$batch_layer)]]), raw$meta, if (nzchar(input$color_by %||% "")) input$color_by else NULL))
    output$pca_before <- renderPlot(pca_before_fn())

    output$confound_ui <- renderUI({
      req(input$batch_col, input$phenotype_col, nzchar(input$batch_col), nzchar(input$phenotype_col), raw$meta)
      cc <- multi_live_confounding_check(raw$meta, input$batch_col, input$phenotype_col)
      if (is.null(cc)) return(NULL)
      if (isTRUE(cc$confounded)) {
        div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
            " The chosen batch column appears confounded with the phenotype column (at least one batch maps to a single phenotype level) - correcting for batch here risks removing real biological signal along with technical variation. Proceed with caution.")
      } else {
        div(class = "empty-note", icon("circle-check"), sprintf(" No strong batch/phenotype confounding detected (chi-square p = %.3f).", cc$p_value %||% NA))
      }
    })

    observeEvent(input$correct_btn, {
      req(proc$scaled_mats, input$batch_layer, input$batch_col, raw$meta)
      m <- proc$scaled_mats[[input$batch_layer]]
      common <- intersect(rownames(m), rownames(raw$meta))
      validate(need(length(common) >= 3, "Not enough samples with both data and batch metadata."))
      batch <- raw$meta[common, input$batch_col]
      res <- multi_live_batch_correct(m[common, , drop = FALSE], batch, method = input$correct_method %||% "combat")
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      proc$batch_corrected <- res$mat
      showNotification("Batch correction applied.", type = "message")
    })
    pca_after_fn <- reactive({
      req(proc$batch_corrected)
      multi_live_pca_plot(multi_live_pca(proc$batch_corrected), raw$meta, if (nzchar(input$color_by %||% "")) input$color_by else NULL)
    })
    output$pca_after <- renderPlot(pca_after_fn())

    output$variance_diagnostic_ui <- renderUI({
      req(proc$batch_corrected, input$batch_col, input$phenotype_col, nzchar(input$batch_col), nzchar(input$phenotype_col))
      pca_b <- multi_live_pca(proc$scaled_mats[[input$batch_layer]])
      pca_a <- multi_live_pca(proc$batch_corrected)
      if (!isTRUE(pca_b$ok) || !isTRUE(pca_a$ok)) return(NULL)
      r2_batch_before <- multi_live_variance_by_group(pca_b$scores, raw$meta, input$batch_col)
      r2_pheno_before <- multi_live_variance_by_group(pca_b$scores, raw$meta, input$phenotype_col)
      r2_batch_after <- multi_live_variance_by_group(pca_a$scores, raw$meta, input$batch_col)
      r2_pheno_after <- multi_live_variance_by_group(pca_a$scores, raw$meta, input$phenotype_col)
      tagList(
        h5("Quantitative diagnostic (R² of PC1/PC2 vs. group, not a \"looks better\" claim)"),
        DT::dataTableOutput(ns("variance_diag_table")),
        p(class = "submodule-desc", "A successful correction typically reduces R² against batch while keeping R² against phenotype roughly stable. Strong correction can remove true biological signal - if phenotype R² also collapses, reconsider the correction method or covariates.")
      )
    })
    output$variance_diag_table <- DT::renderDataTable({
      pca_b <- multi_live_pca(proc$scaled_mats[[req(input$batch_layer)]])
      pca_a <- multi_live_pca(req(proc$batch_corrected))
      df <- data.frame(
        PC = c("PC1", "PC2"),
        `R2 vs batch (before)` = multi_live_variance_by_group(pca_b$scores, raw$meta, input$batch_col),
        `R2 vs batch (after)` = multi_live_variance_by_group(pca_a$scores, raw$meta, input$batch_col),
        `R2 vs phenotype (before)` = multi_live_variance_by_group(pca_b$scores, raw$meta, input$phenotype_col),
        `R2 vs phenotype (after)` = multi_live_variance_by_group(pca_a$scores, raw$meta, input$phenotype_col),
        check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    ## ---- Publish state for the nested MOFA2 sub-module ---------------------
    observe({
      final_mats <- if (!is.null(proc$batch_corrected)) {
        c(proc$scaled_mats[setdiff(names(proc$scaled_mats), input$batch_layer)], setNames(list(proc$batch_corrected), input$batch_layer))
      } else proc$scaled_mats
      live_state$mats <- final_mats
      live_state$meta <- raw$meta
    })

    mod_multi_live_mofa_server("mofa", live_state, multi_results)

    observe({
      if (is.null(multi_results)) return()
      multi_results$live_qc <- list(n_layers = length(raw$mats), validations = raw$validations)
    })
  })
}
