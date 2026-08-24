## R/multiomics/mod_multi_live_mofa.R
## Nested sub-module mounted inside mod_multi_live.R (sections 5-8: MOFA2
## Integration, Factor Results, Cross-Omics Correlation, Export). Reads the
## `live_state` reactiveValues mod_multi_live_server publishes (final
## preprocessed/scaled, matched-sample matrices per layer + uploaded
## metadata) - never touches the precomputed-cohort tabs' own state.
##
## MOFA2 training runs via shiny::ExtendedTask + future::future_promise(),
## gated by ARTHOMIX_ASYNC_AVAILABLE (global.R) exactly like
## mod_methyl_dataset.R's own ~2.1GB preloaded-matrix read - the one other
## genuinely slow operation in this app - so training doesn't freeze the
## rest of the app for every open session. Falls back to a blocking call
## with a clear notice if `future`/`promises` aren't installed.

mod_multi_live_mofa_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("mofa_availability_note")),
    fluidRow(
      column(
        4,
        box(
          width = NULL, title = "5. MOFA2 Integration", status = "primary", solidHeader = FALSE,
          uiOutput(ns("readiness_ui")),
          numericInput(ns("num_factors"), "Number of factors", value = 10, min = 2, max = 25),
          numericInput(ns("seed"), "Random seed", value = 1, min = 1),
          selectInput(ns("convergence"), "Convergence mode", choices = c("Fast" = "fast", "Medium" = "medium", "Slow (most precise)" = "slow")),
          p(class = "submodule-desc", "MOFA2 (Argelaguet et al. 2018/2020) is an unsupervised latent-factor model: X_m ~ Z W_m^T + E_m per view m. Factors capture major sources of variation across the integrated views - they are not automatically disease factors, biomarkers, or causal."),
          actionButton(ns("train_btn"), "Train MOFA2", icon = icon("play"), class = "btn-primary btn-sm", width = "100%"),
          uiOutput(ns("train_status_ui"))
        )
      ),
      column(
        8,
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("train_btn")),
          tabsetPanel(
            id = ns("mofa_tabs"), type = "tabs",
            tabPanel("6. Variance Explained", br(), uiOutput(ns("variance_ui"))),
            tabPanel("6. Factor Scores", br(), uiOutput(ns("scores_ui"))),
            tabPanel("6. Factor Heatmap", br(), uiOutput(ns("heatmap_ui"))),
            tabPanel("6. Feature Loadings", br(), uiOutput(ns("loadings_ui"))),
            tabPanel("7. Cross-Omics Correlation", br(), uiOutput(ns("correlation_ui"))),
            tabPanel("8. Export", br(), uiOutput(ns("export_ui")))
          )
        )
      )
    )
  )
}

mod_multi_live_mofa_server <- function(id, live_state, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    model_state <- reactiveValues(model = NULL, variance_df = NULL, factors_df = NULL, loadings_df = NULL, error = NULL, trained_at = NULL, params = NULL)

    output$mofa_availability_note <- renderUI({
      if (MULTI_MOFA_AVAILABLE) return(NULL)
      div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
          " MOFA2 is not installed in this deployment - the Integration step below is unavailable. This is stated plainly, not worked around with a substitute result.")
    })

    output$readiness_ui <- renderUI({
      mats <- live_state$mats
      if (is.null(mats) || length(mats) < 2) return(div(class = "empty-note", icon("circle-info"), "Complete steps 1-3 (Upload, Matching, Preprocessing) first - at least two preprocessed layers are needed."))
      gr <- multi_live_mofa_guardrails(mats)
      tagList(
        div(class = "empty-note", icon("circle-check"), sprintf(" %d layers ready, %s matched samples, %s total features.", length(mats), format(gr$n_samples, big.mark = ","), format(gr$total_features, big.mark = ","))),
        lapply(gr$warnings, function(w) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"), paste("", w)))
      )
    })

    ## ---- Async training (falls back to blocking if future/promises absent) --
    run_training <- function() {
      multi_live_run_mofa(live_state$mats, num_factors = input$num_factors %||% 10, seed = input$seed %||% 1, convergence_mode = input$convergence %||% "fast")
    }

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      mofa_task <- ExtendedTask$new(function(mats, num_factors, seed, convergence) {
        promises::future_promise(multi_live_run_mofa(mats, num_factors = num_factors, seed = seed, convergence_mode = convergence), seed = TRUE)
      })
      observeEvent(input$train_btn, {
        req(MULTI_MOFA_AVAILABLE, live_state$mats)
        model_state$error <- NULL
        mofa_task$invoke(live_state$mats, input$num_factors %||% 10, input$seed %||% 1, input$convergence %||% "fast")
      })
      observe({
        res <- tryCatch(mofa_task$result(), error = function(e) e)
        if (inherits(res, "shiny.silent.error")) return()
        if (inherits(res, "error")) { model_state$error <- conditionMessage(res); return() }
        if (!isTRUE(res$ok)) { model_state$error <- res$error; return() }
        model_state$model <- res$model
        model_state$variance_df <- multi_live_mofa_variance_df(res$model)
        model_state$factors_df <- multi_live_mofa_factors_df(res$model)
        model_state$loadings_df <- multi_live_mofa_loadings_df(res$model)
        model_state$trained_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        model_state$params <- list(num_factors = res$num_factors, seed = res$seed)
      })
      output$train_status_ui <- renderUI({
        st <- mofa_task$status()
        if (identical(st, "running")) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Training MOFA2 in the background - the rest of the app stays usable while this runs."))
        if (!is.null(model_state$error)) return(div(class = "empty-note", style = "border-color: var(--color-danger, #e34948);", icon("triangle-exclamation"), paste(" ", model_state$error)))
        if (!is.null(model_state$model)) return(div(class = "empty-note", icon("circle-check"), sprintf(" Trained at %s - %d factors.", model_state$trained_at, model_state$params$num_factors %||% NA)))
        NULL
      })
    } else {
      observeEvent(input$train_btn, {
        req(MULTI_MOFA_AVAILABLE, live_state$mats)
        model_state$error <- NULL
        showNotification("Training MOFA2 synchronously (future/promises not installed) - the app will be briefly unresponsive.", type = "message", duration = 5)
        res <- run_training()
        if (!isTRUE(res$ok)) { model_state$error <- res$error; return() }
        model_state$model <- res$model
        model_state$variance_df <- multi_live_mofa_variance_df(res$model)
        model_state$factors_df <- multi_live_mofa_factors_df(res$model)
        model_state$loadings_df <- multi_live_mofa_loadings_df(res$model)
        model_state$trained_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        model_state$params <- list(num_factors = res$num_factors, seed = res$seed)
      })
      output$train_status_ui <- renderUI({
        if (!is.null(model_state$error)) return(div(class = "empty-note", style = "border-color: var(--color-danger, #e34948);", icon("triangle-exclamation"), paste(" ", model_state$error)))
        if (!is.null(model_state$model)) return(div(class = "empty-note", icon("circle-check"), sprintf(" Trained at %s - %d factors.", model_state$trained_at, model_state$params$num_factors %||% NA)))
        NULL
      })
    }

    ## ---- 6. Factor Results --------------------------------------------------
    output$variance_ui <- renderUI({
      if (is.null(model_state$variance_df)) return(multi_empty_state("Train MOFA2 first."))
      tagList(
        multi_plot_or_empty(variance_plot_fn, ns("variance_plot"), height = "380px"),
        div(class = "table-toolbar", downloadButton(ns("dl_variance_png"), "Download plot (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_variance_csv"), "Download data (CSV)", class = "btn-sm")),
        p(class = "submodule-desc", tags$em("How to read this: ", "each bar segment is the fraction of one view's total variance captured by that factor. A factor with high variance explained captures a major source of variation in the data - this alone does not mean it is biologically or clinically meaningful."))
      )
    })
    variance_plot_fn <- reactive(multi_live_mofa_variance_plot(model_state$variance_df))
    output$variance_plot <- renderPlot(variance_plot_fn())
    output$dl_variance_png <- multi_png_download(variance_plot_fn, function() "multiomics_factor_variance_explained.png")
    output$dl_variance_csv <- downloadHandler(function() "multiomics_factor_variance_explained.csv", function(file) utils::write.csv(model_state$variance_df, file, row.names = FALSE))

    output$scores_ui <- renderUI({
      if (is.null(model_state$factors_df)) return(multi_empty_state("Train MOFA2 first."))
      facs <- colnames(model_state$factors_df)
      meta_cols <- if (!is.null(live_state$meta)) colnames(live_state$meta) else character(0)
      tagList(
        fluidRow(
          column(4, selectInput(ns("x_factor"), "X factor", choices = facs, selected = facs[1])),
          column(4, selectInput(ns("y_factor"), "Y factor", choices = facs, selected = facs[min(2, length(facs))])),
          column(4, selectInput(ns("score_color_by"), "Color by", choices = c("(none)" = "", meta_cols)))
        ),
        multi_plot_or_empty(score_plot_fn, ns("score_plot"), height = "420px"),
        div(class = "table-toolbar", downloadButton(ns("dl_score_png"), "Download plot (PNG)", class = "btn-sm")),
        p(class = "submodule-desc", tags$em("How to read this: ", "samples close together have similar scores on the selected factors. Any apparent grouping should be checked against real metadata association (color-by above) before being called meaningful."))
      )
    })
    score_plot_fn <- reactive({
      req(input$x_factor, input$y_factor)
      multi_live_factor_score_plot(model_state$factors_df, input$x_factor, input$y_factor, live_state$meta, if (nzchar(input$score_color_by %||% "")) input$score_color_by else NULL)
    })
    output$score_plot <- renderPlot(score_plot_fn())
    output$dl_score_png <- multi_png_download(score_plot_fn, function() sprintf("multiomics_%s_%s.png", input$x_factor %||% "F1", input$y_factor %||% "F2"))

    output$heatmap_ui <- renderUI({
      if (is.null(model_state$factors_df)) return(multi_empty_state("Train MOFA2 first."))
      tagList(
        multi_plot_or_empty(heatmap_plot_fn, ns("heatmap_plot"), height = "460px"),
        div(class = "table-toolbar", downloadButton(ns("dl_heatmap_png"), "Download plot (PNG)", class = "btn-sm"))
      )
    })
    heatmap_plot_fn <- reactive(multi_live_factor_heatmap(model_state$factors_df))
    output$heatmap_plot <- renderPlot(heatmap_plot_fn())
    output$dl_heatmap_png <- multi_png_download(heatmap_plot_fn, function() "multiomics_factor_heatmap.png")

    output$loadings_ui <- renderUI({
      if (is.null(model_state$loadings_df)) return(multi_empty_state("Train MOFA2 first."))
      facs <- unique(model_state$loadings_df$factor)
      tagList(
        fluidRow(
          column(4, selectInput(ns("loadings_factor"), "Factor", choices = facs)),
          column(4, selectInput(ns("loadings_sign"), "Sign", choices = c("Both" = "both", "Positive only" = "positive", "Negative only" = "negative"))),
          column(4, selectInput(ns("loadings_top_n"), "Show top", choices = c("10" = 10, "20" = 20, "50" = 50), selected = 20))
        ),
        multi_plot_or_empty(loadings_plot_fn, ns("loadings_plot"), height = "460px"),
        div(class = "table-toolbar", downloadButton(ns("dl_loadings_png"), "Download plot (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_loadings_csv"), "Download data (CSV)", class = "btn-sm")),
        p(class = "submodule-desc", tags$em("How to read this: ", "larger absolute loading means a feature contributes more strongly to that factor's direction. A high loading is not by itself evidence of clinical importance."))
      )
    })
    loadings_plot_fn <- reactive({
      req(input$loadings_factor)
      df <- model_state$loadings_df[model_state$loadings_df$factor == input$loadings_factor, , drop = FALSE]
      multi_live_loadings_plot(df, sign = input$loadings_sign %||% "both", top_n = as.integer(input$loadings_top_n %||% 20))
    })
    output$loadings_plot <- renderPlot(loadings_plot_fn())
    output$dl_loadings_png <- multi_png_download(loadings_plot_fn, function() sprintf("multiomics_factor_loadings_%s.png", input$loadings_factor %||% "factor"))
    output$dl_loadings_csv <- downloadHandler(function() "multiomics_factor_loadings.csv", function(file) utils::write.csv(model_state$loadings_df, file, row.names = FALSE))

    ## ---- 7. Cross-Omics Correlation -----------------------------------------
    output$correlation_ui <- renderUI({
      mats <- live_state$mats
      if (is.null(mats) || length(mats) < 2) return(multi_empty_state("Complete preprocessing (steps 1-3) first."))
      tagList(
        h5("Single feature pair"),
        fluidRow(
          column(3, selectInput(ns("corr_omicsA"), "Omics A", choices = names(mats))),
          column(3, uiOutput(ns("featA_ui"))),
          column(3, selectInput(ns("corr_omicsB"), "Omics B", choices = names(mats), selected = names(mats)[min(2, length(mats))])),
          column(3, uiOutput(ns("featB_ui")))
        ),
        radioButtons(ns("corr_method"), "Method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"), inline = TRUE),
        multi_plot_or_empty(scatter_fn, ns("corr_scatter"), height = "360px"),
        hr(),
        h5("Correlation heatmap (top-variance features per side)"),
        fluidRow(
          column(4, sliderInput(ns("heatmap_top_n"), "Top-N features per side", min = 5, max = 60, value = 20)),
          column(4, sliderInput(ns("heatmap_fdr"), "FDR filter (grey out above)", min = 0, max = 1, value = 1, step = 0.01))
        ),
        multi_plot_or_empty(heatmap_corr_fn, ns("corr_heatmap"), height = "480px"),
        div(class = "table-toolbar", downloadButton(ns("dl_corr_png"), "Download heatmap (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_corr_csv"), "Download pairwise data (CSV)", class = "btn-sm"))
      )
    })
    output$featA_ui <- renderUI({ req(input$corr_omicsA); selectInput(ns("corr_featA"), "Feature A", choices = colnames(live_state$mats[[input$corr_omicsA]])) })
    output$featB_ui <- renderUI({ req(input$corr_omicsB); selectInput(ns("corr_featB"), "Feature B", choices = colnames(live_state$mats[[input$corr_omicsB]])) })
    corr_result <- reactive({
      req(input$corr_omicsA, input$corr_omicsB, input$corr_featA, input$corr_featB)
      mA <- live_state$mats[[input$corr_omicsA]]; mB <- live_state$mats[[input$corr_omicsB]]
      common <- intersect(rownames(mA), rownames(mB))
      req(length(common) >= 3)
      multi_live_correlation(mA[common, input$corr_featA], mB[common, input$corr_featB], method = input$corr_method %||% "pearson")
    })
    scatter_fn <- reactive({
      r <- tryCatch(corr_result(), error = function(e) NULL)
      if (is.null(r) || !isTRUE(r$ok)) return(NULL)
      mA <- live_state$mats[[input$corr_omicsA]]; mB <- live_state$mats[[input$corr_omicsB]]
      common <- intersect(rownames(mA), rownames(mB))
      multi_live_correlation_scatter_plot(mA[common, input$corr_featA], mB[common, input$corr_featB],
                                           sprintf("%s: %s", input$corr_omicsA, input$corr_featA), sprintf("%s: %s", input$corr_omicsB, input$corr_featB),
                                           r$r, r$p, r$n)
    })
    output$corr_scatter <- renderPlot(scatter_fn())

    heatmap_corr_data <- reactive({
      req(input$corr_omicsA, input$corr_omicsB)
      multi_live_correlation_heatmap_data(live_state$mats[[input$corr_omicsA]], live_state$mats[[input$corr_omicsB]], top_n = input$heatmap_top_n %||% 20, method = input$corr_method %||% "pearson")
    })
    heatmap_corr_fn <- reactive({
      r <- tryCatch(heatmap_corr_data(), error = function(e) NULL)
      if (is.null(r) || !isTRUE(r$ok)) return(NULL)
      multi_live_correlation_heatmap_plot(r$df, fdr_threshold = input$heatmap_fdr %||% 1)
    })
    output$corr_heatmap <- renderPlot(heatmap_corr_fn())
    output$dl_corr_png <- multi_png_download(heatmap_corr_fn, function() "multiomics_cross_omics_correlation.png")
    output$dl_corr_csv <- downloadHandler(function() "multiomics_cross_omics_correlation.csv", function(file) {
      r <- heatmap_corr_data()
      utils::write.csv(if (isTRUE(r$ok)) r$df else data.frame(), file, row.names = FALSE)
    })

    ## ---- 8. Export -----------------------------------------------------------
    output$export_ui <- renderUI({
      if (is.null(model_state$model)) return(multi_empty_state("Train MOFA2 first."))
      tagList(
        h5("Factor summary"),
        DT::dataTableOutput(ns("factor_summary_table")),
        div(class = "table-toolbar", downloadButton(ns("dl_factor_summary_csv"), "Download factor summary (CSV)", class = "btn-sm")),
        hr(),
        h5("Reproducibility"),
        tags$ul(
          tags$li(sprintf("Trained at: %s", model_state$trained_at %||% "-")),
          tags$li(sprintf("Number of factors: %s", model_state$params$num_factors %||% "-")),
          tags$li(sprintf("Random seed: %s", model_state$params$seed %||% "-")),
          tags$li(sprintf("MOFA2 version: %s", tryCatch(as.character(utils::packageVersion("MOFA2")), error = function(e) "unknown")))
        ),
        div(class = "table-toolbar", downloadButton(ns("dl_bundle"), "Download full results bundle (ZIP)", class = "btn-sm btn-primary"))
      )
    })
    factor_summary_df <- reactive({
      req(model_state$variance_df, model_state$loadings_df)
      facs <- unique(model_state$variance_df$factor)
      do.call(rbind, lapply(facs, function(f) {
        vd <- model_state$variance_df[model_state$variance_df$factor == f, , drop = FALSE]
        main_view <- vd$view[which.max(vd$variance_explained)]
        ld <- model_state$loadings_df[model_state$loadings_df$factor == f, , drop = FALSE]
        ld <- ld[order(-abs(ld$value)), , drop = FALSE]
        data.frame(
          Factor = f,
          `Variance explained (max view)` = sprintf("%.1f%%", 100 * max(vd$variance_explained, na.rm = TRUE)),
          `Main omics contributor` = main_view,
          `Top positive feature` = ld$feature[which.max(ld$value)],
          `Top negative feature` = ld$feature[which.min(ld$value)],
          check.names = FALSE
        )
      }))
    })
    output$factor_summary_table <- DT::renderDataTable({
      DT::datatable(req(factor_summary_df()), rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_factor_summary_csv <- downloadHandler(function() "multiomics_factor_summary.csv", function(file) utils::write.csv(factor_summary_df(), file, row.names = FALSE))

    output$dl_bundle <- downloadHandler(
      filename = function() paste0("multiomics_live_mofa_bundle_", Sys.Date(), ".zip"),
      content = function(file) {
        tmp <- tempfile(); dir.create(tmp)
        if (!is.null(model_state$variance_df)) utils::write.csv(model_state$variance_df, file.path(tmp, "variance_explained.csv"), row.names = FALSE)
        if (!is.null(model_state$factors_df)) utils::write.csv(model_state$factors_df, file.path(tmp, "factor_scores.csv"), row.names = FALSE)
        if (!is.null(model_state$loadings_df)) utils::write.csv(model_state$loadings_df, file.path(tmp, "feature_loadings.csv"), row.names = FALSE)
        fs <- tryCatch(factor_summary_df(), error = function(e) NULL)
        if (!is.null(fs)) utils::write.csv(fs, file.path(tmp, "factor_summary.csv"), row.names = FALSE)
        writeLines(c(
          "# MOFA2 Live Analysis - reproducibility",
          sprintf("Trained at: %s", model_state$trained_at %||% "-"),
          sprintf("Number of factors: %s", model_state$params$num_factors %||% "-"),
          sprintf("Random seed: %s", model_state$params$seed %||% "-"),
          sprintf("MOFA2 version: %s", tryCatch(as.character(utils::packageVersion("MOFA2")), error = function(e) "unknown")),
          "Trained on your own uploaded, matched-sample, preprocessed data - not the precomputed cohort."
        ), file.path(tmp, "reproducibility.md"))
        old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
        setwd(tmp)
        utils::zip(file, files = list.files(tmp))
      }
    )

    observe({
      if (is.null(multi_results) || is.null(model_state$model)) return()
      multi_results$live_mofa <- list(trained_at = model_state$trained_at, params = model_state$params, n_factors = model_state$params$num_factors)
    })
  })
}
