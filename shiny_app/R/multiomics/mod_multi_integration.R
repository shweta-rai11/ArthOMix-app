## R/multiomics/mod_multi_integration.R
## Submodule: Multi-omics Integration (DIABLO & SNF) - for one of the six
## analysis cells, shows the pipeline's own already-computed DIABLO
## (mixOmics::block.splsda, supervised) and SNF (SNFtool::SNF, unsupervised
## patient-similarity fusion) performance, per-patient component scores, and
## panel loadings, plus the master six-part-summary comparison against
## single-omics baselines (spec: "show whether integration actually improves
## performance - do not assume it does"). No model is re-fit, no
## cross-validation is re-run here.

mod_multi_integration_config <- list(
  id = "integration", title = "Multi-omics Integration (DIABLO & SNF)", icon = "layer-group", group = "Data",
  description = "DIABLO (supervised, sparse multi-block discriminant analysis) and SNF (unsupervised patient-similarity fusion) performance, scores, and panels for one analysis cell - and whether integration actually beats single-omics."
)

mod_multi_integration_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Analysis cell", status = "primary", solidHeader = FALSE,
        selectInput(ns("cell"), "Sex / drug / outcome", choices = MULTI_CELL_CHOICES, width = "100%"),
        p(class = "submodule-desc", "Both DIABLO and SNF were run leave-one-out cross-validated, with feature selection strictly confined to each fold's training patients."),
        actionButton(ns("load_btn"), "Load cell", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")
      ),
      box(
        width = NULL, title = "Method notes", status = "primary", solidHeader = FALSE,
        p(tags$strong("DIABLO"), " - supervised: finds a sparse combination of features per omics block that jointly discriminates the response outcome. Appropriate when matched samples and a known phenotype are both available (this cohort's case)."),
        p(tags$strong("SNF"), " - unsupervised: fuses per-omics patient-similarity networks, independent of the outcome label; \"response-informed\" and \"unsupervised\" variants differ only in which features build the per-omics similarity matrix, not in the fusion step itself."),
        p(class = "submodule-desc", "Neither substitutes for the other - DIABLO answers \"which features discriminate response\", SNF answers \"do patients cluster by molecular similarity, and does that line up with response\" (see the Patient Stratification tab).")
      )
    ),
    column(
      8,
      uiOutput(ns("banner_ui")),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("load_btn")),
        tabsetPanel(
          id = ns("tabs"), type = "tabs",
          tabPanel("Performance", br(), uiOutput(ns("performance_ui"))),
          tabPanel("Variance explained", br(), uiOutput(ns("variance_ui"))),
          tabPanel("DIABLO scores", br(), uiOutput(ns("scores_ui"))),
          tabPanel("DIABLO panel", br(), uiOutput(ns("panel_ui"))),
          tabPanel("Integrated vs. single-omics", br(), uiOutput(ns("compare_ui")))
        )
      )
    )
  )
}

mod_multi_integration_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    dat <- reactiveValues(cell = NULL, perf = NULL, snf_perf = NULL, scores = NULL, panel = NULL, cell36 = NULL, baseline36 = NULL, variance_df = NULL)

    observeEvent(input$load_btn, {
      cell <- multi_cell_by_key(input$cell)
      req(cell)
      is_response <- is.na(cell$drug)

      perf_label  <- if (is_response) "DIABLO performance — response (drug-pooled)" else "DIABLO performance — drug x sex (response)"
      score_label <- if (is_response) "DIABLO scores — response (drug-pooled)" else "DIABLO scores — drug x sex (response)"
      panel_label <- if (is_response) "DIABLO panel — response (drug-pooled)" else "DIABLO panel — drug x sex (response)"

      perf  <- multi_read_registry_table(perf_label)
      score <- multi_read_registry_table(score_label)
      panel <- multi_read_registry_table(panel_label)
      s36   <- multi_read_registry_table("Master six-part summary (integrated vs single-omics)")
      snf   <- if (cell$has_snf) multi_read_registry_table("SNF performance — drug x sex") else list(ok = FALSE, df = NULL)

      if (!perf$ok) { showNotification(perf$error, type = "error"); return() }

      dat$cell <- cell
      dat$perf <- multi_filter_cell(perf$df, sex = cell$sex, drug = cell$drug)
      dat$scores <- if (score$ok) multi_filter_cell(score$df, sex = cell$sex, drug = cell$drug) else NULL
      dat$panel <- if (panel$ok) multi_filter_cell(panel$df, sex = cell$sex, drug = cell$drug) else NULL
      dat$snf_perf <- if (snf$ok) multi_filter_cell(snf$df, sex = cell$sex, drug = cell$drug) else NULL
      dat$cell36 <- if (s36$ok) multi_table36_for_cell(s36$df, cell) else NULL
      dat$baseline36 <- if (s36$ok) multi_table36_single_omics_baseline(s36$df, cell) else NULL

      fit_res <- multi_diablo_fit(cell)
      dat$variance_df <- if (fit_res$ok) multi_diablo_variance_df(fit_res$fit) else NULL

      showNotification(sprintf("Loaded %s.", cell$label), type = "message")
    })

    output$banner_ui <- renderUI({
      if (is.null(dat$cell)) return(multi_empty_state("Pick an analysis cell and click \"Load cell\"."))
      p <- dat$perf
      if (is.null(p) || nrow(p) == 0) return(div(class = "empty-note", icon("triangle-exclamation"), "No DIABLO performance row for this cell."))
      ex <- p$excludes_chance %in% TRUE
      div(class = "empty-note", style = if (any(ex)) "border-color: var(--color-warning, #eda100);" else NULL,
          icon(if (any(ex)) "circle-check" else "triangle-exclamation"),
          sprintf(" %s - DIABLO AUROC %.2f [%.2f, %.2f], n=%s. %s",
                  dat$cell$label, p$auroc[1], p$ci_lo[1], p$ci_hi[1], p$n[1],
                  if (any(ex)) "This cell's CI excludes chance." else "This cell's CI includes chance (0.5) - treat as exploratory, not a validated model."))
    })

    output$performance_ui <- renderUI({
      req(dat$cell)
      rows <- list(dat$perf)
      if (!is.null(dat$snf_perf)) rows <- c(rows, list(dat$snf_perf))
      combined <- tryCatch(do.call(function(...) {
        dfs <- list(...)
        dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0, dfs)
        if (length(dfs) == 0) return(NULL)
        common <- Reduce(intersect, lapply(dfs, colnames))
        do.call(rbind, lapply(dfs, function(d) d[, common, drop = FALSE]))
      }, rows), error = function(e) NULL)
      if (is.null(combined)) return(multi_empty_state())
      tagList(
        multi_plot_or_empty(perf_plot_fn, ns("perf_plot"), height = "260px"),
        div(class = "table-toolbar", downloadButton(ns("dl_perf_png"), "Download plot (PNG)", class = "btn-sm")),
        p(class = "submodule-desc", tags$em("How to read this: ", "each point is one model's AUROC with its 95% confidence interval; aqua means the interval sits entirely above chance (0.5), grey means it does not. A statistically significant AUROC is not by itself evidence of clinical utility.")),
        DT::dataTableOutput(ns("perf_table"))
      )
    })
    perf_plot_fn <- reactive({
      req(dat$perf)
      combined <- rbind(
        cbind(dat$perf[, intersect(c("model", "auroc", "ci_lo", "ci_hi", "excludes_chance"), colnames(dat$perf)), drop = FALSE]),
        if (!is.null(dat$snf_perf)) {
          sp <- dat$snf_perf
          data.frame(model = paste0("SNF (", sp$variant, ")"), auroc = sp$auroc, ci_lo = sp$ci_lo, ci_hi = sp$ci_hi, excludes_chance = sp$excludes_chance)
        } else NULL
      )
      multi_performance_ci_plot(combined, label_col = "model", title = dat$cell$label)
    })
    output$perf_plot <- renderPlot(perf_plot_fn())
    output$dl_perf_png <- multi_png_download(perf_plot_fn, function() sprintf("multiomics_integration_performance_%s.png", dat$cell$key))
    output$perf_table <- DT::renderDataTable({
      req(dat$perf)
      DT::datatable(dat$perf, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$scores_ui <- renderUI({
      if (is.null(dat$scores) || nrow(dat$scores) == 0) return(multi_empty_state("No per-patient score table for this cell."))
      tagList(
        multi_plot_or_empty(scores_plot_fn, ns("scores_plot"), height = "360px"),
        div(class = "table-toolbar", downloadButton(ns("dl_scores_png"), "Download plot (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_scores_csv"), "Download data (CSV)", class = "btn-sm"))
      )
    })
    scores_plot_fn <- reactive(multi_diablo_score_plot(dat$scores))
    output$scores_plot <- renderPlot(scores_plot_fn())
    output$dl_scores_png <- multi_png_download(scores_plot_fn, function() sprintf("multiomics_diablo_scores_%s.png", dat$cell$key))
    output$dl_scores_csv <- downloadHandler(function() sprintf("multiomics_diablo_scores_%s.csv", dat$cell$key), function(file) utils::write.csv(dat$scores, file, row.names = FALSE))

    output$panel_ui <- renderUI({
      if (is.null(dat$panel) || nrow(dat$panel) == 0) return(multi_empty_state("No panel/loadings table for this cell."))
      tagList(
        selectInput(ns("panel_top_n"), "Show top", choices = c("10" = 10, "20" = 20, "50" = 50), selected = 20, width = "160px"),
        multi_plot_or_empty(panel_plot_fn, ns("panel_plot"), height = "460px"),
        div(class = "table-toolbar", downloadButton(ns("dl_panel_png"), "Download plot (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_panel_csv"), "Download data (CSV)", class = "btn-sm")),
        p(class = "submodule-desc", tags$em("How to read this: ", "bar length is the feature's loading on DIABLO's first component within its own omics block; a larger absolute value means it contributes more strongly to that component's direction, not that it is individually clinically important.")),
        DT::dataTableOutput(ns("panel_table"))
      )
    })
    panel_plot_fn <- reactive(multi_diablo_panel_plot(dat$panel, top_n = as.integer(input$panel_top_n %||% 20)))
    output$panel_plot <- renderPlot(panel_plot_fn())
    output$dl_panel_png <- multi_png_download(panel_plot_fn, function() sprintf("multiomics_diablo_panel_%s.png", dat$cell$key))
    output$dl_panel_csv <- downloadHandler(function() sprintf("multiomics_diablo_panel_%s.csv", dat$cell$key), function(file) utils::write.csv(dat$panel, file, row.names = FALSE))
    output$panel_table <- DT::renderDataTable({
      req(dat$panel)
      DT::datatable(dat$panel, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    output$variance_ui <- renderUI({
      if (is.null(dat$variance_df)) return(multi_empty_state("No saved DIABLO fit for this cell (variance-explained data unavailable)."))
      tagList(
        p(class = "submodule-desc", "Proportion of each omics block's own variance captured by each DIABLO component - the real, correctly-scoped analog of a \"variance explained\" plot for a supervised sparse discriminant method (mixOmics's own prop_expl_var, not re-derived here). This is not the same quantity as a MOFA2 factor's variance explained across the whole dataset."),
        multi_plot_or_empty(variance_plot_fn, ns("variance_plot"), height = "320px"),
        div(class = "table-toolbar", downloadButton(ns("dl_variance_png"), "Download plot (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_variance_csv"), "Download data (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("variance_table"))
      )
    })
    variance_plot_fn <- reactive(multi_diablo_variance_plot(dat$variance_df))
    output$variance_plot <- renderPlot(variance_plot_fn())
    output$dl_variance_png <- multi_png_download(variance_plot_fn, function() sprintf("multiomics_diablo_variance_%s.png", dat$cell$key))
    output$dl_variance_csv <- downloadHandler(function() sprintf("multiomics_diablo_variance_%s.csv", dat$cell$key), function(file) utils::write.csv(dat$variance_df, file, row.names = FALSE))
    output$variance_table <- DT::renderDataTable({
      req(dat$variance_df)
      DT::datatable(dat$variance_df, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    output$compare_ui <- renderUI({
      if (is.null(dat$cell36) && is.null(dat$baseline36)) return(multi_empty_state("No cross-omics comparison rows for this cell."))
      tagList(
        p(class = "submodule-desc", "Every method run for this cell (integrated) alongside the drug-pooled single-omics baselines - integration is not assumed to win; compare the AUROC/CI directly."),
        DT::dataTableOutput(ns("compare_table"))
      )
    })
    output$compare_table <- DT::renderDataTable({
      combined <- rbind(dat$cell36, dat$baseline36)
      req(combined)
      DT::datatable(combined[, intersect(c("question", "sex", "drug", "omics", "method", "variant", "n", "auroc", "ci_lo", "ci_hi", "excludes_chance"), colnames(combined)), drop = FALSE],
                    rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    observe({
      if (is.null(dat$cell) || is.null(multi_results)) return()
      multi_results$integration <- list(cell = dat$cell, perf = dat$perf, snf_perf = dat$snf_perf, scores = dat$scores, panel = dat$panel)
    })
  })
}
