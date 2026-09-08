## R/multiomics/05_Biomarker_Discovery/mod_multi_biomarker.R
## Biomarker Discovery: live DIABLO (mixOmics::block.splsda) feature selection.

mod_multi_biomarker_config <- list(
  id = "biomarker", title = "Biomarker Discovery", icon = "star", group = "Biomarker modeling",
  description = "Supervised DIABLO feature selection across Transcriptomics + Methylomics."
)

mb_eval_cards <- function(ev) {
  tagList(
    div(style = "display:flex; gap:10px; flex-wrap:wrap;",
        mi_stat_card(sprintf("%.3f", ev$auc), "AUROC"),
        mi_stat_card(sprintf("%.2f-%.2f", ev$ci_lo, ev$ci_hi), "95% CI (DeLong)"),
        mi_stat_card(sprintf("%.3f", ev$ber), "Balanced error rate"),
        mi_stat_card(ev$n, "Samples scored")),
    (function() {
      call_i <- if (!is.null(ev$auroc_call)) ev$auroc_call else multi_auroc_call(ev$auc, ev$ci_lo, ev$ci_hi)
      switch(call_i,
        above_chance = mi_ok("The confidence interval excludes 0.5 in the expected direction (AUROC > 0.5): genuine, statistically supported discrimination."),
        below_chance = mi_stop("The confidence interval is entirely BELOW 0.5: this classifier performs significantly WORSE than random guessing. This is NOT a validation success - it typically indicates a label/orientation error or overfitting to a spuriously anti-correlated small sample, and should be investigated before being reported as a result."),
        mi_warn("The confidence interval includes 0.5 - performance is not distinguishable from chance on these samples.")
      )
    })(),
    p(class = "submodule-desc", sprintf("Classes: %s. Positive class for the score: %s. Prediction distance: %s.",
                                        paste(sprintf("%s = %d", names(ev$n_per_class), as.integer(ev$n_per_class)), collapse = ", "), ev$pos_class, ev$distance))
  )
}

mb_holdout_summary_ui <- function(res) {
  if (is.null(res$holdout)) {
    return(if (isTRUE((res$holdout_frac %||% 0) > 0)) {
      p(class = "submodule-desc", "The hold-out could not be scored.")
    } else {
      mi_warn("No hold-out was set aside for this run (Model tab: \"Hold-out fraction\" = 0). The signature below is association-only - DIABLO's own sparsity-selected features from the analyzed cohort, not an out-of-sample-validated biomarker panel. Raise the hold-out fraction, or run External Validation, before describing it as validated.")
    })
  }
  ev <- res$holdout
  if (!isTRUE(ev$ok)) return(mi_warn(ev$error))
  tagList(
    p(class = "submodule-desc", sprintf("%d training samples were used for selection and cross-validation; %d samples (%.0f%%) were sealed and scored once with the final model.",
                                        length(res$train_ids), length(res$test_ids), 100 * res$holdout_frac)),
    mb_eval_cards(ev)
  )
}

mb_external_summary_ui <- function(ev) {
  if (!isTRUE(ev$ok)) {
    return(tagList(mi_stop(ev$error), if (!is.null(ev$coverage)) mb_coverage_ui(ev$coverage)))
  }
  tagList(
    p(class = "submodule-desc", sprintf("Panel refitted on %d training samples, using the %.0f%% of panel features also measured in the external cohort. Applied to %d external samples (each cohort standardised separately).",
                                        ev$n_train, 100 * ev$overall_coverage, ev$n)),
    if (length(ev$notes) > 0) p(class = "submodule-desc", paste(ev$notes, collapse = " ")),
    mb_coverage_ui(ev$coverage),
    mb_eval_cards(ev)
  )
}

mb_coverage_ui <- function(cov) {
  tags$ul(lapply(names(cov), function(b) tags$li(sprintf("%s: %d of %d panel features measured in the external cohort%s", b, cov[[b]]$n_present, cov[[b]]$n_panel,
                                                        if (length(cov[[b]]$missing) > 0) sprintf(" (missing: %s%s)", paste(head(cov[[b]]$missing, 8), collapse = ", "), if (length(cov[[b]]$missing) > 8) ", ..." else "") else ""))))
}

mod_multi_biomarker_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("active_dataset_banner")),
    tabsetPanel(
      id = ns("tabs"), type = "tabs",
      tabPanel(
        "Setup", br(),
        box(
          width = NULL, title = "1. Data source", status = "primary", solidHeader = FALSE,
          radioButtons(ns("data_source"), "Data source",
                       choices = c("Active Multi-Omics Dataset (Dataset Workspace)" = "active", "Preloaded (RA anti-TNF cohort)" = "preloaded"),
                       selected = "active", inline = TRUE),
          conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
                            selectInput(ns("preloaded_cell"), "Analysis cell", choices = MULTI_CELL_CHOICES, width = "100%")),
          uiOutput(ns("source_note")),
          uiOutput(ns("block_role_ui"))
        ),
        box(
          width = NULL, title = "2. Outcome", status = "primary", solidHeader = FALSE,
          uiOutput(ns("outcome_ui")),
          uiOutput(ns("class_selection_ui"))
        ),
        box(
          width = NULL, title = "3. Feature filtering & scaling", status = "primary", solidHeader = FALSE, collapsible = TRUE,
          p(class = "submodule-desc", "Variance-ranked feature cap per block, applied before cross-validation. Leave at default to keep all features."),
          uiOutput(ns("feature_filter_ui")),
          checkboxInput(ns("scale_blocks"), "Scale each block (mean-center, unit variance)", value = TRUE)
        ),
        box(width = NULL, title = "4. Data check", status = "primary", solidHeader = FALSE, DT::dataTableOutput(ns("data_check_table")))
      ),
      tabPanel("Model", br(), fluidRow(column(4, uiOutput(ns("model_params_ui"))), column(8, uiOutput(ns("planned_run_ui"))))),
      tabPanel("Signature", br(), uiOutput(ns("signature_ui"))),
      tabPanel("Performance", br(), uiOutput(ns("performance_ui"))),
      tabPanel("Stability", br(), uiOutput(ns("stability_ui"))),
      tabPanel("Integration", br(), uiOutput(ns("integration_ui"))),
      tabPanel("Plots", br(), uiOutput(ns("plots_ui"))),
      tabPanel("External Validation", br(), fluidRow(column(4, uiOutput(ns("ext_params_ui"))), column(8, uiOutput(ns("ext_results_ui"))))),
      tabPanel("Sex-Stratified", br(), fluidRow(column(4, uiOutput(ns("ss_params_ui"))), column(8, uiOutput(ns("ss_results_ui")))))
    )
  )
}

mod_multi_biomarker_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$active_dataset_banner <- renderUI(multi_active_dataset_banner(multi_dataset, multi_results))

    mb_dataset <- reactive({
      if (identical(input$data_source, "preloaded")) {
        req(input$preloaded_cell)
        mi_preloaded_cell_dataset(input$preloaded_cell)
      } else {
        if (is.null(multi_dataset) || !isTRUE(multi_dataset$active) || length(multi_dataset$layers %||% list()) < 1) {
          return(list(ok = FALSE, error = "No Active Multi-Omics Dataset yet - build one on the Dataset Workspace tab, or switch to \"Preloaded RA anti-TNF cohort\" above."))
        }
        list(
          ok = TRUE, layers = multi_dataset$layers, sample_meta = multi_dataset$sample_meta, outcome_col = NULL,
          label = sprintf("Active Multi-Omics Dataset (%s)", paste(names(multi_dataset$layers), collapse = " + ")),
          provenance = sprintf("Active Multi-Omics Dataset from the Dataset Workspace tab (source: %s).", multi_dataset$source %||% "unknown"),
          omics_type = stats::setNames(
            vapply(names(multi_dataset$layers), function(nm) multi_dataset$layer_meta[[nm]]$omics_type %||% "other", character(1)),
            names(multi_dataset$layers)
          )
        )
      }
    })

    output$source_note <- renderUI({
      d <- mb_dataset()
      if (!isTRUE(d$ok)) return(mi_warn(d$error))
      tagList(mi_ok(d$provenance),
              if (!identical(input$data_source, "preloaded")) multi_harmonisation_note(multi_results, n_here = tryCatch(mb_val_raw()$n_shared, error = function(e) NULL)))
    })

    output$block_role_ui <- renderUI({
      d <- mb_dataset()
      if (!isTRUE(d$ok)) return(NULL)
      choices <- names(d$layers)
      if (length(choices) < 2) return(mi_stop("DIABLO requires two matched omics blocks (Transcriptomics + Methylomics); only one is available in this dataset."))
      if (!is.null(d$omics_type)) {
        guess_t <- choices[d$omics_type[choices] %in% c("rnaseq", "mirna", "genomics")]
        guess_m <- choices[d$omics_type[choices] == "methylation"]
      } else {
        guess_t <- choices[grepl("transcript|rna|express|gene", choices, ignore.case = TRUE)]
        guess_m <- choices[grepl("methyl|cpg|beta", choices, ignore.case = TRUE)]
      }
      sel_t <- if (length(guess_t) > 0) guess_t[1] else choices[1]
      sel_m <- if (length(guess_m) > 0) guess_m[1] else choices[min(2, length(choices))]
      tagList(
        selectInput(ns("transcript_block"), "Transcriptomics block", choices = choices, selected = sel_t),
        selectInput(ns("methyl_block"), "Methylomics block", choices = choices, selected = sel_m)
      )
    })

    mb_raw_layers <- reactive({
      d <- mb_dataset()
      if (!isTRUE(d$ok)) return(NULL)
      req(input$transcript_block, input$methyl_block)
      if (identical(input$transcript_block, input$methyl_block)) return(NULL)
      mb_select_blocks(d$layers, input$transcript_block, input$methyl_block)
    })

    mb_val_raw <- reactive({
      layers <- mb_raw_layers()
      req(layers)
      mi_validate_dataset(layers, mb_dataset()$sample_meta, outcome_col = NULL)
    })

    outcome_candidates <- reactive({
      d <- mb_dataset()
      if (!isTRUE(d$ok) || is.null(d$sample_meta)) return(character(0))
      if (!is.null(d$outcome_col)) return(d$outcome_col)
      colnames(d$sample_meta)
    })

    output$outcome_ui <- renderUI({
      cands <- outcome_candidates()
      if (length(cands) == 0) return(mi_warn("No sample metadata with a candidate outcome column is available for this dataset."))
      selectInput(ns("outcome_col"), "Outcome variable", choices = cands, selected = multi_harmonisation_outcome_default(multi_results, cands), width = "100%")
    })

    mb_outcome_full <- reactive({
      d <- mb_dataset()
      req(isTRUE(d$ok), !is.null(d$sample_meta), input$outcome_col)
      stats::setNames(d$sample_meta[[input$outcome_col]], rownames(d$sample_meta))
    })

    mb_outcome_summary_raw <- reactive({
      d <- mb_dataset(); v <- mb_val_raw()
      if (!isTRUE(d$ok) || is.null(v) || is.null(input$outcome_col)) return(NULL)
      mi_outcome_summary(d$sample_meta, input$outcome_col, v$shared_ids)
    })

    output$class_selection_ui <- renderUI({
      o <- mb_outcome_summary_raw()
      if (is.null(o) || !identical(o$type, "categorical") || is.null(o$class_counts)) return(NULL)
      cls <- names(o$class_counts)
      labels <- sprintf("%s (n=%d)", mb_friendly_class_label(cls), as.integer(o$class_counts))
      choices <- stats::setNames(cls, labels)
      selectizeInput(ns("classes_selected"), "Classes to include", choices = choices, selected = cls, multiple = TRUE, width = "100%")
    })

    mb_eligible_ids <- reactive({
      v <- mb_val_raw()
      req(v)
      ids <- v$shared_ids
      if (!is.null(input$classes_selected) && length(input$classes_selected) > 0) {
        y <- mb_outcome_full()
        ids <- intersect(ids, names(y)[!is.na(y) & as.character(y) %in% input$classes_selected])
      }
      ids
    })

    output$feature_filter_ui <- renderUI({
      layers <- mb_raw_layers()
      req(layers)
      tagList(
        numericInput(ns("max_features_t"), "Max transcriptomic features (variance-ranked)", value = mb_default_max_features(ncol(layers$Transcriptomics)), min = 10, max = ncol(layers$Transcriptomics), step = 10),
        numericInput(ns("max_features_m"), "Max methylomic features (variance-ranked)", value = mb_default_max_features(ncol(layers$Methylomics)), min = 10, max = ncol(layers$Methylomics), step = 10)
      )
    })

    mb_final_layers <- reactive({
      layers <- mb_raw_layers(); ids <- mb_eligible_ids()
      req(layers, length(ids) > 0)
      max_t <- input$max_features_t %||% mb_default_max_features(ncol(layers$Transcriptomics))
      max_m <- input$max_features_m %||% mb_default_max_features(ncol(layers$Methylomics))
      list(
        Transcriptomics = mb_variance_prefilter(layers$Transcriptomics[ids, , drop = FALSE], max_t),
        Methylomics = mb_variance_prefilter(layers$Methylomics[ids, , drop = FALSE], max_m)
      )
    })

    mb_val <- reactive({
      layers <- tryCatch(mb_final_layers(), error = function(e) NULL)
      if (is.null(layers)) return(NULL)
      mi_validate_dataset(layers, NULL, outcome_col = NULL)
    })
    mb_outcome <- reactive({
      v <- mb_val()
      if (is.null(v) || is.null(input$outcome_col)) return(NULL)
      mi_outcome_summary(mb_dataset()$sample_meta, input$outcome_col, v$shared_ids)
    })
    mb_elig <- reactive({
      v <- mb_val()
      if (is.null(v)) return(list(ok = FALSE, reason = "Select a dataset and both omics blocks above first."))
      mi_diablo_eligibility(v, mb_outcome())
    })

    output$data_check_table <- DT::renderDataTable({
      DT::datatable(mb_data_check_table(mb_val(), mb_outcome(), mb_elig()), rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    output$model_params_ui <- renderUI({
      v <- mb_val(); o <- mb_outcome(); elig <- mb_elig()
      if (is.null(v) || !isTRUE(v$ok)) return(box(width = NULL, title = "Model parameters", status = "primary", solidHeader = FALSE, mi_warn("Complete Setup first (data source, blocks, outcome).")))
      n <- v$n_shared
      loo_ok <- mi_diablo_loo_feasible(n)
      min_class_n <- if (!is.null(o$class_counts)) min(o$class_counts) else max(2, floor(n / 2))

      box(
        width = NULL, title = "Model parameters", status = "primary", solidHeader = FALSE,
        if (!isTRUE(elig$ok)) mi_stop(elig$reason) else mi_ok("DIABLO is available for this dataset."),
        div(style = if (!isTRUE(elig$ok)) "opacity:0.5; pointer-events:none;" else NULL,
          h5("Components"),
          checkboxInput(ns("ncomp_auto"), "Tune the number of components with perf() on a non-sparse block.plsda (mixOmics DIABLO workflow; at least 3 repeats)", value = TRUE),
          numericInput(ns("ncomp"), "Number of components (ncomp)", value = max(mi_diablo_feasible_ncomp(o$n_classes %||% 2, min_class_n)), min = 1, max = max(mi_diablo_feasible_ncomp(o$n_classes %||% 2, min_class_n)), step = 1),
          p(class = "submodule-desc", sprintf("Bounded by outcome classes and the smallest class size for this dataset (feasible: 1-%d).", max(mi_diablo_feasible_ncomp(o$n_classes %||% 2, min_class_n)))),
          hr(),
          h5("Feature selection (keepX)"),
          textInput(ns("keepx_t"), "Transcriptomics keepX per component", value = "20,10"),
          textInput(ns("keepx_m"), "Methylomics keepX per component", value = "20,10"),
          p(class = "submodule-desc", "Comma-separated, one value per component."),
          checkboxInput(ns("keepx_auto"), "Auto-tune keepX instead (mixOmics::tune.block.splsda() grid search - slower)", value = FALSE),
          hr(),
          h5("Block relationship (design)"),
          sliderInput(ns("design_val"), "Transcriptomics <-> Methylomics", min = 0, max = 1, value = 0.1, step = 0.05),
          p(class = "submodule-desc", "Cross-block weight in DIABLO's design matrix - 0.1 is mixOmics's own documented default."),
          hr(),
          h5("Validation"),
          radioButtons(ns("validation_method"), "Method", choices = c("M-fold CV" = "mfold", if (loo_ok) c("Leave-one-out" = "loo")), selected = "mfold", inline = TRUE),
          if (!loo_ok) p(class = "submodule-desc", sprintf("Leave-one-out is disabled for this dataset (n = %d exceeds %d).", n, MI_DIABLO_LOO_MAX_N)),
          conditionalPanel(condition = sprintf("input['%s'] == 'mfold'", ns("validation_method")),
            numericInput(ns("folds"), "Number of folds", value = mi_diablo_feasible_folds(min_class_n), min = 2, max = mi_diablo_max_folds(min_class_n), step = 1),
            numericInput(ns("nrepeat"), "Repeats", value = mi_diablo_feasible_repeats(n), min = 1, max = 50, step = 1)
          ),
          hr(),
          h5("Prediction distance"),
          selectInput(ns("distance"), NULL, choices = c("Automatically selected during tuning" = "automatic", "centroids.dist" = "centroids.dist", "mahalanobis.dist" = "mahalanobis.dist", "max.dist" = "max.dist")),
          numericInput(ns("seed"), "Random seed", value = 1, min = 1, step = 1),
          hr(),
          h5("Sealed hold-out"),
          numericInput(ns("holdout_frac"), "Hold-out fraction (0 = none)", value = 0.2, min = 0, max = 0.5, step = 0.05),
          p(class = "submodule-desc", "A stratified share of the matched samples is set aside before any modelling and scored once with the final model - it never enters feature selection, tuning or cross-validation. Defaults to 20%: without a hold-out (0%), the signature table below is association-only - DIABLO's own sparsity-selected features from the analyzed cohort, not an out-of-sample-validated panel."),
          hr(),
          if (isTRUE(elig$ok)) div(
            if (isTRUE(input$keepx_auto)) mi_warn("Auto-tuning keepX (checked above) may take several minutes for high-dimensional data.") else NULL,
            actionButton(ns("run_btn"), "Run analysis", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")),
          uiOutput(ns("run_status_ui"))
        )
      )
    })

    mb_params <- reactive(list(
      ncomp_mode = if (isTRUE(input$ncomp_auto %||% TRUE)) "tuned" else "manual", ncomp = input$ncomp,
      keepx_mode = if (isTRUE(input$keepx_auto)) "automatic" else "manual",
      keepx_manual = if (!isTRUE(input$keepx_auto)) list(
        Transcriptomics = as.integer(trimws(strsplit(input$keepx_t %||% "10", ",")[[1]])),
        Methylomics = as.integer(trimws(strsplit(input$keepx_m %||% "10", ",")[[1]]))
      ) else NULL,
      design_mode = "custom",
      design_custom = {
        m <- matrix(input$design_val %||% 0.1, 2, 2, dimnames = list(c("Transcriptomics", "Methylomics"), c("Transcriptomics", "Methylomics"))); diag(m) <- 0; m
      },
      validation_mode = "manual", validation_method = input$validation_method %||% "mfold",
      folds = input$folds, nrepeat = input$nrepeat, distance = input$distance %||% "automatic",
      scale = isTRUE(input$scale_blocks %||% TRUE),
      holdout_frac = input$holdout_frac %||% 0
    ))

    output$planned_run_ui <- renderUI({
      v <- mb_val(); req(v)
      p <- mb_params()
      box(width = NULL, title = "Planned run", status = "primary", solidHeader = FALSE,
          tags$ul(
            tags$li(sprintf("Blocks: Transcriptomics + Methylomics | Samples: %d | Outcome: %s", v$n_shared, input$outcome_col %||% "-")),
            tags$li(sprintf("Components: %s", if (identical(p$ncomp_mode, "tuned")) "tuned by perf() (upper bound set above)" else p$ncomp %||% "-")),
            tags$li(sprintf("Feature selection: %s", if (identical(p$keepx_mode, "automatic")) "Auto-tuned (grid search, nested inside every CV fold for the reported performance)" else "Manual (set above)")),
            tags$li(sprintf("Validation: %s, %s-fold x %s repeats", if (identical(p$validation_method, "loo")) "Leave-one-out" else "M-fold CV", p$folds %||% "-", p$nrepeat %||% "-")),
            tags$li(sprintf("Scaling: %s", if (p$scale) "Each block mean-centered, unit variance" else "None")),
            tags$li(sprintf("Sealed hold-out: %s", if (isTRUE((p$holdout_frac %||% 0) > 0)) sprintf("%.0f%% of matched samples, scored once after modelling", 100 * p$holdout_frac) else "none"))
          ),
          multi_empty_state("Click \"Run analysis\" (Model tab, left) to compute the Signature/Performance/Stability/Integration/Plots tabs.")
      )
    })

    mb_run_analysis <- function(layers, outcome, sample_ids, params, seed) {
      split <- mb_holdout_split(outcome, sample_ids, params$holdout_frac %||% 0, seed = seed)
      if (!isTRUE(split$ok)) return(list(ok = FALSE, error = split$error))
      train_ids <- split$train_ids
      set.seed(seed)
      diablo_res <- mi_diablo_run(layers, outcome, train_ids, params)
      if (!isTRUE(diablo_res$ok)) return(list(ok = FALSE, error = diablo_res$error))
      Y <- droplevels(factor(outcome[train_ids]))
      cv_roc <- if (nlevels(Y) == 2) {
        Xm <- lapply(layers, function(m) m[train_ids, , drop = FALSE])
        mb_cv_roc(Xm, Y, diablo_res$params, seed = seed)
      } else NULL
      holdout <- if (length(split$test_ids) > 0) mb_holdout_evaluate(diablo_res, layers, outcome, split$test_ids) else NULL
      list(ok = TRUE, diablo = diablo_res, cv_roc = cv_roc, holdout = holdout,
           train_ids = train_ids, test_ids = split$test_ids, holdout_frac = split$frac)
    }

    mb_state <- reactiveValues(result = NULL, error = NULL, submitted = FALSE, outcome_used = NULL, layers_used = NULL, dataset_label = NULL)

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      mb_task <- ExtendedTask$new(function(layers, outcome, ids, params, seed) {
        promises::future_promise(mb_run_analysis(layers, outcome, ids, params, seed), seed = TRUE)
      })
      observeEvent(input$run_btn, {
        v <- req(mb_val()); layers <- req(mb_final_layers())
        validate(need(isTRUE(mb_elig()$ok), mb_elig()$reason))
        outcome <- mb_outcome_full()
        mb_state$error <- NULL; mb_state$result <- NULL; mb_state$submitted <- TRUE
        mb_state$outcome_used <- outcome; mb_state$layers_used <- layers; mb_state$dataset_label <- mb_dataset()$label
        mb_state$preselected_input <- identical(input$data_source, "preloaded") || identical(multi_dataset$source %||% "", "preloaded")
        p <- mb_params(); ids <- v$shared_ids; seed <- input$seed %||% 1
        session$onFlushed(function() mb_task$invoke(layers, outcome, ids, p, seed), once = TRUE)
      })
      observe({
        res <- tryCatch(mb_task$result(), error = function(e) e)
        if (inherits(res, "shiny.silent.error")) return()
        mb_state$submitted <- FALSE
        if (inherits(res, "error")) { mb_state$error <- conditionMessage(res); return() }
        if (!isTRUE(res$ok)) { mb_state$error <- res$error; return() } else { mb_state$result <- res }
      })
      mb_running <- reactive(identical(mb_task$status(), "running") || isTRUE(mb_state$submitted))
    } else {
      observeEvent(input$run_btn, {
        validate(need(isTRUE(mb_elig()$ok), mb_elig()$reason))
        showNotification("Running DIABLO synchronously - the app will be briefly unresponsive.", type = "message", duration = 5)
        v <- req(mb_val()); layers <- req(mb_final_layers())
        outcome <- mb_outcome_full()
        mb_state$outcome_used <- outcome; mb_state$layers_used <- layers; mb_state$dataset_label <- mb_dataset()$label
        mb_state$preselected_input <- identical(input$data_source, "preloaded") || identical(multi_dataset$source %||% "", "preloaded")
        res <- mb_run_analysis(layers, outcome, v$shared_ids, mb_params(), input$seed %||% 1)
        mb_state$error <- if (!isTRUE(res$ok)) res$error else NULL
        if (isTRUE(res$ok)) mb_state$result <- res
      })
      mb_running <- reactive(FALSE)
    }

    observe({
      if (isTRUE(mb_running())) {
        shinyjs::disable(ns("run_btn")); shinyjs::html(ns("run_btn"), as.character(tagList(icon("spinner", class = "fa-spin"), " Running...")))
      } else {
        shinyjs::enable(ns("run_btn")); shinyjs::html(ns("run_btn"), as.character(tagList(icon("play"), " Run analysis")))
      }
    })

    output$run_status_ui <- renderUI({
      if (!isTRUE(input$run_btn > 0)) return(NULL)
      if (isTRUE(mb_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running - see the Signature/Performance/Stability tabs once this finishes."))
      if (!is.null(mb_state$error)) return(mi_stop(mb_state$error))
      if (!is.null(mb_state$result)) return(mi_ok("Finished - see the Signature/Performance/Stability/Integration/Plots tabs."))
      NULL
    })

    gate_ui <- function(body_fn) {
      if (!isTRUE(input$run_btn > 0)) return(multi_empty_state("Set parameters and click \"Run analysis\" (Model tab) to see results here."))
      if (isTRUE(mb_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running - the rest of the app stays usable while this runs."))
      if (!is.null(mb_state$error)) return(mi_stop(mb_state$error))
      if (is.null(mb_state$result)) return(multi_empty_state())
      body_fn(mb_state$result)
    }

    sig_df <- reactive({ req(mb_state$result); mb_signature_table(mb_state$result$diablo) })

    output$signature_ui <- renderUI(gate_ui(function(res) {
      sig <- sig_df()
      if (is.null(sig)) return(mi_warn("No features were selected by this model."))
      n_t <- length(unique(sig$feature[sig$omics == "Transcriptomics"]))
      n_m <- length(unique(sig$feature[sig$omics == "Methylomics"]))
      tagList(
        div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px;",
            mi_stat_card(n_t, "Transcriptomic features selected", ARTHOMIX_COLORS$violet),
            mi_stat_card(n_m, "Methylomic features selected", ARTHOMIX_COLORS$orange),
            mi_stat_card(res$diablo$params$ncomp, "Components"), mi_stat_card(res$diablo$params$n_samples, "Samples")),
        div(class = "table-toolbar", downloadButton(ns("dl_signature"), "Download signature (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("signature_table"))
      )
    }))
    output$signature_table <- DT::renderDataTable({
      sig <- req(sig_df())
      cols <- intersect(c("omics", "feature", "component", "loading", "rank_within_block", "selection_frequency", "stability_category"), colnames(sig))
      DT::datatable(sig[, cols, drop = FALSE], rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_signature <- downloadHandler(function() "biomarker_discovery_signature.csv", function(file) utils::write.csv(req(sig_df()), file, row.names = FALSE))

    output$performance_ui <- renderUI(gate_ui(function(res) {
      perf_sum <- mi_diablo_performance_summary(res$diablo)
      p <- res$diablo$params
      tagList(
        box(width = NULL, title = "Actual parameters used", status = "primary", solidHeader = FALSE,
            tags$ul(
              tags$li(sprintf("Components: %d (%s)", p$ncomp, if (identical(p$ncomp_mode, "tuned")) p$ncomp_note %||% "tuned" else "set manually")),
              tags$li(sprintf("keepX - Transcriptomics: %s | Methylomics: %s", paste(p$keepX$Transcriptomics, collapse = ","), paste(p$keepX$Methylomics, collapse = ","))),
              tags$li(sprintf("Prediction distance: %s", p$distance)),
              tags$li(sprintf("Validation: %s, %d-fold x %d repeat(s)", p$validation_method, p$folds, p$nrepeat)),
              tags$li(sprintf("Reported error rates and AUROC: %s", perf_sum$error_source %||% "-")),
              if (isTRUE(p$loo_downgraded)) tags$li(mi_warn("Leave-one-out was requested but this dataset is too large for it - M-fold CV was used instead."))
            )),
        box(width = NULL, title = "Cross-validated performance (held-out folds only - never training performance)", status = "primary", solidHeader = FALSE,
            if (is.null(perf_sum)) mi_warn("Performance could not be computed for this configuration.") else tagList(
              div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                  mi_stat_card(sprintf("%.3f", perf_sum$ber), "Balanced error rate (BER)"),
                  mi_stat_card(sprintf("%.3f", perf_sum$overall_error), "Overall error rate")),
              multi_plot_or_empty(function() mi_diablo_error_bar_plot(perf_sum), ns("perf_error_plot"), height = "260px"),
              div(class = "table-toolbar", downloadButton(ns("dl_perf_error_png"), "Download plot (PNG)", class = "btn-sm")),
              if (!is.null(perf_sum$auc)) tagList(h5("AUROC (pooled out-of-fold predictions, per-fold refits; DeLong 95% CI)"), DT::dataTableOutput(ns("perf_auc_table")))
            )),
        box(width = NULL, title = "Cross-validated ROC (pooled out-of-fold predictions)", status = "primary", solidHeader = FALSE,
            if (is.null(res$cv_roc)) mi_warn("Not available - ROC requires a binary outcome (this run has a different number of classes).")
            else tagList(
              multi_plot_or_empty(function() mb_roc_plot(res$cv_roc), ns("perf_roc_plot"), height = "380px"),
              div(class = "table-toolbar", downloadButton(ns("dl_perf_roc_png"), "Download plot (PNG)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Sealed hold-out evaluation (samples never seen during selection, tuning or CV)", status = "primary", solidHeader = FALSE,
            if (isTRUE(mb_state$preselected_input) && !is.null(res$holdout)) mi_warn("The preloaded input is the 50-gene/100-CpG panel DIABLO already selected on ALL samples of this cell. A hold-out drawn from it isn't independent, so its performance is optimistic, not real evidence. Use a full-matrix dataset (Upload or GEO) for a valid hold-out."),
            mb_holdout_summary_ui(res)),
        box(width = NULL, title = "Reproducibility - parameters actually used", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("repro_table")),
            div(class = "table-toolbar",
                downloadButton(ns("dl_params"), "Download model parameters (CSV)", class = "btn-sm"),
                downloadButton(ns("dl_matched_samples"), "Download matched sample list (CSV)", class = "btn-sm"),
                downloadButton(ns("dl_performance"), "Download performance table (CSV)", class = "btn-sm"))),
        box(width = NULL, title = "Software versions", status = "primary", solidHeader = FALSE, DT::dataTableOutput(ns("repro_versions_table")))
      )
    }))
    output$perf_auc_table <- DT::renderDataTable({
      DT::datatable(mi_diablo_performance_summary(req(mb_state$result)$diablo)$auc, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$dl_performance <- downloadHandler(function() "biomarker_discovery_performance.csv", function(file) {
      res <- req(mb_state$result); perf_sum <- mi_diablo_performance_summary(res$diablo)
      df <- data.frame(metric = c("BER", "Overall error", names(perf_sum$per_class_error)), value = c(perf_sum$ber, perf_sum$overall_error, as.numeric(perf_sum$per_class_error)))
      utils::write.csv(df, file, row.names = FALSE)
    })
    output$perf_error_plot <- multi_render_plotly(function() { mi_diablo_error_bar_plot(mi_diablo_performance_summary(req(mb_state$result)$diablo)) })
    output$dl_perf_error_png <- multi_png_download(function() { mi_diablo_error_bar_plot(mi_diablo_performance_summary(req(mb_state$result)$diablo)) }, function() "biomarker_discovery_error_rate.png")
    output$perf_roc_plot <- multi_render_plotly(function() { mb_roc_plot(req(req(mb_state$result)$cv_roc)) })
    output$dl_perf_roc_png <- multi_png_download(function() { mb_roc_plot(req(req(mb_state$result)$cv_roc)) }, function() "biomarker_discovery_cv_roc.png")

    output$stability_ui <- renderUI(gate_ui(function(res) {
      sig <- sig_df()
      if (is.null(sig)) return(mi_warn("No selected features to assess."))
      dist <- as.data.frame(table(stability = sig$stability_category))
      tagList(
        p(class = "submodule-desc", sprintf(
          "Stable: >=%.0f%% of CV repeats. Moderate: %.0f-%.0f%%. Low: <%.0f%%. (Different from SNF Clustering's own chance-corrected >=75%% stability metric.)",
          MB_STABILITY_THRESHOLDS$stable * 100, MB_STABILITY_THRESHOLDS$moderate * 100, MB_STABILITY_THRESHOLDS$stable * 100, MB_STABILITY_THRESHOLDS$moderate * 100
        )),
        div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px;",
            lapply(seq_len(nrow(dist)), function(i) mi_stat_card(dist$Freq[i], as.character(dist$stability[i])))),
        multi_plot_or_empty(function() mb_stability_plot(sig), ns("stability_plot"), height = "460px"),
        div(class = "table-toolbar", downloadButton(ns("dl_stability_png"), "Download plot (PNG)", class = "btn-sm")),
        h5("Stability by omics block"),
        DT::dataTableOutput(ns("stability_by_block_table")),
        div(class = "table-toolbar", downloadButton(ns("dl_stability"), "Download stability table (CSV)", class = "btn-sm"))
      )
    }))
    output$stability_plot <- multi_render_plotly(function() { mb_stability_plot(req(sig_df())) })
    output$dl_stability_png <- multi_png_download(function() { mb_stability_plot(req(sig_df())) }, function() "biomarker_discovery_stability.png")
    output$stability_by_block_table <- DT::renderDataTable({
      sig <- req(sig_df())
      tab <- as.data.frame(table(omics = sig$omics, stability = sig$stability_category))
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$dl_stability <- downloadHandler(function() "biomarker_discovery_stability.csv", function(file) {
      sig <- req(sig_df())
      utils::write.csv(sig[, intersect(c("omics", "feature", "component", "selection_frequency", "stability_category"), colnames(sig)), drop = FALSE], file, row.names = FALSE)
    })

    output$integration_ui <- renderUI(gate_ui(function(res) {
      fit <- res$diablo$fit; p <- res$diablo$params
      cc <- mb_component_correlation(fit)
      sel <- mi_diablo_selected_features_df(fit)
      corr <- if (!is.null(sel)) mi_diablo_selected_correlation_data(mb_state$layers_used, sel, "Transcriptomics", "Methylomics", rownames(fit$variates$Transcriptomics)) else list(ok = FALSE)
      tagList(
        box(width = NULL, title = "DIABLO integration structure", status = "primary", solidHeader = FALSE,
            tags$ul(
              tags$li(sprintf("Design (block relationship weight): %.2f", p$design["Transcriptomics", "Methylomics"])),
              tags$li(sprintf("Design mode: %s", p$design_mode))
            )),
        box(width = NULL, title = "Component correlation (block scores, component 1)", status = "primary", solidHeader = FALSE,
            if (is.null(cc)) mi_warn("Not available.") else tagList(
              mi_stat_card(sprintf("%.2f", cc$r), "Correlation (r)"),
              multi_plot_or_empty(function() mb_component_correlation_plot(fit, mb_state$outcome_used), ns("int_component_plot"), height = "360px"),
              div(class = "table-toolbar", downloadButton(ns("dl_int_component_png"), "Download plot (PNG)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Selected cross-omics feature relationships", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Correlation between selected Transcriptomics and Methylomics features (matched samples only)."),
            if (!isTRUE(corr$ok)) mi_warn(corr$error %||% "Not available for this run.")
            else tagList(
              multi_plot_or_empty(function() multi_live_correlation_heatmap_plot(corr$df), ns("int_corr_plot"), height = "420px"),
              div(class = "table-toolbar", downloadButton(ns("dl_int_corr_png"), "Download plot (PNG)", class = "btn-sm"))
            ))
      )
    }))
    output$int_component_plot <- multi_render_plotly(function() { mb_component_correlation_plot(req(mb_state$result)$diablo$fit, mb_state$outcome_used) })
    output$dl_int_component_png <- multi_png_download(function() { mb_component_correlation_plot(req(mb_state$result)$diablo$fit, mb_state$outcome_used) }, function() "biomarker_discovery_component_correlation.png")
    output$int_corr_plot <- multi_render_plotly(function() {
      res <- req(mb_state$result)
      sel <- req(mi_diablo_selected_features_df(res$diablo$fit))
      corr <- mi_diablo_selected_correlation_data(mb_state$layers_used, sel, "Transcriptomics", "Methylomics", rownames(res$diablo$fit$variates$Transcriptomics))
      req(isTRUE(corr$ok))
      multi_live_correlation_heatmap_plot(corr$df)
    })
    output$dl_int_corr_png <- multi_png_download(function() {
      res <- req(mb_state$result)
      sel <- req(mi_diablo_selected_features_df(res$diablo$fit))
      corr <- mi_diablo_selected_correlation_data(mb_state$layers_used, sel, "Transcriptomics", "Methylomics", rownames(res$diablo$fit$variates$Transcriptomics))
      req(isTRUE(corr$ok))
      multi_live_correlation_heatmap_plot(corr$df)
    }, function() "biomarker_discovery_cross_omics_correlation.png")

    output$plots_ui <- renderUI(gate_ui(function(res) {
      fit <- res$diablo$fit
      sel <- mi_diablo_selected_features_df(fit)
      scores_df <- mi_diablo_sample_scores_df(fit, mb_state$outcome_used)
      ids <- rownames(fit$variates$Transcriptomics)
      tagList(
        box(width = NULL, title = "Sample plot (component scores, colored by outcome)", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() multi_diablo_score_plot(scores_df), ns("plot_sample"), height = "320px"),
            div(class = "table-toolbar", downloadButton(ns("dl_plot_sample_png"), "Download plot (PNG)", class = "btn-sm"))),
        box(width = NULL, title = "Feature loadings (top selected features)", status = "primary", solidHeader = FALSE,
            if (is.null(sel)) mi_warn("No features were selected.") else tagList(
              multi_plot_or_empty(function() multi_diablo_panel_plot(mi_diablo_panel_df_for_plot(sel, 1)), ns("plot_loadings"), height = "380px"),
              div(class = "table-toolbar", downloadButton(ns("dl_plot_loadings_png"), "Download plot (PNG)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Feature-selection stability", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mb_stability_plot(sig_df()), ns("plot_stability"), height = "420px"),
            div(class = "table-toolbar", downloadButton(ns("dl_plot_stability_png"), "Download plot (PNG)", class = "btn-sm"))),
        box(width = NULL, title = "Selected features x samples (annotated by outcome)", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mb_heatmap_plot(mb_state$layers_used, sig_df(), mb_state$outcome_used, ids), ns("plot_heatmap"), height = "460px"),
            div(class = "table-toolbar", downloadButton(ns("dl_plot_heatmap_png"), "Download plot (PNG)", class = "btn-sm"))),
        box(width = NULL, title = "Variance explained (within block)", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() multi_diablo_variance_plot(multi_diablo_variance_df(fit)), ns("plot_variance"), height = "300px"),
            div(class = "table-toolbar", downloadButton(ns("dl_plot_variance_png"), "Download plot (PNG)", class = "btn-sm"))),
        if (!is.null(res$cv_roc)) box(width = NULL, title = "Cross-validated ROC", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mb_roc_plot(res$cv_roc), ns("plot_roc"), height = "380px"),
            div(class = "table-toolbar", downloadButton(ns("dl_plot_roc_png"), "Download plot (PNG)", class = "btn-sm"))),
        p(class = "submodule-desc", "Circos-style plot not implemented - see \"Selected cross-omics feature relationships\" on the Integration tab for the same data as a heatmap.")
      )
    }))
    output$plot_sample <- multi_render_plotly(function() { multi_diablo_score_plot(mi_diablo_sample_scores_df(req(mb_state$result)$diablo$fit, mb_state$outcome_used)) })
    output$dl_plot_sample_png <- multi_png_download(function() { multi_diablo_score_plot(mi_diablo_sample_scores_df(req(mb_state$result)$diablo$fit, mb_state$outcome_used)) }, function() "biomarker_discovery_sample_plot.png")
    output$plot_loadings <- multi_render_plotly(function() { multi_diablo_panel_plot(mi_diablo_panel_df_for_plot(req(mi_diablo_selected_features_df(req(mb_state$result)$diablo$fit)), 1)) })
    output$dl_plot_loadings_png <- multi_png_download(function() { multi_diablo_panel_plot(mi_diablo_panel_df_for_plot(req(mi_diablo_selected_features_df(req(mb_state$result)$diablo$fit)), 1)) }, function() "biomarker_discovery_feature_loadings.png")
    output$plot_stability <- multi_render_plotly(function() { mb_stability_plot(req(sig_df())) })
    output$dl_plot_stability_png <- multi_png_download(function() { mb_stability_plot(req(sig_df())) }, function() "biomarker_discovery_feature_stability.png")
    output$plot_heatmap <- multi_render_plotly(function() {
      res <- req(mb_state$result)
      mb_heatmap_plot(mb_state$layers_used, req(sig_df()), mb_state$outcome_used, rownames(res$diablo$fit$variates$Transcriptomics))
    })
    output$dl_plot_heatmap_png <- multi_png_download(function() {
      res <- req(mb_state$result)
      mb_heatmap_plot(mb_state$layers_used, req(sig_df()), mb_state$outcome_used, rownames(res$diablo$fit$variates$Transcriptomics))
    }, function() "biomarker_discovery_feature_sample_heatmap.png")
    output$plot_variance <- multi_render_plotly(function() { multi_diablo_variance_plot(multi_diablo_variance_df(req(mb_state$result)$diablo$fit)) })
    output$dl_plot_variance_png <- multi_png_download(function() { multi_diablo_variance_plot(multi_diablo_variance_df(req(mb_state$result)$diablo$fit)) }, function() "biomarker_discovery_variance_explained.png")
    output$plot_roc <- multi_render_plotly(function() { mb_roc_plot(req(req(mb_state$result)$cv_roc)) })
    output$dl_plot_roc_png <- multi_png_download(function() { mb_roc_plot(req(req(mb_state$result)$cv_roc)) }, function() "biomarker_discovery_cv_roc.png")

    output$repro_table <- DT::renderDataTable({
      res <- req(mb_state$result)
      preprocessing_note <- sprintf(
        "Variance prefilter before cross-validation (Transcriptomics <= %s features, Methylomics <= %s features); %s",
        format(input$max_features_t %||% NA, big.mark = ","), format(input$max_features_m %||% NA, big.mark = ","),
        if (isTRUE(mb_params()$scale)) "each block mean-centered and scaled to unit variance inside DIABLO." else "no additional scaling applied."
      )
      DT::datatable(mb_reproducibility_table(res$diablo, mb_state$dataset_label, preprocessing_note, input$seed %||% 1), rownames = FALSE, options = list(dom = "t", pageLength = 20), class = "stripe hover compact")
    })
    output$repro_versions_table <- DT::renderDataTable({
      DT::datatable(multi_package_versions(), rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$dl_params <- downloadHandler(function() "biomarker_discovery_parameters.csv", function(file) {
      res <- req(mb_state$result)
      utils::write.csv(mb_reproducibility_table(res$diablo, mb_state$dataset_label, "See app for details.", input$seed %||% 1), file, row.names = FALSE)
    })
    output$dl_matched_samples <- downloadHandler(function() "biomarker_discovery_matched_samples.csv", function(file) {
      res <- req(mb_state$result)
      ids <- rownames(res$diablo$fit$variates$Transcriptomics)
      utils::write.csv(mb_matched_sample_table(mb_val(), mb_state$outcome_used, ids), file, row.names = FALSE)
    })

    ss_state <- reactiveValues(result = NULL, error = NULL, submitted = FALSE)

    ## ---- External validation: transfer the selected panel to an independent cohort ----
    ext_state <- reactiveValues(result = NULL, layers = NULL, meta = NULL, notes = character(0))

    ext_read_layers <- reactive({
      fe <- input$ext_expr_file; fm <- input$ext_meth_file; fs <- input$ext_meta_file
      req(fe, fm, fs)
      e <- mb_read_external_layer(fe$datapath, fe$name); m <- mb_read_external_layer(fm$datapath, fm$name)
      meta <- mo_read_meta_file(fs)
      errs <- c(if (!isTRUE(e$ok)) sprintf("Expression: %s", e$error), if (!isTRUE(m$ok)) sprintf("Methylation: %s", m$error), if (is.null(meta)) "Metadata: could not be read as a table.")
      if (length(errs) > 0) return(list(ok = FALSE, error = paste(errs, collapse = " ")))
      layers <- list(Transcriptomics = e$mat, Methylomics = m$mat)
      ov <- multi_live_sample_overlap(layers)
      chk <- multi_live_same_patient_check(ov)
      if (!isTRUE(chk$ok)) return(list(ok = FALSE, error = paste(chk$message, chk$detail)))
      ids <- intersect(ov$shared_ids, rownames(meta))
      if (length(ids) < 6) return(list(ok = FALSE, error = "Fewer than 6 external samples are present in both layers and the metadata sheet."))
      list(ok = TRUE, layers = lapply(ov$mats, function(x) x[ids, , drop = FALSE]), meta = meta[ids, , drop = FALSE], n = length(ids))
    })

    output$ext_params_ui <- renderUI({
      box(width = NULL, title = "External cohort", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Upload an independent cohort with matched expression and methylation from the same individuals. First column must be sample ID in every file, with the same patient IDs across all three files. The panel is refitted on the discovery training samples and scored here - nothing from this cohort enters feature selection."),
          fileInput(ns("ext_expr_file"), "External expression matrix (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
          fileInput(ns("ext_meth_file"), "External methylation matrix (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
          fileInput(ns("ext_meta_file"), "External sample metadata (CSV, first column = sample ID)", accept = c(".csv")),
          uiOutput(ns("ext_outcome_ui")),
          numericInput(ns("ext_min_coverage"), "Minimum share of panel features that must be measured in the external cohort", value = MB_TRANSFER_MIN_COVERAGE, min = 0.1, max = 1, step = 0.05),
          uiOutput(ns("ext_run_ui")))
    })

    output$ext_outcome_ui <- renderUI({
      rd <- tryCatch(ext_read_layers(), error = function(e) NULL)
      if (is.null(rd)) return(p(class = "submodule-desc", "Upload all three files to continue."))
      if (!isTRUE(rd$ok)) return(mi_stop(rd$error))
      cols <- colnames(rd$meta)
      guess <- cols[grepl("outcome|response|group|status|flare|class", cols, ignore.case = TRUE)]
      tagList(
        mi_ok(sprintf("%d external samples with both layers and metadata.", rd$n)),
        selectInput(ns("ext_outcome_col"), "External outcome column", choices = cols, selected = if (length(guess) > 0) guess[1] else cols[1]),
        uiOutput(ns("ext_level_map_ui")))
    })

    output$ext_level_map_ui <- renderUI({
      rd <- tryCatch(ext_read_layers(), error = function(e) NULL); req(rd, isTRUE(rd$ok), input$ext_outcome_col)
      res <- mb_state$result; req(res)
      classes <- levels(res$diablo$fit$Y)
      lv <- sort(unique(as.character(rd$meta[[input$ext_outcome_col]])))
      lv <- lv[!is.na(lv) & nzchar(lv)]
      tagList(
        p(class = "submodule-desc", sprintf("Map the external levels onto the discovery classes (positive class = \"%s\", negative class = \"%s\"). Unmapped levels are dropped.", classes[2], classes[1])),
        checkboxGroupInput(ns("ext_pos_levels"), sprintf("External levels treated as \"%s\"", classes[2]), choices = lv, selected = if (classes[2] %in% lv) classes[2] else NULL),
        checkboxGroupInput(ns("ext_neg_levels"), sprintf("External levels treated as \"%s\"", classes[1]), choices = lv, selected = if (classes[1] %in% lv) classes[1] else NULL))
    })

    output$ext_run_ui <- renderUI({
      if (is.null(mb_state$result)) return(mi_warn("Run the discovery analysis (Model tab) first - the external cohort is scored with that panel."))
      actionButton(ns("ext_run_btn"), "Validate panel on external cohort", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")
    })

    observeEvent(input$ext_run_btn, {
      res <- mb_state$result
      if (is.null(res)) { showNotification("Run the discovery analysis first.", type = "error"); return() }
      rd <- tryCatch(ext_read_layers(), error = function(e) NULL)
      if (is.null(rd) || !isTRUE(rd$ok)) { showNotification(rd$error %||% "Upload the three external files first.", type = "error"); return() }
      classes <- levels(res$diablo$fit$Y)
      if (length(input$ext_pos_levels) == 0 || length(input$ext_neg_levels) == 0 || length(intersect(input$ext_pos_levels, input$ext_neg_levels)) > 0) {
        showNotification("Assign at least one external level to each discovery class, with no level in both.", type = "error"); return()
      }
      ext_outcome <- mb_map_external_outcome(rd$meta, input$ext_outcome_col, input$ext_pos_levels, input$ext_neg_levels, classes[2], classes[1])
      showNotification("Refitting the panel and scoring the external cohort...", type = "message", duration = 5)
      ext_state$result <- mb_panel_transfer(res$diablo, mb_state$layers_used, res$train_ids, mb_state$outcome_used, rd$layers, ext_outcome,
                                            min_coverage = input$ext_min_coverage %||% MB_TRANSFER_MIN_COVERAGE)
      ext_state$result$external_label <- sprintf("%s + %s", input$ext_expr_file$name, input$ext_meth_file$name)
      multi_results$biomarker_external <- ext_state$result
    })

    observeEvent(mb_state$result, { ext_state$result <- NULL })

    output$ext_results_ui <- renderUI({
      if (is.null(mb_state$result)) return(multi_empty_state("Run the discovery analysis first, then upload an external cohort on the left."))
      if (is.null(ext_state$result)) return(multi_empty_state("Upload the external cohort files, map its outcome levels and click \"Validate panel on external cohort\"."))
      ev <- ext_state$result
      tagList(
        box(width = NULL, title = sprintf("Panel transfer to external cohort: %s", ev$external_label %||% ""), status = "primary", solidHeader = FALSE,
            mb_external_summary_ui(ev)),
        if (isTRUE(ev$ok)) box(width = NULL, title = "External ROC", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mb_roc_plot(list(roc = ev$roc, auc = ev$auc, n_used = ev$n, n_total = ev$n, folds = NA, pos_class = ev$pos_class, neg_class = ev$neg_class)), ns("ext_roc_plot"), height = "360px"),
            div(class = "table-toolbar", downloadButton(ns("dl_external_scores"), "Download external sample scores (CSV)", class = "btn-sm")))
      )
    })
    output$dl_external_scores <- downloadHandler(function() "biomarker_external_validation_scores.csv", function(file) {
      ev <- req(ext_state$result); req(isTRUE(ev$ok)); utils::write.csv(ev$scores, file, row.names = FALSE)
    })

    ss_sex_col <- reactive({
      d <- mb_dataset()
      if (!isTRUE(d$ok) || is.null(d$sample_meta)) return(NULL)
      cands <- multi_sex_candidates(d$sample_meta)
      if (length(cands) == 0) NULL else cands[1]
    })

    ss_covariate_choices <- reactive({
      d <- mb_dataset()
      if (!isTRUE(d$ok) || is.null(d$sample_meta)) return(character(0))
      setdiff(colnames(d$sample_meta), c(input$outcome_col, ss_sex_col()))
    })

    output$ss_params_ui <- renderUI({
      layers <- tryCatch(mb_raw_layers(), error = function(e) NULL)
      if (is.null(layers)) return(box(width = NULL, title = "Sex-Stratified parameters", status = "primary", solidHeader = FALSE, mi_warn("Complete Setup first (data source, blocks, outcome).")))
      sex_col <- ss_sex_col()
      cov_choices <- ss_covariate_choices()
      cov_guess <- cov_choices[grepl("drug|treatment|therapy", cov_choices, ignore.case = TRUE)][1]

      box(
        width = NULL, title = "Sex-Stratified parameters", status = "primary", solidHeader = FALSE,
        if (is.null(sex_col)) mi_warn("No Sex/Gender column detected - \"All (pooled)\" still works.")
        else mi_ok(sprintf("Sex column detected: \"%s\".", sex_col)),
        hr(),
        h5("Engine"),
        radioButtons(ns("ss_engine"), NULL, choices = MSS_ENGINE_CHOICES, selected = "diablo"),
        hr(),
        h5("Sex stratification"),
        radioButtons(ns("ss_sex_mode"), NULL, choices = MSS_SEX_MODE_CHOICES, selected = "pooled"),
        hr(),
        h5("Adjustment covariate"),
        selectInput(ns("ss_covariate"), NULL, choices = c("(none)" = "", stats::setNames(cov_choices, cov_choices)), selected = cov_guess %||% ""),
        p(class = "submodule-desc", "Adjusted for during in-fold feature selection (e.g. drug). Optional."),
        hr(),
        h5("Parameters (defaults = the pipeline's own constants)"),
        fluidRow(
          column(6, numericInput(ns("ss_top_expr"), "Top expression features", value = MSS_DEFAULTS$top_expr, min = 5, max = 2000, step = 5)),
          column(6, numericInput(ns("ss_top_meth"), "Top methylation features", value = MSS_DEFAULTS$top_meth, min = 5, max = 2000, step = 5))
        ),
        numericInput(ns("ss_ncomp"), "DIABLO components (ncomp)", value = MSS_DEFAULTS$ncomp, min = 1, max = 5, step = 1),
        fluidRow(
          column(6, numericInput(ns("ss_folds"), "CV folds", value = MSS_DEFAULTS$folds, min = 2, max = 10, step = 1)),
          column(6, numericInput(ns("ss_repeats"), "CV repeats", value = MSS_DEFAULTS$repeats, min = 1, max = 20, step = 1))
        ),
        fluidRow(
          column(6, numericInput(ns("ss_keepx_min"), "keepX min", value = MSS_DEFAULTS$keepx_min, min = 2, max = 50, step = 1)),
          column(6, numericInput(ns("ss_keepx_max"), "keepX max", value = MSS_DEFAULTS$keepx_max, min = 2, max = 100, step = 1))
        ),
        conditionalPanel(condition = sprintf("input['%s'] == 'rf'", ns("ss_engine")),
                          numericInput(ns("ss_ntree"), "Random Forest ntree", value = MSS_DEFAULTS$ntree, min = 100, max = 2000, step = 100)),
        hr(),
        actionButton(ns("ss_run_btn"), "Run Sex-Stratified Analysis", icon = icon("play"), class = "btn-primary btn-sm", width = "100%"),
        uiOutput(ns("ss_run_status_ui"))
      )
    })

    ss_params <- reactive(list(
      top_expr = input$ss_top_expr %||% MSS_DEFAULTS$top_expr, top_meth = input$ss_top_meth %||% MSS_DEFAULTS$top_meth,
      top_expr_rf = input$ss_top_expr %||% MSS_DEFAULTS$top_expr_rf, top_meth_rf = input$ss_top_meth %||% MSS_DEFAULTS$top_meth_rf,
      ncomp = input$ss_ncomp %||% MSS_DEFAULTS$ncomp,
      folds = input$ss_folds %||% MSS_DEFAULTS$folds, repeats = input$ss_repeats %||% MSS_DEFAULTS$repeats,
      keepx_min = input$ss_keepx_min %||% MSS_DEFAULTS$keepx_min, keepx_max = input$ss_keepx_max %||% MSS_DEFAULTS$keepx_max,
      ntree = input$ss_ntree %||% MSS_DEFAULTS$ntree,
      min_selected_diablo = MSS_DEFAULTS$min_selected_diablo, min_selected_rf = MSS_DEFAULTS$min_selected_rf,
      min_train_rows_rf = MSS_DEFAULTS$min_train_rows_rf
    ))

    run_ss <- function() {
      layers <- req(mb_raw_layers()); d <- req(mb_dataset())
      mss_run_stratified(
        layers$Transcriptomics, layers$Methylomics, d$sample_meta, input$outcome_col,
        covariate_col = if (nzchar(input$ss_covariate %||% "")) input$ss_covariate else NULL,
        sex_mode = input$ss_sex_mode %||% "pooled", engine = input$ss_engine %||% "diablo", params = ss_params()
      )
    }

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      ss_task <- ExtendedTask$new(function(expr, meth, meta, outcome_col, covariate_col, sex_mode, engine, params) {
        promises::future_promise(mss_run_stratified(expr, meth, meta, outcome_col, covariate_col, sex_mode, engine, params), seed = TRUE)
      })
      observeEvent(input$ss_run_btn, {
        layers <- req(mb_raw_layers()); d <- req(mb_dataset())
        req(input$outcome_col)
        sex_mode <- input$ss_sex_mode %||% "pooled"
        validate(need(identical(sex_mode, "pooled") || !is.null(ss_sex_col()), "No Sex/Gender column in this dataset's metadata - switch Sex stratification to \"All (pooled)\"."))
        expr <- layers$Transcriptomics; meth <- layers$Methylomics
        meta <- d$sample_meta; outcome_col <- input$outcome_col
        covariate_col <- if (nzchar(input$ss_covariate %||% "")) input$ss_covariate else NULL
        engine <- input$ss_engine %||% "diablo"
        p <- ss_params()
        ss_state$error <- NULL; ss_state$result <- NULL; ss_state$submitted <- TRUE
        session$onFlushed(function() ss_task$invoke(expr, meth, meta, outcome_col, covariate_col, sex_mode, engine, p), once = TRUE)
      })
      observe({
        res <- tryCatch(ss_task$result(), error = function(e) e)
        if (inherits(res, "shiny.silent.error")) return()
        ss_state$submitted <- FALSE
        if (inherits(res, "error")) { ss_state$error <- conditionMessage(res); return() }
        if (!isTRUE(res$ok)) { ss_state$error <- res$error; return() } else { ss_state$result <- res }
      })
      ss_running <- reactive(identical(ss_task$status(), "running") || isTRUE(ss_state$submitted))
    } else {
      observeEvent(input$ss_run_btn, {
        sex_mode <- input$ss_sex_mode %||% "pooled"
        validate(need(identical(sex_mode, "pooled") || !is.null(ss_sex_col()), "No Sex/Gender column in this dataset's metadata - switch Sex stratification to \"All (pooled)\"."))
        showNotification("Running Sex-Stratified analysis synchronously - the app will be briefly unresponsive.", type = "message", duration = 5)
        res <- run_ss()
        ss_state$error <- if (!isTRUE(res$ok)) res$error else NULL
        if (isTRUE(res$ok)) ss_state$result <- res
      })
      ss_running <- reactive(FALSE)
    }

    observe({
      if (isTRUE(ss_running())) {
        shinyjs::disable(ns("ss_run_btn"))
        shinyjs::html(ns("ss_run_btn"), as.character(tagList(icon("spinner", class = "fa-spin"), " Running...")))
      } else {
        shinyjs::enable(ns("ss_run_btn"))
        shinyjs::html(ns("ss_run_btn"), as.character(tagList(icon("play"), " Run Sex-Stratified Analysis")))
      }
    })

    output$ss_run_status_ui <- renderUI({
      if (!isTRUE(input$ss_run_btn > 0)) return(NULL)
      if (isTRUE(ss_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running - nested CV can take a minute or more."))
      if (!is.null(ss_state$error)) return(mi_stop(ss_state$error))
      if (!is.null(ss_state$result)) return(mi_ok("Finished - see the results panel to the right."))
      NULL
    })

    output$ss_results_ui <- renderUI({
      if (!isTRUE(input$ss_run_btn > 0)) return(multi_empty_state("Set parameters and click \"Run Sex-Stratified Analysis\" to see results here."))
      if (isTRUE(ss_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running - the rest of the app stays usable while this runs."))
      if (!is.null(ss_state$error)) return(mi_stop(ss_state$error))
      res <- ss_state$result
      if (is.null(res)) return(multi_empty_state())
      perf <- res$performance
      panels <- res$panels
      skipped <- Filter(function(r) !isTRUE(r$ok), res$strata)
      ok_strata <- Filter(function(r) isTRUE(r$ok), res$strata)
      freq_expr <- do.call(rbind, lapply(names(ok_strata), function(nm) {
        f <- mss_selection_frequency(ok_strata[[nm]]$cv$fold_selected_features, "expr")
        if (is.null(f)) return(NULL); cbind(stratum = nm, f)
      }))
      freq_meth <- do.call(rbind, lapply(names(ok_strata), function(nm) {
        f <- mss_selection_frequency(ok_strata[[nm]]$cv$fold_selected_features, "meth")
        if (is.null(f)) return(NULL); cbind(stratum = nm, f)
      }))

      tagList(
        box(width = NULL, title = "Performance by stratum (nested cross-validated AUROC)", status = "primary", solidHeader = FALSE,
            div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px;",
                lapply(seq_len(nrow(perf)), function(i) {
                  call_i <- if ("auroc_call" %in% colnames(perf)) perf$auroc_call[i] else multi_auroc_call(perf$auroc[i], perf$ci_lo[i], perf$ci_hi[i])
                  mark <- switch(call_i, above_chance = " *", below_chance = " †", "")
                  col <- switch(call_i, above_chance = ARTHOMIX_COLORS$aqua, below_chance = ARTHOMIX_COLORS$red, ARTHOMIX_COLORS$blue)
                  mi_stat_card(
                    sprintf("%.3f [%.3f, %.3f]", perf$auroc[i], perf$ci_lo[i], perf$ci_hi[i]),
                    sprintf("%s (n=%d)%s", perf$stratum[i], perf$n[i], mark),
                    col)
                })),
            p(class = "submodule-desc", "* genuinely exceeds chance (95% CI entirely above 0.5). ",
              "† significantly BELOW chance (95% CI entirely below 0.5) - this is NOT a validation success; it typically signals a label/orientation error or overfitting to a spuriously anti-correlated small sample, and should be investigated, not read as evidence of discrimination."),
            DT::dataTableOutput(ns("ss_perf_table")),
            div(class = "table-toolbar", downloadButton(ns("ss_dl_perf"), "Download performance (CSV)", class = "btn-sm"))),
        box(width = NULL, title = if (identical(res$engine, "rf")) "Feature panel (Random Forest importance)" else "Feature panel (DIABLO loadings)", status = "primary", solidHeader = FALSE,
            if (!is.null(res$panel_note)) mi_warn(res$panel_note),
            if (is.null(panels) || nrow(panels) == 0) mi_warn("No panel available.") else tagList(
              DT::dataTableOutput(ns("ss_panel_table")),
              div(class = "table-toolbar", downloadButton(ns("ss_dl_panel"), "Download panel (CSV)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Biomarker comparison by sex", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "One row per feature. A filled Female/Male/Pooled cell means the feature was selected in that stratum's panel; blank means it wasn't."),
            if (is.null(res$panels_wide) || nrow(res$panels_wide) == 0) mi_warn("Not available - run \"Female and Male separately\" (or compare separate pooled/female/male runs) to populate this comparison.") else tagList(
              DT::dataTableOutput(ns("ss_panel_wide_table")),
              div(class = "table-toolbar", downloadButton(ns("ss_dl_panel_wide"), "Download comparison (CSV)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Feature-selection stability (frequency across CV folds)", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "How often each feature was selected across the folds x repeats run for that stratum (same concept as the Stability tab)."),
            h5("Expression"),
            if (is.null(freq_expr) || nrow(freq_expr) == 0) mi_warn("Not available.") else DT::dataTableOutput(ns("ss_stability_expr_table")),
            h5("Methylation"),
            if (is.null(freq_meth) || nrow(freq_meth) == 0) mi_warn("Not available.") else DT::dataTableOutput(ns("ss_stability_meth_table"))),
        if (length(skipped) > 0) box(width = NULL, title = "Strata skipped", status = "primary", solidHeader = FALSE,
            tags$ul(lapply(names(skipped), function(nm) tags$li(sprintf("%s: %s", nm, skipped[[nm]]$error)))))
      )
    })

    output$ss_perf_table <- DT::renderDataTable({
      DT::datatable(req(ss_state$result)$performance, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$ss_dl_perf <- downloadHandler(function() "biomarker_sex_stratified_performance.csv", function(file) utils::write.csv(req(ss_state$result)$performance, file, row.names = FALSE))
    output$ss_panel_table <- DT::renderDataTable({
      DT::datatable(req(ss_state$result)$panels, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$ss_dl_panel <- downloadHandler(function() "biomarker_sex_stratified_panel.csv", function(file) utils::write.csv(req(ss_state$result)$panels, file, row.names = FALSE))
    output$ss_panel_wide_table <- DT::renderDataTable({
      DT::datatable(req(ss_state$result)$panels_wide, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$ss_dl_panel_wide <- downloadHandler(function() "biomarker_sex_stratified_comparison.csv", function(file) utils::write.csv(req(ss_state$result)$panels_wide, file, row.names = FALSE))
    output$ss_stability_expr_table <- DT::renderDataTable({
      res <- req(ss_state$result)
      ok_strata <- Filter(function(r) isTRUE(r$ok), res$strata)
      df <- do.call(rbind, lapply(names(ok_strata), function(nm) { f <- mss_selection_frequency(ok_strata[[nm]]$cv$fold_selected_features, "expr"); if (is.null(f)) NULL else cbind(stratum = nm, f) }))
      DT::datatable(req(df), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    output$ss_stability_meth_table <- DT::renderDataTable({
      res <- req(ss_state$result)
      ok_strata <- Filter(function(r) isTRUE(r$ok), res$strata)
      df <- do.call(rbind, lapply(names(ok_strata), function(nm) { f <- mss_selection_frequency(ok_strata[[nm]]$cv$fold_selected_features, "meth"); if (is.null(f)) NULL else cbind(stratum = nm, f) }))
      DT::datatable(req(df), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    observe({
      if (is.null(multi_results) || is.null(mb_state$result)) return()
      multi_results$biomarker <- list(df = sig_df(), cohort = mb_state$dataset_label)
    })
    observe({
      if (is.null(multi_results)) return()
      multi_results$biomarker_stratified <- if (!is.null(ss_state$result)) list(cohort = mb_state$dataset_label, result = ss_state$result) else NULL
    })

    ## ---- provenance records (session-wide Analysis records log) ----
    observeEvent(mb_state$result, {
      res <- mb_state$result; if (is.null(res) || !isTRUE(res$ok)) return()
      r <- res$diablo; ps <- mi_diablo_performance_summary(r)
      arthomix_provenance_push(arthomix_provenance_record(
        module = "mod_multi_biomarker",
        checksum_input = list(train_ids = res$train_ids, test_ids = res$test_ids, features = lapply(r$fit$X, colnames), classes = r$params$classes),
        params = list(n_train = length(res$train_ids), n_holdout = length(res$test_ids), holdout_frac = res$holdout_frac, outcome = input$outcome_col,
                      classes = r$params$classes, ncomp = r$params$ncomp, ncomp_mode = r$params$ncomp_mode, keepX = r$params$keepX, keepx_mode = r$params$keepx_mode,
                      design = r$params$design["Transcriptomics", "Methylomics"], distance = r$params$distance, validation = r$params$validation_method,
                      folds = r$params$folds, nrepeat = r$params$nrepeat, nested = isTRUE(r$nested), prefilter_t = input$max_features_t, prefilter_m = input$max_features_m,
                      ber = ps$ber, oof_auroc = if (!is.null(ps$auc)) ps$auc$AUC[1] else NA_real_,
                      holdout_auroc = if (!is.null(res$holdout) && isTRUE(res$holdout$ok)) res$holdout$auc else NA_real_,
                      preselected_input = isTRUE(mb_state$preselected_input)),
        seed = input$seed %||% 1, packages = c("mixOmics", "pROC", "caret"),
        extra = list(dataset = mb_state$dataset_label)
      ), session = session, dedupe = TRUE)
    })
    observeEvent(ext_state$result, {
      ev <- ext_state$result; if (is.null(ev)) return()
      arthomix_provenance_push(arthomix_provenance_record(
        module = "mod_multi_biomarker_external_validation",
        checksum_input = list(scores = ev$scores, coverage = ev$coverage, error = ev$error),
        params = list(ok = isTRUE(ev$ok), external = ev$external_label, n_external = ev$n, n_train = ev$n_train, coverage = ev$overall_coverage,
                      min_coverage = input$ext_min_coverage, outcome_col = input$ext_outcome_col, pos_levels = input$ext_pos_levels, neg_levels = input$ext_neg_levels,
                      auroc = ev$auc, ci_lo = ev$ci_lo, ci_hi = ev$ci_hi, ber = ev$ber, error = ev$error),
        seed = NULL, packages = c("mixOmics", "pROC", "org.Hs.eg.db"),
        extra = list(dataset = mb_state$dataset_label)
      ), session = session, dedupe = TRUE)
    })
    observeEvent(ss_state$result, {
      r <- ss_state$result; if (is.null(r) || !isTRUE(r$ok)) return()
      arthomix_provenance_push(arthomix_provenance_record(
        module = "mod_multi_biomarker_sexstratified",
        checksum_input = list(perf = r$performance %||% r$per_stratum %||% names(r)),
        params = list(engine = input$ss_engine, sex_mode = input$ss_sex_mode, covariate = input$ss_covariate, top_expr = input$ss_top_expr, top_meth = input$ss_top_meth,
                      ncomp = input$ss_ncomp, folds = input$ss_folds, repeats = input$ss_repeats),
        seed = NULL, packages = c("mixOmics", "randomForest", "limma", "pROC"),
        extra = list(dataset = mb_state$dataset_label)
      ), session = session, dedupe = TRUE)
    })
  })
}
