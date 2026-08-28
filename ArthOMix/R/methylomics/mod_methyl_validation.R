## R/methylomics/mod_methyl_validation.R
## Applies the Diagnostic Classifier's already-trained models (results$diagnostic_models) to an
## independent EXTERNAL cohort - never retrains. Mirrors mod_methyl_diagnostic.R's per-model tab
## layout (one tab per algorithm, Model Comparison, Test External Data, Export) but every "test"
## evaluation here is scored on the external validation cohort loaded on the Datasets tab, not the
## internal train/test split. Reuses mod_methyl_diagnostic.R's dxm_* helpers and DXM_MODEL_SPECS registry.

mod_methyl_validation_config <- list(
  id = "validation", title = "Validation", icon = "flask-vial", group = "Biomarker modeling",
  description = "External validation of the Diagnostic Classifier's trained models on an independent cohort"
)

## vld_* helpers - new to this module only

## Required-vs-available CpG overlap; `ok` gates every downstream step since
## predict() requires the exact trained column set (100% match).
vld_feature_alignment <- function(required_ids, available_ids) {
  matched <- intersect(required_ids, available_ids)
  missing <- setdiff(required_ids, available_ids)
  extra <- setdiff(available_ids, required_ids)
  n_req <- length(required_ids)
  list(required = n_req, available = length(available_ids), matched = length(matched),
       missing = missing, n_missing = length(missing), extra = extra, n_extra = length(extra),
       overlap_pct = if (n_req > 0) 100 * length(matched) / n_req else NA_real_,
       ok = n_req > 0 && length(missing) == 0, matched_ids = matched)
}

## Duplicate-sample-ID check when real IDs exist; otherwise reports
## independence as unverified rather than assuming it.
vld_sample_overlap <- function(train_ids, val_ids) {
  if (length(train_ids) == 0 || length(val_ids) == 0) {
    return(list(status = "unknown", message = "Sample independence could not be verified.", n_overlap = NA_integer_))
  }
  ov <- intersect(train_ids, val_ids)
  if (length(ov) > 0) {
    list(status = "overlap", n_overlap = length(ov), ids = ov,
         message = sprintf("%d sample ID(s) appear in both the training cohort and this validation cohort - independence is NOT confirmed.", length(ov)))
  } else {
    list(status = "independent", n_overlap = 0,
         message = "No shared sample IDs detected between the training cohort and this validation cohort.")
  }
}

## Component / Training Pipeline / Validation Requirement / Status audit table.
vld_compat_table <- function(model, cohort, align) {
  row <- function(component, training, requirement, status) {
    data.frame(Component = component, `Training Pipeline` = training,
               `Validation Requirement / Applied` = requirement, Status = status,
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  do.call(rbind, list(
    row("Data representation",
        "M-value (logit of beta): log2(b / (1-b)), fixed and elementwise",
        sprintf("Uploaded as %s; converted with the identical fixed transform before scoring",
                if (identical(cohort$scale_declared, "m")) "M-value (no conversion needed)" else "beta-value"),
        "Compatible"),
    row("Feature / CpG set",
        sprintf("%d CpG(s) fixed by the trained model (chosen upstream via WGCNA / Feature Selection / manual pick, never selected in this module)", length(model$feature_ids)),
        sprintf("Same %d CpG ID(s), same order, required exactly - %d matched, %d missing (%.1f%% overlap)",
                align$required, align$matched, align$n_missing, align$overlap_pct),
        if (align$ok) "Pass" else "Fail"),
    row("Scaling (center/scale; SVM/kNN only)",
        "Fit once, inside caret::train(), on the training partition only",
        "Reused unchanged from the stored fitted model object via predict() - never refit on validation data",
        "Pass"),
    row("Normalization / batch correction",
        "Noob normalization + granulocyte cell-type adjustment, performed upstream of this app (script09 provenance)",
        cohort$normalization_note, cohort$normalization_status),
    row("Missing values",
        "Detected/reported only, never imputed (dxm_validate_checklist)",
        sprintf("%d missing value(s) at required CpGs in this validation cohort - samples with any missing required CpG are excluded from scoring, never imputed", cohort$n_missing),
        if (cohort$n_missing == 0) "Pass" else "Warning"),
    row("Classification threshold",
        sprintf("%.3f, chosen from the training ROC only via Youden's J / a fixed sensitivity-specificity target / 0.5 (never from test or validation data)", model$threshold),
        "Reused unchanged - never re-optimized on validation labels",
        "Pass"),
    row("Platform / array type",
        "Illumina Infinium methylation array (per upstream provenance)",
        cohort$platform, cohort$platform_status)
  ))
}

## Bootstrap CI for sensitivity/specificity/accuracy (AUC already gets a DeLong CI).
vld_bootstrap_ci <- function(y, prob, threshold, stat_fn, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(y)
  vals <- vapply(seq_len(n_boot), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(y[idx])) < 2) return(NA_real_)
    tryCatch(stat_fn(y[idx], prob[idx], threshold), error = function(e) NA_real_)
  }, numeric(1))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 50) return(c(NA_real_, NA_real_))
  stats::quantile(vals, c(0.025, 0.975), names = FALSE)
}

vld_ci_label <- function(ci) if (is.null(ci) || any(is.na(ci))) "NA" else sprintf("%.3f-%.3f", ci[1], ci[2])

## which.max() coerces character input to numeric (silently producing all-NA and an empty
## index), so it can't be used on "%Y-%m-%d %H:%M:%S" ran_at timestamps - use ordinary
## character comparison instead, which sorts those timestamps correctly.
vld_which_latest <- function(timestamps) which(timestamps == max(timestamps))[1]

## validate()/need() only display when caught by a render/output context - inside an
## observeEvent (button click) there's no output to catch it, so a failed check is otherwise
## completely silent. Use `if (vld_notify_fail(cond, msg)) return()` in observers instead.
vld_notify_fail <- function(cond, msg) {
  if (isTRUE(cond)) return(FALSE)
  showNotification(msg, type = "error", duration = 10)
  TRUE
}

## Scores one trained model (a results$diagnostic_models entry) against the loaded external
## cohort, exactly as trained: no refitting, re-tuning, or threshold re-optimization. Returns
## invisible(FALSE) (after a user-facing notification) on any incompatibility or data problem.
vld_do_run_external <- function(mid, spec, m, cohort, vms, vruns) {
  if (vld_notify_fail(isTRUE(cohort$loaded), "Load an external validation cohort first (Datasets tab).")) return(invisible(FALSE))
  if (vld_notify_fail(identical(m$ref_level, cohort$ref_lab) && identical(m$comp_level, cohort$comp_lab),
        sprintf("This model was trained on different class labels (%s/%s) than the loaded cohort (%s/%s) - reload the cohort or retrain with matching labels.",
                m$ref_level, m$comp_level, cohort$ref_lab, cohort$comp_lab))) return(invisible(FALSE))
  if (vld_notify_fail(!(identical(cohort$source, "preloaded") && identical(m$sex_stratum, "male")),
        "The preloaded external cohort is all-female; cannot validate a male-stratum model. Upload an independent cohort with male samples instead.")) return(invisible(FALSE))
  a <- vld_feature_alignment(m$feature_ids, rownames(cohort$m))
  if (vld_notify_fail(a$ok, sprintf("%d of %d required CpG(s) are missing from this validation cohort. See Compatibility.", a$n_missing, a$required))) return(invisible(FALSE))

  Xv <- as.data.frame(t(cohort$m[m$feature_ids, , drop = FALSE]))
  yv <- cohort$y
  complete <- stats::complete.cases(Xv)
  if (!all(complete)) {
    showNotification(sprintf("%d of %d validation sample(s) excluded: missing value(s) at a required CpG (never imputed).", sum(!complete), nrow(Xv)),
                      type = "warning", duration = 10)
  }
  Xv <- Xv[complete, , drop = FALSE]; yv <- yv[complete]
  if (vld_notify_fail(nrow(Xv) >= 5, "Fewer than 5 validation samples with complete data at all required CpGs - cannot compute reliable performance metrics.")) return(invisible(FALSE))
  if (vld_notify_fail(length(unique(yv)) == 2, "The validation cohort (after filtering) contains only one class - sensitivity/specificity/AUC require both classes.")) return(invisible(FALSE))

  prob <- dxm_predict_prob(m$fit$model, Xv)
  roc <- dxm_roc_bundle(yv, prob)
  if (vld_notify_fail(!is.null(roc), "Could not compute an ROC curve for this validation cohort.")) return(invisible(FALSE))
  met <- dxm_metrics_bundle(yv, prob, m$threshold, roc)
  conf <- dxm_confusion(yv, prob, m$threshold)
  seed <- m$seed %||% 42
  ci_sens <- vld_bootstrap_ci(yv, prob, m$threshold, function(y, p, th) dxm_confusion(y, p, th)$sensitivity, seed = seed)
  ci_spec <- vld_bootstrap_ci(yv, prob, m$threshold, function(y, p, th) dxm_confusion(y, p, th)$specificity, seed = seed + 1)
  ci_acc  <- vld_bootstrap_ci(yv, prob, m$threshold, function(y, p, th) dxm_confusion(y, p, th)$accuracy, seed = seed + 2)

  vms$ext_prob <- prob; vms$ext_roc <- roc; vms$ext_metrics <- met; vms$ext_confusion <- conf; vms$ext_y <- yv
  vms$ci_sens <- ci_sens; vms$ci_spec <- ci_spec; vms$ci_acc <- ci_acc
  vms$n_used <- nrow(Xv); vms$n_dropped <- sum(!complete)
  vms$ran <- TRUE; vms$ran_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  vms$roc_generated <- FALSE; vms$calib_generated <- FALSE; vms$ext_calib <- NULL

  key <- paste(mid, m$analysis_type, paste(m$feature_ids, collapse = "|"))
  vruns[[key]] <- list(model_id = mid, label = spec$label, analysis_type = m$analysis_type, feature_ids = m$feature_ids,
                        threshold = m$threshold, train_metrics = m$train_metrics, train_roc = m$train_roc, cv_roc = m$cv_roc,
                        ext_metrics = met, ext_roc = roc, n_used = nrow(Xv), n_dropped = sum(!complete), ran_at = vms$ran_at)

  showNotification(sprintf("%s: external validation complete (n=%d, external AUC = %.3f).", spec$label, nrow(Xv), roc$auc), type = "message")
  invisible(TRUE)
}

## Per-model tab: same generic pattern as mod_methyl_diagnostic.R's dxm_render_model_panel(), but
## "test" is always the loaded external cohort, and "Results: Training" reuses the metrics already
## published by Diagnostic Classifier rather than recomputing them.
vld_render_model_panel <- function(mid, spec, ns, entries, cohort, vms, m) {
  if (length(entries) == 0) {
    return(p(class = "text-muted", sprintf(
      "No trained %s model available yet. Go to Diagnostic Classifier, train %s on any feature panel, click \"Run Test Evaluation\", then return here.",
      spec$label, spec$label)))
  }

  setup_box <- box(width = 12, status = "primary", solidHeader = TRUE, title = "Setup",
    if (length(entries) > 1) selectInput(ns(paste0(mid, "_run_select")), "Trained run",
      choices = stats::setNames(names(entries), vapply(entries, function(e) sprintf("%s (%d feat., tested %s)", e$analysis_type, length(e$feature_ids), e$ran_at %||% "unknown"), character(1)))),
    tags$ul(
      tags$li(sprintf("Analysis type: %s (%d feature%s)", m$analysis_type, length(m$feature_ids), if (length(m$feature_ids) == 1) "" else "s")),
      tags$li(sprintf("Phenotype: %s (reference) vs %s (comparison/positive class)", m$ref_level, m$comp_level)),
      tags$li(sprintf("Training stratum / mode: %s / %s", tools::toTitleCase(as.character(m$sex_stratum %||% "NA")), m$mode)),
      tags$li(sprintf("Trained on %d sample(s); published by Diagnostic Classifier %s", m$train_n, m$ran_at %||% "unknown")),
      tags$li(sprintf("External validation cohort: %s", if (isTRUE(cohort$loaded)) cohort$label else "not loaded yet (see Datasets tab)")),
      tags$li(sprintf("External test: %s", if (isTRUE(vms$ran)) "evaluated" else "not yet evaluated"))
    ))

  train_box <- box(width = 12, status = "success", solidHeader = TRUE, title = "Results: Training",
    fluidRow(
      valueBox(sprintf("%.3f", m$train_metrics$auc), "Training AUC", icon = icon("chart-line"), color = "blue", width = 3),
      valueBox(sprintf("%.3f", m$cv_roc$mean_auc), sprintf("Mean CV AUC (+/- %.3f SD, %d folds)", m$cv_roc$sd_auc, m$cv_roc$n_folds), icon = icon("layer-group"), color = "blue", width = 3),
      valueBox(sprintf("%.3f", m$threshold), "Classification threshold", icon = icon("ruler"), color = "light-blue", width = 3),
      valueBox(length(m$feature_ids), "Feature(s) used", icon = icon("dna"), color = "light-blue", width = 3)
    ),
    DT::dataTableOutput(ns(paste0(mid, "_train_metrics_table"))),
    br(), actionButton(ns(paste0(mid, "_run_btn")), "Run External Validation", icon = icon("play"), class = "btn-primary"))

  out <- list(setup_box, train_box)

  if (isTRUE(vms$ran)) {
    out <- c(out, list(box(width = 12, status = "success", solidHeader = TRUE, title = "Results: External Test",
      p(class = "submodule-desc", sprintf("Scored on %s (n = %d after filtering).", cohort$label %||% "the loaded external cohort", vms$n_used)),
      fluidRow(
        valueBox(sprintf("%.3f", vms$ext_metrics$auc), sprintf("External AUC (95%% CI %.3f-%.3f)", vms$ext_metrics$auc_ci_lo, vms$ext_metrics$auc_ci_hi), icon = icon("chart-area"), color = "blue", width = 3),
        valueBox(sprintf("%.3f", vms$ext_metrics$sensitivity), sprintf("Sensitivity (95%% CI %s)", vld_ci_label(vms$ci_sens)), icon = icon("check"), color = "light-blue", width = 3),
        valueBox(sprintf("%.3f", vms$ext_metrics$specificity), sprintf("Specificity (95%% CI %s)", vld_ci_label(vms$ci_spec)), icon = icon("shield"), color = "light-blue", width = 3),
        valueBox(vms$n_used, sprintf("Validated N (%d ref, %d comp)", sum(vms$ext_y == DXM_NEG), sum(vms$ext_y == DXM_POS)), icon = icon("users"), color = "teal", width = 3)
      ),
      DT::dataTableOutput(ns(paste0(mid, "_ext_metrics_table"))),
      hr(), p(class = "submodule-desc", dxm_overfitting_note(m$train_metrics$auc, m$cv_roc$mean_auc, vms$ext_metrics$auc)))))

    out <- c(out, list(box(width = 12, status = "primary", solidHeader = TRUE, title = "ROC / AUC",
      actionButton(ns(paste0(mid, "_roc_btn")), "Generate ROC/AUC", icon = icon("chart-area"), class = "btn-primary"),
      downloadButton(ns(paste0(mid, "_roc_png")), "Download ROC plot (PNG)", class = "btn-sm"),
      br(), br(),
      if (isTRUE(vms$roc_generated)) fluidRow(
        column(6, plotOutput(ns(paste0(mid, "_roc_train_plot")))),
        column(6, plotOutput(ns(paste0(mid, "_roc_cv_plot"))))
      ),
      if (isTRUE(vms$roc_generated)) fluidRow(
        column(6, plotOutput(ns(paste0(mid, "_roc_ext_plot"))))
      ))))

    out <- c(out, list(box(width = 12, status = "primary", solidHeader = TRUE, title = "Diagnostics",
      fluidRow(
        column(6, h5("Training Confusion Matrix"), DT::dataTableOutput(ns(paste0(mid, "_confusion_train_table")))),
        column(6, h5("External Test Confusion Matrix"), DT::dataTableOutput(ns(paste0(mid, "_confusion_ext_table"))))
      ),
      hr(),
      actionButton(ns(paste0(mid, "_calib_btn")), "Generate Calibration", icon = icon("chart-line"), class = "btn-primary btn-sm"),
      if (isTRUE(vms$calib_generated)) plotOutput(ns(paste0(mid, "_calib_plot"))))))
  }
  tagList(out)
}

## Registers one model tab's server logic - mirrors dxm_register_model_server()'s shape.
vld_register_model_server <- function(mid, spec, input, output, session, ns, avail_models, cohort, vms, vruns) {

  entries_r <- reactive({ Filter(function(m) identical(m$model_id, mid), avail_models()) })

  selected_key <- reactive({
    entries <- entries_r(); if (length(entries) == 0) return(NULL)
    sel <- input[[paste0(mid, "_run_select")]]
    if (!is.null(sel) && sel %in% names(entries)) return(sel)
    names(entries)[vld_which_latest(vapply(entries, function(m) m$ran_at %||% "", character(1)))]
  })

  observeEvent(input[[paste0(mid, "_run_btn")]], {
    key <- selected_key(); req(key)
    m <- entries_r()[[key]]
    vld_do_run_external(mid, spec, m, cohort, vms, vruns)
  })

  observeEvent(input[[paste0(mid, "_roc_btn")]], { req(vms$ran); vms$roc_generated <- TRUE })

  observeEvent(input[[paste0(mid, "_calib_btn")]], {
    req(vms$ran)
    vms$ext_calib <- dxm_calibration(vms$ext_y, vms$ext_prob)
    vms$calib_generated <- TRUE
  })

  output[[paste0(mid, "_train_metrics_table")]] <- DT::renderDataTable({
    key <- selected_key(); req(key)
    DT::datatable(dxm_metrics_display(entries_r()[[key]]$train_metrics), rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  output[[paste0(mid, "_ext_metrics_table")]] <- DT::renderDataTable({
    req(vms$ran); DT::datatable(dxm_metrics_display(vms$ext_metrics), rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  output[[paste0(mid, "_confusion_train_table")]] <- DT::renderDataTable({
    key <- selected_key(); req(key); m <- entries_r()[[key]]; req(m$confusion_train)
    DT::datatable(as.data.frame.matrix(m$confusion_train$table), options = list(dom = "t", paging = FALSE))
  })
  output[[paste0(mid, "_confusion_ext_table")]] <- DT::renderDataTable({
    req(vms$ran); DT::datatable(as.data.frame.matrix(vms$ext_confusion$table), options = list(dom = "t", paging = FALSE))
  })

  output[[paste0(mid, "_roc_train_plot")]] <- renderPlot({
    key <- selected_key(); req(key, vms$roc_generated); m <- entries_r()[[key]]; req(m$train_roc)
    p <- dxm_plot_roc(list(Training = m$train_roc), sprintf("%s - Training ROC", spec$label)); vms$last_roc_plot <- p; p
  })
  output[[paste0(mid, "_roc_cv_plot")]] <- renderPlot({
    key <- selected_key(); req(key, vms$roc_generated); m <- entries_r()[[key]]; req(m$cv_roc)
    dxm_plot_cv_roc(m$cv_roc, sprintf("%s - Cross-Validated ROC", spec$label))
  })
  output[[paste0(mid, "_roc_ext_plot")]] <- renderPlot({
    req(vms$roc_generated, vms$ext_roc); dxm_plot_roc(list(`External Test` = vms$ext_roc), sprintf("%s - External Test ROC", spec$label))
  })
  output[[paste0(mid, "_calib_plot")]] <- renderPlot({
    req(vms$calib_generated, vms$ext_calib); p <- dxm_plot_calibration(vms$ext_calib, sprintf("%s - External Test Calibration", spec$label)); vms$last_calib_plot <- p; p
  })

  output[[paste0(mid, "_roc_png")]] <- downloadHandler(
    filename = function() sprintf("methylomics_validation_%s_roc.png", mid),
    content = function(file) { req(vms$last_roc_plot); ggplot2::ggsave(file, vms$last_roc_plot, width = 7, height = 6, dpi = 150) }
  )

  output[[paste0("panel_", mid)]] <- renderUI({
    key <- selected_key(); entries <- entries_r()
    vld_render_model_panel(mid, spec, ns, entries, cohort, vms, if (!is.null(key)) entries[[key]] else NULL)
  })

  invisible(NULL)
}

## UI

mod_methyl_validation_ui <- function(id) {
  ns <- NS(id)
  model_tabs <- lapply(DXM_MODEL_SPECS, function(spec) {
    tabPanel(spec$label, br(), withSpinner(uiOutput(ns(paste0("panel_", spec$id))), color = "#2563EB", type = 6))
  })

  do.call(tabsetPanel, c(
    list(id = ns("main_tabs"), type = "tabs",

      tabPanel("Datasets", br(),
        p(class = "submodule-desc",
          "Identifies every trained model in the Diagnostic Classifier and loads external validation data."),
        fluidRow(
          column(4,
            box(width = NULL, title = "Trained models (from Diagnostic Classifier)", status = "primary", solidHeader = TRUE,
              uiOutput(ns("model_summary_ui"))),
            box(width = NULL, title = "Validation cohort source", status = "primary", solidHeader = FALSE,
              radioButtons(ns("cohort_source"), NULL,
                choices = c("Preloaded external validation cohort" = "preloaded", "Upload an independent cohort" = "upload"),
                selected = "preloaded"),
              conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("cohort_source")),
                p(class = "submodule-desc", "43 samples, all female, Noob-renormalized and granulocyte-adjusted."),
                radioButtons(ns("cohort_sex_stratum"), "Sex stratum", choices = c("All samples" = "all", "Female" = "female", "Male" = "male"), selected = "all", inline = TRUE),
                actionButton(ns("load_preloaded_btn"), "Load Preloaded External Cohort", icon = icon("download"), class = "btn-primary")),
              conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("cohort_source")),
                fileInput(ns("val_upload_matrix"), "Methylation matrix (CSV/TSV, probes x samples)"),
                radioButtons(ns("val_upload_scale"), "Matrix scale", choices = c("Beta-value" = "beta", "M-value" = "m"), selected = "beta", inline = TRUE),
                fileInput(ns("val_upload_sheet"), "Sample sheet (CSV/TSV)"),
                uiOutput(ns("val_upload_col_ui")),
                radioButtons(ns("val_upload_sex_stratum"), "Sex stratum", choices = c("All samples" = "all", "Female" = "female", "Male" = "male"), selected = "all", inline = TRUE),
                uiOutput(ns("val_upload_sex_col_ui")),
                actionButton(ns("load_upload_btn"), "Load Uploaded Cohort", icon = icon("upload"), class = "btn-primary"))
            )
          ),
          column(8,
            uiOutput(ns("cohort_compare_gate_ui")),
            DT::dataTableOutput(ns("cohort_compare_table")),
            uiOutput(ns("cohort_overlap_ui")),
            hr(),
            h5("Validation cohort history (this session)"),
            DT::dataTableOutput(ns("cohort_history_table"))
          )
        )
      ),

      tabPanel("Compatibility", br(),
        p(class = "submodule-desc", "Checks each model's required CpGs are present in the external cohort."),
        uiOutput(ns("compat_gate_ui")),
        h5("All trained models"),
        DT::dataTableOutput(ns("compat_overview_table")),
        hr(),
        h5("Detailed compatibility audit"),
        selectInput(ns("compat_model_select"), "Model", choices = NULL, width = "100%"),
        uiOutput(ns("compat_kpi_ui")),
        DT::dataTableOutput(ns("compat_table"))
      )
    ),
    unname(model_tabs),
    list(
      tabPanel("Model Comparison", br(), withSpinner(uiOutput(ns("compare_ui")), color = "#2563EB", type = 6)),
      tabPanel("Test External Data", br(), withSpinner(uiOutput(ns("testdata_ui")), color = "#2563EB", type = 6)),
      tabPanel("Export", br(), withSpinner(uiOutput(ns("export_ui")), color = "#2563EB", type = 6))
    )
  ))
}

## Server

mod_methyl_validation_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cohort <- reactiveValues(
      loaded = FALSE, source = NULL, label = NULL, m = NULL, y = NULL,
      sample_ids = NULL, had_id_col = FALSE, n_missing = 0,
      normalization_note = NULL, normalization_status = NULL,
      platform = NULL, platform_status = NULL, accession = NULL,
      scale_declared = NULL, n_samples = NULL, ref_lab = NULL, comp_lab = NULL
    )
    vhist <- reactiveValues(rows = list())
    vruns <- reactiveValues()
    compare_state <- reactiveValues(generated = FALSE, bundles = NULL)

    ## All models published by the Diagnostic Classifier (one per model/feature-panel
    ## combination that was run and test-evaluated there).
    avail_models <- reactive({ results$diagnostic_models %||% list() })

    ## Reference model used only to seed the validation cohort's class labels at load
    ## time (the most recently tested model, across all algorithms) - each model is still
    ## checked individually for label/feature compatibility before being scored.
    ref_model <- reactive({
      models <- avail_models(); req(length(models) > 0)
      models[[vld_which_latest(vapply(models, function(m) m$ran_at %||% "", character(1)))]]
    })

    observe({
      models <- avail_models()
      ch <- if (length(models) == 0) character(0) else stats::setNames(names(models),
        vapply(models, function(m) sprintf("%s - %s (%d feat., tested %s)", m$label, m$analysis_type, length(m$feature_ids), m$ran_at %||% "unknown"), character(1)))
      updateSelectInput(session, "compat_model_select", choices = ch)
    })

    ## ---- Cohort --------------------------------------------------------

    output$model_summary_ui <- renderUI({
      models <- avail_models()
      if (length(models) == 0) {
        return(p(class = "text-danger", icon("triangle-exclamation"),
          " No trained model available yet. Go to Diagnostic Classifier, train a model on any tab, click \"Run Test Evaluation\", then return here. Every model you test-evaluate there becomes available here."))
      }
      tagList(
        p(class = "submodule-desc", sprintf("%d trained model(s) available.", length(models))),
        DT::dataTableOutput(ns("model_list_table"))
      )
    })

    output$model_list_table <- DT::renderDataTable({
      models <- avail_models(); req(length(models) > 0)
      df <- do.call(rbind, lapply(models, function(m) data.frame(
        Algorithm = m$label, `Analysis type` = m$analysis_type, Features = length(m$feature_ids),
        Contrast = sprintf("%s vs %s", m$ref_level, m$comp_level), Stratum = tools::toTitleCase(as.character(m$sex_stratum %||% "NA")),
        `Train AUC` = round(m$train_metrics$auc %||% NA_real_, 3), `CV AUC` = round(m$cv_roc$mean_auc %||% NA_real_, 3),
        `Internal-test AUC` = round(m$test_internal_metrics$auc %||% NA_real_, 3), `Tested at` = m$ran_at %||% "unknown",
        check.names = FALSE, stringsAsFactors = FALSE)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    output$val_upload_col_ui <- renderUI({
      req(input$val_upload_sheet)
      ps <- methyl_parse_sample_sheet(input$val_upload_sheet$datapath, input$val_upload_sheet$name)
      if (!isTRUE(ps$ok)) return(p(class = "text-danger", ps$error))
      cols <- colnames(ps$df)
      guess <- cols[which(grepl("class|phenotype|group|status|disease", cols, ignore.case = TRUE))[1]]
      selectInput(ns("val_upload_pheno_col"), "Phenotype / class column", choices = cols, selected = guess %||% cols[1])
    })

    output$val_upload_sex_col_ui <- renderUI({
      req(input$val_upload_sheet)
      ps <- methyl_parse_sample_sheet(input$val_upload_sheet$datapath, input$val_upload_sheet$name)
      if (!isTRUE(ps$ok)) return(NULL)
      cols <- colnames(ps$df)
      guess <- cols[which(grepl("sex|gender", cols, ignore.case = TRUE))[1]]
      selectInput(ns("val_upload_sex_col"), "Sex column (used when Sex stratum is Female/Male)", choices = cols, selected = guess %||% cols[1])
    })

    observeEvent(input$load_preloaded_btn, {
      models <- avail_models()
      if (vld_notify_fail(length(models) > 0, "No trained model available yet - train and test-evaluate at least one model in Diagnostic Classifier first.")) return()
      if (vld_notify_fail(METH_DIAG_DATA_AVAILABLE, "The preloaded external cohort isn't available in this deployment.")) return()
      rm <- ref_model()
      dd <- load_default_diagnostic_train_test()
      if (vld_notify_fail(!is.null(dd), "Could not load the preloaded external cohort.")) return()
      ext <- dd$external
      ref_lab <- rm$ref_level; comp_lab <- rm$comp_level
      sex_sel <- input$cohort_sex_stratum %||% "all"
      sex_code <- switch(sex_sel, male = "M", female = "F", NA_character_)
      sex_keep <- if (is.na(sex_code)) rep(TRUE, nrow(ext$pheno)) else ext$pheno$sex == sex_code
      keep <- sex_keep & ext$pheno$group %in% c(ref_lab, comp_lab)
      if (vld_notify_fail(sum(keep) > 5, sprintf(
        "Fewer than 5 preloaded external cohort samples match sex = %s and the most recently trained model's class labels (%s / %s).%s",
        dxm_sex_label(sex_sel), ref_lab, comp_lab,
        if (identical(sex_sel, "male")) " This published external cohort is all-female (43/43) - select \"All samples\" or \"Female\", or upload an independent cohort with male samples." else ""))) return()
      beta_sub <- ext$beta[, keep, drop = FALSE]
      pheno_sub <- ext$pheno[keep, , drop = FALSE]
      Mv <- dxm_beta_to_m(beta_sub)
      colnames(Mv) <- pheno_sub$gsm
      y <- factor(ifelse(pheno_sub$group == comp_lab, DXM_POS, DXM_NEG), levels = c(DXM_NEG, DXM_POS))

      cohort$loaded <- TRUE; cohort$source <- "preloaded"; cohort$label <- sprintf("External validation cohort (preloaded, %s)", dxm_sex_label(sex_sel))
      cohort$m <- Mv; cohort$y <- y; cohort$sample_ids <- colnames(Mv); cohort$had_id_col <- TRUE
      cohort$n_missing <- sum(is.na(Mv))
      cohort$normalization_note <- "Also Noob-renormalized and granulocyte cell-type-adjusted specifically to match the training cohort's own preprocessing (script09 provenance) - same processing family as training"
      cohort$normalization_status <- "Pass"
      cohort$platform <- "Illumina 450K, same upstream reprocessing pipeline as the training cohort"
      cohort$platform_status <- "Compatible"
      cohort$accession <- "Preloaded external cohort"
      cohort$scale_declared <- "beta"
      cohort$n_samples <- ncol(Mv)
      cohort$ref_lab <- ref_lab; cohort$comp_lab <- comp_lab

      vhist$rows[[length(vhist$rows) + 1]] <- list(label = cohort$label, n = ncol(Mv), source = "Preloaded (external cohort)", loaded_at = format(Sys.time(), "%H:%M:%S"))
      showNotification(sprintf("Loaded external validation cohort: %d samples (%s=%d, %s=%d).", ncol(Mv), ref_lab, sum(y == DXM_NEG), comp_lab, sum(y == DXM_POS)), type = "message", duration = 10)
    })

    observeEvent(input$load_upload_btn, {
      models <- avail_models()
      if (vld_notify_fail(length(models) > 0, "No trained model available yet - train and test-evaluate at least one model in Diagnostic Classifier first.")) return()
      rm <- ref_model()
      req(input$val_upload_matrix, input$val_upload_sheet)
      pm <- methyl_parse_matrix(input$val_upload_matrix$datapath, input$val_upload_matrix$name)
      if (vld_notify_fail(isTRUE(pm$ok), pm$error %||% "Could not parse the uploaded methylation matrix.")) return()
      ps <- methyl_parse_sample_sheet(input$val_upload_sheet$datapath, input$val_upload_sheet$name)
      if (vld_notify_fail(isTRUE(ps$ok), ps$error %||% "Could not parse the uploaded sample sheet.")) return()
      mat <- pm$mat; sheet <- ps$df
      had_id_col <- length(intersect(c("sample", "Sample", "sample_id", "Sample_ID"), colnames(sheet))) > 0
      sample_ids <- methyl_sheet_sample_ids(sheet, colnames(mat))
      common <- intersect(colnames(mat), sample_ids)
      if (vld_notify_fail(length(common) >= 3, "Fewer than 3 samples matched between the uploaded matrix and sample sheet.")) return()
      mat <- mat[, common, drop = FALSE]; sheet <- sheet[match(common, sample_ids), , drop = FALSE]
      pheno_col <- input$val_upload_pheno_col
      if (vld_notify_fail(!is.null(pheno_col) && nzchar(pheno_col) && pheno_col %in% colnames(sheet), "Select a phenotype/class column from the uploaded sample sheet.")) return()
      grp_raw <- trimws(as.character(sheet[[pheno_col]]))
      ref_lab <- rm$ref_level; comp_lab <- rm$comp_level
      if (vld_notify_fail(all(c(ref_lab, comp_lab) %in% grp_raw),
        sprintf("The most recently trained model's class labels (%s / %s) were not both found in the chosen phenotype column - relabel the sample sheet to match.", ref_lab, comp_lab))) return()
      sex_sel <- input$val_upload_sex_stratum %||% "all"
      sex_keep <- rep(TRUE, nrow(sheet))
      if (!identical(sex_sel, "all")) {
        sex_col <- input$val_upload_sex_col
        if (vld_notify_fail(!is.null(sex_col) && nzchar(sex_col) && sex_col %in% colnames(sheet),
              "Select a sex column from the uploaded sample sheet to filter by Female/Male, or choose \"All samples\".")) return()
        sex_norm <- dxm_normalize_sex(sheet[[sex_col]])
        target <- if (identical(sex_sel, "male")) "M" else "F"
        sex_keep <- !is.na(sex_norm) & sex_norm == target
        if (vld_notify_fail(sum(sex_keep) > 0, sprintf("No samples matched sex = %s in the selected sex column.", dxm_sex_label(sex_sel)))) return()
      }
      keep <- sex_keep & (grp_raw %in% c(ref_lab, comp_lab))
      if (vld_notify_fail(sum(keep) >= 3, "Fewer than 3 samples match the trained model's class labels (and sex stratum, if selected).")) return()
      mat <- mat[, keep, drop = FALSE]; grp_raw <- grp_raw[keep]
      if (vld_notify_fail(sum(duplicated(rownames(mat))) == 0, "Uploaded matrix has duplicated CpG IDs.")) return()
      Mv <- if (identical(input$val_upload_scale, "m")) mat else dxm_beta_to_m(mat)
      y <- factor(ifelse(grp_raw == comp_lab, DXM_POS, DXM_NEG), levels = c(DXM_NEG, DXM_POS))

      cohort$loaded <- TRUE; cohort$source <- "upload"; cohort$label <- sprintf("Uploaded: %s", input$val_upload_matrix$name)
      cohort$m <- Mv; cohort$y <- y; cohort$sample_ids <- colnames(Mv); cohort$had_id_col <- had_id_col
      cohort$n_missing <- sum(is.na(Mv))
      cohort$normalization_note <- "Unknown - user-uploaded data; normalization method cannot be verified by this app"
      cohort$normalization_status <- "Warning"
      cohort$platform <- "Not specified by the uploaded data - platform compatibility could not be verified"
      cohort$platform_status <- "Warning"
      cohort$accession <- input$val_upload_matrix$name
      cohort$scale_declared <- input$val_upload_scale %||% "beta"
      cohort$n_samples <- ncol(Mv)
      cohort$ref_lab <- ref_lab; cohort$comp_lab <- comp_lab

      vhist$rows[[length(vhist$rows) + 1]] <- list(label = cohort$label, n = ncol(Mv), source = "Upload", loaded_at = format(Sys.time(), "%H:%M:%S"))
      showNotification(sprintf("Loaded uploaded cohort: %d samples.", ncol(Mv)), type = "message")
    })

    cohort_compare_table_data <- reactive({
      req(cohort$loaded, length(avail_models()) > 0)
      m <- ref_model()
      data.frame(
        Field = c("Role", "Source / accession", "N samples", sprintf("%s (reference)", m$ref_level), sprintf("%s (comparison)", m$comp_level), "CpG rows available"),
        `Training Cohort` = c("Training (model development)",
                               if (identical(m$mode, "preloaded")) "Training cohort (internal, preloaded)" else "User upload",
                               m$train_n, unname(m$train_class_table[DXM_NEG] %||% NA), unname(m$train_class_table[DXM_POS] %||% NA), NA),
        `Independent Validation Cohort` = c("Independent validation - never used for training or tuning",
                               cohort$accession %||% cohort$label, cohort$n_samples,
                               sum(cohort$y == DXM_NEG), sum(cohort$y == DXM_POS), nrow(cohort$m)),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })

    output$cohort_compare_gate_ui <- renderUI({
      if (length(avail_models()) == 0) return(NULL)
      if (!isTRUE(cohort$loaded)) return(p(class = "text-muted", "Load a validation cohort to see the training-vs-validation comparison (shown for the most recently trained model)."))
      NULL
    })

    output$cohort_compare_table <- DT::renderDataTable({
      DT::datatable(cohort_compare_table_data(), rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    output$cohort_overlap_ui <- renderUI({
      req(cohort$loaded, length(avail_models()) > 0)
      m <- ref_model()
      ov <- vld_sample_overlap(if (isTRUE(cohort$had_id_col)) m$train_sample_ids else character(0),
                                if (isTRUE(cohort$had_id_col)) cohort$sample_ids else character(0))
      p(class = if (identical(ov$status, "overlap")) "text-danger" else "submodule-desc",
        icon(if (identical(ov$status, "overlap")) "triangle-exclamation" else "circle-info"), " ", ov$message)
    })

    output$cohort_history_table <- DT::renderDataTable({
      req(length(vhist$rows) > 0)
      df <- do.call(rbind, lapply(seq_along(vhist$rows), function(i) {
        r <- vhist$rows[[i]]
        data.frame(Cohort = sprintf("Cohort %d", i), Label = r$label, N = r$n, Source = r$source, `Loaded at` = r$loaded_at, check.names = FALSE)
      }))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    ## ---- Compatibility ---------------------------------------------------

    compat_overview <- reactive({
      models <- avail_models(); req(cohort$loaded, length(models) > 0)
      do.call(rbind, lapply(models, function(m) {
        a <- vld_feature_alignment(m$feature_ids, rownames(cohort$m))
        data.frame(Model = m$label, `Analysis type` = m$analysis_type, `Required CpGs` = a$required,
                   Matched = a$matched, Missing = a$n_missing, `Overlap %` = round(a$overlap_pct, 1),
                   Status = if (a$ok) "Pass" else "Fail", check.names = FALSE, stringsAsFactors = FALSE)
      }))
    })

    output$compat_gate_ui <- renderUI({
      models <- avail_models()
      if (length(models) == 0 || !isTRUE(cohort$loaded)) {
        return(p(class = "text-muted", "Load at least one trained model and a validation cohort on the Datasets tab first."))
      }
      ov <- compat_overview()
      n_fail <- sum(ov$Status == "Fail")
      if (n_fail == nrow(ov)) {
        return(div(class = "alert alert-danger", icon("ban"),
          " Validation is BLOCKED for every trained model: none of them have all their required CpGs present in this validation cohort. Load a cohort containing every required CpG, or choose a different validation cohort."))
      }
      if (n_fail > 0) {
        return(div(class = "alert alert-warning", icon("triangle-exclamation"),
          sprintf(" %d of %d model(s) are blocked due to missing CpGs (see table below); the remaining %d can be validated.", n_fail, nrow(ov), nrow(ov) - n_fail)))
      }
      div(class = "alert alert-success", icon("circle-check"),
          sprintf(" All required CpGs are present for all %d model(s). Safe to proceed to any model tab.", nrow(ov)))
    })

    output$compat_overview_table <- DT::renderDataTable({
      DT::datatable(compat_overview(), rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    compat_selected_model <- reactive({ req(input$compat_model_select); avail_models()[[input$compat_model_select]] })

    output$compat_kpi_ui <- renderUI({
      m <- compat_selected_model(); req(m, cohort$loaded)
      a <- vld_feature_alignment(m$feature_ids, rownames(cohort$m))
      fluidRow(
        valueBox(a$required, "Required CpGs", icon = icon("dna"), color = "blue", width = 2),
        valueBox(a$available, "Available in validation cohort", icon = icon("table"), color = "light-blue", width = 3),
        valueBox(a$matched, "Matched CpGs", icon = icon("check"), color = if (a$ok) "green" else "yellow", width = 2),
        valueBox(a$n_missing, "Missing CpGs", icon = icon("triangle-exclamation"), color = if (a$n_missing == 0) "green" else "red", width = 2),
        valueBox(sprintf("%.1f%%", a$overlap_pct), "Feature overlap", icon = icon("percent"), color = if (a$ok) "green" else "yellow", width = 3)
      )
    })

    output$compat_table <- DT::renderDataTable({
      m <- compat_selected_model(); req(m, cohort$loaded)
      DT::datatable(vld_compat_table(m, cohort, vld_feature_alignment(m$feature_ids, rownames(cohort$m))), rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    ## ---- Per-model tabs ----------------------------------------------

    model_states <- lapply(DXM_MODEL_SPECS, function(spec) {
      vms <- reactiveValues(ran = FALSE, ran_at = NULL,
                             ext_prob = NULL, ext_roc = NULL, ext_metrics = NULL, ext_confusion = NULL, ext_calib = NULL, ext_y = NULL,
                             ci_sens = NULL, ci_spec = NULL, ci_acc = NULL, n_used = NULL, n_dropped = NULL,
                             roc_generated = FALSE, calib_generated = FALSE, last_roc_plot = NULL, last_calib_plot = NULL)
      vld_register_model_server(spec$id, spec, input, output, session, ns, avail_models, cohort, vms, vruns)
      vms
    })
    names(model_states) <- names(DXM_MODEL_SPECS)

    ## ---- Model Comparison --------------------------------------------

    output$compare_ui <- renderUI({
      keys <- names(shiny::reactiveValuesToList(vruns))
      if (length(keys) == 0) return(p(class = "text-muted", "No completed external validation runs yet - run external validation on any model tab first."))
      tagList(
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Model Comparison",
          p(class = "submodule-desc", icon("circle-info"),
            " \"External AUC\" is measured on the independent external validation cohort loaded on the Datasets tab - never the internal train/test split."),
          DT::dataTableOutput(ns("compare_table")),
          downloadButton(ns("compare_download"), "Download comparison (CSV)", class = "btn-sm")
        ),
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Compare ROC Curves",
          selectizeInput(ns("compare_select"), "Select runs to compare",
            choices = stats::setNames(keys, vapply(keys, function(k) sprintf("%s - %s", vruns[[k]]$label, vruns[[k]]$analysis_type), character(1))), multiple = TRUE),
          selectInput(ns("compare_curve"), "Curve to compare", choices = c("External Test" = "ext", "Training" = "train", "Cross-Validated" = "cv")),
          actionButton(ns("compare_roc_btn"), "Generate ROC Comparison", icon = icon("chart-area"), class = "btn-primary"),
          br(), br(),
          if (isTRUE(compare_state$generated)) plotOutput(ns("compare_roc_plot"))
        )
      )
    })

    output$compare_table <- DT::renderDataTable({
      all_runs <- shiny::reactiveValuesToList(vruns)
      req(length(all_runs) > 0)
      tbl <- do.call(rbind, lapply(all_runs, function(r) data.frame(
        Model = r$label, `Feature set` = sprintf("%s (%d)", r$analysis_type, length(r$feature_ids)),
        `Train AUC` = round(r$train_metrics$auc %||% NA_real_, 3), `CV AUC` = round(r$cv_roc$mean_auc %||% NA_real_, 3),
        `External AUC` = round(r$ext_metrics$auc, 3), `N validated` = r$n_used,
        Threshold = round(r$threshold, 3), `Ran at` = r$ran_at, check.names = FALSE)))
      DT::datatable(tbl, rownames = FALSE, filter = "top", options = list(pageLength = 10))
    })

    output$compare_download <- downloadHandler(
      filename = function() "methylomics_validation_model_comparison.csv",
      content = function(file) {
        all_runs <- shiny::reactiveValuesToList(vruns)
        tbl <- do.call(rbind, lapply(all_runs, function(r) data.frame(
          model = r$label, analysis_type = r$analysis_type, n_features = length(r$feature_ids),
          features = paste(r$feature_ids, collapse = ";"), threshold = r$threshold,
          train_auc = r$train_metrics$auc %||% NA_real_, cv_auc = r$cv_roc$mean_auc %||% NA_real_,
          external_auc = r$ext_metrics$auc, n_validated = r$n_used, ran_at = r$ran_at)))
        utils::write.csv(tbl, file, row.names = FALSE)
      }
    )

    observeEvent(input$compare_roc_btn, {
      sel <- input$compare_select
      validate(need(length(sel) > 0, "Select at least one run to compare."))
      all_runs <- shiny::reactiveValuesToList(vruns)
      bundles <- stats::setNames(lapply(sel, function(k) switch(input$compare_curve,
                    ext = all_runs[[k]]$ext_roc, train = all_runs[[k]]$train_roc, cv = all_runs[[k]]$cv_roc$overall)),
                    vapply(sel, function(k) sprintf("%s (%s)", all_runs[[k]]$label, all_runs[[k]]$analysis_type), character(1)))
      compare_state$bundles <- bundles; compare_state$generated <- TRUE
    })
    output$compare_roc_plot <- renderPlot({
      req(compare_state$generated)
      dxm_plot_roc_compare(compare_state$bundles, sprintf("ROC Comparison (%s)",
        switch(input$compare_curve, ext = "External Test", train = "Training", cv = "Cross-Validated")))
    })

    ## ---- Test External Data (external cohort snapshot) ------------------------

    output$testdata_ui <- renderUI({
      if (!isTRUE(cohort$loaded)) return(p(class = "text-muted", "Load an external validation cohort on the Datasets tab first."))
      box(width = 12, status = "primary", solidHeader = TRUE, title = "Test External Data",
        p(sprintf("%d external validation samples (%s=%d, %s=%d).", cohort$n_samples, cohort$ref_lab %||% "reference",
                   sum(cohort$y == DXM_NEG), cohort$comp_lab %||% "comparison", sum(cohort$y == DXM_POS))),
        tags$ul(
          tags$li(sprintf("CpG rows available in the validation cohort matrix: %d", nrow(cohort$m))),
          tags$li(sprintf("Missing values in the validation cohort matrix: %d", cohort$n_missing))
        ),
        hr(),
        h5("Feature availability per trained model"),
        DT::dataTableOutput(ns("testdata_feature_table"))
      )
    })

    output$testdata_feature_table <- DT::renderDataTable({
      models <- avail_models(); req(cohort$loaded, length(models) > 0)
      df <- do.call(rbind, lapply(models, function(m) {
        a <- vld_feature_alignment(m$feature_ids, rownames(cohort$m))
        data.frame(Model = m$label, `Analysis type` = m$analysis_type, `Training features` = a$required,
                   `Available in external cohort` = a$available, `Shared features` = a$matched,
                   `Unmatched (dropped) features` = a$n_missing, check.names = FALSE, stringsAsFactors = FALSE)
      }))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    ## ---- Export ------------------------------------------------------

    output$export_ui <- renderUI({
      keys <- names(shiny::reactiveValuesToList(vruns))
      if (length(keys) == 0) return(p(class = "text-muted", "No completed external validation runs yet - run external validation on any model tab first."))
      box(width = 12, status = "primary", solidHeader = TRUE, title = "Export",
        downloadButton(ns("export_metrics_csv"), "Download all metrics (CSV)"),
        downloadButton(ns("export_featavail_csv"), "Download feature availability per model (CSV)"))
    })
    output$export_metrics_csv <- downloadHandler(
      filename = function() "methylomics_validation_metrics.csv",
      content = function(file) {
        all_runs <- shiny::reactiveValuesToList(vruns)
        tbl <- do.call(rbind, lapply(all_runs, function(r) data.frame(
          model = r$label, analysis_type = r$analysis_type, n_features = length(r$feature_ids),
          features = paste(r$feature_ids, collapse = ";"), threshold = r$threshold,
          train_auc = r$train_metrics$auc %||% NA_real_, cv_auc = r$cv_roc$mean_auc %||% NA_real_,
          external_auc = r$ext_metrics$auc, n_validated = r$n_used, ran_at = r$ran_at)))
        utils::write.csv(tbl, file, row.names = FALSE)
      }
    )
    output$export_featavail_csv <- downloadHandler(
      filename = function() "methylomics_validation_feature_availability.csv",
      content = function(file) {
        models <- avail_models(); req(cohort$loaded, length(models) > 0)
        df <- do.call(rbind, lapply(models, function(m) {
          a <- vld_feature_alignment(m$feature_ids, rownames(cohort$m))
          data.frame(model = m$label, analysis_type = m$analysis_type, features_required = a$required,
                     features_available = a$available, features_matched = a$matched, features_missing = a$n_missing)
        }))
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    invisible(NULL)
  })
}
