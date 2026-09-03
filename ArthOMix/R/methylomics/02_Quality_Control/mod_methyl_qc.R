## R/methylomics/02_Quality_Control/mod_methyl_qc.R - Methylomics Quality Control sub-module (UI + server).
## Each tab is an independent, button-driven eventReactive over the shared
## methyl_dataset; nothing here mutates methyl_dataset itself.
mod_methyl_qc_config <- list(
  id = "qc", title = "Quality Control", icon = "magnifying-glass-chart", group = "Data",
  description = "Performs probe and sample QC"
)

mod_methyl_qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    withSpinner(uiOutput(ns("body_ui")), color = "#2563EB", type = 6),
    tags$hr(),
    uiOutput(ns("default_ui_wrap"))
  )
}

qc_tab_title <- function(ic, label) tagList(icon(ic), " ", label)

mod_methyl_qc_server <- function(id, methyl_dataset, methyl_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    svg_download_link <- function(plot_id) {
      actionLink(ns(paste0(plot_id, "_svgdl")),
                 label = tagList(icon("file-image"), " Download SVG"),
                 style = "font-size:12px; margin-left:12px; color:#2563EB;")
    }
    wire_svg_download <- function(plot_id, filename) {
      observeEvent(input[[paste0(plot_id, "_svgdl")]], {
        shinyjs::runjs(sprintf(
          "if (window.Plotly) { var gd = document.getElementById('%s'); if (gd) Plotly.downloadImage(gd, {format: 'svg', filename: '%s'}); }",
          ns(plot_id), filename
        ))
      }, ignoreInit = TRUE)
    }

    register_has_run_gate <- function(gate_id, has_run_flag_fn, result_output_id, not_run_message) {
      output[[gate_id]] <- renderUI({
        if (isTRUE(has_run_flag_fn())) uiOutput(ns(result_output_id))
        else div(class = "card", p(class = "empty-note", icon("circle-info"), not_run_message))
      })
    }

    plot_shown <- reactiveValues()
    lazy_plot_ids <- c("viz_pca_plot", "viz_pca3d_plot", "mds_plot", "density_plot", "boxplot_plot",
                        "violin_plot", "corr_heatmap", "mean_sd_plot",
                        "detp_heatmap", "beadcount_dist", "control_heatmap",
                        "pca_outlier_plot", "outlier_diagnostic_plot", "dendro_plot")
    lapply(lazy_plot_ids, function(pid) {
      plot_shown[[pid]] <- FALSE
      observeEvent(input[[paste0(pid, "_gen_btn")]], {
        plot_shown[[pid]] <- TRUE
        shinyjs::hide(id = paste0(pid, "_gen_wrap"))
        shinyjs::show(id = paste0(pid, "_out_wrap"))
        shinyjs::delay(200, shinyjs::runjs(sprintf(
          "if (window.Plotly) { var gd = document.getElementById('%s'); if (gd && gd.data) Plotly.Plots.resize(gd); }",
          ns(pid)
        )))
      }, ignoreInit = TRUE)
    })

    lazy_plot_ui <- function(plot_id, height, label = "Generate plot", with_svg = TRUE, kind = c("plotly", "plot")) {
      kind <- match.arg(kind)
      out_widget <- if (kind == "plotly") plotly::plotlyOutput(ns(plot_id), height = height)
                    else plotOutput(ns(plot_id), height = height)
      tagList(
        div(id = ns(paste0(plot_id, "_gen_wrap")),
            actionButton(ns(paste0(plot_id, "_gen_btn")), label, icon = icon("chart-line"), class = "btn-primary btn-sm")),
        shinyjs::hidden(div(id = ns(paste0(plot_id, "_out_wrap")),
            if (with_svg) div(svg_download_link(plot_id)),
            withSpinner(out_widget, color = "#2563EB", type = 6)
        ))
      )
    }

    run_info_line <- function(run_at, extra = NULL) {
      p(class = "empty-note", icon("clock"),
        sprintf("Run at %s.%s", format(run_at, "%Y-%m-%d %H:%M:%S"), if (!is.null(extra)) paste0(" ", extra) else ""))
    }

    output$default_ui_wrap <- renderUI({
      if (!isTRUE(methyl_dataset$preloaded)) return(NULL)
      tagList(
        actionLink(ns("toggle_historical_btn"),
                   label = tagList(icon("clock-rotate-left"), " Historical pipeline reference (reproduced from the completed run, not recomputed) - click to show/hide"),
                   class = "empty-note", style = "cursor:pointer;"),
        shinyjs::hidden(div(id = ns("historical_wrap"), style = "margin-top:8px;",
          withSpinner(uiOutput(ns("default_ui")), color = "#2563EB", type = 6)
        ))
      )
    })
    observeEvent(input$toggle_historical_btn, {
      shinyjs::toggle(id = "historical_wrap")
    })

    output$default_ui <- renderUI({
      if (!isTRUE(methyl_dataset$preloaded)) {
        return(div(class = "card",
          div(class = "card-title", icon("circle-info"), "No default analysis loaded"),
          p(class = "submodule-desc", "Load the preloaded whole-blood dataset on the Methylomics Dataset tab to see its default, already-completed QC here.")
        ))
      }
      pheno <- load_default_meth_pheno()
      sexcheck <- load_default_meth_qc_sexcheck()
      req(pheno, sexcheck)
      tagList(
        div(class = "empty-note", icon("flask"),
            "Historical reference: sample-level QC and probe-level filtering for the preloaded whole-blood dataset, reproduced from the completed pipeline run - nothing here is recomputed. The live, interactive tool below runs the same kinds of checks against the real matrix instead, one method at a time as you run each tab."),
        div(class = "card",
            div(class = "card-title", icon("users"), "Cohort composition"),
            radioButtons(ns("qc_sex"), "Stratum", inline = TRUE,
                         choices = c("Female" = "F", "Male" = "M", "Both" = "both"), selected = "both"),
            DT::dataTableOutput(ns("qc_cohort_table"))
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-line"), "PCA outliers & chrY sex-check"),
            fluidRow(
              valueBox(sum(sexcheck$outlier), "PCA outliers (>5 MAD, within-sex)", icon = icon("triangle-exclamation"), color = if (sum(sexcheck$outlier) > 0) "yellow" else "green", width = 4),
              valueBox(sum(sexcheck$sex_mismatch), "Sex mismatches (chrY vs reported)", icon = icon("venus-mars"), color = if (sum(sexcheck$sex_mismatch) > 0) "red" else "green", width = 4),
              valueBox(nrow(sexcheck), "Samples checked", icon = icon("vial"), color = "light-blue", width = 4)
            ),
            p(class = "submodule-desc", "Outlier detection is run separately within each sex (whole-blood methylation separates strongly by sex on early PCs, so pooling would misflag one sex as \"outliers\" relative to the other)."),
            DT::dataTableOutput(ns("qc_outlier_table"))
        ),
        div(class = "card",
            div(class = "card-title", icon("filter"), "Probe-filtering cascade"),
            p(class = "submodule-desc", "Applied as five sequential, independently logged criteria (ChAMP-style: cg-prefix restriction, Zhou et al. 2017 MASK_general, multi-hit removal, sex-chromosome removal, >5% missingness) - a cohort-wide, probe-level QC pass run once on all 689 samples before the sex-stratified analysis, so it is intentionally the same for every Stratum selection above (probe filtering doesn't depend on which sex is being viewed; only the two sample-level tables above it do)."),
            DT::dataTableOutput(ns("qc_cascade_table"))
        )
      )
    })

    output$qc_cohort_table <- DT::renderDataTable({
      pheno <- load_default_meth_pheno()
      df <- if (identical(input$qc_sex, "both")) pheno else pheno[pheno$sex == input$qc_sex, , drop = FALSE]
      tbl <- as.data.frame.matrix(table(df$group, df$sex))
      tbl <- cbind(group = rownames(tbl), tbl)
      DT::datatable(tbl, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$qc_outlier_table <- DT::renderDataTable({
      sexcheck <- load_default_meth_qc_sexcheck()
      df <- if (identical(input$qc_sex, "both")) sexcheck else sexcheck[sexcheck$sex == input$qc_sex, , drop = FALSE]
      df <- df[df$outlier | df$sex_mismatch, , drop = FALSE]
      if (nrow(df) == 0) df <- data.frame(message = "No flagged samples in this stratum.")
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$qc_cascade_table <- DT::renderDataTable({
      DT::datatable(METH_QC_PROBE_CASCADE, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatCurrency(columns = c("retained", "removed"), currency = "", interval = 3, mark = ",", digits = 0)
    })

    anno_result <- reactive({
      req(methyl_dataset$array_type)
      methyl_get_annotation(methyl_dataset$array_type)
    })

    is_illumina_array <- reactive({
      req(methyl_dataset$array_type)
      methyl_dataset$array_type %in% METHYL_ARRAY_TYPES_ILLUMINA
    })

    manual_exclude <- reactiveVal(character(0))
    overview_has_run <- reactiveVal(FALSE)
    sample_qc_has_run <- reactiveVal(FALSE)
    probe_qc_has_run <- reactiveVal(FALSE)
    sex_qc_has_run <- reactiveVal(FALSE)
    batch_qc_has_run <- reactiveVal(FALSE)
    outlier_qc_has_run <- reactiveVal(FALSE)

    observeEvent(methyl_dataset$beta, {
      manual_exclude(character(0))
      overview_has_run(FALSE); sample_qc_has_run(FALSE); probe_qc_has_run(FALSE)
      sex_qc_has_run(FALSE); batch_qc_has_run(FALSE); outlier_qc_has_run(FALSE)
    }, ignoreNULL = TRUE)

    register_has_run_gate("overview_gate", overview_has_run, "overview_result_ui",
      "Not run yet - click \"Run Overview QC\" above to see it.")
    register_has_run_gate("sample_qc_gate", sample_qc_has_run, "sample_qc_result_ui",
      "Not run yet - set options above and click \"Run Sample QC\".")
    register_has_run_gate("probe_qc_gate", probe_qc_has_run, "probe_qc_result_ui",
      "Not run yet - set filters above and click \"Run Probe QC\". Changing a filter alone never recomputes anything - only clicking this button does.")
    register_has_run_gate("sex_qc_gate", sex_qc_has_run, "sex_qc_result_ui",
      "Not run yet - click \"Run Sex QC\" above.")
    register_has_run_gate("batch_qc_gate", batch_qc_has_run, "batch_result_ui",
      "Not run yet - choose a method above and click \"Run Batch QC\".")
    register_has_run_gate("outlier_qc_gate", outlier_qc_has_run, "outlier_result_ui",
      "Not run yet - choose method(s) above and click \"Run Outlier Detection\".")

    output$live_stratum_ui <- renderUI({
      req(methyl_dataset$beta, methyl_dataset$sample_sheet)
      all_ids <- colnames(methyl_dataset$beta)
      if (is.null(input$live_group_col) || !nzchar(input$live_group_col)) {
        return(radioButtons(ns("live_stratum"), "Stratum",
                             choices = stats::setNames("__all__", sprintf("All (n=%d)", length(all_ids))),
                             selected = "__all__", inline = TRUE))
      }
      sheet <- methyl_dataset$sample_sheet
      id_col <- intersect(c("sample", "Sample", "sample_id", "Sample_ID"), colnames(sheet))
      sample_ids <- if (length(id_col) > 0) as.character(sheet[[id_col[1]]]) else rownames(sheet)
      grp <- stats::setNames(as.character(sheet[[input$live_group_col]]), sample_ids)[all_ids]
      counts <- table(grp, useNA = "no")
      lvl_names <- names(counts) %||% character(0)
      lvl_choices <- stats::setNames(lvl_names, sprintf("%s (n=%d)", lvl_names, as.integer(counts)))
      choices <- c(stats::setNames("__all__", sprintf("All (n=%d)", length(all_ids))), lvl_choices)
      radioButtons(ns("live_stratum"), "Stratum", choices = choices, selected = "__all__", inline = TRUE)
    })

    current_subgroup <- reactive({
      req(methyl_dataset$beta)
      sg <- methyl_qc_subgroup_filter(methyl_dataset$beta, methyl_dataset$sample_sheet, input$live_group_col, input$live_stratum)
      methyl_apply_manual_exclude(sg, manual_exclude())
    })
    current_rg_subset <- reactive({
      sg <- current_subgroup()
      if (is.null(methyl_dataset$rg_set)) return(NULL)
      tryCatch(methyl_dataset$rg_set[, sg$included], error = function(e) NULL)
    })
    current_mset_subset <- reactive({
      sg <- current_subgroup()
      if (is.null(methyl_dataset$mset)) return(NULL)
      tryCatch(methyl_dataset$mset[, sg$included], error = function(e) NULL)
    })

    stratum_all_samples <- reactive({
      req(methyl_dataset$beta)
      sg <- methyl_qc_subgroup_filter(methyl_dataset$beta, methyl_dataset$sample_sheet, input$live_group_col, input$live_stratum)
      data.frame(sample = sg$included, call_rate = methyl_sample_call_rate(sg$mat), stringsAsFactors = FALSE)
    })

    output$manual_table <- DT::renderDataTable({
      df <- stratum_all_samples()
      df$manually_excluded <- df$sample %in% manual_exclude()
      sel <- which(df$sample %in% manual_exclude())
      DT::datatable(df, rownames = FALSE, selection = list(mode = "multiple", selected = sel),
                    options = list(scrollX = TRUE, pageLength = 10),
                    class = "stripe hover compact") %>%
        DT::formatRound(columns = "call_rate", digits = 4)
    })

    observeEvent(input$manual_apply_btn, {
      df <- stratum_all_samples()
      manual_exclude(df$sample[input$manual_table_rows_selected %||% integer(0)])
    })
    observeEvent(input$manual_clear_btn, {
      manual_exclude(character(0))
    })

    output$manual_exclude_summary <- renderUI({
      ex <- manual_exclude()
      if (length(ex) == 0) return(p(class = "empty-note", icon("circle-check"), "No samples are manually excluded."))
      p(class = "empty-note", icon("user-slash"),
        sprintf("%d sample(s) manually excluded: %s", length(ex), paste(ex, collapse = ", ")))
    })

    ruvm_available <- reactive({
      !is.null(methyl_dataset$rg_set) && requireNamespace("missMethyl", quietly = TRUE)
    })

    output$body_ui <- renderUI({
      if (is.null(methyl_dataset$beta)) {
        return(div(class = "card",
          div(class = "card-title", icon("upload"), "Live QC dashboard"),
          p(class = "submodule-desc", "Upload a beta/M-value matrix or raw IDAT files on the Methylomics Dataset tab - or load the preloaded whole-blood dataset, if this deployment has the live matrix available - to run this interactively.")
        ))
      }
      div(class = "tx-menu-wrap",
        tabsetPanel(
          id = ns("tabs"), type = "tabs", header = tagList(tags$hr()),
          tabPanel(value = "Overview", title = qc_tab_title("gauge-high", "Overview"),
                   br(), withSpinner(uiOutput(ns("overview_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Sample QC", title = qc_tab_title("users", "Sample QC"),
                   br(), withSpinner(uiOutput(ns("sample_qc_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Probe QC", title = qc_tab_title("filter", "Probe QC"),
                   br(), withSpinner(uiOutput(ns("probe_qc_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Sex QC", title = qc_tab_title("venus-mars", "Sex QC"),
                   br(), withSpinner(uiOutput(ns("sex_qc_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Batch QC", title = qc_tab_title("wand-magic-sparkles", "Batch QC"),
                   br(), withSpinner(uiOutput(ns("batch_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Outlier QC", title = qc_tab_title("magnifying-glass", "Outlier QC"),
                   br(), withSpinner(uiOutput(ns("outlier_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Visualizations", title = qc_tab_title("chart-line", "Visualizations"),
                   br(), withSpinner(uiOutput(ns("viz_ui")), color = "#2563EB", type = 6)),
          tabPanel(value = "Reports & Export", title = qc_tab_title("file-export", "Reports & Export"),
                   br(), withSpinner(uiOutput(ns("reports_ui")), color = "#2563EB", type = 6))
        )
      )
    })

    output$overview_controls_ui <- renderUI({
      req(methyl_dataset$beta)
      sheet <- methyl_dataset$sample_sheet
      tagList(
        div(class = "empty-note", icon("circle-info"),
            sprintf("Dataset loaded - %s probe(s) x %s sample(s) (%s). Nothing has been analyzed yet: open a tab below, set its parameters, and click that tab's own Run button.",
                    format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta), methyl_dataset$source %||% "unlabeled dataset")),
        div(class = "card", style = "padding:10px 16px;",
            fluidRow(
              column(6,
                if (!is.null(sheet)) {
                  sheet_cols <- colnames(sheet)
                  default_col <- intersect(METHYL_SEX_COL_CANDIDATES, sheet_cols)
                  selectInput(ns("live_group_col"), "Subgroup column", choices = c("All samples (no subgroup)" = "", sheet_cols),
                              selected = if (length(default_col) > 0) default_col[1] else "", width = "100%")
                } else p(class = "empty-note", style = "margin:6px 0;", icon("circle-info"), "No sample sheet - every QC method below runs on all samples.")
              ),
              column(6, withSpinner(uiOutput(ns("live_stratum_ui")), color = "#2563EB", type = 6, proxy.height = "40px"))
            )
        )
      )
    })

    output$overview_summary_ui <- renderUI({
      req(methyl_dataset$beta)
      sheet <- methyl_dataset$sample_sheet
      batch_cols <- methyl_batch_columns(sheet)
      sex_cols <- if (!is.null(sheet)) intersect(METHYL_SEX_COL_CANDIDATES, colnames(sheet)) else character(0)
      sg <- current_subgroup()
      tagList(
        div(class = "card",
            div(class = "card-title", icon("database"), "Loaded dataset"),
            fluidRow(
              valueBox(format(length(sg$included), big.mark = ","), sprintf("Samples in scope (%s)", sg$label), icon = icon("vial"), color = "purple", width = 4),
              valueBox(format(nrow(methyl_dataset$beta), big.mark = ","), "CpGs / probes", icon = icon("dna"), color = "aqua", width = 4),
              valueBox(if (identical(methyl_dataset$input_scale, "beta")) "Beta values" else "M-values", "Scale", icon = icon("ruler"), color = "light-blue", width = 4)
            ),
            p(class = "submodule-desc", methyl_dataset$source %||% "Unlabeled dataset."),
            p(class = "submodule-desc",
              sprintf("Metadata columns available: %s.", if (!is.null(sheet)) paste(colnames(sheet), collapse = ", ") else "none - no sample sheet uploaded.")),
            p(class = "submodule-desc",
              sprintf("Batch/chip column detected: %s. Sex/gender column detected: %s.",
                      if (length(batch_cols) > 0) paste(batch_cols, collapse = ", ") else "none",
                      if (length(sex_cols) > 0) paste(sex_cols, collapse = ", ") else "none"))
        ),
        div(class = "card",
            div(class = "card-title", icon("gauge-high"), "Basic QC pass"),
            p(class = "submodule-desc", "A lightweight, self-contained call-rate/missingness summary for the current stratum - independent of every other tab below. Nothing is calculated until you click Run."),
            actionButton(ns("run_overview_btn"), "Run Overview QC", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("overview_gate"))
      )
    })

    output$overview_ui <- renderUI({
      req(methyl_dataset$beta)
      tagList(
        uiOutput(ns("overview_controls_ui")),
        withSpinner(uiOutput(ns("overview_summary_ui")), color = "#2563EB", type = 6)
      )
    })

    overview_result <- eventReactive(input$run_overview_btn, {
      sg <- current_subgroup()
      validate(need(length(sg$included) >= 1, "No samples in the current stratum."))
      mat <- sg$mat
      call_rate <- methyl_sample_call_rate(mat)
      list(n_samples = ncol(mat), n_probes = nrow(mat),
           overall_missing_pct = 100 * mean(is.na(mat)),
           median_call_rate = stats::median(call_rate, na.rm = TRUE),
           range_check = methyl_check_beta_range(mat, methyl_dataset$input_scale %||% "beta"),
           subgroup = sg, run_at = Sys.time())
    })
    observeEvent(overview_result(), overview_has_run(TRUE))

    output$overview_result_ui <- renderUI({
      req(overview_result())
      o <- overview_result()
      status <- methyl_qc_status_badge(o)
      status_color <- c(pass = "green", warning = "yellow", fail = "red")[status$status]
      status_label <- c(pass = "Pass", warning = "Warning", fail = "Fail")[status$status]
      tagList(
        run_info_line(o$run_at, sprintf("%d sample(s), %d probe(s) analyzed.", o$n_samples, o$n_probes)),
        div(class = "card",
            div(class = "card-title", icon("gauge-high"), "Basic QC status"),
            fluidRow(
              valueBox(status_label, "Status", icon = icon(switch(status$status, pass = "circle-check", warning = "triangle-exclamation", fail = "circle-xmark")), color = status_color, width = 4),
              valueBox(sprintf("%.4f", o$median_call_rate), "Median sample call rate", icon = icon("vial"), color = "light-blue", width = 4),
              valueBox(sprintf("%.2f%%", o$overall_missing_pct), "Overall missingness", icon = icon("percent"), color = "yellow", width = 4)
            ),
            tags$ul(lapply(status$reasons, tags$li))
        )
      )
    })

    output$sample_qc_ui <- renderUI({
      req(methyl_dataset$beta)
      has_idat <- !is.null(methyl_dataset$rg_set) && !is.null(methyl_dataset$detp)
      sg <- current_subgroup()
      tagList(
        div(class = "card",
            div(class = "card-title", icon("venus-mars"), "Selected subgroup"),
            p(class = "empty-note", icon("filter"), sprintf("Every method below runs on: %s.", sg$label)),
            if (isTRUE(sg$low_n)) p(class = "empty-note", icon("triangle-exclamation"),
              sprintf("Only %d sample(s) in this stratum - PCA, clustering, and outlier detection will be unreliable with this few samples.", length(sg$included)))
        ),
        div(class = "card",
            div(class = "card-title", icon("users"), "Sample QC options"),
            fluidRow(
              column(6, numericInput(ns("call_rate_min"), "Minimum sample call rate", value = 0.95, min = 0, max = 1, step = 0.01)),
              column(6, numericInput(ns("sample_detp_thresh"), "Detection failure threshold (per-probe p-value)", value = 0.01, min = 0, max = 1, step = 0.01))
            ),
            fluidRow(
              column(6,
                checkboxInput(ns("f_failed_probe_pct"), "Maximum failed-probe percentage", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_failed_probe_pct")),
                  numericInput(ns("failed_probe_pct_max"), "Max failed probes per sample (%)", value = 5, min = 0, max = 100, step = 1))
              ),
              column(6,
                checkboxInput(ns("f_min_intensity"), "Minimum signal intensity", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_min_intensity")),
                  numericInput(ns("min_intensity_thresh"), "Minimum median log2 intensity", value = 10, step = 0.5))
              )
            ),
            if (!has_idat) p(class = "empty-note", icon("circle-info"),
              "This dataset has no raw IDAT input, so detection p-values and per-sample intensity aren't available - \"Maximum failed-probe percentage\", \"Minimum signal intensity\", bisulfite conversion efficiency, and median intensity below all show as unavailable. Call-rate filtering still works from the beta/M-value matrix alone."),
            actionButton(ns("run_sample_qc_btn"), "Run Sample QC", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        div(class = "card",
            div(class = "card-title", icon("table-list"), "Manual sample inclusion/exclusion"),
            p(class = "submodule-desc", "Select rows to exclude, then Apply - excluded samples are removed from the sample scope every QC method reads from, starting the next time each method's own Run button is clicked."),
            DT::dataTableOutput(ns("manual_table")),
            div(style = "margin-top:8px;",
                actionButton(ns("manual_apply_btn"), "Apply exclusions", icon = icon("check"), class = "btn-primary btn-sm"),
                actionButton(ns("manual_clear_btn"), "Clear all exclusions", icon = icon("rotate-left"), class = "btn-default btn-sm", style = "margin-left:8px;")
            ),
            uiOutput(ns("manual_exclude_summary"))
        ),
        uiOutput(ns("sample_qc_gate"))
      )
    })

    sample_qc_result <- eventReactive(input$run_sample_qc_btn, {
      sg <- current_subgroup()
      validate(need(length(sg$included) >= 1, "No samples in the current stratum."))
      mat <- sg$mat

      withProgress(message = "Running Sample QC", value = 0.1, {
        call_rate <- methyl_sample_call_rate(mat)
        sample_qc <- data.frame(sample = colnames(mat), call_rate = call_rate,
                                 call_rate_flag = call_rate < (input$call_rate_min %||% 0.95))

        incProgress(0.3, detail = "Raw-intensity metrics")
        rg <- current_rg_subset(); mset <- current_mset_subset()
        bisulfite <- if (!is.null(rg)) methyl_bisulfite_conversion(rg) else list(ok = FALSE, reason = "No IDAT input.")
        median_int <- if (!is.null(mset)) methyl_median_intensity(mset) else list(ok = FALSE, reason = "No IDAT input.")

        incProgress(0.3, detail = "Failed-probe % / intensity")
        failed_pct <- methyl_sample_failed_probe_pct(mat, methyl_dataset$detp, input$sample_detp_thresh %||% 0.01)
        if (isTRUE(failed_pct$ok)) sample_qc$failed_probe_pct <- failed_pct$pct[sample_qc$sample]
        if (isTRUE(input$f_failed_probe_pct) && isTRUE(failed_pct$ok)) {
          sample_qc$failed_probe_pct_flag <- sample_qc$failed_probe_pct > (input$failed_probe_pct_max %||% 5)
        }
        low_intensity <- methyl_sample_low_intensity(median_int, input$min_intensity_thresh %||% 10)
        if (isTRUE(low_intensity$ok)) sample_qc$intensity_score <- low_intensity$score[sample_qc$sample]
        if (isTRUE(input$f_min_intensity) && isTRUE(low_intensity$ok)) {
          sample_qc$low_intensity_flag <- low_intensity$low[sample_qc$sample]
        }
      })

      list(mat = mat, subgroup = sg, sample_qc = sample_qc, bisulfite = bisulfite, median_int = median_int,
           settings = list(call_rate_min = input$call_rate_min, sample_detp_thresh = input$sample_detp_thresh,
                            f_failed_probe_pct = input$f_failed_probe_pct, failed_probe_pct_max = input$failed_probe_pct_max,
                            f_min_intensity = input$f_min_intensity, min_intensity_thresh = input$min_intensity_thresh,
                            stratum_label = sg$label),
           run_at = Sys.time())
    })
    observeEvent(sample_qc_result(), sample_qc_has_run(TRUE))

    output$sample_qc_result_ui <- renderUI({
      req(sample_qc_result())
      r <- sample_qc_result()
      tagList(
        run_info_line(r$run_at, sprintf("%d sample(s) analyzed, stratum: %s.", ncol(r$mat), r$subgroup$label)),
        div(class = "card",
            div(class = "card-title", icon("users"), "Sample call-rate QC"),
            DT::dataTableOutput(ns("sample_qc_table")),
            if (!isTRUE(r$bisulfite$ok) && !isTRUE(r$median_int$ok) && identical(r$bisulfite$reason, r$median_int$reason)) {
              p(class = "empty-note", icon("circle-info"),
                sprintf("Bisulfite conversion efficiency and median intensity: %s", r$bisulfite$reason))
            } else tagList(
              if (isTRUE(r$bisulfite$ok)) tagList(h5("Bisulfite conversion efficiency"), DT::dataTableOutput(ns("bisulfite_table")))
              else p(class = "empty-note", icon("circle-info"), sprintf("Bisulfite conversion efficiency: %s", r$bisulfite$reason)),
              if (isTRUE(r$median_int$ok)) tagList(h5("Median intensity"), DT::dataTableOutput(ns("median_int_table")))
              else p(class = "empty-note", icon("circle-info"), sprintf("Median intensity: %s", r$median_int$reason))
            )
        )
      )
    })

    output$sample_qc_table <- DT::renderDataTable({
      req(sample_qc_result())
      df <- sample_qc_result()$sample_qc
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatRound(columns = "call_rate", digits = 4)
    })

    output$bisulfite_table <- DT::renderDataTable({
      req(sample_qc_result())
      pct <- sample_qc_result()$bisulfite$pct
      df <- data.frame(sample = names(pct), bisulfite_conversion_pct = as.numeric(pct))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatRound(columns = "bisulfite_conversion_pct", digits = 2)
    })

    output$median_int_table <- DT::renderDataTable({
      req(sample_qc_result())
      df <- sample_qc_result()$median_int$detail
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatRound(columns = setdiff(colnames(df), "sample"), digits = 2)
    })

    output$probe_qc_ui <- renderUI({
      req(methyl_dataset$beta)
      has_idat <- !is.null(methyl_dataset$rg_set)
      scale_is_beta <- identical(methyl_dataset$input_scale, "beta")
      tagList(
        div(class = "card",
            div(class = "card-title", icon("filter"), "Probe filtering options"),
            fluidRow(
              column(6,
                checkboxInput(ns("f_detp"), "Filter by detection p-value", value = has_idat),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_detp")),
                  numericInput(ns("detp_thresh"), "Detection p-value threshold", value = 0.01, min = 0, max = 1, step = 0.01)),
                if (!has_idat) p(class = "empty-note", icon("circle-info"), "Requires raw IDAT input - upload IDAT files on the Dataset tab to enable this."),
                checkboxInput(ns("f_beadcount"), "Filter by bead count", value = has_idat && !is.null(methyl_dataset$beadcount)),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_beadcount")),
                  numericInput(ns("beadcount_thresh"), "Minimum bead count", value = 3, min = 1, step = 1)),
                checkboxInput(ns("f_missing"), "Remove probes with missing values", value = TRUE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_missing")),
                  numericInput(ns("missing_max"), "Max missing fraction allowed per probe", value = 0, min = 0, max = 1, step = 0.05))
              ),
              column(6,
                if (isTRUE(is_illumina_array())) tagList(
                  checkboxInput(ns("f_snp"), "Remove SNP-overlapping probes", value = TRUE),
                  checkboxInput(ns("f_noncpg"), "Remove non-CpG (CpH) probes", value = TRUE),
                  radioButtons(ns("sexchr_mode"), "Sex-chromosome handling", inline = TRUE,
                               choices = c("Keep all" = "keep", "Remove X and Y" = "remove_xy", "Remove Y only" = "remove_y_only"),
                               selected = "keep")
                ) else p(class = "empty-note", icon("circle-info"),
                         "SNP / non-CpG / sex-chromosome probe filters only apply to Illumina array data and are hidden for this dataset type."),
                checkboxInput(ns("f_crossreactive"), "Remove cross-reactive probes (uploaded list)", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_crossreactive")),
                  fileInput(ns("crossreactive_file"), "Probe exclusion list (one ID per line)", accept = c(".csv", ".txt", ".tsv"))),
                checkboxInput(ns("f_maf"), "Minor allele frequency (MAF) filter (uploaded list)", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_maf")),
                  fileInput(ns("maf_file"), "Probe MAF table (probe_id, maf columns)", accept = c(".csv", ".txt", ".tsv")),
                  numericInput(ns("maf_max"), "Maximum allowed MAF", value = 0.05, min = 0, max = 0.5, step = 0.01)),
                checkboxInput(ns("f_variance"), "Variance filter", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_variance")),
                  numericInput(ns("variance_min"), "Minimum variance", value = 0, min = 0, step = 0.001)),
                checkboxInput(ns("f_sd"), "Standard deviation filter (low-variance CpGs)", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_sd")),
                  numericInput(ns("sd_min"), "Minimum standard deviation", value = 0, min = 0, step = 0.01)),
                checkboxInput(ns("f_meanrange"), "Mean value range filter", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_meanrange")),
                  fluidRow(
                    column(6, numericInput(ns("mean_lo"), "Min", value = if (scale_is_beta) 0.01 else -6, step = 0.01)),
                    column(6, numericInput(ns("mean_hi"), "Max", value = if (scale_is_beta) 0.99 else 6, step = 0.01))
                  ))
              )
            ),
            actionButton(ns("run_probe_qc_btn"), "Run Probe QC", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("probe_qc_gate"))
      )
    })

    probe_qc_result <- eventReactive(input$run_probe_qc_btn, {
      sg <- current_subgroup()
      validate(need(length(sg$included) >= 3, sprintf(
        "Only %d sample(s) in this stratum (%s) - need at least 3 for probe QC.", length(sg$included), sg$label)))
      mat <- sg$mat
      anno <- anno_result()

      crossreactive_ids <- NULL
      if (isTRUE(input$f_crossreactive) && !is.null(input$crossreactive_file)) {
        pl <- methyl_parse_probe_list(input$crossreactive_file$datapath, input$crossreactive_file$name)
        if (isTRUE(pl$ok)) crossreactive_ids <- pl$ids
      }
      maf_table <- NULL
      if (isTRUE(input$f_maf) && !is.null(input$maf_file)) {
        pl <- methyl_parse_maf_list(input$maf_file$datapath, input$maf_file$name)
        if (isTRUE(pl$ok)) maf_table <- pl$maf
      }

      withProgress(message = "Running Probe QC", value = 0, {
        filters <- list()
        incProgress(0.3, detail = "Applying filters")
        if (isTRUE(input$f_detp)) filters$detection_p <- methyl_filter_detection_p(mat, methyl_dataset$detp, input$detp_thresh %||% 0.01)
        if (isTRUE(input$f_beadcount)) filters$beadcount <- methyl_filter_beadcount(mat, methyl_dataset$beadcount, input$beadcount_thresh %||% 3)
        if (isTRUE(input$f_snp) && isTRUE(is_illumina_array())) filters$snp <- methyl_filter_snp(mat, anno)
        if (isTRUE(input$f_noncpg) && isTRUE(is_illumina_array())) filters$non_cpg <- methyl_filter_non_cpg(mat)
        if (!identical(input$sexchr_mode %||% "keep", "keep") && isTRUE(is_illumina_array()))
          filters$sex_chr <- methyl_filter_sex_chr(mat, anno, mode = input$sexchr_mode)
        if (isTRUE(input$f_crossreactive)) filters$cross_reactive <- methyl_filter_cross_reactive(mat, crossreactive_ids)
        if (isTRUE(input$f_maf)) filters$maf <- methyl_filter_maf(mat, maf_table, input$maf_max %||% 0.05)
        if (isTRUE(input$f_missing)) filters$missing <- methyl_filter_missing(mat, input$missing_max %||% 0)
        if (isTRUE(input$f_variance)) filters$variance <- methyl_filter_variance(mat, input$variance_min %||% 0)
        if (isTRUE(input$f_sd)) filters$sd <- methyl_filter_sd(mat, input$sd_min %||% 0)
        if (isTRUE(input$f_meanrange)) filters$mean_range <- methyl_filter_mean_range(mat, input$mean_lo %||% 0, input$mean_hi %||% 1)

        keep <- rep(TRUE, nrow(mat))
        for (f in filters) keep <- keep & f$keep
        filtered <- mat[keep, , drop = FALSE]

        incProgress(0.5, detail = "Retention cascade")
        cascade <- methyl_probe_retention_cascade(nrow(mat), filters)

        settings <- list(
          f_detp = input$f_detp, detp_thresh = input$detp_thresh,
          f_beadcount = input$f_beadcount, beadcount_thresh = input$beadcount_thresh,
          f_snp = input$f_snp, f_noncpg = input$f_noncpg, sexchr_mode = input$sexchr_mode %||% "keep",
          f_crossreactive = input$f_crossreactive,
          f_maf = input$f_maf, maf_max = input$maf_max,
          f_missing = input$f_missing, missing_max = input$missing_max,
          f_variance = input$f_variance, variance_min = input$variance_min,
          f_sd = input$f_sd, sd_min = input$sd_min,
          f_meanrange = input$f_meanrange, mean_lo = input$mean_lo, mean_hi = input$mean_hi,
          stratum_label = sg$label
        )
      })

      list(mat = mat, filters = filters, keep = keep, filtered = filtered, cascade = cascade,
           subgroup = sg, anno = anno, settings = settings, run_at = Sys.time())
    })
    observeEvent(probe_qc_result(), probe_qc_has_run(TRUE))

    output$probe_qc_result_ui <- renderUI({
      req(probe_qc_result())
      r <- probe_qc_result()
      tagList(
        run_info_line(r$run_at, sprintf("Stratum: %s.", r$subgroup$label)),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Probe filtering summary"),
            fluidRow(
              valueBox(format(nrow(r$mat), big.mark = ","), "Probes in", icon = icon("dna"), color = "light-blue", width = 3),
              valueBox(format(nrow(r$filtered), big.mark = ","), "Probes kept", icon = icon("check"), color = "green", width = 3),
              valueBox(format(nrow(r$mat) - nrow(r$filtered), big.mark = ","), "Probes removed", icon = icon("xmark"), color = if (nrow(r$mat) > nrow(r$filtered)) "yellow" else "green", width = 3),
              valueBox(ncol(r$mat), "Samples", icon = icon("vial"), color = "purple", width = 3)
            ),
            if (length(r$filters) > 0) DT::dataTableOutput(ns("filter_table")) else p(class = "empty-note", icon("circle-info"), "No probe filters were selected.")
        ),
        div(class = "card",
            div(class = "card-title", icon("filter-circle-dollar"), "Probe retention flowchart", svg_download_link("cascade_plot")),
            p(class = "submodule-desc", "Probes remaining after each selected filter is applied in sequence (cumulative, not independent)."),
            withSpinner(plotly::plotlyOutput(ns("cascade_plot"), height = 320), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("download"), "Download filtered matrix"),
            downloadButton(ns("download_filtered_probeqc"), "Filtered matrix (CSV)", class = "btn-default btn-sm")
        )
      )
    })

    output$filter_table <- DT::renderDataTable({
      req(probe_qc_result())
      r <- probe_qc_result()
      df <- data.frame(
        filter = names(r$filters),
        probes_removed = vapply(r$filters, function(f) sum(!f$keep), integer(1)),
        note = vapply(r$filters, `[[`, character(1), "note")
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$cascade_plot <- plotly::renderPlotly({
      req(probe_qc_result())
      plotly::ggplotly(methyl_plot_cascade(probe_qc_result()$cascade))
    })
    wire_svg_download("cascade_plot", "probe_filtering_cascade.svg")

    output$download_filtered_probeqc <- downloadHandler(
      filename = function() sprintf("methylomics_qc_filtered_%s.csv", methyl_dataset$array_type %||% "matrix"),
      content = function(file) {
        req(probe_qc_result())
        m <- probe_qc_result()$filtered
        utils::write.csv(data.frame(probe_id = rownames(m), m, check.names = FALSE), file, row.names = FALSE)
      }
    )

    output$sex_qc_ui <- renderUI({
      req(methyl_dataset$beta)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("venus-mars"), "Sex check"),
            p(class = "submodule-desc", "Predicts each sample's sex from its methylation profile (minfi::getSex() when raw IDAT is available, else a chrX/chrY beta heuristic) and compares it against the sample sheet's reported sex, if any. Independent of every other tab - nothing here changes what any other method does."),
            actionButton(ns("run_sex_qc_btn"), "Run Sex QC", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("sex_qc_gate"))
      )
    })

    sex_qc_result <- eventReactive(input$run_sex_qc_btn, {
      sg <- current_subgroup()
      validate(need(length(sg$included) >= 1, "No samples in the current stratum."))
      mat <- sg$mat
      anno <- anno_result()
      rg <- current_rg_subset()
      reported_sex <- {
        sheet <- methyl_dataset$sample_sheet
        if (is.null(sheet)) NULL else {
          sample_ids <- methyl_sheet_sample_ids(sheet, colnames(methyl_dataset$beta))
          sex_col <- intersect(METHYL_SEX_COL_CANDIDATES, colnames(sheet))
          if (length(sex_col) == 0) NULL else stats::setNames(as.character(sheet[[sex_col[1]]]), sample_ids)
        }
      }
      sex <- methyl_sex_check(mat, anno, rg, reported_sex)
      list(sex = sex, subgroup = sg, run_at = Sys.time())
    })
    observeEvent(sex_qc_result(), sex_qc_has_run(TRUE))

    output$sex_qc_result_ui <- renderUI({
      req(sex_qc_result())
      r <- sex_qc_result()
      tagList(
        run_info_line(r$run_at, sprintf("Stratum: %s.", r$subgroup$label)),
        div(class = "card",
            div(class = "card-title", icon("venus-mars"), "Sex check results"),
            if (isTRUE(r$sex$ok)) tagList(
              p(class = "empty-note", icon("circle-info"), r$sex$method),
              if (!is.na(r$sex$n_mismatch)) {
                if (r$sex$n_mismatch > 0) p(class = "empty-note", icon("triangle-exclamation"),
                  sprintf("%d of %d sample(s) have a predicted sex that disagrees with their reported sex - possible mislabeled samples.", r$sex$n_mismatch, nrow(r$sex$detail)))
                else p(class = "empty-note", icon("circle-check"), "Predicted sex agrees with reported sex for every sample.")
              } else p(class = "empty-note", icon("circle-info"),
                       "No sex/gender column found in the sample sheet - showing predicted sex only, with nothing to compare it against."),
              DT::dataTableOutput(ns("sex_table"))
            ) else tagList(
              p(class = "empty-note", icon("triangle-exclamation"), r$sex$reason),
              if (isTRUE(methyl_dataset$preloaded)) p(class = "empty-note", icon("circle-info"),
                "The preloaded matrix already had sex-chromosome probes removed - upload your own data (with sex-chromosome probes retained) to run this check live.")
            )
        ),
        if (isTRUE(r$sex$ok)) div(class = "card",
            div(class = "card-title", icon("chart-line"), "X vs. Y methylation", svg_download_link("sex_scatter")),
            withSpinner(plotly::plotlyOutput(ns("sex_scatter"), height = 380), color = "#2563EB", type = 6)
        ),
        if (isTRUE(r$sex$ok) && !is.na(r$sex$n_mismatch) && r$sex$n_mismatch > 0) div(class = "card",
            div(class = "card-title", icon("triangle-exclamation"), "Discordant samples"),
            p(class = "submodule-desc", "Select rows and exclude them via the manual-exclusion mechanism (same one used on the Sample QC tab) if these are believed to be mislabeled."),
            DT::dataTableOutput(ns("discordant_table")),
            actionButton(ns("sex_exclude_btn"), "Exclude selected discordant samples", icon = icon("user-slash"), class = "btn-warning btn-sm", style = "margin-top:8px;")
        )
      )
    })

    output$sex_table <- DT::renderDataTable({
      req(sex_qc_result())
      DT::datatable(sex_qc_result()$sex$detail, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })

    output$sex_scatter <- plotly::renderPlotly({
      req(sex_qc_result())
      d <- sex_qc_result()$sex$detail
      xcol <- intersect(c("mean_chrX", "chrX_median"), colnames(d))[1]
      ycol <- intersect(c("mean_chrY", "chrY_median"), colnames(d))[1]
      validate(need(!is.na(xcol) && !is.na(ycol), "X/Y methylation values are not available for this sex-check method."))
      d$mismatch_label <- if ("sex_mismatch" %in% colnames(d)) ifelse(d$sex_mismatch, "Mismatch", "Concordant") else "Predicted"
      gg <- ggplot(d, aes(x = .data[[xcol]], y = .data[[ycol]], color = predicted_sex, shape = mismatch_label, text = sample)) +
        geom_point(size = 2.5, alpha = 0.85) +
        scale_color_manual(values = arthomix_pair(d$predicted_sex)) +
        labs(x = "Mean/median chrX methylation", y = "Mean/median chrY methylation", color = "Predicted sex", shape = NULL) +
        theme_arthomix()
      plotly::ggplotly(gg, tooltip = "text")
    })
    wire_svg_download("sex_scatter", "sex_check_x_vs_y.svg")

    output$discordant_table <- DT::renderDataTable({
      req(sex_qc_result())
      d <- sex_qc_result()$sex$detail
      d <- d[isTRUE(d$sex_mismatch) | (!is.na(d$sex_mismatch) & d$sex_mismatch), , drop = FALSE]
      DT::datatable(d, rownames = FALSE, selection = "multiple", options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })

    observeEvent(input$sex_exclude_btn, {
      req(sex_qc_result())
      d <- sex_qc_result()$sex$detail
      d <- d[isTRUE(d$sex_mismatch) | (!is.na(d$sex_mismatch) & d$sex_mismatch), , drop = FALSE]
      sel <- d$sample[input$discordant_table_rows_selected %||% integer(0)]
      if (length(sel) > 0) manual_exclude(union(manual_exclude(), sel))
    })

    output$batch_ui <- renderUI({
      req(methyl_dataset$beta)
      sheet <- methyl_dataset$sample_sheet
      batch_cols <- methyl_batch_columns(sheet)
      ruvm_ok <- isTRUE(ruvm_available())
      if (length(batch_cols) == 0 && !ruvm_ok) {
        return(div(class = "card",
          div(class = "card-title", icon("wand-magic-sparkles"), "Batch QC"),
          p(class = "empty-note", icon("triangle-exclamation"), "No batch variable detected in uploaded metadata."),
          p(class = "submodule-desc", "ComBat needs a batch/chip/plate/slide column in the sample sheet. RUVm instead needs raw IDAT input (control probes) plus the missMethyl package - neither is available for this dataset.")
        ))
      }
      tagList(
        div(class = "card",
            div(class = "card-title", icon("wand-magic-sparkles"), "Batch QC"),
            if (length(batch_cols) == 0) p(class = "empty-note", icon("triangle-exclamation"), "No batch variable detected in uploaded metadata - ComBat is unavailable below; RUVm doesn't need one."),
            selectInput(ns("batch_method"), "Method",
                        choices = c(if (length(batch_cols) > 0) c("ComBat" = "combat"), if (ruvm_ok) c("RUVm" = "ruvm")),
                        selected = if (length(batch_cols) > 0) "combat" else "ruvm"),
            conditionalPanel(condition = sprintf("input['%s'] == 'combat'", ns("batch_method")),
              if (length(batch_cols) > 0) selectInput(ns("batch_col"), "Batch / chip column", choices = batch_cols, selected = batch_cols[1])),
            conditionalPanel(condition = sprintf("input['%s'] == 'ruvm'", ns("batch_method")),
              if (!is.null(sheet)) selectInput(ns("ruvm_group_col"), "Factor of interest to protect (e.g. disease/control)", choices = colnames(sheet), selected = colnames(sheet)[1])
              else p(class = "empty-note", icon("triangle-exclamation"), "No sample sheet - RUVm needs a factor-of-interest column."),
              numericInput(ns("ruvm_k"), "Unwanted-variation factors (k)", value = 1, min = 1, max = 10, step = 1),
              p(class = "submodule-desc", "RUVm (Maksimovic et al. 2015) estimates unwanted variation from the array's own internal negative-control probes, conditioned on the factor above - it doesn't need a batch label the way ComBat does.")),
            actionButton(ns("run_batch_btn"), "Run Batch QC", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("batch_qc_gate"))
      )
    })

    batch_qc_result <- eventReactive(input$run_batch_btn, {
      validate(need(input$batch_method %in% c("combat", "ruvm"), "Choose a method before running Batch QC."))
      sg <- current_subgroup()
      validate(need(length(sg$included) >= 1, "No samples in the current stratum."))
      mat <- sg$mat
      sheet <- methyl_dataset$sample_sheet
      before <- methyl_pca_scores(mat)

      scale <- methyl_dataset$input_scale %||% "beta"
      if (identical(input$batch_method, "combat")) {
        sample_ids <- methyl_sheet_sample_ids(sheet, colnames(mat))
        batch <- stats::setNames(as.character(sheet[[input$batch_col]]), sample_ids)[colnames(mat)]
        validate(need(!any(is.na(batch)), "Some samples in this run have no value in the chosen batch column."))
        out <- methyl_batch_correct_combat(mat, batch, input_scale = scale)
      } else {
        sample_ids <- methyl_sheet_sample_ids(sheet, colnames(mat))
        batch <- stats::setNames(as.character(sheet[[input$ruvm_group_col]]), sample_ids)[colnames(mat)]
        out <- methyl_batch_correct_ruvm(mat, current_rg_subset(), batch, k = input$ruvm_k %||% 1, input_scale = scale)
      }
      after <- if (isTRUE(out$ok)) methyl_pca_scores(if (!is.null(out$corrected)) out$corrected else mat) else NULL
      list(out = out, before = before, after = after, batch = batch, method = input$batch_method, subgroup = sg, run_at = Sys.time())
    })
    observeEvent(batch_qc_result(), batch_qc_has_run(TRUE))

    output$batch_result_ui <- renderUI({
      req(batch_qc_result())
      cr <- batch_qc_result()
      if (!isTRUE(cr$out$ok)) return(tagList(
        run_info_line(cr$run_at),
        div(class = "card", p(class = "empty-note", icon("triangle-exclamation"), cr$out$reason))
      ))
      tagList(
        run_info_line(cr$run_at, sprintf("Method: %s. Stratum: %s.", if (identical(cr$method, "ruvm")) "RUVm" else "ComBat", cr$subgroup$label)),
        div(class = "card",
            div(class = "card-title", icon("chart-line"), sprintf("PCA before vs. after %s", if (identical(cr$method, "ruvm")) "RUVm" else "ComBat")),
            fluidRow(
              column(6, h5("Before", svg_download_link("batch_pca_before")), withSpinner(plotly::plotlyOutput(ns("batch_pca_before"), height = 340), color = "#2563EB", type = 6)),
              column(6, h5("After", svg_download_link("batch_pca_after")), withSpinner(plotly::plotlyOutput(ns("batch_pca_after"), height = 340), color = "#2563EB", type = 6))
            )
        ),
        div(class = "card",
            div(class = "card-title", icon("table"), "Variance explained (top PCs)"),
            DT::dataTableOutput(ns("batch_variance_table"))
        )
      )
    })

    .batch_pca_plot <- function(scores_obj, batch, color_lab = "Batch") {
      validate(need(isTRUE(scores_obj$ok), if (isTRUE(scores_obj$ok)) "" else scores_obj$reason))
      df <- data.frame(sample = rownames(scores_obj$scores), PC1 = scores_obj$scores[, 1], PC2 = scores_obj$scores[, 2],
                        batch = as.character(batch)[rownames(scores_obj$scores)])
      gg <- ggplot(df, aes(x = PC1, y = PC2, color = batch, text = sample)) +
        geom_point(size = 2.5, alpha = 0.85) +
        scale_color_manual(values = arthomix_pair(df$batch)) +
        labs(color = color_lab) + theme_arthomix()
      plotly::ggplotly(gg, tooltip = "text")
    }
    output$batch_pca_before <- plotly::renderPlotly({
      cr <- batch_qc_result()
      .batch_pca_plot(cr$before, cr$batch, if (identical(cr$method, "ruvm")) "Factor of interest" else "Batch")
    })
    output$batch_pca_after  <- plotly::renderPlotly({
      cr <- batch_qc_result()
      .batch_pca_plot(cr$after, cr$batch, if (identical(cr$method, "ruvm")) "Factor of interest" else "Batch")
    })
    wire_svg_download("batch_pca_before", "batch_correction_pca_before.svg")
    wire_svg_download("batch_pca_after", "batch_correction_pca_after.svg")

    output$batch_variance_table <- DT::renderDataTable({
      cr <- batch_qc_result()
      df <- data.frame(
        PC = paste0("PC", seq_along(cr$before$var_explained)),
        variance_explained_before = round(cr$before$var_explained * 100, 2),
        variance_explained_after = round((cr$after$var_explained[seq_along(cr$before$var_explained)] %||% NA), 2)
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$outlier_ui <- renderUI({
      req(methyl_dataset$beta)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("magnifying-glass"), "Outlier detection settings"),
            fluidRow(
              column(6,
                checkboxGroupInput(ns("outlier_methods"), "Outlier detection method(s)",
                                   choices = c("PCA-based detection" = "pca", "Hierarchical clustering" = "hclust",
                                               "Correlation-based detection" = "correlation", "Distance-based detection (Mahalanobis)" = "mahalanobis"),
                                   selected = c("pca", "hclust"))
              ),
              column(6,
                conditionalPanel(condition = sprintf("input['%s'].includes('pca')", ns("outlier_methods")),
                  numericInput(ns("pca_sd"), "PCA distance threshold (SD)", value = 3, min = 1, step = 0.5)),
                conditionalPanel(condition = sprintf("input['%s'].includes('hclust')", ns("outlier_methods")),
                  numericInput(ns("hclust_height_frac"), "Hierarchical-clustering singleton height fraction", value = 0.5, min = 0.05, max = 1, step = 0.05)),
                conditionalPanel(condition = sprintf("input['%s'].includes('correlation')", ns("outlier_methods")),
                  numericInput(ns("corr_k"), "Correlation MAD multiplier (k)", value = 3, min = 1, step = 0.5)),
                conditionalPanel(condition = sprintf("input['%s'].includes('mahalanobis')", ns("outlier_methods")),
                  numericInput(ns("mahal_alpha"), "Mahalanobis distance alpha", value = 0.01, min = 0.001, max = 0.2, step = 0.001))
              )
            ),
            actionButton(ns("run_outlier_btn"), "Run Outlier Detection", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("outlier_qc_gate"))
      )
    })

    outlier_qc_result <- eventReactive(input$run_outlier_btn, {
      sg <- current_subgroup()
      validate(need(length(sg$included) >= 4, sprintf(
        "Only %d sample(s) in this stratum - need at least 4 for outlier detection.", length(sg$included))))
      mat <- sg$mat
      methods_sel <- input$outlier_methods %||% character(0)
      validate(need(length(methods_sel) > 0, "Select at least one outlier-detection method."))

      not_selected <- function(label) list(ok = FALSE, reason = sprintf("%s was not selected for this run.", label))

      withProgress(message = "Running Outlier Detection", value = 0.2, {
        sample_qc <- data.frame(sample = colnames(mat), stringsAsFactors = FALSE)
        pca_detail <- if ("pca" %in% methods_sel) methyl_sample_outliers_pca(mat, sd_threshold = input$pca_sd %||% 3) else not_selected("PCA-based detection")
        if ("pca" %in% methods_sel) sample_qc$pca_outlier <- if (isTRUE(pca_detail$ok)) pca_detail$outlier[colnames(mat)] else NA

        incProgress(0.3, detail = "Clustering / correlation / distance")
        hc_res <- if ("hclust" %in% methods_sel) methyl_sample_outliers_hclust(mat, height_frac = input$hclust_height_frac %||% 0.5) else not_selected("Hierarchical clustering")
        if ("hclust" %in% methods_sel) sample_qc$hclust_outlier <- if (isTRUE(hc_res$ok)) hc_res$outlier[colnames(mat)] else NA

        corr_detail <- if ("correlation" %in% methods_sel) methyl_sample_outliers_correlation(mat, k = input$corr_k %||% 3) else not_selected("Correlation-based detection")
        if ("correlation" %in% methods_sel) sample_qc$correlation_outlier <- if (isTRUE(corr_detail$ok)) corr_detail$outlier[colnames(mat)] else NA

        mahal_detail <- if ("mahalanobis" %in% methods_sel) methyl_sample_outliers_mahalanobis(mat, alpha = input$mahal_alpha %||% 0.01) else not_selected("Distance-based detection (Mahalanobis)")
        if ("mahalanobis" %in% methods_sel) sample_qc$mahalanobis_outlier <- if (isTRUE(mahal_detail$ok)) mahal_detail$outlier[colnames(mat)] else NA

        outlier_scores <- methyl_outlier_score_table(sample_qc)
      })

      list(mat = mat, subgroup = sg, sample_qc = sample_qc, outlier_scores = outlier_scores,
           pca_detail = pca_detail, mahal_detail = mahal_detail, corr_detail = corr_detail, hc_res = hc_res,
           methods = methods_sel, run_at = Sys.time())
    })
    observeEvent(outlier_qc_result(), outlier_qc_has_run(TRUE))

    output$outlier_result_ui <- renderUI({
      req(outlier_qc_result())
      r <- outlier_qc_result()
      tagList(
        run_info_line(r$run_at, sprintf("Method(s): %s. Stratum: %s.", paste(r$methods, collapse = ", "), r$subgroup$label)),
        div(class = "card",
            div(class = "card-title", icon("ranking-star"), "Outlier score table"),
            p(class = "submodule-desc", "Number of the selected methods that flagged each sample. Nothing is excluded automatically - select rows below and click \"Apply Sample Exclusions\" to remove them from the sample scope every QC method reads from."),
            DT::dataTableOutput(ns("outlier_score_table")),
            div(style = "margin-top:8px;", actionButton(ns("apply_outlier_exclusions_btn"), "Apply Sample Exclusions", icon = icon("user-slash"), class = "btn-warning btn-sm")),
            uiOutput(ns("manual_exclude_summary"))
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-line"), "PCA sample structure"),
            if (isTRUE(r$pca_detail$ok))
              lazy_plot_ui("pca_outlier_plot", height = 380, label = "Generate PCA plot")
            else p(class = "empty-note", icon("triangle-exclamation"), r$pca_detail$reason)
        ),
        div(class = "card",
            div(class = "card-title", icon("ranking-star"), "Outlier diagnostic: PCA distance from centroid"),
            if (isTRUE(r$pca_detail$ok))
              lazy_plot_ui("outlier_diagnostic_plot", height = 340, label = "Generate diagnostic plot")
            else p(class = "empty-note", icon("triangle-exclamation"), r$pca_detail$reason)
        ),
        div(class = "card",
            div(class = "card-title", icon("code-branch"), "Dendrogram"),
            p(class = "submodule-desc", "Shown when \"Hierarchical clustering\" was selected for this run."),
            if ("hclust" %in% r$methods)
              lazy_plot_ui("dendro_plot", height = 420, label = "Generate dendrogram", with_svg = FALSE, kind = "plot")
            else p(class = "empty-note", icon("circle-info"), "Select \"Hierarchical clustering\" and re-run to see the dendrogram.")
        )
      )
    })

    output$outlier_score_table <- DT::renderDataTable({
      req(outlier_qc_result())
      DT::datatable(outlier_qc_result()$outlier_scores, rownames = FALSE, selection = "multiple",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })

    observeEvent(input$apply_outlier_exclusions_btn, {
      req(outlier_qc_result())
      df <- outlier_qc_result()$outlier_scores
      sel <- df$sample[input$outlier_score_table_rows_selected %||% integer(0)]
      if (length(sel) > 0) manual_exclude(union(manual_exclude(), sel))
    })

    output$pca_outlier_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["pca_outlier_plot"]]))
      req(outlier_qc_result())
      p <- outlier_qc_result()$pca_detail
      df <- data.frame(sample = colnames(outlier_qc_result()$mat), PC1 = p$scores[, 1], PC2 = p$scores[, 2],
                        outlier = ifelse(p$outlier, "Flagged", "Within range"))
      gg <- ggplot(df, aes(x = PC1, y = PC2, color = outlier, text = sample)) +
        geom_point(size = 2.5, alpha = 0.85) +
        scale_color_manual(values = c("Within range" = ARTHOMIX_COLORS$blue, "Flagged" = ARTHOMIX_STATUS$critical)) +
        labs(color = NULL) + theme_arthomix()
      plotly::ggplotly(gg, tooltip = "text")
    })
    wire_svg_download("pca_outlier_plot", "outlier_pca_scatter.svg")

    output$outlier_diagnostic_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["outlier_diagnostic_plot"]]))
      req(outlier_qc_result())
      r <- outlier_qc_result()
      diag <- methyl_outlier_diagnostic_table(r$sample_qc, r$pca_detail, r$mahal_detail)
      validate(need("pca_distance" %in% colnames(diag), "PCA distance is not available for this run."))
      plotly::ggplotly(methyl_plot_outlier_diagnostic(diag, "pca_distance", "PCA distance from centroid"))
    })
    wire_svg_download("outlier_diagnostic_plot", "outlier_diagnostic.svg")

    output$dendro_plot <- renderPlot({
      req(isTRUE(plot_shown[["dendro_plot"]]))
      req(outlier_qc_result())
      hc_res <- outlier_qc_result()$hc_res
      validate(need(isTRUE(hc_res$ok), if (isTRUE(hc_res$ok)) "" else hc_res$reason))
      dd <- ggdendro::dendro_data(hc_res$hc)
      labs_df <- dd$labels
      labs_df$outlier <- hc_res$outlier[as.character(labs_df$label)]
      ggplot() +
        geom_segment(data = dd$segments, aes(x = x, y = y, xend = xend, yend = yend), color = ARTHOMIX_COLORS$axis) +
        geom_text(data = labs_df, aes(x = x, y = 0, label = label, color = outlier), angle = 90, hjust = 1, size = 2.6) +
        scale_color_manual(values = c(`TRUE` = ARTHOMIX_STATUS$critical, `FALSE` = ARTHOMIX_COLORS$ink_muted), guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0.35, 0.05))) +
        labs(x = NULL, y = "Height") + theme_arthomix() +
        theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
    })

    output$viz_ui <- renderUI({
      req(methyl_dataset$beta)
      has_ctrl <- !is.null(current_rg_subset())
      has_detp <- !is.null(methyl_dataset$detp)
      has_beadcount <- !is.null(methyl_dataset$beadcount)

      tagList(
        uiOutput(ns("viz_probe_qc_gate")),
        div(class = "card",
            div(class = "card-title", icon("chart-area"), "Detection p-value distribution (raw IDAT only)"),
            if (has_detp) lazy_plot_ui("detp_heatmap", height = 320, label = "Generate detection p-value plot")
            else p(class = "empty-note", icon("circle-info"), "Requires raw IDAT input - upload IDAT files on the Dataset tab.")
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Bead count distribution (raw IDAT only)"),
            if (has_beadcount) lazy_plot_ui("beadcount_dist", height = 280, label = "Generate bead count plot")
            else p(class = "empty-note", icon("circle-info"), "Requires raw IDAT input - upload IDAT files on the Dataset tab.")
        ),
        div(class = "card",
            div(class = "card-title", icon("border-all"), "Control-probe heatmap (raw IDAT only)"),
            if (has_ctrl) lazy_plot_ui("control_heatmap", height = 420, label = "Generate control-probe heatmap")
            else p(class = "empty-note", icon("circle-info"), "Requires raw IDAT input - control-probe intensities aren't present in a beta/M-value matrix upload or the preloaded dataset's bundled matrix.")
        )
      )
    })

    output$viz_probe_qc_gate <- renderUI({
      if (!isTRUE(probe_qc_has_run())) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "No figures to show yet"),
          p(class = "empty-note", icon("circle-info"),
            "PCA, MDS, density, boxplot, violin, correlation-heatmap, and mean-SD plots below all need Probe QC's filtered matrix - nothing here has run yet."),
          actionButton(ns("viz_goto_probeqc_btn"), "Go to Probe QC", icon = icon("arrow-right"), class = "btn-primary btn-sm")
        ))
      }
      sheet <- methyl_dataset$sample_sheet
      color_choices <- c("Sample (no coloring)" = "__none__")
      if (!is.null(sheet)) color_choices <- c(color_choices, stats::setNames(colnames(sheet), colnames(sheet)))
      tagList(
        div(class = "card", div(class = "card-title", icon("palette"), "Color by (PCA/MDS)"),
            selectInput(ns("viz_color_by"), NULL, choices = color_choices, selected = "__none__", width = "300px")),
        div(class = "card", div(class = "card-title", icon("chart-line"), "PCA (2D)"),
            lazy_plot_ui("viz_pca_plot", height = 400, label = "Generate PCA (2D)")),
        div(class = "card", div(class = "card-title", icon("cube"), "PCA (3D, PC1-PC3)"),
            p(class = "submodule-desc", "3D scatter is rendered via WebGL - drag to rotate, scroll to zoom. SVG export isn't supported for 3D plots (PNG only, via the camera icon)."),
            lazy_plot_ui("viz_pca3d_plot", height = 460, label = "Generate PCA (3D)", with_svg = FALSE)),
        div(class = "card", div(class = "card-title", icon("braille"), "MDS"),
            p(class = "submodule-desc", "Classical multidimensional scaling on pairwise sample distances - a distinct view of sample structure from PCA."),
            lazy_plot_ui("mds_plot", height = 380, label = "Generate MDS")),
        div(class = "card", div(class = "card-title", icon("chart-area"), "Beta-value density: before vs. after Probe QC"),
            lazy_plot_ui("density_plot", height = 340, label = "Generate density plot")),
        div(class = "card", div(class = "card-title", icon("chart-column"), "Beta-value boxplot (per sample)"),
            lazy_plot_ui("boxplot_plot", height = 380, label = "Generate boxplot")),
        div(class = "card", div(class = "card-title", icon("chart-simple"), "Beta-value violin plot (per sample)"),
            lazy_plot_ui("violin_plot", height = 380, label = "Generate violin plot")),
        div(class = "card", div(class = "card-title", icon("table-cells"), "Sample correlation heatmap"),
            p(class = "submodule-desc", "Top 5,000 most-variable probes; samples ordered by hierarchical clustering."),
            lazy_plot_ui("corr_heatmap", height = 420, label = "Generate correlation heatmap")),
        div(class = "card", div(class = "card-title", icon("chart-line"), "Mean-SD plot"),
            p(class = "submodule-desc", "Per-probe mean vs. standard deviation - a flat trend line means variance isn't confounded with mean methylation."),
            lazy_plot_ui("mean_sd_plot", height = 360, label = "Generate mean-SD plot"))
      )
    })

    observeEvent(input$viz_goto_probeqc_btn, {
      updateTabsetPanel(session, "tabs", selected = "Probe QC")
    })

    .viz_color_vec <- function(mat) {
      cb <- input$viz_color_by %||% "__none__"
      if (identical(cb, "__none__") || is.null(methyl_dataset$sample_sheet)) {
        list(vals = stats::setNames(rep("Sample", ncol(mat)), colnames(mat)), pal = c("Sample" = ARTHOMIX_COLORS$blue))
      } else {
        sheet <- methyl_dataset$sample_sheet
        sample_ids <- methyl_sheet_sample_ids(sheet, colnames(mat))
        vals <- stats::setNames(as.character(sheet[[cb]]), sample_ids)
        list(vals = vals, pal = arthomix_pair(stats::na.omit(vals)))
      }
    }

    output$viz_pca_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["viz_pca_plot"]]))
      req(probe_qc_result())
      r <- probe_qc_result()
      p <- methyl_pca_scores(r$filtered)
      validate(need(isTRUE(p$ok), if (isTRUE(p$ok)) "" else p$reason))
      cv <- .viz_color_vec(r$filtered)
      df <- data.frame(x = p$scores[, 1], y = p$scores[, 2], color = cv$vals[rownames(p$scores)], text = rownames(p$scores))
      plotly::ggplotly(methyl_plot_scatter2d(df, "PC1", "PC2", NULL, cv$pal), tooltip = "text")
    })
    wire_svg_download("viz_pca_plot", "pca_2d.svg")

    output$viz_pca3d_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["viz_pca3d_plot"]]))
      req(probe_qc_result())
      r <- probe_qc_result()
      p <- methyl_pca_scores(r$filtered, n_pcs = 10)
      validate(need(isTRUE(p$ok) && ncol(p$scores) >= 3, "Fewer than 3 principal components are available for this run."))
      cv <- .viz_color_vec(r$filtered)
      df <- data.frame(PC1 = p$scores[, 1], PC2 = p$scores[, 2], PC3 = p$scores[, 3],
                        color = cv$vals[rownames(p$scores)], sample = rownames(p$scores))
      plotly::plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~color, colors = unname(cv$pal),
                       text = ~sample, type = "scatter3d", mode = "markers", marker = list(size = 4)) %>%
        plotly::layout(scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")))
    })

    output$mds_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["mds_plot"]]))
      req(probe_qc_result())
      m <- methyl_mds_scores(probe_qc_result()$filtered)
      validate(need(isTRUE(m$ok), if (isTRUE(m$ok)) "" else m$reason))
      cv <- .viz_color_vec(probe_qc_result()$filtered)
      df <- data.frame(x = m$scores[, 1], y = m$scores[, 2], color = cv$vals[rownames(m$scores)], text = rownames(m$scores))
      plotly::ggplotly(methyl_plot_scatter2d(df, "Dim1", "Dim2", NULL, cv$pal), tooltip = "text")
    })
    wire_svg_download("mds_plot", "mds.svg")

    viz_density_df <- reactive({
      req(probe_qc_result())
      r <- probe_qc_result()
      probe_ids <- rownames(r$mat)[order(stats::runif(nrow(r$mat)))[seq_len(min(5000, nrow(r$mat)))]]
      before <- methyl_beta_density_sample(r$mat, seed_probes = probe_ids)
      before$stage <- "Before filtering"
      after_ids <- intersect(probe_ids, rownames(r$filtered))
      if (length(after_ids) > 0) {
        after <- methyl_beta_density_sample(r$filtered, seed_probes = after_ids); after$stage <- "After filtering"
        rbind(before, after)
      } else before
    })

    output$density_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["density_plot"]]))
      x_lab <- if (identical(methyl_dataset$input_scale, "beta")) "Beta value" else "M-value"
      plotly::ggplotly(methyl_plot_density(viz_density_df(), x_lab))
    })
    wire_svg_download("density_plot", "beta_density.svg")

    output$boxplot_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["boxplot_plot"]]))
      plotly::ggplotly(methyl_plot_boxplot(viz_density_df()))
    })
    wire_svg_download("boxplot_plot", "beta_boxplot.svg")

    output$violin_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["violin_plot"]]))
      plotly::ggplotly(methyl_plot_violin(viz_density_df()))
    })
    wire_svg_download("violin_plot", "beta_violin.svg")

    output$corr_heatmap <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["corr_heatmap"]]))
      req(probe_qc_result())
      cr <- methyl_sample_correlation(probe_qc_result()$filtered)
      validate(need(isTRUE(cr$ok), if (isTRUE(cr$ok)) "" else cr$reason))
      plotly::ggplotly(methyl_plot_corr_heatmap(cr$cor))
    })
    wire_svg_download("corr_heatmap", "sample_correlation_heatmap.svg")

    output$mean_sd_plot <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["mean_sd_plot"]]))
      req(probe_qc_result())
      plotly::ggplotly(methyl_plot_mean_sd(methyl_mean_sd_table(probe_qc_result()$filtered)))
    })
    wire_svg_download("mean_sd_plot", "mean_sd_plot.svg")

    output$detp_heatmap <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["detp_heatmap"]]))
      req(methyl_dataset$detp)
      mat <- current_subgroup()$mat
      common <- intersect(rownames(methyl_dataset$detp), rownames(mat))
      validate(need(length(common) > 0, "No overlapping probes between the detection p-value matrix and the loaded matrix."))
      probes <- common[order(stats::runif(length(common)))[seq_len(min(150, length(common)))]]
      sub <- methyl_dataset$detp[probes, colnames(mat), drop = FALSE]
      df <- as.data.frame(as.table(sub)); colnames(df) <- c("probe", "sample", "detection_p")
      plotly::ggplotly(methyl_plot_detp_heatmap(df))
    })
    wire_svg_download("detp_heatmap", "detection_pvalue_distribution.svg")

    output$beadcount_dist <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["beadcount_dist"]]))
      req(methyl_dataset$beadcount)
      mat <- current_subgroup()$mat
      common <- intersect(rownames(methyl_dataset$beadcount), rownames(mat))
      validate(need(length(common) > 0, "No overlapping probes between the bead-count matrix and the loaded matrix."))
      vals <- as.numeric(methyl_dataset$beadcount[common, colnames(mat), drop = FALSE])
      plotly::ggplotly(methyl_plot_beadcount_dist(vals, input$beadcount_thresh %||% 3))
    })
    wire_svg_download("beadcount_dist", "bead_count_distribution.svg")

    output$control_heatmap <- plotly::renderPlotly({
      req(isTRUE(plot_shown[["control_heatmap"]]))
      ctrl <- methyl_control_probe_matrix(current_rg_subset())
      validate(need(isTRUE(ctrl$ok), if (isTRUE(ctrl$ok)) "" else ctrl$reason))
      m <- ctrl$mat
      if (nrow(m) > 60) m <- m[order(methyl_row_vars(m), decreasing = TRUE)[seq_len(60)], , drop = FALSE]
      df <- as.data.frame(as.table(m)); colnames(df) <- c("control_probe", "sample", "log2_intensity")
      gg <- ggplot(df, aes(x = sample, y = control_probe, fill = log2_intensity)) +
        geom_tile() +
        scale_fill_gradient(low = ARTHOMIX_COLORS$blue, high = ARTHOMIX_COLORS$orange) +
        labs(x = NULL, y = NULL, fill = "log2 intensity") + theme_arthomix() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), axis.text.y = element_text(size = 5))
      plotly::ggplotly(gg)
    })
    wire_svg_download("control_heatmap", "control_probe_heatmap.svg")

    pdf_report_available <- reactive({
      has_pandoc <- requireNamespace("rmarkdown", quietly = TRUE) &&
        isTRUE(tryCatch(rmarkdown::pandoc_available(), error = function(e) FALSE))
      has_latex <- any(nzchar(Sys.which(c("pdflatex", "xelatex", "lualatex")))) ||
        (requireNamespace("tinytex", quietly = TRUE) && isTRUE(tryCatch(tinytex::is_tinytex(), error = function(e) FALSE)))
      has_pandoc && has_latex
    })

    current_qc_pieces <- reactive(list(
      overview = if (isTRUE(overview_has_run())) overview_result() else NULL,
      sample_qc = if (isTRUE(sample_qc_has_run())) sample_qc_result() else NULL,
      probe_qc = if (isTRUE(probe_qc_has_run())) probe_qc_result() else NULL,
      sex_qc = if (isTRUE(sex_qc_has_run())) sex_qc_result() else NULL,
      outlier_qc = if (isTRUE(outlier_qc_has_run())) outlier_qc_result() else NULL,
      batch_qc = if (isTRUE(batch_qc_has_run())) batch_qc_result() else NULL
    ))

    output$reports_ui <- renderUI({
      req(methyl_dataset$beta)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("file-export"), "Downloads"),
            p(class = "submodule-desc", "Each export uses whichever QC tabs have actually been run this session - run a tab first if its export below is unavailable."),
            if (isTRUE(probe_qc_has_run())) tagList(
              downloadButton(ns("dl_filtered_beta"), "Filtered beta matrix (CSV)", class = "btn-default btn-sm"),
              downloadButton(ns("dl_filtered_mvalue"), "Filtered M-value matrix (CSV)", class = "btn-default btn-sm", style = "margin-left:8px;"),
              downloadButton(ns("dl_filter_summary"), "Probe filter summary (CSV)", class = "btn-default btn-sm", style = "margin-left:8px;")
            ) else p(class = "empty-note", icon("circle-info"), "Filtered matrix / probe filter summary: run Probe QC first."),
            tags$div(style = "margin-top:10px;",
              if (isTRUE(sample_qc_has_run())) downloadButton(ns("dl_sample_qc"), "Sample QC table (CSV)", class = "btn-default btn-sm")
              else p(class = "empty-note", icon("circle-info"), "Sample QC table: run Sample QC first.")
            ),
            tags$div(style = "margin-top:10px;",
              if (isTRUE(batch_qc_has_run()) && isTRUE(batch_qc_result()$out$ok))
                downloadButton(ns("dl_batch_corrected"), "Batch-corrected matrix (CSV)", class = "btn-default btn-sm")
              else p(class = "empty-note", icon("circle-info"), "Batch-corrected matrix: run Batch QC first.")
            ),
            tags$div(style = "margin-top:10px;",
              downloadButton(ns("dl_qc_summary"), "QC summary (CSV)", class = "btn-default btn-sm")
            )
        ),
        div(class = "card",
            div(class = "card-title", icon("file-lines"), "QC report & figures"),
            p(class = "submodule-desc", "A self-contained report (every figure baked in as an image, opens standalone in any browser) plus every figure as a separate PNG - built from whichever tabs have been run; sections not yet run are listed as skipped."),
            downloadButton(ns("dl_report_html"), "QC report (HTML)", class = "btn-primary btn-sm"),
            if (isTRUE(pdf_report_available())) downloadButton(ns("dl_report_pdf"), "QC report (PDF)", class = "btn-default btn-sm", style = "margin-left:8px;")
            else p(class = "empty-note", style = "margin-top:6px;", icon("circle-info"),
                   "PDF export needs a LaTeX toolchain (rmarkdown + tinytex) that isn't installed in this deployment - the HTML report above opens in any browser and can be printed to PDF from there."),
            downloadButton(ns("dl_figures_zip"), "All figures (ZIP)", class = "btn-default btn-sm", style = "margin-left:8px;")
        ),
        if (isTRUE(probe_qc_has_run())) div(class = "card",
            div(class = "card-title", icon("clipboard-list"), "Probe QC reproducibility"),
            p(class = "submodule-desc", sprintf("Exact Probe QC settings used - stratum: %s.", probe_qc_result()$settings$stratum_label)),
            DT::dataTableOutput(ns("repro_table"))
        ),
        if (isTRUE(probe_qc_has_run())) div(class = "card",
            div(class = "card-title", icon("code"), "Copy R code (Probe QC)"),
            p(class = "submodule-desc", "An equivalent minfi/ChAMP-style Bioconductor snippet for the Probe QC filters above."),
            actionButton(ns("copy_code_btn"), "Copy to clipboard", icon = icon("copy"), class = "btn-default btn-sm"),
            tags$pre(id = ns("r_code_pre"), style = "margin-top:8px; max-height:320px; overflow:auto;", verbatimTextOutput(ns("r_code")))
        )
      )
    })

    output$dl_filtered_beta <- downloadHandler(
      filename = function() "methylomics_qc_filtered_beta.csv",
      content = function(file) {
        req(probe_qc_result())
        m <- probe_qc_result()$filtered
        if (!identical(methyl_dataset$input_scale, "beta")) m <- methyl_mvalue_to_beta(m)
        utils::write.csv(data.frame(probe_id = rownames(m), m, check.names = FALSE), file, row.names = FALSE)
      }
    )
    output$dl_filtered_mvalue <- downloadHandler(
      filename = function() "methylomics_qc_filtered_mvalue.csv",
      content = function(file) {
        req(probe_qc_result())
        m <- probe_qc_result()$filtered
        if (identical(methyl_dataset$input_scale, "beta")) m <- methyl_beta_to_mvalue(m)
        utils::write.csv(data.frame(probe_id = rownames(m), m, check.names = FALSE), file, row.names = FALSE)
      }
    )
    output$dl_filter_summary <- downloadHandler(
      filename = function() "methylomics_probe_filter_summary.csv",
      content = function(file) {
        req(probe_qc_result())
        r <- probe_qc_result()
        df <- data.frame(filter = names(r$filters),
                          probes_removed = vapply(r$filters, function(f) sum(!f$keep), integer(1)),
                          note = vapply(r$filters, `[[`, character(1), "note"))
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_sample_qc <- downloadHandler(
      filename = function() "methylomics_sample_qc.csv",
      content = function(file) {
        req(sample_qc_result())
        utils::write.csv(sample_qc_result()$sample_qc, file, row.names = FALSE)
      }
    )
    output$dl_batch_corrected <- downloadHandler(
      filename = function() sprintf("methylomics_qc_batch_corrected_%s.csv", batch_qc_result()$method %||% "corrected"),
      content = function(file) {
        req(batch_qc_result())
        r <- batch_qc_result()
        validate(need(isTRUE(r$out$ok), "Batch correction did not complete successfully for this run."))
        m <- r$out$corrected
        utils::write.csv(data.frame(probe_id = rownames(m), m, check.names = FALSE), file, row.names = FALSE)
      }
    )
    output$dl_qc_summary <- downloadHandler(
      filename = function() "methylomics_qc_summary.csv",
      content = function(file) {
        p <- current_qc_pieces()
        df <- methyl_qc_summary_table(overview = p$overview, sample_qc = p$sample_qc, probe_qc = p$probe_qc,
                                       sex_qc = p$sex_qc, outlier_qc = p$outlier_qc, batch_qc = p$batch_qc)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    report_plots <- reactive({
      p <- current_qc_pieces()
      methyl_qc_report_plots(methyl_dataset, probe_qc = p$probe_qc, outlier_qc = p$outlier_qc)
    })

    output$dl_report_html <- downloadHandler(
      filename = function() "methylomics_qc_report.html",
      content = function(file) {
        rp <- report_plots()
        p <- current_qc_pieces()
        summary_df <- methyl_qc_summary_table(overview = p$overview, sample_qc = p$sample_qc, probe_qc = p$probe_qc,
                                               sex_qc = p$sex_qc, outlier_qc = p$outlier_qc, batch_qc = p$batch_qc)
        subtitle <- sprintf("Source: %s.", methyl_dataset$source %||% "n/a")
        out <- methyl_qc_report_html(methyl_dataset, summary_df, rp$plots, rp$skipped, subtitle)
        validate(need(isTRUE(out$ok), if (isTRUE(out$ok)) "" else out$reason))
        file.copy(out$path, file, overwrite = TRUE)
      }
    )

    output$dl_report_pdf <- downloadHandler(
      filename = function() "methylomics_qc_report.pdf",
      content = function(file) {
        validate(need(isTRUE(pdf_report_available()), "PDF export requires rmarkdown + a LaTeX toolchain, not available in this deployment."))
        rp <- report_plots()
        p <- current_qc_pieces()
        summary_df <- methyl_qc_summary_table(overview = p$overview, sample_qc = p$sample_qc, probe_qc = p$probe_qc,
                                               sex_qc = p$sex_qc, outlier_qc = p$outlier_qc, batch_qc = p$batch_qc)
        rmd <- tempfile(fileext = ".Rmd")
        fig_paths <- vapply(names(rp$plots), function(nm) {
          f <- tempfile(fileext = ".png")
          ggplot2::ggsave(f, rp$plots[[nm]], width = 8, height = 5, dpi = 110)
          f
        }, character(1))
        lines <- c(
          "---", "title: \"Methylomics QC report\"", "output: pdf_document", "---",
          sprintf("Source: %s.", methyl_rmd_safe_text(methyl_dataset$source %||% "n/a")),
          "", "## Summary", "",
          "| Metric | Value |", "|---|---|",
          sprintf("| %s | %s |", methyl_rmd_safe_text(summary_df$metric), methyl_rmd_safe_text(summary_df$value)),
          "", "## Figures", "",
          unlist(lapply(seq_along(fig_paths), function(i) c(sprintf("### %s", methyl_rmd_safe_text(names(fig_paths)[i])), sprintf("![](%s)", fig_paths[i]), "")))
        )
        writeLines(lines, rmd)
        out <- tryCatch(rmarkdown::render(rmd, output_file = file, output_format = "pdf_document", quiet = TRUE, envir = new.env()), error = function(e) e)
        unlink(fig_paths)
        validate(need(!inherits(out, "error"), paste("PDF rendering failed:", if (inherits(out, "error")) conditionMessage(out) else "")))
      }
    )

    output$dl_figures_zip <- downloadHandler(
      filename = function() "methylomics_qc_figures.zip",
      content = function(file) {
        rp <- report_plots()
        out <- methyl_qc_report_zip(rp$plots)
        validate(need(isTRUE(out$ok), if (isTRUE(out$ok)) "" else out$reason))
        file.copy(out$path, file, overwrite = TRUE)
      }
    )

    output$repro_table <- DT::renderDataTable({
      req(probe_qc_result())
      s <- probe_qc_result()$settings
      flat <- s[!vapply(s, is.list, logical(1))]
      flat <- flat[!names(flat) %in% "stratum_label"]
      df <- data.frame(setting = names(flat), value = vapply(flat, function(v) paste(v, collapse = ", "), character(1)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 30), class = "stripe hover compact")
    })

    output$r_code <- renderText({
      req(probe_qc_result())
      methyl_qc_r_code(probe_qc_result()$settings)
    })

    observeEvent(input$copy_code_btn, {
      shinyjs::runjs(sprintf(
        "navigator.clipboard.writeText(document.getElementById('%s').innerText).catch(function(){});",
        ns("r_code")
      ))
    })
  })
}
