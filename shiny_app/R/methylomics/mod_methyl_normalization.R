## R/methylomics/mod_methyl_normalization.R
## Methylomics sub-module: Normalization. Reads the shared `methyl_dataset`
## reactiveValues (R/methylomics/normalization.R has the actual method
## implementations, diagnostics, status-detection, and validation helpers -
## see its own header for what each one calls and why).
##
## ---- Two data pathways, one shared diagnostics/status engine -------------
## Both the preloaded dataset and an uploaded one run through the SAME
## automatic diagnostics (methyl_norm_diagnostics()) and normalization-status
## detection (methyl_norm_status()) below - there is no separate "preloaded
## logic" for that part. What differs is what happens next: the preloaded
## whole-blood dataset was analyzed from its original author-normalized data
## (see the Dataset tab), a deliberate choice that preserves direct
## comparability with the originally published beta values rather than
## risking systematic differences from an independently re-implemented
## normalization - so for preloaded data this tab stops at "here is what the
## data looks like and here is its status", it never offers a live method
## picker, filters, or a Run button that would reprocess it. An uploaded
## dataset gets the full live workflow: filters, a method picker whose
## choices depend on what was actually uploaded (raw IDAT unlocks Noob/
## Functional normalization/SWAN/Dasen/Stratified quantile/Noob+SWAN, all of
## which need an RGChannelSet; BMIQ/PBC/Noob+BMIQ need Type I/II manifest
## annotation - 450K/EPIC only - but the matrix-based ones work from a beta
## matrix alone; plain quantile normalization and "no normalization" always
## work), a Compare Methods panel, and progressively-revealed before/after
## results. Nothing is guessed or offered when its input isn't actually
## available.
##
## Follows the same "compute a candidate result, then an explicit button
## promotes it to the shared dataset" pattern mod_preprocessing.R already
## uses for transcriptomics ("Use this as the active dataset") - the run
## button alone never overwrites methyl_dataset$beta, so a normalization
## can be tried, inspected via the before/after diagnostics, and either
## kept or discarded.

mod_methyl_normalization_config <- list(
  id = "normalization", title = "Normalization", icon = "wave-square", group = "Data",
  description = "Automatic data diagnostics and normalization-status detection, methylation-specific filters, and Noob, Functional normalization, SWAN, Dasen, Stratified quantile normalization (IDAT only), BMIQ and PBC (450K/EPIC beta values), Noob+BMIQ / Noob+SWAN sequential workflows, and plain quantile normalization (any input) - with method comparison and before/after QC."
)

METHYL_NORM_METHODS_IDAT <- c(
  "Noob" = "noob", "Functional normalization" = "funnorm", "SWAN" = "swan",
  "Dasen" = "dasen", "Stratified quantile normalization" = "stratified_quantile"
)
METHYL_NORM_METHODS_MATRIX <- c("BMIQ" = "bmiq", "PBC" = "pbc")
METHYL_NORM_METHODS_UNIVERSAL <- c("No normalization / keep current data" = "none", "Quantile normalization (plain)" = "quantile")
METHYL_NORM_METHODS_COMBO_SWAN <- c("Noob + SWAN" = "noob_swan")
METHYL_NORM_METHODS_COMBO_BMIQ <- c("Noob + BMIQ" = "noob_bmiq")

## Methods that operate directly on a beta/M-value matrix (filters are
## applied BEFORE the method runs); every other method runs on raw IDAT
## input and filters are applied to its resulting beta AFTER the method
## runs, since minfi/wateRmelon's raw-intensity methods work on the array's
## full probe set from the RGChannelSet, not a pre-filtered matrix.
METHYL_NORM_MATRIX_METHOD_KEYS <- c("bmiq", "pbc", "quantile", "none")

mod_methyl_normalization_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    withSpinner(uiOutput(ns("body_ui")), color = "#2563EB", type = 6)
  )
}

## Small helper: an actionLink that toggles a hidden wrapper div - the same
## "nothing renders until expanded" affordance mod_methyl_qc.R's historical-
## reference toggle already uses, applied here to every progressive-reveal
## section in the results area.
.methyl_norm_toggle_ui <- function(ns, id, label) {
  tagList(
    actionLink(ns(paste0(id, "_toggle")),
               label = tagList(icon("chevron-down"), " ", label),
               class = "empty-note", style = "cursor:pointer; font-weight:600;"),
    shinyjs::hidden(div(id = ns(paste0(id, "_wrap")), style = "margin-top:8px; margin-bottom:12px;",
      uiOutput(ns(paste0(id, "_body")))
    ))
  )
}

## Sub-tab title (icon + label) - same convention as mod_methyl_qc.R's own
## qc_tab_title() for its Overview/Sample QC/Probe QC/... tabs.
norm_tab_title <- function(ic, label) tagList(icon(ic), " ", label)

mod_methyl_normalization_server <- function(id, methyl_dataset, methyl_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ==== Shared: annotation, diagnostics, status (both pathways) =========

    anno_result <- reactive({
      req(methyl_dataset$array_type)
      methyl_get_annotation(methyl_dataset$array_type)
    })
    norm_anno_result <- reactive({
      req(methyl_dataset$array_type)
      methyl_get_norm_annotation(methyl_dataset$array_type)
    })
    manifest_available <- reactive(isTRUE(anno_result()$ok))
    has_idat <- reactive(!is.null(methyl_dataset$rg_set))
    is_beta_scale <- reactive(identical(methyl_dataset$input_scale, "beta"))
    is_illumina <- reactive({
      !is.null(methyl_dataset$array_type) && methyl_dataset$array_type %in% METHYL_ARRAY_TYPES_ILLUMINA
    })

    diag_result <- reactive({
      req(methyl_dataset$beta)
      methyl_norm_diagnostics(methyl_dataset$beta, methyl_dataset, anno_result())
    })
    status_result <- reactive({
      req(methyl_dataset$beta)
      methyl_norm_status(methyl_dataset$beta, methyl_dataset, anno_result())
    })

    diagnostics_card <- function(d, extra_note = NULL) {
      tagList(
        div(class = "card",
            div(class = "card-title", icon("magnifying-glass-chart"), "Automatic data diagnostics"),
            if (!is.null(extra_note)) p(class = "submodule-desc", extra_note),
            fluidRow(
              valueBox(format(d$n_probes, big.mark = ","), "Probes / features", icon = icon("dna"), color = "aqua", width = 3),
              valueBox(d$n_samples, "Samples", icon = icon("vial"), color = "purple", width = 3),
              valueBox(sprintf("%.2f%%", d$pct_missing), "Missing values", icon = icon("circle-question"), color = if (d$pct_missing > 5) "yellow" else "green", width = 3),
              valueBox(d$representation_short, "Representation", icon = icon("ruler"), color = "light-blue", width = 3)
            ),
            .methyl_norm_toggle_ui(ns, "diag_details", "Show details")
        )
      )
    }

    output$diag_details_body <- renderUI({
      req(methyl_dataset$beta)
      d <- diag_result()
      tagList(
        tags$table(class = "table table-condensed", style = "width:100%; font-size:13px;",
          tags$tbody(
            tags$tr(tags$td("Representation"), tags$td(d$representation)),
            tags$tr(tags$td("Value range"), tags$td(sprintf("%.4f - %.4f", d$value_min, d$value_max))),
            tags$tr(tags$td("Median"), tags$td(sprintf("%.4f", d$median))),
            tags$tr(tags$td("Mean"), tags$td(sprintf("%.4f", d$mean))),
            tags$tr(tags$td("Standard deviation"), tags$td(sprintf("%.4f", d$sd))),
            tags$tr(tags$td("Infinite values"), tags$td(format(d$n_infinite, big.mark = ","))),
            tags$tr(tags$td("Detected platform"), tags$td(d$platform)),
            tags$tr(tags$td("Raw methylated/unmethylated intensities available"), tags$td(if (d$has_raw_intensity) "Yes" else "No")),
            tags$tr(tags$td("Beta values available"), tags$td(if (d$has_beta) "Yes" else "No")),
            tags$tr(tags$td("M-values available"), tags$td(if (d$has_mvalue) "Yes" else "No")),
            tags$tr(tags$td("Type I/II probe-design annotation available"), tags$td(if (d$has_probe_type_info) "Yes" else "No"))
          )
        )
      )
    })

    status_card <- function(st, extra, title = "Normalization status") {
      badge <- switch(st$status,
        raw = list(color = "yellow", icon = "triangle-exclamation", label = "Raw, unprocessed IDAT data"),
        bias_detected = list(color = "yellow", icon = "triangle-exclamation", label = "Type I/II probe-design bias detected"),
        no_bias_detected = list(color = "green", icon = "circle-check", label = "No Type I/II probe-design bias detected"),
        list(color = "light-blue", icon = "circle-question", label = "Unable to determine automatically"))
      tagList(
        div(class = "card",
            div(class = "card-title", icon(badge$icon), title),
            fluidRow(valueBox(badge$label, "Status", icon = icon(badge$icon), color = badge$color, width = 12)),
            p(class = "submodule-desc", st$message),
            extra
        )
      )
    }

    ## ==== Branch: no dataset loaded at all ==================================

    output$body_ui <- renderUI({
      if (isTRUE(methyl_dataset$preloaded)) return(preloaded_ui())
      if (is.null(methyl_dataset$beta)) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "No methylation dataset loaded"),
          p(class = "submodule-desc", "Upload a beta/M-value matrix or raw IDAT files on the Methylomics Dataset tab first.")
        ))
      }
      uploaded_ui()
    })

    ## ==== Branch: preloaded dataset - diagnostics/status only, no live =====
    ## reprocessing (see this file's header comment for why).

    preloaded_ui <- function() {
      if (is.null(methyl_dataset$beta)) {
        return(div(class = "card",
          div(class = "card-title", icon("circle-info"), "No separate normalization step for this dataset"),
          p(class = "submodule-desc",
            "The preloaded whole-blood dataset's live beta matrix isn't available in this deployment, so there is nothing to run diagnostics against. The dataset was analyzed from its original author-normalized data - see the Dataset tab.")
        ))
      }
      d <- diag_result(); st <- status_result()
      ## Deliberately NOT status_card(st, ...): that helper's badge reflects
      ## methyl_norm_status()'s narrow Type I/II probe-design-bias signal
      ## (see normalization.R's header on that function), which is a
      ## DIFFERENT question from "was this dataset normalized" - for this
      ## specific cohort the documented, ground-truth answer to the latter
      ## is already known (author-normalized before GEO deposition; see the
      ## Dataset tab) and is not something a heuristic should second-guess
      ## or appear to contradict. The heuristic's actual reading is still
      ## shown, just demoted to a secondary technical note rather than the
      ## headline badge, so it stays informative without reading as "this
      ## dataset needs normalizing" for a dataset this tab will never offer
      ## to reprocess anyway.
      tagList(
        diagnostics_card(d, "Computed automatically from the preloaded whole-blood cohort's bundled beta matrix."),
        div(class = "card",
            div(class = "card-title", icon("circle-check"), "Normalization status"),
            fluidRow(valueBox("Using existing (author) normalization", "Status", icon = icon("circle-check"), color = "green", width = 12)),
            p(class = "submodule-desc",
              "This is the app's own reference dataset, analyzed from its original author-normalized data (Liu et al. 2013, deposited on GEO) rather than reprocessed here - a deliberate choice that preserves direct comparability with the originally published beta values (see the Dataset tab). \"Re-normalize\" and \"Compare normalization methods\" are intentionally not offered for this dataset."),
            if (!is.null(st$bias)) div(class = "empty-note", icon("flask"), sprintf(
              "Technical note: this tab's automatic Type I/II probe-design-bias check (a narrower signal than \"was this normalized at all\" - see the Filters/Method sections' own explanations below for the distinction) reads %s (Kolmogorov-Smirnov statistic = %.3f) for this dataset. GEO-deposited processed beta values from this era typically received background/dye/batch correction but not necessarily a probe-design-aware step (BMIQ/SWAN); either way, this dataset's own author-normalized values are kept as-is rather than reprocessed here.",
              if (identical(st$status, "bias_detected")) "substantial residual Type I/II difference" else if (identical(st$status, "no_bias_detected")) "little residual Type I/II difference" else "an inconclusive result",
              st$bias$ks_stat))
        )
      )
    }

    ## ==== Branch: uploaded dataset - full live workflow ====================

    available_methods <- reactive({
      c(
        if (isTRUE(has_idat())) METHYL_NORM_METHODS_IDAT,
        if (isTRUE(has_idat())) METHYL_NORM_METHODS_COMBO_SWAN,
        if (isTRUE(has_idat()) && isTRUE(manifest_available())) METHYL_NORM_METHODS_COMBO_BMIQ,
        if (isTRUE(manifest_available()) && isTRUE(is_beta_scale())) METHYL_NORM_METHODS_MATRIX,
        METHYL_NORM_METHODS_UNIVERSAL
      )
    })

    default_method <- reactive({
      am <- available_methods(); st <- status_result()
      if (identical(st$status, "no_bias_detected") && "none" %in% am) return("none")
      if ("noob" %in% am) return("noob")
      if ("bmiq" %in% am) return("bmiq")
      unname(utils::tail(am, 1))
    })

    recommendation_text <- reactive({
      methyl_norm_recommendation(methyl_dataset, status_result(), available_methods())
    })

    ## Section "Already-normalized data": when the status card says
    ## "no_bias_detected", the filters/method-picker/compare-methods sections
    ## below stay hidden until the user actively picks "Re-normalize" or
    ## "Compare normalization methods" - "Keep current normalization" just
    ## confirms the status quo (the currently loaded matrix already IS the
    ## active dataset; nothing here needs to change it). Reset whenever a
    ## genuinely new dataset is loaded.
    norm_choice <- reactiveVal(NULL)
    observeEvent(methyl_dataset$beta, norm_choice(NULL), ignoreNULL = TRUE)
    observeEvent(input$choice_keep_btn, norm_choice("keep"), ignoreInit = TRUE)
    observeEvent(input$choice_renorm_btn, norm_choice("renormalize"), ignoreInit = TRUE)
    observeEvent(input$choice_compare_btn, norm_choice("compare"), ignoreInit = TRUE)
    observeEvent(input$choice_reset_btn, norm_choice(NULL), ignoreInit = TRUE)

    show_live_workflow <- reactive({
      st <- status_result()
      if (!identical(st$status, "no_bias_detected")) return(TRUE)
      identical(norm_choice(), "renormalize") || identical(norm_choice(), "compare")
    })

    uploaded_ui <- function() {
      d <- diag_result(); st <- status_result()

      already_normalized_choice_ui <- if (identical(st$status, "no_bias_detected")) {
        if (is.null(norm_choice())) {
          div(style = "display:flex; gap:8px; flex-wrap:wrap; margin-top:8px;",
              actionButton(ns("choice_keep_btn"), "Keep current normalization", icon = icon("check"), class = "btn-primary btn-sm"),
              actionButton(ns("choice_renorm_btn"), "Re-normalize", icon = icon("rotate"), class = "btn-default btn-sm"),
              actionButton(ns("choice_compare_btn"), "Compare normalization methods", icon = icon("scale-balanced"), class = "btn-default btn-sm")
          )
        } else if (identical(norm_choice(), "keep")) {
          tagList(
            div(class = "empty-note", icon("circle-check"),
                "No changes made - the currently loaded matrix remains the active Methylomics dataset for every sub-module below."),
            actionLink(ns("choice_reset_btn"), "Change this choice", style = "font-size:12px;")
          )
        } else {
          actionLink(ns("choice_reset_btn"), tagList(icon("arrow-left"), " Back to normalization-status choices"), style = "font-size:12px;")
        }
      } else NULL

      tagList(
        diagnostics_card(d),
        status_card(st, already_normalized_choice_ui),
        if (isTRUE(show_live_workflow())) tagList(
          div(class = "card",
              div(class = "card-title", icon("lightbulb"), "Recommended approach"),
              p(class = "submodule-desc", recommendation_text()),
              p(class = "empty-note", icon("circle-info"), "This is advisory only - every compatible method below remains selectable regardless of this recommendation.")
          ),
          ## Filters / Method & Run / Compare Methods as sub-tabs (same
          ## tabsetPanel convention mod_methyl_qc.R already uses for its own
          ## Overview/Sample QC/Probe QC/... tabs) rather than one long
          ## vertical stack of cards - a tabsetPanel renders every tabPanel's
          ## content up front (tab switching is CSS visibility, not
          ## conditional rendering), so filters set on the Filters tab stay
          ## live and are read by both Method & Run and Compare Methods
          ## regardless of which tab is currently open.
          div(class = "tx-menu-wrap",
              tabsetPanel(
                id = ns("config_tabs"), type = "tabs", header = tagList(tags$hr()),
                tabPanel(value = "Filters", title = norm_tab_title("filter", "Filters"),
                         br(), filters_ui()),
                tabPanel(value = "Method & Run", title = norm_tab_title("sliders", "Method & Run"),
                         br(), method_ui(), withSpinner(uiOutput(ns("results_ui")), color = "#2563EB", type = 6)),
                tabPanel(value = "Compare Methods", title = norm_tab_title("scale-balanced", "Compare Methods"),
                         br(), compare_methods_ui())
              )
          )
        )
      )
    }

    ## ---- Filters --------------------------------------------------------

    filters_ui <- function() {
      sheet <- methyl_dataset$sample_sheet
      island_choices <- if ("island_relation" %in% colnames(norm_anno_result()$anno %||% data.frame())) {
        sort(unique(stats::na.omit(as.character(norm_anno_result()$anno$island_relation))))
      } else character(0)
      region_choices <- if ("gene_region" %in% colnames(norm_anno_result()$anno %||% data.frame())) {
        sort(unique(stats::na.omit(as.character(norm_anno_result()$anno$gene_region))))
      } else character(0)
      chr_choices <- if (isTRUE(manifest_available())) sort(unique(stats::na.omit(as.character(anno_result()$anno$chr)))) else character(0)

      div(class = "card",
          div(class = "card-title", icon("filter"), "Filters"),
          p(class = "submodule-desc", "Every filter below is optional and only applied when you click Run Normalization or Compare Methods - nothing here changes the loaded dataset by itself."),
          h5("General"),
          fluidRow(
            column(6,
              checkboxInput(ns("f_probe_missing"), "Probe missingness threshold", value = FALSE),
              conditionalPanel(condition = sprintf("input['%s']", ns("f_probe_missing")),
                numericInput(ns("probe_missing_max"), "Max missing fraction per probe", value = 0.05, min = 0, max = 1, step = 0.01)),
              checkboxInput(ns("f_sample_missing"), "Sample missingness threshold", value = FALSE),
              conditionalPanel(condition = sprintf("input['%s']", ns("f_sample_missing")),
                numericInput(ns("sample_missing_max"), "Max missing fraction per sample", value = 0.1, min = 0, max = 1, step = 0.01))
            ),
            column(6,
              checkboxInput(ns("f_sample_detrate"), "Minimum sample detection rate", value = FALSE),
              conditionalPanel(condition = sprintf("input['%s']", ns("f_sample_detrate")),
                fluidRow(
                  column(6, numericInput(ns("detp_thresh"), "Detection p-value threshold", value = 0.01, min = 0, max = 1, step = 0.01)),
                  column(6, numericInput(ns("min_det_rate"), "Min sample detection rate (%)", value = 95, min = 0, max = 100, step = 1))
                )),
              if (!isTRUE(has_idat())) p(class = "empty-note", icon("circle-info"), "Detection-rate filtering requires raw IDAT input - unavailable for this dataset.")
            )
          ),
          if (isTRUE(is_illumina())) tagList(
            h5("Methylation-specific"),
            fluidRow(
              column(6,
                checkboxInput(ns("f_snp"), "Remove SNP-overlapping probes", value = FALSE),
                checkboxInput(ns("f_crossreactive"), "Remove cross-reactive probes (uploaded list)", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_crossreactive")),
                  fileInput(ns("crossreactive_file"), "Probe exclusion list (one ID per line)", accept = c(".csv", ".txt", ".tsv"))),
                checkboxInput(ns("f_sexchr"), "Remove sex-chromosome (chrX/chrY) probes", value = FALSE),
                checkboxInput(ns("f_chr"), "Exclude specific chromosome(s)", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("f_chr")),
                  selectizeInput(ns("exclude_chr"), NULL, choices = chr_choices, multiple = TRUE))
              ),
              column(6,
                if (length(island_choices) > 0) tagList(
                  checkboxInput(ns("f_island"), "Filter by CpG island relation", value = FALSE),
                  conditionalPanel(condition = sprintf("input['%s']", ns("f_island")),
                    checkboxGroupInput(ns("island_categories"), "Keep probes in:", choices = island_choices, selected = island_choices, inline = TRUE))
                ) else p(class = "empty-note", icon("circle-info"), "CpG island-relation annotation is unavailable for this array type."),
                if (length(region_choices) > 0) tagList(
                  checkboxInput(ns("f_generegion"), "Filter by gene region", value = FALSE),
                  conditionalPanel(condition = sprintf("input['%s']", ns("f_generegion")),
                    checkboxGroupInput(ns("gene_regions"), "Keep probes in:", choices = region_choices, selected = region_choices, inline = TRUE))
                ) else p(class = "empty-note", icon("circle-info"), "Gene-region annotation is unavailable for this array type.")
              )
            )
          ) else p(class = "empty-note", icon("circle-info"), "SNP / cross-reactive / chromosome / island / gene-region filters only apply to Illumina array data and are hidden for this dataset type."),
          if (!is.null(sheet) && ncol(sheet) > 0) tagList(
            h5("Biological-signal-preservation check (optional)"),
            selectInput(ns("group_col_check"), "Group column to check PCA separation against before vs. after",
                        choices = c("(none)" = "", colnames(sheet)), width = "100%")
          )
      )
    }

    cross_reactive_ids <- reactive({
      if (!isTRUE(input$f_crossreactive) || is.null(input$crossreactive_file)) return(NULL)
      pl <- methyl_parse_probe_list(input$crossreactive_file$datapath, input$crossreactive_file$name)
      if (isTRUE(pl$ok)) pl$ids else NULL
    })

    ## Every filter's keep-vector, applied by the caller (matrix methods:
    ## before running; raw-intensity methods: on the resulting beta) - see
    ## this file's header comment on why the order differs by method.
    build_probe_filters <- function(m) {
      filters <- list()
      if (isTRUE(input$f_probe_missing)) filters$probe_missing <- methyl_filter_missing(m, input$probe_missing_max %||% 0.05)
      if (isTRUE(input$f_snp) && isTRUE(is_illumina())) filters$snp <- methyl_filter_snp(m, anno_result())
      if (isTRUE(input$f_sexchr) && isTRUE(is_illumina())) filters$sex_chr <- methyl_filter_sex_chr(m, anno_result(), mode = "remove_xy")
      if (isTRUE(input$f_chr) && isTRUE(is_illumina())) filters$chromosome <- methyl_filter_chromosome(m, anno_result(), input$exclude_chr %||% character(0))
      if (isTRUE(input$f_crossreactive) && isTRUE(is_illumina())) filters$cross_reactive <- methyl_filter_cross_reactive(m, cross_reactive_ids())
      if (isTRUE(input$f_island)) filters$island <- methyl_filter_island_relation(m, norm_anno_result(), input$island_categories %||% character(0))
      if (isTRUE(input$f_generegion)) filters$gene_region <- methyl_filter_gene_region(m, norm_anno_result(), input$gene_regions %||% character(0))
      filters
    }

    ## ---- Method picker ----------------------------------------------------

    cond_any <- function(input_id, values) paste(sprintf("input['%s'] == '%s'", ns(input_id), values), collapse = " || ")

    method_ui <- function() {
      am <- available_methods()
      tagList(
        div(class = "card",
            div(class = "card-title", icon("sliders"), "Normalization method"),
            if (!isTRUE(has_idat())) p(class = "empty-note", icon("circle-info"),
              "No raw IDAT input for this dataset - Noob, Functional normalization, SWAN, Dasen, Stratified quantile normalization, Noob+SWAN, and Noob+BMIQ all require it and are hidden below."),
            if (!isTRUE(manifest_available())) p(class = "empty-note", icon("circle-info"),
              sprintf("No manifest annotation for %s - BMIQ, PBC, and Noob+BMIQ (all need Type I/II probe design) are hidden below.", methyl_dataset$array_type %||% "this array type")),
            if (isTRUE(manifest_available()) && !isTRUE(is_beta_scale())) p(class = "empty-note", icon("circle-info"),
              "Uploaded matrix is on the M-value scale - BMIQ and PBC both require actual beta values (0-1) and are hidden below. Re-upload as beta values, or use plain quantile normalization instead."),
            radioButtons(ns("method"), NULL, choices = am, selected = default_method()),
            uiOutput(ns("method_info_ui")),
            conditionalPanel(condition = cond_any("method", c("noob", "noob_bmiq", "noob_swan")),
              fluidRow(
                column(6, numericInput(ns("noob_offset"), "Background offset", value = 15, min = 0, step = 1)),
                column(6, radioButtons(ns("noob_dye"), "Dye-bias correction", inline = TRUE,
                                       choices = c("Single-sample (ssNoob)" = "single", "Reference array" = "reference"), selected = "single"))
              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'funnorm'", ns("method")),
              numericInput(ns("funnorm_npcs"), "Number of control-probe PCs", value = 2, min = 1, max = 10, step = 1)),
            conditionalPanel(condition = cond_any("method", c("bmiq", "noob_bmiq")),
              fluidRow(
                column(4, numericInput(ns("bmiq_nfit"), "Probes used per type (nfit)", value = 50000, min = 1000, step = 1000)),
                column(4, numericInput(ns("bmiq_nl"), "Number of states/components (nL)", value = 3, min = 2, max = 4, step = 1)),
                column(4, numericInput(ns("bmiq_tol"), "Convergence tolerance", value = 0.001, min = 0.0001, max = 0.01, step = 0.0001))
              )),
            conditionalPanel(condition = cond_any("method", c("swan", "dasen", "stratified_quantile", "pbc", "quantile", "none")),
              p(class = "empty-note", icon("circle-info"), "No additional tunable parameters for this method.")),
            actionButton(ns("run_btn"), "Run Normalization", icon = icon("play"), class = "btn-primary")
        )
      )
    }

    output$method_info_ui <- renderUI({
      req(input$method)
      info <- METHYL_NORM_METHOD_INFO[[input$method]]
      if (is.null(info)) return(NULL)
      div(class = "empty-note", style = "margin-bottom:10px;",
          icon("circle-info"), strong(sprintf(" %s: ", info$category)), info$text)
    })

    ## ---- Compare Methods --------------------------------------------------

    compare_methods_ui <- function() {
      am <- available_methods()
      div(class = "card",
          div(class = "card-title", icon("scale-balanced"), "Compare methods"),
          p(class = "submodule-desc", "Select two or more compatible methods to run side by side, using the same filters configured above. Nothing runs until you click Compare Methods."),
          checkboxGroupInput(ns("compare_methods_sel"), NULL, choices = am, inline = TRUE),
          actionButton(ns("compare_btn"), "Compare Methods", icon = icon("scale-balanced"), class = "btn-default btn-sm"),
          withSpinner(uiOutput(ns("compare_result_ui")), color = "#2563EB", type = 6)
      )
    }

    ## ---- Core: run one method on one (already sample-filtered) input ------
    ## Returns list(ok, beta, note) - shared by both the single Run button
    ## and the Compare Methods panel so both call the exact same code path.
    run_one_method <- function(method_key, mat_in, rg_in, mset_in) {
      switch(method_key,
        noob = methyl_norm_noob(rg_in, input$noob_offset %||% 15, input$noob_dye %||% "single"),
        funnorm = methyl_norm_funnorm(rg_in, input$funnorm_npcs %||% 2),
        swan = methyl_norm_swan(rg_in, mset_in),
        dasen = methyl_norm_dasen(mset_in),
        stratified_quantile = methyl_norm_stratified_quantile(rg_in),
        noob_swan = methyl_norm_noob_swan(rg_in, input$noob_offset %||% 15, input$noob_dye %||% "single"),
        noob_bmiq = methyl_norm_noob_bmiq(rg_in, anno_result(), input$noob_offset %||% 15, input$noob_dye %||% "single",
                                           input$bmiq_nfit %||% 50000, input$bmiq_nl %||% 3, input$bmiq_tol %||% 0.001),
        bmiq = methyl_norm_bmiq(mat_in, anno_result(), input$bmiq_nfit %||% 50000, input$bmiq_nl %||% 3, input$bmiq_tol %||% 0.001),
        pbc = methyl_norm_pbc(mat_in, anno_result()),
        quantile = methyl_norm_quantile(mat_in),
        none = list(ok = TRUE, beta = mat_in, note = "No normalization applied - matrix kept exactly as loaded."),
        list(ok = FALSE, reason = "Unknown method.")
      )
    }

    ## Sample-level filters + subsetting, shared by Run and Compare Methods.
    sample_scope <- reactive({
      mat <- methyl_dataset$beta
      keep <- rep(TRUE, ncol(mat)); notes <- character(0)
      if (isTRUE(input$f_sample_missing)) {
        f <- methyl_filter_samples_missingness(mat, input$sample_missing_max %||% 0.1)
        keep <- keep & f$keep; notes <- c(notes, f$note)
      }
      if (isTRUE(input$f_sample_detrate) && !is.null(methyl_dataset$detp)) {
        fp <- methyl_sample_failed_probe_pct(mat, methyl_dataset$detp, input$detp_thresh %||% 0.01)
        if (isTRUE(fp$ok)) {
          det_rate <- 100 - fp$pct[colnames(mat)]
          keep2 <- !is.na(det_rate) & det_rate >= (input$min_det_rate %||% 95)
          keep <- keep & keep2
          notes <- c(notes, sprintf("%d sample(s) below %.0f%% minimum detection rate removed.", sum(!keep2), input$min_det_rate %||% 95))
        }
      }
      list(mat = mat[, keep, drop = FALSE],
           rg = if (!is.null(methyl_dataset$rg_set)) tryCatch(methyl_dataset$rg_set[, keep], error = function(e) NULL) else NULL,
           mset = if (!is.null(methyl_dataset$mset)) tryCatch(methyl_dataset$mset[, keep], error = function(e) NULL) else NULL,
           removed_samples = colnames(mat)[!keep], notes = notes)
    })

    group_labels_for_check <- function(sample_ids) {
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet) || is.null(input$group_col_check) || !nzchar(input$group_col_check)) return(NULL)
      ids <- methyl_sheet_sample_ids(sheet, sample_ids)
      stats::setNames(as.character(sheet[[input$group_col_check]]), ids)[sample_ids]
    }

    ## One end-to-end run: sample filters -> method -> probe filters
    ## (ordered per method type, see this file's header) -> before/after
    ## comparison on the shared probe set. Returns NULL (with a `reason`)
    ## on failure instead of throwing.
    run_full <- function(method_key) {
      sc <- sample_scope()
      if (ncol(sc$mat) < 2) return(list(ok = FALSE, reason = "Fewer than 2 samples remain after the selected sample filters - relax the thresholds above."))
      is_matrix_method <- method_key %in% METHYL_NORM_MATRIX_METHOD_KEYS

      if (is_matrix_method) {
        filters <- build_probe_filters(sc$mat)
        keep <- rep(TRUE, nrow(sc$mat)); for (f in filters) keep <- keep & f$keep
        mat_in <- sc$mat[keep, , drop = FALSE]
        res <- run_one_method(method_key, mat_in, sc$rg, sc$mset)
      } else {
        filters <- list()
        res <- run_one_method(method_key, sc$mat, sc$rg, sc$mset)
        if (isTRUE(res$ok)) {
          filters <- build_probe_filters(res$beta)
          keep <- rep(TRUE, nrow(res$beta)); for (f in filters) keep <- keep & f$keep
          res$beta <- res$beta[keep, , drop = FALSE]
        }
      }
      if (!isTRUE(res$ok)) return(res)
      if (nrow(res$beta) < 1) return(list(ok = FALSE, reason = "All probes were removed by the selected filters - relax them and try again."))

      common_probes <- intersect(rownames(sc$mat), rownames(res$beta))
      before <- sc$mat[common_probes, , drop = FALSE]
      after <- res$beta[common_probes, , drop = FALSE]
      filter_notes <- c(sc$notes, vapply(filters, `[[`, character(1), "note"))
      removed_probes <- setdiff(rownames(sc$mat), common_probes)

      group_labels <- group_labels_for_check(colnames(after))
      validation <- methyl_norm_validation(before, after, anno_result(), group_labels)

      list(ok = TRUE, method = method_key, method_label = names(available_methods())[available_methods() == method_key],
           before = before, after = after, note = res$note,
           removed_probes = removed_probes, removed_samples = sc$removed_samples, filter_notes = filter_notes,
           n_probes_before = nrow(sc$mat), n_samples_before = ncol(sc$mat),
           validation = validation, run_at = Sys.time())
    }

    ## ---- Run Normalization --------------------------------------------

    norm_has_run <- reactiveVal(FALSE)
    observeEvent(methyl_dataset$beta, norm_has_run(FALSE), ignoreNULL = TRUE)

    norm_result <- eventReactive(input$run_btn, {
      req(input$method)
      validate(need(!is.null(methyl_dataset$beta), "Upload a dataset on the Dataset tab first."))
      r <- withProgress(message = sprintf("Running %s", names(available_methods())[available_methods() == input$method]), value = 0.3, {
        run_full(input$method)
      })
      validate(need(isTRUE(r$ok), r$reason %||% "Normalization failed."))
      r
    })

    observeEvent(norm_result(), {
      r <- norm_result()
      methyl_results$normalization <- list(method = r$method_label, note = r$note, n_probes = nrow(r$after), n_samples = ncol(r$after))
      norm_has_run(TRUE)
    })

    observeEvent(input$promote_btn, {
      r <- norm_result()
      methyl_dataset$beta <- r$after
      ## Every method except plain quantile normalization / "no
      ## normalization" (both scale-agnostic per normalization.R's header)
      ## derives or requires actual beta values, so it's safe to stamp
      ## "beta" for those - but quantile-normalizing (or leaving unchanged)
      ## an uploaded M-value matrix does not change its scale, and
      ## relabeling it "beta" here would corrupt input_scale for every
      ## downstream module that trusts it.
      methyl_dataset$input_scale <- if (r$method %in% c("quantile", "none")) methyl_dataset$input_scale else "beta"
      methyl_dataset$source <- sprintf("%s, normalized with %s", methyl_dataset$source %||% "Methylation dataset", r$method_label)
      showNotification(sprintf("The %s-normalized matrix is now the active Methylomics dataset for every sub-module below.", r$method_label), type = "message")
    }, ignoreInit = TRUE)

    output$download_normalized <- downloadHandler(
      filename = function() sprintf("methylomics_normalized_%s.csv", norm_result()$method),
      content = function(file) {
        m <- norm_result()$after
        utils::write.csv(data.frame(probe_id = rownames(m), m, check.names = FALSE), file, row.names = FALSE)
      }
    )
    output$download_removed_probes <- downloadHandler(
      filename = function() sprintf("methylomics_removed_probes_%s.csv", norm_result()$method),
      content = function(file) utils::write.csv(data.frame(probe_id = norm_result()$removed_probes), file, row.names = FALSE)
    )
    output$download_removed_samples <- downloadHandler(
      filename = function() sprintf("methylomics_removed_samples_%s.csv", norm_result()$method),
      content = function(file) utils::write.csv(data.frame(sample_id = norm_result()$removed_samples), file, row.names = FALSE)
    )
    output$download_processing_record <- downloadHandler(
      filename = function() sprintf("methylomics_normalization_record_%s.txt", norm_result()$method),
      content = function(file) {
        r <- norm_result()
        writeLines(methyl_norm_processing_record(list(
          method_label = r$method_label, dataset_source = methyl_dataset$source, run_at = r$run_at,
          n_samples_before = r$n_samples_before, n_samples_after = ncol(r$after),
          n_probes_before = r$n_probes_before, n_probes_after = nrow(r$after),
          filter_notes = r$filter_notes, note = r$note
        )), file)
      }
    )

    ## ---- Progressive results -----------------------------------------

    output$results_ui <- renderUI({
      req(norm_has_run(), norm_result())
      r <- norm_result()
      interp <- methyl_norm_interpretation(r$validation)
      status_color <- c(pass = "#e6f4ea", warning = "#fff4e5", neutral = "#eef2f7")[interp$status]
      tagList(
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Normalization summary"),
            div(style = sprintf("background:%s; padding:10px 14px; border-radius:6px; margin-bottom:10px;", status_color),
                icon(switch(interp$status, pass = "circle-check", warning = "triangle-exclamation", "circle-info")), " ", interp$text),
            fluidRow(
              valueBox(r$method_label, "Method", icon = icon("wave-square"), color = "light-blue", width = 3),
              valueBox(sprintf("%s -> %s", format(r$n_samples_before, big.mark = ","), ncol(r$after)), "Samples", icon = icon("vial"), color = "purple", width = 3),
              valueBox(sprintf("%s -> %s", format(r$n_probes_before, big.mark = ","), format(nrow(r$after), big.mark = ",")), "Probes", icon = icon("dna"), color = "green", width = 3),
              valueBox(sprintf("%.2f%% -> %.2f%%", r$validation$missing_before, r$validation$missing_after), "Missingness", icon = icon("circle-question"), color = "yellow", width = 3)
            ),
            fluidRow(
              valueBox(sprintf("%.4f", stats::median(r$before, na.rm = TRUE)), "Median beta before", icon = icon("chart-line"), color = "light-blue", width = 3),
              valueBox(sprintf("%.4f", stats::median(r$after, na.rm = TRUE)), "Median beta after", icon = icon("chart-line"), color = "green", width = 3),
              valueBox(sprintf("%.4f", mean(r$before, na.rm = TRUE)), "Mean beta before", icon = icon("chart-bar"), color = "light-blue", width = 3),
              valueBox(sprintf("%.4f", mean(r$after, na.rm = TRUE)), "Mean beta after", icon = icon("chart-bar"), color = "green", width = 3)
            ),
            p(class = "empty-note", icon("clock"), sprintf("Run at %s. %d probe(s) and %d sample(s) removed by the selected filters.",
              format(r$run_at, "%Y-%m-%d %H:%M:%S"), length(r$removed_probes), length(r$removed_samples))),
            div(style = "display:flex; gap:8px; flex-wrap:wrap;",
                actionButton(ns("promote_btn"), "Use this as the active Methylomics dataset", icon = icon("check"), class = "btn-primary btn-sm"),
                downloadButton(ns("download_normalized"), "Normalized matrix (CSV)", class = "btn-default btn-sm"),
                downloadButton(ns("download_removed_probes"), "Filtered probe list (CSV)", class = "btn-default btn-sm"),
                downloadButton(ns("download_removed_samples"), "Removed sample list (CSV)", class = "btn-default btn-sm"),
                downloadButton(ns("download_processing_record"), "Processing details (TXT)", class = "btn-default btn-sm")
            )
        ),
        div(class = "card",
            .methyl_norm_toggle_ui(ns, "bvsa", "Show before vs. after"),
            .methyl_norm_toggle_ui(ns, "stats", "Show normalization statistics"),
            .methyl_norm_toggle_ui(ns, "filtered", "Show filtered probes / removed samples"),
            .methyl_norm_toggle_ui(ns, "details", "Show processing details")
        )
      )
    })

    lapply(c("diag_details", "bvsa", "stats", "filtered", "details"), function(sec) {
      observeEvent(input[[paste0(sec, "_toggle")]], shinyjs::toggle(id = paste0(sec, "_wrap")), ignoreInit = TRUE)
    })

    output$bvsa_body <- renderUI({
      req(norm_result())
      tagList(
        p(class = "submodule-desc", "Type I and Type II probes have systematically different beta-value distributions before normalization; a well-normalized dataset shows the two largely converge. Each line/box is one sample."),
        withSpinner(plotOutput(ns("density_plot"), height = 340), color = "#2563EB", type = 6),
        withSpinner(plotOutput(ns("boxplot_plot"), height = 340), color = "#2563EB", type = 6),
        fluidRow(
          column(6,
            actionButton(ns("bvsa_pca_btn"), "Generate before/after PCA", icon = icon("chart-line"), class = "btn-default btn-sm"),
            shinyjs::hidden(div(id = ns("bvsa_pca_wrap"), style = "margin-top:8px;", withSpinner(plotOutput(ns("bvsa_pca_plot"), height = 340), color = "#2563EB", type = 6)))
          ),
          column(6,
            actionButton(ns("bvsa_corr_btn"), "Generate before/after correlation", icon = icon("table-cells"), class = "btn-default btn-sm"),
            shinyjs::hidden(div(id = ns("bvsa_corr_wrap"), style = "margin-top:8px;",
              fluidRow(
                column(6, withSpinner(plotOutput(ns("bvsa_corr_before_plot"), height = 320), color = "#2563EB", type = 6)),
                column(6, withSpinner(plotOutput(ns("bvsa_corr_after_plot"), height = 320), color = "#2563EB", type = 6))
              )))
          )
        )
      )
    })
    observeEvent(input$bvsa_pca_btn, shinyjs::show(id = "bvsa_pca_wrap"), ignoreInit = TRUE)
    observeEvent(input$bvsa_corr_btn, shinyjs::show(id = "bvsa_corr_wrap"), ignoreInit = TRUE)

    beta_long_before_after <- function(r) {
      to_long <- function(m, label) {
        df <- as.data.frame(m); df$probe <- rownames(m)
        long <- tidyr::pivot_longer(df, -probe, names_to = "sample", values_to = "value")
        long$stage <- label; long
      }
      n_sub <- min(nrow(r$before), 20000)
      idx <- if (nrow(r$before) > n_sub) sample(nrow(r$before), n_sub) else seq_len(nrow(r$before))
      df <- rbind(to_long(r$before[idx, , drop = FALSE], "Before normalization"), to_long(r$after[idx, , drop = FALSE], "After normalization"))
      df$stage <- factor(df$stage, levels = c("Before normalization", "After normalization"))
      colnames(df)[colnames(df) == "value"] <- "beta"
      df
    }

    output$density_plot <- renderPlot({
      r <- norm_result()
      x_lab <- if (identical(methyl_dataset$input_scale, "beta") && !identical(r$method, "quantile")) "Beta value" else if (identical(methyl_dataset$input_scale, "m")) "M-value" else "Beta value"
      methyl_plot_density(beta_long_before_after(r), x_lab)
    })
    output$boxplot_plot <- renderPlot({
      req(norm_result())
      methyl_plot_boxplot(beta_long_before_after(norm_result()))
    })
    output$bvsa_pca_plot <- renderPlot({
      req(norm_result())
      r <- norm_result()
      pb <- methyl_pca_scores(r$before); pa <- methyl_pca_scores(r$after)
      validate(need(isTRUE(pb$ok) && isTRUE(pa$ok), "Not enough samples/probes for a before/after PCA."))
      mk <- function(pca, label) data.frame(x = pca$scores[, 1], y = pca$scores[, 2], color = label, text = rownames(pca$scores))
      df <- rbind(mk(pb, "Before normalization"), mk(pa, "After normalization"))
      methyl_plot_scatter2d(df, "PC1", "PC2", NULL, c("Before normalization" = ARTHOMIX_COLORS$blue, "After normalization" = ARTHOMIX_COLORS$orange)) + facet_wrap(~color)
    })
    output$bvsa_corr_before_plot <- renderPlot({
      req(norm_result())
      cb <- methyl_sample_correlation(norm_result()$before)
      validate(need(isTRUE(cb$ok), "Not enough samples/probes for a correlation heatmap."))
      methyl_plot_corr_heatmap(cb$cor) + labs(title = "Before normalization")
    })
    output$bvsa_corr_after_plot <- renderPlot({
      req(norm_result())
      ca <- methyl_sample_correlation(norm_result()$after)
      validate(need(isTRUE(ca$ok), "Not enough samples/probes for a correlation heatmap."))
      methyl_plot_corr_heatmap(ca$cor) + labs(title = "After normalization")
    })

    output$stats_body <- renderUI({
      req(norm_result())
      v <- norm_result()$validation
      rows <- list(
        c("Mean row variance", sprintf("%.5f", v$var_before), sprintf("%.5f", v$var_after)),
        c("Mean sample-sample correlation", sprintf("%.4f", v$mean_cor_before), sprintf("%.4f", v$mean_cor_after)),
        c("Type I/II distribution KS statistic", if (is.na(v$ks_before)) "n/a" else sprintf("%.4f", v$ks_before), if (is.na(v$ks_after)) "n/a" else sprintf("%.4f", v$ks_after)),
        c("Missingness (%)", sprintf("%.2f", v$missing_before), sprintf("%.2f", v$missing_after))
      )
      if (!is.null(v$signal_check)) {
        rows[[length(rows) + 1]] <- c("PC1 ~ group R² (signal-preservation check)", sprintf("%.4f", v$signal_check$r2_before), sprintf("%.4f", v$signal_check$r2_after))
      }
      df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
      colnames(df) <- c("Metric", "Before", "After")
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$filtered_body <- renderUI({
      req(norm_result())
      r <- norm_result()
      tagList(
        h5(sprintf("Removed probes (%d)", length(r$removed_probes))),
        if (length(r$removed_probes) > 0) DT::dataTableOutput(ns("removed_probes_table")) else p(class = "empty-note", icon("circle-check"), "No probes were removed."),
        h5(sprintf("Removed samples (%d)", length(r$removed_samples))),
        if (length(r$removed_samples) > 0) DT::dataTableOutput(ns("removed_samples_table")) else p(class = "empty-note", icon("circle-check"), "No samples were removed.")
      )
    })
    output$removed_probes_table <- DT::renderDataTable({
      req(norm_result())
      DT::datatable(data.frame(probe_id = norm_result()$removed_probes), rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })
    output$removed_samples_table <- DT::renderDataTable({
      req(norm_result())
      DT::datatable(data.frame(sample_id = norm_result()$removed_samples), rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    output$details_body <- renderUI({
      req(norm_result())
      r <- norm_result()
      tags$pre(style = "white-space:pre-wrap; font-size:12px;",
        methyl_norm_processing_record(list(
          method_label = r$method_label, dataset_source = methyl_dataset$source, run_at = r$run_at,
          n_samples_before = r$n_samples_before, n_samples_after = ncol(r$after),
          n_probes_before = r$n_probes_before, n_probes_after = nrow(r$after),
          filter_notes = r$filter_notes, note = r$note
        ))
      )
    })

    ## ---- Compare Methods ---------------------------------------------

    compare_result <- eventReactive(input$compare_btn, {
      sel <- input$compare_methods_sel
      validate(need(length(sel) >= 2, "Select at least two methods to compare."))
      out <- withProgress(message = "Comparing methods", value = 0, {
        lapply(seq_along(sel), function(i) {
          incProgress(1 / length(sel), detail = names(available_methods())[available_methods() == sel[i]])
          r <- run_full(sel[i])
          r$key <- sel[i]
          r
        })
      })
      names(out) <- sel
      out
    })

    output$compare_result_ui <- renderUI({
      req(compare_result())
      res <- compare_result()
      ok <- Filter(function(r) isTRUE(r$ok), res)
      failed <- Filter(function(r) !isTRUE(r$ok), res)
      tagList(
        if (length(failed) > 0) div(class = "empty-note", icon("triangle-exclamation"),
          sprintf("%d method(s) could not be compared: %s.", length(failed),
            paste(sprintf("%s (%s)", names(available_methods())[match(names(failed), available_methods())], vapply(failed, `[[`, character(1), "reason")), collapse = "; "))),
        if (length(ok) >= 1) tagList(
          DT::dataTableOutput(ns("compare_table")),
          withSpinner(plotOutput(ns("compare_density_plot"), height = 380), color = "#2563EB", type = 6)
        ) else p(class = "empty-note", icon("circle-info"), "No selected method could be run on this dataset with the current filters.")
      )
    })

    output$compare_table <- DT::renderDataTable({
      res <- Filter(function(r) isTRUE(r$ok), compare_result())
      df <- do.call(rbind, lapply(res, function(r) data.frame(
        method = r$method_label, samples = ncol(r$after), probes = nrow(r$after),
        mean_variance = r$validation$var_after, mean_correlation = r$validation$mean_cor_after,
        type_ks_stat = r$validation$ks_after, stringsAsFactors = FALSE
      )))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatRound(columns = c("mean_variance", "mean_correlation", "type_ks_stat"), digits = 4)
    })

    output$compare_density_plot <- renderPlot({
      res <- Filter(function(r) isTRUE(r$ok), compare_result())
      validate(need(length(res) > 0, "No result to plot."))
      base_before <- res[[1]]$before
      n_sub <- min(nrow(base_before), 15000)
      idx <- if (nrow(base_before) > n_sub) sample(nrow(base_before), n_sub) else seq_len(nrow(base_before))
      to_long <- function(m, label) {
        df <- as.data.frame(m[idx, , drop = FALSE]); df$probe <- rownames(df)
        long <- tidyr::pivot_longer(df, -probe, names_to = "sample", values_to = "beta")
        long$method <- label; long
      }
      parts <- c(list(to_long(base_before, "Before normalization")), lapply(res, function(r) to_long(r$after, r$method_label)))
      df <- do.call(rbind, parts)
      pal <- arthomix_pair(df$method)
      ggplot(df, aes(x = beta, color = method)) +
        geom_density(linewidth = 0.5) +
        scale_color_manual(values = pal) +
        labs(x = "Value", y = "Density", color = NULL) + theme_arthomix()
    })
  })
}
