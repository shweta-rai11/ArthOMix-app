## R/multiomics/03_DIABLO_SNF_Integration/mod_multi_integration.R
## Submodule: Multi-omics Integration - a live, data-adaptive DIABLO
## (mixOmics::block.splsda, supervised) and SNF (SNFtool::SNF, unsupervised)
## engine, plus a head-to-head Compare tab against single-omics performance.
## Runs on whichever data is actually selected below - the Active
## Multi-Omics Dataset built on the Dataset Workspace tab, or one preloaded
## RA anti-TNF analysis cell (recomputed live from that cell's own saved
## fit - multiomics_integration_helpers.R::mi_preloaded_cell_dataset())
## - never a fixed template. Every parameter range, tuning grid, and CV
## fold/repeat count is derived from the data in front of it
## (multiomics_integration_helpers.R); nothing below renders a result,
## score, plot, or table until its own blue "Run" button is clicked.
##
## DIABLO and SNF keep entirely separate parameter panels and result
## sections (never mixed) since they answer different questions - "which
## features discriminate a known outcome" vs. "do patients cluster by
## molecular similarity, independent of any outcome". "Outcome variable" is
## the one exception kept in the shared Data Selection panel rather than
## duplicated per tab: it's a property of the dataset (which metadata
## column is a usable label), not a DIABLO- or SNF-specific model setting -
## SNF's own panel never reads it.

mod_multi_integration_config <- list(
  id = "integration", title = "Multi-omics Integration (DIABLO & SNF)", icon = "layer-group", group = "Data",
  description = "DIABLO (supervised) and SNF (unsupervised) integration, with single-omics comparison and sex-stratified analysis."
)

## ---------------------------------------------------------------------------
## Small shared UI pieces
## ---------------------------------------------------------------------------

mi_stat_card <- function(value, label, color = ARTHOMIX_COLORS$blue) {
  div(class = "card", style = "flex:1 1 140px; text-align:center; padding:10px;",
      div(style = sprintf("font-size:1.3em; font-weight:600; color:%s;", color), value),
      div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label))
}

mi_warn <- function(...) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"), " ", ...)
mi_ok   <- function(...) div(class = "empty-note", icon("circle-check"), " ", ...)
mi_stop <- function(...) div(class = "empty-note", style = "border-color: var(--color-danger, #e34948);", icon("circle-xmark"), " ", ...)

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_multi_integration_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("active_dataset_banner")),

    tabsetPanel(
      id = ns("tabs"), type = "tabs",
      tabPanel(
        "Data", br(),
        box(
          width = NULL, title = "1. Data selection", status = "primary", solidHeader = FALSE,
          radioButtons(ns("data_source"), "Data source",
                       choices = c("Active Multi-Omics Dataset (Dataset Workspace)" = "active", "Reference / Example Dataset (RA anti-TNF cohort)" = "preloaded"),
                       selected = "active", inline = TRUE),
          conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
                            selectInput(ns("preloaded_cell"), "Analysis cell", choices = MULTI_CELL_CHOICES, width = "100%")),
          uiOutput(ns("outcome_ui")),
          uiOutput(ns("source_note"))
        ),
        uiOutput(ns("validation_ui"))
      ),
      tabPanel("DIABLO", br(), fluidRow(column(4, uiOutput(ns("diablo_params_ui"))), column(8, uiOutput(ns("diablo_results_ui"))))),
      tabPanel("SNF", br(), fluidRow(column(4, uiOutput(ns("snf_params_ui"))), column(8, uiOutput(ns("snf_results_ui"))))),
      tabPanel("Compare", br(), uiOutput(ns("compare_ui"))),
      tabPanel("Sex-Stratified", br(), fluidRow(column(4, uiOutput(ns("ss_params_ui"))), column(8, uiOutput(ns("ss_results_ui")))))
    )
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_multi_integration_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$active_dataset_banner <- renderUI(multi_active_dataset_banner(multi_dataset))

    ## =========================================================================
    ## 1. Data selection - one reactive dataset object feeds everything below,
    ## whether it came from the preloaded-cell adapter or the shared
    ## multi_dataset (Dataset Workspace's Active Multi-Omics Dataset, already
    ## N-omics-generic). Never merges the two.
    ## =========================================================================
    mi_dataset <- reactive({
      if (identical(input$data_source, "preloaded")) {
        req(input$preloaded_cell)
        mi_preloaded_cell_dataset(input$preloaded_cell)
      } else {
        if (is.null(multi_dataset) || !isTRUE(multi_dataset$active) || length(multi_dataset$layers %||% list()) < 1) {
          return(list(ok = FALSE, error = "No Active Multi-Omics Dataset yet. Build one on the Dataset Workspace tab, or switch to \"Preloaded RA anti-TNF cohort\" above."))
        }
        list(
          ok = TRUE, layers = multi_dataset$layers, sample_meta = multi_dataset$sample_meta, outcome_col = NULL,
          label = sprintf("Active Multi-Omics Dataset (%s)", paste(names(multi_dataset$layers), collapse = " + ")),
          provenance = sprintf("Active Multi-Omics Dataset from the Dataset Workspace tab (source: %s).", multi_dataset$source %||% "unknown")
        )
      }
    })

    output$source_note <- renderUI({
      d <- mi_dataset()
      if (!isTRUE(d$ok)) return(mi_warn(d$error))
      mi_ok(d$provenance)
    })

    outcome_candidates <- reactive({
      d <- mi_dataset()
      if (!isTRUE(d$ok) || is.null(d$sample_meta)) return(character(0))
      if (!is.null(d$outcome_col)) return(d$outcome_col)
      colnames(d$sample_meta)
    })

    output$outcome_ui <- renderUI({
      cands <- outcome_candidates()
      if (length(cands) == 0) return(mi_warn("No sample metadata with a candidate outcome column is available for this dataset."))
      selectInput(ns("outcome_col"), "Outcome variable", choices = cands, selected = cands[1], width = "100%")
    })

    ## =========================================================================
    ## 2. Common validation (spec section 5) - always visible once a dataset
    ## is selected, regardless of which subtab is open. Sample matching is
    ## always by matched sample ID (mi_validate_dataset() -> the shared
    ## multi_live_sample_overlap()), never by row position.
    ## =========================================================================
    ## mi_val() deliberately does NOT depend on input$outcome_col - block
    ## selection/tuning-mode widgets below are rendered from mi_val() alone,
    ## so switching the outcome variable never resets a user's in-progress
    ## parameter choices. The outcome summary is a separate, smaller reactive
    ## that DOES depend on it.
    mi_val <- reactive({
      d <- mi_dataset()
      if (!isTRUE(d$ok)) return(NULL)
      mi_validate_dataset(d$layers, d$sample_meta, outcome_col = NULL)
    })
    ## Block-scoped counterpart of mi_val() - eligibility/overlap/"planned
    ## run" figures for DIABLO and SNF must be computed from only the blocks
    ## the user actually checked (input$d_blocks/input$s_blocks), not every
    ## block in the active dataset. Using the whole-dataset mi_val() for
    ## these meant deselecting a poorly-matched third block never recovered
    ## the larger sample overlap that the two remaining blocks alone would
    ## allow, and an unselected block's own missing values could block SNF
    ## even though the selected blocks were complete. mi_val() itself is
    ## still used for the dataset-wide "2. Data validation" cards, which are
    ## deliberately meant to describe the whole active dataset.
    val_for_blocks <- function(blocks) {
      d <- mi_dataset()
      if (!isTRUE(d$ok)) return(NULL)
      ## NULL covers both "checkboxGroupInput hasn't rendered yet" and "user
      ## unchecked every block" (Shiny reports both as NULL) - defaulting to
      ## every block in that case reproduces the old whole-dataset behavior
      ## for the ambiguous case while still narrowing correctly the moment a
      ## real (non-empty) selection exists.
      blocks <- if (is.null(blocks)) names(d$layers) else intersect(blocks, names(d$layers))
      if (length(blocks) < 1) return(NULL)
      mi_validate_dataset(d$layers[blocks], d$sample_meta, outcome_col = NULL)
    }
    mi_outcome <- reactive({
      d <- mi_dataset(); v <- mi_val()
      if (!isTRUE(d$ok) || is.null(v) || is.null(input$outcome_col)) return(NULL)
      mi_outcome_summary(d$sample_meta, input$outcome_col, v$shared_ids)
    })

    output$validation_ui <- renderUI({
      d <- mi_dataset(); v <- mi_val(); o <- mi_outcome()
      if (!isTRUE(d$ok) || is.null(v)) return(NULL)
      if (!isTRUE(v$ok)) return(box(width = NULL, title = "2. Data validation", status = "primary", solidHeader = FALSE, mi_warn(v$error)))

      block_cards <- lapply(names(v$per_block), function(nm) {
        b <- v$per_block[[nm]]
        if (!isTRUE(b$ok)) return(mi_stat_card("-", nm, ARTHOMIX_COLORS$red))
        div(class = "card", style = "flex:1 1 200px; padding:12px;",
            div(style = "font-weight:600; margin-bottom:4px;", nm),
            div(style = "font-size:0.85em;", sprintf("Samples: %s | Features: %s", format(b$n_samples, big.mark = ","), format(b$n_features, big.mark = ","))),
            div(style = "font-size:0.85em;", sprintf("Missing: %.1f%% | Zero-variance: %d | Dup. samples: %d | Dup. features: %d", b$pct_missing, b$n_zero_variance, b$n_duplicate_samples, b$n_duplicate_features)))
      })

      outcome_note <- if (!is.null(o)) {
        if (identical(o$type, "categorical")) {
          counts_txt <- paste(sprintf("%s = %d", names(o$class_counts), as.integer(o$class_counts)), collapse = ", ")
          tagList(
            div(sprintf("Outcome \"%s\": %d classes (%s), n = %d matched.", o$column, o$n_classes, counts_txt, o$n)),
            if (isTRUE(o$imbalanced)) mi_warn("Class sizes are substantially imbalanced - prefer BER over overall error rate.")
          )
        } else div(sprintf("Outcome \"%s\" looks continuous (%d values) - DIABLO needs a categorical outcome.", o$column, o$n))
      } else div(class = "submodule-desc", "Select an outcome variable above to see its class breakdown.")

      box(
        width = NULL, title = "2. Data validation", status = "primary", solidHeader = FALSE,
        div(style = "display:flex; gap:12px; flex-wrap:wrap; margin-bottom:10px;", block_cards),
        div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px;",
            mi_stat_card(v$n_blocks, "Omics blocks"),
            mi_stat_card(v$n_shared, "Matched samples across all blocks", ARTHOMIX_COLORS$aqua)),
        outcome_note,
        if (!isTRUE(v$reliable_matching)) mi_stop(v$mismatch_message) else mi_ok(sprintf("%d matched samples will be used for integration.", v$n_shared))
      )
    })

    ## =========================================================================
    ## 3. DIABLO subtab
    ## =========================================================================
    ## `submitted` is set eagerly, synchronously, in the same observer that
    ## calls diablo_task$invoke() - `future::multisession` worker cold-start
    ## (spinning up a background R process the first time, loading mixOmics
    ## in it) can itself take a few real seconds, during which
    ## diablo_task$status() can still read "initial" rather than "running".
    ## Without this flag, every render branch below fell through to nothing
    ## visible during that gap (confirmed directly: an empty/unchanged panel
    ## for several seconds after clicking "Run DIABLO", not just a briefly
    ## slow spinner) - `submitted` closes that gap so the spinner appears
    ## the instant the button is clicked, not once the async machinery
    ## catches up.
    diablo_state <- reactiveValues(result = NULL, error = NULL, submitted = FALSE)

    diablo_elig <- reactive({
      v <- val_for_blocks(input$d_blocks)
      if (is.null(v)) return(list(ok = FALSE, reason = "Select a dataset and at least one omics block above first."))
      o <- if (is.null(input$outcome_col)) NULL else mi_outcome_summary(mi_dataset()$sample_meta, input$outcome_col, v$shared_ids)
      mi_diablo_eligibility(v, o)
    })

    output$diablo_params_ui <- renderUI({
      d <- mi_dataset(); v <- mi_val()
      if (!isTRUE(d$ok) || is.null(v) || !isTRUE(v$ok)) return(box(width = NULL, title = "DIABLO parameters", status = "primary", solidHeader = FALSE, mi_warn("Select a valid dataset above first.")))
      elig <- diablo_elig()
      blocks <- v$block_labels
      v_sel <- val_for_blocks(input$d_blocks) %||% v
      n <- v_sel$n_shared
      loo_ok <- mi_diablo_loo_feasible(n)
      o_sel <- if (is.null(input$outcome_col)) NULL else mi_outcome_summary(mi_dataset()$sample_meta, input$outcome_col, v_sel$shared_ids)
      min_class_n <- if (!is.null(o_sel$class_counts)) min(o_sel$class_counts) else max(2, floor(n / 2))

      tagList(
        box(
          width = NULL, title = "DIABLO parameters", status = "primary", solidHeader = FALSE,
          if (!isTRUE(elig$ok)) mi_stop(elig$reason) else mi_ok("DIABLO is available for this dataset."),
          div(style = if (!isTRUE(elig$ok)) "opacity:0.5; pointer-events:none;" else NULL,
            h5("Data"),
            checkboxGroupInput(ns("d_blocks"), "Omics blocks", choices = blocks, selected = blocks),
            p(class = "submodule-desc", "Outcome variable is set once, above (Data selection) - it applies to both DIABLO and Compare."),
            hr(),
            h5("Components"),
            numericInput(ns("d_ncomp"), "Number of components (ncomp)", value = max(mi_diablo_feasible_ncomp(o_sel$n_classes %||% 2, min_class_n)), min = 1, max = max(mi_diablo_feasible_ncomp(o_sel$n_classes %||% 2, min_class_n)), step = 1),
            p(class = "submodule-desc", sprintf("Feasible range for this dataset: 1-%d.", max(mi_diablo_feasible_ncomp(o_sel$n_classes %||% 2, min_class_n)))),
            hr(),
            h5("Feature selection"),
            uiOutput(ns("d_keepx_manual_ui")),
            p(class = "submodule-desc", "Features retained per omics block, per component (comma-separated)."),
            checkboxInput(ns("d_keepx_auto"), "Auto-tune keepX instead (mixOmics::tune.block.splsda() grid search - slower)", value = FALSE),
            hr(),
            h5("Block relationship (design)"),
            uiOutput(ns("d_design_ui")),
            p(class = "submodule-desc", "Cross-block weight in DIABLO's design matrix - mixOmics default is 0.1."),
            hr(),
            h5("Validation"),
            radioButtons(ns("d_validation_method"), "Method",
                         choices = c("M-fold CV" = "mfold", if (loo_ok) c("Leave-one-out" = "loo")), selected = "mfold", inline = TRUE),
            if (!loo_ok) p(class = "submodule-desc", sprintf("Leave-one-out is disabled (n = %d exceeds %d).", n, MI_DIABLO_LOO_MAX_N)),
            conditionalPanel(condition = sprintf("input['%s'] == 'mfold'", ns("d_validation_method")),
              numericInput(ns("d_folds"), "Number of folds", value = mi_diablo_feasible_folds(min_class_n), min = 2, max = mi_diablo_max_folds(min_class_n), step = 1),
              numericInput(ns("d_nrepeat"), "Repeats", value = mi_diablo_feasible_repeats(n), min = 1, max = 50, step = 1)
            ),
            hr(),
            h5("Prediction distance"),
            selectInput(ns("d_distance"), NULL, choices = c("Automatically selected during tuning" = "automatic", "centroids.dist" = "centroids.dist", "mahalanobis.dist" = "mahalanobis.dist", "max.dist" = "max.dist")),
            tags$details(tags$summary("Advanced parameters"),
              numericInput(ns("d_tol"), "Convergence tolerance", value = 1e-06, min = 1e-9, max = 1e-2),
              numericInput(ns("d_maxiter"), "Maximum iterations", value = 100, min = 10, max = 1000, step = 10),
              numericInput(ns("d_seed"), "Random seed", value = 1, min = 1)
            ),
            hr(),
            uiOutput(ns("d_summary_pre")),
            if (isTRUE(elig$ok)) div(
              if (isTRUE(input$d_keepx_auto)) mi_warn("Auto-tuning keepX (checked above) may take several minutes for high-dimensional datasets.") else NULL,
              actionButton(ns("d_run_btn"), "Run DIABLO", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")),
            uiOutput(ns("d_run_status_ui"))
          )
        )
      )
    })

    output$d_keepx_manual_ui <- renderUI({
      v <- mi_val(); req(v)
      tagList(lapply(v$block_labels, function(b) textInput(ns(paste0("d_keepx_", make.names(b))), sprintf("%s keepX per component", b), value = "20,10")))
    })
    output$d_design_ui <- renderUI({
      blocks <- req(input$d_blocks)
      if (length(blocks) < 2) return(mi_warn("Select at least two blocks."))
      pairs <- utils::combn(blocks, 2, simplify = FALSE)
      tagList(lapply(pairs, function(p) sliderInput(ns(paste0("d_design_", make.names(p[1]), "_", make.names(p[2]))), sprintf("%s <-> %s", p[1], p[2]), min = 0, max = 1, value = 0.1, step = 0.05)))
    })

    d_params <- reactive({
      ## The design (block-relationship) and validation sliders/fields are
      ## now always visible, always real numbers - so this always builds a
      ## concrete design_custom/validation config rather than choosing
      ## between an "Automatic" constant and a hidden "Custom" panel. Only
      ## keepX keeps a real Automatic-vs-Manual distinction, since Automatic
      ## keepX is a genuine (slow) grid search, not a fixed default like the
      ## others - driven by the single "Auto-tune keepX instead" checkbox.
      blocks <- input$d_blocks %||% character(0)
      design_custom <- NULL
      if (length(blocks) >= 2) {
        design_custom <- matrix(0.1, length(blocks), length(blocks), dimnames = list(blocks, blocks))
        pairs <- utils::combn(blocks, 2, simplify = FALSE)
        for (p in pairs) {
          val <- input[[paste0("d_design_", make.names(p[1]), "_", make.names(p[2]))]] %||% 0.1
          design_custom[p[1], p[2]] <- val; design_custom[p[2], p[1]] <- val
        }
      }
      keepx_manual <- NULL
      if (!isTRUE(input$d_keepx_auto)) {
        keepx_manual <- stats::setNames(lapply(blocks, function(b) {
          txt <- input[[paste0("d_keepx_", make.names(b))]] %||% "10"
          as.integer(trimws(strsplit(txt, ",")[[1]]))
        }), blocks)
      }
      list(
        design_mode = "custom", design_custom = design_custom,
        ncomp_mode = "manual", ncomp = input$d_ncomp,
        keepx_mode = if (isTRUE(input$d_keepx_auto)) "automatic" else "manual", keepx_manual = keepx_manual,
        validation_mode = "manual", validation_method = input$d_validation_method %||% "mfold",
        folds = input$d_folds, nrepeat = input$d_nrepeat, distance = input$d_distance %||% "automatic",
        ## Threaded straight into mixOmics::tune.block.splsda()/perf()'s own
        ## `seed=` argument inside mi_diablo_run() - an external set.seed()
        ## call around this reactive would NOT work, since both of those
        ## mixOmics calls unconditionally reseed internally from their own
        ## `seed` argument (default NULL, which re-randomizes every call).
        seed = input$d_seed %||% 1
      )
    })

    output$d_summary_pre <- renderUI({
      v <- req(val_for_blocks(input$d_blocks)); p <- d_params()
      tags$div(class = "submodule-desc",
        tags$strong("Planned run: "), sprintf(
          "Blocks: %s | Samples: %d | Outcome: %s | Components: %s | keepX: %s | Validation: %s, %s-fold x %s repeats | Distance: %s",
          paste(input$d_blocks, collapse = " + "), v$n_shared, input$outcome_col %||% "-",
          p$ncomp %||% "-",
          if (identical(p$keepx_mode, "automatic")) "Auto-tuned (grid search)" else "Manual (set above)",
          if (identical(p$validation_method, "loo")) "Leave-one-out" else "M-fold CV", p$folds %||% "-", p$nrepeat %||% "-",
          if (identical(p$distance, "automatic")) "Auto-selected during tuning" else p$distance
        ))
    })

    run_diablo <- function() {
      d <- mi_dataset(); v <- mi_val()
      layers <- d$layers[input$d_blocks]
      outcome <- stats::setNames(d$sample_meta[[input$outcome_col]], rownames(d$sample_meta))
      mi_diablo_run(layers, outcome, v$shared_ids, params = d_params())
    }

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      ## The "Random seed" field (d_seed) rides along inside `params` itself
      ## (d_params() above) straight into mi_diablo_run()'s own
      ## mixOmics::tune.block.splsda()/perf() `seed=` arguments - note
      ## future_promise()'s own `seed = TRUE` here is unrelated: it only
      ## asks `future` for a parallel-safe RNG stream for whatever runs
      ## inside the worker process, not the user-chosen integer in the UI.
      diablo_task <- ExtendedTask$new(function(layers, outcome, ids, params) {
        promises::future_promise(mi_diablo_run(layers, outcome, ids, params), seed = TRUE)
      })
      observeEvent(input$d_run_btn, {
        v <- req(val_for_blocks(input$d_blocks)); d <- req(mi_dataset())
        validate(need(isTRUE(diablo_elig()$ok), diablo_elig()$reason))
        layers <- d$layers[input$d_blocks]
        outcome <- stats::setNames(d$sample_meta[[input$outcome_col]], rownames(d$sample_meta))
        diablo_state$error <- NULL
        diablo_state$result <- NULL
        diablo_state$submitted <- TRUE
        ## Pinned at invoke time - diablo_results_ui reads this, not
        ## input$outcome_col live, so changing the outcome dropdown after a
        ## run never relabels the sample-score plot against the wrong column.
        diablo_state$outcome_used <- outcome
        ## Pinned alongside outcome_used so the Compare tab can tell whether
        ## input$outcome_col has since been changed by the user - without
        ## this, c_supervised() below would silently compare DIABLO's fit
        ## against a different outcome than it was actually trained on.
        diablo_state$outcome_col_used <- input$outcome_col
        p <- d_params()
        ids <- v$shared_ids
        ## diablo_task$invoke() calls its ExtendedTask function body
        ## SYNCHRONOUSLY before returning (confirmed directly by timing it:
        ## tens of seconds in this app's actual reactive-domain context,
        ## not the near-instant dispatch its docs describe in isolation) -
        ## so calling it inline here would block THIS SAME observer, and
        ## nothing (including the `submitted` flag just set above) can
        ## flush to the browser until a running observer/render finishes.
        ## session$onFlushed(..., once = TRUE) defers the actual invoke()
        ## to the next flush cycle, so the "Running..." spinner state set
        ## above reaches the browser first. Every argument is snapshotted
        ## here (not read live inside the deferred callback), since that
        ## callback runs outside a normal reactive-consumer context.
        session$onFlushed(function() diablo_task$invoke(layers, outcome, ids, p), once = TRUE)
      })
      observe({
        res <- tryCatch(diablo_task$result(), error = function(e) e)
        if (inherits(res, "shiny.silent.error")) return()
        diablo_state$submitted <- FALSE
        if (inherits(res, "error")) { diablo_state$error <- conditionMessage(res); return() }
        if (!isTRUE(res$ok)) { diablo_state$error <- res$error; return() } else { diablo_state$result <- res }
      })
      diablo_running <- reactive(identical(diablo_task$status(), "running") || isTRUE(diablo_state$submitted))
    } else {
      observeEvent(input$d_run_btn, {
        validate(need(isTRUE(diablo_elig()$ok), diablo_elig()$reason))
        showNotification("Running DIABLO synchronously - the app will be briefly unresponsive.", type = "message", duration = 5)
        diablo_state$outcome_used <- stats::setNames(mi_dataset()$sample_meta[[input$outcome_col]], rownames(mi_dataset()$sample_meta))
        diablo_state$outcome_col_used <- input$outcome_col
        res <- run_diablo()
        diablo_state$error <- if (!isTRUE(res$ok)) res$error else NULL
        if (isTRUE(res$ok)) diablo_state$result <- res
      })
      diablo_running <- reactive(FALSE)
    }

    ## Mutates the existing "Run DIABLO" button in place via shinyjs (already
    ## active app-wide, ui.R::useShinyjs()) rather than recreating it through
    ## renderUI - recreating an actionButton resets its click counter to 0,
    ## which would break the "has this been run at least once" gating used
    ## throughout this tab (input$d_run_btn > 0).
    observe({
      if (isTRUE(diablo_running())) {
        shinyjs::disable(ns("d_run_btn"))
        shinyjs::html(ns("d_run_btn"), as.character(tagList(icon("spinner", class = "fa-spin"), " Running DIABLO...")))
      } else {
        shinyjs::enable(ns("d_run_btn"))
        shinyjs::html(ns("d_run_btn"), as.character(tagList(icon("play"), " Run DIABLO")))
      }
    })

    ## Shown right next to the "Run DIABLO" button itself (not just in the
    ## results column, which can be out of view or easy to miss while
    ## automatic tuning genuinely takes real time) - unmistakable feedback at
    ## the exact point the user just clicked.
    output$d_run_status_ui <- renderUI({
      if (!isTRUE(input$d_run_btn > 0)) return(NULL)
      if (isTRUE(diablo_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running DIABLO - see the results panel once this finishes."))
      if (!is.null(diablo_state$error)) return(mi_stop(diablo_state$error))
      if (!is.null(diablo_state$result)) return(mi_ok("Finished - see the results panel to the right."))
      NULL
    })

    output$diablo_results_ui <- renderUI({
      if (!isTRUE(input$d_run_btn > 0)) return(multi_empty_state("Set parameters and click \"Run DIABLO\" to see results here."))
      if (isTRUE(diablo_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running DIABLO - the rest of the app stays usable while this runs."))
      if (!is.null(diablo_state$error)) return(mi_stop(diablo_state$error))
      res <- diablo_state$result
      if (is.null(res)) return(multi_empty_state())
      perf_sum <- mi_diablo_performance_summary(res)
      sel <- mi_diablo_selected_features_df(res$fit)
      scores_df <- mi_diablo_sample_scores_df(res$fit, diablo_state$outcome_used)
      p <- res$params
      keepx_txt <- paste(sprintf("%s: %s", names(p$keepX), vapply(p$keepX, paste, character(1), collapse = ",")), collapse = " | ")

      tagList(
        box(width = NULL, title = "Actual parameters used", status = "primary", solidHeader = FALSE,
            tags$ul(
              tags$li(sprintf("Components: %d", p$ncomp)),
              tags$li(sprintf("keepX per block: %s", keepx_txt)),
              tags$li(sprintf("Block relationship: %s", p$design_mode)),
              tags$li(sprintf("Prediction distance: %s%s", p$distance, if (identical(p$distance_mode, "automatic")) " (auto-selected)" else "")),
              tags$li(sprintf("Validation: %s, %d-fold x %d repeat(s)", p$validation_method, p$folds, p$nrepeat)),
              if (isTRUE(p$loo_downgraded)) tags$li(mi_warn("Leave-one-out was requested but this dataset is too large for it - M-fold CV was used instead."))
            )),
        box(width = NULL, title = "Model", status = "primary", solidHeader = FALSE,
            div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                mi_stat_card("Valid", "Model status", ARTHOMIX_COLORS$aqua),
                mi_stat_card(p$n_samples, "Samples"), mi_stat_card(length(p$blocks), "Blocks"), mi_stat_card(p$ncomp, "Components"))),
        box(width = NULL, title = "Performance (cross-validated)", status = "primary", solidHeader = FALSE,
            if (is.null(perf_sum)) mi_warn("Performance could not be computed for this configuration.") else tagList(
              div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                  mi_stat_card(sprintf("%.3f", perf_sum$ber), "Balanced error rate (BER)"),
                  mi_stat_card(sprintf("%.3f", perf_sum$overall_error), "Overall error rate")),
              multi_plot_or_empty(function() mi_diablo_error_bar_plot(perf_sum), ns("d_error_plot"), height = "260px"),
              if (!is.null(perf_sum$auc)) tagList(h5("AUC"), DT::dataTableOutput(ns("d_auc_table")))
            )),
        box(width = NULL, title = "Selected multi-omics features", status = "primary", solidHeader = FALSE,
            if (is.null(sel)) mi_warn("No features were selected.") else tagList(
              multi_plot_or_empty(function() multi_diablo_panel_plot(mi_diablo_panel_df_for_plot(sel, 1)), ns("d_panel_plot"), height = "380px"),
              DT::dataTableOutput(ns("d_selected_table")),
              div(class = "table-toolbar", downloadButton(ns("d_dl_selected"), "Download selected features (CSV)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Sample plot", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() multi_diablo_score_plot(scores_df), ns("d_score_plot"), height = "320px")),
        box(width = NULL, title = "Variance explained (within block)", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() multi_diablo_variance_plot(multi_diablo_variance_df(res$fit)), ns("d_variance_plot"), height = "300px")),
        if (length(p$blocks) >= 2) box(width = NULL, title = "Correlation between selected variables across blocks", status = "primary", solidHeader = FALSE,
            fluidRow(column(6, selectInput(ns("d_corr_a"), "Block A", choices = p$blocks, selected = p$blocks[1])),
                     column(6, selectInput(ns("d_corr_b"), "Block B", choices = p$blocks, selected = p$blocks[min(2, length(p$blocks))]))),
            uiOutput(ns("d_corr_plot_ui"))),
        box(width = NULL, title = "Feature-selection stability", status = "primary", solidHeader = FALSE, uiOutput(ns("d_stability_ui")))
      )
    })

    output$d_auc_table <- DT::renderDataTable({
      req(diablo_state$result)
      DT::datatable(mi_diablo_performance_summary(diablo_state$result)$auc, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$d_selected_table <- DT::renderDataTable({
      req(diablo_state$result)
      DT::datatable(mi_diablo_selected_features_df(diablo_state$result$fit), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$d_dl_selected <- downloadHandler(function() "diablo_selected_features.csv", function(file) utils::write.csv(mi_diablo_selected_features_df(diablo_state$result$fit), file, row.names = FALSE))

    ## multi_plot_or_empty() only decides whether to draw the plotOutput
    ## placeholder or an empty-state message - it does not register the
    ## render logic. Each placeholder above needs its own renderPlot binding.
    output$d_error_plot <- renderPlot({
      res <- req(diablo_state$result)
      perf_sum <- mi_diablo_performance_summary(res)
      req(perf_sum)
      mi_diablo_error_bar_plot(perf_sum)
    })
    output$d_panel_plot <- renderPlot({
      res <- req(diablo_state$result)
      sel <- req(mi_diablo_selected_features_df(res$fit))
      multi_diablo_panel_plot(mi_diablo_panel_df_for_plot(sel, 1))
    })
    output$d_score_plot <- renderPlot({
      res <- req(diablo_state$result)
      multi_diablo_score_plot(mi_diablo_sample_scores_df(res$fit, diablo_state$outcome_used))
    })
    output$d_variance_plot <- renderPlot({
      res <- req(diablo_state$result)
      multi_diablo_variance_plot(multi_diablo_variance_df(res$fit))
    })

    output$d_corr_plot_ui <- renderUI({
      req(diablo_state$result, input$d_corr_a, input$d_corr_b, input$d_corr_a != input$d_corr_b)
      v <- mi_val(); d <- mi_dataset()
      sel <- mi_diablo_selected_features_df(diablo_state$result$fit)
      corr <- mi_diablo_selected_correlation_data(d$layers, sel, input$d_corr_a, input$d_corr_b, v$shared_ids)
      if (!isTRUE(corr$ok)) return(mi_warn(corr$error))
      multi_plot_or_empty(function() multi_live_correlation_heatmap_plot(corr$df), ns("d_corr_plot"), height = "420px")
    })
    output$d_corr_plot <- renderPlot({
      req(diablo_state$result, input$d_corr_a, input$d_corr_b, input$d_corr_a != input$d_corr_b)
      v <- mi_val(); d <- mi_dataset()
      sel <- mi_diablo_selected_features_df(diablo_state$result$fit)
      corr <- mi_diablo_selected_correlation_data(d$layers, sel, input$d_corr_a, input$d_corr_b, v$shared_ids)
      req(isTRUE(corr$ok))
      multi_live_correlation_heatmap_plot(corr$df)
    })
    output$d_stability_ui <- renderUI({
      res <- req(diablo_state$result)
      stab <- mi_diablo_stability_df(res, block = res$params$blocks[1], comp = 1)
      if (is.null(stab)) return(mi_warn("Feature-selection stability is not available for this configuration (needs more than one CV repeat)."))
      DT::dataTableOutput(ns("d_stability_table"))
    })
    output$d_stability_table <- DT::renderDataTable({
      res <- req(diablo_state$result)
      DT::datatable(mi_diablo_stability_df(res, block = res$params$blocks[1], comp = 1), rownames = FALSE, options = list(pageLength = 10), class = "stripe hover compact")
    })

    ## =========================================================================
    ## 4. SNF subtab
    ## =========================================================================
    ## `submitted` - same rationale as diablo_state$submitted above.
    snf_state <- reactiveValues(result = NULL, error = NULL, submitted = FALSE)

    snf_elig <- reactive({
      v <- val_for_blocks(input$s_blocks)
      if (is.null(v)) return(list(ok = FALSE, reason = "Select a dataset and at least one omics block above first."))
      mi_snf_eligibility(v)
    })

    output$snf_params_ui <- renderUI({
      d <- mi_dataset(); v <- mi_val()
      if (!isTRUE(d$ok) || is.null(v) || !isTRUE(v$ok)) return(box(width = NULL, title = "SNF parameters", status = "primary", solidHeader = FALSE, mi_warn("Select a valid dataset above first.")))
      elig <- snf_elig()
      blocks <- v$block_labels
      v_sel <- val_for_blocks(input$s_blocks) %||% v
      n <- v_sel$n_shared
      k_range <- mi_snf_feasible_k_range(n)

      box(
        width = NULL, title = "SNF parameters", status = "primary", solidHeader = FALSE,
        if (!isTRUE(elig$ok)) mi_stop(elig$reason) else mi_ok("SNF is available for this dataset."),
        h5("Data"),
        checkboxGroupInput(ns("s_blocks"), "Omics blocks", choices = blocks, selected = blocks),
        hr(),
        h5("Preprocessing"),
        checkboxInput(ns("s_standardize"), "Standardize each block (z-score)", value = TRUE),
        p(class = "submodule-desc", "Missing values (per block, matched samples):"),
        uiOutput(ns("s_missing_ui")),
        hr(),
        h5("K - Nearest neighbors"),
        numericInput(ns("s_k"), NULL, value = k_range$default, min = 2, max = k_range$max, step = 1),
        p(class = "submodule-desc", sprintf("Must be smaller than the number of matched samples (%d). Feasible range: %d–%d.", n, k_range$min, k_range$max)),
        checkboxInput(ns("s_k_auto"), "Auto-tune K instead (grid search, scored by fused-network silhouette - slower)", value = FALSE),
        hr(),
        h5("Alpha - Similarity scaling"),
        sliderInput(ns("s_alpha"), NULL, min = MI_SNF_ALPHA_RANGE$min, max = MI_SNF_ALPHA_RANGE$max, value = MI_SNF_ALPHA_RANGE$default, step = 0.05),
        p(class = "submodule-desc", "Scaling parameter for similarity-network construction."),
        checkboxInput(ns("s_alpha_auto"), "Auto-tune Alpha instead (grid search - slower)", value = FALSE),
        hr(),
        h5("T - Fusion iterations"),
        numericInput(ns("s_t"), NULL, value = MI_SNF_T_CANDIDATES[2], min = 5, max = 100, step = 5),
        p(class = "submodule-desc", "Number of network-fusion iterations."),
        checkboxInput(ns("s_t_auto"), sprintf("Auto-converge T instead (stops early once the fused network stabilizes, searched %d–%d)", MI_SNF_T_CANDIDATES[1], max(MI_SNF_T_CANDIDATES)), value = FALSE),
        hr(),
        h5("Number of clusters"),
        numericInput(ns("s_n_clusters"), "Number of clusters", value = 2, min = 2, max = min(6, n - 1), step = 1),
        checkboxInput(ns("s_cluster_auto"), sprintf("Auto-estimate instead (eigengap on the fused network, searched 2–%d)", min(6, n - 1)), value = FALSE),
        hr(),
        h5("Clustering technique"),
        selectInput(ns("s_cluster_method"), NULL, choices = MI_SNF_CLUSTER_METHODS, selected = "spectral"),
        p(class = "submodule-desc", "Partitions the same fused network; does not affect SNF fusion itself."),
        hr(),
        h5("Reproducibility"),
        ## SNFtool::spectralClustering() (and cluster::pam for the "PAM"
        ## technique) are k-means-based internally, so identical K/Alpha/T
        ## settings can still yield different cluster labels run to run
        ## without a fixed seed - the same reproducibility gap DIABLO's own
        ## "Random seed" field (d_seed above) already closes for
        ## mixOmics::tune.block.splsda()/perf().
        numericInput(ns("s_seed"), "Random seed", value = 1, min = 1),
        tags$details(tags$summary("Advanced: network diagnostics"),
                      checkboxInput(ns("s_show_diagnostics"), "Show individual/fused network diagnostics after running", value = TRUE)),
        hr(),
        uiOutput(ns("s_summary_pre")),
        if (isTRUE(elig$ok)) div(
          if (isTRUE(input$s_k_auto) || isTRUE(input$s_alpha_auto) || isTRUE(input$s_t_auto)) mi_warn("Auto-tuning K/alpha/T (checked above) may take a little longer on high-dimensional datasets.") else NULL,
          actionButton(ns("s_run_btn"), "Run SNF", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")),
        uiOutput(ns("s_run_status_ui"))
      )
    })

    output$s_missing_ui <- renderUI({
      v <- req(val_for_blocks(input$s_blocks))
      rows <- lapply(names(v$per_block), function(nm) {
        b <- v$per_block[[nm]]
        sprintf("%s: %.1f%%", nm, if (isTRUE(b$ok)) b$pct_missing else NA)
      })
      any_missing <- any(vapply(v$per_block, function(b) isTRUE(b$ok) && b$n_missing > 0, logical(1)))
      tagList(
        tags$ul(lapply(rows, tags$li)),
        if (any_missing) mi_stop("SNF requires complete data - resolve missing values on the Dataset Workspace tab first.")
        else mi_ok("No missing values detected in the selected blocks.")
      )
    })

    output$s_summary_pre <- renderUI({
      v <- req(val_for_blocks(input$s_blocks))
      tags$div(class = "submodule-desc", tags$strong("Planned run: "), sprintf(
        "Blocks: %s | Samples: %d | K: %s | Alpha: %s | T: %s | Clusters: %s | Technique: %s | Seed: %s",
        paste(input$s_blocks, collapse = " + "), v$n_shared,
        if (isTRUE(input$s_k_auto)) sprintf("Auto-tuned (starting near %s)", input$s_k) else input$s_k,
        if (isTRUE(input$s_alpha_auto)) sprintf("Auto-tuned (starting near %.2f)", input$s_alpha %||% MI_SNF_ALPHA_RANGE$default) else input$s_alpha,
        if (isTRUE(input$s_t_auto)) sprintf("Auto-converged (up to %d)", max(MI_SNF_T_CANDIDATES)) else input$s_t,
        if (isTRUE(input$s_cluster_auto)) "Auto-estimated (eigengap)" else input$s_n_clusters,
        names(MI_SNF_CLUSTER_METHODS)[MI_SNF_CLUSTER_METHODS == (input$s_cluster_method %||% "spectral")],
        input$s_seed %||% 1
      ))
    })

    s_params <- reactive(list(
      standardize = isTRUE(input$s_standardize %||% TRUE),
      k_mode = if (isTRUE(input$s_k_auto)) "automatic" else "manual",
      alpha_mode = if (isTRUE(input$s_alpha_auto)) "automatic" else "manual",
      t_mode = if (isTRUE(input$s_t_auto)) "automatic" else "manual",
      k = input$s_k, alpha = input$s_alpha, t = input$s_t,
      cluster_mode = if (isTRUE(input$s_cluster_auto)) "automatic" else "manual", n_clusters = input$s_n_clusters,
      cluster_method = input$s_cluster_method %||% "spectral",
      ## Threaded straight into a set.seed() call inside mi_snf_run() itself
      ## (multiomics_integration_helpers.R), right before the
      ## kmeans-based SNFtool::spectralClustering()/cluster::pam() calls -
      ## same rationale as DIABLO's own d_seed above.
      seed = input$s_seed %||% 1
    ))

    run_snf <- function() {
      d <- mi_dataset()
      v <- val_for_blocks(input$s_blocks)
      layers <- lapply(d$layers[input$s_blocks], function(m) m[v$shared_ids, , drop = FALSE])
      mi_snf_run(layers, params = s_params())
    }

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      snf_task <- ExtendedTask$new(function(layers, params) promises::future_promise(mi_snf_run(layers, params), seed = TRUE))
      observeEvent(input$s_run_btn, {
        validate(need(isTRUE(snf_elig()$ok), snf_elig()$reason))
        d <- req(mi_dataset()); v <- req(val_for_blocks(input$s_blocks))
        layers <- lapply(d$layers[input$s_blocks], function(m) m[v$shared_ids, , drop = FALSE])
        snf_state$error <- NULL
        snf_state$result <- NULL
        snf_state$submitted <- TRUE
        ## Deferred to the next flush cycle - same rationale as the DIABLO
        ## observer above (snf_task$invoke() blocks this observer
        ## synchronously in this app's actual reactive-domain context).
        p <- s_params()
        session$onFlushed(function() snf_task$invoke(layers, p), once = TRUE)
      })
      observe({
        res <- tryCatch(snf_task$result(), error = function(e) e)
        if (inherits(res, "shiny.silent.error")) return()
        snf_state$submitted <- FALSE
        if (inherits(res, "error")) { snf_state$error <- conditionMessage(res); return() }
        if (!isTRUE(res$ok)) { snf_state$error <- res$error; return() } else { snf_state$result <- res }
      })
      snf_running <- reactive(identical(snf_task$status(), "running") || isTRUE(snf_state$submitted))
    } else {
      observeEvent(input$s_run_btn, {
        validate(need(isTRUE(snf_elig()$ok), snf_elig()$reason))
        showNotification("Running SNF synchronously (future/promises not installed).", type = "message", duration = 5)
        res <- run_snf()
        snf_state$error <- if (!isTRUE(res$ok)) res$error else NULL
        if (isTRUE(res$ok)) snf_state$result <- res
      })
      snf_running <- reactive(FALSE)
    }

    ## Same rationale as the DIABLO button observer above - mutate in place,
    ## never recreate via renderUI (that would reset the click counter).
    observe({
      if (isTRUE(snf_running())) {
        shinyjs::disable(ns("s_run_btn"))
        shinyjs::html(ns("s_run_btn"), as.character(tagList(icon("spinner", class = "fa-spin"), " Running SNF...")))
      } else {
        shinyjs::enable(ns("s_run_btn"))
        shinyjs::html(ns("s_run_btn"), as.character(tagList(icon("play"), " Run SNF")))
      }
    })

    ## Same rationale as d_run_status_ui above - visible right next to the
    ## "Run SNF" button itself.
    output$s_run_status_ui <- renderUI({
      if (!isTRUE(input$s_run_btn > 0)) return(NULL)
      if (isTRUE(snf_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running SNF - see the results panel once this finishes."))
      if (!is.null(snf_state$error)) return(mi_stop(snf_state$error))
      if (!is.null(snf_state$result)) return(mi_ok("Finished - see the results panel to the right."))
      NULL
    })

    output$snf_results_ui <- renderUI({
      if (!isTRUE(input$s_run_btn > 0)) return(multi_empty_state("Set parameters and click \"Run SNF\" to see results here."))
      if (isTRUE(snf_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running SNF - the rest of the app stays usable while this runs."))
      if (!is.null(snf_state$error)) return(mi_stop(snf_state$error))
      res <- snf_state$result
      if (is.null(res)) return(multi_empty_state())
      p <- res$params
      cl_tab <- table(res$clusters)

      tagList(
        box(width = NULL, title = "Actual parameters used", status = "primary", solidHeader = FALSE,
            tags$ul(
              tags$li(sprintf("K: %d (%s)", p$k, p$k_mode)), tags$li(sprintf("Alpha: %.2f (%s)", p$alpha, p$alpha_mode)), tags$li(sprintf("T: %d (%s)", p$t, p$t_mode)),
              tags$li(sprintf("Clusters: %d (%s)", p$n_clusters, p$cluster_mode)),
              tags$li(sprintf("Clustering technique: %s", names(MI_SNF_CLUSTER_METHODS)[MI_SNF_CLUSTER_METHODS == (p$cluster_method %||% "spectral")])),
              tags$li(sprintf("Standardized: %s", p$standardize)),
              tags$li(sprintf("Random seed: %s", p$seed %||% "-"))
            )),
        box(width = NULL, title = "Fused network", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mi_snf_fused_heatmap(res$W, res$clusters), ns("s_heatmap"), height = "420px")),
        box(width = NULL, title = "Clusters", status = "primary", solidHeader = FALSE,
            div(style = "display:flex; gap:10px; flex-wrap:wrap;", lapply(names(cl_tab), function(cl) mi_stat_card(cl_tab[[cl]], sprintf("Cluster %s", cl)))),
            DT::dataTableOutput(ns("s_assign_table")),
            div(class = "table-toolbar", downloadButton(ns("s_dl_assign"), "Download cluster assignments (CSV)", class = "btn-sm"))),
        box(width = NULL, title = "Cluster visualization", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mi_snf_pca_cluster_plot(mi_dataset()$layers[input$s_blocks], res$clusters), ns("s_pca_plot"), height = "380px")),
        box(width = NULL, title = "Cluster quality", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mi_snf_cluster_estimate_plot(res$cluster_estimate), ns("s_estimate_plot"), height = "260px")),
        if (isTRUE(input$s_show_diagnostics)) box(width = NULL, title = "Network concordance (with fused network)", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Network-agreement metric between each block's clustering and the fused clustering, not predictive accuracy."),
            DT::dataTableOutput(ns("s_concordance_table"))),
        uiOutput(ns("s_posthoc_ui"))
      )
    })

    ## multi_plot_or_empty() (in snf_results_ui above) only decides whether to
    ## draw the plotOutput placeholder or an empty-state message - it does not
    ## register the render logic itself (same caveat noted above the d_* DIABLO
    ## renderPlot bindings). These three were missing their bindings, which is
    ## why Fused network / Cluster visualization / Cluster quality rendered as
    ## blank boxes even after a successful SNF run.
    output$s_heatmap <- renderPlot({
      res <- req(snf_state$result)
      mi_snf_fused_heatmap(res$W, res$clusters)
    })
    output$s_pca_plot <- renderPlot({
      res <- req(snf_state$result)
      mi_snf_pca_cluster_plot(mi_dataset()$layers[input$s_blocks], res$clusters)
    })
    output$s_estimate_plot <- renderPlot({
      res <- req(snf_state$result)
      mi_snf_cluster_estimate_plot(res$cluster_estimate)
    })

    output$s_assign_table <- DT::renderDataTable({
      res <- req(snf_state$result)
      DT::datatable(data.frame(sample = names(res$clusters), cluster = as.integer(res$clusters)), rownames = FALSE, options = list(pageLength = 15), class = "stripe hover compact")
    })
    output$s_dl_assign <- downloadHandler(function() "snf_cluster_assignments.csv", function(file) {
      res <- req(snf_state$result)
      utils::write.csv(data.frame(sample = names(res$clusters), cluster = as.integer(res$clusters)), file, row.names = FALSE)
    })
    output$s_concordance_table <- DT::renderDataTable({
      res <- req(snf_state$result)
      conc <- mi_snf_concordance(res)
      if (is.null(conc)) return(DT::datatable(data.frame(Note = "Not available for this configuration.")))
      DT::datatable(conc, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    snf_posthoc <- reactive({
      res <- req(snf_state$result); d <- mi_dataset(); v <- mi_val()
      req(!is.null(d$sample_meta), !is.null(input$outcome_col), input$outcome_col %in% colnames(d$sample_meta))
      outcome <- stats::setNames(d$sample_meta[[input$outcome_col]], rownames(d$sample_meta))
      mi_snf_posthoc_outcome(res, outcome, v$shared_ids)
    })

    output$s_posthoc_ui <- renderUI({
      res <- req(snf_state$result); d <- mi_dataset()
      if (is.null(d$sample_meta) || is.null(input$outcome_col) || !input$outcome_col %in% colnames(d$sample_meta)) {
        return(box(width = NULL, title = "Post-hoc outcome evaluation", status = "primary", solidHeader = FALSE, multi_empty_state("No outcome variable selected - post-hoc evaluation is optional and skipped.")))
      }
      ph <- tryCatch(snf_posthoc(), error = function(e) NULL)
      box(width = NULL, title = "Post-hoc outcome evaluation", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", tags$strong("Post-hoc evaluation, not SNF prediction accuracy"), " - clusters were formed without the outcome; this checks whether they line up with it."),
          if (is.null(ph)) mi_warn("Not enough matched samples with a valid outcome for this comparison.") else tagList(
            verbatimTextOutput(ns("s_posthoc_table")),
            div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                mi_stat_card(sprintf("%.3f", ph$nmi), "NMI (cluster vs. outcome)"),
                mi_stat_card(sprintf("%.3f", ph$ari), "Adjusted Rand Index"),
                mi_stat_card(sprintf("%.3g", ph$fisher_p), "Fisher's exact p"))
          ))
    })
    output$s_posthoc_table <- renderPrint(print(req(snf_posthoc())$table))

    ## =========================================================================
    ## 5. Compare subtab (spec sections 29-30) - reads whatever DIABLO/SNF
    ## results have already been computed above; never recomputes them.
    ## =========================================================================
    output$compare_ui <- renderUI({
      tagList(
        box(width = NULL, title = "Supervised: DIABLO vs. single-omics", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Single-omics baselines: elastic-net per view, same samples and outcome as the DIABLO run."),
            if (is.null(diablo_state$result)) multi_empty_state("Run DIABLO first (DIABLO tab).")
            else tagList(actionButton(ns("c_run_supervised"), "Run Comparison", icon = icon("play"), class = "btn-primary btn-sm"), uiOutput(ns("c_supervised_ui")))
        ),
        box(width = NULL, title = "Unsupervised: SNF vs. single-omics clustering", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Per-block clustering vs. fused SNF clustering, compared by NMI/ARI - agreement metrics, not accuracy."),
            if (is.null(snf_state$result)) multi_empty_state("Run SNF first (SNF tab).")
            else tagList(actionButton(ns("c_run_unsupervised"), "Run Comparison", icon = icon("play"), class = "btn-primary btn-sm"), uiOutput(ns("c_unsupervised_ui")))
        )
      )
    })

    c_supervised <- eventReactive(input$c_run_supervised, {
      ## The outcome selector is read live at click time, but the fitted
      ## DIABLO result was trained on whatever outcome was selected back
      ## when "Run DIABLO" was clicked (diablo_state$outcome_col_used) - if
      ## the user has since switched the outcome dropdown, comparing against
      ## the CURRENT selection would silently pit DIABLO's real performance
      ## against a single-omics baseline fit on a different label, and
      ## present the mismatch as a head-to-head result. Returned as the same
      ## fail-soft list(ok, error) shape mi_compare_supervised() itself uses
      ## (not validate()) - output$c_supervised_ui below wraps this call in
      ## tryCatch(..., error = function(e) NULL), which would silently
      ## swallow a validate() condition into a blank panel instead of
      ## showing the warning.
      if (!identical(input$outcome_col, diablo_state$outcome_col_used)) {
        return(list(ok = FALSE, error = sprintf(
          "The outcome variable has changed since DIABLO was run (fit on \"%s\", now \"%s\" is selected). Re-run DIABLO on the current outcome before comparing.",
          diablo_state$outcome_col_used %||% "-", input$outcome_col %||% "-")))
      }
      d <- mi_dataset(); v <- mi_val()
      outcome <- stats::setNames(d$sample_meta[[input$outcome_col]], rownames(d$sample_meta))
      mi_compare_supervised(d$layers[diablo_state$result$params$blocks], outcome, v$shared_ids, diablo_state$result)
    }, ignoreInit = TRUE)
    output$c_supervised_ui <- renderUI({
      res <- tryCatch(c_supervised(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      if (!isTRUE(res$ok)) return(mi_stop(res$error))
      tagList(
        multi_plot_or_empty(function() mi_compare_bar_plot(res$table), ns("c_sup_plot"), height = "260px"),
        DT::dataTableOutput(ns("c_sup_table")),
        p(class = "submodule-desc", tags$em(res$note))
      )
    })
    ## Same missing-renderPlot bug as the SNF plots above - the plotOutput
    ## placeholder existed (260px tall) but nothing ever drew into it, which
    ## is exactly the "long blank space" between the Run Comparison button
    ## and the results table.
    output$c_sup_plot <- renderPlot({
      mi_compare_bar_plot(req(c_supervised())$table)
    })
    output$c_sup_table <- DT::renderDataTable({
      DT::datatable(req(c_supervised())$table, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    c_unsupervised <- eventReactive(input$c_run_unsupervised, {
      mi_compare_unsupervised(snf_state$result)
    }, ignoreInit = TRUE)
    output$c_unsupervised_ui <- renderUI({
      res <- tryCatch(c_unsupervised(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      if (!isTRUE(res$ok)) return(mi_stop(res$error))
      DT::dataTableOutput(ns("c_unsup_table"))
    })
    output$c_unsup_table <- DT::renderDataTable({
      DT::datatable(req(c_unsupervised())$table, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    ## =========================================================================
    ## 6. Sex-Stratified subtab - exact-match port of
    ## Research_05_multiomics_sexstratified's own nested-CV pipeline
    ## (multiomics_sexstratified_engine.R::mss_run_stratified()), covering
    ## both DIABLO and Random Forest, run per sex stratum on whichever
    ## dataset is selected above (Active or preloaded cell). Kept separate
    ## from the DIABLO tab above: this engine's feature-selection/CV design
    ## is fixed (in-fold covariate adjustment, specific top-K/keepX), not the
    ## freeform tune.block.splsda()-based engine DIABLO users already rely
    ## on - this is the one place a user's own uploaded data gets the same
    ## leakage-safe analysis that produced Table34/35/37/39.
    ## =========================================================================
    ss_state <- reactiveValues(result = NULL, error = NULL, submitted = FALSE)

    ss_sex_col <- reactive({
      d <- mi_dataset()
      if (!isTRUE(d$ok) || is.null(d$sample_meta)) return(NULL)
      cands <- multi_sex_candidates(d$sample_meta)
      if (length(cands) == 0) NULL else cands[1]
    })

    ss_covariate_choices <- reactive({
      d <- mi_dataset()
      if (!isTRUE(d$ok) || is.null(d$sample_meta)) return(character(0))
      setdiff(colnames(d$sample_meta), c(input$outcome_col, ss_sex_col()))
    })

    output$ss_block_role_ui <- renderUI({
      v <- req(mi_val())
      blocks <- v$block_labels
      if (length(blocks) < 2) return(mi_warn("This engine needs at least two omics blocks (expression + methylation)."))
      guess_e <- Find(function(nm) grepl("transcript|express|rna", nm, ignore.case = TRUE), blocks) %||% blocks[1]
      guess_m <- Find(function(nm) grepl("methyl", nm, ignore.case = TRUE), blocks) %||% blocks[min(2, length(blocks))]
      tagList(
        selectInput(ns("ss_expr_block"), "Expression/transcriptomics block", choices = blocks, selected = guess_e),
        selectInput(ns("ss_meth_block"), "Methylation block", choices = blocks, selected = guess_m)
      )
    })

    output$ss_params_ui <- renderUI({
      d <- mi_dataset(); v <- mi_val()
      if (!isTRUE(d$ok) || is.null(v) || !isTRUE(v$ok)) return(box(width = NULL, title = "Sex-Stratified parameters", status = "primary", solidHeader = FALSE, mi_warn("Select a valid dataset above first.")))
      sex_col <- ss_sex_col()
      cov_choices <- ss_covariate_choices()
      cov_guess <- cov_choices[grepl("drug|treatment|therapy", cov_choices, ignore.case = TRUE)][1]

      box(
        width = NULL, title = "Sex-Stratified parameters", status = "primary", solidHeader = FALSE,
        if (is.null(sex_col)) mi_warn("No Sex/Gender column detected - \"All (pooled)\" still works; Female/Male splitting is unavailable.")
        else mi_ok(sprintf("Sex column detected: \"%s\".", sex_col)),
        hr(),
        h5("Omics blocks"),
        uiOutput(ns("ss_block_role_ui")),
        hr(),
        h5("Engine"),
        radioButtons(ns("ss_engine"), NULL, choices = MSS_ENGINE_CHOICES, selected = "diablo"),
        hr(),
        h5("Sex stratification"),
        radioButtons(ns("ss_sex_mode"), NULL, choices = MSS_SEX_MODE_CHOICES, selected = "pooled"),
        hr(),
        h5("Adjustment covariate"),
        selectInput(ns("ss_covariate"), NULL, choices = c("(none)" = "", stats::setNames(cov_choices, cov_choices)), selected = cov_guess %||% ""),
        p(class = "submodule-desc", "Adjusted for during in-fold feature selection (e.g. drug), matching the pipeline's own design - optional."),
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
      d <- req(mi_dataset())
      mss_run_stratified(
        d$layers[[input$ss_expr_block]], d$layers[[input$ss_meth_block]], d$sample_meta, input$outcome_col,
        covariate_col = if (nzchar(input$ss_covariate %||% "")) input$ss_covariate else NULL,
        sex_mode = input$ss_sex_mode %||% "pooled", engine = input$ss_engine %||% "diablo", params = ss_params()
      )
    }

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      ss_task <- ExtendedTask$new(function(expr, meth, meta, outcome_col, covariate_col, sex_mode, engine, params) {
        promises::future_promise(mss_run_stratified(expr, meth, meta, outcome_col, covariate_col, sex_mode, engine, params), seed = TRUE)
      })
      observeEvent(input$ss_run_btn, {
        d <- req(mi_dataset())
        req(input$outcome_col, input$ss_expr_block, input$ss_meth_block)
        sex_mode <- input$ss_sex_mode %||% "pooled"
        ## A sex column is only required when actually splitting by sex -
        ## "All (pooled)" (the default, and the only sensible mode for an
        ## already-single-sex preloaded cell) runs with no sex column at all.
        validate(need(identical(sex_mode, "pooled") || !is.null(ss_sex_col()), "No Sex/Gender column detected in this dataset's metadata - switch Sex stratification to \"All (pooled)\", or use a dataset with a Sex/Gender column."))
        expr <- d$layers[[input$ss_expr_block]]; meth <- d$layers[[input$ss_meth_block]]
        meta <- d$sample_meta; outcome_col <- input$outcome_col
        covariate_col <- if (nzchar(input$ss_covariate %||% "")) input$ss_covariate else NULL
        engine <- input$ss_engine %||% "diablo"
        p <- ss_params()
        ss_state$error <- NULL; ss_state$result <- NULL; ss_state$submitted <- TRUE
        ## Deferred to the next flush cycle - same rationale as the DIABLO/SNF
        ## observers above (invoke() blocks this observer synchronously).
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
        validate(need(identical(sex_mode, "pooled") || !is.null(ss_sex_col()), "No Sex/Gender column detected in this dataset's metadata - switch Sex stratification to \"All (pooled)\", or use a dataset with a Sex/Gender column."))
        showNotification("Running Sex-Stratified analysis synchronously (future/promises not installed) - the app will be briefly unresponsive.", type = "message", duration = 5)
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

      tagList(
        box(width = NULL, title = "Performance by stratum (nested cross-validated AUROC)", status = "primary", solidHeader = FALSE,
            div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px;",
                lapply(seq_len(nrow(perf)), function(i) mi_stat_card(
                  sprintf("%.3f [%.3f, %.3f]", perf$auroc[i], perf$ci_lo[i], perf$ci_hi[i]),
                  sprintf("%s (n=%d)%s", perf$stratum[i], perf$n[i], if (isTRUE(perf$excludes_chance[i])) " *" else ""),
                  if (isTRUE(perf$excludes_chance[i])) ARTHOMIX_COLORS$aqua else ARTHOMIX_COLORS$blue))),
            p(class = "submodule-desc", "* excludes chance (95% CI does not cross 0.5) - the same criterion Table34/37/39 use."),
            DT::dataTableOutput(ns("ss_perf_table")),
            div(class = "table-toolbar", downloadButton(ns("ss_dl_perf"), "Download performance (CSV)", class = "btn-sm"))),
        box(width = NULL, title = if (identical(res$engine, "rf")) "Feature panel (Random Forest importance)" else "Feature panel (DIABLO loadings)", status = "primary", solidHeader = FALSE,
            if (!is.null(res$panel_note)) mi_warn(res$panel_note),
            if (is.null(panels) || nrow(panels) == 0) mi_warn("No panel available.") else tagList(
              DT::dataTableOutput(ns("ss_panel_table")),
              div(class = "table-toolbar", downloadButton(ns("ss_dl_panel"), "Download panel (CSV)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Biomarker comparison by sex", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "One row per candidate feature - a value in the Female/Male/Pooled column means that feature was selected in that stratum's own full-cohort panel; a blank cell means it wasn't."),
            if (is.null(res$panels_wide) || nrow(res$panels_wide) == 0) mi_warn("Not available - run \"Female and Male separately\" (or compare separate pooled/female/male runs) to populate this comparison.") else tagList(
              DT::dataTableOutput(ns("ss_panel_wide_table")),
              div(class = "table-toolbar", downloadButton(ns("ss_dl_panel_wide"), "Download comparison (CSV)", class = "btn-sm"))
            )),
        if (length(skipped) > 0) box(width = NULL, title = "Strata skipped", status = "primary", solidHeader = FALSE,
            tags$ul(lapply(names(skipped), function(nm) tags$li(sprintf("%s: %s", nm, skipped[[nm]]$error)))))
      )
    })

    output$ss_perf_table <- DT::renderDataTable({
      DT::datatable(req(ss_state$result)$performance, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$ss_dl_perf <- downloadHandler(function() "sex_stratified_performance.csv", function(file) utils::write.csv(req(ss_state$result)$performance, file, row.names = FALSE))
    output$ss_panel_table <- DT::renderDataTable({
      DT::datatable(req(ss_state$result)$panels, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$ss_dl_panel <- downloadHandler(function() "sex_stratified_panel.csv", function(file) utils::write.csv(req(ss_state$result)$panels, file, row.names = FALSE))
    output$ss_panel_wide_table <- DT::renderDataTable({
      DT::datatable(req(ss_state$result)$panels_wide, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$ss_dl_panel_wide <- downloadHandler(function() "sex_stratified_biomarker_comparison.csv", function(file) utils::write.csv(req(ss_state$result)$panels_wide, file, row.names = FALSE))

    ## =========================================================================
    ## Publish - kept cell$label-shaped so multi_qc_scorecard()/
    ## multi_analysis_summary_table() (multiomics_helpers.R) keep working
    ## unmodified for either data source.
    ## =========================================================================
    observe({
      if (is.null(multi_results)) return()
      d <- tryCatch(mi_dataset(), error = function(e) NULL)
      if (is.null(d) || !isTRUE(d$ok)) return()
      multi_results$integration <- list(
        cell = list(label = d$label),
        perf = if (!is.null(diablo_state$result)) mi_diablo_performance_summary(diablo_state$result) else NULL,
        snf_perf = if (!is.null(snf_state$result)) snf_state$result$params else NULL
      )
      multi_results$integration_stratified <- if (!is.null(ss_state$result)) list(cell = list(label = d$label), result = ss_state$result) else NULL
    })
  })
}
