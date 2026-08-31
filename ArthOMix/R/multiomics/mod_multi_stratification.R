## R/multiomics/mod_multi_stratification.R
## Submodule: SNF Clustering - a live, data-adaptive unsupervised patient
## stratification workflow (Similarity Network Fusion, Wang et al. 2014,
## SNFtool::SNF) on either the Active Multi-Omics Dataset (Dataset Workspace
## tab) or a preloaded RA anti-TNF analysis cell. Reuses the tested SNF
## engine already built for the Multi-omics Integration submodule
## (mi_snf_*()/mi_ari()/mi_snf_concordance() - multiomics_integration_live_
## helpers.R) rather than re-implementing it; this file's own helpers
## (snf_clustering_helpers.R) add only what that engine does not already
## cover for a dedicated stratification workflow: preprocessing, single-omics
## fallback, clinical association testing, resampling-based stability,
## parameter sensitivity, and per-feature cluster-association ranking.
##
## Workflow (spec): Data -> Validate -> Configure SNF -> Run -> Clusters ->
## Stability -> Clinical -> Features. Nothing under Clusters/Stability/
## Clinical/Features renders until the blue "Run SNF Clustering" button is
## clicked - every ..._ui below gates on `input$run_btn > 0` first, exactly
## like Multi-omics Integration's own diablo_results_ui/snf_results_ui.
## Clusters are always unsupervised molecular groups; any clinical/outcome
## comparison is explicitly post-hoc and never used to choose K/alpha/T or
## the cluster count.

mod_multi_stratification_config <- list(
  id = "stratification", title = "SNF Clustering", icon = "diagram-project", group = "Data",
  description = "Unsupervised patient clustering via Similarity Network Fusion (SNF)."
)

## ---------------------------------------------------------------------------
## Small shared UI pieces (mi_stat_card()/mi_warn()/mi_ok()/mi_stop() are
## already defined, globally, in mod_multi_integration.R - reused as-is
## rather than redefined, so both submodules read identically).
## ---------------------------------------------------------------------------

mod_multi_stratification_ui <- function(id) {
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
                       choices = c("Active Multi-Omics Dataset (Dataset Workspace)" = "active", "Preloaded RA anti-TNF cohort" = "preloaded"),
                       selected = "active", inline = TRUE),
          conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
                            selectInput(ns("preloaded_cell"), "Analysis cell", choices = MULTI_CELL_CHOICES, width = "100%")),
          uiOutput(ns("source_note"))
        ),
        uiOutput(ns("blocks_ui")),
        uiOutput(ns("validation_ui")),
        uiOutput(ns("preproc_ui")),
        uiOutput(ns("ready_ui"))
      ),
      tabPanel("SNF Setup", br(), uiOutput(ns("setup_ui"))),
      tabPanel("Clusters", br(), uiOutput(ns("clusters_ui"))),
      tabPanel("Stability", br(), uiOutput(ns("stability_ui"))),
      tabPanel("Clinical", br(), uiOutput(ns("clinical_ui"))),
      tabPanel("Features", br(), uiOutput(ns("features_ui")))
    )
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_multi_stratification_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    blk_id <- function(prefix, b) paste0(prefix, "_", make.names(b))

    output$active_dataset_banner <- renderUI(multi_active_dataset_banner(multi_dataset))

    ## =========================================================================
    ## 1. Data selection (spec section 25) - one reactive dataset object,
    ## converging on the same shape mi_dataset() uses in Multi-omics
    ## Integration, never merging preloaded and Active Multi-Omics Dataset.
    ## =========================================================================
    sc_dataset <- reactive({
      if (identical(input$data_source, "preloaded")) {
        req(input$preloaded_cell)
        sfc_preloaded_dataset(input$preloaded_cell)
      } else {
        sfc_active_dataset(multi_dataset)
      }
    })

    output$source_note <- renderUI({
      d <- sc_dataset()
      if (!isTRUE(d$ok)) return(mi_warn(d$error))
      if (!MULTI_SNF_LIVE_AVAILABLE) return(mi_stop("SNFtool is not installed in this deployment - SNF Clustering is unavailable."))
      mi_ok(d$provenance)
    })

    output$blocks_ui <- renderUI({
      d <- sc_dataset()
      if (!isTRUE(d$ok) || length(d$layers) == 0) return(NULL)
      choices <- names(d$layers)
      box(width = NULL, title = "2. Modalities", status = "primary", solidHeader = FALSE,
          checkboxGroupInput(ns("blocks"), "Available modalities (select at least one; at least two for true multi-omics SNF)", choices = choices, selected = choices, inline = TRUE))
    })

    ## =========================================================================
    ## 2. Validation (spec sections 6-8) - raw block QC + sample matching,
    ## before any preprocessing is applied. sfc_validate_dataset() handles
    ## both the single- and multi-block case (never fails just because one
    ## modality was deselected).
    ## =========================================================================
    sc_val_raw <- reactive({
      d <- sc_dataset()
      if (!isTRUE(d$ok)) return(NULL)
      blocks <- intersect(input$blocks %||% names(d$layers), names(d$layers))
      if (length(blocks) == 0) return(NULL)
      sfc_validate_dataset(d$layers[blocks], d$sample_meta)
    })

    output$validation_ui <- renderUI({
      v <- sc_val_raw()
      if (is.null(v)) return(NULL)
      if (!isTRUE(v$ok)) return(box(width = NULL, title = "3. Data validation", status = "primary", solidHeader = FALSE, mi_warn(v$error)))
      block_cards <- lapply(names(v$per_block), function(nm) {
        b <- v$per_block[[nm]]
        if (!isTRUE(b$ok)) return(mi_stat_card("-", nm, ARTHOMIX_COLORS$red))
        div(class = "card", style = "flex:1 1 220px; padding:12px;",
            div(style = "font-weight:600; margin-bottom:4px;", nm),
            div(style = "font-size:0.85em;", sprintf("Samples: %s | Features: %s", format(b$n_samples, big.mark = ","), format(b$n_features, big.mark = ","))),
            div(style = "font-size:0.85em;", sprintf("Missing: %.1f%% | Zero-variance: %d | Dup. samples: %d | Dup. features: %d", b$pct_missing, b$n_zero_variance, b$n_duplicate_samples, b$n_duplicate_features)))
      })
      n_excluded <- if (length(v$per_block) >= 1) max(vapply(v$per_block, function(b) if (isTRUE(b$ok)) b$n_samples else 0, integer(1))) - v$n_shared else 0
      box(
        width = NULL, title = "3. Data validation", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Patient/sample identifiers are matched by rowname intersection, not row position."),
        div(style = "display:flex; gap:12px; flex-wrap:wrap; margin-bottom:10px;", block_cards),
        div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px;",
            mi_stat_card(v$n_blocks, "Modalities selected"),
            mi_stat_card(v$n_shared, "Matched patients", ARTHOMIX_COLORS$aqua),
            mi_stat_card(max(0, n_excluded), "Excluded (not present in every selected modality)", ARTHOMIX_COLORS$orange)),
        if (!isTRUE(v$reliable_matching)) mi_stop(v$mismatch_message) else mi_ok(sprintf("%d matched patients will be used.", v$n_shared))
      )
    })

    ## =========================================================================
    ## 3. Preprocessing (spec sections 9-11) - per selected block: data-type
    ## detection drives which transforms are even offered; missing-value
    ## handling is never automatic; feature filtering is optional and never
    ## forces blocks to share a feature count.
    ## =========================================================================
    output$preproc_ui <- renderUI({
      d <- sc_dataset(); v <- sc_val_raw()
      if (!isTRUE(d$ok) || is.null(v) || !isTRUE(v$ok)) return(NULL)
      box(
        width = NULL, title = "4. Preprocessing", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Order: missing values, then transform, then feature filter - applied per modality."),
        tagList(lapply(v$block_labels, function(b) {
          mat <- d$layers[[b]]
          mk <- make.names(b)
          omics_type <- tryCatch(d$layer_meta[[b]]$omics_type, error = function(e) NULL)
          tc <- sfc_transform_choices(mat, omics_type)
          n_missing <- sum(is.na(mat))
          tags$details(
            open = TRUE,
            tags$summary(sprintf("%s (%d x %d, detected: %s)", b, nrow(mat), ncol(mat), tc$type)),
            p(class = "submodule-desc", tc$note),
            fluidRow(
              column(4, selectInput(ns(blk_id("transform", b)), "Transformation", choices = tc$choices, selected = "none")),
              column(4, selectInput(ns(blk_id("missing", b)), "Missing values",
                                     choices = c("None (fails if any remain)" = "none", "Remove features above threshold" = "remove_cols",
                                                 "Remove samples above threshold" = "remove_rows", "Impute mean" = "mean", "Impute median" = "median"),
                                     selected = "none")),
              column(4, numericInput(ns(blk_id("misspct", b)), "Missingness threshold (%)", value = 50, min = 0, max = 100, step = 5))
            ),
            if (n_missing > 0) mi_warn(sprintf("%d missing value(s) (%.1f%%) in this block - choose a missing-value option above.", n_missing, 100 * n_missing / (nrow(mat) * ncol(mat)))) else mi_ok("No missing values detected."),
            fluidRow(
              column(6, radioButtons(ns(blk_id("filtcrit", b)), "Feature filter", choices = c("None (use all features)" = "none", "Top-N by variance" = "variance", "Top-N by MAD" = "mad"), selected = "none", inline = TRUE)),
              column(6, numericInput(ns(blk_id("filtn", b)), "Features to keep (if filtering)", value = min(2000, ncol(mat)), min = 10, max = ncol(mat), step = 10))
            ),
            if (ncol(mat) > 3000) mi_warn(sprintf("%s features is high-dimensional - consider a variance/MAD filter for computational feasibility.", format(ncol(mat), big.mark = ","))),
            hr()
          )
        }))
      )
    })

    ## Applies the chosen preprocessing to each selected block's own full
    ## sample set (never restricted to the matched set first) - matching is
    ## re-evaluated afterward (sc_val2()) since removing high-missingness
    ## samples can change which patients remain.
    sc_ready <- reactive({
      d <- req(sc_dataset()); v <- req(sc_val_raw())
      out <- list(); errors <- character(0)
      for (b in v$block_labels) {
        mat <- d$layers[[b]]
        pp <- sfc_preprocess_block(
          mat,
          transform = input[[blk_id("transform", b)]] %||% "none",
          missing_method = input[[blk_id("missing", b)]] %||% "none",
          missing_threshold = input[[blk_id("misspct", b)]] %||% 50,
          filter_criterion = input[[blk_id("filtcrit", b)]] %||% "none",
          filter_top_n = input[[blk_id("filtn", b)]]
        )
        if (!isTRUE(pp$ok)) { errors <- c(errors, sprintf("%s: %s", b, pp$error)); next }
        out[[b]] <- pp$mat
      }
      list(layers = out, errors = errors)
    })

    sc_val2 <- reactive({
      r <- sc_ready()
      if (length(r$layers) == 0) return(NULL)
      sfc_validate_dataset(r$layers, sc_dataset()$sample_meta)
    })

    sc_elig <- reactive({
      r <- sc_ready()
      if (length(r$errors) > 0) return(list(ok = FALSE, reason = paste(r$errors, collapse = " "), mode = NULL))
      if (!MULTI_SNF_LIVE_AVAILABLE) return(list(ok = FALSE, reason = "SNFtool is not installed in this deployment.", mode = NULL))
      sfc_eligibility(sc_val2())
    })

    output$ready_ui <- renderUI({
      v2 <- sc_val2(); e <- sc_elig()
      if (is.null(v2)) return(NULL)
      box(width = NULL, title = "5. Ready-to-run status", status = "primary", solidHeader = FALSE,
          if (isTRUE(e$ok)) mi_ok(sprintf("Ready - %d matched patients after preprocessing (%s).", v2$n_shared, if (identical(e$mode, "single_omics")) "Single-Omics Clustering" else "multi-omics SNF")) else mi_stop(e$reason))
    })

    ## =========================================================================
    ## 4. SNF Setup (spec sections 12-15) - K/Alpha/T/cluster-count, each
    ## independently Automatic or Manual, ranges always derived from the
    ## matched-sample count in front of them (never one universal value).
    ## =========================================================================
    output$setup_ui <- renderUI({
      v2 <- sc_val2(); e <- sc_elig()
      if (is.null(v2)) return(box(width = NULL, title = "SNF parameters", status = "primary", solidHeader = FALSE, mi_warn("Select and validate a dataset on the Data tab first.")))
      n <- v2$n_shared
      k_range <- mi_snf_feasible_k_range(max(n, 4))
      multi_mode <- identical(e$mode, "multi_omics_snf")

      box(
        width = NULL, title = "SNF parameters", status = "primary", solidHeader = FALSE,
        if (!isTRUE(e$ok)) mi_stop(e$reason) else mi_ok(sprintf("%s is available for this dataset (%d matched patients).", if (multi_mode) "SNF" else "Single-Omics Clustering", n)),
        div(style = if (!isTRUE(e$ok)) "opacity:0.5; pointer-events:none;" else NULL,
          h5("Scaling"),
          checkboxInput(ns("standardize"), "Standardize each block (z-score) before building the similarity network", value = TRUE),
          p(class = "submodule-desc", "SNFtool::standardNormalization(). Disable if already Autoscaled under Data > Preprocessing."),
          hr(),
          h5("Distance / similarity"),
          selectInput(ns("dist_method"), NULL, choices = c("Squared Euclidean distance (SNFtool::dist2)" = "dist2")),
          p(class = "submodule-desc", "Feeds SNFtool::affinityMatrix() - the only distance metric this engine implements."),
          hr(),
          h5("K - Nearest neighbors"),
          numericInput(ns("k"), "K", value = k_range$default, min = 2, max = k_range$max, step = 1),
          p(class = "submodule-desc", sprintf("Must be smaller than matched patients (%d). Feasible range: %d-%d.", n, k_range$min, k_range$max)),
          checkboxInput(ns("k_auto"), "Auto-tune K instead (grid search, scored by fused-network silhouette - slower)", value = FALSE),
          hr(),
          h5("Alpha - Affinity scaling"),
          sliderInput(ns("alpha"), "Alpha", min = MI_SNF_ALPHA_RANGE$min, max = MI_SNF_ALPHA_RANGE$max, value = MI_SNF_ALPHA_RANGE$default, step = 0.05),
          p(class = "submodule-desc", "Controls affinity scaling."),
          checkboxInput(ns("alpha_auto"), "Auto-tune Alpha instead (grid search - slower)", value = FALSE),
          hr(),
          if (multi_mode) tagList(
            h5("T - Fusion iterations"),
            numericInput(ns("t"), "T", value = MI_SNF_T_CANDIDATES[2], min = 5, max = 100, step = 5),
            p(class = "submodule-desc", "Number of network-fusion diffusion iterations (commonly 20-50)."),
            checkboxInput(ns("t_auto"), sprintf("Auto-converge T instead (stops early once the fused network stabilizes, searched %d-%d)", MI_SNF_T_CANDIDATES[1], max(MI_SNF_T_CANDIDATES)), value = FALSE),
            hr()
          ) else p(class = "submodule-desc", "T does not apply - only one modality selected, no fusion step."),
          h5("Number of clusters"),
          numericInput(ns("n_clusters"), "Number of clusters", value = 2, min = 2, max = min(6, max(n - 1, 2)), step = 1),
          checkboxInput(ns("cluster_auto"), sprintf("Auto-estimate instead (eigengap search, 2-%d)", min(6, max(n - 1, 2))), value = FALSE),
          hr(),
          h5("Clustering technique"),
          selectInput(ns("cluster_method"), NULL, choices = MI_SNF_CLUSTER_METHODS, selected = "spectral"),
          p(class = "submodule-desc", "Partitions the same fused network regardless of technique."),
          hr(),
          h5("Reproducibility"),
          ## SNFtool::spectralClustering() (and cluster::pam for the "PAM"
          ## technique) are k-means-based internally, so identical
          ## K/Alpha/T settings can still yield different cluster
          ## assignments run to run without a fixed seed - same
          ## reproducibility gap Multi-omics Integration's own SNF tab
          ## closes with its "Random seed" field.
          numericInput(ns("seed"), "Random seed", value = 1, min = 1),
          hr(),
          uiOutput(ns("summary_pre")),
          if (isTRUE(e$ok)) div(
            if (isTRUE(input$k_auto) || isTRUE(input$alpha_auto) || isTRUE(input$t_auto)) mi_warn("Auto-tuning plus the default stability check may take longer for larger cohorts.") else NULL,
            actionButton(ns("run_btn"), "Run SNF Clustering", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")),
          uiOutput(ns("run_status_ui"))
        )
      )
    })

    output$summary_pre <- renderUI({
      v2 <- req(sc_val2())
      tags$div(class = "submodule-desc", tags$strong("Planned run: "), sprintf(
        "Modalities: %s | Matched patients: %d | K: %s | Alpha: %s | T: %s | Clusters: %s | Technique: %s | Seed: %s",
        paste(names(sc_ready()$layers), collapse = " + "), v2$n_shared,
        if (isTRUE(input$k_auto)) sprintf("Auto-tuned (starting near %s)", input$k) else input$k,
        if (isTRUE(input$alpha_auto)) sprintf("Auto-tuned (starting near %.2f)", input$alpha %||% MI_SNF_ALPHA_RANGE$default) else input$alpha,
        if (isTRUE(input$t_auto)) sprintf("Auto-converged (up to %d)", max(MI_SNF_T_CANDIDATES)) else (input$t %||% "-"),
        if (isTRUE(input$cluster_auto)) "Auto-estimated (eigengap)" else input$n_clusters,
        names(MI_SNF_CLUSTER_METHODS)[MI_SNF_CLUSTER_METHODS == (input$cluster_method %||% "spectral")],
        input$seed %||% 1
      ))
    })

    sc_params <- reactive(list(
      standardize = isTRUE(input$standardize %||% TRUE),
      k_mode = if (isTRUE(input$k_auto)) "automatic" else "manual",
      alpha_mode = if (isTRUE(input$alpha_auto)) "automatic" else "manual",
      t_mode = if (isTRUE(input$t_auto)) "automatic" else "manual",
      k = input$k, alpha = input$alpha, t = input$t,
      cluster_mode = if (isTRUE(input$cluster_auto)) "automatic" else "manual", n_clusters = input$n_clusters,
      cluster_method = input$cluster_method %||% "spectral",
      ## Threaded into set.seed() inside sfc_snf_run()/mi_snf_run()
      ## (snf_clustering_helpers.R / multiomics_integration_live_helpers.R),
      ## right before their kmeans-based SNFtool::spectralClustering()/
      ## cluster::pam() calls.
      seed = input$seed %||% 1
    ))

    ## =========================================================================
    ## 5. Run (spec section 5, "nothing before the blue button") - the main
    ## run bundles the headline clustering AND a default stability check
    ## (fixed, sensible defaults) in one async task, since clusters must
    ## never be shown/labeled without stability evidence (spec rule 8). The
    ## Stability tab additionally offers its own "Recompute" controls for a
    ## custom resample count/fraction, and a separate parameter-sensitivity
    ## check (both secondary, explicit actions - the heavier of the two,
    ## sensitivity, is optional/best-effort per spec section 19).
    ## =========================================================================
    ## `submitted` is set eagerly, synchronously, in the same observer that
    ## calls run_task$invoke() - future::multisession worker cold-start
    ## (spinning up a background R process the first time, loading SNFtool
    ## in it) can itself take a few real seconds, during which
    ## run_task$status() can still read something other than "running".
    ## Without this flag, every render branch below (gate(), run_status_ui,
    ## the button itself) fell through to nothing visible during that gap -
    ## confirmed directly on the sibling Multi-omics Integration submodule's
    ## identical pattern (DIABLO/SNF tabs, mod_multi_integration.R): an
    ## empty/unchanged panel for several seconds after clicking, not just a
    ## briefly slow spinner. `submitted` closes that gap so the spinner
    ## appears the instant the button is clicked.
    state <- reactiveValues(result = NULL, stability = NULL, error = NULL, submitted = FALSE, layers_used = NULL, sample_meta = NULL, dataset_label = NULL)

    snapshot_inputs <- function() {
      r <- sc_ready(); v2 <- req(sc_val2())
      layers <- lapply(r$layers[input$blocks], function(m) m[v2$shared_ids, , drop = FALSE])
      state$error <- NULL
      state$result <- NULL
      state$submitted <- TRUE
      state$layers_used <- layers
      state$dataset_label <- sc_dataset()$label
      state$sample_meta <- sc_dataset()$sample_meta
      layers
    }

    if (isTRUE(ARTHOMIX_ASYNC_AVAILABLE)) {
      run_task <- ExtendedTask$new(function(layers, params) promises::future_promise(sfc_snf_run_with_stability(layers, params), seed = TRUE))
      observeEvent(input$run_btn, {
        validate(need(isTRUE(sc_elig()$ok), sc_elig()$reason))
        layers <- snapshot_inputs()
        p <- sc_params()
        ## run_task$invoke() calls its ExtendedTask function body
        ## SYNCHRONOUSLY before returning (confirmed directly on the
        ## sibling Multi-omics Integration submodule's identical pattern -
        ## tens of seconds in this app's actual reactive-domain context,
        ## not the near-instant dispatch its docs describe in isolation) -
        ## so calling it inline here would block THIS SAME observer, and
        ## nothing (including snapshot_inputs()'s `submitted` flag) can
        ## flush to the browser until a running observer/render finishes.
        ## session$onFlushed(..., once = TRUE) defers the actual invoke()
        ## to the next flush cycle, so the "Running..." spinner reaches
        ## the browser first. `layers`/`p` are snapshotted here, not read
        ## live inside the deferred callback, since that callback runs
        ## outside a normal reactive-consumer context.
        session$onFlushed(function() run_task$invoke(layers, p), once = TRUE)
      })
      observe({
        out <- tryCatch(run_task$result(), error = function(e) e)
        if (inherits(out, "shiny.silent.error")) return()
        state$submitted <- FALSE
        if (inherits(out, "error")) { state$error <- conditionMessage(out); return() }
        if (!isTRUE(out$ok)) { state$error <- out$error; return() }
        state$result <- out$res
        state$stability <- out$stability
      })
      run_running <- reactive(identical(run_task$status(), "running") || isTRUE(state$submitted))
    } else {
      observeEvent(input$run_btn, {
        validate(need(isTRUE(sc_elig()$ok), sc_elig()$reason))
        showNotification("Running SNF Clustering synchronously - the app will be briefly unresponsive.", type = "message", duration = 5)
        layers <- snapshot_inputs()
        out <- sfc_snf_run_with_stability(layers, sc_params())
        state$submitted <- FALSE
        if (!isTRUE(out$ok)) { state$error <- out$error } else { state$result <- out$res; state$stability <- out$stability }
      })
      run_running <- reactive(FALSE)
    }

    ## Mutates the existing "Run SNF Clustering" button in place via shinyjs
    ## (already active app-wide, ui.R::useShinyjs()) rather than recreating
    ## it through renderUI - recreating an actionButton resets its click
    ## counter to 0, which would break the "has this been run at least
    ## once" gating used throughout this module (input$run_btn > 0).
    observe({
      if (isTRUE(run_running())) {
        shinyjs::disable(ns("run_btn"))
        shinyjs::html(ns("run_btn"), as.character(tagList(icon("spinner", class = "fa-spin"), " Running...")))
      } else {
        shinyjs::enable(ns("run_btn"))
        shinyjs::html(ns("run_btn"), as.character(tagList(icon("play"), " Run SNF Clustering")))
      }
    })

    output$run_status_ui <- renderUI({
      if (!isTRUE(input$run_btn > 0)) return(NULL)
      if (isTRUE(run_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running - see the Clusters tab once this finishes."))
      if (!is.null(state$error)) return(mi_stop(state$error))
      if (!is.null(state$result)) return(mi_ok("Finished - see Clusters / Stability / Clinical / Features."))
      NULL
    })

    gate <- function(body_fn) {
      if (!isTRUE(input$run_btn > 0)) return(multi_empty_state("Set parameters (Data / SNF Setup) and click \"Run SNF Clustering\" to see results here."))
      if (isTRUE(run_running())) return(div(class = "empty-note", icon("spinner", class = "fa-spin"), " Running - the rest of the app stays usable while this runs."))
      if (!is.null(state$error)) return(mi_stop(state$error))
      if (is.null(state$result)) return(multi_empty_state())
      body_fn()
    }

    ## =========================================================================
    ## 6. Clusters (spec section 16) + Analysis Summary / downloads (27-28).
    ## =========================================================================
    output$clusters_ui <- renderUI(gate(function() {
      res <- state$result; p <- res$params
      multi_mode <- identical(p$mode, "multi_omics_snf")
      cl_tab <- table(res$clusters)

      tagList(
        box(width = NULL, title = "SNF Summary", status = "primary", solidHeader = FALSE,
            tags$ul(lapply(sfc_summary_lines(res, state$stability), tags$li)),
            div(class = "table-toolbar", downloadButton(ns("dl_bundle"), "Download full results bundle (ZIP)", class = "btn-sm"))),
        box(width = NULL, title = "Actual parameters used", status = "primary", solidHeader = FALSE,
            tags$ul(
              tags$li(sprintf("Mode: %s", if (multi_mode) "Multi-omics SNF" else "Single-Omics Clustering (only one modality selected)")),
              tags$li(sprintf("K: %d (%s)", p$k, p$k_mode)), tags$li(sprintf("Alpha: %.2f (%s)", p$alpha, p$alpha_mode)),
              tags$li(sprintf("T: %s", if (is.na(p$t)) "not applicable" else sprintf("%d (%s)", p$t, p$t_mode))),
              tags$li(sprintf("Clusters: %d (%s)", p$n_clusters, p$cluster_mode)),
              tags$li(sprintf("Clustering technique: %s", names(MI_SNF_CLUSTER_METHODS)[MI_SNF_CLUSTER_METHODS == (p$cluster_method %||% "spectral")])),
              tags$li(sprintf("Standardized: %s", p$standardize))
            )),
        box(width = NULL, title = "Cluster summary", status = "primary", solidHeader = FALSE,
            div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                mi_stat_card(p$n_samples, "Matched patients"), mi_stat_card(p$n_clusters, "Clusters"),
                lapply(names(cl_tab), function(cl) mi_stat_card(sprintf("%d (%.0f%%)", cl_tab[[cl]], 100 * cl_tab[[cl]] / sum(cl_tab)), sprintf("Cluster %s", cl))))),
        box(width = NULL, title = "Patient cluster plot", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() sfc_spectral_embedding_plot(res$W, res$clusters), ns("cl_embed_plot"), height = "380px"),
            p(class = "submodule-desc", "Spectral embedding of the fused patient similarity network.")),
        box(width = NULL, title = "Cluster heatmap", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() sfc_feature_heatmap(state$layers_used, res$clusters), ns("cl_heatmap"), height = "420px"),
            p(class = "submodule-desc", "Top-variance features per modality, used to build the similarity network, ordered by cluster.")),
        box(width = NULL, title = "Fused similarity network", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mi_snf_fused_heatmap(res$W, res$clusters), ns("cl_fused"), height = "380px"),
            if (length(res$Wall) >= 1) tagList(
              selectInput(ns("cl_inspect_block"), "Inspect one modality's own network", choices = names(res$Wall)),
              multi_plot_or_empty(function() { req(input$cl_inspect_block); mi_snf_fused_heatmap(res$Wall[[input$cl_inspect_block]], res$clusters) }, ns("cl_block_net"), height = "340px")
            )),
        box(width = NULL, title = "Modality contribution", status = "primary", solidHeader = FALSE,
            if (!multi_mode) multi_empty_state("Not applicable - only one modality was used (Single-Omics Clustering).") else tagList(
              p(class = "submodule-desc", "Normalized Mutual Information between each modality's network and the fused network (SNFtool::concordanceNetworkNMI)."),
              multi_plot_or_empty(function() sfc_concordance_bar_plot(mi_snf_concordance(res)), ns("cl_conc_plot"), height = "260px"),
              DT::dataTableOutput(ns("cl_conc_table"))
            )),
        box(width = NULL, title = "Cluster quality (candidate cluster counts)", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mi_snf_cluster_estimate_plot(res$cluster_estimate), ns("cl_estimate_plot"), height = "260px")),
        box(width = NULL, title = "Patient-to-cluster assignments", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("cl_assign_table")),
            div(class = "table-toolbar", downloadButton(ns("dl_assign"), "Download cluster assignments (CSV)", class = "btn-sm")))
      )
    }))

    output$cl_conc_table <- DT::renderDataTable({
      req(state$result)
      conc <- mi_snf_concordance(state$result)
      if (is.null(conc)) return(DT::datatable(data.frame(Note = "Not available for this configuration.")))
      DT::datatable(conc, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$cl_assign_table <- DT::renderDataTable({
      req(state$result)
      DT::datatable(data.frame(patient_id = names(state$result$clusters), cluster = as.integer(state$result$clusters)), rownames = FALSE, options = list(pageLength = 15), class = "stripe hover compact")
    })
    output$dl_assign <- downloadHandler(function() "snf_cluster_assignments.csv", function(file) {
      res <- req(state$result)
      utils::write.csv(data.frame(patient_id = names(res$clusters), cluster = as.integer(res$clusters)), file, row.names = FALSE)
    })

    ## multi_plot_or_empty() only decides whether to draw the plotOutput
    ## placeholder or an empty-state message - it does not register the
    ## render logic. Each placeholder used above needs its own renderPlot
    ## binding (mirrors mod_multi_integration.R's own d_error_plot/etc.).
    output$cl_embed_plot <- renderPlot({ res <- req(state$result); sfc_spectral_embedding_plot(res$W, res$clusters) })
    output$cl_heatmap <- renderPlot({ res <- req(state$result); sfc_feature_heatmap(state$layers_used, res$clusters) })
    output$cl_fused <- renderPlot({ res <- req(state$result); mi_snf_fused_heatmap(res$W, res$clusters) })
    output$cl_block_net <- renderPlot({
      res <- req(state$result); req(input$cl_inspect_block)
      mi_snf_fused_heatmap(res$Wall[[input$cl_inspect_block]], res$clusters)
    })
    output$cl_conc_plot <- renderPlot({ res <- req(state$result); sfc_concordance_bar_plot(mi_snf_concordance(res)) })
    output$cl_estimate_plot <- renderPlot({ res <- req(state$result); mi_snf_cluster_estimate_plot(res$cluster_estimate) })

    ## =========================================================================
    ## 7. Stability (spec section 18, REQUIRED) + Sensitivity (section 19).
    ## =========================================================================
    output$stability_ui <- renderUI(gate(function() {
      stab <- state$stability
      tagList(
        box(width = NULL, title = "Cluster stability (default: 20 resamples, 80% subsample)", status = "primary", solidHeader = FALSE,
            if (!isTRUE(stab$ok)) mi_warn(stab$error %||% "Stability could not be computed.") else tagList(
              div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                  mi_stat_card(stab$verdict, "Verdict", if (identical(stab$verdict, "Stable")) ARTHOMIX_STATUS$good else if (identical(stab$verdict, "Moderately stable")) ARTHOMIX_STATUS$warning else ARTHOMIX_STATUS$critical),
                  mi_stat_card(sprintf("%.2f", stab$mean_ari), "Mean ARI"), mi_stat_card(sprintf("%.2f", stab$sd_ari), "SD ARI"),
                  mi_stat_card(stab$n_resamples, "Successful resamples")),
              multi_plot_or_empty(function() sfc_stability_plot(stab), ns("st_plot"), height = "300px"),
              p(class = "submodule-desc", sprintf("Verdict thresholds: mean ARI >= %.2f = Stable, >= %.2f = Moderately stable, else Unstable. (This is a per-sample cluster-membership-agreement stability, distinct from Biomarker Discovery's feature-selection-frequency stability, which uses its own >=%.0f%% “Stable” cutoff on a different, uncorrected-proportion scale.)", SFC_STABILITY_THRESHOLDS$stable, SFC_STABILITY_THRESHOLDS$moderate, MB_STABILITY_THRESHOLDS$stable * 100)),
              div(class = "table-toolbar", downloadButton(ns("dl_stability"), "Download stability metrics (CSV)", class = "btn-sm"))
            )),
        box(width = NULL, title = "Recompute stability with custom settings", status = "primary", solidHeader = FALSE,
            fluidRow(
              column(4, numericInput(ns("st_n_resamples"), "Number of resamples", value = 20, min = 5, max = 100, step = 5)),
              column(4, sliderInput(ns("st_subsample_frac"), "Subsample fraction", min = 0.5, max = 0.95, value = 0.8, step = 0.05)),
              column(4, numericInput(ns("st_seed"), "Random seed", value = 1, min = 1))
            ),
            actionButton(ns("st_recompute_btn"), "Recompute stability", icon = icon("rotate"), class = "btn-primary btn-sm"),
            uiOutput(ns("st_custom_ui"))),
        box(width = NULL, title = "Parameter sensitivity (K / Alpha / T / cluster count)", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Reruns clustering at the low/high end of each parameter and compares to the reference by ARI."),
            actionButton(ns("sens_run_btn"), "Run parameter sensitivity check", icon = icon("play"), class = "btn-primary btn-sm"),
            uiOutput(ns("sens_ui")))
      )
    }))

    st_custom <- eventReactive(input$st_recompute_btn, {
      req(state$result, state$layers_used)
      sfc_stability_run(state$layers_used, state$result$clusters, state$result$params,
                         n_resamples = input$st_n_resamples %||% 20, subsample_frac = input$st_subsample_frac %||% 0.8, seed = input$st_seed %||% 1)
    }, ignoreInit = TRUE)
    output$st_custom_ui <- renderUI({
      stab <- tryCatch(st_custom(), error = function(e) NULL)
      if (is.null(stab)) return(NULL)
      if (!isTRUE(stab$ok)) return(mi_warn(stab$error))
      tagList(
        div(style = "display:flex; gap:10px; flex-wrap:wrap;",
            mi_stat_card(stab$verdict, "Verdict"), mi_stat_card(sprintf("%.2f", stab$mean_ari), "Mean ARI"), mi_stat_card(stab$n_resamples, "Successful resamples")),
        multi_plot_or_empty(function() sfc_stability_plot(stab), ns("st_custom_plot"), height = "280px")
      )
    })
    output$dl_stability <- downloadHandler(function() "snf_stability_metrics.csv", function(file) {
      stab <- req(state$stability)
      utils::write.csv(data.frame(resample = seq_along(stab$ari), ari = stab$ari), file, row.names = FALSE)
    })
    output$st_plot <- renderPlot({ req(isTRUE(state$stability$ok)); sfc_stability_plot(state$stability) })
    output$st_custom_plot <- renderPlot({ stab <- req(st_custom()); req(isTRUE(stab$ok)); sfc_stability_plot(stab) })

    sens <- eventReactive(input$sens_run_btn, {
      req(state$result, state$layers_used)
      sfc_sensitivity_run(state$layers_used, state$result$params, state$result$clusters, seed = 1)
    }, ignoreInit = TRUE)
    output$sens_ui <- renderUI({
      s <- tryCatch(sens(), error = function(e) NULL)
      if (is.null(s)) return(NULL)
      tagList(
        multi_plot_or_empty(function() sfc_sensitivity_plot(s), ns("sens_plot"), height = "280px"),
        DT::dataTableOutput(ns("sens_table")),
        p(class = "submodule-desc", "Low sensitivity = stable assignments across the range; high = they change substantially.")
      )
    })
    output$sens_table <- DT::renderDataTable({
      DT::datatable(req(sens())$summary, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$sens_plot <- renderPlot({ sfc_sensitivity_plot(req(sens())) })

    ## =========================================================================
    ## 8. Clinical (spec sections 20-23) - fully optional, only ever offers
    ## variables actually detected in this dataset's own metadata; never
    ## fabricates a field. Always post-hoc, never used to pick K/alpha/T/n_clusters.
    ## =========================================================================
    sc_clinical <- reactive({
      req(state$result)
      sfc_detect_clinical(state$sample_meta)
    })

    output$clinical_ui <- renderUI(gate(function() {
      det <- sc_clinical()
      if (length(det$categorical) == 0 && length(det$continuous) == 0 && is.null(det$survival)) {
        return(box(width = NULL, title = "Clinical association", status = "primary", solidHeader = FALSE, multi_empty_state("Clinical outcome unavailable for this dataset.")))
      }
      tagList(
        p(class = "submodule-desc", tags$strong("This is post-hoc evaluation, not model validation"), " - clusters were formed without clinical variables; this checks alignment afterward only."),
        if (length(det$categorical) > 0) box(width = NULL, title = "Categorical variables (Fisher's exact test)", status = "primary", solidHeader = FALSE,
            checkboxGroupInput(ns("clin_cat_vars"), "Select variable(s) to test", choices = det$categorical, inline = TRUE),
            uiOutput(ns("clin_cat_ui"))),
        if (length(det$continuous) > 0) box(width = NULL, title = "Continuous variables (Kruskal-Wallis test)", status = "primary", solidHeader = FALSE,
            checkboxGroupInput(ns("clin_cont_vars"), "Select variable(s) to test", choices = det$continuous, inline = TRUE),
            uiOutput(ns("clin_cont_ui"))),
        box(width = NULL, title = "Survival", status = "primary", solidHeader = FALSE,
            if (is.null(det$survival) || !requireNamespace("survival", quietly = TRUE)) multi_empty_state("Survival data were not detected for this dataset.")
            else uiOutput(ns("clin_surv_ui")))
      )
    }))

    ## One combined plot (faceted by variable) + one results table per kind,
    ## rather than a plot per selected checkbox - avoids dynamically-named
    ## renderPlot bindings (see the comment above cl_embed_plot/etc.).
    clin_cat_results <- reactive({
      req(length(input$clin_cat_vars) > 0)
      sfc_clinical_run(state$result$clusters, state$sample_meta, input$clin_cat_vars, kind = "categorical")
    })
    output$clin_cat_ui <- renderUI({
      if (length(input$clin_cat_vars %||% character(0)) == 0) return(multi_empty_state("Select at least one variable above."))
      res <- clin_cat_results()
      tagList(
        if (length(res) > 1) p(class = "submodule-desc", "Multiple variables selected - p-values BH-corrected across the selection shown."),
        multi_plot_or_empty(function() sfc_categorical_multi_plot(state$result$clusters, state$sample_meta, input$clin_cat_vars), ns("clin_cat_plot"), height = "300px"),
        DT::dataTableOutput(ns("clin_cat_table"))
      )
    })
    output$clin_cat_plot <- renderPlot({ req(length(input$clin_cat_vars) > 0); sfc_categorical_multi_plot(state$result$clusters, state$sample_meta, input$clin_cat_vars) })
    output$clin_cat_table <- DT::renderDataTable({
      res <- req(clin_cat_results())
      df <- do.call(rbind, lapply(names(res), function(v) {
        r <- res[[v]]
        if (!isTRUE(r$ok)) return(data.frame(variable = v, test = NA, p_value = NA, p_fdr = NA, effect = NA, n = NA, note = r$error))
        data.frame(variable = v, test = r$test, p_value = signif(r$p_value, 3), p_fdr = signif(r$p_fdr, 3), effect = signif(r$effect, 3), n = r$n, note = "")
      }))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    clin_cont_results <- reactive({
      req(length(input$clin_cont_vars) > 0)
      sfc_clinical_run(state$result$clusters, state$sample_meta, input$clin_cont_vars, kind = "continuous")
    })
    output$clin_cont_ui <- renderUI({
      if (length(input$clin_cont_vars %||% character(0)) == 0) return(multi_empty_state("Select at least one variable above."))
      res <- clin_cont_results()
      tagList(
        if (length(res) > 1) p(class = "submodule-desc", "Multiple variables selected - p-values BH-corrected across the selection shown."),
        multi_plot_or_empty(function() sfc_continuous_multi_plot(state$result$clusters, state$sample_meta, input$clin_cont_vars), ns("clin_cont_plot"), height = "300px"),
        DT::dataTableOutput(ns("clin_cont_table"))
      )
    })
    output$clin_cont_plot <- renderPlot({ req(length(input$clin_cont_vars) > 0); sfc_continuous_multi_plot(state$result$clusters, state$sample_meta, input$clin_cont_vars) })
    output$clin_cont_table <- DT::renderDataTable({
      res <- req(clin_cont_results())
      df <- do.call(rbind, lapply(names(res), function(v) {
        r <- res[[v]]
        if (!isTRUE(r$ok)) return(data.frame(variable = v, test = NA, p_value = NA, p_fdr = NA, n = NA, note = r$error))
        data.frame(variable = v, test = r$test, p_value = signif(r$p_value, 3), p_fdr = signif(r$p_fdr, 3), n = r$n, note = "")
      }))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    clin_surv <- reactive({
      det <- sc_clinical(); req(det$survival)
      time <- stats::setNames(state$sample_meta[[det$survival$time_col]], rownames(state$sample_meta))
      event <- stats::setNames(state$sample_meta[[det$survival$event_col]], rownames(state$sample_meta))
      sfc_test_survival(state$result$clusters, time, event)
    })
    output$clin_surv_ui <- renderUI({
      surv <- tryCatch(clin_surv(), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
      if (!isTRUE(surv$ok)) return(mi_warn(surv$error))
      tagList(
        multi_plot_or_empty(function() sfc_km_plot(surv), ns("clin_km_plot"), height = "360px"),
        div(style = "display:flex; gap:10px; flex-wrap:wrap;",
            mi_stat_card(sprintf("%.3g", surv$logrank_p), "Log-rank p"),
            if (!is.null(surv$hr)) mi_stat_card(sprintf("%.2f (%.2f-%.2f)", surv$hr$hr, surv$hr$lo, surv$hr$hi), "Hazard ratio (95% CI)"),
            mi_stat_card(surv$n, "Patients with survival data")),
        h5("Number at risk"),
        DT::dataTableOutput(ns("clin_km_risk_table")),
        p(class = "submodule-desc", "Log-rank test compares survival across clusters; hazard ratio shown only for 2 clusters (Cox).")
      )
    })
    output$clin_km_risk_table <- DT::renderDataTable({
      surv <- req(clin_surv())
      DT::datatable(sfc_km_risk_table(surv), rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })
    output$clin_km_plot <- renderPlot({ surv <- req(clin_surv()); req(isTRUE(surv$ok)); sfc_km_plot(surv) })

    ## =========================================================================
    ## 9. Features (spec section 24) - association with the already-computed
    ## clusters, distinct from "features used to construct the network"
    ## (every feature in a selected block was used for that).
    ## =========================================================================
    output$features_ui <- renderUI(gate(function() {
      res <- state$result
      tagList(
        p(class = "submodule-desc", "Features ranked here correlate with the unsupervised clusters (Kruskal-Wallis test), not necessarily the same features used to build the network."),
        fluidRow(column(6, selectInput(ns("feat_block"), "Modality", choices = names(state$layers_used))),
                 column(6, numericInput(ns("feat_top_n"), "Features to show", value = 25, min = 5, max = 200, step = 5))),
        uiOutput(ns("feat_result_ui"))
      )
    }))

    feat_rank <- reactive({
      req(input$feat_block, state$layers_used, state$result)
      sfc_feature_ranking(state$layers_used[[input$feat_block]], state$result$clusters, input$feat_block, top_n = input$feat_top_n %||% 25)
    })
    output$feat_result_ui <- renderUI({
      r <- tryCatch(feat_rank(), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
      if (!isTRUE(r$ok)) return(mi_warn(r$error))
      tagList(
        multi_plot_or_empty(function() sfc_feature_rank_plot(r$table), ns("feat_plot"), height = "420px"),
        DT::dataTableOutput(ns("feat_table")),
        div(class = "table-toolbar", downloadButton(ns("dl_feat"), "Download feature ranking (CSV)", class = "btn-sm")),
        p(class = "submodule-desc", sprintf("%d features tested in this modality; BH-FDR corrected.", r$n_features_tested))
      )
    })
    output$feat_table <- DT::renderDataTable({
      DT::datatable(req(feat_rank())$table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$feat_plot <- renderPlot({ r <- req(feat_rank()); req(isTRUE(r$ok)); sfc_feature_rank_plot(r$table) })
    output$dl_feat <- downloadHandler(function() sprintf("snf_feature_ranking_%s.csv", make.names(input$feat_block %||% "block")), function(file) {
      utils::write.csv(req(feat_rank())$table, file, row.names = FALSE)
    })

    ## =========================================================================
    ## 10. Full results bundle (spec section 28) - everything computed so
    ## far this session, whatever subset that is; nothing fabricated for
    ## tabs the user never opened.
    ## =========================================================================
    output$dl_bundle <- downloadHandler(
      filename = function() paste0("snf_clustering_bundle_", Sys.Date(), ".zip"),
      content = function(file) {
        res <- req(state$result)
        tmp <- tempfile(); dir.create(tmp)
        utils::write.csv(data.frame(patient_id = names(res$clusters), cluster = as.integer(res$clusters)), file.path(tmp, "cluster_assignments.csv"), row.names = FALSE)
        if (!is.null(state$stability) && isTRUE(state$stability$ok)) utils::write.csv(data.frame(resample = seq_along(state$stability$ari), ari = state$stability$ari), file.path(tmp, "stability_metrics.csv"), row.names = FALSE)
        conc <- tryCatch(mi_snf_concordance(res), error = function(e) NULL)
        if (!is.null(conc)) utils::write.csv(conc, file.path(tmp, "modality_contribution.csv"), row.names = FALSE)
        feat <- tryCatch(feat_rank(), error = function(e) NULL)
        if (!is.null(feat) && isTRUE(feat$ok)) utils::write.csv(feat$table, file.path(tmp, "feature_ranking.csv"), row.names = FALSE)
        writeLines(c(
          "# SNF Clustering - reproducibility / analysis summary",
          sprintf("Dataset: %s", state$dataset_label %||% "-"),
          sfc_summary_lines(res, state$stability),
          sprintf("SNFtool version: %s", tryCatch(as.character(utils::packageVersion("SNFtool")), error = function(e) "unknown")),
          "Clusters are unsupervised molecular groups; any clinical association reported alongside is post-hoc and not used to choose the cluster count."
        ), file.path(tmp, "analysis_summary.md"))
        old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
        setwd(tmp)
        utils::zip(file, files = list.files(tmp))
      }
    )

    ## =========================================================================
    ## Publish - kept in the same shape multi_qc_scorecard()/
    ## multi_analysis_summary_table() (multiomics_helpers.R) already read
    ## (`r$stratification$drug` as a display label) so those shared,
    ## untouched cross-submodule summaries keep working unmodified.
    ## =========================================================================
    observe({
      if (is.null(multi_results)) return()
      d <- tryCatch(sc_dataset(), error = function(e) NULL)
      if (is.null(d) || !isTRUE(d$ok)) return()
      multi_results$stratification <- list(
        drug = d$label,
        params = if (!is.null(state$result)) state$result$params else NULL,
        stability = if (!is.null(state$stability) && isTRUE(state$stability$ok)) state$stability$verdict else NULL,
        clusters = if (!is.null(state$result)) state$result$clusters else NULL
      )
    })
  })
}
