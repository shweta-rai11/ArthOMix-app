## server.R
## ArthOMix Explorer

function(input, output, session) {

  dataset <- local({
    d <- load_default_dataset()
    d$source_type <- "preloaded"
    d$is_bundled_reference <- TRUE
    d$geo_ids <- MERGED_DEFAULT_GEO_IDS
    d$declared_data_type <- NA_character_
    d$load_id <- 0L
    do.call(reactiveValues, d)
  })

  results <- reactiveValues()

  ## Session-wide analysis-record + diagnostics stores (R/provenance.R), shown by the header "Analysis records" button.
  session$userData$arthomix_provenance <- reactiveVal(list())
  session$userData$arthomix_diagnostics <- reactiveVal(list())

  output$analysis_records_badge <- renderUI({
    n <- length(arthomix_provenance_records(session))
    if (n == 0) return(NULL)
    span(class = "app-header-badge", n)
  })

  observeEvent(input$analysis_records_btn, {
    showModal(modalDialog(
      title = tagList(icon("clipboard-list"), " Analysis records - this session"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      p(class = "submodule-desc",
        "One row per analysis run in this session, in the order they ran. Each record holds the exact parameters, seed, R and package versions, and a checksum of the input data. Download the full JSON to attach to a report or a manuscript supplement. Records are kept when the dataset changes; results are not."),
      p(class = "submodule-desc", icon("circle-info"),
        " Coverage: every analysis run button in the Transcriptomics module writes a record here (Overview normalisation, Preprocessing, Differential Expression, WGCNA, Candidate Genes, MR, Colocalisation, Feature Selection, Diagnostic Model incl. External Validation, Sex Interaction, Cross-Tissue Replication, Cross-Ancestry MR Replication, Functional Enrichment, Immune Deconvolution, Biomarker Card); Cross-Omics (Expression/Methylation Integration, Biomarker Convergence, MR Evidence); and Multi-Omics (Cohort Harmonization, DIABLO/SNF Integration, SNF Clustering, Biomarker Discovery, Gene-CpG Mapping, Pathways). Methylomics is not yet wired to this shared log - its DMP tab writes its own downloadable manifest, as does Cross-Omics' Expression/Methylation Integration tab (in addition to writing here)."),
      div(class = "table-toolbar",
          downloadButton("download_analysis_records", "Download all records (.json)", class = "btn-sm btn-default"),
          actionButton("clear_analysis_records_btn", "Clear records", icon = icon("trash"), class = "btn-sm btn-default")),
      DT::dataTableOutput("analysis_records_table"),
      tags$hr(),
      h4(icon("stethoscope"), " Diagnostics: absorbed errors and warnings"),
      p(class = "submodule-desc",
        "When an optional step fails, the app shows \"Not available\" rather than crashing the tab. The reason is kept here so you can see why. Identical messages are collapsed with a count, newest first."),
      div(class = "table-toolbar",
          actionButton("clear_diagnostics_btn", "Clear diagnostics", icon = icon("trash"), class = "btn-sm btn-default")),
      DT::dataTableOutput("analysis_diagnostics_table")
    ))
  }, ignoreInit = TRUE)

  output$analysis_records_table <- DT::renderDataTable({
    tbl <- arthomix_provenance_summary_table(arthomix_provenance_records(session))
    DT::datatable(tbl, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE, order = list(list(0, "desc"))),
                  class = "stripe hover compact")
  })
  output$analysis_diagnostics_table <- DT::renderDataTable({
    tbl <- arthomix_diag_summary_table(arthomix_diag_entries(session))
    DT::datatable(tbl, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  })
  output$download_analysis_records <- downloadHandler(
    filename = function() sprintf("arthomix_analysis_records_%s.json", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      recs <- arthomix_provenance_json_safe(arthomix_provenance_records(session))
      writeLines(jsonlite::toJSON(recs, pretty = TRUE, auto_unbox = TRUE, force = TRUE, null = "null"), file)
    }
  )
  observeEvent(input$clear_analysis_records_btn, arthomix_provenance_clear(session), ignoreInit = TRUE)
  observeEvent(input$clear_diagnostics_btn, arthomix_diag_clear(session), ignoreInit = TRUE)

  observeEvent(dataset$source, {
    for (nm in names(results)) results[[nm]] <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$home_browse_modules, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "modules")
  }, ignoreInit = TRUE)
  observeEvent(input$home_cta_browse_modules, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "modules")
  }, ignoreInit = TRUE)
  lapply(MODULE_REGISTRY, function(m) {
    if (identical(m$status, "available") && !identical(m$id, "arthochat")) {
      observeEvent(input[[paste0("home_card_open_", m$id)]], {
        updateTabsetPanel(session, "sidebar_tabs", selected = m$tab)
      }, ignoreInit = TRUE)
    }
  })

  lapply(MODULE_REGISTRY, function(m) {
    if (identical(m$status, "available") && !identical(m$id, "arthochat")) {
      observeEvent(input[[paste0("open_", m$id)]], {
        updateTabsetPanel(session, "sidebar_tabs", selected = m$tab)
      }, ignoreInit = TRUE)
    }
  })

  mod_dataset_server("tx_dataset", dataset)

  methyl_dataset <- reactiveValues(
    beta = NULL, input_scale = NULL, array_type = NULL, sample_sheet = NULL,
    rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL, source = NULL,
    preloaded = FALSE,
    source_type = NULL
  )
  methyl_results <- reactiveValues()
  mod_methyl_dataset_server("mx_dataset", methyl_dataset)
  lapply(MX_MODULES, function(m) m$server(paste0("mx_", m$config$id), methyl_dataset, methyl_results))
  observeEvent(methyl_dataset$source, {
    for (nm in names(methyl_results)) methyl_results[[nm]] <- NULL
  }, ignoreInit = TRUE)

  cross_dataset <- reactiveValues(
    user_expr_df = NULL, user_expr_source = NULL, user_expr_wide = NULL, user_expr_mapping = NULL, user_expr_sample_cols = character(0),
    user_meth_df = NULL, user_meth_source = NULL, user_meth_wide = NULL, user_meth_mapping = NULL, user_meth_sample_cols = character(0)
  )
  cross_results <- reactiveValues()
  mod_cross_dataset_server("cx_dataset", cross_dataset, results, methyl_results)
  lapply(CX_MODULES, function(m) {
    if (identical(m$config$id, "integration")) {
      m$server(paste0("cx_", m$config$id), cross_dataset, cross_results, dataset, results, methyl_dataset, methyl_results)
    } else if (identical(m$config$id, "mrstage")) {
      m$server(paste0("cx_", m$config$id), cross_dataset, cross_results, app_session = session)
    } else {
      m$server(paste0("cx_", m$config$id), cross_dataset, cross_results)
    }
  })

  multi_dataset <- reactiveValues(
    table_label = NULL, df = NULL, source = NULL,
    layers = list(), layer_meta = list(), sample_meta = NULL,
    overlap = NULL, active = FALSE, loaded_at = NULL
  )
  multi_results <- reactiveValues()
  mod_multi_dataset_server("mo_dataset", multi_dataset, multi_results)
  lapply(MULTI_MODULES, function(m) m$server(paste0("mo_", m$config$id), multi_dataset, multi_results))
  observeEvent(input[["mo_dataset-goto_crossomics"]], {
    updateTabsetPanel(session, "sidebar_tabs", selected = "crossomics")
  })
  observeEvent(multi_dataset$source, {
    for (nm in names(multi_results)) multi_results[[nm]] <- NULL
  }, ignoreInit = TRUE)

  agent_run_hooks <- new.env(parent = emptyenv())
  agent_run_hooks$transcriptomics <- list()
  title_to_module_id <- function(modules_list, title) {
    hit <- Find(function(m) identical(m$config$title, title), modules_list)
    if (is.null(hit)) NULL else hit$config$id
  }

  current_module_context <- reactive({
    top <- input$sidebar_tabs %||% "home"
    switch(top,
      transcriptomics = list(
        module = "transcriptomics",
        submodule_id = title_to_module_id(TX_MODULES, input$tx_menu),
        view_label = if (is.null(input$tx_menu)) "Transcriptomics" else sprintf("Transcriptomics / %s", input$tx_menu)
      ),
      methylomics = list(
        module = "methylomics",
        submodule_id = title_to_module_id(MX_MODULES, input$mx_menu),
        view_label = if (is.null(input$mx_menu)) "Methylomics" else sprintf("Methylomics / %s", input$mx_menu)
      ),
      crossomics = list(
        module = "crossomics",
        submodule_id = title_to_module_id(CX_MODULES, input$cx_menu),
        view_label = if (is.null(input$cx_menu)) "Cross-Omics" else sprintf("Cross-Omics / %s", input$cx_menu)
      ),
      multiomics = list(
        module = "multiomics",
        submodule_id = title_to_module_id(MULTI_MODULES, input$mo_menu),
        view_label = if (is.null(input$mo_menu)) "Multi-Omics" else sprintf("Multi-Omics / %s", input$mo_menu)
      ),
      list(module = "app", submodule_id = NULL, view_label = "ArthOMix (no specific module open)")
    )
  })

  output$cx_arthochat_hint <- renderUI({
    hint <- if (identical(input$cx_menu, "MR Evidence")) {
      "Need help selecting MR data - which sex's evidence, or preloaded vs. your own upload? Ask ArthOChat."
    } else {
      "Questions about how the panels converge, or which sex/data source to select? Ask ArthOChat."
    }
    arthochat_shortcut_ui(hint, compact = TRUE)
  })

  if (ARTHOMIX_CHAT_ENABLED) {
    mod_arthochat_server(
      "arthochat", dataset, results,
      methyl_dataset, methyl_results,
      cross_dataset, cross_results,
      multi_dataset, multi_results,
      current_context = current_module_context,
      run_hooks = agent_run_hooks
    )
  }

  tx_server_hooks <- setNames(
    lapply(TX_MODULES, function(m) m$server(paste0("tx_", m$config$id), dataset, results)),
    vapply(TX_MODULES, function(m) m$config$id, character(1))
  )
  agent_run_hooks$transcriptomics$dge <- tx_server_hooks$dge$run

  added <- reactiveValues(ids = character(0))

  ## Same steps as pressing a sub-module card's "Add" button; also used by the guided workflow.
  ## Inserted by pipeline rank (step, then catalogue order) rather than in click order, so the tabs
  ## behind the hidden strip always read Dataset -> step 2 -> ... -> step 9 -> Sub-modules.
  tx_insert_submodule <- function(m) {
    hid <- m$config$id
    r <- tx_rank(hid)
    later <- Filter(function(x) tx_rank(x) > r, added$ids)
    target <- if (length(later)) {
      TX_MODULES_BY_ID[[later[[which.min(vapply(later, tx_rank, integer(1)))]]]]$config$title
    } else "Sub-modules"
    insertTab(
      session = session, inputId = "tx_menu",
      tabPanel(m$config$title, br(), m$ui(paste0("tx_", hid))),
      target = target, position = "before", select = TRUE
    )
    added$ids <- union(added$ids, hid)
    shinyjs::addClass(id = paste0("smcard_", hid), class = "sm-card-active")
    shinyjs::html(id = paste0("smstate_", hid), html = "Added")
  }

  lapply(TX_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% added$ids) {
        removeTab(session = session, inputId = "tx_menu", target = m$config$title)
        added$ids <- setdiff(added$ids, hid)
        shinyjs::removeClass(id = paste0("smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("smstate_", hid), html = "Add")
      } else {
        tx_insert_submodule(m)
      }
    }, ignoreInit = TRUE)
  })

  output$sm_active_count <- renderUI({
    n <- length(added$ids)
    span(
      class = "sm-count-badge",
      sprintf("%d of %d sub-modules added", n, length(TX_MODULES))
    )
  })

  observeEvent(input$sm_search, {
    q <- tolower(trimws(input$sm_search %||% ""))
    lapply(TX_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  mx_added <- reactiveValues(ids = character(0))

  ## Same steps as pressing a sub-module card's "Add" button; also used by the guided workflow.
  mx_insert_submodule <- function(m) {
    hid <- m$config$id
    insertTab(
      session = session, inputId = "mx_menu",
      tabPanel(m$config$title, br(), m$ui(paste0("mx_", hid))),
      target = "Sub-modules", position = "before", select = TRUE
    )
    mx_added$ids <- union(mx_added$ids, hid)
    shinyjs::addClass(id = paste0("mx_smcard_", hid), class = "sm-card-active")
    shinyjs::html(id = paste0("mx_smstate_", hid), html = "Added")
  }

  lapply(MX_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("mx_sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% mx_added$ids) {
        removeTab(session = session, inputId = "mx_menu", target = m$config$title)
        mx_added$ids <- setdiff(mx_added$ids, hid)
        shinyjs::removeClass(id = paste0("mx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mx_smstate_", hid), html = "Add")
      } else {
        mx_insert_submodule(m)
      }
    }, ignoreInit = TRUE)
  })

  output$mx_sm_active_count <- renderUI({
    n <- length(mx_added$ids)
    span(class = "sm-count-badge", sprintf("%d of %d sub-modules added", n, length(MX_MODULES)))
  })

  observeEvent(input$mx_sm_search, {
    q <- tolower(trimws(input$mx_sm_search %||% ""))
    lapply(MX_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("mx_smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  cx_added <- reactiveValues(ids = character(0))

  lapply(CX_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("cx_sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% cx_added$ids) {
        removeTab(session = session, inputId = "cx_menu", target = m$config$title)
        cx_added$ids <- setdiff(cx_added$ids, hid)
        shinyjs::removeClass(id = paste0("cx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("cx_smstate_", hid), html = "Add")
      } else {
        insertTab(
          session = session, inputId = "cx_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("cx_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        cx_added$ids <- union(cx_added$ids, hid)
        shinyjs::addClass(id = paste0("cx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("cx_smstate_", hid), html = "Added")
      }
    }, ignoreInit = TRUE)
  })

  output$cx_sm_active_count <- renderUI({
    n <- length(cx_added$ids)
    span(class = "sm-count-badge", sprintf("%d of %d sub-modules added", n, length(CX_MODULES)))
  })

  observeEvent(input$cx_sm_search, {
    q <- tolower(trimws(input$cx_sm_search %||% ""))
    lapply(CX_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("cx_smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  mo_added <- reactiveValues(ids = character(0))

  lapply(MULTI_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("mo_sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% mo_added$ids) {
        removeTab(session = session, inputId = "mo_menu", target = m$config$title)
        mo_added$ids <- setdiff(mo_added$ids, hid)
        shinyjs::removeClass(id = paste0("mo_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mo_smstate_", hid), html = "Add")
      } else {
        insertTab(
          session = session, inputId = "mo_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("mo_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        mo_added$ids <- union(mo_added$ids, hid)
        shinyjs::addClass(id = paste0("mo_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mo_smstate_", hid), html = "Added")
      }
    }, ignoreInit = TRUE)
  })

  output$mo_sm_active_count <- renderUI({
    n <- length(mo_added$ids)
    span(class = "sm-count-badge", sprintf("%d of %d sub-modules added", n, length(MULTI_MODULES)))
  })

  observeEvent(input$mo_sm_search, {
    q <- tolower(trimws(input$mo_sm_search %||% ""))
    lapply(MULTI_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("mo_smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  jump_to_submodule <- function(mod_id, inner_tab = NULL, sm_filter = NULL, inner_tabset = "tabs") {
    cfg <- TX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% added$ids) {
      updateTabsetPanel(session, "tx_menu", selected = cfg$title)
      if (!is.null(inner_tab)) {
        updateTabsetPanel(session, paste0("tx_", mod_id, "-", inner_tabset), selected = inner_tab)
      }
    } else {
      updateTabsetPanel(session, "tx_menu", selected = "Sub-modules")
      updateTextInput(session, "sm_search", value = sm_filter %||% cfg$title)
    }
  }

  ## ---------------------------------------------------------------------------------------------
  ## Transcriptomics workflow shell (R/workflow_guide.R). Navigation only: it reads the modules'
  ## existing results and opens their existing tabs, and computes no statistics of its own.
  ## ---------------------------------------------------------------------------------------------
  wf_tx <- reactiveValues(
    contrast = NULL, sex_mode = NULL, selected = character(0),
    stamps = list(), fingerprints = list(), seq = 0L, active_key = "__dataset__",
    configured = FALSE
  )

  tx_dataset_state <- reactive({
    meta <- dataset$meta
    list(
      ok = !is.null(dataset$expr) && NCOL(dataset$expr) > 0 && is.data.frame(meta) && nrow(meta) > 0,
      sex_ok = length(wf_sex_levels(meta, wf_sex_col(meta))) >= 2,
      load_id = dataset$load_id %||% 0L
    )
  })

  ## The one comparison every module inherits: the guide's, if it has set one; otherwise whatever the
  ## Differential Expression page itself is showing; otherwise the default that page would pick from
  ## the loaded metadata. Read-only everywhere else, which is what the banner says.
  tx_effective_contrast <- reactive({
    if (!is.null(wf_tx$contrast)) return(c(wf_tx$contrast, list(source = "guide")))
    col <- input[["tx_dge-contrast_col"]]
    ref <- input[["tx_dge-ref_group"]]
    comp <- input[["tx_dge-comp_group"]]
    if (!is.null(col) && nzchar(ref %||% "") && nzchar(comp %||% "") && !identical(ref, comp)) {
      return(list(col = col, ref = ref, comp = comp, source = "page"))
    }
    meta <- dataset$meta
    gc <- wf_group_col(meta)
    if (is.null(gc)) return(NULL)
    dc <- wf_default_contrast(meta[[gc]])
    if (is.null(dc)) return(NULL)
    list(col = gc, ref = dc$ref, comp = dc$comp, source = "default")
  })

  tx_wf_settings <- function() {
    ct <- tx_effective_contrast()
    list(contrast = if (is.null(ct)) NULL else ct[c("col", "ref", "comp")],
         sex_mode = wf_tx$sex_mode,
         stamps = wf_tx$stamps, fingerprints = wf_tx$fingerprints)
  }

  ## Every module's status in one place, so the sidebar, the cards and the breadcrumb agree.
  tx_statuses <- reactive({
    ds <- tx_dataset_state()
    set <- tx_wf_settings()
    titles <- wf_step_titles()
    stats::setNames(lapply(names(TX_STEP_MAP), function(k) {
      wf_module_status(k, ds, results, set, titles = titles)
    }), names(TX_STEP_MAP))
  })

  ## Records a fingerprint the moment a module stores a result, so a later change to the dataset, the
  ## contrast, the sex design or an upstream run shows up as "Stale - re-run" here and downstream.
  ## The read of results is inside observe(), never in a module body at setup time.
  lapply(names(TX_STEP_MAP), function(k) {
    entry <- TX_STEP_MAP[[k]]
    observe({
      has <- wf_has_result(entry, results)
      isolate({
        if (!has) {
          if (!is.null(wf_tx$stamps[[k]])) {
            st <- wf_tx$stamps; st[[k]] <- NULL; wf_tx$stamps <- st
            fp <- wf_tx$fingerprints; fp[[k]] <- NULL; wf_tx$fingerprints <- fp
          }
          return()
        }
        wf_tx$seq <- (wf_tx$seq %||% 0L) + 1L
        st <- wf_tx$stamps; st[[k]] <- wf_tx$seq; wf_tx$stamps <- st
        fp <- wf_tx$fingerprints
        set <- isolate(tx_wf_settings())
        set$stamps <- st
        fp[[k]] <- wf_fingerprint(k, isolate(tx_dataset_state()), set)
        wf_tx$fingerprints <- fp
      })
    })
  })

  ## Pipeline rank: position in TX_STEP_MAP_ORDER, which runs step by step and, inside a step, in the
  ## catalogue's own order. Used for tab insert position and for sidebar ordering.
  tx_primary_key <- function(hid) {
    hit <- Find(function(e) identical(e$id, hid) && !e$alias, TX_STEP_MAP)
    if (is.null(hit)) NULL else hit$key
  }
  tx_rank <- function(hid) {
    k <- tx_primary_key(hid)
    if (is.null(k)) length(TX_STEP_MAP_ORDER) + 1L else match(k, TX_STEP_MAP_ORDER)
  }

  ## Opens a module's inner tab once that tab exists; a module inserted in this same flush has not
  ## rendered its own tabset yet, so a plain updateTabsetPanel would be dropped.
  tx_open_inner_tab <- function(mod_id, tabset, tab) {
    sel <- sprintf("#tx_%s-%s li a[data-value='%s']", mod_id, tabset, tab)
    shinyjs::runjs(sprintf(
      "(function retry(n){ var a = $('%s'); if (a.length) { a.tab('show'); } else if (n > 0) { setTimeout(function(){ retry(n - 1); }, 100); } })(40);",
      gsub("'", "\\\\'", sel)
    ))
  }

  ## Single click route for the step sidebar. Event-priority, so re-rendering the list with new
  ## statuses never fires it (per-row action buttons would, every time their counters reset).
  observeEvent(input$tx_nav_go, {
    key <- as.character(input$tx_nav_go)
    wf_tx$active_key <- key
    if (identical(key, "__dataset__")) {
      updateTabsetPanel(session, "tx_menu", selected = "Dataset")
      return()
    }
    if (identical(key, "__catalogue__")) {
      updateTabsetPanel(session, "tx_menu", selected = "Sub-modules")
      return()
    }
    if (identical(key, "__design__")) {
      updateTabsetPanel(session, "tx_menu", selected = "Custom Analysis Design")
      wf_guide$bar_hidden <- FALSE
      return()
    }
    e <- TX_STEP_MAP[[key]]
    req(e)
    m <- TX_MODULES_BY_ID[[e$id]]
    if (e$id %in% added$ids) {
      updateTabsetPanel(session, "tx_menu", selected = m$config$title)
    } else {
      tx_insert_submodule(m)
    }
    if (!is.null(e$tabset)) tx_open_inner_tab(e$id, e$tabset, e$tab)
  }, ignoreInit = TRUE)

  observeEvent(input$tx_crumb_home, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "home")
  }, ignoreInit = TRUE)

  ## A tab opened from anywhere else (search, ArthOChat, a catalogue card, the guide) moves the
  ## highlight too - unless the row already showing is one of that same tab's rows, which is how the
  ## step 8 External Validation row stays highlighted instead of snapping back to step 7.
  observeEvent(input$tx_menu, {
    sel <- input$tx_menu
    hid <- title_to_module_id(TX_MODULES, sel)
    cur <- TX_STEP_MAP[[wf_tx$active_key %||% ""]]
    if (is.null(hid)) {
      wf_tx$active_key <- switch(sel,
                                 "Sub-modules" = "__catalogue__",
                                 "Custom Analysis Design" = "__design__",
                                 "__dataset__")
      return()
    }
    if (!is.null(cur) && identical(cur$id, hid)) return()
    wf_tx$active_key <- tx_primary_key(hid) %||% wf_tx$active_key
  })

  tx_step_keys <- function(n) {
    names(TX_STEP_MAP)[vapply(TX_STEP_MAP, `[[`, integer(1), "step") == n]
  }

  output$tx_sidebar_dynamic_nav <- renderUI({
    sts <- tx_statuses()
    titles <- wf_step_titles()
    ds <- tx_dataset_state()
    active <- wf_tx$active_key %||% "__dataset__"
    shown_ids <- added$ids

    ## The entry point, always first; the tracker under it only once there is progress to track.
    design_btn <- tags$ul(class = "sidebar-nav sidebar-nav-design", tags$li(tags$a(
      href = "#", class = paste("sidebar-nav-item sidebar-design-item",
                                if (identical(active, "__design__")) "active"),
      title = "Customise the analysis - pick the dataset, the comparison, the sex design and which analyses to run",
      onclick = "Shiny.setInputValue('tx_nav_go', '__design__', {priority: 'event'}); return false;",
      icon("compass-drafting"), tags$span(class = "sidebar-step-label", "Custom Analysis Design")
    )))
    add_btn <- tags$ul(class = "sidebar-nav sidebar-nav-foot", tags$li(tags$a(
      href = "#", class = "sidebar-nav-item sidebar-add-item",
      onclick = "Shiny.setInputValue('tx_nav_go', '__catalogue__', {priority: 'event'}); return false;",
      icon("plus"), "Add analyses"
    )))
    open_row <- function(k) {
      e <- TX_STEP_MAP[[k]]
      omics_sidebar_step_row(k, titles[[k]], TX_MODULES_BY_ID[[e$id]]$config$icon,
                             sts[[k]]$status, sts[[k]]$reason, active = identical(active, k))
    }
    ## Loading a dataset is not a run either: the row says what is loaded, but stays neutral.
    dataset_row <- omics_sidebar_step_row(
      "__dataset__", "Dataset", "database", "ready",
      dataset$source %||% "No dataset loaded yet.",
      active = identical(active, "__dataset__"))

    ## Before an analysis is set up there is no progress to track, so nine empty step headers would
    ## be noise. Until then this is a plain list of what is open; it becomes the step tracker the
    ## moment the guide starts.
    if (!isTRUE(wf_tx$configured)) {
      return(tagList(
        design_btn,
        tags$ul(
          class = "sidebar-nav sidebar-step-group",
          tags$li(class = "sidebar-step-header", tags$span(class = "sidebar-step-name", "On this page")),
          dataset_row,
          lapply(Filter(function(k) TX_STEP_MAP[[k]]$id %in% shown_ids, names(TX_STEP_MAP)), open_row)
        ),
        add_btn
      ))
    }

    tagList(
      design_btn,
      lapply(WF_STEPS, function(st) {
        n <- st$n
        ## Only the analyses the design actually chose - an empty step is noise, not progress.
        keys <- intersect(tx_step_keys(n), wf_tx$selected)
        ## Steps 3 and 4 are setup, not analyses - once the design is fixed they are nothing to
        ## track. Step 1 stays because the Dataset tab has to remain reachable.
        if (!length(keys) && n != 1L) return(NULL)
        ## A step only goes green when the analyses in it have actually produced results; setup
        ## steps have no analyses, so they stay neutral.
        step_status <- wf_step_status(vapply(keys, function(k) sts[[k]]$status, character(1)))
        rows <- if (n == 1L) list(dataset_row)
                else lapply(Filter(function(k) TX_STEP_MAP[[k]]$id %in% shown_ids, keys), open_row)
        hint <- if (length(rows)) NULL else tags$li(class = "sidebar-step-empty", switch(as.character(n),
          "3" = "Set in the guide, or on the Differential Expression page.",
          "4" = "Set in the guide.",
          "Nothing added yet - use + Add analyses."))
        tags$ul(class = "sidebar-nav sidebar-step-group",
                omics_sidebar_step_header(n, st$label, step_status), rows, hint)
      }),
      add_btn
    )
  })
  outputOptions(output, "tx_sidebar_dynamic_nav", suspendWhenHidden = FALSE)

  ## Catalogue cards carry the same status, using the existing smstate_<id> label and sm-card-active
  ## class. A locked card still opens its tab; it only says what it is waiting for.
  observe({
    sts <- tx_statuses()
    ids <- added$ids
    all_states <- paste(paste0("sm-state-", names(WF_STATUS_LABEL)), collapse = " ")
    lapply(Filter(function(e) !e$alias, TX_STEP_MAP), function(e) {
      s <- sts[[e$key]]
      shinyjs::html(id = paste0("smstate_", e$id), html = if (e$id %in% ids) "Added" else "Add")
      shinyjs::html(id = paste0("smreason_", e$id), html = htmltools::htmlEscape(
        if (identical(s$status, "done")) "Done" else
        if (identical(s$status, "ready")) "Ready" else s$reason))
      shinyjs::removeClass(id = paste0("smcard_", e$id), class = all_states)
      shinyjs::addClass(id = paste0("smcard_", e$id), class = paste0("sm-state-", s$status))
    })
  })

  ## Which of the nine guide steps the page currently on screen belongs to - NULL when the open page
  ## (Custom Analysis Design, Sub-modules) spans more than one step and only the guide itself knows
  ## which one it is on. Shared by the breadcrumb and the header step strip so both agree with what
  ## is actually visible instead of the strip trailing behind wherever the guide was last advanced.
  tx_active_step <- reactive({
    key <- wf_tx$active_key %||% "__dataset__"
    if (identical(key, "__catalogue__") || identical(key, "__design__")) NULL
    else if (identical(key, "__dataset__") || is.null(TX_STEP_MAP[[key]])) 1L
    else TX_STEP_MAP[[key]]$step
  })

  ## Breadcrumb in the page header, in the slot the page subtitle used to fill.
  output$tx_page_subtitle <- renderUI({
    key <- wf_tx$active_key %||% "__dataset__"
    titles <- wf_step_titles()
    step <- tx_active_step()
    leaf <- if (identical(key, "__catalogue__")) "Add analyses"
            else if (identical(key, "__design__")) "Custom Analysis Design"
            else if (identical(key, "__dataset__") || is.null(TX_STEP_MAP[[key]])) "Dataset"
            else titles[[key]]
    crumb <- list(step = step, leaf = leaf)
    sep <- tags$span(class = "crumb-sep", HTML("&rsaquo;"))
    tags$nav(
      class = "page-breadcrumb",
      tags$a(href = "#", class = "crumb-link",
             onclick = "Shiny.setInputValue('tx_crumb_home', Math.random(), {priority: 'event'}); return false;",
             icon("house"), " Home"),
      sep, tags$span(class = "crumb-here", "Transcriptomics"),
      if (!is.null(crumb$step)) tagList(sep, tags$span(class = "crumb-step",
        sprintf("Step %d · %s", crumb$step, WF_STEP_LABEL[[crumb$step]]))),
      sep, tags$strong(class = "crumb-leaf", crumb$leaf)
    )
  })

  ## The header carries the step strip while the guide is running, so where you are is always on
  ## screen next to the breadcrumb. (This slot used to hold the read-only contrast line; the
  ## breadcrumb and the design page already state the contrast, so it was saying it a third time.)
  output$tx_settings_banner <- renderUI({
    if (!isTRUE(wf_guide$active)) return(NULL)
    ## The page on screen wins when it maps to one specific step (e.g. the Dataset tab is always
    ## step 1); on a page that spans several steps (Custom Analysis Design, Sub-modules) there is
    ## nothing to disagree with, so the guide's own position stands.
    step <- tx_active_step() %||% wf_guide$step
    div(class = "wf-header-steps", wf_bar_steps_ui(step, wf_guide$done))
  })
  outputOptions(output, "tx_settings_banner", suspendWhenHidden = FALSE)

  output$tx_sidebar_arthochat_hint <- renderUI({
    mod_id <- title_to_module_id(TX_MODULES, input$tx_menu)
    hint <- if (is.null(mod_id)) {
      "Ask ArthOChat about any module title or your dataset."
    } else {
      cfg <- TX_MODULES_BY_ID[[mod_id]]$config
      paste0("Ask ArthOChat about ", cfg$title, ".")
    }
    arthochat_shortcut_ui(hint, compact = TRUE)
  })
  outputOptions(output, "tx_sidebar_arthochat_hint", suspendWhenHidden = FALSE)

  output$mx_sidebar_arthochat_hint <- renderUI({
    mod_id <- title_to_module_id(MX_MODULES, input$mx_menu)
    hint <- if (is.null(mod_id)) {
      "Ask ArthOChat about any module title or your dataset."
    } else {
      cfg <- MX_MODULES_BY_ID[[mod_id]]$config
      paste0("Ask ArthOChat about ", cfg$title, ".")
    }
    arthochat_shortcut_ui(hint, compact = TRUE)
  })
  outputOptions(output, "mx_sidebar_arthochat_hint", suspendWhenHidden = FALSE)

  highlight_mx_sidebar <- function() {
    req(input$mx_menu)
    shinyjs::runjs(sprintf(
      "(function retry(n){
         var col = $('.omics-sidebar-col:visible');
         if (col.length) {
           col.find('.sidebar-nav-item').removeClass('active');
           col.find('.sidebar-nav-item[data-match=\"%s\"]').addClass('active');
         } else if (n > 0) {
           setTimeout(function(){ retry(n - 1); }, 50);
         }
       })(20);",
      gsub('(["\\\\])', "\\\\\\1", input$mx_menu)
    ))
  }
  observeEvent(input$mx_menu, highlight_mx_sidebar())
  observeEvent(input$sidebar_tabs, {
    if (identical(input$sidebar_tabs, "methylomics")) highlight_mx_sidebar()
  }, ignoreInit = TRUE)

  jump_to_mx_submodule <- function(mod_id, sm_filter = NULL) {
    cfg <- MX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% mx_added$ids) {
      updateTabsetPanel(session, "mx_menu", selected = cfg$title)
    } else {
      updateTabsetPanel(session, "mx_menu", selected = "Sub-modules")
      updateTextInput(session, "mx_sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_methylomics_dataset, {
    updateTabsetPanel(session, "mx_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_methylomics_celltype, {
    jump_to_mx_submodule("celltype", sm_filter = "Cell-Type Deconvolution")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_methylomics_submodules, {
    updateTabsetPanel(session, "mx_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  MX_DYNAMIC_NAV_IDS <- c("qc", "normalization", "dmp")
  lapply(MX_DYNAMIC_NAV_IDS, function(hid) {
    m <- MX_MODULES_BY_ID[[hid]]
    observeEvent(input[[paste0("sidebar_nav_methylomics_dyn_", hid)]], {
      updateTabsetPanel(session, "mx_menu", selected = m$config$title)
    }, ignoreInit = TRUE)
  })

  output$mx_sidebar_dynamic_nav <- renderUI({
    tagList(lapply(MX_DYNAMIC_NAV_IDS, function(hid) {
      if (!hid %in% mx_added$ids) return(NULL)
      m <- MX_MODULES_BY_ID[[hid]]
      tags$li(
        tags$a(
          id = paste0("sidebar_nav_methylomics_dyn_", hid), href = "#",
          class = "sidebar-nav-item action-button",
          `data-match` = m$config$title,
          icon(m$config$icon), m$config$title
        )
      )
    }))
  })
  outputOptions(output, "mx_sidebar_dynamic_nav", suspendWhenHidden = FALSE)

  jump_to_cx_submodule <- function(mod_id, sm_filter = NULL) {
    cfg <- CX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% cx_added$ids) {
      updateTabsetPanel(session, "cx_menu", selected = cfg$title)
    } else {
      updateTabsetPanel(session, "cx_menu", selected = "Sub-modules")
      updateTextInput(session, "cx_sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_crossomics_dataset, {
    updateTabsetPanel(session, "cx_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_crossomics_submodules, {
    updateTabsetPanel(session, "cx_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  lapply(CX_MODULES, function(m) {
    hid <- m$config$id
    observeEvent(input[[paste0("sidebar_nav_crossomics_dyn_", hid)]], {
      updateTabsetPanel(session, "cx_menu", selected = m$config$title)
    }, ignoreInit = TRUE)
  })

  output$cx_sidebar_dynamic_nav <- renderUI({
    tagList(lapply(CX_MODULES, function(m) {
      hid <- m$config$id
      if (!hid %in% cx_added$ids) return(NULL)
      tags$li(
        tags$a(
          id = paste0("sidebar_nav_crossomics_dyn_", hid), href = "#",
          class = "sidebar-nav-item action-button",
          `data-match` = m$config$title,
          icon(m$config$icon), m$config$title
        )
      )
    }))
  })
  outputOptions(output, "cx_sidebar_dynamic_nav", suspendWhenHidden = FALSE)

  jump_to_mo_submodule <- function(mod_id, sm_filter = NULL) {
    cfg <- MULTI_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% mo_added$ids) {
      updateTabsetPanel(session, "mo_menu", selected = cfg$title)
    } else {
      updateTabsetPanel(session, "mo_menu", selected = "Sub-modules")
      updateTextInput(session, "mo_sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_multiomics_dataset, {
    updateTabsetPanel(session, "mo_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_multiomics_submodules, {
    updateTabsetPanel(session, "mo_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  lapply(MULTI_MODULES, function(m) {
    hid <- m$config$id
    observeEvent(input[[paste0("sidebar_nav_multiomics_dyn_", hid)]], {
      updateTabsetPanel(session, "mo_menu", selected = m$config$title)
    }, ignoreInit = TRUE)
  })

  output$mo_sidebar_dynamic_nav <- renderUI({
    tagList(lapply(MULTI_MODULES, function(m) {
      hid <- m$config$id
      if (!hid %in% mo_added$ids) return(NULL)
      tags$li(
        tags$a(
          id = paste0("sidebar_nav_multiomics_dyn_", hid), href = "#",
          class = "sidebar-nav-item action-button",
          `data-match` = m$config$title,
          icon(m$config$icon), m$config$title
        )
      )
    }))
  })
  outputOptions(output, "mo_sidebar_dynamic_nav", suspendWhenHidden = FALSE)

  ## Guided "Start an analysis" workflow (R/workflow_guide.R): opens existing tabs the same way their
  ## sub-module "Add" cards and sidebar links do.
  wf_guide <- workflow_guide_server(
    input, output, session, dataset, results, methyl_dataset, methyl_results,
    nav = list(
      dataset = function(layer) {
        updateTabsetPanel(session, "sidebar_tabs", selected = layer)
        updateTabsetPanel(session, if (identical(layer, "transcriptomics")) "tx_menu" else "mx_menu", selected = "Dataset")
      },
      ## Cross-Omics / Multi-Omics: no guided pipeline, so the picker just opens the section.
      section = function(layer) updateTabsetPanel(session, "sidebar_tabs", selected = layer),
      apply = function(layer, keys) {
        if (!identical(layer, "transcriptomics")) return(invisible())
        ids <- unique(vapply(intersect(keys, names(TX_STEP_MAP)),
                             function(k) TX_STEP_MAP[[k]]$id, character(1)))
        for (hid in setdiff(ids, added$ids)) tx_insert_submodule(TX_MODULES_BY_ID[[hid]])
        first <- intersect(keys, names(TX_STEP_MAP))
        if (length(first)) {
          e <- TX_STEP_MAP[[first[[1]]]]
          updateTabsetPanel(session, "tx_menu", selected = TX_MODULES_BY_ID[[e$id]]$config$title)
          wf_tx$active_key <- first[[1]]
        }
      },
      design = function(layer) {
        updateTabsetPanel(session, "sidebar_tabs", selected = layer)
        updateTabsetPanel(session, "tx_menu", selected = "Custom Analysis Design")
        wf_tx$active_key <- "__design__"
      },
      ## Inserts a module's tab without switching to it, so its inputs exist for the guide to read.
      mount = function(layer, hid) {
        if (!identical(layer, "transcriptomics") || hid %in% added$ids) return(invisible())
        tx_insert_submodule(TX_MODULES_BY_ID[[hid]])
        updateTabsetPanel(session, "tx_menu", selected = "Custom Analysis Design")
      },
      open = function(layer, hid) {
        updateTabsetPanel(session, "sidebar_tabs", selected = layer)
        if (identical(layer, "transcriptomics")) {
          m <- TX_MODULES_BY_ID[[hid]]
          if (hid %in% added$ids) updateTabsetPanel(session, "tx_menu", selected = m$config$title) else tx_insert_submodule(m)
        } else {
          m <- MX_MODULES_BY_ID[[hid]]
          if (hid %in% mx_added$ids) updateTabsetPanel(session, "mx_menu", selected = m$config$title) else mx_insert_submodule(m)
        }
      }
    )
  )

  ## The guide owns the contrast, the sex design and the chosen analyses; the shell mirrors them so
  ## the sidebar, the cards and the read-only banner all use the one setting (transcriptomics only).
  observe({
    if (!identical(wf_guide$layer, "transcriptomics")) return()
    wf_tx$contrast <- wf_guide$contrast
    wf_tx$sex_mode <- wf_guide$sex_mode
    wf_tx$selected <- wf_guide$selected %||% character(0)
    wf_tx$configured <- isTRUE(wf_guide$configured)
  })

  observeEvent(input$header_search_submit, {
    q <- tolower(trimws(input$header_search_submit %||% ""))
    req(nzchar(q))

    mod_hit <- Find(function(m) identical(m$status, "available") && grepl(q, tolower(m$title), fixed = TRUE), MODULE_REGISTRY)
    if (!is.null(mod_hit)) {
      if (identical(mod_hit$id, "arthochat")) {
        shinyjs::runjs(ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT)
      } else {
        updateTabsetPanel(session, "sidebar_tabs", selected = mod_hit$tab)
      }
      return()
    }

    sm_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), TX_MODULES)
    if (!is.null(sm_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "transcriptomics")
      jump_to_submodule(sm_hit$config$id, sm_filter = sm_hit$config$title)
      return()
    }

    mx_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), MX_MODULES)
    if (!is.null(mx_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "methylomics")
      jump_to_mx_submodule(mx_hit$config$id, sm_filter = mx_hit$config$title)
      return()
    }

    cx_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), CX_MODULES)
    if (!is.null(cx_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "crossomics")
      jump_to_cx_submodule(cx_hit$config$id, sm_filter = cx_hit$config$title)
      return()
    }

    mo_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), MULTI_MODULES)
    if (!is.null(mo_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "multiomics")
      jump_to_mo_submodule(mo_hit$config$id, sm_filter = mo_hit$config$title)
      return()
    }

    showNotification(sprintf('No module or sub-module matched "%s".', input$header_search_submit), type = "warning")
  }, ignoreInit = TRUE)

  observeEvent(input$theme_toggle_btn, {
    shinyjs::runjs("
      var html = document.documentElement;
      var next = html.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
      html.setAttribute('data-theme', next);
      try { localStorage.setItem('arthomix-theme', next); } catch (e) {}
    ")
  }, ignoreInit = TRUE)

  output$mx_page_subtitle <- renderUI({
    sel <- input$mx_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = "Load a methylation dataset",
      "Sub-modules" = "Add or remove sub-modules.",
      sel
    )
    p(txt)
  })

  output$cx_page_subtitle <- renderUI({
    sel <- input$cx_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = "Load a DEG/DMP data",
      "Sub-modules" = "Add or remove sub-modules.",
      sel
    )
    p(txt)
  })

  output$mo_page_subtitle <- renderUI({
    sel <- input$mo_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = NULL,
      "Sub-modules" = "Add or remove sub-modules.",
      sel
    )
    if (is.null(txt)) return(NULL)
    p(txt)
  })
}
