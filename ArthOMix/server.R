## server.R
## ArthOMix Explorer

function(input, output, session) {

  dataset <- local({
    d <- load_default_dataset()
    d$source_type <- "preloaded"
    d$is_bundled_reference <- TRUE
    d$geo_ids <- MERGED_DEFAULT_GEO_IDS
    d$declared_data_type <- NA_character_
    do.call(reactiveValues, d)
  })

  results <- reactiveValues()

  ## Session-wide analysis-record + diagnostics stores (R/provenance.R). Every analysis run button in the
  ## Transcriptomics module pushes its provenance record here (Methylomics / Cross-Omics / Multi-Omics are not
  ## yet wired; the methylomics DMP tab writes its own downloadable manifest); the header "Analysis records"
  ## button lists and exports them, and the same panel shows the reason behind every error/warning the app
  ## absorbed into a "Not available".
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
        " Coverage: every analysis run button in the Transcriptomics module writes a record here (Overview normalisation, Preprocessing, Differential Expression, WGCNA, Candidate Genes, MR, Colocalisation, Feature Selection, Diagnostic Model incl. External Validation, Sex Interaction, Cross-Tissue Replication, Cross-Ancestry MR Replication, Functional Enrichment, Immune Deconvolution, Nomogram, Biomarker Card). The Methylomics, Cross-Omics and Multi-Omics modules are not yet wired to this log; the methylomics DMP tab writes its own downloadable manifest."),
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
    hint <- if (identical(input$cx_menu, "Cross-Omics MR")) {
      "Need help selecting MR data - which sex's evidence, or preloaded vs. your own upload? Ask ArthOChat."
    } else {
      "Questions about how the panels converge, or which sex/data source to select? Ask ArthOChat."
    }
    arthochat_shortcut_ui(hint, compact = TRUE)
  })

  mod_arthochat_server(
    "arthochat", dataset, results,
    methyl_dataset, methyl_results,
    cross_dataset, cross_results,
    multi_dataset, multi_results,
    current_context = current_module_context,
    run_hooks = agent_run_hooks
  )

  tx_server_hooks <- setNames(
    lapply(TX_MODULES, function(m) m$server(paste0("tx_", m$config$id), dataset, results)),
    vapply(TX_MODULES, function(m) m$config$id, character(1))
  )
  agent_run_hooks$transcriptomics$dge <- tx_server_hooks$dge$run

  added <- reactiveValues(ids = character(0))

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
        insertTab(
          session = session, inputId = "tx_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("tx_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        added$ids <- union(added$ids, hid)
        shinyjs::addClass(id = paste0("smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("smstate_", hid), html = "Added")
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
        insertTab(
          session = session, inputId = "mx_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("mx_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        mx_added$ids <- union(mx_added$ids, hid)
        shinyjs::addClass(id = paste0("mx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mx_smstate_", hid), html = "Added")
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

  jump_to_submodule <- function(mod_id, inner_tab = NULL, sm_filter = NULL) {
    cfg <- TX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% added$ids) {
      updateTabsetPanel(session, "tx_menu", selected = cfg$title)
      if (!is.null(inner_tab)) {
        updateTabsetPanel(session, paste0("tx_", mod_id, "-tabs"), selected = inner_tab)
      }
    } else {
      updateTabsetPanel(session, "tx_menu", selected = "Sub-modules")
      updateTextInput(session, "sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_transcriptomics_dataset, {
    updateTabsetPanel(session, "tx_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_transcriptomics_submodules, {
    updateTabsetPanel(session, "tx_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  highlight_tx_sidebar <- function() {
    req(input$tx_menu)
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
      gsub('(["\\\\])', "\\\\\\1", input$tx_menu)
    ))
  }
  observeEvent(input$tx_menu, highlight_tx_sidebar())
  observeEvent(input$sidebar_tabs, {
    if (identical(input$sidebar_tabs, "transcriptomics")) highlight_tx_sidebar()
  }, ignoreInit = TRUE)

  lapply(TX_MODULES, function(m) {
    hid <- m$config$id
    observeEvent(input[[paste0("sidebar_nav_transcriptomics_dyn_", hid)]], {
      updateTabsetPanel(session, "tx_menu", selected = m$config$title)
    }, ignoreInit = TRUE)
  })

  output$tx_sidebar_dynamic_nav <- renderUI({
    tagList(lapply(TX_MODULES, function(m) {
      hid <- m$config$id
      if (!hid %in% added$ids) return(NULL)
      tags$li(
        tags$a(
          id = paste0("sidebar_nav_transcriptomics_dyn_", hid), href = "#",
          class = "sidebar-nav-item action-button",
          `data-match` = m$config$title,
          icon(m$config$icon), m$config$title
        )
      )
    }))
  })
  outputOptions(output, "tx_sidebar_dynamic_nav", suspendWhenHidden = FALSE)

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

  output$tx_page_subtitle <- renderUI({
    sel <- input$tx_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = "Load a transcriptomics dataset",
      "Sub-modules" = "Add or remove sub-modules.",
      sel
    )
    p(txt)
  })

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
