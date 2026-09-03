## R/multiomics/01_Data_Workspace/mod_multi_dataset.R
## Multi-Omics Dataset Workspace - the front door for the Multi-Omics module.
## Lets the user pick a preloaded dataset, upload their own, or retrieve one

mod_multi_dataset_config <- list(
  id = "dataset", title = "Dataset Workspace", icon = "database",
  description = "Select, upload, or retrieve datasets for multi-omics analysis."
)

MO_MAX_BLOCKS <- 8
MO_BLOCK_IDS <- paste0("block", seq_len(MO_MAX_BLOCKS))

MO_SOURCE_CHOICES <- c("Reference / Example Dataset" = "preloaded", "Upload Dataset" = "upload", "Retrieve from GEO" = "geo")

mo_status_badge <- function(status) {
  color <- switch(status$level, ready = ARTHOMIX_COLORS$aqua, review = ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)
  icn <- switch(status$level, ready = "circle-check", review = "triangle-exclamation", "circle-xmark")
  tagList(
    span(style = sprintf("color:%s; font-weight:600;", color), icon(icn), " ", status$label),
    if (length(status$reasons) > 0) tags$ul(style = "margin:4px 0 0 0; padding-left:18px; font-size:0.85em; color:var(--color-ink-muted, #898781);",
                                              lapply(status$reasons, tags$li))
  )
}

mo_load_first_msg <- function(dataset_source) {
  if (identical(dataset_source, "preloaded")) "Click \"Load Reference Dataset\" above first." else "Add at least two datasets and click \"Validate Datasets\" first."
}

mo_dataset_block_card <- function(label, n_samples, n_features, status) {
  color <- switch(status$level, ready = ARTHOMIX_COLORS$aqua, review = ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)
  div(class = "card", style = sprintf("flex:1 1 220px; padding:14px; border-left:4px solid %s;", color),
      div(style = "font-weight:600; margin-bottom:6px;", label),
      div(style = "font-size:0.85em; color:var(--color-ink-muted, #898781);",
          sprintf("Samples: %s", format(n_samples, big.mark = ","))),
      div(style = "font-size:0.85em; color:var(--color-ink-muted, #898781);",
          sprintf("Features: %s", format(n_features, big.mark = ","))),
      div(style = "margin-top:6px;", mo_status_badge(status))
  )
}

mo_summary_table <- function(layer_meta) {
  if (length(layer_meta) == 0) return(NULL)
  do.call(rbind, lapply(names(layer_meta), function(nm) {
    lm <- layer_meta[[nm]]
    v <- lm$validation
    data.frame(
      Dataset = nm,
      `Omics Type` = names(MULTI_LIVE_OMICS_TYPES)[match(lm$omics_type, MULTI_LIVE_OMICS_TYPES)] %||% lm$omics_type,
      Samples = if (!is.null(v)) v$n_samples else NA,
      Features = if (!is.null(v)) v$n_features else NA,
      Processing = lm$processing %||% "Not processed",
      Status = lm$status$label %||% "Unknown",
      check.names = FALSE
    )
  }))
}

mo_provenance_ui <- function(layer_meta) {
  if (length(layer_meta) == 0) return(multi_empty_state("No datasets selected yet."))
  tagList(lapply(names(layer_meta), function(nm) {
    p <- layer_meta[[nm]]$provenance %||% list()
    box(width = NULL, title = nm, status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
        tags$table(class = "table table-condensed", style = "font-size:0.88em;",
                    tags$tbody(
                      tags$tr(tags$td(tags$strong("Source")), tags$td(p$source %||% "Unknown")),
                      tags$tr(tags$td(tags$strong("Detail")), tags$td(p$detail %||% "")),
                      tags$tr(tags$td(tags$strong("Omics")), tags$td(names(MULTI_LIVE_OMICS_TYPES)[match(layer_meta[[nm]]$omics_type, MULTI_LIVE_OMICS_TYPES)] %||% "")),
                      tags$tr(tags$td(tags$strong("Imported")), tags$td(p$imported_at %||% ""))
                    )))
  }))
}

mo_preloaded_blocks_ui <- function() {
  if (!MULTI_DATA_AVAILABLE) {
    return(div(class = "empty-note", icon("triangle-exclamation"),
               "The preloaded multi-omics dataset isn't available in this deployment."))
  }
  rna <- multi_read_registry_table("RNA-seq QC summary")
  meth <- multi_read_registry_table("Methylation QC summary")
  matching <- multi_read_registry_table("Patient sample matching (all 80 patients)")

  rna_pbmc <- if (rna$ok && all(c("cell_type", "n_samples", "n_genes_retained") %in% colnames(rna$df))) rna$df[rna$df$cell_type == "PBMC", , drop = FALSE] else NULL
  rna_n_samples <- if (!is.null(rna_pbmc) && nrow(rna_pbmc) > 0) rna_pbmc$n_samples[1] else NA
  rna_n_features <- if (!is.null(rna_pbmc) && nrow(rna_pbmc) > 0) rna_pbmc$n_genes_retained[1] else NA
  meth_n_samples <- if (meth$ok && "n_samples_retained" %in% colnames(meth$df)) meth$df$n_samples_retained[1] else NA
  meth_n_features <- if (meth$ok && "n_probes_retained" %in% colnames(meth$df)) meth$df$n_probes_retained[1] else NA
  n_total <- if (matching$ok) nrow(matching$df) else NA

  ready_status <- list(level = "ready", label = "Ready", reasons = character(0))
  tagList(
    div(style = "display:flex; gap:14px; flex-wrap:wrap;",
        if (rna$ok) mo_dataset_block_card("Transcriptomics", rna_n_samples, rna_n_features, ready_status),
        if (meth$ok) mo_dataset_block_card("Methylomics", meth_n_samples, meth_n_features, ready_status)
    ),
    box(width = NULL, title = "Study Context", status = "primary", solidHeader = FALSE, style = "margin-top:12px;",
        tags$table(class = "table table-condensed", style = "font-size:0.9em;",
                    tags$tbody(
                      tags$tr(tags$td(tags$strong("Study")), tags$td("Rheumatoid arthritis anti-TNF treatment study (Tao et al. 2021)")),
                      tags$tr(tags$td(tags$strong("Samples")), tags$td(if (!is.na(n_total)) format(n_total, big.mark = ",") else "Unavailable")),
                      tags$tr(tags$td(tags$strong("Disease")), tags$td("Rheumatoid arthritis")),
                      tags$tr(tags$td(tags$strong("Time point")), tags$td("Baseline, pre-treatment"))
                    )))
  )
}

mo_block_id <- function(i, mode) paste0(if (identical(mode, "geo")) "g" else "u", "block", i)

mo_file_input_id <- function(bid, gen) paste0(bid, "_file_g", gen)
mo_meta_file_input_id <- function(bid, gen) paste0(bid, "_meta_file_g", gen)

mo_block_ui <- function(ns, i, mode = c("upload", "geo"), gen = 0) {
  mode <- match.arg(mode)
  bid <- mo_block_id(i, mode)
  box(
    width = NULL, title = sprintf("Dataset %d%s", i, if (i <= 2) " (required)" else ""),
    status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
    selectInput(ns(paste0(bid, "_type")), "Omics type", choices = MULTI_LIVE_OMICS_TYPES,
                selected = if (i == 1) "rnaseq" else if (i == 2) "methylation" else "other"),
    textInput(ns(paste0(bid, "_label")), "Display label", value = sprintf("Dataset %d", i)),
    uiOutput(ns(paste0(bid, "_feature_id_note"))),
    if (identical(mode, "upload")) tagList(
      fileInput(ns(mo_file_input_id(bid, gen)), "File (CSV, TSV, TXT, XLSX, or RDS)",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".Rds")),
      uiOutput(ns(paste0(bid, "_omics_type_note"))),
      fileInput(ns(mo_meta_file_input_id(bid, gen)), "Sample metadata for this dataset (optional, CSV, first column = sample ID)",
                accept = c(".csv")),
      uiOutput(ns(paste0(bid, "_orient_note"))),
      uiOutput(ns(paste0(bid, "_shape_ui")))
    ) else tagList(
      textInput(ns(paste0(bid, "_geo_acc")), "GEO accession", placeholder = "GSE12345"),
      actionButton(ns(paste0(bid, "_geo_fetch")), "Fetch from GEO", icon = icon("cloud-arrow-down"), class = "btn-primary btn-sm"),
      fileInput(ns(mo_meta_file_input_id(bid, gen)), "Additional sample metadata for this dataset (optional, CSV, first column = sample ID)",
                accept = c(".csv")),
      uiOutput(ns(paste0(bid, "_geo_platform_ui"))),
      uiOutput(ns(paste0(bid, "_geo_status")))
    ),
    p(class = "submodule-desc", "Rows must match the sample IDs in your metadata.")
  )
}

mod_multi_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      width = NULL, title = "Select Dataset Source", status = "primary", solidHeader = FALSE,
      radioButtons(ns("dataset_source"), NULL, choices = MO_SOURCE_CHOICES, selected = "preloaded", inline = TRUE)
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'preloaded'", ns("dataset_source")),
      box(width = NULL, title = "Reference / Example Dataset", status = "primary", solidHeader = FALSE,
          selectInput(ns("preloaded_pick"), "Select a reference dataset",
                      choices = c("RA anti-TNF Multi-Omics Dataset" = "ra_antitnf"), width = "100%"),
          selectInput(ns("preloaded_cell"), "Analysis cell (matched sex x drug/outcome subset)", choices = MULTI_CELL_CHOICES, width = "100%"),
          actionButton(ns("load_preloaded_btn"), "Load Reference Dataset", icon = icon("database"), class = "btn-primary btn-sm"),
          div(style = "margin-top:10px;", uiOutput(ns("preloaded_blocks_ui"))),
          div(style = "margin-top:10px;", uiOutput(ns("preloaded_activate_ui")))
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("load_preloaded_btn")),
        box(width = NULL, title = "Table", status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
            if (!MULTI_DATA_AVAILABLE) div(class = "empty-note", icon("triangle-exclamation"), "Not available in this deployment.")
            else tagList(
              selectInput(ns("table_pick"), NULL, choices = names(MULTI_TABLE_REGISTRY), width = "100%"),
              actionButton(ns("load_table_btn"), "Load", icon = icon("upload"), class = "btn-primary btn-sm")
            ))
      )
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload' || input['%s'] == 'geo' || input['%s'] == 'preloaded'", ns("dataset_source"), ns("dataset_source"), ns("dataset_source")),
      fluidRow(
        column(
          4,
          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload'", ns("dataset_source")),
            uiOutput(ns("upload_blocks_ui")),
            uiOutput(ns("add_upload_block_ui"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'geo'", ns("dataset_source")),
            div(class = "empty-note", icon("circle-info"),
                "Enter one GEO Series accession per dataset below. Inspect what's fetched before adding it - not every layer will be compatible."),
            uiOutput(ns("geo_blocks_ui")),
            uiOutput(ns("add_geo_block_ui")),
            box(width = NULL, title = "Sample Metadata", status = "primary", solidHeader = FALSE,
                p(class = "submodule-desc", "Imported automatically from each fetched GEO series' own sample metadata."))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload' || input['%s'] == 'geo'", ns("dataset_source"), ns("dataset_source")),
            actionButton(ns("validate_btn"), "Validate Datasets", icon = icon("check-double"), class = "btn-primary btn-sm", width = "100%")
          )
        ),
        column(8, uiOutput(ns("pipeline_ui")))
      ),

      conditionalPanel(
        condition = sprintf("input['%s'] > 0 || input['%s'] > 0", ns("validate_btn"), ns("load_preloaded_btn")),
        hr(),
        box(width = NULL, title = "Data Provenance", status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
            uiOutput(ns("provenance_ui")))
      ),

      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("activate_btn")),
        hr(),
        box(width = NULL, title = mod_multi_mofa_config$title, status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", mod_multi_mofa_config$description),
            mod_multi_mofa_ui(ns("integrated")))
      )
    )
  )
}

mo_label_omics_type <- function(label, input, n_upload, n_geo, mode) {
  if (identical(mode, "preloaded")) {
    return(switch(label, Transcriptomics = "rnaseq", Methylomics = "methylation", "other"))
  }
  ids <- if (identical(mode, "upload")) seq_len(n_upload) else seq_len(n_geo)
  bmode <- if (identical(mode, "upload")) "upload" else "geo"
  for (i in ids) {
    bid <- mo_block_id(i, bmode)
    if (identical(input[[paste0(bid, "_label")]], label)) return(input[[paste0(bid, "_type")]] %||% "other")
  }
  "other"
}

mo_apply_matching <- function(mats, method, meta = NULL, patient_col = NULL, mapping_df = NULL) {
  if (identical(method, "patient_id") && !is.null(meta) && !is.null(patient_col) && patient_col %in% colnames(meta)) {
    meta_ids <- rownames(meta)
    meta_norm <- ch_normalize_id(meta_ids)
    meta_norm_dup <- duplicated(meta_norm) | duplicated(meta_norm, fromLast = TRUE)
    dropped <- list()
    mats <- Map(function(m, label) {
      rn <- rownames(m)
      exact_idx <- match(rn, meta_ids)
      rn_norm <- ch_normalize_id(rn)
      norm_idx <- match(rn_norm, meta_norm)
      use_norm <- is.na(exact_idx) & !is.na(norm_idx) & !meta_norm_dup[ifelse(is.na(norm_idx), 1L, norm_idx)]
      final_idx <- ifelse(!is.na(exact_idx), exact_idx, ifelse(use_norm, norm_idx, NA_integer_))
      keep <- !is.na(final_idx)
      dropped[[label]] <<- sum(!keep)
      out <- m[keep, , drop = FALSE]
      rownames(out) <- as.character(meta[[patient_col]][final_idx[keep]])
      out
    }, mats, names(mats))
    return(list(mats = mats, dropped = dropped))
  }
  if (identical(method, "mapping") && !is.null(mapping_df) && ncol(mapping_df) >= 2) {
    canonical <- as.character(mapping_df[[1]])
    dropped <- list()
    mats <- Map(function(m, label) {
      col <- colnames(mapping_df)[tolower(colnames(mapping_df)) %in% tolower(c(label, gsub("[^A-Za-z0-9]", "_", label)))]
      if (length(col) == 0) { dropped[[label]] <<- nrow(m); return(m) }
      map_ids <- as.character(mapping_df[[col[1]]])
      map_norm <- ch_normalize_id(map_ids)
      map_norm_dup <- duplicated(map_norm) | duplicated(map_norm, fromLast = TRUE)
      rn <- rownames(m)
      exact_idx <- match(rn, map_ids)
      rn_norm <- ch_normalize_id(rn)
      norm_idx <- match(rn_norm, map_norm)
      use_norm <- is.na(exact_idx) & !is.na(norm_idx) & !map_norm_dup[ifelse(is.na(norm_idx), 1L, norm_idx)]
      final_idx <- ifelse(!is.na(exact_idx), exact_idx, ifelse(use_norm, norm_idx, NA_integer_))
      keep <- !is.na(final_idx)
      dropped[[label]] <<- sum(!keep)
      out <- m[keep, , drop = FALSE]
      rownames(out) <- canonical[final_idx[keep]]
      out
    }, mats, names(mats))
    return(list(mats = mats, dropped = dropped))
  }
  list(mats = mats, dropped = list())
}

mo_merge_sample_meta <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b) || nrow(b) == 0) return(a)
  all_ids <- union(rownames(a), rownames(b))
  out <- data.frame(row.names = all_ids)
  for (cn in colnames(a)) out[[cn]] <- a[match(all_ids, rownames(a)), cn]
  for (cn in colnames(b)) out[[cn]] <- b[match(all_ids, rownames(b)), cn]
  out
}

mo_read_meta_file <- function(fi) {
  if (is.null(fi)) return(NULL)
  m <- tryCatch(as.data.frame(data.table::fread(fi$datapath, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(m) || ncol(m) < 1) return(NULL)
  rownames(m) <- as.character(m[[1]])
  m
}

mod_multi_dataset_server <- function(id, multi_dataset, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$preloaded_blocks_ui <- renderUI({
      if (!isTRUE(input$load_preloaded_btn > 0)) return(multi_empty_state("Click \"Load Dataset\" to preview this dataset."))
      mo_preloaded_blocks_ui()
    })

    loaded_table <- eventReactive(input$load_table_btn, {
      req(input$table_pick)
      multi_read_registry_table(input$table_pick)
    }, ignoreInit = TRUE)

    observeEvent(input$load_table_btn, {
      res <- loaded_table()
      if (isTRUE(res$ok)) multi_dataset$table_label <- input$table_pick
    }, ignoreInit = TRUE)

    preloaded_load_attempted <- reactiveVal(FALSE)

    output$preloaded_activate_ui <- renderUI({
      req(preloaded_load_attempted(), identical(input$dataset_source, "preloaded"))
      if (length(raw$mats) == 0) {
        return(div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
                    " Could not load this analysis cell - see the notification for details."))
      }
      div(class = "empty-note", icon("circle-check"),
          tags$strong(" Reference / Example Dataset loaded: "), "RA anti-TNF Multi-Omics Dataset (Transcriptomics + Methylomics). ",
          "Use the pipeline below (1. Preview and Validate through 5. Compatibility and Activate) - the same steps used for an uploaded dataset.")
    })

    observeEvent(input$load_preloaded_btn, {
      req(identical(input$dataset_source, "preloaded"))
      preloaded_load_attempted(TRUE)
      multi_dataset$source <- "preloaded"
      multi_dataset$active <- TRUE
      multi_dataset$layers <- list()
      multi_dataset$layer_meta <- list()
      multi_dataset$sample_meta <- NULL
      multi_dataset$overlap <- NULL
      multi_dataset$loaded_at <- Sys.time()

      res <- mi_preloaded_cell_dataset(input$preloaded_cell)
      if (!isTRUE(res$ok)) {
        showNotification(res$error, type = "error")
        raw$mats <- list(); raw$validations <- list(); raw$labels <- list(); raw$provenance <- list(); raw$meta <- NULL
        return()
      }
      raw$mats <- res$layers
      raw$validations <- stats::setNames(lapply(names(res$layers), function(nm) multi_live_validate_matrix(res$layers[[nm]], layer_label = nm)), names(res$layers))
      raw$labels <- stats::setNames(as.list(names(res$layers)), names(res$layers))
      raw$provenance <- stats::setNames(lapply(names(res$layers), function(nm) list(source = "Reference / Example Dataset", detail = res$provenance, imported_at = format(Sys.time(), "%d %b %Y %H:%M"))), names(res$layers))
      raw$meta <- res$sample_meta
    }, ignoreInit = TRUE)

    observeEvent(input$dataset_source, {
      multi_dataset$source <- NULL
      multi_dataset$active <- FALSE
      multi_dataset$layers <- list()
      multi_dataset$layer_meta <- list()
      multi_dataset$sample_meta <- NULL
      multi_dataset$overlap <- NULL
      multi_dataset$loaded_at <- NULL

      raw$mats <- list(); raw$validations <- list(); raw$labels <- list()
      raw$provenance <- list(); raw$meta <- NULL
      proc$filtered_mats <- NULL; proc$scaled_mats <- NULL; proc$batch_corrected <- NULL

      for (nm in names(geo_raw_platforms)) geo_raw_platforms[[nm]] <- NULL
      for (nm in names(geo_fetched)) geo_fetched[[nm]] <- NULL
      for (nm in names(block_shape)) block_shape[[nm]] <- NULL

      for (i in seq_len(MO_MAX_BLOCKS)) {
        ubid <- mo_block_id(i, "upload")
        output[[paste0(ubid, "_orient_note")]] <- renderUI(NULL)
        output[[paste0(ubid, "_omics_type_note")]] <- renderUI(NULL)
      }

      n_upload_blocks(2)
      n_geo_blocks(2)
      block_reset_gen(block_reset_gen() + 1)

      preloaded_load_attempted(FALSE)
    })

    n_upload_blocks <- reactiveVal(2)
    n_geo_blocks <- reactiveVal(2)
    block_reset_gen <- reactiveVal(0)
    geo_raw_platforms <- reactiveValues()   # gblockN -> list(accession, platforms)
    geo_fetched <- reactiveValues()          # gblockN -> list(mat, meta, platform, accession, collapsed)
    raw <- reactiveValues(mats = list(), validations = list(), labels = list(), provenance = list(), meta = NULL)
    proc <- reactiveValues(filtered_mats = NULL, scaled_mats = NULL, batch_corrected = NULL)
    block_shape <- reactiveValues()

    output$pipeline_ui <- renderUI({
      if (length(raw$validations) == 0) return(multi_empty_state(mo_load_first_msg(input$dataset_source)))
      gen <- block_reset_gen()
      tabsetPanel(
        id = ns("pipeline_tabs"), type = "tabs",
        tabPanel("1. Preview and Validate", br(), uiOutput(ns("validate_ui"))),
        tabPanel(
          "2. Sample Matching", br(),
          box(width = NULL, title = "Sample Matching Method", status = "primary", solidHeader = FALSE, collapsible = TRUE,
              radioButtons(ns("matching_method"), "Sample matching method",
                           choices = c("Exact Sample ID" = "exact", "Patient ID (from metadata)" = "patient_id", "User-defined mapping file" = "mapping"),
                           selected = "exact"),
              conditionalPanel(condition = sprintf("input['%s'] == 'patient_id'", ns("matching_method")),
                                uiOutput(ns("patient_id_col_ui"))),
              conditionalPanel(condition = sprintf("input['%s'] == 'mapping'", ns("matching_method")),
                                fileInput(ns(paste0("mapping_file_g", gen)), "Sample mapping file (CSV, one column per dataset's own sample IDs)", accept = c(".csv")))
          ),
          uiOutput(ns("matching_ui"))
        ),
        tabPanel("3. Preprocessing", br(), uiOutput(ns("preprocess_ui"))),
        tabPanel("4. Batch Diagnostics", br(), uiOutput(ns("batch_ui"))),
        tabPanel("5. Compatibility and Activate", br(), uiOutput(ns("compat_ui")))
      )
    })

    output$upload_blocks_ui <- renderUI({
      gen <- block_reset_gen()
      tagList(lapply(seq_len(n_upload_blocks()), function(i) mo_block_ui(ns, i, mode = "upload", gen = gen)))
    })
    output$geo_blocks_ui <- renderUI({
      gen <- block_reset_gen()
      tagList(lapply(seq_len(n_geo_blocks()), function(i) mo_block_ui(ns, i, mode = "geo", gen = gen)))
    })

    output$add_upload_block_ui <- renderUI({
      if (n_upload_blocks() >= MO_MAX_BLOCKS) return(div(class = "submodule-desc", sprintf("Maximum of %d datasets per session.", MO_MAX_BLOCKS)))
      actionButton(ns("add_upload_block_btn"), "+ Add Another Dataset", icon = icon("plus"), class = "btn-primary btn-sm")
    })
    observeEvent(input$add_upload_block_btn, n_upload_blocks(min(MO_MAX_BLOCKS, n_upload_blocks() + 1)))

    output$add_geo_block_ui <- renderUI({
      if (n_geo_blocks() >= MO_MAX_BLOCKS) return(div(class = "submodule-desc", sprintf("Maximum of %d datasets per session.", MO_MAX_BLOCKS)))
      actionButton(ns("add_geo_block_btn"), "+ Add Another Dataset", icon = icon("plus"), class = "btn-primary btn-sm")
    })
    observeEvent(input$add_geo_block_btn, n_geo_blocks(min(MO_MAX_BLOCKS, n_geo_blocks() + 1)))

    lapply(seq_len(MO_MAX_BLOCKS), function(i) {
      ubid <- mo_block_id(i, "upload"); gbid <- mo_block_id(i, "geo")

      output[[paste0(ubid, "_feature_id_note")]] <- renderUI({
        otype <- input[[paste0(ubid, "_type")]] %||% "other"
        p(class = "submodule-desc", sprintf("Feature identifier: %s", MULTI_LIVE_FEATURE_ID_LABELS[[otype]] %||% "Feature ID"))
      })
      output[[paste0(gbid, "_feature_id_note")]] <- renderUI({
        otype <- input[[paste0(gbid, "_type")]] %||% "other"
        p(class = "submodule-desc", sprintf("Feature identifier: %s", MULTI_LIVE_FEATURE_ID_LABELS[[otype]] %||% "Feature ID"))
      })

      observeEvent(input[[mo_file_input_id(ubid, block_reset_gen())]], {
        fi <- input[[mo_file_input_id(ubid, block_reset_gen())]]
        req(fi)
        ext <- tolower(tools::file_ext(fi$name))
        df <- if (ext %in% c("csv", "tsv", "txt")) tryCatch(as.data.frame(data.table::fread(fi$datapath, showProgress = FALSE, nrows = 200)), error = function(e) NULL) else NULL

        orient_det <- multi_live_detect_orientation(df)
        shape_det <- multi_live_detect_table_shape(df)
        long_det <- if (!is.null(df) && identical(shape_det$shape, "long")) multi_live_detect_long_columns(df) else NULL

        if (isTRUE(orient_det$confident)) {
          output[[paste0(ubid, "_orient_note")]] <- renderUI(
            div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
                sprintf(" %s Orientation below has been set accordingly - confirm or change it before validating.", orient_det$reason)))
        } else {
          output[[paste0(ubid, "_orient_note")]] <- renderUI(NULL)
        }

        otype_sel <- input[[paste0(ubid, "_type")]] %||% "other"
        preview_mat <- tryCatch({
          if (is.null(df) || ncol(df) < 2) NULL
          else if (identical(orient_det$suggested %||% "samples_rows", "features_rows")) {
            m <- t(as.matrix(vapply(df[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x)), numeric(nrow(df)))))
            colnames(m) <- as.character(df[[1]]); m
          } else {
            m <- as.matrix(vapply(df[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x)), numeric(nrow(df))))
            colnames(m) <- colnames(df)[-1]; m
          }
        }, error = function(e) NULL)
        otype_det <- if (!is.null(preview_mat)) multi_live_detect_omics_type(preview_mat) else NULL
        output[[paste0(ubid, "_omics_type_note")]] <- renderUI({
          if (is.null(otype_det) || !otype_sel %in% c("rnaseq", "methylation")) return(NULL)
          mismatch <- (identical(otype_sel, "rnaseq") && identical(otype_det$detected, "methylation")) ||
            (identical(otype_sel, "methylation") && identical(otype_det$detected, "rnaseq")) ||
            (otype_sel %in% c("rnaseq", "methylation") && identical(otype_det$detected, "unclassifiable"))
          if (!mismatch) return(NULL)
          div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
              sprintf(" Structural check: %s Selected omics type above is %s - please confirm this is correct.", otype_det$reason,
                      names(MULTI_LIVE_OMICS_TYPES)[match(otype_sel, MULTI_LIVE_OMICS_TYPES)] %||% otype_sel))
        })

        block_shape[[ubid]] <- list(
          shape = shape_det$shape, shape_reason = shape_det$reason,
          orient_selected = orient_det$suggested, cols = if (!is.null(df)) colnames(df) else character(0),
          long_det = long_det
        )
      }, ignoreInit = TRUE)

      output[[paste0(ubid, "_shape_ui")]] <- renderUI({
        bs <- block_shape[[ubid]]
        wide_ui <- radioButtons(ns(paste0(ubid, "_orient")), "Data orientation",
                                 choices = c("Sample and Feature (first column = sample ID)" = "samples_rows",
                                             "Feature and Sample (first column = feature ID)" = "features_rows"),
                                 selected = bs$orient_selected %||% "samples_rows")
        if (is.null(bs)) return(wide_ui)
        tagList(
          radioButtons(ns(paste0(ubid, "_table_shape")), "Table shape",
                       choices = c("Wide (samples/features matrix)" = "wide", "Long/tidy (one row per measurement)" = "long"),
                       selected = bs$shape %||% "wide", inline = TRUE),
          if (!is.null(bs$shape_reason)) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"), sprintf(" %s", bs$shape_reason)) else NULL,
          conditionalPanel(condition = sprintf("input['%s'] == 'wide'", ns(paste0(ubid, "_table_shape"))), wide_ui),
          conditionalPanel(condition = sprintf("input['%s'] == 'long'", ns(paste0(ubid, "_table_shape"))), {
            det <- bs$long_det
            col_choices <- c("(none)", bs$cols %||% character(0))
            pick <- function(field_id, label, sel) selectInput(ns(paste0(ubid, field_id)), label, choices = col_choices, selected = sel %||% "(none)")
            tagList(
              pick("_long_feature", "Feature ID column (gene/protein/probe/compound)", det$feature_col),
              pick("_long_sample", "Sample ID column", det$sample_col),
              pick("_long_value", "Measurement value column", det$value_col),
              pick("_long_group", "Group/Condition column (optional)", det$group_col),
              if (length(det$warnings %||% character(0)) > 0) div(class = "empty-note", icon("circle-info"), tags$ul(lapply(det$warnings, tags$li))) else NULL
            )
          })
        )
      })

      observeEvent(input[[paste0(gbid, "_geo_fetch")]], {
        acc <- input[[paste0(gbid, "_geo_acc")]]
        res <- multi_geo_layer_fetch(acc)
        if (!res$ok) {
          geo_fetched[[gbid]] <- NULL
          geo_raw_platforms[[gbid]] <- NULL
          output[[paste0(gbid, "_geo_status")]] <- renderUI(div(class = "empty-note", icon("triangle-exclamation"), res$error))
          output[[paste0(gbid, "_geo_platform_ui")]] <- renderUI(NULL)
          return()
        }
        geo_raw_platforms[[gbid]] <- list(accession = res$accession, platforms = res$platforms)
        plat_names <- names(res$platforms)
        output[[paste0(gbid, "_geo_platform_ui")]] <- renderUI(
          if (length(plat_names) > 1) selectInput(ns(paste0(gbid, "_geo_platform")), "This series spans multiple platforms - pick one", choices = plat_names, width = "100%")
          else NULL
        )
      })

      observe({
        gp <- geo_raw_platforms[[gbid]]
        if (is.null(gp)) return()
        chosen <- if (length(gp$platforms) > 1) input[[paste0(gbid, "_geo_platform")]] else names(gp$platforms)[1]
        req(chosen)
        eset <- gp$platforms[[chosen]]
        pm <- multi_geo_platform_matrix(eset, collapse_genes = identical(input[[paste0(gbid, "_type")]] %||% "other", "rnaseq"))
        if (!pm$ok) {
          geo_fetched[[gbid]] <- NULL
          output[[paste0(gbid, "_geo_status")]] <- renderUI(div(class = "empty-note", icon("triangle-exclamation"), pm$error))
          return()
        }
        geo_fetched[[gbid]] <- list(mat = pm$mat, meta = pm$meta, platform = pm$platform, accession = gp$accession, collapsed = pm$collapsed)
        output[[paste0(gbid, "_geo_status")]] <- renderUI(
          div(class = "empty-note", icon("circle-check"),
              sprintf("Fetched %s (platform %s): %s samples x %s features.%s", gp$accession, pm$platform,
                      format(nrow(pm$mat), big.mark = ","), format(ncol(pm$mat), big.mark = ","),
                      if (isTRUE(pm$collapsed)) " Probes collapsed to genes." else "")))
      })
    })

    output$patient_id_col_ui <- renderUI({
      cols <- if (!is.null(raw$meta)) colnames(raw$meta) else character(0)
      if (length(cols) == 0) return(div(class = "empty-note", icon("triangle-exclamation"), "Upload sample metadata first."))
      guess <- cols[grepl("patient", cols, ignore.case = TRUE)]
      selectInput(ns("patient_id_col"), "Patient ID column", choices = cols, selected = if (length(guess) > 0) guess[1] else cols[1])
    })

    observeEvent(input$validate_btn, {
      mode <- input$dataset_source
      gen <- block_reset_gen()
      mats <- list(); validations <- list(); labels <- list(); provenance <- list()

      long_group_dfs <- list()
      explicit_meta_dfs <- list()
      geo_meta_dfs <- list()

      if (identical(mode, "upload")) {
        for (i in seq_len(n_upload_blocks())) {
          ubid <- mo_block_id(i, "upload")
          fi <- input[[mo_file_input_id(ubid, gen)]]
          if (is.null(fi)) next
          label <- input[[paste0(ubid, "_label")]] %||% sprintf("Dataset %d", i)
          otype <- input[[paste0(ubid, "_type")]] %||% "other"
          shape <- input[[paste0(ubid, "_table_shape")]] %||% "wide"

          if (identical(shape, "long")) {
            full_df <- tryCatch(as.data.frame(data.table::fread(fi$datapath, showProgress = FALSE)), error = function(e) NULL)
            if (is.null(full_df)) { showNotification(sprintf("%s: could not read this file as a table.", label), type = "error"); next }
            none_to_na <- function(x) if (is.null(x) || identical(x, "(none)")) NA_character_ else x
            fcol <- none_to_na(input[[paste0(ubid, "_long_feature")]])
            scol <- none_to_na(input[[paste0(ubid, "_long_sample")]])
            vcol <- none_to_na(input[[paste0(ubid, "_long_value")]])
            gcol <- input[[paste0(ubid, "_long_group")]]
            piv <- multi_live_pivot_long(full_df, fcol, scol, vcol, gcol)
            if (!isTRUE(piv$ok)) { showNotification(sprintf("%s: %s", label, piv$error), type = "error"); next }
            mats[[label]] <- piv$mat
            provenance[[label]] <- list(source = "User Upload", detail = sprintf("%s (long/tidy format - %s)", fi$name, piv$note), imported_at = format(Sys.time(), "%d %b %Y %H:%M"))
            if (!is.null(piv$group_df)) {
              if (!is.null(piv$group_df$warning)) showNotification(sprintf("%s: %s", label, piv$group_df$warning), type = "warning")
              if (!is.null(piv$group_df$df)) long_group_dfs[[label]] <- piv$group_df$df
            }
          } else {
            orient <- input[[paste0(ubid, "_orient")]] %||% "samples_rows"
            res <- multi_live_read_matrix(fi$datapath, orientation = orient, filename = fi$name)
            if (!res$ok) { showNotification(sprintf("%s: %s", label, res$error), type = "error"); next }
            if (isTRUE(res$n_coerced_na > 0)) {
              showNotification(sprintf("%s: %d value(s) could not be read as numbers and became missing. If this is a long/tidy table, switch \"Table shape\" above to Long/tidy.", label, res$n_coerced_na), type = "warning")
            }
            mats[[label]] <- res$mat
            provenance[[label]] <- list(source = "User Upload", detail = fi$name, imported_at = format(Sys.time(), "%d %b %Y %H:%M"))
          }
          if (otype %in% c("rnaseq", "methylation")) {
            det <- multi_live_detect_omics_type(mats[[label]])
            mismatch <- identical(det$detected, "unclassifiable") ||
              (identical(otype, "rnaseq") && identical(det$detected, "methylation")) ||
              (identical(otype, "methylation") && identical(det$detected, "rnaseq"))
            if (mismatch) {
              showNotification(sprintf(
                "%s: rejected. %s This dataset was selected as %s, but its structure %s. Supported omics types here are Transcriptomics and DNA Methylomics - if this is one of those, check the feature-ID column and value scale; otherwise pick its correct omics type above.",
                label, det$reason,
                names(MULTI_LIVE_OMICS_TYPES)[match(otype, MULTI_LIVE_OMICS_TYPES)] %||% otype,
                if (identical(det$detected, "unclassifiable")) "could not be confidently classified as either" else sprintf("looks like %s instead", if (identical(det$detected, "methylation")) "DNA Methylomics" else "Transcriptomics")
              ), type = "error", duration = 15)
              mats[[label]] <- NULL
              next
            }
          }
          validations[[label]] <- multi_live_validate_matrix(mats[[label]], layer_label = label)
          labels[[ubid]] <- label
          mfi <- input[[mo_meta_file_input_id(ubid, gen)]]
          m <- mo_read_meta_file(mfi)
          if (!is.null(m)) {
            explicit_meta_dfs[[label]] <- m
          } else if (!is.null(mfi)) {
            showNotification(sprintf("%s: metadata file could not be read as a table and was skipped.", label), type = "warning")
          }
        }
      } else if (identical(mode, "geo")) {
        for (i in seq_len(n_geo_blocks())) {
          gbid <- mo_block_id(i, "geo")
          gf <- geo_fetched[[gbid]]
          if (is.null(gf)) next
          label <- input[[paste0(gbid, "_label")]] %||% sprintf("Dataset %d", i)
          otype <- input[[paste0(gbid, "_type")]] %||% "other"
          if (otype %in% c("rnaseq", "methylation")) {
            det <- multi_live_detect_omics_type(gf$mat)
            mismatch <- identical(det$detected, "unclassifiable") ||
              (identical(otype, "rnaseq") && identical(det$detected, "methylation")) ||
              (identical(otype, "methylation") && identical(det$detected, "rnaseq"))
            if (mismatch) {
              showNotification(sprintf(
                "%s: rejected. %s This dataset was selected as %s, but its structure %s.", label, det$reason,
                names(MULTI_LIVE_OMICS_TYPES)[match(otype, MULTI_LIVE_OMICS_TYPES)] %||% otype,
                if (identical(det$detected, "unclassifiable")) "could not be confidently classified as either Transcriptomics or DNA Methylomics" else sprintf("looks like %s instead", if (identical(det$detected, "methylation")) "DNA Methylomics" else "Transcriptomics")
              ), type = "error", duration = 15)
              next
            }
          }
          v <- multi_live_validate_matrix(gf$mat, layer_label = label)
          mats[[label]] <- gf$mat
          validations[[label]] <- v
          labels[[gbid]] <- label
          provenance[[label]] <- list(source = "NCBI GEO", detail = sprintf("%s (platform %s)", gf$accession, gf$platform), imported_at = format(Sys.time(), "%d %b %Y %H:%M"))
          if (!is.null(gf$meta)) geo_meta_dfs[[label]] <- gf$meta
          mfi <- input[[mo_meta_file_input_id(gbid, gen)]]
          m <- mo_read_meta_file(mfi)
          if (!is.null(m)) {
            explicit_meta_dfs[[label]] <- m
          } else if (!is.null(mfi)) {
            showNotification(sprintf("%s: metadata file could not be read as a table and was skipped.", label), type = "warning")
          }
        }
      }

      if (length(mats) < 2) showNotification("Add at least two datasets to proceed.", type = "warning")

      raw$mats <- mats
      raw$validations <- validations
      raw$labels <- labels
      raw$provenance <- provenance
      proc$filtered_mats <- NULL; proc$scaled_mats <- NULL; proc$batch_corrected <- NULL

      if (length(long_group_dfs) > 0) {
        merged_groups <- Reduce(mo_merge_sample_meta, long_group_dfs)
        raw$meta <- mo_merge_sample_meta(raw$meta, merged_groups)
      }
      if (length(geo_meta_dfs) > 0) {
        merged_geo <- Reduce(mo_merge_sample_meta, geo_meta_dfs)
        raw$meta <- mo_merge_sample_meta(raw$meta, merged_geo)
      }
      if (length(explicit_meta_dfs) > 0) {
        merged_explicit <- Reduce(mo_merge_sample_meta, explicit_meta_dfs)
        raw$meta <- mo_merge_sample_meta(raw$meta, merged_explicit)
      }

      showNotification(sprintf("Validated %d dataset(s).", length(mats)), type = "message")
    })

    output$validate_ui <- renderUI({
      if (length(raw$validations) == 0) return(multi_empty_state(mo_load_first_msg(input$dataset_source)))
      tagList(
        div(style = "display:flex; gap:14px; flex-wrap:wrap; margin-bottom:14px;",
            lapply(names(raw$validations), function(nm) {
              v <- raw$validations[[nm]]
              st <- multi_dataset_status(v)
              mo_dataset_block_card(nm, v$n_samples %||% NA, v$n_features %||% NA, st)
            })
        ),
        DT::dataTableOutput(ns("qc_table")),
        p(class = "submodule-desc", tags$em("Every value above is computed from your actual datasets.")),
        lapply(raw$validations, function(v) {
          if (!isTRUE(v$ok)) return(NULL)
          issues <- c(
            if (v$n_duplicate_samples > 0) sprintf("%d duplicate sample ID(s)", v$n_duplicate_samples),
            if (v$n_duplicate_features > 0) sprintf("%d duplicate feature ID(s)", v$n_duplicate_features),
            if (v$n_zero_variance > 0) sprintf("%d zero-variance feature(s)", v$n_zero_variance),
            if (v$n_non_finite > 0) sprintf("%d non-finite (Inf/NaN) value(s)", v$n_non_finite),
            if (v$pct_missing > 20) sprintf("%.0f%% missing values", v$pct_missing)
          )
          if (length(issues) == 0) return(NULL)
          div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
              sprintf(" %s: %s.", v$layer, paste(issues, collapse = "; ")))
        })
      )
    })
    output$qc_table <- DT::renderDataTable({
      tbl <- req(multi_live_qc_summary_table(raw$validations))
      DT::datatable(tbl, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    overlap <- reactive({
      req(length(raw$mats) >= 2)
      mats <- raw$mats
      if (identical(input$matching_method, "patient_id")) {
        res <- mo_apply_matching(mats, "patient_id", meta = raw$meta, patient_col = input$patient_id_col)
        mats <- res$mats
      } else if (identical(input$matching_method, "mapping")) {
        mf <- input[[paste0("mapping_file_g", block_reset_gen())]]
        mapping_df <- if (!is.null(mf)) tryCatch(as.data.frame(data.table::fread(mf$datapath, showProgress = FALSE)), error = function(e) NULL) else NULL
        if (!is.null(mf) && is.null(mapping_df)) {
          showNotification("Sample mapping file could not be read as a table and was skipped.", type = "warning")
        }
        res <- mo_apply_matching(mats, "mapping", mapping_df = mapping_df)
        mats <- res$mats
      }
      multi_live_sample_overlap(mats)
    })

    output$matching_ui <- renderUI({
      ov <- tryCatch(overlap(), error = function(e) NULL)
      if (is.null(ov) || !isTRUE(ov$ok)) return(multi_empty_state(mo_load_first_msg(input$dataset_source)))
      tagList(
        div(style = "display:flex; gap:10px; flex-wrap:wrap;",
            lapply(names(ov$per_layer), function(nm) div(class = "card", style = "flex:1 1 140px; text-align:center; padding:10px;",
                                                            div(style = sprintf("font-size:1.3em; font-weight:600; color:%s;", ARTHOMIX_COLORS$blue), ov$per_layer[[nm]]),
                                                            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", nm))),
            div(class = "card", style = "flex:1 1 140px; text-align:center; padding:10px;",
                div(style = sprintf("font-size:1.3em; font-weight:600; color:%s;", ARTHOMIX_COLORS$aqua), ov$n_shared),
                div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", "Common across all"))
        ),
        if (ov$n_shared < 3) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
                                   " Fewer than 3 samples shared across selected datasets. Try a different matching method or provide a mapping file.")
        else div(class = "empty-note", icon("circle-check"), sprintf(" Matched-sample integration below will use exactly these %d samples.", ov$n_shared)),
        hr(),
        h5("Missing-data QC"),
        selectInput(ns("miss_layer"), "Dataset", choices = names(raw$mats)),
        multi_plot_or_empty(function() multi_live_missingness_by_omics_plot(raw$validations), ns("miss_by_omics"), height = "260px"),
        fluidRow(
          column(6, sliderInput(ns("max_sample_missing"), "Max sample missingness (%)", min = 0, max = 100, value = 50)),
          column(6, sliderInput(ns("max_feature_missing"), "Max feature missingness (%)", min = 0, max = 100, value = 50))
        ),
        multi_plot_or_empty(sample_miss_plot_fn, ns("sample_miss_plot"), height = "300px"),
        multi_plot_or_empty(feature_miss_plot_fn, ns("feature_miss_plot"), height = "260px"),
        selectInput(ns("impute_method"), "Missing-value handling (applied per dataset before normalization)",
                    choices = c("Leave as-is" = "none", "Mean imputation" = "mean", "Median imputation" = "median",
                                "Remove samples/features exceeding thresholds" = "remove_rows"), selected = "none"),
        p(class = "submodule-desc", tags$strong("Missing-value handling affects downstream results."), " Nothing is imputed until you pick a method and apply it in step 3.")
      )
    })
    output$miss_by_omics <- renderPlot(multi_live_missingness_by_omics_plot(raw$validations))
    sample_miss_plot_fn <- reactive(multi_live_sample_missingness_plot(multi_live_missingness(raw$mats[[req(input$miss_layer)]]), threshold = input$max_sample_missing))
    output$sample_miss_plot <- renderPlot(sample_miss_plot_fn())
    feature_miss_plot_fn <- reactive(multi_live_feature_missingness_plot(multi_live_missingness(raw$mats[[req(input$miss_layer)]])))
    output$feature_miss_plot <- renderPlot(feature_miss_plot_fn())

    output$preprocess_ui <- renderUI({
      if (length(raw$mats) < 2) return(multi_empty_state(mo_load_first_msg(input$dataset_source)))
      tagList(
        p(class = "submodule-desc", "Normalization choices are restricted to what's appropriate for each dataset's omics type."),
        uiOutput(ns("norm_controls")),
        actionButton(ns("preprocess_btn"), "Apply normalization, filtering, and scaling", icon = icon("gears"), class = "btn-primary btn-sm"),
        p(class = "submodule-desc", tags$em("Normalization is never applied silently - choose it based on assay type and analysis goals.")),
        uiOutput(ns("preprocess_result_ui"))
      )
    })

    output$norm_controls <- renderUI({
      req(length(raw$mats) > 0)
      tagList(lapply(names(raw$mats), function(label) {
        otype <- mo_label_omics_type(label, input, n_upload_blocks(), n_geo_blocks(), input$dataset_source)
        box(width = NULL, title = label, status = "primary", solidHeader = FALSE, collapsible = TRUE,
            selectInput(ns(paste0("norm_", make.names(label))), "Normalization", choices = MULTI_LIVE_NORM_CHOICES[[otype]] %||% MULTI_LIVE_NORM_CHOICES$other),
            numericInput(ns(paste0("topvar_", make.names(label))), "Feature filtering: keep top-N most variable features (blank = keep all)", value = NA, min = 10))
      }))
    })

    observeEvent(input$preprocess_btn, {
      ov <- tryCatch(overlap(), error = function(e) NULL)
      validate(need(!is.null(ov) && isTRUE(ov$ok) && ov$n_shared >= 3, "Need at least 3 matched samples across datasets before preprocessing."))
      matched <- ov$shared_ids
      use_mats <- ov$mats %||% raw$mats

      out_mats <- list()
      for (label in names(use_mats)) {
        m <- use_mats[[label]][matched, , drop = FALSE]
        imp <- multi_live_handle_missing(m, method = input$impute_method %||% "none",
                                          max_sample_missing_pct = input$max_sample_missing %||% 100,
                                          max_feature_missing_pct = input$max_feature_missing %||% 100)
        if (!imp$ok) { showNotification(sprintf("%s: %s", label, imp$error), type = "error"); next }
        m <- imp$mat
        otype <- mo_label_omics_type(label, input, n_upload_blocks(), n_geo_blocks(), input$dataset_source)
        norm <- multi_live_normalize(m, otype, input[[paste0("norm_", make.names(label))]] %||% "none")
        m <- norm$mat
        top_n <- input[[paste0("topvar_", make.names(label))]]
        if (!is.null(top_n) && !is.na(top_n) && top_n > 0) {
          filt <- multi_live_filter_features(m, criterion = "variance", keep_top_n = top_n)
          m <- filt$mat
        }
        out_mats[[label]] <- m
      }
      req(length(out_mats) >= 2)
      proc$filtered_mats <- out_mats
      proc$scaled_mats <- lapply(out_mats, multi_live_scale)
      showNotification("Preprocessing applied.", type = "message")
    })

    output$preprocess_result_ui <- renderUI({
      if (is.null(proc$filtered_mats)) return(NULL)
      tagList(
        h5("Feature retention"),
        fluidRow(lapply(names(proc$filtered_mats), function(nm) column(6,
          multi_plot_or_empty(function() multi_live_retention_plot(ncol(raw$mats[[nm]]), ncol(proc$filtered_mats[[nm]])), ns(paste0("retention_", make.names(nm))), height = "220px")
        ))),
        h5("Distribution before / after normalization"),
        selectInput(ns("dist_layer"), "Dataset", choices = names(proc$filtered_mats)),
        radioButtons(ns("dist_kind"), NULL, choices = c("Boxplot" = "box", "Density" = "density"), inline = TRUE),
        fluidRow(
          column(6, p(tags$strong("Before")), multi_plot_or_empty(before_dist_fn, ns("dist_before"), height = "280px")),
          column(6, p(tags$strong("After")), multi_plot_or_empty(after_dist_fn, ns("dist_after"), height = "280px"))
        ),
        h5("Cross-dataset scale comparison (after z-score standardization)"),
        multi_plot_or_empty(scale_compare_fn, ns("scale_compare"), height = "300px")
      )
    })
    before_dist_fn <- reactive(multi_live_distribution_plot(raw$mats[[req(input$dist_layer)]], kind = input$dist_kind %||% "box"))
    output$dist_before <- renderPlot(before_dist_fn())
    after_dist_fn <- reactive(multi_live_distribution_plot(proc$filtered_mats[[req(input$dist_layer)]], kind = input$dist_kind %||% "box"))
    output$dist_after <- renderPlot(after_dist_fn())
    scale_compare_fn <- reactive({
      req(proc$scaled_mats)
      multi_live_scale_comparison_plot(proc$scaled_mats, names(proc$scaled_mats))
    })
    output$scale_compare <- renderPlot(scale_compare_fn())

    observe({
      req(proc$filtered_mats)
      for (nm in names(proc$filtered_mats)) {
        local({
          nm_local <- nm
          output[[paste0("retention_", make.names(nm_local))]] <- renderPlot(multi_live_retention_plot(ncol(raw$mats[[nm_local]]), ncol(proc$filtered_mats[[nm_local]])))
        })
      }
    })

    output$batch_ui <- renderUI({
      if (is.null(proc$scaled_mats)) return(multi_empty_state("Apply preprocessing (step 3) first."))
      meta_cols <- if (!is.null(raw$meta)) colnames(raw$meta) else character(0)
      tagList(
        selectInput(ns("batch_layer"), "Dataset", choices = names(proc$scaled_mats)),
        if (length(meta_cols) == 0) div(class = "empty-note", icon("circle-info"), "No metadata uploaded - batch-effect assessment isn't possible. PCA below still shows structure without a color-by variable.")
        else tagList(
          selectInput(ns("color_by"), "Color PCA by", choices = c("(none)" = "", meta_cols)),
          selectInput(ns("batch_col"), "Batch column", choices = c("(none)" = "", meta_cols)),
          selectInput(ns("phenotype_col"), "Phenotype/group column", choices = c("(none)" = "", meta_cols))
        ),
        h5("Before correction"),
        fluidRow(
          column(6, p(tags$strong("PCA")), multi_plot_or_empty(pca_before_fn, ns("pca_before"), height = "340px")),
          column(6, p(tags$strong("Sample correlation")), multi_plot_or_empty(corr_before_fn, ns("corr_before"), height = "340px"))
        ),
        uiOutput(ns("confound_ui")),
        conditionalPanel(condition = sprintf("input['%s'] != '' && input['%s'] != ''", ns("batch_col"), ns("phenotype_col")), tagList(
          selectInput(ns("correct_method"), "Correction method", choices = c("ComBat (empirical Bayes)" = "combat", "limma::removeBatchEffect" = "limma")),
          uiOutput(ns("confound_override_ui")),
          actionButton(ns("correct_btn"), "Apply batch correction", icon = icon("play"), class = "btn-primary btn-sm"),
          uiOutput(ns("batch_success_panel")),
          uiOutput(ns("batch_after_ui"))
        ))
      )
    })
    pca_before_fn <- reactive(multi_live_pca_plot(multi_live_pca(proc$scaled_mats[[req(input$batch_layer)]]), raw$meta, if (nzchar(input$color_by %||% "")) input$color_by else NULL))
    output$pca_before <- renderPlot(pca_before_fn())
    corr_before_fn <- reactive({
      d <- multi_live_sample_correlation_data(proc$scaled_mats[[req(input$batch_layer)]])
      req(isTRUE(d$ok))
      multi_live_correlation_heatmap_plot(d$df)
    })
    output$corr_before <- renderPlot(corr_before_fn())

    output$confound_ui <- renderUI({
      req(input$batch_col, input$phenotype_col, nzchar(input$batch_col), nzchar(input$phenotype_col), raw$meta)
      cc <- multi_live_confounding_check(raw$meta, input$batch_col, input$phenotype_col)
      if (is.null(cc)) return(NULL)
      if (isTRUE(cc$confounded)) {
        div(class = "empty-note", style = "border-color: var(--color-danger, #d9534f);", icon("triangle-exclamation"),
            " Potential confounding detected. The batch column and the phenotype column are strongly associated - batch correction may remove genuine biological signal and cannot reliably separate batch from phenotype. Correction is blocked below unless you explicitly override this.")
      } else {
        div(class = "empty-note", icon("circle-check"), sprintf(" No strong batch/phenotype confounding detected (chi-square p = %.3f).", cc$p_value %||% NA))
      }
    })
    confounded_now <- reactive({
      req(input$batch_col, input$phenotype_col, nzchar(input$batch_col), nzchar(input$phenotype_col), raw$meta)
      cc <- multi_live_confounding_check(raw$meta, input$batch_col, input$phenotype_col)
      isTRUE(cc$confounded)
    })
    output$confound_override_ui <- renderUI({
      if (!isTRUE(tryCatch(confounded_now(), error = function(e) FALSE))) return(NULL)
      checkboxInput(ns("confound_override"), "I understand batch and phenotype appear confounded and want to proceed anyway.", value = FALSE)
    })

    observeEvent(input$correct_btn, {
      req(proc$scaled_mats, input$batch_layer, input$batch_col, raw$meta)
      confounded <- isTRUE(tryCatch(confounded_now(), error = function(e) FALSE))
      if (confounded && !isTRUE(input$confound_override)) {
        showNotification("Batch correction is blocked because the batch variable appears confounded with the phenotype variable. Check the box above to proceed anyway.", type = "error")
        return()
      }
      m <- proc$scaled_mats[[input$batch_layer]]
      common <- intersect(rownames(m), rownames(raw$meta))
      if (length(common) < 3) {
        showNotification("Not enough samples with both data and batch metadata.", type = "error")
        return()
      }
      batch <- raw$meta[common, input$batch_col]
      res <- multi_live_batch_correct(m[common, , drop = FALSE], batch, method = input$correct_method %||% "combat")
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      proc$batch_corrected <- res$mat
      proc$batch_summary <- list(
        dimensions = dim(res$mat), batch_col = input$batch_col, phenotype_col = input$phenotype_col,
        method = res$method, n_batches = length(unique(batch)), layer = input$batch_layer,
        timestamp = Sys.time()
      )
      showNotification("Batch correction applied.", type = "message")
    })
    pca_after_fn <- reactive({
      req(proc$batch_corrected)
      multi_live_pca_plot(multi_live_pca(proc$batch_corrected), raw$meta, if (nzchar(input$color_by %||% "")) input$color_by else NULL)
    })
    output$pca_after <- renderPlot(pca_after_fn())
    corr_after_fn <- reactive({
      d <- multi_live_sample_correlation_data(req(proc$batch_corrected))
      req(isTRUE(d$ok))
      multi_live_correlation_heatmap_plot(d$df)
    })
    output$corr_after <- renderPlot(corr_after_fn())

    output$batch_after_ui <- renderUI({
      tagList(
        h5("After correction"),
        fluidRow(
          column(6, p(tags$strong("PCA")), multi_plot_or_empty(pca_after_fn, ns("pca_after"), height = "340px")),
          column(6, p(tags$strong("Sample correlation")), multi_plot_or_empty(corr_after_fn, ns("corr_after"), height = "340px"))
        ),
        uiOutput(ns("variance_diagnostic_ui")),
        if (!is.null(proc$batch_corrected)) div(class = "table-toolbar", downloadButton(ns("dl_batch_corrected_csv"), "Download corrected data (CSV)", class = "btn-sm"))
      )
    })

    output$batch_success_panel <- renderUI({
      s <- proc$batch_summary
      if (is.null(s) || !identical(s$layer, input$batch_layer)) return(NULL)
      div(class = "empty-note", icon("circle-check"),
          tags$strong(" Batch correction completed successfully. "),
          sprintf("Dimensions: %s samples x %s features. Batch variable: \"%s\" (%d level%s). Biological variable preserved: \"%s\". Method: %s.",
                  format(s$dimensions[1], big.mark = ","), format(s$dimensions[2], big.mark = ","),
                  s$batch_col, s$n_batches, if (s$n_batches == 1) "" else "s", s$phenotype_col,
                  if (identical(s$method, "combat")) "ComBat (empirical Bayes)" else "limma::removeBatchEffect"))
    })
    output$dl_batch_corrected_csv <- downloadHandler(
      filename = function() sprintf("%s_batch_corrected.csv", make.names(input$batch_layer %||% "dataset")),
      content = function(file) {
        m <- req(proc$batch_corrected)
        utils::write.csv(m, file, row.names = TRUE)
      }
    )

    output$variance_diagnostic_ui <- renderUI({
      req(proc$batch_corrected, input$batch_col, input$phenotype_col, nzchar(input$batch_col), nzchar(input$phenotype_col))
      pca_b <- multi_live_pca(proc$scaled_mats[[input$batch_layer]])
      pca_a <- multi_live_pca(proc$batch_corrected)
      if (!isTRUE(pca_b$ok) || !isTRUE(pca_a$ok)) return(NULL)
      tagList(
        h5("Quantitative diagnostic (R² of PC1/PC2 vs. group, not a \"looks better\" claim)"),
        DT::dataTableOutput(ns("variance_diag_table")),
        p(class = "submodule-desc", "A good correction reduces R² against batch while keeping R² against phenotype stable. If phenotype R² also collapses, reconsider the method or covariates.")
      )
    })
    output$variance_diag_table <- DT::renderDataTable({
      pca_b <- multi_live_pca(proc$scaled_mats[[req(input$batch_layer)]])
      pca_a <- multi_live_pca(req(proc$batch_corrected))
      df <- data.frame(
        PC = c("PC1", "PC2"),
        `R2 vs batch (before)` = multi_live_variance_by_group(pca_b$scores, raw$meta, input$batch_col),
        `R2 vs batch (after)` = multi_live_variance_by_group(pca_a$scores, raw$meta, input$batch_col),
        `R2 vs phenotype (before)` = multi_live_variance_by_group(pca_b$scores, raw$meta, input$phenotype_col),
        `R2 vs phenotype (after)` = multi_live_variance_by_group(pca_a$scores, raw$meta, input$phenotype_col),
        check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t"), class = "stripe hover compact")
    })

    compat <- reactive({
      req(length(raw$validations) > 0)
      ov <- tryCatch(overlap(), error = function(e) NULL)
      multi_dataset_compatibility(raw$validations, ov, has_metadata = !is.null(raw$meta))
    })

    output$compat_ui <- renderUI({
      cmp <- tryCatch(compat(), error = function(e) NULL)
      if (is.null(cmp) || length(cmp$per_layer) == 0) return(multi_empty_state(mo_load_first_msg(input$dataset_source)))
      overall_color <- if (identical(cmp$overall_label, "READY")) ARTHOMIX_COLORS$aqua
        else if (grepl("READY WITH REVIEW", cmp$overall_label)) ARTHOMIX_COLORS$yellow
        else ARTHOMIX_COLORS$red
      tagList(
        box(width = NULL, title = "Dataset Compatibility", status = "primary", solidHeader = FALSE,
            tags$table(class = "table table-condensed",
                        tags$tbody(
                          lapply(cmp$per_layer, function(p) tags$tr(tags$td(p$label), tags$td(mo_status_badge(p$status)))),
                          tags$tr(tags$td(tags$strong("Sample matching")), tags$td(
                            if (cmp$sample_matching_ok) span(style = sprintf("color:%s;", ARTHOMIX_COLORS$aqua), icon("circle-check"), " OK")
                            else span(style = sprintf("color:%s;", ARTHOMIX_COLORS$red), icon("circle-xmark"), " Insufficient"))),
                          tags$tr(tags$td(tags$strong("Metadata")), tags$td(if (cmp$has_metadata) "Available" else "Not provided"))
                        )),
            div(style = sprintf("margin-top:10px; font-weight:700; color:%s;", overall_color), sprintf("Overall status: %s", cmp$overall_label))
        ),
        box(width = NULL, title = "Datasets Available", status = "primary", solidHeader = FALSE,
            uiOutput(ns("active_checkbox_ui")),
            actionButton(ns("activate_btn"), "Use Selected Datasets for Multi-Omics Analysis", icon = icon("play"), class = "btn-primary btn-sm"),
            uiOutput(ns("activate_message_ui")))
      )
    })

    output$active_checkbox_ui <- renderUI({
      cmp <- tryCatch(compat(), error = function(e) NULL)
      req(cmp)
      eligible <- names(Filter(function(p) p$status$level != "not_compatible", cmp$per_layer))
      checkboxGroupInput(ns("active_layers"), NULL, choices = names(cmp$per_layer), selected = eligible)
    })

    observeEvent(input$activate_btn, {
      chosen <- input$active_layers
      validate(need(length(chosen) >= 2, "Select at least two compatible datasets."))
      cmp <- compat()
      bad <- intersect(chosen, names(Filter(function(p) p$status$level == "not_compatible", cmp$per_layer)))
      validate(need(length(bad) == 0, sprintf("\"%s\" is marked Not Compatible and cannot be used for analysis.", paste(bad, collapse = ", "))))

      validate(need(!is.null(proc$scaled_mats), "Preprocess the selected datasets first (step 3: \"Apply normalization, filtering, and scaling\") before activating them."))
      final_mats <- if (!is.null(proc$batch_corrected) && !is.null(input$batch_layer)) {
        c(proc$scaled_mats[setdiff(names(proc$scaled_mats), input$batch_layer)], setNames(list(proc$batch_corrected), input$batch_layer))
      } else proc$scaled_mats
      final_mats <- final_mats[intersect(names(final_mats), chosen)]
      validate(need(length(final_mats) >= 2, "At least two of the selected datasets must have been preprocessed (step 3) before activating them."))

      ov_now <- tryCatch(overlap(), error = function(e) NULL)
      layer_meta <- lapply(names(final_mats), function(nm) list(
        omics_type = mo_label_omics_type(nm, input, n_upload_blocks(), n_geo_blocks(), input$dataset_source),
        validation = raw$validations[[nm]],
        status = cmp$per_layer[[nm]]$status,
        processing = if (!is.null(proc$filtered_mats)) "Normalized, filtered, scaled" else "Not processed",
        provenance = raw$provenance[[nm]]
      ))
      names(layer_meta) <- names(final_mats)

      multi_dataset$source <- input$dataset_source
      multi_dataset$active <- TRUE
      multi_dataset$layers <- final_mats
      multi_dataset$layer_meta <- layer_meta
      multi_dataset$sample_meta <- raw$meta
      multi_dataset$overlap <- ov_now
      multi_dataset$loaded_at <- Sys.time()

      output$activate_message_ui <- renderUI(
        div(class = "empty-note", icon("circle-check"),
            sprintf(" Active analysis dataset: %s (%d datasets, %s matched samples).", paste(names(final_mats), collapse = " + "), length(final_mats),
                    format(if (!is.null(ov_now)) ov_now$n_shared else NA, big.mark = ","))))
    })

    output$provenance_ui <- renderUI(mo_provenance_ui(multi_dataset$layer_meta %||% list()))

    mod_multi_mofa_server("integrated", multi_dataset, multi_results)
  })
}
