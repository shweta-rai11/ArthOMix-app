## R/transcriptomics/04_Differential_Expression/mod_dge.R - Differential Expression submodule.
## Fits a live limma or DESeq2 contrast between two levels of any metadata
## column, with an optional covariate/filter column; method must match the

mod_dge_config <- list(
  id = "dge", group = "Data",
  title = "Differential Expression",
  description = "Fit a limma or DESeq2 model comparing two levels of a metadata column, optionally stratified by sex (Female/Male).",
  icon = "chart-column"
)

dge_read_expr_upload <- function(datapath, filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("rds")) {
    loaded <- safe_read_rds(datapath)
    validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
    obj <- loaded$value
    m <- if (is.matrix(obj)) obj else as.matrix(obj)
    storage.mode(m) <- "double"
    return(m)
  }
  res <- cx_read_table(datapath, filename)
  validate(need(res$ok, res$error))
  df <- res$df
  validate(need(ncol(df) >= 2, "Expression file needs a feature-ID column plus at least one sample column."))
  rn <- as.character(df[[1]])
  raw <- df[, -1, drop = FALSE]
  m <- suppressWarnings({ mm <- as.matrix(raw); storage.mode(mm) <- "double"; mm })
  rownames(m) <- rn
  colnames(m) <- colnames(raw)
  m
}

dge_read_table_upload <- function(datapath, filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("rds")) {
    loaded <- safe_read_rds(datapath)
    validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
    obj <- loaded$value
    validate(need(is.data.frame(obj), "The uploaded RDS file must contain a data frame."))
    return(as.data.frame(obj))
  }
  res <- cx_read_table(datapath, filename)
  validate(need(res$ok, res$error))
  res$df
}

dge_clean_expr_matrix <- function(m) {
  notes <- character(0)
  m <- as.matrix(m)
  storage.mode(m) <- "double"

  if (anyDuplicated(rownames(m))) {
    n_dup <- sum(duplicated(rownames(m)))
    rn <- rownames(m)
    sum_m <- rowsum(m, group = rn)
    cnt <- as.numeric(table(rn)[rownames(sum_m)])
    m <- sum_m / cnt
    notes <- c(notes, sprintf("%d duplicate feature ID(s) collapsed to their mean expression.", n_dup))
  }

  row_all_na <- rowSums(!is.na(m)) == 0
  if (any(row_all_na)) {
    notes <- c(notes, sprintf("%d feature(s) with no data (all missing) removed.", sum(row_all_na)))
    m <- m[!row_all_na, , drop = FALSE]
  }
  col_all_na <- colSums(!is.na(m)) == 0
  if (any(col_all_na)) {
    notes <- c(notes, sprintf("%d sample(s) with no data (all missing) removed.", sum(col_all_na)))
    m <- m[, !col_all_na, drop = FALSE]
  }

  vars <- apply(m, 1, stats::var, na.rm = TRUE)
  zero_var <- is.na(vars) | vars == 0
  if (any(zero_var)) {
    notes <- c(notes, sprintf("%d zero-variance feature(s) removed (identical value in every sample).", sum(zero_var)))
    m <- m[!zero_var, , drop = FALSE]
  }

  list(mat = m, notes = notes)
}

mod_dge_ui <- function(id) {
  ns <- NS(id)
  tagList(
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = "Data source", status = "primary", solidHeader = FALSE,
            radioButtons(
              ns("data_source"), "Data source",
              choiceNames = list(
                "Use Dataset Pipeline (the dataset already loaded via the Dataset tab)",
                "Upload your own data"
              ),
              choiceValues = list("pipeline", "upload"), selected = "pipeline"
            ),
            uiOutput(ns("upload_ui"))
          ),
          box(
            width = NULL, title = "Contrast", status = "primary", solidHeader = FALSE,
            uiOutput(ns("method_and_columns_ui")),
            uiOutput(ns("ref_comp_ui")),
            uiOutput(ns("covariate_controls_ui")),
            numericInput(ns("padj_cut"), "Adjusted p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput(ns("lfc_cut"), "Absolute log2 fold-change cutoff", value = 0.1, min = 0, step = 0.1),
            actionButton(ns("run_btn"), "Run differential expression", icon = icon("play"), class = "btn-primary btn-sm"),
            div(style = "margin-top:10px;", uiOutput(ns("saved_runs_ui")))
          )
        ),
        column(
          8,
          conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_btn")),
            box(
              width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
              withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6),
              withSpinner(plotOutput(ns("volcano"), height = 460), color = "#2563EB", type = 6),
              div(class = "table-toolbar", downloadButton(ns("download_volcano_png"), "Download volcano plot (PNG, 7×6in @ 300dpi)", class = "btn-sm"))
            )
          )
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("run_btn")),
        box(
          width = 12, title = "Heatmap of top differentially expressed genes", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Per-gene z-scored expression, clustered by gene and sample, for the most significant genes in this contrast."),
          fluidRow(
            column(4, numericInput(ns("heatmap_n"), "Top significant genes to show", value = 30, min = 2, max = 200, step = 5)),
            column(8, div(style = "padding-top: 25px;", downloadButton(ns("download_heatmap_png"), "Download heatmap (PNG, 300dpi)", class = "btn-sm")))
          ),
          withSpinner(plotOutput(ns("heatmap"), height = 520), color = "#2563EB", type = 6)
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("run_btn")),
        box(
          width = 12, title = "Result table", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Every tested gene, with an up/down/not-significant “direction” column. Use the column filter below the header to show only Up or only Down genes."),
          div(class = "table-toolbar",
              downloadButton(ns("download_dge"), "Download CSV", class = "btn-sm"),
              downloadButton(ns("download_provenance"), "Download analysis record (.json)", class = "btn-sm btn-default")),
          DT::dataTableOutput(ns("dge_table"))
        )
      )
  )
}

mod_dge_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    meta_upload_raw <- reactive({
      req(input$dge_meta_file)
      dge_read_table_upload(input$dge_meta_file$datapath, input$dge_meta_file$name)
    })
    annot_upload_raw <- reactive({
      req(input$dge_annot_file)
      dge_read_table_upload(input$dge_annot_file$datapath, input$dge_annot_file$name)
    })

    upload_type_check <- reactive({
      req(input$dge_expr_file)
      expr <- tryCatch(dge_read_expr_upload(input$dge_expr_file$datapath, input$dge_expr_file$name), error = function(e) NULL)
      req(expr)
      tx_validate_expr_upload(expr, input$dge_declared_data_type)
    })

    output$upload_type_warning_ui <- renderUI({
      checked <- tryCatch(upload_type_check(), error = function(e) NULL)
      req(checked)
      if (!isTRUE(checked$ok)) {
        div(class = "empty-note", icon("triangle-exclamation"), checked$error)
      } else if (!is.null(checked$note)) {
        div(class = "empty-note", icon("triangle-exclamation"), checked$note)
      } else {
        NULL
      }
    })

    output$upload_ui <- renderUI({
      req(input$data_source)
      if (!identical(input$data_source, "upload")) return(NULL)
      tagList(
        tags$hr(),
        div(class = "upload-step-label", "STEP 1 · Upload your files"),
        p(strong("Expression matrix"), " - CSV, TSV, TXT, XLSX, or RDS. Features in rows, samples in columns; for delimited/XLSX files, the first column is the feature ID."),
        radioButtons(ns("dge_declared_data_type"), "Data type", inline = TRUE,
                     choices = c("Raw counts" = "raw", "Normalized (TPM/FPKM/CPM)" = "normalized", "Already log-transformed" = "logtransformed"),
                     selected = "normalized"),
        fileInput(ns("dge_expr_file"), "Expression matrix", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".Rds")),
        uiOutput(ns("upload_type_warning_ui")),
        p(strong("Sample metadata"), " - CSV, TSV, TXT, XLSX, or RDS; one row per sample."),
        fileInput(ns("dge_meta_file"), "Sample metadata", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".Rds")),
        p(strong("Feature annotation"), " (optional) - only needed if the expression matrix's row IDs aren't already gene symbols (e.g. probe IDs)."),
        fileInput(ns("dge_annot_file"), "Feature annotation (optional)", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".Rds")),
        uiOutput(ns("column_mapping_ui"))
      )
    })

    output$column_mapping_ui <- renderUI({
      req(input$data_source == "upload")
      out <- tagList()

      meta_df <- tryCatch(meta_upload_raw(), error = function(e) NULL)
      if (!is.null(meta_df)) {
        cols <- colnames(meta_df)
        out <- tagAppendChildren(
          out,
          div(class = "upload-step-label", "STEP 2 · Map metadata columns"),
          selectInput(ns("map_sample_id"), "Sample ID column (must match expression column names)", choices = cols,
                      selected = guess_column_by_name(cols, c("sample", "sample_id", "id", "geo_accession", "accession")),
                      selectize = FALSE)
        )
      }

      annot_df <- tryCatch(annot_upload_raw(), error = function(e) NULL)
      if (!is.null(annot_df)) {
        acols <- colnames(annot_df)
        out <- tagAppendChildren(
          out,
          div(class = "upload-step-label", "STEP 3 · Map annotation columns"),
          selectInput(ns("map_annot_feature_id"), "Feature/probe ID column (must match expression row IDs)", choices = acols,
                      selected = guess_column_by_name(acols, c("probe", "probe_id", "feature", "feature_id", "id")),
                      selectize = FALSE),
          selectInput(ns("map_annot_symbol"), "Gene symbol column", choices = acols,
                      selected = guess_column_by_name(acols, c("gene_symbol", "symbol", "gene", "genesymbol")),
                      selectize = FALSE)
        )
      }
      out
    })

    active_upload_input <- reactive({
      req(input$dge_expr_file, input$dge_meta_file, input$map_sample_id)
      expr <- dge_read_expr_upload(input$dge_expr_file$datapath, input$dge_expr_file$name)
      meta <- meta_upload_raw()
      meta$sample <- as.character(meta[[input$map_sample_id]])

      if (!is.null(input$dge_annot_file) && !is.null(input$map_annot_feature_id) && !is.null(input$map_annot_symbol)) {
        annot <- annot_upload_raw()
        fid <- as.character(annot[[input$map_annot_feature_id]])
        sym <- as.character(annot[[input$map_annot_symbol]])
        ok <- !is.na(fid) & nzchar(fid) & !is.na(sym) & nzchar(sym)
        fid <- fid[ok]; sym <- sym[ok]
        map <- setNames(sym[!duplicated(fid)], fid[!duplicated(fid)])
        matched <- rownames(expr) %in% names(map)
        new_rn <- rownames(expr)
        new_rn[matched] <- unname(map[new_rn[matched]])
        rownames(expr) <- new_rn
      }

      checked <- tx_validate_expr_upload(expr, input$dge_declared_data_type)
      validate(need(isTRUE(checked$ok), checked$error))
      expr <- checked$mat

      cleaned <- dge_clean_expr_matrix(expr)
      list(
        expr = cleaned$mat, meta = meta,
        source_label = sprintf(
          "Uploaded dataset: %s + %s%s", input$dge_expr_file$name, input$dge_meta_file$name,
          if (!is.null(input$dge_annot_file)) paste0(" + ", input$dge_annot_file$name) else ""
        ),
        mode = "upload", clean_notes = c(cleaned$notes, checked$note),
        declared_type = input$dge_declared_data_type %||% NA_character_
      )
    })

    cur_source <- reactive({
      if (identical(input$data_source, "upload")) {
        up <- tryCatch(active_upload_input(), error = function(e) NULL)
        if (!is.null(up)) return(up)
      }
      list(expr = dataset$expr, meta = dataset$meta, source_label = dataset$source,
           mode = "pipeline", clean_notes = character(0),
           declared_type = dataset$declared_data_type %||% NA_character_)
    })

    candidate_columns <- reactive({
      meta <- cur_source()$meta
      req(meta)
      cols <- setdiff(colnames(meta), "sample")
      keep <- vapply(cols, function(cl) {
        nu <- length(unique(stats::na.omit(meta[[cl]])))
        nu >= 2 && nu <= 20
      }, logical(1))
      cols[keep]
    })

    output$method_and_columns_ui <- renderUI({
      cols <- candidate_columns()
      validate(need(length(cols) > 0, "The loaded metadata has no column with 2-20 distinct values to contrast on."))
      default_contrast <- if ("group" %in% cols) "group" else cols[1]
      tagList(
        radioButtons(
          ns("method"), "Method",
          choiceNames = list(
            "limma (microarray, log-scale, or already-normalised data)",
            "DESeq2 (RNA-seq raw counts)"
          ),
          choiceValues = list("limma", "deseq2"), selected = "limma"
        ),
        selectInput(ns("contrast_col"), "Contrast column", choices = cols, selected = default_contrast, selectize = FALSE),
        selectInput(ns("covariate_col"), "Second column (optional)", choices = c("(none)", cols), selected = "(none)", selectize = FALSE)
      )
    })

    output$ref_comp_ui <- renderUI({
      req(input$contrast_col)
      meta <- cur_source()$meta
      req(input$contrast_col %in% colnames(meta))
      lvls <- sort(unique(na.omit(as.character(meta[[input$contrast_col]]))))
      validate(need(length(lvls) >= 2, "This column has fewer than two distinct values."))
      tagList(
        selectInput(ns("ref_group"), "Reference level", choices = lvls, selected = lvls[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison level", choices = lvls, selected = lvls[min(2, length(lvls))], selectize = FALSE)
      )
    })

    output$covariate_controls_ui <- renderUI({
      req(input$covariate_col)
      if (identical(input$covariate_col, "(none)")) return(NULL)
      tagList(
        radioButtons(
          ns("covariate_mode"), NULL,
          choiceNames = list("Filter to one level of this column", "Adjust for this column as a covariate"),
          choiceValues = list("filter", "adjust"), selected = "filter"
        ),
        uiOutput(ns("covariate_level_ui"))
      )
    })

    output$covariate_level_ui <- renderUI({
      req(input$covariate_col, input$covariate_mode)
      if (identical(input$covariate_col, "(none)") || !identical(input$covariate_mode, "filter")) return(NULL)
      meta <- cur_source()$meta
      req(input$covariate_col %in% colnames(meta))
      lvls <- sort(unique(na.omit(as.character(meta[[input$covariate_col]]))))
      selectInput(ns("covariate_level"), "Restrict to", choices = lvls, selected = lvls[1], selectize = FALSE)
    })

    compute_dge_fit <- function(contrast_col, ref_group, comp_group, covariate_col, covariate_mode, covariate_level, method) {
      cs <- cur_source()
      meta <- cs$meta

      ## Case-insensitive lookup against real metadata column names/values, for the ArthOChat tool-invocation path.
      .ci_resolve <- function(value, choices) {
        if (is.null(value) || !length(choices)) return(value)
        hit <- choices[tolower(choices) == tolower(value)]
        if (length(hit)) hit[[1]] else value
      }
      contrast_col <- .ci_resolve(contrast_col, colnames(meta))
      validate(need(!is.null(contrast_col) && contrast_col %in% colnames(meta),
                    sprintf("Column \"%s\" was not found in the sample metadata.", contrast_col %||% "")))
      actual_levels <- unique(as.character(stats::na.omit(meta[[contrast_col]])))
      ref_group <- .ci_resolve(ref_group, actual_levels)
      comp_group <- .ci_resolve(comp_group, actual_levels)

      validate(need(ref_group != comp_group, "Reference and comparison level must be different."))

      meta <- meta[!is.na(meta[[contrast_col]]) & as.character(meta[[contrast_col]]) %in% c(ref_group, comp_group), , drop = FALSE]

      use_covariate <- !is.null(covariate_col) && !identical(covariate_col, "(none)")
      covariate_label <- ""
      if (use_covariate) {
        meta <- meta[!is.na(meta[[covariate_col]]), , drop = FALSE]
        if (identical(covariate_mode, "filter")) {
          req(covariate_level)
          meta <- meta[as.character(meta[[covariate_col]]) == covariate_level, , drop = FALSE]
          covariate_label <- sprintf(" (%s = %s)", covariate_col, covariate_level)
        }
      }

      common <- intersect(colnames(cs$expr), meta$sample)
      validate(need(length(common) >= 6, "Fewer than 6 samples match this contrast (and filter, if set); pick a different combination."))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      expr <- cs$expr[, common, drop = FALSE]

      grp <- factor(as.character(meta[[contrast_col]]), levels = c(ref_group, comp_group))
      validate(need(all(table(grp) >= 3), "Each level needs at least 3 samples in this contrast to fit a model."))

      adjust_for_covariate <- use_covariate && identical(covariate_mode, "adjust")
      covar <- NULL
      if (adjust_for_covariate) {
        covar <- factor(as.character(meta[[covariate_col]]))
        validate(need(length(unique(covar)) >= 2, "The covariate column has only one level left in this contrast. Pick \"Filter\" instead, or a different column."))
        covariate_label <- sprintf(" adjusted for %s", covariate_col)
      }

      used_method <- method
      declared_type <- cs$declared_type
      if (!is.null(declared_type) && !is.na(declared_type) && nzchar(declared_type)) {
        is_counts <- identical(declared_type, "raw")
        is_normalized_totals <- identical(declared_type, "normalized")
      } else {
        is_counts <- looks_like_raw_counts(expr)
        is_normalized_totals <- looks_like_normalized_totals(expr)
      }
      if (identical(used_method, "deseq2")) {
        validate(need(is_counts, "DESeq2 needs raw, non-negative integer counts, but this data has negative or non-integer values (it looks already normalised/log-scale). Pick limma instead, or load raw counts directly via Dataset → Upload your own data (Preprocessing → Batch Correction always outputs normalised, log-scale data, even with Preprocessing's own log2 set to \"Skip\")."))
        validate(need(!is_normalized_totals, "This data's per-sample totals are tightly pinned near a fixed value (e.g. ~1e6) - the signature of TPM/FPKM/CPM-normalised expression, not raw sequencing counts. DESeq2 requires raw counts; pick limma instead, or load a raw count matrix via Dataset → Upload your own data."))
      } else if (identical(used_method, "limma")) {
        validate(need(!(is_counts && !is_normalized_totals),
          "This data looks like raw, non-negative sequencing counts (wide value range, not library-size-normalised) - limma assumes continuous, roughly-normal data and can give misleading results on raw counts. Pick DESeq2 instead, or load/normalise this to continuous, log-scale data first."))
      }

      tt <- if (identical(used_method, "deseq2")) {
        counts <- round(as.matrix(expr))
        storage.mode(counts) <- "integer"
        col_data <- data.frame(grp = grp, row.names = colnames(counts))
        design <- ~grp
        if (adjust_for_covariate) { col_data$covar <- covar; design <- ~covar + grp }
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = design)
        dds <- tryCatch(suppressMessages(DESeq2::DESeq(dds)),
                         error = function(e) validate(need(FALSE, paste("DESeq2 could not fit this model:", conditionMessage(e)))))
        res <- DESeq2::results(dds, contrast = c("grp", levels(grp)[2], levels(grp)[1]))
        out <- as.data.frame(res)
        out$gene <- rownames(out)
        rownames(out) <- NULL
        out <- dplyr::rename(out, logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj)
        out <- out[!is.na(out$P.Value), , drop = FALSE]
        out[order(out$P.Value), c("gene", setdiff(colnames(out), "gene"))]
      } else {
        design <- tryCatch(
          if (adjust_for_covariate) model.matrix(~0 + grp + covar) else model.matrix(~0 + grp),
          error = function(e) validate(need(FALSE, paste("Could not build a design matrix for this contrast/covariate combination:", conditionMessage(e))))
        )
        safe_levels <- make.names(levels(grp), unique = TRUE)
        colnames(design)[seq_len(nlevels(grp))] <- safe_levels
        aw <- limma::arrayWeights(expr, design)
        fit <- limma::lmFit(expr, design, weights = aw)
        cm <- limma::makeContrasts(contrasts = paste0(safe_levels[2], "-", safe_levels[1]), levels = design)
        fit2 <- tryCatch(limma::eBayes(limma::contrasts.fit(fit, cm)),
                          error = function(e) validate(need(FALSE, paste("limma could not fit this contrast:", conditionMessage(e)))))
        out <- limma::topTable(fit2, number = Inf, sort.by = "P")
        out$gene <- rownames(out)
        rownames(out) <- NULL
        out[, c("gene", setdiff(colnames(out), "gene"))]
      }

      design_formula <- if (identical(used_method, "deseq2")) {
        if (adjust_for_covariate) sprintf("~ %s + %s", covariate_col, contrast_col) else sprintf("~ %s", contrast_col)
      } else {
        if (adjust_for_covariate) sprintf("~0 + %s + %s", contrast_col, covariate_col) else sprintf("~0 + %s", contrast_col)
      }
      test_label <- switch(used_method,
        limma = "moderated t-test (limma eBayes, array-quality-weighted)",
        deseq2 = "Wald test (DESeq2)",
        used_method
      )

      list(
        table = tt, method = used_method,
        expr = expr, grp = grp,
        design_formula = design_formula, test_label = test_label,
        n_ref = sum(grp == levels(grp)[1]), n_comp = sum(grp == levels(grp)[2]),
        label = paste0(comp_group, " vs ", ref_group, " (", contrast_col, ")", covariate_label)
      )
    }

    fit_result <- eventReactive(input$run_btn, {
      req(isTRUE(input$run_btn > 0))
      req(input$contrast_col, input$ref_group, input$comp_group)
      compute_dge_fit(
        input$contrast_col, input$ref_group, input$comp_group,
        input$covariate_col %||% "(none)", input$covariate_mode %||% "filter", input$covariate_level,
        input$method
      )
    })

    dge_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, dge_has_run(TRUE), ignoreInit = TRUE)

    observeEvent(cur_source(), {
      dge_has_run(FALSE)
    }, ignoreInit = TRUE)

    compute_dge_significance <- function(fit, padj_cut, lfc_cut) {
      fit$table %>%
        mutate(
          significant = !is.na(adj.P.Val) & adj.P.Val < padj_cut & abs(logFC) > lfc_cut,
          direction = case_when(
            !significant ~ "Not significant",
            logFC > 0 ~ "Up",
            TRUE ~ "Down"
          )
        )
    }

    sig_table <- reactive({
      req(fit_result())
      compute_dge_significance(fit_result(), input$padj_cut, input$lfc_cut)
    })

    save_dge_result <- function(res, df) {
      out <- list(
        contrast = res$label,
        method = res$method,
        n_samples = res$n_ref + res$n_comp,
        n_tested = nrow(df),
        n_significant = sum(df$significant),
        n_up = sum(df$direction == "Up"),
        n_down = sum(df$direction == "Down"),
        top_hits = head(df$gene[order(df$adj.P.Val)], 10)
      )
      results$dge <- out

      runs <- results$dge_runs %||% list()
      run_id <- paste0("run", length(runs) + 1L)
      runs[[run_id]] <- list(
        contrast = res$label, method = res$method,
        n_samples = res$n_ref + res$n_comp,
        table = df[, c("gene", "logFC", "adj.P.Val", "direction")]
      )
      if (length(runs) > 8) runs <- utils::tail(runs, 8)
      results$dge_runs <- runs

      showNotification(
        sprintf("Saved \"%s\" as run %d of %d this session - pick it up in Candidate Gene Identification's DEG contrast pickers.",
                res$label, length(runs), length(runs)),
        type = "message", duration = 6
      )
      out
    }

    observeEvent(input$run_btn, {
      df <- tryCatch(sig_table(), error = function(e) NULL)
      req(df)
      save_dge_result(fit_result(), df)
    })

    run_dge_now <- function(contrast_col, ref_group, comp_group, method,
                             covariate_col = "(none)", covariate_mode = "filter", covariate_level = NULL,
                             padj_cut = 0.05, lfc_cut = 0.1) {
      fit <- compute_dge_fit(contrast_col, ref_group, comp_group, covariate_col, covariate_mode, covariate_level, method)
      df <- compute_dge_significance(fit, padj_cut, lfc_cut)
      save_dge_result(fit, df)
    }

    latest_run_matching <- function(pattern, negate = FALSE) {
      runs <- results$dge_runs %||% list()
      if (length(runs) == 0) return(NULL)
      labels <- vapply(runs, function(r) r$contrast, character(1))
      matches <- grepl(pattern, labels, ignore.case = TRUE, perl = TRUE)
      hit <- which(if (negate) !matches else matches)
      if (length(hit) == 0) return(NULL)
      labels[[utils::tail(hit, 1)]]
    }

    output$saved_runs_ui <- renderUI({
      sex_pattern <- "\\b(female|male)\\b|\\bF\\b|\\bM\\b"
      all_label <- latest_run_matching(sex_pattern, negate = TRUE)
      female_label <- latest_run_matching("\\bfemale\\b|\\bF\\b")
      male_label <- latest_run_matching("\\bmale\\b|\\bM\\b")

      status_row <- function(sex, label) {
        if (is.null(label)) {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s DEG - not run yet", sex))
        } else {
          tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s DEG completed: ", sex)), label)
        }
      }

      tagList(
        p(class = "submodule-desc", style = "margin-bottom: 4px;",
          "Run status (Female/Male are used by Candidate Gene Identification):"),
        tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
                status_row("All samples", all_label), status_row("Female", female_label), status_row("Male", male_label))
      )
    })

    output$summary_ui <- renderUI({
      if (!dge_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Set a contrast on the left, then click \"Run differential expression\"."))
      }
      res <- fit_result()
      df <- sig_table()
      n_sig <- sum(df$significant)
      n_up <- sum(df$direction == "Up")
      n_down <- sum(df$direction == "Down")
      method_label <- switch(res$method, limma = "limma", deseq2 = "DESeq2", res$method)
      tagList(
        p(strong("Contrast: "), res$label, sprintf(" (n = %d vs %d, %s)", res$n_comp, res$n_ref, method_label)),
        p(strong("Model: "), tags$code(res$design_formula), sprintf(", %s.", res$test_label)),
        p(strong(format(nrow(df), big.mark = ",")), " genes tested, ",
          strong(format(n_sig, big.mark = ",")), " significant at the current cutoffs."),
        div(class = "pipeline-status-strip",
            span(class = "badge-up", sprintf("%s upregulated", format(n_up, big.mark = ","))),
            span(class = "badge-down", sprintf("%s downregulated", format(n_down, big.mark = ","))))
      )
    })
    outputOptions(output, "summary_ui", suspendWhenHidden = FALSE)

    volcano_plot_obj <- reactive({
      df <- sig_table()
      top_labels <- df %>% dplyr::filter(significant) %>% dplyr::arrange(adj.P.Val) %>% head(15)
      n_up <- sum(df$direction == "Up")
      n_down <- sum(df$direction == "Down")
      sig_p <- suppressWarnings(max(df$P.Value[df$significant], na.rm = TRUE))
      xr <- max(c(abs(df$logFC), input$lfc_cut), na.rm = TRUE) * 1.08

      p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = direction, size = direction)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(
          values = c(Up = ARTHOMIX_COLORS$red, Down = ARTHOMIX_COLORS$blue, `Not significant` = "#B8BEC8"),
          name = NULL
        ) +
        scale_size_manual(values = c(Up = 1.7, Down = 1.7, `Not significant` = 1.0), guide = "none") +
        geom_vline(xintercept = c(-input$lfc_cut, input$lfc_cut), linetype = "dashed", linewidth = 0.4, color = "#8A929C") +
        scale_x_continuous(limits = c(-xr, xr)) +
        labs(x = expression(log[2] ~ "fold-change"), y = expression(-log[10] ~ "(" * italic(P) * "-value)")) +
        theme_classic(base_size = 14) +
        theme(
          axis.title = element_text(face = "bold", color = "#1A1A1A"),
          axis.text = element_text(color = "#1A1A1A"),
          axis.line = element_line(color = "#1A1A1A", linewidth = 0.4),
          axis.ticks = element_line(color = "#1A1A1A", linewidth = 0.4),
          legend.position = "top", legend.text = element_text(size = 11),
          plot.margin = ggplot2::margin(16, 18, 10, 10)
        )
      if (is.finite(sig_p)) {
        p <- p + geom_hline(yintercept = -log10(sig_p), linetype = "dashed", linewidth = 0.4, color = "#8A929C")
      }
      p +
        ggrepel::geom_text_repel(data = top_labels, mapping = aes(x = logFC, y = -log10(P.Value), label = gene),
                                  inherit.aes = FALSE, size = 3.2, fontface = "italic", color = "#1A1A1A",
                                  segment.size = 0.3, max.overlaps = 20, show.legend = FALSE) +
        annotate("text", x = -xr, y = Inf, label = sprintf("%s down", format(n_down, big.mark = ",")),
                 hjust = 0, vjust = 1.6, size = 4.2, fontface = "bold", color = ARTHOMIX_COLORS$blue) +
        annotate("text", x = xr, y = Inf, label = sprintf("%s up", format(n_up, big.mark = ",")),
                 hjust = 1, vjust = 1.6, size = 4.2, fontface = "bold", color = ARTHOMIX_COLORS$red)
    })

    output$volcano <- renderPlot({
      if (!dge_has_run()) return(NULL)
      p <- tryCatch(volcano_plot_obj(), error = function(e) NULL)
      req(p)
      p
    })
    outputOptions(output, "volcano", suspendWhenHidden = FALSE)

    output$download_volcano_png <- downloadHandler(
      filename = function() "volcano_plot.png",
      content = function(file) {
        if (!dge_has_run()) stop("No differential expression run yet in this session - click \"Run differential expression\" first.")
        ggsave(file, plot = volcano_plot_obj(), width = 7, height = 6, dpi = 300, bg = "white")
      }
    )

    top_de_genes <- reactive({
      df <- sig_table()
      n <- input$heatmap_n %||% 30
      df %>% dplyr::filter(significant) %>% dplyr::arrange(adj.P.Val) %>% head(n) %>% dplyr::pull(gene)
    })

    heatmap_matrix <- reactive({
      res <- fit_result()
      genes <- top_de_genes()
      validate(need(length(genes) >= 2, "Fewer than 2 significant genes at the current cutoffs to draw a heatmap. Lower the cutoffs, or pick a contrast with more signal."))
      m <- res$expr[genes, , drop = FALSE]
      if (identical(res$method, "deseq2")) m <- log2(m + 1)
      m
    })

    heatmap_args <- reactive({
      res <- fit_result()
      mat <- heatmap_matrix()
      mat_z <- t(scale(t(mat)))
      mat_z[!is.finite(mat_z)] <- 0
      ann_col <- data.frame(Group = res$grp, row.names = colnames(mat_z))
      grp_colors <- setNames(c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$red), levels(res$grp))
      list(
        mat = mat_z, annotation_col = ann_col, annotation_colors = list(Group = grp_colors),
        color = grDevices::colorRampPalette(c(ARTHOMIX_COLORS$blue, "#FFFFFF", ARTHOMIX_COLORS$red))(100),
        show_rownames = nrow(mat_z) <= 60, show_colnames = FALSE,
        clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
        clustering_method = "complete", fontsize = 10, fontsize_row = 7, border_color = NA
      )
    })

    output$heatmap <- renderPlot({
      if (!dge_has_run()) return(NULL)
      a <- tryCatch(heatmap_args(), error = function(e) NULL)
      req(a)
      ph <- pheatmap::pheatmap(
        a$mat, annotation_col = a$annotation_col, annotation_colors = a$annotation_colors,
        color = a$color, show_rownames = a$show_rownames, show_colnames = a$show_colnames,
        clustering_distance_rows = a$clustering_distance_rows, clustering_distance_cols = a$clustering_distance_cols,
        clustering_method = a$clustering_method, fontsize = a$fontsize, fontsize_row = a$fontsize_row,
        border_color = a$border_color, silent = TRUE
      )
      grid::grid.newpage()
      grid::grid.draw(ph$gtable)
    })
    outputOptions(output, "heatmap", suspendWhenHidden = FALSE)

    output$download_heatmap_png <- downloadHandler(
      filename = function() "dge_heatmap.png",
      content = function(file) {
        if (!dge_has_run()) stop("No differential expression run yet in this session - click \"Run differential expression\" first.")
        a <- heatmap_args()
        pheatmap::pheatmap(
          a$mat, annotation_col = a$annotation_col, annotation_colors = a$annotation_colors,
          color = a$color, show_rownames = a$show_rownames, show_colnames = a$show_colnames,
          clustering_distance_rows = a$clustering_distance_rows, clustering_distance_cols = a$clustering_distance_cols,
          clustering_method = a$clustering_method, fontsize = a$fontsize, fontsize_row = a$fontsize_row,
          border_color = a$border_color,
          filename = file, width = 7, height = max(5, 0.15 * nrow(a$mat) + 3), res = 300
        )
      }
    )

    output$dge_table <- DT::renderDataTable({
      req(dge_has_run())
      df <- tryCatch(sig_table(), error = function(e) NULL)
      req(df)
      DT::datatable(df, rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "dge_table", suspendWhenHidden = FALSE)

    output$download_dge <- downloadHandler(
      filename = function() "differential_expression.csv",
      content = function(file) {
        if (!dge_has_run()) stop("No differential expression run yet in this session - click \"Run differential expression\" first.")
        write.csv(sig_table(), file, row.names = FALSE)
      }
    )

    dge_provenance_record <- reactive({
      if (!dge_has_run()) stop("No differential expression run yet in this session - click \"Run differential expression\" first.")
      res <- fit_result()
      declared_type <- cur_source()$declared_type
      arthomix_provenance_record(
        module = "mod_dge",
        checksum_input = list(expr = res$expr, grp = as.character(res$grp)),
        params = list(
          method = res$method,
          data_source = input$data_source,
          contrast_col = input$contrast_col,
          reference_level = input$ref_group,
          comparison_level = input$comp_group,
          covariate_col = input$covariate_col %||% "(none)",
          covariate_mode = if (!identical(input$covariate_col %||% "(none)", "(none)")) input$covariate_mode %||% NA_character_ else NA_character_,
          covariate_level = if (identical(input$covariate_mode %||% "", "filter")) input$covariate_level %||% NA_character_ else NA_character_,
          padj_cut = input$padj_cut, lfc_cut = input$lfc_cut,
          design_formula = res$design_formula,
          n_reference = res$n_ref, n_comparison = res$n_comp
        ),
        seed = NULL,
        packages = if (identical(res$method, "deseq2")) "DESeq2" else "limma",
        extra = if (!is.null(declared_type) && !is.na(declared_type) && nzchar(declared_type))
          list(declared_data_type = declared_type) else list()
      )
    })

    output$download_provenance <- arthomix_provenance_download_handler(dge_provenance_record, "mod_dge_provenance")

    list(run = run_dge_now)
  })
}
