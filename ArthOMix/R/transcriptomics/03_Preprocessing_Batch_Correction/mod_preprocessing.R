## R/transcriptomics/03_Preprocessing_Batch_Correction/mod_preprocessing.R - Preprocessing and Batch Correction (Section 2.2).

mod_preprocessing_config <- list(
  id = "preprocessing", group = "Data",
  title = "Preprocessing and Batch Correction",
  description = "Preprocess and merge one or more datasets based on shared genes and probes, then normalise and batch-correct on the uploaded or pre-loaded data.",
  icon = "filter"
)


## Hover info icon revealing `text` in a floating card (CSS only, see .field-hint* in www/custom.css).
mod_pp_field_hint <- function(text) {
  tags$span(class = "field-hint", tabindex = "0",
            icon("circle-info"),
            tags$span(class = "field-hint-box", text))
}

## Log2-transform radio choices, shared by every preprocessing path in this module -
## each option carries a hover hint explaining when to use it.
pp_log2_choice_names <- function(auto_label = "Auto-detect (recommended)") {
  list(
    span(class = "field-label-with-hint", auto_label,
         mod_pp_field_hint("Log2-transforms the data only if it looks untransformed (e.g. raw intensities or counts). Use this unless you have a reason not to.")),
    span(class = "field-label-with-hint", "Force log2 transform",
         mod_pp_field_hint("Always log2-transforms the data, even if it may already be transformed. Use only if you're sure the data is raw and auto-detect got it wrong.")),
    span(class = "field-label-with-hint", "Skip log2 transform",
         mod_pp_field_hint("Never log2-transforms the data. Use for raw RNA-seq counts that will be normalised downstream (e.g. TMM, DESeq2), or data that's already log-scaled."))
  )
}

pp_guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
  for (term in exact) {
    hit <- cols[tolower(cols) == tolower(term)]
    if (length(hit) > 0) return(hit[1])
  }
  for (term in contains) {
    hit <- cols[grepl(term, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  fallback
}


pp_collapse_probes_to_genes <- function(expr, annot, method = c("median", "maxmean", "mean")) {
  method <- match.arg(method)
  cols <- colnames(annot)
  probe_col <- pp_guess_col(cols, c("probe", "probe_id", "probeid", "probeset", "probesetid", "id"))
  gene_col  <- pp_guess_col(cols, c("gene_symbol", "genesymbol", "symbol", "gene"),
                             fallback = if (length(cols) >= 2) cols[2] else cols[1])
  map <- stats::setNames(as.character(annot[[gene_col]]), as.character(annot[[probe_col]]))
  sym <- unname(map[rownames(expr)])
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym, fixed = TRUE)
  validate(need(any(keep), "None of this expression matrix's row IDs matched the annotation file's probe-ID column, or every match was to more than one gene. Check the annotation file's columns."))
  ex <- expr[keep, , drop = FALSE]; sym <- sym[keep]

  if (method == "maxmean") {
    row_mean <- rowMeans(ex, na.rm = TRUE)
    keep_idx <- tapply(seq_along(sym), sym, function(idx) idx[which.max(row_mean[idx])])
    out <- ex[unlist(keep_idx), , drop = FALSE]
    rownames(out) <- names(keep_idx)
    return(out[order(rownames(out)), , drop = FALSE])
  }
  agg_fn <- if (method == "median") stats::median else base::mean
  apply(ex, 2, function(col_vals) tapply(col_vals, sym, agg_fn, na.rm = TRUE))
}


PP_COHORT_LABELS <- c(
  "GSE93272"  = "Whole Blood Training Cohort A",
  "GSE110169" = "Whole Blood Training Cohort B",
  "GSE15573"  = "PBMC Validation Cohort",
  "GSE89408"  = "Synovial Tissue Validation Cohort"
)
PP_MERGED_COHORT_LABEL <- "Whole Blood Training Cohort (Merged)"

pp_cohort_label <- function(id) {
  if (identical(id, default_dataset_entry$id)) return(PP_MERGED_COHORT_LABEL)
  if (id %in% names(PP_COHORT_LABELS)) unname(PP_COHORT_LABELS[[id]]) else id
}

## Same choices as mod_dataset.R's preloaded_choices(), with this tab's display names.
pp_cohort_choices <- function() {
  ids <- unname(preloaded_choices())
  stats::setNames(ids, vapply(ids, pp_cohort_label, character(1)))
}


pp_preloaded_read <- function(choice_id, log2_choice, dataset = NULL) {
  if (identical(choice_id, "__current__")) {
    ## Prefers the Dataset tab's staged preview, falls back to the active dataset.
    use_expr <- dataset$staged_expr %||% dataset$expr
    use_meta <- dataset$staged_meta %||% dataset$meta
    use_label <- dataset$staged_source %||% dataset$source
    validate(need(!is.null(dataset) && !is.null(use_expr),
                  "No dataset is currently loaded. Preview one on the Dataset tab first."))
    expr <- use_expr
    meta <- use_meta
    label <- use_label %||% "Currently Loaded Dataset"
  } else {
    gse <- choice_id
    if (identical(gse, default_dataset_entry$id)) {
      ## Already merged and batch-corrected - no raw file to read for this one.
      d <- load_default_dataset()
      expr <- d$expr; meta <- d$meta
    } else if (identical(gse, "GSE89408")) {
      d <- load_individual_dataset(gse)
      validate(need(!is.null(d), paste("Raw data for", gse, "was not found on disk.")))
      expr <- d$expr; meta <- d$meta
    } else {
      eset <- get_raw_eset(gse)
      validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
      expr <- get_collapsed_genes(gse)
      meta <- eset_harmonize_meta(eset, gse)
      keep <- !is.na(meta$group)
      meta <- meta[keep, , drop = FALSE]
      expr <- expr[, meta$sample, drop = FALSE]
    }
    label <- pp_cohort_label(gse)
  }
  if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_

  n_samples_before <- ncol(expr); n_genes_before <- nrow(expr)
  q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
  needs_log <- if (identical(log2_choice, "force")) TRUE
               else if (identical(log2_choice, "skip")) FALSE
               else isTRUE(!is.na(q99) && q99 > 100)
  if (needs_log) {
    expr[expr <= 0] <- NA
    expr <- log2(expr)
    expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  }
  ## Median-impute per gene so Batch Correction never sees NAs (same approach as filter_and_transform_expr()).
  if (anyNA(expr)) {
    row_med <- apply(expr, 1, stats::median, na.rm = TRUE)
    na_idx <- which(is.na(expr), arr.ind = TRUE)
    expr[na_idx] <- row_med[na_idx[, 1]]
  }

  list(label = label, expr = as.matrix(expr), meta = meta,
       n_samples_before = n_samples_before, n_samples_after = ncol(expr),
       n_genes_before = n_genes_before, n_genes_after = nrow(expr),
       log2_applied = needs_log)
}

## UI


pp_tab_title <- function(ic, label) {
  tagList(icon(ic), " ", label)
}

pp_study_wise_limma <- function(expr, meta, study_col, ref_group, comp_group, min_n = 3) {
  studies <- sort(unique(stats::na.omit(as.character(meta[[study_col]]))))
  per_study <- list(); skipped <- character(0)
  for (s in studies) {
    idx <- which(as.character(meta[[study_col]]) == s & as.character(meta$group) %in% c(ref_group, comp_group))
    y <- factor(as.character(meta$group[idx]), levels = c(ref_group, comp_group))
    if (length(idx) < 2 * min_n || any(table(y) < min_n)) { skipped <- c(skipped, s); next }
    x <- as.matrix(expr[, meta$sample[idx], drop = FALSE])
    x <- x[stats::complete.cases(x), , drop = FALSE]
    if (nrow(x) < 10) { skipped <- c(skipped, s); next }
    design <- stats::model.matrix(~ y)
    pos <- x[x > 0]
    q99 <- if (length(pos)) stats::quantile(pos, 0.99, na.rm = TRUE) else NA_real_
    looks_like_counts <- is.finite(q99) && q99 > 100 && all(x >= 0, na.rm = TRUE) &&
      mean(abs(x - round(x)) < 1e-6, na.rm = TRUE) > 0.99
    if (looks_like_counts) {
      dge <- edgeR::DGEList(counts = round(x))
      keep <- edgeR::filterByExpr(dge, group = y)
      dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
      v <- limma::voom(dge, design)
      fit <- limma::eBayes(limma::lmFit(v, design))
      method <- "limma-voom (TMM, raw counts)"
    } else {
      if (needs_quantile_norm(summarize_norm_diagnostics(x))) {
        x <- limma::normalizeBetweenArrays(x, method = "quantile")
        method <- "limma (quantile-normalised within study)"
      } else {
        method <- "limma"
      }
      fit <- limma::eBayes(limma::lmFit(x, design))
    }
    se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, 2]
    per_study[[s]] <- data.frame(
      gene = rownames(fit$coefficients), logFC = unname(fit$coefficients[, 2]), se = unname(se),
      p = unname(fit$p.value[, 2]), n_ref = sum(y == ref_group), n_comp = sum(y == comp_group),
      method = method, stringsAsFactors = FALSE
    )
  }
  list(per_study = per_study, skipped = skipped)
}

pp_dl_meta <- function(per_study) {
  genes <- Reduce(union, lapply(per_study, `[[`, "gene"))
  k <- length(per_study)
  Y  <- vapply(per_study, function(d) d$logFC[match(genes, d$gene)], numeric(length(genes)))
  SE <- vapply(per_study, function(d) d$se[match(genes, d$gene)], numeric(length(genes)))
  if (is.null(dim(Y))) { Y <- matrix(Y, ncol = k); SE <- matrix(SE, ncol = k) }
  colnames(Y) <- colnames(SE) <- names(per_study)
  W <- 1 / SE^2
  n_st <- rowSums(!is.na(Y) & is.finite(SE))
  sw <- rowSums(W, na.rm = TRUE)
  theta_fe <- rowSums(W * Y, na.rm = TRUE) / sw
  Q <- rowSums(W * (Y - theta_fe)^2, na.rm = TRUE)
  df <- pmax(n_st - 1, 0)
  C <- sw - rowSums(W^2, na.rm = TRUE) / sw
  tau2 <- (Q - df) / C
  tau2[!is.finite(tau2) | tau2 < 0] <- 0
  Wr <- 1 / (SE^2 + tau2)
  swr <- rowSums(Wr, na.rm = TRUE)
  theta_re <- rowSums(Wr * Y, na.rm = TRUE) / swr
  se_re <- sqrt(1 / swr)
  z <- theta_re / se_re
  p <- 2 * stats::pnorm(-abs(z))
  I2 <- ifelse(Q > 0 & df > 0, pmax(0, (Q - df) / Q) * 100, 0)
  Qp <- ifelse(df > 0, stats::pchisq(Q, df, lower.tail = FALSE), NA_real_)
  res <- data.frame(
    gene = genes, n_studies = n_st,
    logFC_pooled = theta_re, se_pooled = se_re,
    ci_low = theta_re - 1.96 * se_re, ci_high = theta_re + 1.96 * se_re,
    z = z, p_value = p, adj_p = NA_real_,
    logFC_fixed = theta_fe, tau2 = tau2, Q = Q, Q_pvalue = Qp, I2 = I2,
    stringsAsFactors = FALSE
  )
  colnames(Y) <- paste0("logFC_", names(per_study))
  res <- cbind(res, as.data.frame(Y))
  res <- res[res$n_studies >= 2, , drop = FALSE]
  if (nrow(res)) res$adj_p <- stats::p.adjust(res$p_value, method = "BH")
  res[order(res$p_value), , drop = FALSE]
}

mod_preprocessing_ui <- function(id) {
  ns <- NS(id)

  div(
    ## Same nav-tabs styling class as the outer tx_menu tabset (www/custom.css .tx-menu-wrap).
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("tabs"), type = "tabs",
      header = tagList(
        tags$hr()
      ),
      tabPanel(
        value = "Preprocessing", title = pp_tab_title("broom", "Preprocessing"),
        br(), uiOutput(ns("preprocessing_tab_ui"))
      ),
      tabPanel(
        value = "Merge datasets", title = pp_tab_title("code-merge", "Merge Datasets"),
        br(), uiOutput(ns("merge_tab_ui"))
      ),
      tabPanel(
        value = "Batch correction", title = pp_tab_title("wand-magic-sparkles", "Batch Correction"),
        br(), uiOutput(ns("batch_tab_ui"))
      ),
      tabPanel(
        value = "Study-wise meta-analysis", title = pp_tab_title("layer-group", "Study-wise Meta-analysis"),
        br(), uiOutput(ns("meta_tab_ui"))
      ),
      tabPanel(
        value = "Explore", title = tagList(icon("magnifying-glass-chart"), " Data Exploration"),
        br(), mod_data_exploration_ui(ns("eda"))
      )
    )
  )
}

mod_preprocessing_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    
    observeEvent(dataset$source_type, {
      hide <- dataset$source_type %in% c("uploaded", "geo")
      shinyjs::runjs(sprintf(
        "(function(){
           var sel = '#%s';
           var hide = %s;
           function applyVisibility(tabset) {
             var targets = ['Merge datasets', 'Batch correction', 'Explore'];
             var lis = targets.map(function(v){ return tabset.find('a[data-value=\"' + v + '\"]').closest('li'); });
             if (!lis.some(function(li){ return li.length; })) return;
             if (hide) {
               var activeHidden = lis.some(function(li){ return li.hasClass('active'); });
               if (activeHidden) tabset.find('a[data-value=\"Preprocessing\"]').tab('show');
               lis.forEach(function(li){ li.hide(); });
             } else {
               lis.forEach(function(li){ li.show(); });
             }
           }
           var existing = $(sel);
           if (existing.length) { applyVisibility(existing); return; }
           var observer = new MutationObserver(function(){
             var tabset = $(sel);
             if (tabset.length) { observer.disconnect(); applyVisibility(tabset); }
           });
           observer.observe(document.body, { childList: true, subtree: true });
           setTimeout(function(){ observer.disconnect(); }, 600000);
         })();",
        ns("tabs"), if (hide) "true" else "false"
      ))
    }, ignoreNULL = FALSE)

    pp_progress <- reactive({
      merged_ok <- !is.null(tryCatch(merged(), error = function(e) NULL))
      batch_ok  <- !is.null(tryCatch(result(), error = function(e) NULL))
      pl_res <- preloaded_results()
      n_pl <- length(input$preloaded_selected %||% character(0))
      n <- max(1, n_pl)
      n_ready <- sum(vapply(pl_res, function(r) isTRUE(r$ok), logical(1)))
      list(n = n, n_ready = n_ready, merged_ok = merged_ok, batch_ok = batch_ok)
    })


    output$pipeline_summary <- renderUI({
      pr <- pp_progress()
      step_state <- function(done, current) if (done) "done" else if (current) "current" else "future"

      steps <- list(
        list(number = 1, label = "Preprocessing", sublabel = "Filter low-quality samples & genes",
             state = step_state(pr$n_ready == pr$n, pr$n_ready < pr$n)),
        list(number = 2, label = "Merge Datasets", sublabel = "Harmonize and merge cohorts",
             state = step_state(pr$merged_ok, pr$n_ready == pr$n && !pr$merged_ok)),
        list(number = 3, label = "Batch Correction", sublabel = "Remove batch effects (ComBat)",
             state = step_state(pr$batch_ok, pr$merged_ok && !pr$batch_ok))
      )
      pipeline_summary_ui(steps)
    })


    PP_TRAINING_GEO_IDS <- c("GSE93272", "GSE110169")
    PP_TRAINING_COHORT_LABEL <- "Whole Blood Training Cohorts A and B"

    
    available_example_groups <- reactive({
      grps <- unlist(lapply(PP_TRAINING_GEO_IDS, function(gse) {
        eset <- get_raw_eset(gse)
        if (is.null(eset)) return(character(0))
        unique(stats::na.omit(eset_harmonize_meta(eset, gse)$group))
      }))
      sort(unique(grps))
    })

    preloaded_results_val <- reactiveVal(NULL)
    observeEvent(input$preloaded_run, {
      req(length(input$preloaded_selected) > 0)
      choices <- input$preloaded_selected
      out <- lapply(choices, function(choice_id) {
        tryCatch(
          list(ok = TRUE, value = pp_preloaded_read(choice_id, input$preloaded_log2, dataset)),
          error = function(e) list(ok = FALSE,
                                    label = if (identical(choice_id, "__current__")) "Currently Loaded Dataset" else pp_cohort_label(choice_id),
                                    error = conditionMessage(e))
        )
      })
      names(out) <- vapply(seq_along(out), function(i) {
        if (isTRUE(out[[i]]$ok)) out[[i]]$value$label else out[[i]]$label
      }, character(1))
      preloaded_results_val(out)
    }, ignoreInit = TRUE)
    preloaded_results <- function() preloaded_results_val()

    output$preloaded_status_ui <- renderUI({
      if (is.null(input$preloaded_run) || input$preloaded_run == 0) {
        msg <- if (dataset$source_type %in% c("geo", "uploaded")) {
          "Click \"Preprocess\" to run it on the currently loaded dataset."
        } else {
          "Select one or more cohorts above, then click \"Load and Preprocess Selected Cohorts\"."
        }
        return(p(class = "empty-note", icon("circle-info"), msg))
      }
      res <- preloaded_results()
      if (is.null(res)) {
        return(div(class = "empty-note", icon("triangle-exclamation"), "Select at least one cohort above first."))
      }
      tagList(lapply(res, function(r) {
        if (isTRUE(r$ok)) {
          v <- r$value
          div(class = "empty-note", icon("check"),
              sprintf("%s: %s of %s samples kept, %s of %s features kept%s.",
                      v$label, v$n_samples_after, v$n_samples_before,
                      format(v$n_genes_after, big.mark = ","), format(v$n_genes_before, big.mark = ","),
                      if (v$log2_applied) ", log2-transformed" else ""))
        } else {
          div(class = "empty-note", icon("triangle-exclamation"),
              sprintf("%s: %s", r$label, r$error))
        }
      }))
    })

    
    output$activate_current_ui <- renderUI({
      res <- preloaded_results()
      req(length(res) >= 1, isTRUE(res[[1]]$ok))
      tagList(
        actionButton(ns("activate_current_btn"), "Use this preprocessed data as the active dataset",
                     icon = icon("check"), class = "btn-primary btn-sm"),
        uiOutput(ns("activate_current_status_ui"))
      )
    })

    observeEvent(input$activate_current_btn, {
      res <- preloaded_results()
      req(length(res) >= 1, isTRUE(res[[1]]$ok))
      v <- res[[1]]$value
      dataset$expr <- v$expr
      dataset$meta <- v$meta
      dataset$source <- paste0(dataset$source %||% "Currently loaded dataset", " (preprocessed)")
      output$activate_current_status_ui <- renderUI(
        div(class = "empty-note", icon("check"), "This is now the active dataset. Every other sub-module will use it.")
      )
    }, ignoreInit = TRUE)

    
    output$single_dataset_note <- renderUI({
      if (is.null(dataset$expr)) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    if (identical(dataset$source_type, "uploaded"))
                      "Nothing is currently loaded. Upload a dataset on the Dataset tab first."
                    else
                      "Nothing is currently loaded. Fetch a GEO series on the Dataset tab first."))
      }
      div(class = "empty-note", icon("check"),
          sprintf("Using %s: %s genes x %s samples.", dataset$source %||% "Currently Loaded Dataset",
                   format(nrow(dataset$expr), big.mark = ","), ncol(dataset$expr)))
    })

    ## Feeds merge_inputs() below: pick bundled cohorts or the currently-loaded dataset, no upload here.
    output$preprocessing_tab_ui <- renderUI({
      if (dataset$source_type %in% c("geo", "uploaded")) {
        is_upload <- identical(dataset$source_type, "uploaded")
        return(box(
          width = 12,
          title = tagList(icon("broom"), if (is_upload) " Preprocess Uploaded Data" else " Preprocess Fetched GEO Data"),
          status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc",
            if (is_upload)
              "Preprocesses the dataset uploaded on the Dataset tab. A single uploaded dataset is already one dataset, so merging and batch correction don't apply here."
            else
              "Preprocesses the dataset fetched from NCBI GEO on the Dataset tab. A single fetched series is already one dataset, so merging and batch correction don't apply here."),
          uiOutput(ns("single_dataset_note")),
          shinyjs::hidden(checkboxGroupInput(ns("preloaded_selected"), "Cohorts",
                              choices = c("Currently Loaded Dataset" = "__current__"), selected = "__current__")),
          radioButtons(ns("preloaded_log2"), "Log2 Transform",
                       choiceNames = pp_log2_choice_names(), choiceValues = list("auto", "force", "skip"),
                       selected = "auto", inline = TRUE),
          actionButton(ns("preloaded_run"), "Preprocess", icon = icon("play"), class = "btn-primary btn-sm"),
          div(style = "margin-top:8px;", uiOutput(ns("preloaded_status_ui"))),
          uiOutput(ns("activate_current_ui"))
        ))
      }
      box(
        width = 12, title = tagList(icon("database"), " Preloaded Data"), status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc",
          "Select one or more bundled blood or synovium cohorts, or use whatever is currently loaded on the Dataset tab. No upload required - each selection is preprocessed independently."),
        checkboxGroupInput(ns("preloaded_selected"), "Cohorts",
                            choices = c(pp_cohort_choices(), "Currently Loaded Dataset" = "__current__"),
                            inline = TRUE),
        p(class = "empty-note", icon("triangle-exclamation"),
          "These cohorts are raw and platform-specific. Normalization and merging are handled in the next two steps."),
        radioButtons(ns("preloaded_log2"), "Log2 Transform (applied to every selected cohort)",
                     choiceNames = pp_log2_choice_names("Auto-detect per cohort (recommended)"),
                     choiceValues = list("auto", "force", "skip"), selected = "auto", inline = TRUE),
        actionButton(ns("preloaded_run"), "Load and Preprocess Selected Cohorts", icon = icon("play"), class = "btn-primary btn-sm"),
        div(style = "margin-top:8px;", uiOutput(ns("preloaded_status_ui")))
      )
    })

    merge_inputs <- reactive({
      res <- preloaded_results()
      failed <- Filter(Negate(function(r) isTRUE(r$ok)), res %||% list())
      validate(need(length(failed) == 0,
                    sprintf("%s failed to load: %s. Fix and re-run before merging.",
                            paste(vapply(failed, `[[`, character(1), "label"), collapse = ", "),
                            paste(vapply(failed, `[[`, character(1), "error"), collapse = "; "))))
      preloaded_vals <- lapply(Filter(function(r) isTRUE(r$ok), res %||% list()), `[[`, "value")
      all_vals <- preloaded_vals
      validate(need(length(all_vals) > 0,
                    "Load at least one preloaded cohort above before merging."))
      all_vals
    })

    output$merge_tab_ui <- renderUI({
      tagList(
        div(
          class = "card",
          div(class = "card-title", icon("code-merge"), "Merge data for this step"),
          radioButtons(
            ns("merge_mode"), NULL,
            choiceNames = list(
              tagList(
                div(class = "data-source-option-title", "Merge the example pipeline's training datasets"),
                div(class = "data-source-option-desc", "Merges the two raw training datasets on shared genes, before batch correction.")
              ),
              tagList(
                div(class = "data-source-option-title", "Merge your own data"),
                div(class = "data-source-option-desc", "Merges the datasets from the Preprocessing tab.")
              )
            ),
            choiceValues = list("example", "own"), selected = "example"
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'example'", ns("merge_mode")),
          uiOutput(ns("merge_example_ui"))
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'own'", ns("merge_mode")),
          box(
            width = 12, title = tagList(icon("dna"), " Probe-to-gene collapsing (optional)"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc",
              "Turn this on if the loaded datasets are still at probe level (e.g. raw Affymetrix IDs) rather than one row per gene. Applies to every selected dataset using the same annotation file - for different platforms, collapse each dataset separately before uploading."),
            checkboxInput(ns("collapse_probes"), "My selected data is at probe level - collapse to one row per gene before merging", value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s']", ns("collapse_probes")),
              div(class = "field-label-with-hint", span("Probe → gene annotation file"),
                  mod_pp_field_hint("CSV or RDS with at least two columns: the probe/feature ID (matching the expression matrix's row IDs) and the gene symbol. Column names are auto-guessed.")),
              fileInput(ns("collapse_annot_file"), NULL, accept = c(".csv", ".rds", ".Rds")),
              radioButtons(
                ns("collapse_method"), "Collapsing method",
                choiceNames = list(
                  "Median across probes per gene (recommended - matches most published methods, e.g. Zhu et al. 2021)",
                  "Highest-mean probe represents the gene (this app's own preloaded-dataset convention)",
                  "Mean across probes per gene"
                ),
                choiceValues = list("median", "maxmean", "mean"), selected = "median"
              ),
              p(class = "empty-note", icon("circle-info"),
                "Probes with no match in the annotation file, or matched to more than one gene, are dropped - never split or duplicated.")
            )
          ),
          uiOutput(ns("merge_select_ui")),
          uiOutput(ns("merge_venn_ui")),
          box(width = 12, title = tagList(icon("code-merge"), " Merge"), status = "primary", solidHeader = FALSE,
              actionButton(ns("merge_btn"), "Merge datasets", icon = icon("code-merge"), class = "btn-primary btn-sm"),
              withSpinner(uiOutput(ns("merge_summary_ui")), color = "#2563EB", type = 6))
        )
      )
    })

    
    example_default_groups <- reactive({
      groups <- available_example_groups()
      if (length(intersect(c("HC", "RA"), groups)) > 0) intersect(c("HC", "RA"), groups) else groups
    })

    
    output$merge_example_ui <- renderUI({
      groups <- available_example_groups()
      default_groups <- example_default_groups()
      ## Excludes by group identity (not a hardcoded "SLE" label), so this stays accurate if sources change.
      excluded_groups <- setdiff(groups, default_groups)
      excluded_note <- if (length(excluded_groups) > 0) {
        sprintf(" %s also present (%s) - tick above to include for a different comparison.",
                if (length(excluded_groups) == 1) "One other group is" else "Other groups are",
                paste(excluded_groups, collapse = ", "))
      } else ""
      tagList(
        div(
          class = "card",
          if (example_merge_from_raw()) div(class = "empty-note", icon("circle-info"),
              sprintf("Merges %s on shared genes, before batch correction.", PP_TRAINING_COHORT_LABEL)),
          if (length(groups) > 0) tagList(
            div(style = "margin-top:12px;",
                strong("Diagnosis groups to include"),
                p(class = "submodule-desc",
                  sprintf("Defaults to control (HC) + RA only, matching this project's own training cohort.%s", excluded_note)),
                checkboxGroupInput(ns("example_groups"), NULL, choices = groups, selected = default_groups, inline = TRUE))
          ),
          div(style = "margin-top:12px;",
              actionButton(ns("merge_use_example_btn"), "Merge these datasets", icon = icon("code-merge"), class = "btn-primary btn-sm"))
        ),
        uiOutput(ns("merge_venn_example_ui"))
      )
    })

    output$merge_venn_example_ui <- renderUI({
      if (!isTRUE((input$merge_use_example_btn %||% 0) > 0)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Click “Merge these datasets” above to merge and see the result."))
      }
      if (!example_merge_from_raw()) {
        m <- example_live_merge()
        tagList(
          box(width = 12, title = tagList(icon("circle-info"), " Using the merged, batch-corrected cohort"), status = "primary", solidHeader = FALSE,
              p(class = "submodule-desc",
                "Uses the real merged, already-batch-corrected training cohort directly (fast), so Batch Correction below will find little or no residual batch effect left to remove."),
              fluidRow(
                valueBox(format(ncol(m$expr), big.mark = ","), "Samples", icon = icon("users"), color = "light-blue", width = 4),
                valueBox(format(nrow(m$expr), big.mark = ","), "Genes", icon = icon("dna"), color = "green", width = 4),
                valueBox(length(unique(m$meta$dataset)), "Source datasets", icon = icon("layer-group"), color = "purple", width = 4)
              ),
              DT::dataTableOutput(ns("merge_example_composition_table")))
        )
      } else {
        tagList(
          box(width = 12, title = tagList(icon("diagram-project"), " Feature overlap across datasets"), status = "primary", solidHeader = FALSE,
              p(class = "submodule-desc", sprintf(
                "%s, collapsed to gene symbol.",
                PP_TRAINING_COHORT_LABEL
              )),
              fluidRow(
                column(7,
                  withSpinner(plotOutput(ns("venn_plot_example"), height = 440), color = "#2563EB", type = 6),
                  div(class = "table-toolbar", downloadButton(ns("download_venn_example_png"), "Download diagram (PNG)", class = "btn-primary btn-sm"))
                ),
                column(5, DT::dataTableOutput(ns("venn_table_example")))
              )),
          box(width = 12, title = tagList(icon("table-list"), " Region breakdown"), status = "primary", solidHeader = FALSE,
              p(class = "submodule-desc", "Every exact combination of datasets a feature belongs to, largest first. These are the same regions the diagram above draws, shown here as a full table."),
              div(class = "table-toolbar", downloadButton(ns("download_venn_example"), "Download region counts (CSV)", class = "btn-sm")),
              DT::dataTableOutput(ns("venn_region_table_example")))
        )
      }
    })

    output$merge_example_composition_table <- DT::renderDataTable({
      m <- example_live_merge()
      tbl <- as.data.frame(table(Dataset = m$meta$dataset, Group = m$meta$group))
      tbl <- tbl[tbl$Freq > 0, , drop = FALSE]
      colnames(tbl)[colnames(tbl) == "Freq"] <- "Samples"
      DT::datatable(tbl, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    ## Gene-symbol-collapsed feature sets for the two training GEO sources (raw probe IDs don't overlap across platforms).
    example_overlap_sets <- reactive({
      sets <- lapply(PP_TRAINING_GEO_IDS, function(gse) {
        collapsed <- get_collapsed_genes(gse)
        if (is.null(collapsed)) return(NULL)
        rownames(collapsed)
      })
      names(sets) <- PP_TRAINING_GEO_IDS
      sets <- Filter(Negate(is.null), sets)
      validate(need(length(sets) >= 2, "Raw source data for the example cohort's training datasets was not found on disk."))
      sets
    })

    
    example_merge_from_raw <- reactive({
      sel <- input$example_groups %||% example_default_groups()
      !setequal(sel, example_default_groups())
    })

    
    example_live_merge <- reactive({
      if (!example_merge_from_raw()) {
        d <- load_default_dataset()
        meta <- d$meta
        if (!"batch" %in% colnames(meta) || all(is.na(meta$batch))) meta$batch <- meta$dataset
        return(list(expr = d$expr, meta = meta, sources = PP_TRAINING_COHORT_LABEL, n_dup_features = 0L))
      }
      parts <- lapply(PP_TRAINING_GEO_IDS, function(gse) {
        eset <- get_raw_eset(gse)
        validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
        expr <- get_collapsed_genes(gse)
        meta <- eset_harmonize_meta(eset, gse)
        ## Restricts to the checked groups (defaults to HC+RA), falling back to all non-NA groups pre-render.
        wanted_groups <- input$example_groups %||% available_example_groups()
        keep <- !is.na(meta$group) & meta$group %in% wanted_groups
        meta <- meta[keep, , drop = FALSE]
        expr <- expr[, meta$sample, drop = FALSE]

        ## Same auto-detect log2 rule as every other preprocessing path in this module.
        q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
        if (isTRUE(!is.na(q99) && q99 > 100)) {
          expr[expr <= 0] <- NA
          expr <- log2(expr)
          expr <- expr[stats::complete.cases(expr), , drop = FALSE]
        }
        if (anyNA(expr)) {
          row_med <- apply(expr, 1, stats::median, na.rm = TRUE)
          na_idx <- which(is.na(expr), arr.ind = TRUE)
          expr[na_idx] <- row_med[na_idx[, 1]]
        }
        if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_
        list(expr = expr, meta = meta, label = pp_cohort_label(gse))
      })

      common <- Reduce(intersect, lapply(parts, function(x) rownames(x$expr)))
      validate(need(length(common) >= 20, "Fewer than 20 common genes between the two training datasets."))
      n_dup_features <- sum(vapply(parts, function(x) expr_raw_health(x$expr)$n_duplicated_features, integer(1)))
      merged_expr <- do.call(cbind, lapply(parts, function(x) x$expr[common, , drop = FALSE]))
      metas <- lapply(parts, function(x) { m <- x$meta; m$dataset <- x$label; m })
      all_cols <- unique(unlist(lapply(metas, colnames)))
      metas <- lapply(metas, function(m) {
        missing_cols <- setdiff(all_cols, colnames(m))
        for (cl in missing_cols) m[[cl]] <- NA
        m[, all_cols, drop = FALSE]
      })
      merged_meta <- do.call(rbind, metas)
      if (!"batch" %in% colnames(merged_meta) || all(is.na(merged_meta$batch))) merged_meta$batch <- merged_meta$dataset
      validate(need(identical(colnames(merged_expr), merged_meta$sample),
                    "Internal error: merged expression columns and metadata sample order do not match. Please report this as a bug."))
      list(expr = merged_expr, meta = merged_meta, sources = PP_TRAINING_COHORT_LABEL, n_dup_features = n_dup_features)
    })

    venn_plot_example_obj <- reactive({
      draw_overlap_venn(example_overlap_sets())
    })

    output$venn_plot_example <- renderPlot({
      venn_plot_example_obj()
    })

    output$download_venn_example_png <- downloadHandler(
      filename = function() "example_feature_overlap_venn.png",
      content = function(file) ggsave(file, plot = venn_plot_example_obj(), width = 7.5, height = 6.5, dpi = 300, bg = "white")
    )

    output$venn_table_example <- DT::renderDataTable({
      sets <- example_overlap_sets()
      common <- Reduce(intersect, sets)
      df <- data.frame(dataset = c(names(sets), "Common to all"), n_features = c(lengths(sets), length(common)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    venn_regions_example <- reactive(overlap_region_sizes(example_overlap_sets()))

    output$venn_region_table_example <- DT::renderDataTable({
      DT::datatable(venn_regions_example(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_venn_example <- downloadHandler(
      filename = function() "example_feature_overlap_regions.csv",
      content = function(file) write.csv(venn_regions_example(), file, row.names = FALSE)
    )

    
    output$merge_select_ui <- renderUI({
      
      lst <- tryCatch(merge_inputs(), error = function(e) NULL)
      req(length(lst) >= 2)
      labels <- vapply(lst, `[[`, character(1), "label")
      box(width = 12, title = tagList(icon("check-double"), " Choose which datasets to merge"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "All preprocessed datasets are selected by default. Uncheck any dataset you do not want included in this merge."),
          checkboxGroupInput(ns("merge_selected"), NULL, choices = labels, selected = labels, inline = TRUE))
    })

    ## Shared annotation table for optional probe-collapsing, reused across every selected dataset.
    collapse_annot <- reactive({
      validate(need(!is.null(input$collapse_annot_file),
                    "Probe-to-gene collapsing is turned on, but no annotation file has been uploaded yet - add one above, or turn the checkbox off if your data is already at gene level."))
      path <- input$collapse_annot_file$datapath
      if (grepl("\\.rds$", input$collapse_annot_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        d <- loaded$value
        validate(need(is.data.frame(d), "The annotation RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    selected_lst <- reactive({
      lst <- merge_inputs()
      labels <- vapply(lst, `[[`, character(1), "label")
      sel <- if (length(lst) < 2) labels else (input$merge_selected %||% labels)
      validate(need(length(sel) >= 1, "Select at least one dataset to merge."))
      lst <- lst[labels %in% sel]
      if (isTRUE(input$collapse_probes)) {
        annot <- collapse_annot()
        lst <- lapply(lst, function(x) {
          x$expr <- pp_collapse_probes_to_genes(x$expr, annot, input$collapse_method %||% "median")
          x
        })
      }
      lst
    })

    output$merge_venn_ui <- renderUI({
      
      lst <- tryCatch(selected_lst(), error = function(e) e)
      if (inherits(lst, "error")) {
        return(div(class = "empty-note", icon("circle-info"), conditionMessage(lst)))
      }
      if (length(lst) < 2) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Fewer than two datasets are selected, so there is nothing to compare. Click \"Merge datasets\" below to continue with just this one, or select more datasets above."))
      }
      tagList(
        box(width = 12, title = tagList(icon("diagram-project"), " Feature overlap across datasets"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Only features (genes/probes) present in every selected dataset are kept in the merge."),
            fluidRow(
              column(7,
                withSpinner(plotOutput(ns("venn_plot_custom"), height = 440), color = "#2563EB", type = 6),
                div(class = "table-toolbar", downloadButton(ns("download_venn_custom_png"), "Download diagram (PNG)", class = "btn-primary btn-sm"))
              ),
              column(5, DT::dataTableOutput(ns("venn_table_custom")))
            )),
        box(width = 12, title = tagList(icon("table-list"), " Region breakdown"), status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Every exact combination of datasets a feature belongs to, largest first. These are the same regions the diagram above draws, shown here as a full table."),
            div(class = "table-toolbar", downloadButton(ns("download_venn_custom"), "Download region counts (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("venn_region_table_custom")))
      )
    })

    overlap_sets <- reactive({
      lst <- selected_lst()
      validate(need(length(lst) >= 2, "Fewer than two datasets are selected, so there is nothing to compare."))
      setNames(lapply(lst, function(x) rownames(x$expr)), vapply(lst, `[[`, character(1), "label"))
    })

    venn_plot_custom_obj <- reactive({
      draw_overlap_venn(overlap_sets())
    })

    output$venn_plot_custom <- renderPlot({
      venn_plot_custom_obj()
    })

    output$download_venn_custom_png <- downloadHandler(
      filename = function() "feature_overlap_venn.png",
      content = function(file) ggsave(file, plot = venn_plot_custom_obj(), width = 7.5, height = 6.5, dpi = 300, bg = "white")
    )

    output$venn_table_custom <- DT::renderDataTable({
      sets <- overlap_sets()
      common <- Reduce(intersect, sets)
      df <- data.frame(dataset = c(names(sets), "Common to all"), n_features = c(lengths(sets), length(common)))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    venn_regions_custom <- reactive(overlap_region_sizes(overlap_sets()))

    output$venn_region_table_custom <- DT::renderDataTable({
      DT::datatable(venn_regions_custom(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_venn_custom <- downloadHandler(
      filename = function() "feature_overlap_regions.csv",
      content = function(file) write.csv(venn_regions_custom(), file, row.names = FALSE)
    )

    ## Triggered by either "Merge datasets" or the example-path button.
    merged <- eventReactive(list(input$merge_btn, input$merge_use_example_btn), {
      if (identical(input$merge_mode, "example")) {
        return(example_live_merge())
      }

      lst <- selected_lst()
      if (length(lst) == 1) {
        x <- lst[[1]]
        meta <- x$meta
        if (!"dataset" %in% colnames(meta)) meta$dataset <- x$label
        n_dup_features <- expr_raw_health(x$expr)$n_duplicated_features
        
        expr_dedup <- x$expr[!duplicated(rownames(x$expr)), , drop = FALSE]
        return(list(expr = expr_dedup, meta = meta, sources = x$label, n_dup_features = n_dup_features))
      }
      sets <- lapply(lst, function(x) rownames(x$expr))
      common <- Reduce(intersect, sets)
      validate(need(length(common) >= 20,
                    "Fewer than 20 features are in common across the selected datasets. Check that every uploaded dataset uses the same type of row name, for example all gene symbols or all the same probe IDs."))
      
      n_dup_features <- sum(vapply(lst, function(x) expr_raw_health(x$expr)$n_duplicated_features, integer(1)))
      merged_expr <- do.call(cbind, lapply(lst, function(x) x$expr[common, , drop = FALSE]))
      metas <- lapply(lst, function(x) { m <- x$meta; m$dataset <- x$label; m })
      all_cols <- unique(unlist(lapply(metas, colnames)))
      metas <- lapply(metas, function(m) {
        missing_cols <- setdiff(all_cols, colnames(m))
        for (cl in missing_cols) m[[cl]] <- NA
        m[, all_cols, drop = FALSE]
      })
      merged_meta <- tryCatch(do.call(rbind, metas), error = function(e) {
        validate("Could not combine metadata across datasets. A column with the same name has a different type in different datasets, for example numeric in one and text in another. Rename or fix that column, then preprocess again.")
      })
      if (!"batch" %in% colnames(merged_meta) || all(is.na(merged_meta$batch))) merged_meta$batch <- merged_meta$dataset
      validate(need(identical(colnames(merged_expr), merged_meta$sample),
                    "Internal error: merged expression columns and metadata sample order do not match. Please report this as a bug."))
      list(expr = merged_expr, meta = merged_meta, sources = paste(vapply(lst, `[[`, character(1), "label"), collapse = " + "), n_dup_features = n_dup_features)
    }, ignoreInit = TRUE)

    output$merge_summary_ui <- renderUI({
      if (!isTRUE((input$merge_btn %||% 0) > 0) && !isTRUE((input$merge_use_example_btn %||% 0) > 0)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Select datasets above, then click \"Merge datasets\"."))
      }
      m <- tryCatch(merged(), error = function(e) e)
      if (inherits(m, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), conditionMessage(m)))
      }
      tagList(
        div(class = "empty-note", icon("check"),
            sprintf("Merged: %s samples x %s common features, from %s.", ncol(m$expr), format(nrow(m$expr), big.mark = ","), m$sources)),
        if ((m$n_dup_features %||% 0) > 0) div(class = "empty-note", icon("triangle-exclamation"),
            sprintf("%d duplicated feature identifier(s) were detected across the merged dataset(s); only the first occurrence of each was kept.", m$n_dup_features)),
        div(class = "table-toolbar",
            downloadButton(ns("download_merged_expr"), "Download expression matrix (CSV)", class = "btn-sm"),
            downloadButton(ns("download_merged_meta"), "Download sample metadata (CSV)", class = "btn-sm"),
            downloadButton(ns("download_merged_rds"), "Download both (RDS)", class = "btn-sm")),
        p(class = "submodule-desc", "Composition by dataset (and group, if mapped):"),
        DT::dataTableOutput(ns("merge_composition_table"))
      )
    })

    output$download_merged_expr <- downloadHandler(
      filename = function() "merged_expression.csv",
      content = function(file) {
        m <- merged()
        write.csv(data.frame(feature = rownames(m$expr), m$expr, check.names = FALSE), file, row.names = FALSE)
      }
    )
    output$download_merged_meta <- downloadHandler(
      filename = function() "merged_metadata.csv",
      content = function(file) write.csv(merged()$meta, file, row.names = FALSE)
    )
    output$download_merged_rds <- downloadHandler(
      filename = function() "merged_dataset.rds",
      content = function(file) {
        m <- merged()
        saveRDS(list(expr = m$expr, meta = m$meta, sources = m$sources), file)
      }
    )

    output$merge_composition_table <- DT::renderDataTable({
      m <- merged()
      grp <- if ("group" %in% colnames(m$meta)) m$meta$group else rep("(no group column)", nrow(m$meta))
      df <- as.data.frame(table(dataset = m$meta$dataset, group = grp))
      colnames(df)[3] <- "n_samples"
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    

    active_meta_df <- reactive({
      m <- merged()
      m$meta
    })

    output$settings_ui <- renderUI({
      m <- tryCatch(merged(), error = function(e) NULL)
      if (is.null(m)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Finish the Merge datasets tab first (or preprocess and merge just one dataset there) before configuring batch correction."))
      }
      cols <- colnames(m$meta)
      batch_default <- intersect(c("batch", "batch_full", "dataset"), cols)
      protect_default <- intersect(c("group", "sex"), cols)
      
      max_pc <- max(2, min(5, ncol(m$expr)))
      pc_choices <- setNames(seq_len(max_pc), paste0("PC", seq_len(max_pc)))
      tagList(
        selectInput(ns("color_by"), "Color PCA by", choices = cols,
                    selected = if ("group" %in% cols) "group" else cols[1], selectize = FALSE),
        selectInput(ns("pc_x"), "X axis", choices = pc_choices, selected = 1, selectize = FALSE),
        selectInput(ns("pc_y"), "Y axis", choices = pc_choices, selected = min(2, max_pc), selectize = FALSE),
        checkboxInput(ns("show_ellipse"), "Show group confidence ellipses", value = TRUE),
        checkboxInput(ns("show_labels"), "Label points with sample ID", value = FALSE),
        selectInput(ns("batch_col"), "Batch column to correct for",
                    choices = cols,
                    selected = if (length(batch_default) > 0) batch_default[1] else cols[1], selectize = FALSE),
        uiOutput(ns("confound_ui")),
        {
          tagList(
            radioButtons(
              ns("norm_method"), "Normalisation",
              choiceNames = list(
                "Auto-detect (recommended)",
                "Already normalised, skip",
                "Quantile normalisation (microarray or log-scale data)",
                "TMM plus log2-CPM (RNA-seq raw counts)"
              ),
              choiceValues = list("auto", "skip", "quantile", "tmm"),
              selected = "auto"
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] != 'tmm'", ns("norm_method")),
              sliderInput(ns("min_pct"), "Drop genes below this percentile of mean expression", min = 0, max = 90, value = 0, step = 5)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'tmm'", ns("norm_method")),
              p(class = "empty-note", icon("circle-info"),
                "Low-count genes are filtered with edgeR::filterByExpr() by group instead of the percentile slider. Set log2 to \"Skip\" for this dataset in the Preprocessing tab first."),
              radioButtons(
                ns("tmm_correction_stage"), "Batch-correct raw counts, or the normalised log2-CPM?",
                choiceNames = list(
                  "After TMM: ComBat or limma on log2-CPM below (typical batch effects)",
                  "Before TMM: ComBat-seq on raw counts (Zhang, Parmigiani and Johnson 2020, better for large, count-driven batch effects)"
                ),
                choiceValues = list("post", "pre"), selected = "post"
              ),
              p(class = "empty-note", icon("circle-info"),
                "ComBat-seq (before TMM) ignores the correction method, prior, reference batch and exclude-outliers options below. It always protects the group column directly, and uses the batch column and biological covariates chosen here.")
            ),
            selectInput(ns("protect_cols"), "Biological covariates to protect (won't be treated as batch)",
                        choices = cols, selected = protect_default, multiple = TRUE),
            checkboxInput(ns("skip_combat"), "Skip batch correction (normalise only)", value = FALSE),
            tags$hr(),
            checkboxInput(ns("show_advanced"), strong("Show advanced batch-correction options"), value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s']", ns("show_advanced")),
              selectInput(ns("correction_method"), "Correction method",
                          choices = c("ComBat (empirical Bayes, recommended)" = "combat",
                                      "limma::removeBatchEffect (simple linear adjustment)" = "limma",
                                      "Surrogate Variable Analysis - SVA (unknown/hidden sources, no batch column needed)" = "sva"),
                          selected = "combat", selectize = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'combat'", ns("correction_method")),
                radioButtons(ns("combat_prior"), "ComBat empirical Bayes prior",
                             choiceNames = list("Parametric (faster, default)", "Non-parametric (slower, more robust for small or uneven batches)"),
                             choiceValues = list("param", "nonparam"), selected = "param"),
                checkboxInput(ns("combat_mean_only"), "Adjust batch mean only (skip variance adjustment)", value = FALSE),
                uiOutput(ns("ref_batch_ui"))
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'sva'", ns("correction_method")),
                p(class = "empty-note", icon("circle-info"),
                  "SVA estimates unwanted variation directly from the data instead of using the batch column above - useful when the real source of batch effects is unknown or only partly captured by a column you have. It still protects the biological covariates chosen below. Leek JT et al., Bioinformatics 2012;28(6):882-883."),
                numericInput(ns("sva_n_sv"), "Number of surrogate variables (0 = auto-estimate)", value = 0, min = 0, max = 20, step = 1)
              ),
              selectInput(ns("batch_col2"), "Combine with a second column into an interaction batch (optional, for example dataset by scan batch)",
                          choices = c("(none)", cols), selectize = FALSE),
              sliderInput(ns("variance_pct"), "Also drop genes below this percentile of variance", min = 0, max = 90, value = 0, step = 5),
              checkboxInput(ns("exclude_outliers"), "Exclude samples flagged as QC outliers before correcting (uses the outlier sensitivity below)", value = FALSE)
            )
          )
        },
        sliderInput(ns("mad_k"), "Outlier sensitivity (MADs from the cohort median)", min = 2, max = 6, value = 3, step = 0.5),
        uiOutput(ns("confound_override_ui")),
        actionButton(ns("run_btn"), "Run normalisation and batch correction",
                      icon = icon("play"), class = "btn-primary btn-sm")
      )
    })

    confounded_now <- reactive({
      meta <- tryCatch(active_meta_df(), error = function(e) NULL)
      req(meta, input$batch_col, input$batch_col %in% colnames(meta), "group" %in% colnames(meta))
      cc <- multi_live_confounding_check(meta, input$batch_col, "group")
      isTRUE(cc$confounded)
    })

    output$confound_ui <- renderUI({
      meta <- tryCatch(active_meta_df(), error = function(e) NULL)
      req(meta, input$batch_col, input$batch_col %in% colnames(meta))
      if (!"group" %in% colnames(meta)) return(NULL)
      cc <- multi_live_confounding_check(meta, input$batch_col, "group")
      if (is.null(cc)) return(NULL)
      if (isTRUE(cc$confounded)) {
        div(class = "empty-note", style = "border-color: var(--color-danger, #d9534f);", icon("triangle-exclamation"),
            " Potential confounding detected: the batch column and the group (phenotype) column are strongly associated - every batch level maps to essentially one group. Batch correction may remove genuine biological signal and cannot reliably separate batch from phenotype. Correction is blocked below unless you explicitly override this.")
      } else {
        div(class = "empty-note", icon("circle-check"), sprintf(" No strong batch/phenotype confounding detected (chi-square p = %.3f).", cc$p_value %||% NA))
      }
    })

    output$confound_override_ui <- renderUI({
      if (!isTRUE(tryCatch(confounded_now(), error = function(e) FALSE))) return(NULL)
      checkboxInput(ns("confound_override"), "I understand batch and the group column appear confounded and want to proceed anyway.", value = FALSE)
    })

    output$ref_batch_ui <- renderUI({
      req(input$batch_col)
      meta <- tryCatch(active_meta_df(), error = function(e) NULL)
      req(meta)
      lvls <- sort(unique(stats::na.omit(as.character(meta[[input$batch_col]]))))
      selectInput(ns("ref_batch"), "Reference batch (optional). Other batches are shifted to match this one instead of a pooled average.",
                  choices = c("(none)", lvls), selected = "(none)", selectize = FALSE)
    })

    
    observeEvent(list(input$batch_col, input$batch_col2), {
      meta <- tryCatch(active_meta_df(), error = function(e) NULL)
      req(meta)
      cols <- colnames(meta)
      exclude <- c(input$batch_col, if (!identical(input$batch_col2 %||% "(none)", "(none)")) input$batch_col2 else NULL)
      choices <- setdiff(cols, exclude)
      current <- intersect(input$protect_cols %||% character(0), choices)
      updateSelectInput(session, "protect_cols", choices = choices, selected = current)
    }, ignoreInit = TRUE)

    result <- eventReactive(input$run_btn, {
      req(input$batch_col, input$color_by)

      {
        skip_combat_early <- isTRUE(input$skip_combat)
        correction_method_early <- input$correction_method %||% "combat"
        using_batch_col_for_correction <- !skip_combat_early && (
          identical(correction_method_early, "combat") ||
          identical(correction_method_early, "limma") ||
          (identical(input$norm_method %||% "auto", "tmm") && identical(input$tmm_correction_stage %||% "post", "pre"))
        )
        if (using_batch_col_for_correction) {
          meta_cc <- tryCatch(active_meta_df(), error = function(e) NULL)
          if (!is.null(meta_cc) && input$batch_col %in% colnames(meta_cc) && "group" %in% colnames(meta_cc)) {
            cc <- multi_live_confounding_check(meta_cc, input$batch_col, "group")
            validate(need(
              is.null(cc) || !isTRUE(cc$confounded) || isTRUE(input$confound_override),
              "Batch correction is blocked: the batch column and the group (phenotype) column appear confounded (every batch level maps to a single group). This cannot reliably separate batch from phenotype and correction could remove genuine biological signal. Check the override box above \"Run normalisation and batch correction\" to proceed anyway, or choose SVA or a different batch column."
            ))
          }
        }
      }

      combat_fallback_note <- NULL

      {
        m <- merged()
        expr <- m$expr
        meta <- m$meta
        sources <- m$sources

        norm_method <- input$norm_method %||% "auto"
        skip_combat <- isTRUE(input$skip_combat)
        already_corrected <- FALSE

        if (identical(norm_method, "tmm")) {
          ## TMM + log2-CPM for raw RNA-seq counts: edgeR::filterByExpr() by group, then calcNormFactors(method="TMM").
          validate(need(all(expr >= 0, na.rm = TRUE),
                        "TMM normalisation expects raw, non-negative counts, but this data has negative values, which suggests it is already log-transformed. Preprocess this dataset again with log2 set to \"Skip\"."))
          
          non_integer_frac <- mean(abs(as.matrix(expr) - round(as.matrix(expr))) > 1e-6, na.rm = TRUE)
          validate(need(non_integer_frac < 0.01,
                        "TMM normalisation expects raw integer counts, but most values in this data are non-integer, which suggests it has already been normalised (e.g. CPM/RPKM/TPM, or quantile-normalised microarray intensities). Preprocess the original raw count matrix again with log2 set to \"Skip\"."))
          validate(need("group" %in% colnames(meta), "TMM normalisation needs a group column to filter low-count genes by."))
          grp <- factor(meta$group)
          validate(need(length(unique(na.omit(grp))) >= 2, "TMM normalisation needs at least two group levels."))
          counts <- round(as.matrix(expr))
          storage.mode(counts) <- "integer"
          dge0 <- edgeR::DGEList(counts = counts)
          keepg <- edgeR::filterByExpr(dge0, group = grp)
          n_before <- nrow(expr)
          validate(need(sum(keepg) >= 50, "Fewer than 50 genes pass edgeR's expression filter for this group split."))
          counts_f <- counts[keepg, , drop = FALSE]

          tmm_stage <- input$tmm_correction_stage %||% "post"
          if (!skip_combat && identical(tmm_stage, "pre")) {
            ## ComBat-seq: negative-binomial batch correction on raw counts before TMM (Zhang, Parmigiani & Johnson 2020).
            cs_batch_primary <- as.character(meta[[input$batch_col]])
            cs_use_batch2 <- !identical(input$batch_col2 %||% "(none)", "(none)") && (input$batch_col2 %in% colnames(meta))
            cs_batch <- if (cs_use_batch2) paste(cs_batch_primary, as.character(meta[[input$batch_col2]]), sep = "_") else cs_batch_primary
            validate(need(length(unique(na.omit(cs_batch))) >= 2, "The chosen batch column (or combination) needs at least two levels for ComBat-seq."))
            validate(need(all(table(cs_batch) >= 2), "Every level of the chosen batch column (or combination) needs at least 2 samples for ComBat-seq."))
            
            validate(need(!anyNA(cs_batch), sprintf(
              "%d sample(s) have no value in the chosen batch column (or combination) - every sample needs a batch value. Fix the metadata (e.g. map a batch column for every merged dataset) or choose a different batch column.",
              sum(is.na(cs_batch)))))

            cs_batch_cols_used <- c(input$batch_col, if (cs_use_batch2) input$batch_col2 else NULL)
            cs_protect <- intersect(input$protect_cols %||% character(0), colnames(meta))
            protect_dropped_for_batch <- intersect(setdiff(cs_protect, "group"), cs_batch_cols_used)
            cs_covar_cols <- setdiff(cs_protect, c("group", cs_batch_cols_used))
            cs_covar_cols <- cs_covar_cols[vapply(cs_covar_cols, function(cl) length(unique(na.omit(meta[[cl]]))) >= 2, logical(1))]
            cs_covar_mod <- if (length(cs_covar_cols) > 0) {
              meta_mod <- meta
              for (cl in cs_covar_cols) meta_mod[[cl]] <- ifelse(is.na(meta_mod[[cl]]), "Unknown", meta_mod[[cl]])
              
              covar_terms <- paste0("`", cs_covar_cols, "`")
              stats::model.matrix(stats::as.formula(paste("~", paste(covar_terms, collapse = " + "))), data = meta_mod)
            } else {
              NULL
            }

            cs_fell_back <- FALSE
            counts_adj <- tryCatch(
              sva::ComBat_seq(counts = counts_f, batch = cs_batch, group = as.character(grp), covar_mod = cs_covar_mod),
              error = function(e) {
                cs_fell_back <<- TRUE
                sva::ComBat_seq(counts = counts_f, batch = cs_batch_primary, group = as.character(grp))
              }
            )
            if (cs_fell_back) {
              cs_use_batch2 <- FALSE
              cs_covar_cols <- character(0)
              combat_fallback_note <- "ComBat-seq failed with the requested interaction batch column and/or protected covariates and was retried using only the primary batch column and the group label."
            }
            dge_before <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_f), method = "TMM")
            dge_after  <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_adj), method = "TMM")
            expr_prenorm <- counts_f
            expr_qnorm  <- edgeR::cpm(dge_before, log = TRUE, prior.count = 1)
            expr_combat <- edgeR::cpm(dge_after,  log = TRUE, prior.count = 1)
            already_corrected <- TRUE
            correction_method <- "combat_seq"
            protect <- intersect(c("group", cs_covar_cols), c("group", cs_protect))
            use_batch2 <- cs_use_batch2; ref_batch <- NULL
            combat_prior <- NA_character_; combat_mean_only <- FALSE
            norm_label <- "TMM plus log2-CPM, batch-corrected on raw counts with ComBat-seq before normalisation"
          } else {
            dge <- edgeR::calcNormFactors(dge0[keepg, , keep.lib.sizes = FALSE], method = "TMM")
            expr_prenorm <- counts_f
            expr_qnorm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
            norm_label <- "TMM (edgeR::calcNormFactors) plus log2-CPM"
          }
          needs_log <- FALSE; q99 <- NA_real_
          apply_qnorm <- TRUE
        } else {
          ## low-expression and low-variance gene filter, on the merged matrix
          n_before <- nrow(expr)
          gene_mean <- rowMeans(expr, na.rm = TRUE)
          gene_var  <- apply(expr, 1, stats::var, na.rm = TRUE)
          mean_cutoff <- stats::quantile(gene_mean, input$min_pct / 100, na.rm = TRUE)
          var_cutoff  <- stats::quantile(gene_var, (input$variance_pct %||% 0) / 100, na.rm = TRUE)
          keep <- !is.na(gene_mean) & !is.na(gene_var) & gene_mean >= mean_cutoff & gene_var > 0 & gene_var >= var_cutoff
          expr <- expr[keep, , drop = FALSE]
          validate(need(nrow(expr) >= 50, "Fewer than 50 genes remain after filtering. Lower the expression/variance percentile cutoffs."))
          expr_prenorm <- expr
          needs_log <- FALSE; q99 <- NA_real_

          diag_before <- summarize_norm_diagnostics(expr_prenorm)
          apply_qnorm <- switch(norm_method,
            skip = FALSE,
            quantile = TRUE,
            needs_quantile_norm(diag_before) # "auto"
          )
          expr_qnorm <- if (apply_qnorm) {
            mtx <- limma::normalizeBetweenArrays(as.matrix(expr_prenorm), method = "quantile")
            rownames(mtx) <- rownames(expr_prenorm); colnames(mtx) <- colnames(expr_prenorm)
            mtx
          } else {
            as.matrix(expr_prenorm)
          }
          norm_label <- switch(norm_method,
            skip = "None, used as loaded",
            quantile = "Quantile normalisation (forced)",
            if (apply_qnorm) "Quantile normalisation (auto-detected as needed)" else "None, auto-detected as already normalised"
          )
        }

        
        n_excluded_outliers <- 0L
        if (!already_corrected) {
          ## Optionally drop QC-flagged samples before correcting, not just flag them afterwards.
          if (isTRUE(input$exclude_outliers)) {
            qc_pre <- compute_sample_qc(expr_qnorm, mad_k = input$mad_k)
            flagged <- qc_pre$sample[qc_pre$flag_signal | qc_pre$flag_detected | qc_pre$flag_cor]
            if (length(flagged) > 0) {
              validate(need(ncol(expr_qnorm) - length(flagged) >= 6,
                            "Excluding QC-flagged samples would leave fewer than 6 samples. Lower the outlier sensitivity, or turn this option off."))
              keep_samples <- setdiff(colnames(expr_qnorm), flagged)
              expr_prenorm <- expr_prenorm[, keep_samples, drop = FALSE]
              expr_qnorm   <- expr_qnorm[, keep_samples, drop = FALSE]
              meta <- meta[match(keep_samples, meta$sample), , drop = FALSE]
              n_excluded_outliers <- length(flagged)
            }
          }

          batch_primary <- as.character(meta[[input$batch_col]])
          use_batch2 <- !identical(input$batch_col2 %||% "(none)", "(none)") && (input$batch_col2 %in% colnames(meta))
          batch <- if (use_batch2) paste(batch_primary, as.character(meta[[input$batch_col2]]), sep = "_") else batch_primary
          if (!skip_combat) {
            validate(need(length(unique(na.omit(batch))) >= 2, "The chosen batch column (or combination) needs at least two levels. If you don't need batch correction, tick \"Skip batch correction\" above."))
            validate(need(all(table(batch) >= 2), "Every level of the chosen batch column (or combination) needs at least 2 samples for correction."))
            
            validate(need(!anyNA(batch), sprintf(
              "%d sample(s) have no value in the chosen batch column (or combination) - every sample needs a batch value. Fix the metadata (e.g. map a batch column for every merged dataset) or choose a different batch column.",
              sum(is.na(batch)))))
          }

          protect <- intersect(input$protect_cols %||% character(0), colnames(meta))
          protect <- protect[vapply(protect, function(cl) length(unique(na.omit(meta[[cl]]))) >= 2, logical(1))]
          ## Drop batch columns from protect up front - protecting the corrected-for column makes batch/mod collinear.
          batch_cols_used <- c(input$batch_col, if (use_batch2) input$batch_col2 else NULL)
          protect_dropped_for_batch <- intersect(protect, batch_cols_used)
          protect <- setdiff(protect, batch_cols_used)
          mod <- if (length(protect) > 0) {
            meta_mod <- meta
            for (cl in protect) meta_mod[[cl]] <- ifelse(is.na(meta_mod[[cl]]), "Unknown", meta_mod[[cl]])
            
            protect_terms <- paste0("`", protect, "`")
            stats::model.matrix(stats::as.formula(paste("~", paste(protect_terms, collapse = " + "))), data = meta_mod)
          } else {
            NULL
          }

          combat_prior <- input$combat_prior %||% "param"
          combat_mean_only <- isTRUE(input$combat_mean_only)
          correction_method <- input$correction_method %||% "combat"
          ref_batch <- if (!use_batch2 && !identical(input$ref_batch %||% "(none)", "(none)")) input$ref_batch else NULL

          run_combat <- function(b, use_mod = TRUE, use_ref = TRUE) {
            sva::ComBat(dat = expr_qnorm, batch = b, mod = if (use_mod) mod else NULL,
                        par.prior = identical(combat_prior, "param"), mean.only = combat_mean_only,
                        ref.batch = if (use_ref) ref_batch else NULL)
          }
          run_limma <- function(b) {
            design <- if (!is.null(mod)) mod else matrix(1, ncol(expr_qnorm), 1)
            limma::removeBatchEffect(expr_qnorm, batch = b, design = design)
          }
          
          run_sva <- function() {
            mod_full <- if (!is.null(mod)) mod else matrix(1, ncol(expr_qnorm), 1)
            mod0 <- matrix(1, ncol(expr_qnorm), 1)
            ## vfilter restricts num.sv/sva() to the most variable genes so estimation stays fast on large matrices.
            n_genes <- nrow(expr_qnorm)
            vfilt <- if (n_genes > 2000) 2000L else NULL
            n_sv <- as.integer(input$sva_n_sv %||% 0)
            if (n_sv <= 0) {
              
              set.seed(ARTHOMIX_TX_ML_SEED)
              n_sv <- tryCatch(sva::num.sv(as.matrix(expr_qnorm), mod_full, method = "be", vfilter = vfilt),
                                error = function(e) NA_integer_)
            }
            n_sv <- if (is.na(n_sv)) 1L else max(1L, min(n_sv, ncol(expr_qnorm) - ncol(mod_full) - 1L, 20L))
            
            set.seed(ARTHOMIX_TX_ML_SEED)
            sv_obj <- sva::sva(as.matrix(expr_qnorm), mod_full, mod0, n.sv = n_sv, vfilter = vfilt)
            validate(need(sv_obj$n.sv >= 1, "SVA did not find any significant hidden sources of variation to correct for - try ComBat or limma instead, or set the number of surrogate variables manually."))
            limma::removeBatchEffect(expr_qnorm, covariates = sv_obj$sv, design = mod_full)
          }

          expr_combat <- if (skip_combat) {
            expr_qnorm
          } else if (identical(correction_method, "limma")) {
            tryCatch(run_limma(batch), error = function(e) {
              use_batch2 <<- FALSE
              combat_fallback_note <<- "limma::removeBatchEffect failed with the interaction of the two chosen batch columns and was retried using only the primary batch column."
              run_limma(batch_primary)
            })
          } else if (identical(correction_method, "sva")) {
            run_sva()
          } else {
            
            tryCatch(run_combat(batch, use_mod = TRUE, use_ref = TRUE),
              error = function(e) tryCatch({
                use_batch2 <<- FALSE
                ref_batch <<- NULL
                combat_fallback_note <<- "ComBat failed with the full model (interaction batch column plus reference batch) and was retried using only the primary batch column, without a reference batch."
                run_combat(batch_primary, use_mod = TRUE, use_ref = FALSE)
              }, error = function(e2) {
                protect <<- character(0)
                combat_fallback_note <<- "ComBat failed again without a reference batch and was retried using only the primary batch column, with no reference batch and no protected covariates."
                run_combat(batch_primary, use_mod = FALSE, use_ref = FALSE)
              }))
          }
        }
        n_after <- nrow(expr_prenorm)
      }

      diag_before <- summarize_norm_diagnostics(expr_prenorm)
      diag_after  <- summarize_norm_diagnostics(expr_qnorm)
      norm_diag <- rbind(cbind(stage = "before normalisation", diag_before),
                          cbind(stage = "after normalisation", diag_after))

      qc <- compute_sample_qc(expr_combat, mad_k = input$mad_k)
      qc <- merge(qc, meta[, intersect(c("sample", "group"), colnames(meta))], by = "sample", all.x = TRUE)

      before <- pca_of(expr_qnorm)
      after  <- pca_of(expr_combat)

      list(
        expr_prenorm = expr_prenorm, expr_qnorm = expr_qnorm, expr_combat = expr_combat,
        before = before, after = after, meta = meta, qc = qc,
        n_before = n_before, n_after = n_after,
        needs_log = needs_log, q99 = q99, apply_qnorm = apply_qnorm, norm_label = norm_label, norm_diag = norm_diag,
        protect = protect, protect_dropped_for_batch = protect_dropped_for_batch,
        batch_col = input$batch_col, color_by = input$color_by,
        skip_combat = skip_combat, combat_prior = combat_prior,
        combat_mean_only = combat_mean_only, sources = sources, correction_method = correction_method,
        n_excluded_outliers = n_excluded_outliers, use_batch2 = use_batch2, ref_batch = ref_batch,
        combat_fallback_note = combat_fallback_note
      )
    })

    ##  Value boxes 

    output$vb_samples <- renderValueBox({
      res <- result()
      valueBox(nrow(res$meta), "Samples", icon = icon("users"), color = "light-blue")
    })
    output$vb_genes_kept <- renderValueBox({
      res <- result()
      valueBox(format(res$n_after, big.mark = ","), "Genes retained", icon = icon("check"), color = "green")
    })
    output$vb_genes_dropped <- renderValueBox({
      res <- result()
      valueBox(format(res$n_before - res$n_after, big.mark = ","), "Genes filtered out", icon = icon("filter"), color = "purple")
    })
    output$vb_flagged <- renderValueBox({
      res <- result()
      n_flagged <- sum(res$qc$flag_signal | res$qc$flag_detected | res$qc$flag_cor)
      valueBox(n_flagged, "Samples flagged", icon = icon("triangle-exclamation"),
                color = if (n_flagged > 0) "red" else "green")
    })

    output$decisions_ui <- renderUI({
      res <- result()
      tagList(
        p(icon("circle-info"), " Normalisation: ", strong(res$norm_label), "."),
        p(icon("circle-info"), " Batch correction: ",
          if (isTRUE(res$skip_combat)) strong("skipped, normalisation only") else tagList(
            strong(switch(res$correction_method %||% "combat",
              combat_seq = "ComBat-seq, on raw counts before TMM normalisation",
              limma = "limma::removeBatchEffect",
              sva = "Surrogate Variable Analysis (unknown/hidden sources), regressed out via limma::removeBatchEffect",
              if (identical(res$combat_prior, "param")) "ComBat, parametric prior" else "ComBat, non-parametric prior"
            )),
            if (identical(res$correction_method, "combat") && isTRUE(res$combat_mean_only)) ", mean only" else "",
            if (isTRUE(res$use_batch2)) " on an interaction of the two chosen batch columns" else "",
            if (!is.null(res$ref_batch)) paste0(", referenced to batch \"", res$ref_batch, "\"") else ""
          ), "."),
        if (!isTRUE(res$skip_combat)) p(icon("circle-info"), " Protected ",
          if (length(res$protect) > 0) paste(res$protect, collapse = " and ") else "no biological covariates, since none were selected, or none had two or more levels in the current metadata",
          ", so that signal is not removed as if it were batch."),
        if (length(res$protect_dropped_for_batch) > 0) p(icon("triangle-exclamation"), " ",
          strong(paste(res$protect_dropped_for_batch, collapse = " and ")),
          if (length(res$protect_dropped_for_batch) > 1) " were" else " was",
          " also selected as the batch column (or part of it), so ", if (length(res$protect_dropped_for_batch) > 1) "they were" else "it was",
          " dropped from the protected covariates instead - a column can't be both corrected for and protected from correction at the same time."),
        if (isTRUE(res$n_excluded_outliers > 0)) p(icon("circle-info"), " Excluded ", strong(res$n_excluded_outliers),
          " sample(s) flagged as QC outliers before correcting."),
        if (!is.null(res$combat_fallback_note)) p(icon("triangle-exclamation"), " ", res$combat_fallback_note,
          " The settings above still show what was requested - this line reflects what actually ran.")
      )
    })

    ## QC plots 

    output$signal_plot <- renderPlot({
      res <- result()
      qc_bar_plot(res$qc, "signal", "flag_signal", "Total signal")
    })
    output$detected_plot <- renderPlot({
      res <- result()
      qc_bar_plot(res$qc, "detected", "flag_detected", "Detected features")
    })
    output$cor_plot <- renderPlot({
      res <- result()
      qc_bar_plot(res$qc, "mean_cor", "flag_cor", "Mean correlation")
    })

    dist_summary <- function(expr, meta) {
      qs <- apply(expr, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      df <- data.frame(sample = colnames(expr), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ],
                        upper = qs[4, ], ymax = qs[5, ])
      merge(df, meta[, intersect(c("sample", "group"), colnames(meta))], by = "sample", all.x = TRUE)
    }

    output$dist_plot <- renderPlot({
      res <- result()
      before_df <- dist_summary(res$expr_prenorm, res$meta); before_df$stage <- "Before normalisation"
      after_df  <- dist_summary(res$expr_qnorm, res$meta);   after_df$stage  <- "After normalisation"
      box_df <- rbind(before_df, after_df)
      box_df$stage <- factor(box_df$stage, levels = c("Before normalisation", "After normalisation"))
      if (!"group" %in% colnames(box_df)) box_df$group <- "all samples"

      ggplot(box_df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle,
                           upper = upper, ymax = ymax, fill = group)) +
        geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15) +
        scale_fill_manual(values = arthomix_pair(box_df$group)) +
        facet_wrap(~stage, ncol = 1, scales = "free_x") +
        labs(x = NULL, y = "Expression", fill = NULL) +
        theme_arthomix() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
    })

    output$pca_before <- renderPlot({
      res <- result()
      plot_pca_advanced(res$before, res$meta, res$color_by,
                          pc_x = as.integer(input$pc_x), pc_y = as.integer(input$pc_y),
                          title_suffix = "Before batch correction",
                          show_ellipse = input$show_ellipse, show_labels = input$show_labels)
    })
    output$pca_after <- renderPlot({
      res <- result()
      suffix <- if (isTRUE(res$skip_combat)) "After batch correction (skipped - identical to \"before\")"
                else paste0("After batch correction (vs. ", res$batch_col, ")")
      plot_pca_advanced(res$after, res$meta, res$color_by,
                          pc_x = as.integer(input$pc_x), pc_y = as.integer(input$pc_y),
                          title_suffix = suffix,
                          show_ellipse = input$show_ellipse, show_labels = input$show_labels)
    })
    output$scree_plot <- renderPlot({
      res <- result()
      scree_plot(res$after$var_exp)
    })

    assoc_pvalue <- function(pc1, batch) {
      df <- data.frame(pc1 = pc1, b = batch)
      df <- df[!is.na(df$b), ]
      tryCatch(summary(aov(pc1 ~ b, data = df))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    }

    output$summary_ui <- renderUI({
      res <- result()
      before_p <- assoc_pvalue(res$before$df$PC1, res$meta[[res$batch_col]])
      after_p  <- assoc_pvalue(res$after$df$PC1,  res$meta[[res$batch_col]])
      if (isTRUE(res$skip_combat)) {
        return(tagList(
          p("PC1 versus ", strong(res$batch_col), " association (one-way ANOVA p-value): ",
            sprintf("%.2g", before_p), " before, and ", sprintf("%.2g", after_p), " after."),
          p(class = "empty-note", icon("triangle-exclamation"),
            "These two numbers are identical because \"Skip batch correction\" is ticked in Settings, so \"after\" is the same uncorrected data as \"before\". Untick it to see whether correction actually reduces this column's effect on PC1.")
        ))
      }
      tagList(
        p("PC1 versus ", strong(res$batch_col), " association (one-way ANOVA p-value): ",
          sprintf("%.2g", before_p), " before correction, and ", sprintf("%.2g", after_p), " after."),
        p(class = "empty-note", icon("circle-info"),
          "A larger p-value after correction means PC1 is less explained by that column: the correction reduced its effect.")
      )
    })

    ##  Normalisation diagnostics table

    output$norm_table <- DT::renderDataTable({
      res <- result()
      DT::datatable(res$norm_diag, rownames = FALSE,
                     options = list(pageLength = 4, dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatRound(columns = c("max_value", "min_value", "median_sd", "iqr_sd", "median_range", "iqr_range"), digits = 3)
    })

    output$download_norm <- downloadHandler(
      filename = function() "normalization_diagnostics.csv",
      content = function(file) write.csv(result()$norm_diag, file, row.names = FALSE)
    )

    ##  QC table 

    qc_table_display <- reactive({
      res <- result()
      df <- res$qc
      df$reason <- apply(df[, c("flag_signal", "flag_detected", "flag_cor")], 1, function(r) {
        reasons <- c("low or high signal", "low or high detected features", "low cohort correlation")[r]
        if (length(reasons) == 0) "" else paste(reasons, collapse = "; ")
      })
      df$signal <- round(df$signal, 1)
      df$mean_cor <- round(df$mean_cor, 3)
      df[, c("sample", intersect("group", colnames(df)), "signal", "detected", "mean_cor", "reason")]
    })

    output$qc_table <- DT::renderDataTable({
      df <- qc_table_display()
      flagged_only <- df[df$reason != "", , drop = FALSE]
      DT::datatable(
        if (nrow(flagged_only) > 0) flagged_only else df[0, ],
        rownames = FALSE, filter = "top",
        options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact"
      )
    })

    output$download_qc <- downloadHandler(
      filename = function() "qc_metrics.csv",
      content = function(file) write.csv(qc_table_display(), file, row.names = FALSE)
    )

    ## PCA table 

    pca_table <- reactive({
      res <- result()
      rename_pc <- function(df, suffix) {
        pc_cols <- setdiff(colnames(df), "sample")
        colnames(df)[match(pc_cols, colnames(df))] <- paste0(pc_cols, "_", suffix)
        df
      }
      merge(rename_pc(res$before$df, "before"), rename_pc(res$after$df, "after"), by = "sample")
    })

    output$pca_table <- DT::renderDataTable({
      DT::datatable(pca_table(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_pca <- downloadHandler(
      filename = function() "pca_before_after.csv",
      content = function(file) write.csv(pca_table(), file, row.names = FALSE)
    )

    

    output$activate_ui <- renderUI({
      box(width = 12, title = "Use this dataset app-wide", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Once you're happy with the result above, make it the active dataset for every other sub-module (WGCNA, differential expression, feature selection, etc.)."),
          actionButton(ns("activate_btn"), "Use this as the active dataset", icon = icon("check"), class = "btn-primary btn-sm"),
          uiOutput(ns("activate_status_ui")))
    })

    observeEvent(input$activate_btn, {
      res <- result()
      dataset$expr <- res$expr_combat
      dataset$meta <- res$meta
      
      was_uploaded <- grepl("Uploaded dataset:", res$sources %||% "")
      was_geo <- !was_uploaded && grepl("NCBI GEO:", res$sources %||% "")
      dataset$source <- paste0(if (was_uploaded || was_geo) "Uploaded dataset (preprocessed + " else "Preloaded dataset (preprocessed + ",
                                if (isTRUE(res$skip_combat)) "normalised" else "batch-corrected", "): ",
                                res$sources %||% "your data")
      
      dataset$source_type <- if (was_uploaded) "uploaded" else if (was_geo) "geo" else "preloaded"
      dataset$is_bundled_reference <- FALSE
      output$activate_status_ui <- renderUI(
        div(class = "empty-note", icon("check"), "This is now the active dataset. Every other sub-module will use it.")
      )
    }, ignoreInit = TRUE)


    bc_section <- function(icon_name, title, ..., desc = NULL) {
      tagList(
        tags$h4(
          style = "margin: 20px 0 4px 0; padding-top: 14px; border-top: 1px solid var(--color-border); font-size: 15px; font-weight: 600; color: var(--color-ink); display: flex; align-items: center; gap: 8px;",
          icon(icon_name), title
        ),
        if (!is.null(desc)) p(class = "submodule-desc", desc),
        ...
      )
    }

    output$results_top_ui <- renderUI({
      res <- tryCatch(result(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"),
            "Set the options on the left, then click \"Run normalisation and batch correction\" to see results here."))
      }
      tagList(

        fluidRow(
          valueBoxOutput(ns("vb_samples"), width = 6),
          valueBoxOutput(ns("vb_genes_kept"), width = 6)
        ),
        fluidRow(
          valueBoxOutput(ns("vb_genes_dropped"), width = 6),
          valueBoxOutput(ns("vb_flagged"), width = 6)
        ),
        bc_section("clipboard-list", "Pipeline decisions",
          withSpinner(uiOutput(ns("decisions_ui")), color = "#2563EB", type = 6)
        )
      )
    })

    output$results_rest_ui <- renderUI({
      res <- tryCatch(result(), error = function(e) NULL)
      req(res)
      tagList(
        
        bc_section("chart-simple", "Per-sample signal",
          withSpinner(plotOutput(ns("signal_plot"), height = 210), color = "#2563EB", type = 6)
        ),
        bc_section("chart-simple", "Detected features per sample",
          withSpinner(plotOutput(ns("detected_plot"), height = 210), color = "#2563EB", type = 6)
        ),
        bc_section("chart-simple", "Mean correlation to cohort",
          withSpinner(plotOutput(ns("cor_plot"), height = 210), color = "#2563EB", type = 6)
        ),
        bc_section("chart-area", "Expression distribution: before and after normalisation",
          withSpinner(plotOutput(ns("dist_plot"), height = 300), color = "#2563EB", type = 6)
        ),
        bc_section("chart-line", "Scree plot",
          desc = "Percentage of variance explained per component, after batch correction. A sharp drop after PC2 or PC3 means most of the structure is captured by the axes plotted below.",
          withSpinner(plotOutput(ns("scree_plot"), height = 220), color = "#2563EB", type = 6)
        ),
        if (isTRUE(res$skip_combat)) div(class = "empty-note", icon("triangle-exclamation"),
          strong(" \"Skip batch correction\" is ticked in Settings: "),
          "no correction was run, so the panels below are the same normalised data twice, and every before/after comparison is necessarily identical. Untick it and re-run to see an actual before/after."),
        bc_section("braille", "PCA before batch correction",
          withSpinner(plotOutput(ns("pca_before"), height = 340), color = "#2563EB", type = 6)
        ),
        bc_section("braille", "PCA after batch correction",
          withSpinner(plotOutput(ns("pca_after"), height = 340), color = "#2563EB", type = 6)
        ),
        bc_section("scale-balanced", "Batch effect, before and after",
          withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6)
        ),
        bc_section("table-list", "Normalisation diagnostics",
          desc = "Per-sample median and IQR agreement across the cohort. Large 'after' values mean normalisation didn't fully align the samples.",
          div(class = "table-toolbar", downloadButton(ns("download_norm"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("norm_table"))
        ),
        bc_section("triangle-exclamation", "Flagged samples",
          desc = "Samples flagged by signal, detected-feature count, or cohort correlation, with the reason.",
          div(class = "table-toolbar", downloadButton(ns("download_qc"), "Download full QC table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("qc_table"))
        ),
        bc_section("table-cells", "PCA coordinates",
          div(class = "table-toolbar", downloadButton(ns("download_pca"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("pca_table"))
        ),
        uiOutput(ns("activate_ui"))
      )
    })

    
    batch_content <- tagList(
      fluidRow(
        column(
          6,
          bc_section("sliders", "Settings", uiOutput(ns("settings_ui")))
        ),
        column(
          6,
          withSpinner(uiOutput(ns("results_top_ui")), color = "#2563EB", type = 6)
        )
      ),
      withSpinner(uiOutput(ns("results_rest_ui")), color = "#2563EB", type = 6)
    )

    output$batch_tab_ui <- renderUI({
      tagList(
        div(class = "empty-note", icon("circle-info"),
            "Runs on whatever you merged in the Merge datasets tab (or your single preprocessed dataset, if you only configured one)."),
        batch_content
      )
    })

     
    for (bc_out in c("batch_tab_ui", "settings_ui", "ref_batch_ui", "results_top_ui", "results_rest_ui")) {
      outputOptions(output, bc_out, suspendWhenHidden = FALSE)
    }

    output$meta_tab_ui <- renderUI({
      m <- tryCatch(merged(), error = function(e) NULL)
      intro <- div(class = "empty-note", icon("circle-info"),
        "For studies too different to pool - analyse each separately, then combine effect sizes via random-effects meta-analysis.")
      if (is.null(m)) {
        return(tagList(intro, div(class = "empty-note", icon("triangle-exclamation"),
                                  "Finish the Merge datasets tab first - this tab needs two or more merged studies.")))
      }
      cols <- colnames(m$meta)
      if (!"group" %in% cols) return(tagList(intro, div(class = "empty-note", icon("triangle-exclamation"), "The merged metadata needs a group column.")))
      groups <- sort(unique(stats::na.omit(as.character(m$meta$group))))
      study_default <- intersect(c("dataset", "study", "cohort", "batch"), cols)[1]
      tagList(
        intro,
        fluidRow(
          column(6, bc_section("sliders", "Settings",
            selectInput(ns("meta_study_col"), "Study column (one level per independent study)", choices = cols,
                        selected = if (is.na(study_default)) cols[1] else study_default, selectize = FALSE),
            selectInput(ns("meta_ref_group"), "Reference group", choices = groups, selected = groups[1], selectize = FALSE),
            selectInput(ns("meta_comp_group"), "Comparison group", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE),
            numericInput(ns("meta_min_n"), "Minimum samples per group within a study", value = 3, min = 2, step = 1),
            actionButton(ns("meta_run_btn"), "Run study-wise meta-analysis", icon = icon("play"), class = "btn-primary btn-sm")
          )),
          column(6, withSpinner(uiOutput(ns("meta_summary_ui")), color = "#2563EB", type = 6))
        ),
        withSpinner(uiOutput(ns("meta_results_ui")), color = "#2563EB", type = 6)
      )
    })

    meta_result <- eventReactive(input$meta_run_btn, {
      m <- merged()
      req(input$meta_study_col, input$meta_ref_group, input$meta_comp_group)
      validate(need("group" %in% colnames(m$meta), "The merged metadata needs a group column."))
      validate(need(input$meta_ref_group != input$meta_comp_group, "Reference and comparison group must be different."))
      validate(need(input$meta_study_col %in% colnames(m$meta), "Choose a study column."))
      n_studies <- length(unique(stats::na.omit(as.character(m$meta[[input$meta_study_col]]))))
      validate(need(n_studies >= 2, "Study-wise meta-analysis needs at least two studies in the chosen study column. With a single study, use the Batch correction tab (tick \"Skip batch correction\" if there is no batch)."))
      min_n <- max(2L, as.integer(input$meta_min_n %||% 3))
      sw <- withProgress(message = "Fitting one limma model per study...", value = 0.3,
        pp_study_wise_limma(m$expr, m$meta, input$meta_study_col, input$meta_ref_group, input$meta_comp_group, min_n = min_n))
      validate(need(length(sw$per_study) >= 2, sprintf(
        "Fewer than two studies have at least %d samples in both groups (skipped: %s).", min_n,
        if (length(sw$skipped)) paste(sw$skipped, collapse = ", ") else "none")))
      meta_tbl <- pp_dl_meta(sw$per_study)
      validate(need(nrow(meta_tbl) > 0, "No gene was measured in at least two studies."))
      res <- list(per_study = sw$per_study, skipped = sw$skipped, table = meta_tbl,
                  ref_group = input$meta_ref_group, comp_group = input$meta_comp_group,
                  study_col = input$meta_study_col, run_at = Sys.time())
      if (!is.null(results)) {
        results$study_meta <- list(
          contrast = sprintf("%s vs %s", input$meta_comp_group, input$meta_ref_group),
          studies = names(sw$per_study), n_genes = nrow(meta_tbl),
          n_significant = sum(meta_tbl$adj_p < 0.05, na.rm = TRUE), timestamp = Sys.time())
      }
      res
    }, ignoreInit = TRUE)

    output$meta_summary_ui <- renderUI({
      r <- tryCatch(meta_result(), error = function(e) NULL)
      if (is.null(r)) return(div(class = "empty-note", icon("circle-info"), "Not run yet. Choose the study column and groups, then click \"Run study-wise meta-analysis\"."))
      n_sig <- sum(r$table$adj_p < 0.05, na.rm = TRUE)
      n_het <- sum(r$table$I2 > 50, na.rm = TRUE)
      tagList(
        p(strong("Studies combined: "), paste(names(r$per_study), collapse = ", "),
          if (length(r$skipped)) sprintf(" (skipped for too few samples: %s)", paste(r$skipped, collapse = ", ")) else ""),
        p(strong(format(nrow(r$table), big.mark = ",")), " genes measured in at least two studies; ",
          strong(format(n_sig, big.mark = ",")), sprintf(" with pooled FDR < 0.05 (%s vs %s).", r$comp_group, r$ref_group)),
        p(strong(format(n_het, big.mark = ",")), " genes show substantial between-study heterogeneity (I2 > 50%) - interpret their pooled effect with caution."),
        p(class = "submodule-desc", "No batch correction was applied: each study is modelled separately, so study-level technical differences cannot be mistaken for biology the way they can after pooling. The cost is lower power than a compatible pooled analysis.")
      )
    })

    output$meta_study_table <- DT::renderDataTable({
      r <- meta_result()
      df <- do.call(rbind, lapply(names(r$per_study), function(s) {
        d <- r$per_study[[s]]
        data.frame(study = s, n_ref = d$n_ref[1], n_comp = d$n_comp[1], n_genes = nrow(d), model = d$method[1], stringsAsFactors = FALSE)
      }))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$meta_table <- DT::renderDataTable({
      r <- meta_result()
      num <- setdiff(names(r$table)[vapply(r$table, is.numeric, logical(1))], "n_studies")
      DT::datatable(r$table, rownames = FALSE, filter = "top",
                    options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatSignif(columns = num, digits = 4)
    })

    output$meta_volcano <- renderPlot({
      r <- meta_result()
      df <- r$table
      df$het <- ifelse(is.na(df$I2), "NA", ifelse(df$I2 > 50, "I2 > 50%", "I2 <= 50%"))
      ggplot(df, aes(x = logFC_pooled, y = -log10(pmax(p_value, 1e-300)), color = het)) +
        geom_point(alpha = 0.6, size = 1.4) +
        scale_color_manual(values = c("I2 <= 50%" = ARTHOMIX_COLORS$blue, "I2 > 50%" = ARTHOMIX_COLORS$orange, "NA" = "grey70"), name = "Heterogeneity") +
        labs(x = sprintf("Pooled log2 fold-change (%s vs %s, random effects)", r$comp_group, r$ref_group), y = "-log10 pooled p-value") +
        theme_arthomix(base_size = 12)
    })

    output$meta_download <- downloadHandler(
      filename = function() "study_wise_meta_analysis.csv",
      content = function(file) write.csv(meta_result()$table, file, row.names = FALSE)
    )

    output$meta_results_ui <- renderUI({
      r <- tryCatch(meta_result(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      tagList(
        bc_section("table", "Per-study models", DT::dataTableOutput(ns("meta_study_table"))),
        bc_section("chart-simple", "Pooled effects", plotOutput(ns("meta_volcano"), height = 380)),
        bc_section("table-list", "Meta-analysis table",
                   div(class = "table-toolbar", downloadButton(ns("meta_download"), "Download CSV", class = "btn-sm")),
                   DT::dataTableOutput(ns("meta_table")))
      )
    })
    outputOptions(output, "meta_tab_ui", suspendWhenHidden = FALSE)

  
    mod_data_exploration_server("eda")
  })
}
