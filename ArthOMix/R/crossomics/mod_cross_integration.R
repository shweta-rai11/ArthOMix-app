## R/crossomics/mod_cross_integration.R
## Cross-Omics sub-module: "Expression x Methylation" - integrates the
## Transcriptomics DGE output (gene, log2FC, FDR) and the Methylomics DMP
## output (CpG, gene, Δβ, FDR) at gene level to answer "how does DNA
## methylation relate to gene expression" - the four regulatory quadrants
## (hyper+down, hypo+up, hyper+up, hypo+down) plus a not-significant/
## discordant bucket. Fully self-contained: Preloaded mode reads this
## project's own precomputed DEG/DMP tables directly (cx_load_default_deg()/
## cx_load_default_methylation(), crossomics_integration_helpers.R) rather
## than depending on whatever is currently loaded in the live Transcriptomics/
## Methylomics tabs, so this module never touches those tabs' reactive state.
## Upload mode is fully independent CSV/TSV/TXT/XLSX with auto-detected +
## manually overridable column mapping (crossomics_integration_upload.R).
##
## Every statistical/biological framing follows the "association, not
## causality" requirement throughout: category labels say "potential
## methylation-associated repression/activation", never "causes silencing",
## and sample-level correlation is only ever computed when a real, paired
## sample-level match is detected (see cx_detect_sample_pairing()) -
## otherwise the UI shows the mandated "Unpaired datasets detected" banner
## rather than quietly faking it.

mod_cross_integration_config <- list(
  id = "integration", title = "Expression x Methylation", icon = "dna", group = "Data",
  description = "Gene-level integration of Transcriptomics differential expression and Methylomics differential methylation - quadrant plot, heatmap, correlation, pathways, and downloadable results."
)

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_cross_integration_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Data Input", status = "primary", solidHeader = FALSE,
        radioButtons(ns("sex_stratum"), "Analysis group",
                     choices = c("ALL" = "all", "FEMALE" = "female", "MALE" = "male"),
                     selected = "female", inline = TRUE),
        radioButtons(ns("mode"), "Input type",
                     choices = c("Preloaded data" = "preloaded", "From Dataset tab" = "fromdataset", "Load from app" = "fromapp", "Upload your data" = "upload"),
                     selected = "preloaded"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'preloaded'", ns("mode")),
          p(class = "empty-note", icon("circle-info"),
            "Reads this project's own precomputed, sex-stratified DEG (Transcriptomics) and DMP (Methylomics) tables - independent of whatever is currently loaded on the Transcriptomics/Methylomics tabs."),
          actionButton(ns("load_preloaded"), "Load preloaded data", icon = icon("database"), class = "btn-primary btn-sm")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'fromdataset'", ns("mode")),
          uiOutput(ns("fromdataset_ui"))
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'fromapp'", ns("mode")),
          uiOutput(ns("fromapp_ui"))
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload'", ns("mode")),
          fileInput(ns("expr_file"), "Transcriptomics file", accept = c(".csv", ".tsv", ".txt", ".xlsx"),
                    placeholder = "CSV / TSV / TXT / XLSX"),
          p(class = "empty-note", icon("circle-info"), "Differentially Expressed Genes (DEG) format - one row per gene, with a gene symbol/ID and a log2 fold-change column."),
          fileInput(ns("meth_file"), "Methylomics file", accept = c(".csv", ".tsv", ".txt", ".xlsx"),
                     placeholder = "CSV / TSV / TXT / XLSX"),
          p(class = "empty-note", icon("circle-info"), "Differentially Methylated Position/Region (DMP/DMR) format - one row per CpG or region, with a gene symbol/ID and a Δβ (methylation change) column.")
        )
      ),
      box(
        width = NULL, title = "2. Map Your Data", status = "primary", solidHeader = FALSE,
        uiOutput(ns("mapping_ui"))
      ),
      box(
        width = NULL, title = "3. Integration Setup", status = "primary", solidHeader = FALSE,
        fluidRow(
          column(6, numericInput(ns("expr_thresh"), "Min |log2FC|", value = 1, min = 0, step = 0.1)),
          column(6, numericInput(ns("expr_fdr_thresh"), "Max expression FDR", value = 0.05, min = 0, max = 1, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput(ns("meth_thresh"), "Min |Δβ|", value = 0.10, min = 0, max = 1, step = 0.01)),
          column(6, numericInput(ns("meth_fdr_thresh"), "Max methylation FDR", value = 0.05, min = 0, max = 1, step = 0.01))
        ),
        selectInput(ns("agg_method"), "Methylation aggregation (multiple CpGs/gene)",
                    choices = setNames(names(CX_AGGREGATION_METHODS), CX_AGGREGATION_METHODS), selected = "mean"),
        radioButtons(ns("cor_method"), "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"), inline = TRUE),
        radioButtons(ns("padj_method"), "Multiple-testing adjustment", choices = c("Benjamini-Hochberg" = "BH", "Bonferroni" = "bonferroni"), inline = TRUE),
        tags$hr(),
        checkboxInput(ns("show_sig_only"), "Show significant genes only (plots)", value = FALSE),
        checkboxInput(ns("show_labels"), "Show gene labels", value = FALSE),
        checkboxInput(ns("show_nonsig"), "Show non-significant genes", value = TRUE),
        checkboxInput(ns("show_quadrant_lines"), "Show quadrant boundaries", value = TRUE),
        tags$hr(),
        actionButton(ns("run_integration"), "Run Integration", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")
      ),
      box(
        width = NULL, title = "4. Advanced Filters", status = "primary", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
        selectizeInput(ns("filter_region"), "Genomic region", choices = CX_REGION_FINE_VOCAB, multiple = TRUE, options = list(placeholder = "Any region")),
        selectizeInput(ns("filter_island"), "CpG island status", choices = CX_ISLAND_VOCAB, multiple = TRUE, options = list(placeholder = "Any island status")),
        selectizeInput(ns("filter_evidence"), "Evidence level", choices = CX_EVIDENCE_LEVELS, multiple = TRUE, options = list(placeholder = "Any evidence level")),
        numericInput(ns("filter_min_cpg"), "Minimum CpG count", value = 0, min = 0, step = 1),
        radioButtons(ns("filter_cor_direction"), "Correlation direction (where computed)", choices = c("Any" = "any", "Positive" = "pos", "Negative" = "neg"), inline = TRUE)
      ),
      box(
        width = NULL, title = "Analysis Settings / Reproducibility", status = "primary", solidHeader = FALSE,
        uiOutput(ns("provenance_ui"))
      )
    ),
    column(
      8,
      uiOutput(ns("status_bar")),
      uiOutput(ns("summary_cards")),
      tabsetPanel(
        id = ns("result_tabs"), type = "tabs",
        tabPanel("Validation", br(), uiOutput(ns("validation_ui"))),
        tabPanel("Results Table", br(), uiOutput(ns("results_table_ui"))),
        tabPanel("CpG-Level", br(), uiOutput(ns("cpg_level_ui"))),
        tabPanel("Quadrant Plot", br(), uiOutput(ns("quadrant_ui"))),
        tabPanel("Heatmap", br(), uiOutput(ns("heatmap_ui"))),
        tabPanel("Volcano Plots", br(), uiOutput(ns("volcano_ui"))),
        tabPanel("Correlation", br(), uiOutput(ns("correlation_ui"))),
        tabPanel("Overlap", br(), uiOutput(ns("overlap_ui"))),
        tabPanel("Pathways", br(), uiOutput(ns("pathway_ui"))),
        tabPanel("Network", br(), uiOutput(ns("network_ui"))),
        tabPanel("Genomic View", br(), uiOutput(ns("genomic_ui"))),
        tabPanel("Sex Comparison", br(), uiOutput(ns("sex_comparison_ui"))),
        tabPanel("Downloads", br(), uiOutput(ns("downloads_ui"))),
        tabPanel("Methodology & References", br(), cx_methodology_references_ui())
      )
    )
  )
}

## ---------------------------------------------------------------------------
## Small UI-building helpers (module-local, not shared elsewhere)
## ---------------------------------------------------------------------------

.cx_field_label <- function(field) {
  labels <- c(gene = "Gene ID / Gene Symbol", cpg = "CpG ID / Probe ID", log2fc = "log2 Fold Change",
              dbeta = "Δβ (methylation change)", beta = "Beta value", pvalue = "P-value", fdr = "Adjusted P-value (FDR)",
              chr = "Chromosome", pos = "Genomic position", region = "Region / Feature",
              island = "CpG Island / Relation to Island", sample_id = "Sample ID")
  labels[[field]] %||% field
}

.cx_mapping_selects <- function(ns, prefix, df, mapping, fields) {
  cols <- colnames(df)
  lapply(fields, function(f) {
    selectInput(ns(paste0(prefix, f)), .cx_field_label(f),
                choices = c("(none)" = "", setNames(cols, cols)),
                selected = if (!is.na(mapping[[f]])) mapping[[f]] else "")
  })
}

.cx_read_mapping <- function(input, prefix, fields) {
  vals <- vapply(fields, function(f) {
    v <- input[[paste0(prefix, f)]]
    if (is.null(v) || !nzchar(v)) NA_character_ else v
  }, character(1))
  setNames(vals, fields)
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_cross_integration_server <- function(id, cross_dataset, cross_results,
                                          dataset = NULL, results = NULL,
                                          methyl_dataset = NULL, methyl_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw <- reactiveValues(
      expr_df = NULL, expr_source = NULL, expr_sample_cols = character(0), expr_wide = NULL, expr_mapping = NULL,
      meth_df = NULL, meth_source = NULL, meth_sample_cols = character(0), meth_wide = NULL, meth_mapping = NULL,
      meth_unavailable_reason = NULL
    )
    integ <- reactiveValues(df = NULL, params = NULL, pairing = NULL, provenance = NULL, run_at = NULL,
                             universe = NULL, cpg_level = NULL, id_harmonization = NULL, validation = NULL)
    ## One completed run per stratum (spec sections 13/14/41) - kept
    ## independent of `integ` above (which always reflects the most recently
    ## run stratum) so "Sex Comparison" can compare across runs without
    ## forcing a re-run or a tripled UI.
    integ_by_sex <- reactiveValues(all = NULL, female = NULL, male = NULL)
    active_filter <- reactiveVal("All")
    pathway_result <- reactiveVal(NULL)
    network_hint <- reactiveVal(NULL)

    expr_upload <- reactiveVal(NULL)
    meth_upload <- reactiveVal(NULL)

    ## ---- Data Input: Preloaded ------------------------------------------

    observeEvent(input$load_preloaded, {
      sex <- input$sex_stratum
      deg <- cx_load_default_deg(sex = sex)
      if (is.null(deg)) {
        showNotification("Could not read the preloaded Transcriptomics (DEG) table for this sex stratum.", type = "error")
        return()
      }
      std <- cx_standardize_expression(deg, mapping = c(gene = "gene", log2fc = "logFC", pvalue = "P.Value", fdr = "adj.P.Val"))
      if (!std$ok) { showNotification(std$error, type = "error"); return() }
      raw$expr_df <- std$df
      raw$expr_source <- sprintf("Preloaded DEG table (%s, sex-stratified)", toupper(sex))
      raw$expr_sample_cols <- character(0)
      raw$expr_wide <- NULL

      meth <- cx_load_default_methylation(sex = sex)
      if (!meth$ok) {
        showNotification(meth$error, type = "warning", duration = 10)
        raw$meth_df <- NULL
        raw$meth_source <- NULL
        raw$meth_unavailable_reason <- meth$error
      } else {
        raw$meth_df <- meth$df
        raw$meth_source <- sprintf("Preloaded DMP table (%s, sex-stratified, SVA/bacon-adjusted)", toupper(sex))
        raw$meth_unavailable_reason <- NULL
      }
      raw$meth_sample_cols <- character(0)
      raw$meth_wide <- NULL
      showNotification("Preloaded data loaded.", type = "message")
    })

    ## ---- Data Input: From Dataset tab --------------------------------------
    ## Reads whatever the Dataset tab's "Upload your own data" mode published
    ## into the shared `cross_dataset` store - already standardized and
    ## auto-mapped there, so this is a direct hand-off, not a re-parse.

    output$fromdataset_ui <- renderUI({
      has_expr <- !is.null(cross_dataset$user_expr_df)
      has_meth <- !is.null(cross_dataset$user_meth_df)
      tagList(
        if (has_expr) div(class = "empty-note", icon("check"),
          sprintf("Transcriptomics ready: %s (%s genes).", cross_dataset$user_expr_source, format(nrow(cross_dataset$user_expr_df), big.mark = ","))
        ) else div(class = "empty-note", icon("circle-info"), "No Transcriptomics file uploaded yet."),
        if (has_meth) div(class = "empty-note", icon("check"),
          sprintf("Methylomics ready: %s (%s records).", cross_dataset$user_meth_source, format(nrow(cross_dataset$user_meth_df), big.mark = ","))
        ) else div(class = "empty-note", icon("circle-info"), "No Methylomics file uploaded yet."),
        if (!has_expr || !has_meth) div(class = "empty-note", icon("triangle-exclamation"),
          "Go to the Dataset tab (above) and upload the missing file(s) under \"Upload your own data\", then come back here."),
        actionButton(ns("load_from_dataset"), "Use this data", icon = icon("share-from-square"), class = "btn-primary btn-sm")
      )
    })

    observeEvent(input$load_from_dataset, {
      validate(need(!is.null(cross_dataset$user_expr_df) || !is.null(cross_dataset$user_meth_df),
                    "Upload a file on the Dataset tab first."))
      if (!is.null(cross_dataset$user_expr_df)) {
        raw$expr_df <- cross_dataset$user_expr_df
        raw$expr_source <- cross_dataset$user_expr_source
        raw$expr_wide <- cross_dataset$user_expr_wide
        raw$expr_mapping <- cross_dataset$user_expr_mapping
        raw$expr_sample_cols <- cross_dataset$user_expr_sample_cols %||% character(0)
      }
      if (!is.null(cross_dataset$user_meth_df)) {
        raw$meth_df <- cross_dataset$user_meth_df
        raw$meth_source <- cross_dataset$user_meth_source
        raw$meth_wide <- cross_dataset$user_meth_wide
        raw$meth_mapping <- cross_dataset$user_meth_mapping
        raw$meth_sample_cols <- cross_dataset$user_meth_sample_cols %||% character(0)
        raw$meth_unavailable_reason <- NULL
      }
      showNotification("Loaded your data from the Dataset tab.", type = "message")
    })

    ## ---- Data Input: Load from app (spec section 4, Option A) -------------
    ## "Load from Transcriptomics" reads results$dge_runs, which DOES carry a
    ## full per-gene table (gene, logFC, adj.P.Val, direction - mod_dge.R) so
    ## this fully works. "Load from Methylomics" is honestly limited: every
    ## place Methylomics publishes into methyl_results/results
    ## (mod_methyl_dmp.R, mod_methyl_candidates.R) exposes only summary
    ## counts (comparison label, n_probes, n_sig) - never the full per-CpG
    ## table, which stays local to that module. Per the non-negotiable rule
    ## that Methylomics module files are read-only, this module cannot
    ## reach that table, so the button surfaces the summary honestly and
    ## points to Preloaded/Upload instead of fabricating rows.

    output$fromapp_ui <- renderUI({
      runs <- results$dge_runs
      tx_choices <- if (!is.null(results) && length(runs) > 0) {
        setNames(names(runs), vapply(runs, function(r) sprintf("%s (%s, n=%s)", r$contrast, r$method, r$n_samples), character(1)))
      } else character(0)
      tagList(
        h5("Transcriptomics"),
        if (length(tx_choices) == 0) {
          div(class = "empty-note", icon("circle-info"), "No Transcriptomics run available yet - run a contrast in the Transcriptomics → Differential Expression tab first.")
        } else tagList(
          selectInput(ns("tx_run_pick"), "Contrast run", choices = tx_choices),
          actionButton(ns("load_from_tx"), "Load from Transcriptomics", icon = icon("share-from-square"), class = "btn-primary btn-sm")
        ),
        tags$hr(),
        h5("Methylomics"),
        if (is.null(methyl_results) || is.null(methyl_results$dmp)) {
          div(class = "empty-note", icon("circle-info"), "No Methylomics run available yet - run a comparison in the Methylomics → Differential Methylation tab first.")
        } else tagList(
          p(class = "submodule-desc", sprintf(
            "Live Methylomics session: %s (%s probes tested, %s significant).",
            methyl_results$dmp$comparison %||% "(unnamed)",
            format(methyl_results$dmp$n_probes %||% 0, big.mark = ","), format(methyl_results$dmp$n_sig %||% 0, big.mark = ","))),
          div(class = "empty-note", icon("triangle-exclamation"),
              "The live Methylomics session only shares summary statistics (comparison, probes tested, significant count) with other tabs - not the full per-CpG table, so it cannot be imported here for gene/CpG-level integration. Use \"Preloaded data\" (this project's own precomputed DMP tables) or \"Upload your data\" (export the Methylomics tab's own results table as CSV) instead.")
        )
      )
    })

    observeEvent(input$load_from_tx, {
      req(results, results$dge_runs, input$tx_run_pick)
      run <- results$dge_runs[[input$tx_run_pick]]
      validate(need(!is.null(run), "Selected Transcriptomics run is no longer available."))
      tbl <- run$table
      mapping <- c(gene = "gene", log2fc = "logFC", pvalue = NA_character_, fdr = "adj.P.Val", sample_id = NA_character_)
      std <- cx_standardize_expression(tbl, mapping)
      if (!std$ok) { showNotification(std$error, type = "error"); return() }
      raw$expr_df <- std$df
      raw$expr_source <- sprintf("Live Transcriptomics run: %s (%s)", run$contrast, run$method)
      raw$expr_sample_cols <- character(0)
      raw$expr_wide <- NULL
      showNotification(sprintf("Loaded %s genes from the Transcriptomics tab's \"%s\" run.", format(nrow(std$df), big.mark = ","), run$contrast), type = "message")
    })

    ## ---- Data Input: Upload ------------------------------------------

    observeEvent(input$expr_file, {
      res <- cx_read_and_detect(input$expr_file$datapath, input$expr_file$name, kind = "expression")
      if (!res$ok) { showNotification(res$error, type = "error"); expr_upload(NULL); return() }
      expr_upload(list(df = res$df, mapping = res$mapping, filename = input$expr_file$name))
    })
    observeEvent(input$meth_file, {
      res <- cx_read_and_detect(input$meth_file$datapath, input$meth_file$name, kind = "methylation")
      if (!res$ok) { showNotification(res$error, type = "error"); meth_upload(NULL); return() }
      meth_upload(list(df = res$df, mapping = res$mapping, filename = input$meth_file$name))
    })

    output$mapping_ui <- renderUI({
      if (!identical(input$mode, "upload")) {
        return(div(class = "empty-note", icon("circle-info"), "Column mapping applies to Upload mode - switch \"Input type\" above to Upload."))
      }
      if (is.null(expr_upload()) && is.null(meth_upload())) {
        return(div(class = "empty-note", icon("circle-info"), "Upload both files in step 1 to map their columns here."))
      }
      tagList(
        if (!is.null(expr_upload())) tagList(
          h5("Transcriptomics"),
          p(class = "submodule-desc", sprintf("%s - %s rows, %s columns.", expr_upload()$filename, format(nrow(expr_upload()$df), big.mark = ","), ncol(expr_upload()$df))),
          .cx_mapping_selects(ns, "map_expr_", expr_upload()$df, expr_upload()$mapping, CX_EXPRESSION_FIELDS)
        ),
        if (!is.null(meth_upload())) tagList(
          h5("Methylomics"),
          p(class = "submodule-desc", sprintf("%s - %s rows, %s columns.", meth_upload()$filename, format(nrow(meth_upload()$df), big.mark = ","), ncol(meth_upload()$df))),
          .cx_mapping_selects(ns, "map_meth_", meth_upload()$df, meth_upload()$mapping, CX_METHYLATION_FIELDS)
        ),
        actionButton(ns("confirm_mapping"), "Confirm mapping & standardize", icon = icon("check"), class = "btn-primary btn-sm")
      )
    })

    observeEvent(input$confirm_mapping, {
      eu <- expr_upload()
      if (!is.null(eu)) {
        mapping <- .cx_read_mapping(input, "map_expr_", CX_EXPRESSION_FIELDS)
        std <- cx_standardize_expression(eu$df, mapping)
        if (!std$ok) {
          showNotification(paste("Transcriptomics:", std$error), type = "error")
        } else {
          raw$expr_df <- std$df
          raw$expr_source <- sprintf("Uploaded: %s", eu$filename)
          raw$expr_mapping <- mapping
          raw$expr_wide <- eu$df
          raw$expr_sample_cols <- cx_detect_sample_columns(eu$df, mapping)
        }
      }
      mu <- meth_upload()
      if (!is.null(mu)) {
        mapping <- .cx_read_mapping(input, "map_meth_", CX_METHYLATION_FIELDS)
        std <- cx_standardize_methylation(mu$df, mapping)
        if (!std$ok) {
          showNotification(paste("Methylomics:", std$error), type = "error")
        } else {
          raw$meth_df <- std$df
          raw$meth_source <- sprintf("Uploaded: %s", mu$filename)
          raw$meth_mapping <- mapping
          raw$meth_wide <- mu$df
          raw$meth_sample_cols <- cx_detect_sample_columns(mu$df, mapping)
        }
      }
      if (!is.null(raw$expr_df) && !is.null(raw$meth_df)) {
        n_overlap <- length(intersect(raw$expr_df$gene, unique(raw$meth_df$gene)))
        showNotification(sprintf(
          "%s genes detected - %s methylation-associated genes detected - %s overlapping genes - integration ready.",
          format(nrow(raw$expr_df), big.mark = ","), format(length(unique(raw$meth_df$gene)), big.mark = ","),
          format(n_overlap, big.mark = ",")
        ), type = "message", duration = 8)
      }
    })

    ## ---- Status bar (spec section 30) ------------------------------------

    output$status_bar <- renderUI({
      tx_ok <- !is.null(raw$expr_df)
      mx_ok <- !is.null(raw$meth_df)
      map_ok <- tx_ok && mx_ok
      pairing <- integ$pairing
      sample_txt <- if (is.null(pairing)) "Not available" else if (isTRUE(pairing$paired)) sprintf("✓ (%d matched samples)", pairing$n_common) else "Not available (unpaired)"
      integ_ok <- !is.null(integ$df)
      item <- function(label, ok, extra = NULL) {
        span(class = if (isTRUE(ok)) "cx-status-ok" else "cx-status-pending",
             icon(if (isTRUE(ok)) "check" else "circle-notch"), " ", label, if (!is.null(extra)) paste0(" ", extra) else NULL)
      }
      tagList(
        div(class = "card", style = "padding: 10px 14px; margin-bottom: 10px; display:flex; gap:18px; flex-wrap:wrap; font-size: 0.92em;",
            item("Transcriptomics", tx_ok), item("Methylomics", mx_ok),
            item("Gene mapping", map_ok), item("Sample matching", isTRUE(pairing$paired), sample_txt),
            item("Integration", integ_ok)
        ),
        if (!is.null(raw$meth_unavailable_reason)) div(class = "empty-note", icon("triangle-exclamation"), raw$meth_unavailable_reason)
      )
    })

    ## ---- Run Integration --------------------------------------------------

    observeEvent(input$run_integration, {
      if (is.null(raw$expr_df)) { showNotification("Load or upload a Transcriptomics dataset first.", type = "error"); return() }
      if (is.null(raw$meth_df)) { showNotification("Load or upload a Methylomics dataset first.", type = "error"); return() }

      ## Gene identifier harmonization (spec section 7) - applied to both
      ## sides before the join, so e.g. an Ensembl ID on one side and its
      ## HGNC symbol on the other still match. Only exact/alias-resolved
      ## identifiers are rewritten; ambiguous/unmatched ones keep their
      ## original text (falling back to the previous exact-text join
      ## behavior for those specific genes rather than dropping them).
      id_harm <- cx_harmonize_gene_ids(c(raw$expr_df$gene, unique(raw$meth_df$gene)))
      integ$id_harmonization <- id_harm
      expr_h <- raw$expr_df
      meth_h <- raw$meth_df
      if (isTRUE(id_harm$ok)) {
        expr_h$gene <- cx_apply_harmonization(expr_h$gene, id_harm$df)
        expr_h <- cx_dedup_by_gene(expr_h)
        meth_h$gene <- cx_apply_harmonization(meth_h$gene, id_harm$df)
      }

      agg <- cx_aggregate_methylation(meth_h, method = input$agg_method, meth_thresh = input$meth_thresh, meth_fdr_thresh = input$meth_fdr_thresh)
      if (!agg$ok) { showNotification(agg$error, type = "error"); return() }
      integ$cpg_level <- agg$cpg_level

      expr_j <- data.frame(gene = expr_h$gene, log2fc = expr_h$log2fc,
                            expr_pvalue = expr_h$pvalue, expr_fdr = expr_h$fdr, stringsAsFactors = FALSE)
      meth_j <- agg$df
      colnames(meth_j)[colnames(meth_j) == "pvalue"] <- "meth_pvalue"
      colnames(meth_j)[colnames(meth_j) == "fdr"] <- "meth_fdr"
      joined <- merge(expr_j, meth_j, by = "gene", all = TRUE)

      if (all(is.na(joined$expr_fdr)) && any(!is.na(joined$expr_pvalue))) {
        joined$expr_fdr <- cx_adjust_p(joined$expr_pvalue, input$padj_method)
      }
      if (all(is.na(joined$meth_fdr)) && any(!is.na(joined$meth_pvalue))) {
        joined$meth_fdr <- cx_adjust_p(joined$meth_pvalue, input$padj_method)
      }

      classified <- cx_classify(joined, input$expr_thresh, input$expr_fdr_thresh, input$meth_thresh, input$meth_fdr_thresh)

      pairing <- cx_detect_sample_pairing(raw$expr_sample_cols, raw$meth_sample_cols)
      if (isTRUE(pairing$paired)) {
        expr_mat <- cx_build_gene_sample_matrix(raw$expr_wide, mapping_gene_col(raw$expr_mapping, raw$expr_wide), raw$expr_sample_cols)
        meth_mat <- cx_build_gene_sample_matrix(raw$meth_wide, mapping_gene_col(raw$meth_mapping, raw$meth_wide), raw$meth_sample_cols)
        cor_res <- cx_gene_correlation(expr_mat, meth_mat, pairing$common_samples, method = input$cor_method)
        if (isTRUE(cor_res$ok)) {
          cor_res$df$correlation_fdr <- cx_adjust_p(cor_res$df$p, input$padj_method)
          colnames(cor_res$df) <- c("gene", "correlation_r", "correlation_p", "correlation_n", "correlation_fdr")
          classified <- merge(classified, cor_res$df, by = "gene", all.x = TRUE)
        }
      }
      classified$evidence_level <- cx_classify_evidence(classified, has_correlation = isTRUE(pairing$paired))

      integ$df <- classified
      integ$pairing <- pairing
      integ$validation <- cx_validate_dataset(raw$expr_df, raw$meth_df, id_harm)
      integ_by_sex[[input$sex_stratum]] <- classified
      integ$universe <- classified$gene
      integ$run_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      integ$params <- list(
        sex_stratum = toupper(input$sex_stratum),
        input_mode = switch(input$mode, preloaded = sprintf("Preloaded (%s)", toupper(input$sex_stratum)), fromdataset = "From Dataset tab", fromapp = "Load from app", "Upload"),
        expr_source = raw$expr_source, meth_source = raw$meth_source,
        expr_thresh = input$expr_thresh, expr_fdr_thresh = input$expr_fdr_thresh,
        meth_thresh = input$meth_thresh, meth_fdr_thresh = input$meth_fdr_thresh,
        agg_method = input$agg_method, cor_method = if (isTRUE(pairing$paired)) input$cor_method else NULL,
        padj_method = input$padj_method,
        sample_matching = if (isTRUE(pairing$paired)) sprintf("Paired (%d common samples)", pairing$n_common) else "Not available (unpaired datasets)",
        gene_annotation_source = if (isTRUE(id_harm$ok)) "org.Hs.eg.db (Bioconductor) - exact ID/alias lookup only, no fuzzy matching" else "Not available - matched on exact provided text only",
        ## Keyed off the methylation source's own text, not `input$mode` -
        ## "From Dataset tab" can carry this exact same cx_load_default_
        ## methylation() output (its "Example data" mode calls the identical
        ## loader), so a mode-only check mislabeled that path as "uploaded".
        methylation_platform = if (grepl("bacon-adjusted", raw$meth_source %||% "", fixed = TRUE)) "Illumina 450K (IlluminaHumanMethylation450kanno.ilmn12.hg19)" else "Not specified by uploaded data",
        run_at = integ$run_at
      )
      integ$provenance <- cx_build_provenance(integ$params)
      active_filter("All")

      counts <- table(classified$category)
      cross_results$integration <- list(
        df = classified,
        summary = list(
          n_genes = nrow(classified),
          n_deg = sum(classified$sig_expression, na.rm = TRUE),
          n_dmg = sum(classified$sig_methylation, na.rm = TRUE),
          n_integrated = sum(classified$sig_expression & classified$sig_methylation, na.rm = TRUE),
          counts = as.list(counts)
        ),
        provenance = integ$provenance, params = integ$params, run_at = integ$run_at
      )
      showNotification("Integration complete.", type = "message")
    })

    ## `mapping` (a named char vector, possibly NULL for preloaded data which
    ## has no upload mapping) -> the gene-ID column name to use when building
    ## a genes x samples matrix from the original wide upload.
    mapping_gene_col <- function(mapping, wide_df) {
      if (!is.null(mapping) && !is.na(mapping[["gene"]])) return(mapping[["gene"]])
      "gene"
    }

    ## ---- Summary cards (spec section 27) ----------------------------------

    output$summary_cards <- renderUI({
      req(integ$df)
      df <- integ$df
      counts <- table(df$category)
      card <- function(label, value, key, color = "blue") {
        div(
          class = "card", style = "flex: 1 1 140px; cursor:pointer; text-align:center; padding:10px;",
          onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("card_click"), key),
          div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink),
              format(value, big.mark = ",")),
          div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label)
        )
      }
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          card("Genes analyzed", nrow(df), "All"),
          card("Significant DEGs", sum(df$sig_expression, na.rm = TRUE), "sig_expr_only", "aqua"),
          card("Significant DMGs", sum(df$sig_methylation, na.rm = TRUE), "sig_meth_only", "aqua"),
          card("Integrated genes", sum(df$sig_expression & df$sig_methylation, na.rm = TRUE), "sig_both", "violet"),
          card("Hyper + Down", counts[["Hyper + Down"]] %||% 0L, "Hyper + Down", "red"),
          card("Hypo + Up", counts[["Hypo + Up"]] %||% 0L, "Hypo + Up", "blue"),
          card("Hyper + Up", counts[["Hyper + Up"]] %||% 0L, "Hyper + Up", "orange"),
          card("Hypo + Down", counts[["Hypo + Down"]] %||% 0L, "Hypo + Down", "violet")
      )
    })
    observeEvent(input$card_click, active_filter(input$card_click))

    ## Advanced Filters (spec section 20) applied on top of the category
    ## quick-filter above - each is a cheap in-memory subset of the already-
    ## computed integ$df, so these are live/reactive rather than gated behind
    ## a separate "Apply Filters" button (unlike the genuinely expensive
    ## Pathway/Network analyses, which do stay behind their own Run buttons).
    filtered_df <- reactive({
      req(integ$df)
      df <- cx_filter_by_category(integ$df, active_filter())
      if (isTRUE(input$show_sig_only)) df <- df[df$sig_expression | df$sig_methylation, , drop = FALSE]
      df <- cx_filter_by_region(df, input$filter_region)
      df <- cx_filter_by_island(df, input$filter_island)
      df <- cx_filter_by_evidence(df, input$filter_evidence)
      df <- cx_filter_by_min_cpg(df, input$filter_min_cpg)
      df <- cx_filter_by_correlation_direction(df, input$filter_cor_direction)
      df
    })

    ## ---- Provenance --------------------------------------------------------

    output$provenance_ui <- renderUI({
      if (is.null(integ$provenance)) return(div(class = "empty-note", icon("circle-info"), "Run an integration to see its parameters here."))
      tags$ul(style = "padding-left: 18px; font-size: 0.85em;", lapply(integ$provenance, tags$li))
    })

    ## ---- Validation panel (spec section 6) --------------------------------

    output$validation_ui <- renderUI(cx_validation_checklist_ui(integ$validation))

    ## ---- CpG-Level table (Level 2, spec section 9) ------------------------

    output$cpg_level_ui <- renderUI({
      if (is.null(integ$cpg_level)) return(cx_empty_state())
      tagList(
        p(class = "submodule-desc", "Every individual CpG mapped to an analyzed gene - not collapsed to the gene level. Joined here with that gene's expression direction for context."),
        fluidRow(
          ## A plain text filter, not a selectizeInput with all ~19k gene
          ## names as choices - a client-side selectize list that large is a
          ## well-known browser-freezing anti-pattern (confirmed directly:
          ## an earlier version of this UI made screenshotting/interacting
          ## with this tab time out on this app's own real preloaded data).
          column(6, textInput(ns("cpg_level_gene_pick"), "Filter to gene(s)", placeholder = "Comma-separated symbols, e.g. TP53, BRCA1")),
          column(6, checkboxInput(ns("cpg_level_sig_only"), "Significant CpGs only", value = FALSE))
        ),
        DT::dataTableOutput(ns("cpg_level_table"))
      )
    })

    cpg_level_display <- reactive({
      req(integ$cpg_level)
      df <- integ$cpg_level
      if (!is.null(integ$df)) {
        exp_lookup <- stats::setNames(integ$df$expression_direction, integ$df$gene)
        df$expression_direction <- unname(exp_lookup[df$gene])
      }
      gene_filter <- trimws(strsplit(input$cpg_level_gene_pick %||% "", ",")[[1]])
      gene_filter <- gene_filter[nzchar(gene_filter)]
      if (length(gene_filter) > 0) df <- df[toupper(df$gene) %in% toupper(gene_filter), , drop = FALSE]
      if (isTRUE(input$cpg_level_sig_only)) df <- df[df$sig_cpg %in% TRUE, , drop = FALSE]
      df
    })

    ## Deliberately NOT suspendWhenHidden = FALSE, unlike results_table below -
    ## this table can have hundreds of thousands of rows (spec section 39:
    ## don't render/serialize large tables that aren't even being looked at),
    ## and nothing outside this tab reads its row-selection state, so there's
    ## no reason to force it to compute before the user opens the CpG-Level tab.
    output$cpg_level_table <- DT::renderDataTable({
      req(cpg_level_display())
      DT::datatable(cpg_level_display(), rownames = FALSE,
                     options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    }, server = TRUE)

    ## ---- Sex Comparison (spec section 14) ---------------------------------

    sex_comparison_result <- reactive({
      runs <- list(all = integ_by_sex$all, female = integ_by_sex$female, male = integ_by_sex$male)
      cx_compare_sexes(runs)
    })

    output$sex_comparison_ui <- renderUI({
      cmp <- sex_comparison_result()
      tagList(
        cx_sex_comparison_summary_ui(cmp),
        if (isTRUE(cmp$ok)) DT::dataTableOutput(ns("sex_comparison_table"))
      )
    })
    output$sex_comparison_table <- DT::renderDataTable({
      cmp <- sex_comparison_result()
      req(isTRUE(cmp$ok))
      DT::datatable(cmp$table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
    })

    ## ---- Results Table ------------------------------------------------------

    output$results_table_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      display_cols <- c("gene", "log2fc", "expr_fdr", "expression_direction",
                         "n_cpg_total", "n_cpg_significant", "n_cpg_hyper_sig", "n_cpg_hypo_sig",
                         "dbeta_mean", "dbeta_median", "primary_region", "region",
                         "n_island", "methylation_direction", "meth_fdr",
                         "correlation_r", "correlation_fdr", "category", "evidence_level", "cpg")
      display_cols <- intersect(display_cols, colnames(integ$df))
      tagList(
        fluidRow(
          column(6, selectInput(ns("filter_category"), "Filter by category", choices = CX_FILTER_CHOICES, selected = "All")),
          column(6, selectizeInput(ns("visible_cols"), "Visible columns", choices = display_cols, selected = display_cols, multiple = TRUE))
        ),
        div(class = "table-toolbar",
            downloadButton(ns("dl_table_csv"), "CSV", class = "btn-sm"),
            downloadButton(ns("dl_table_tsv"), "TSV", class = "btn-sm"),
            downloadButton(ns("dl_table_xlsx"), "XLSX", class = "btn-sm"),
            actionButton(ns("view_gene_btn"), "View selected gene", icon = icon("magnifying-glass"), class = "btn-sm")),
        DT::dataTableOutput(ns("results_table"))
      )
    })
    observeEvent(input$filter_category, active_filter(input$filter_category), ignoreInit = TRUE)

    table_df <- reactive({
      req(integ$df)
      df <- filtered_df()
      cols <- intersect(input$visible_cols %||% colnames(df), colnames(df))
      if (length(cols) == 0) cols <- colnames(df)
      df[, cols, drop = FALSE]
    })

    output$results_table <- DT::renderDataTable({
      req(table_df())
      DT::datatable(table_df(), rownames = FALSE, selection = "single",
                     options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    outputOptions(output, "results_table", suspendWhenHidden = FALSE)

    ## Individual CpG rows for one gene, for the detail modal (Level 2 data,
    ## never collapsed - spec section 17).
    gene_cpg_rows <- function(gene) {
      if (is.null(integ$cpg_level)) return(NULL)
      integ$cpg_level[integ$cpg_level$gene == gene, , drop = FALSE]
    }

    observeEvent(input$view_gene_btn, {
      sel <- input$results_table_rows_selected
      validate(need(length(sel) == 1, "Select exactly one row in the table first."))
      ## Indexes into filtered_df() (not the possibly column-subsetted
      ## table_df()) so the modal always has every field regardless of which
      ## "Visible columns" the user picked - row order/count is identical
      ## between the two since column selection never reorders/drops rows.
      row <- filtered_df()[sel, , drop = FALSE]
      showModal(cx_gene_detail_modal(row, integ$pairing, gene_cpg_rows(row$gene[1])))
    })

    ## ---- Quadrant Plot -------------------------------------------------------

    output$quadrant_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        fluidRow(
          ## Plain text, not a selectizeInput with every gene as a choice -
          ## see the CpG-Level tab's comment: a client-side choice list this
          ## large (~19k genes on this app's own real data) reproducibly
          ## hangs the browser.
          column(6, textInput(ns("quadrant_search"), "Search / label gene(s)", placeholder = "Comma-separated symbols, e.g. TP53, BRCA1")),
          column(6, div(style = "padding-top: 24px;",
                        downloadButton(ns("dl_quadrant_png"), "PNG", class = "btn-sm"),
                        downloadButton(ns("dl_quadrant_pdf"), "PDF", class = "btn-sm"),
                        downloadButton(ns("dl_quadrant_svg"), "SVG", class = "btn-sm")))
        ),
        p(class = "submodule-desc", "Each point is a gene. This shows a statistical association between methylation and expression change - it does not establish that one causes the other."),
        plotly::plotlyOutput(ns("quadrant_plot"), height = "520px")
      )
    })

    quadrant_source_df <- reactive({
      req(integ$df)
      df <- integ$df
      if (isTRUE(input$show_nonsig) == FALSE) df <- df[df$category != "Not significant", , drop = FALSE]
      if (isTRUE(input$show_sig_only)) df <- df[df$sig_expression | df$sig_methylation, , drop = FALSE]
      df[!is.na(df$dbeta) & !is.na(df$log2fc), , drop = FALSE]
    })

    quadrant_search_genes <- reactive({
      g <- trimws(strsplit(input$quadrant_search %||% "", ",")[[1]])
      g[nzchar(g)]
    })

    output$quadrant_plot <- plotly::renderPlotly({
      df <- quadrant_source_df()
      validate(need(nrow(df) > 0, "No genes to plot at the current filters/thresholds."))
      cx_quadrant_plotly(df, input$expr_thresh, input$meth_thresh, input$show_quadrant_lines,
                         highlight = quadrant_search_genes(), show_labels = isTRUE(input$show_labels), source_id = ns("quadrant_src"))
    })

    observeEvent(plotly::event_data("plotly_click", source = ns("quadrant_src")), {
      ed <- plotly::event_data("plotly_click", source = ns("quadrant_src"))
      req(ed$key)
      row <- integ$df[integ$df$gene == ed$key[1], , drop = FALSE]
      validate(need(nrow(row) >= 1, "Gene not found."))
      showModal(cx_gene_detail_modal(row[1, , drop = FALSE], integ$pairing, gene_cpg_rows(row$gene[1])))
    })

    ## ---- Heatmap ---------------------------------------------------------

    output$heatmap_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        fluidRow(
          column(4, selectInput(ns("heatmap_selection"), "Genes to show",
                                  choices = c("Top N by |log2FC|" = "topn", "All significant" = "sig", "By category" = "category"))),
          column(4, conditionalPanel(condition = sprintf("input['%s'] == 'topn'", ns("heatmap_selection")),
                                       numericInput(ns("heatmap_topn"), "N", value = 30, min = 5, max = 200))),
          column(4, conditionalPanel(condition = sprintf("input['%s'] == 'category'", ns("heatmap_selection")),
                                       selectInput(ns("heatmap_category"), "Category", choices = CX_CATEGORY_ORDER)))
        ),
        fluidRow(
          column(6, checkboxInput(ns("heatmap_cluster"), "Hierarchical clustering", value = TRUE)),
          column(6, radioButtons(ns("heatmap_scale"), "Scaling", choices = c("None" = "none", "Row" = "row", "Column" = "column"), inline = TRUE))
        ),
        downloadButton(ns("dl_heatmap_png"), "PNG", class = "btn-sm"),
        downloadButton(ns("dl_heatmap_pdf"), "PDF", class = "btn-sm"),
        plotOutput(ns("heatmap_plot"), height = "560px")
      )
    })

    heatmap_genes <- reactive({
      req(integ$df)
      df <- integ$df
      df <- switch(input$heatmap_selection %||% "topn",
        topn = df[order(-abs(df$log2fc %||% 0)), , drop = FALSE][seq_len(min(input$heatmap_topn %||% 30, nrow(df))), , drop = FALSE],
        sig = df[df$sig_expression & df$sig_methylation, , drop = FALSE],
        category = df[df$category == input$heatmap_category, , drop = FALSE]
      )
      df
    })

    cx_heatmap_matrix <- function() {
      df <- heatmap_genes()
      validate(need(nrow(df) >= 2, "Not enough genes selected to build a heatmap (need at least 2)."))
      m <- as.matrix(df[, c("log2fc", "dbeta")])
      rownames(m) <- df$gene
      colnames(m) <- c("Expression (log2FC)", "Methylation (Δβ)")
      m
    }

    output$heatmap_plot <- renderPlot({
      m <- cx_heatmap_matrix()
      ## silent = FALSE (the default) is required here - pheatmap only draws
      ## to the active device when silent is FALSE, and renderPlot's
      ## recording device needs that draw call to happen inside this block.
      pheatmap::pheatmap(m, cluster_rows = isTRUE(input$heatmap_cluster), cluster_cols = FALSE,
                          scale = input$heatmap_scale %||% "none",
                          color = grDevices::colorRampPalette(c(ARTHOMIX_COLORS$blue, "white", ARTHOMIX_COLORS$red))(100),
                          fontsize_row = if (nrow(m) > 60) 5 else 8)
    })

    ## ---- Volcano Plots -----------------------------------------------------

    output$volcano_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        checkboxInput(ns("highlight_hits"), "Highlight cross-omics hits", value = TRUE),
        fluidRow(
          column(6, downloadButton(ns("dl_volcano_expr"), "Download expression volcano (PNG)", class = "btn-sm")),
          column(6, downloadButton(ns("dl_volcano_meth"), "Download methylation volcano (PNG)", class = "btn-sm"))
        ),
        fluidRow(
          column(6, plotOutput(ns("volcano_expr"), height = "440px")),
          column(6, plotOutput(ns("volcano_meth"), height = "440px"))
        )
      )
    })
    output$volcano_expr <- renderPlot(cx_volcano_expr_plot(integ$df, input$expr_thresh, input$expr_fdr_thresh, isTRUE(input$highlight_hits)))
    output$volcano_meth <- renderPlot(cx_volcano_meth_plot(integ$df, input$meth_thresh, input$meth_fdr_thresh, isTRUE(input$highlight_hits)))

    ## ---- Correlation -------------------------------------------------------

    output$correlation_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      if (!isTRUE(integ$pairing$paired)) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
          "Unpaired datasets detected - gene-level differential-change integration is available above, but sample-level correlation is disabled unless matched samples are provided. Upload per-sample expression and methylation matrices that share at least 3 sample identifiers to enable this."))
      }
      tagList(
        p(class = "submodule-desc", sprintf("%d samples matched between the two uploaded datasets.", integ$pairing$n_common)),
        fluidRow(
          column(4, numericInput(ns("cor_min_r"), "Min |r|", value = 0, min = 0, max = 1, step = 0.05)),
          column(4, numericInput(ns("cor_max_fdr"), "Max FDR", value = 1, min = 0, max = 1, step = 0.05)),
          column(4, radioButtons(ns("cor_direction"), "Direction", choices = c("Any" = "any", "Positive" = "pos", "Negative" = "neg"), inline = TRUE))
        ),
        DT::dataTableOutput(ns("correlation_table")),
        tags$hr(),
        ## Plain text, not a selectizeInput with every gene as a choice - see
        ## the CpG-Level tab's comment for why.
        textInput(ns("cor_gene_pick"), "Gene-level correlation plot", placeholder = "Gene symbol, e.g. TP53"),
        plotOutput(ns("correlation_scatter"), height = "380px")
      )
    })

    correlation_filtered <- reactive({
      req(integ$df, "correlation_r" %in% colnames(integ$df))
      df <- integ$df[!is.na(integ$df$correlation_r), c("gene", "correlation_r", "correlation_p", "correlation_fdr"), drop = FALSE]
      df <- df[abs(df$correlation_r) >= (input$cor_min_r %||% 0) & df$correlation_fdr <= (input$cor_max_fdr %||% 1), , drop = FALSE]
      if (identical(input$cor_direction, "pos")) df <- df[df$correlation_r > 0, , drop = FALSE]
      if (identical(input$cor_direction, "neg")) df <- df[df$correlation_r < 0, , drop = FALSE]
      df[order(-abs(df$correlation_r)), , drop = FALSE]
    })
    output$correlation_table <- DT::renderDataTable({
      req(correlation_filtered())
      DT::datatable(correlation_filtered(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })

    output$correlation_scatter <- renderPlot({
      gene <- trimws(input$cor_gene_pick %||% "")
      req(nzchar(gene), integ$pairing$paired)
      cx_gene_correlation_plot(raw$expr_wide, raw$meth_wide, raw$expr_mapping, raw$meth_mapping,
                                gene, integ$pairing$common_samples, input$cor_method %||% "pearson")
    })

    ## ---- Overlap (Venn) ----------------------------------------------------

    output$overlap_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        fluidRow(column(6, plotOutput(ns("venn_plot"), height = "380px")),
                 column(6, plotOutput(ns("overlap_bars"), height = "380px")))
      )
    })
    output$venn_plot <- renderPlot({
      df <- integ$df
      sets <- list(DEGs = df$gene[df$sig_expression], DMGs = df$gene[df$sig_methylation])
      validate(need(length(sets$DEGs) > 0 || length(sets$DMGs) > 0, "No significant genes to compare at the current thresholds."))
      ggVennDiagram::ggVennDiagram(sets, label = "count") +
        ggplot2::scale_fill_gradient(low = "white", high = ARTHOMIX_COLORS$blue) + theme_arthomix() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
    })
    output$overlap_bars <- renderPlot({
      df <- integ$df
      counts <- data.frame(category = CX_CATEGORY_ORDER, n = as.integer(table(factor(df$category, levels = CX_CATEGORY_ORDER))))
      ggplot2::ggplot(counts, ggplot2::aes(x = stats::reorder(category, n), y = n, fill = category)) +
        ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = CX_CATEGORY_COLORS, guide = "none") +
        ggplot2::labs(x = NULL, y = "Genes", title = "Category counts") + theme_arthomix()
    })

    ## ---- Pathways ----------------------------------------------------------

    output$pathway_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        fluidRow(
          column(4, selectInput(ns("pathway_category"), "Gene set", choices = c("All significant (both)" = "sig_both", CX_CATEGORY_ORDER[CX_CATEGORY_ORDER != "Not significant"]))),
          column(4, selectInput(ns("pathway_ontology"), "Ontology", choices = c("GO Biological Process" = "BP", "GO Molecular Function" = "MF", "GO Cellular Component" = "CC", "KEGG" = "KEGG"))),
          column(4, div(style = "padding-top: 24px;", actionButton(ns("run_pathway"), "Run pathway analysis", icon = icon("diagram-project"), class = "btn-primary btn-sm")))
        ),
        p(class = "submodule-desc", "Reactome is not available in this deployment (ReactomePA is not installed); GO and KEGG use clusterProfiler against the current integrated gene list as the universe."),
        DT::dataTableOutput(ns("pathway_table")),
        plotOutput(ns("pathway_dotplot"), height = "460px")
      )
    })

    observeEvent(input$run_pathway, {
      df <- integ$df
      genes <- if (identical(input$pathway_category, "sig_both")) df$gene[df$sig_expression & df$sig_methylation] else df$gene[df$category == input$pathway_category]
      genes <- unique(stats::na.omit(genes))
      if (length(genes) < 3) { showNotification("Need at least 3 genes in the selected set to run pathway enrichment.", type = "error"); return() }
      withProgress(message = "Running pathway enrichment...", {
        res <- tryCatch(cx_run_pathway_enrichment(genes, integ$universe, input$pathway_ontology), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
        if (!isTRUE(res$ok)) { showNotification(res$error, type = "error"); pathway_result(NULL); return() }
        pathway_result(res)
      })
    })
    output$pathway_table <- DT::renderDataTable({
      req(pathway_result())
      DT::datatable(pathway_result()$table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    output$pathway_dotplot <- renderPlot({
      req(pathway_result())
      validate(need(!is.null(pathway_result()$table) && nrow(pathway_result()$table) > 0, "No enriched terms at the default significance cutoff."))
      enrichplot::dotplot(pathway_result()$enrich_obj, showCategory = 15) + theme_arthomix()
    })

    ## ---- Network (static, best-effort) -------------------------------------

    output$network_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      tagList(
        p(class = "submodule-desc", "A static gene-CpG association network (this deployment doesn't have an interactive network graphing package installed)."),
        fluidRow(
          column(6, selectInput(ns("network_category"), "Gene set", choices = c("All significant (both)" = "sig_both", CX_CATEGORY_ORDER[CX_CATEGORY_ORDER != "Not significant"]))),
          column(6, numericInput(ns("network_topn"), "Max genes shown", value = 15, min = 3, max = 40))
        ),
        plotOutput(ns("network_plot"), height = "520px")
      )
    })
    output$network_plot <- renderPlot({
      df <- integ$df
      sub <- if (identical(input$network_category, "sig_both")) df[df$sig_expression & df$sig_methylation, , drop = FALSE] else df[df$category == input$network_category, , drop = FALSE]
      sub <- sub[order(-abs(sub$log2fc %||% 0)), , drop = FALSE]
      sub <- sub[seq_len(min(input$network_topn %||% 15, nrow(sub))), , drop = FALSE]
      validate(need(nrow(sub) >= 2, "Not enough genes in this set to draw a network."))
      cx_gene_cpg_network_plot(sub)
    })

    ## ---- Genomic View --------------------------------------------------------

    output$genomic_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      if (!any(!is.na(integ$df$chr))) {
        return(div(class = "empty-note", icon("circle-info"), "No chromosome/position information is available for this dataset - genomic view is optional and skipped when coordinates aren't present."))
      }
      plotOutput(ns("genomic_plot"), height = "460px")
    })
    output$genomic_plot <- renderPlot(cx_genomic_view_plot(integ$df))

    ## ---- Downloads -----------------------------------------------------------

    output$downloads_ui <- renderUI({
      if (is.null(integ$df)) return(cx_empty_state())
      strata_available <- names(Filter(Negate(is.null), list(all = integ_by_sex$all, female = integ_by_sex$female, male = integ_by_sex$male)))
      tagList(
        h5("Gene-level results (current run)"),
        downloadButton(ns("dl_table_csv"), "Results (CSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_tsv"), "Results (TSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_xlsx"), "Results (XLSX)", class = "btn-sm"),
        tags$hr(),
        h5("CpG-level results (Level 2, current run)"),
        downloadButton(ns("dl_cpg_csv"), "CpG-level (CSV)", class = "btn-sm"),
        tags$hr(),
        h5("By analysis group (spec section 28)"),
        if (length(strata_available) == 0) div(class = "empty-note", icon("circle-info"), "Run at least one of ALL/FEMALE/MALE to enable per-group downloads.") else
          tagList(
            if ("all" %in% strata_available) downloadButton(ns("dl_all_stratum_csv"), "All-sample results (CSV)", class = "btn-sm"),
            if ("female" %in% strata_available) downloadButton(ns("dl_female_csv"), "Female results (CSV)", class = "btn-sm"),
            if ("male" %in% strata_available) downloadButton(ns("dl_male_csv"), "Male results (CSV)", class = "btn-sm")
          ),
        tags$hr(),
        h5("Report"),
        downloadButton(ns("dl_report"), "Analysis report (Markdown)", class = "btn-sm"),
        tags$hr(),
        h5("Everything"),
        downloadButton(ns("dl_all"), "Download All Results (.zip)", class = "btn-primary btn-sm")
      )
    })

    output$dl_table_csv <- downloadHandler(paste0("cross_omics_integration_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ$df, file, row.names = FALSE))
    output$dl_table_tsv <- downloadHandler(paste0("cross_omics_integration_", Sys.Date(), ".tsv"), function(file) utils::write.table(integ$df, file, row.names = FALSE, sep = "\t", quote = FALSE))
    output$dl_table_xlsx <- downloadHandler(paste0("cross_omics_integration_", Sys.Date(), ".xlsx"), function(file) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) { utils::write.csv(integ$df, file, row.names = FALSE); return() }
      openxlsx::write.xlsx(integ$df, file)
    })
    output$dl_cpg_csv <- downloadHandler(paste0("cross_omics_cpg_level_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ$cpg_level, file, row.names = FALSE))
    output$dl_all_stratum_csv <- downloadHandler(paste0("cross_omics_ALL_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ_by_sex$all, file, row.names = FALSE))
    output$dl_female_csv <- downloadHandler(paste0("cross_omics_FEMALE_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ_by_sex$female, file, row.names = FALSE))
    output$dl_male_csv <- downloadHandler(paste0("cross_omics_MALE_", Sys.Date(), ".csv"), function(file) utils::write.csv(integ_by_sex$male, file, row.names = FALSE))

    output$dl_quadrant_png <- downloadHandler(paste0("quadrant_plot_", Sys.Date(), ".png"), function(file) ggplot2::ggsave(file, cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6, dpi = 300))
    output$dl_quadrant_pdf <- downloadHandler(paste0("quadrant_plot_", Sys.Date(), ".pdf"), function(file) ggplot2::ggsave(file, cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6))
    output$dl_quadrant_svg <- downloadHandler(paste0("quadrant_plot_", Sys.Date(), ".svg"), function(file) ggplot2::ggsave(file, cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6))

    output$dl_heatmap_png <- downloadHandler(paste0("heatmap_", Sys.Date(), ".png"), function(file) { grDevices::png(file, width = 7, height = 7, units = "in", res = 300); pheatmap::pheatmap(cx_heatmap_matrix(), cluster_rows = isTRUE(input$heatmap_cluster), cluster_cols = FALSE, scale = input$heatmap_scale %||% "none"); grDevices::dev.off() })
    output$dl_heatmap_pdf <- downloadHandler(paste0("heatmap_", Sys.Date(), ".pdf"), function(file) { grDevices::pdf(file, width = 7, height = 7); pheatmap::pheatmap(cx_heatmap_matrix(), cluster_rows = isTRUE(input$heatmap_cluster), cluster_cols = FALSE, scale = input$heatmap_scale %||% "none"); grDevices::dev.off() })

    output$dl_volcano_expr <- downloadHandler(paste0("volcano_expression_", Sys.Date(), ".png"), function(file) ggplot2::ggsave(file, cx_volcano_expr_plot(integ$df, input$expr_thresh, input$expr_fdr_thresh, isTRUE(input$highlight_hits)), width = 7, height = 6, dpi = 300))
    output$dl_volcano_meth <- downloadHandler(paste0("volcano_methylation_", Sys.Date(), ".png"), function(file) ggplot2::ggsave(file, cx_volcano_meth_plot(integ$df, input$meth_thresh, input$meth_fdr_thresh, isTRUE(input$highlight_hits)), width = 7, height = 6, dpi = 300))

    output$dl_report <- downloadHandler(paste0("cross_omics_report_", Sys.Date(), ".md"), function(file) writeLines(cx_build_report(integ$df, integ$provenance), file))

    output$dl_all <- downloadHandler(
      filename = paste0("cross_omics_all_results_", Sys.Date(), ".zip"),
      content = function(file) {
        tmp <- tempfile(); dir.create(tmp)
        utils::write.csv(integ$df, file.path(tmp, "results.csv"), row.names = FALSE)
        writeLines(cx_build_report(integ$df, integ$provenance), file.path(tmp, "report.md"))
        tryCatch(ggplot2::ggsave(file.path(tmp, "quadrant_plot.png"), cx_quadrant_ggplot(quadrant_source_df(), input$expr_thresh, input$meth_thresh, input$show_quadrant_lines, quadrant_search_genes(), isTRUE(input$show_labels)), width = 7, height = 6, dpi = 300), error = function(e) NULL)
        old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
        setwd(tmp)
        utils::zip(file, files = list.files(tmp))
      }
    )
  })
}
