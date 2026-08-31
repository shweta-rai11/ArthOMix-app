## R/multiomics/mod_multi_biomarkercard.R
## Multi-Omics sub-module: Biomarker Card - integrated per-biomarker
## interpretation, mirroring the "Select Biomarker" -> "Biomarker Card"
## pattern of the Transcriptomics/Methylomics Biomarker Cards
## (mod_biomarkercard.R / mod_methyl_biomarkercard.R), but scoped to what is
## actually multiomics-specific: does a candidate gene–CpG pair have
## Transcriptomics-only, Methylomics-only, or Multiomics-supported evidence?
##
## Data source: exclusively multi_results$concordance$df (Gene–CpG
## Concordance, mod_multi_concordance.R) - the one table in this app that
## already joins a candidate gene's expression evidence to its CpG's
## methylation evidence, with DIABLO/SNF/Joint cross-omics support flags.
## Nothing here recomputes statistics or invents a new score; the evidence
## tier shown is `cx_classify_evidence()` (R/crossomics/
## crossomics_integration_helpers.R), the exact same classifier Cross-Omics
## Integration itself uses, applied read-only for display. Patient Evidence
## (per-sample expression/methylation values) is read directly from
## multi_dataset$layers when the active dataset is live (not preloaded);
## preloaded cohorts only have summary statistics, so that tab reports a
## short "not available" state rather than fabricating sample-level values.
##
## This module never writes to multi_results - it is read-only interpretation
## on top of what Gene–CpG Concordance (and, for DIABLO cross-reference,
## Biomarker Discovery) has already produced.

mod_multi_biomarkercard_config <- list(
  id = "biomarkercard", title = "Biomarker Card", icon = "id-card", group = "Interpretation",
  description = "Integrated gene–CpG biomarker interpretation: expression, methylation, and cross-omics evidence for one candidate at a time."
)

mod_multi_biomarkercard_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("mbc_subtabs"), type = "tabs",
      tabPanel("Select Biomarker", br(), uiOutput(ns("select_ui"))),
      tabPanel("Biomarker Card", br(), uiOutput(ns("card_ui")))
    )
  )
}

## ---------------------------------------------------------------------------
## Small display helpers - shared by every section below. Never silently
## print "NA"; every missing value reads "Not available".
## ---------------------------------------------------------------------------

mbc_val <- function(x) {
  if (length(x) == 0 || is.null(x) || (length(x) == 1 && is.na(x))) "Not available" else as.character(x)
}
mbc_num <- function(x, digits = 3) {
  if (length(x) == 0 || is.null(x) || (length(x) == 1 && is.na(x))) "Not available" else format(round(as.numeric(x), digits), nsmall = digits, big.mark = ",")
}
mbc_row <- function(label, value) tags$tr(tags$td(tags$b(label)), tags$td(value))

mbc_html_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  tags$table(
    class = "table table-condensed table-bordered",
    tags$thead(tags$tr(lapply(colnames(df), function(cn) tags$th(cn)))),
    tags$tbody(lapply(seq_len(nrow(df)), function(i) {
      tags$tr(lapply(df[i, , drop = TRUE], function(v) tags$td(if (length(v) == 0 || is.na(v)) "-" else as.character(v))))
    }))
  )
}

mbc_chip <- function(label, ok) {
  span(class = paste("pipeline-status-chip", if (isTRUE(ok)) "status-done" else "status-neutral"),
       icon(if (isTRUE(ok)) "circle-check" else "circle-minus"), label)
}

## Evidence tier -> chip color. "Strong"/"Moderate candidate" (both layers
## significant) get the green "done" treatment; single-layer-only evidence is
## amber ("pending" - real signal, but only from one omics layer); Discordant/
## Insufficient are neutral grey - never colored as if they were negative
## findings, since "no signal at these thresholds" is not an error state.
MBC_TIER_CLASS <- c(
  "Strong candidate" = "status-done", "Moderate candidate" = "status-done",
  "Expression-only" = "status-pending", "Methylation-only" = "status-pending",
  "Discordant" = "status-neutral", "Insufficient evidence" = "status-neutral"
)

## Short, direct reason why nothing is analyzed yet - never a long paragraph.
mbc_missing_data_note <- function(multi_dataset) {
  md <- multi_dataset
  if (!isTRUE(md$active %||% FALSE)) return("No active Multi-Omics dataset. Load one on the Dataset tab.")
  if (!identical(md$source, "preloaded")) {
    layers <- names(md$layers %||% list())
    if (!"Transcriptomics" %in% layers) return("Missing Transcriptomics data in the active dataset.")
    if (!"Methylomics" %in% layers) return("Missing Methylomics data in the active dataset.")
  }
  "Not analyzed yet. Run Gene–CpG Concordance (Biomarker modeling) first."
}

## ---- Server ----------------------------------------------------------------

mod_multi_biomarkercard_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Real Gene–CpG Concordance table, or NULL if it hasn't been run this
    ## session. Every read of multi_results/multi_dataset below happens
    ## inside a reactive/render/observe consumer, never at moduleServer
    ## setup time.
    conc_df <- reactive({
      conc <- multi_results$concordance
      if (is.null(conc) || is.null(conc$df) || nrow(conc$df) == 0) return(NULL)
      conc$df
    })

    ## Adds `evidence_tier` via the app's own cx_classify_evidence()
    ## (R/crossomics/crossomics_integration_helpers.R) - the same
    ## Transcriptomics-only / Methylomics-only / Moderate+Strong candidate /
    ## Discordant / Insufficient evidence classifier Cross-Omics Integration
    ## already uses. Read-only re-derivation for display; never stored back.
    evidence_df <- reactive({
      df <- conc_df()
      if (is.null(df)) return(NULL)
      has_cor <- "correlation_r" %in% colnames(df) && any(!is.na(df$correlation_r))
      df$evidence_tier <- as.character(cx_classify_evidence(df, has_correlation = has_cor))
      df
    })

    mbc_selected <- reactiveVal(NULL)  # list(gene_symbol=, cpg=)
    has_card <- reactiveVal(FALSE)

    ## Stale-selection guard: a new dataset, or a re-run of Gene–CpG
    ## Concordance, can change or remove the pair currently on the card -
    ## never leave a generated card pointing at results that no longer exist.
    observeEvent(multi_results$concordance, { mbc_selected(NULL); has_card(FALSE) }, ignoreInit = TRUE)
    observeEvent(multi_dataset$source, { mbc_selected(NULL); has_card(FALSE) }, ignoreInit = TRUE)
    observeEvent(input$mbc_search_mode, { mbc_selected(NULL) }, ignoreInit = TRUE)

    ## ==========================================================================
    ## Select Biomarker tab
    ## ==========================================================================

    output$select_ui <- renderUI({
      df <- evidence_df()
      if (is.null(df)) {
        return(tagList(
          multi_active_dataset_banner(multi_dataset),
          div(class = "empty-note", icon("circle-info"), tags$b("Not analyzed yet."), " ", mbc_missing_data_note(multi_dataset))
        ))
      }
      n_meaningful <- sum(df$evidence_tier != "Insufficient evidence", na.rm = TRUE)
      tagList(
        multi_active_dataset_banner(multi_dataset),
        if (n_meaningful == 0) div(class = "empty-note", icon("circle-info"), tags$b("No integrated biomarkers found."), " All gene–CpG pairs are below the significance thresholds used by Gene–CpG Concordance."),
        div(
          class = "card",
          div(class = "card-title", icon("magnifying-glass"), "Find a Biomarker"),
          radioButtons(ns("mbc_search_mode"), NULL, inline = TRUE,
                       choices = c("Browse Integrated Biomarkers" = "browse", "Type a gene or CpG" = "type")),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'browse'", ns("mbc_search_mode")),
            p(class = "submodule-desc", "Click a row to select that gene–CpG pair, then click \"Generate Biomarker Card\"."),
            DT::dataTableOutput(ns("mbc_browse_table"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'type'", ns("mbc_search_mode")),
            textInput(ns("mbc_text_input"), "Gene symbol or CpG ID", placeholder = "e.g. TNF or cg00000029"),
            uiOutput(ns("mbc_text_matches_ui"))
          )
        ),
        div(
          class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap;",
          uiOutput(ns("mbc_selection_status_ui"), inline = TRUE),
          actionButton(ns("mbc_generate_btn"), "Generate Biomarker Card", icon = icon("id-card"), class = "btn-primary btn-sm")
        )
      )
    })

    mbc_browse_disp <- reactive({
      df <- req(evidence_df())
      data.frame(
        Gene = df$gene_symbol, CpG = df$cpg, Evidence = df$evidence_tier,
        Direction = as.character(df$direction_classification),
        DIABLO = ifelse(df$diablo %in% TRUE, "Yes", "-"),
        SNF = ifelse(df$snf %in% TRUE, "Yes", "-"),
        `Priority score` = df$priority_score,
        stringsAsFactors = FALSE, check.names = FALSE
      )
    })

    output$mbc_browse_table <- DT::renderDataTable({
      DT::datatable(mbc_browse_disp(), rownames = FALSE, selection = "single",
                    options = list(pageLength = 10, scrollX = TRUE, order = list(list(6, "desc"))),
                    class = "stripe hover compact")
    })
    outputOptions(output, "mbc_browse_table", suspendWhenHidden = FALSE)

    observeEvent(input$mbc_browse_table_rows_selected, {
      df <- req(evidence_df())
      idx <- input$mbc_browse_table_rows_selected
      if (length(idx) == 1) mbc_selected(list(gene_symbol = df$gene_symbol[idx], cpg = df$cpg[idx]))
    })

    mbc_text_hits <- reactive({
      q <- trimws(input$mbc_text_input %||% "")
      df <- evidence_df()
      if (!nzchar(q) || is.null(df)) return(NULL)
      hit <- toupper(df$gene_symbol) == toupper(q) | toupper(df$cpg) == toupper(q)
      if (!any(hit)) return(NULL)
      df[hit, , drop = FALSE]
    })

    output$mbc_text_matches_ui <- renderUI({
      q <- trimws(input$mbc_text_input %||% "")
      if (!nzchar(q)) return(div(class = "empty-note", icon("circle-info"), "Type a gene symbol or CpG ID above."))
      sub <- mbc_text_hits()
      if (is.null(sub)) return(div(class = "empty-note", icon("circle-info"), sprintf("No match for \"%s\" in the current Gene–CpG Concordance results.", q)))
      tagList(
        p(class = "submodule-desc", sprintf("%d matching pair(s) - click a row to select.", nrow(sub))),
        DT::dataTableOutput(ns("mbc_text_table"))
      )
    })
    output$mbc_text_table <- DT::renderDataTable({
      sub <- req(mbc_text_hits())
      disp <- data.frame(Gene = sub$gene_symbol, CpG = sub$cpg, Evidence = sub$evidence_tier,
                         Direction = as.character(sub$direction_classification), stringsAsFactors = FALSE)
      DT::datatable(disp, rownames = FALSE, selection = "single", options = list(pageLength = 5, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "mbc_text_table", suspendWhenHidden = FALSE)
    observeEvent(input$mbc_text_table_rows_selected, {
      sub <- req(mbc_text_hits())
      idx <- input$mbc_text_table_rows_selected
      if (length(idx) == 1) mbc_selected(list(gene_symbol = sub$gene_symbol[idx], cpg = sub$cpg[idx]))
    })

    output$mbc_selection_status_ui <- renderUI({
      sel <- mbc_selected()
      if (is.null(sel)) tagList(icon("circle-info"), "Select a biomarker above before clicking Generate.")
      else tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("Selected: %s – %s", sel$gene_symbol, sel$cpg)), " - click Generate to build the card.")
    })

    observeEvent(input$mbc_generate_btn, { if (!is.null(mbc_selected())) has_card(TRUE) }, ignoreInit = TRUE)

    ## ==========================================================================
    ## Biomarker Card tab
    ## ==========================================================================

    sel_rows <- reactive({
      sel <- mbc_selected(); df <- evidence_df()
      if (is.null(sel) || is.null(df)) return(NULL)
      out <- df[df$gene_symbol == sel$gene_symbol & df$cpg == sel$cpg, , drop = FALSE]
      if (nrow(out) == 0) NULL else out
    })

    output$card_ui <- renderUI({
      if (!isTRUE(has_card())) return(div(class = "empty-note", icon("circle-info"), "No biomarker generated yet. Select one on \"Select Biomarker\", then click \"Generate Biomarker Card\"."))
      rows <- sel_rows()
      if (is.null(rows)) return(div(class = "empty-note", style = "border-color: var(--color-danger, #e34948);", icon("circle-xmark"), "This biomarker is no longer present in the current results. Go back and re-select one."))
      r <- rows[1, ]
      tagList(
        h4(icon("dna"), sprintf(" %s – %s", mbc_val(r$gene_symbol), mbc_val(r$cpg))),
        tabsetPanel(
          id = ns("mbc_card_subtabs"), type = "tabs",
          tabPanel("Overview", br(), mbc_section_overview(r)),
          tabPanel("Biomarker Status", br(), mbc_section_status(r)),
          tabPanel("Dataset", br(), tagList(multi_active_dataset_banner(multi_dataset), mbc_section_dataset(rows))),
          tabPanel("Expression", br(), mbc_section_expression(r)),
          tabPanel("Methylation", br(), mbc_section_methylation(r)),
          tabPanel("Integrated Evidence", br(), mbc_section_integrated(r)),
          tabPanel("Patient Evidence", br(), box(width = NULL, title = "Patient Evidence", status = "primary", solidHeader = FALSE, uiOutput(ns("mbc_patient_ui")))),
          tabPanel("Download", br(), box(width = NULL, title = "Download", status = "primary", solidHeader = FALSE,
                                          p(class = "submodule-desc", "Full annotated row(s) for this gene–CpG pair, as computed by Gene–CpG Concordance."),
                                          downloadButton(ns("mbc_dl"), "Download CSV", class = "btn-sm")))
        )
      )
    })

    output$mbc_dl <- downloadHandler(
      filename = function() "multiomics_biomarker_card.csv",
      content = function(file) utils::write.csv(req(sel_rows()), file, row.names = FALSE)
    )

    ## ---- Patient Evidence: real per-sample values from the live active
    ## dataset only - preloaded cohorts have no sample-level matrix to read,
    ## and that is disclosed rather than approximated. ----
    patient_data <- reactive({
      rows <- sel_rows()
      if (is.null(rows)) return(list(ok = FALSE, reason = "No biomarker selected."))
      r <- rows[1, ]
      md <- multi_dataset
      if (!isTRUE(md$active %||% FALSE) || identical(md$source, "preloaded"))
        return(list(ok = FALSE, reason = "Patient-level values are not available for preloaded cohort tables (summary statistics only)."))
      layers <- md$layers %||% list()
      if (!all(c("Transcriptomics", "Methylomics") %in% names(layers)))
        return(list(ok = FALSE, reason = "Both Transcriptomics and Methylomics layers are required for patient-level evidence."))
      expr_mat <- layers[["Transcriptomics"]]; meth_mat <- layers[["Methylomics"]]
      gene <- r$gene_symbol; cpg <- r$cpg
      if (is.na(gene) || is.na(cpg) || !(gene %in% colnames(expr_mat)) || !(cpg %in% colnames(meth_mat)))
        return(list(ok = FALSE, reason = sprintf("\"%s\" or \"%s\" is not present in the active dataset's layers.", mbc_val(gene), mbc_val(cpg))))
      match_res <- mcc_match_samples(expr_mat, meth_mat)
      if (!isTRUE(match_res$ok) || length(match_res$common_samples) < 3)
        return(list(ok = FALSE, reason = "Not enough matched samples between the expression and methylation layers."))
      ids <- match_res$common_samples
      tbl <- data.frame(Sample = ids, Expression = round(expr_mat[ids, gene], 3), Methylation = round(meth_mat[ids, cpg], 4), stringsAsFactors = FALSE)
      list(ok = TRUE, tbl = tbl)
    })

    output$mbc_patient_ui <- renderUI({
      pd <- patient_data()
      if (!isTRUE(pd$ok)) return(div(class = "empty-note", icon("circle-info"), pd$reason))
      tagList(
        p(class = "submodule-desc", sprintf("%d matched patient(s) with both expression and methylation values for this pair.", nrow(pd$tbl))),
        DT::dataTableOutput(ns("mbc_patient_table"))
      )
    })
    output$mbc_patient_table <- DT::renderDataTable({
      pd <- patient_data(); req(pd$ok)
      DT::datatable(pd$tbl, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "mbc_patient_table", suspendWhenHidden = FALSE)
  })
}

## ---------------------------------------------------------------------------
## Biomarker Card sub-tab sections - pure functions over one already-resolved
## concordance row (or, for Dataset, all rows for this pair, e.g. per-sex
## strata). No reactivity needed here; the card only rebuilds when the
## selection or the underlying results change (see card_ui above).
## ---------------------------------------------------------------------------

mbc_section_overview <- function(r) {
  box(width = NULL, title = "Identity", status = "primary", solidHeader = FALSE,
      tags$table(class = "table table-condensed",
                 mbc_row("Gene", mbc_val(r$gene_symbol)),
                 mbc_row("Ensembl ID", mbc_val(r$gene_id)),
                 mbc_row("CpG", mbc_val(r$cpg)),
                 mbc_row("Location", if (!is.na(r$chr %||% NA) && !is.na(r$pos %||% NA)) sprintf("%s:%s", r$chr, format(r$pos, big.mark = ",")) else "Not available"),
                 mbc_row("Region", mbc_val(r$region_fine)),
                 mbc_row("Island context", mbc_val(r$island_context)),
                 mbc_row("Dataset / cohort", mbc_val(r$dataset))
      ))
}

mbc_section_status <- function(r) {
  tier <- r$evidence_tier %||% "Insufficient evidence"
  tier_class <- MBC_TIER_CLASS[[tier]] %||% "status-neutral"
  desc <- CX_EVIDENCE_DESCRIPTIONS[[tier]] %||% ""
  tagList(
    box(width = NULL, title = "Evidence Summary", status = "primary", solidHeader = FALSE,
        span(class = paste("pipeline-status-chip", tier_class), icon("certificate"), sprintf("Evidence: %s", tier)),
        p(style = "margin-top:8px;", desc),
        div(style = "display:flex; gap:8px; flex-wrap:wrap; margin-top:8px;",
            mbc_chip("DIABLO", r$diablo %in% TRUE), mbc_chip("SNF", r$snf %in% TRUE), mbc_chip("Joint", r$joint %in% TRUE))
    ),
    box(width = NULL, title = "Biomarker Score", status = "primary", solidHeader = FALSE,
        tags$table(class = "table table-condensed",
                   mbc_row("Priority score", mbc_num(r$priority_score, 1)),
                   mbc_row("Label", mbc_val(r$evidence_label)),
                   mbc_row("DIABLO status", mbc_val(r$diablo_status)),
                   mbc_row("DIABLO loading", mbc_num(r$diablo_loading, 3))
        ))
  )
}

mbc_section_dataset <- function(rows) {
  r <- rows[1, ]
  strata <- if ("sex" %in% colnames(rows) && length(unique(rows$sex)) > 0) paste(unique(rows$sex), collapse = ", ") else "Pooled"
  tagList(
    box(width = NULL, title = "Cohort", status = "primary", solidHeader = FALSE,
        tags$table(class = "table table-condensed",
                   mbc_row("Dataset / cohort", mbc_val(r$dataset)),
                   mbc_row("Strata in this result", strata)
        ),
        if (nrow(rows) > 1) tagList(
          p(class = "submodule-desc", "This gene–CpG pair has more than one row below (e.g. per-sex strata):"),
          mbc_html_table(data.frame(
            Stratum = if ("sex" %in% colnames(rows)) rows$sex else rep("Pooled", nrow(rows)),
            Log2FC = round(rows$log2fc, 3), `Expr FDR` = signif(rows$expr_fdr, 3),
            dBeta = round(rows$dbeta, 3), `Meth FDR` = signif(rows$meth_fdr, 3),
            Direction = as.character(rows$direction_classification),
            check.names = FALSE
          ))
        )
    )
  )
}

mbc_section_expression <- function(r) {
  box(width = NULL, title = "Expression Evidence (Transcriptomics)", status = "primary", solidHeader = FALSE,
      tags$table(class = "table table-condensed",
                 mbc_row("Log2 fold-change", mbc_num(r$log2fc, 3)),
                 mbc_row("FDR", mbc_num(r$expr_fdr, 4)),
                 mbc_row("Direction", mbc_val(r$expression_direction)),
                 mbc_row("Significant", if (isTRUE(r$sig_expression)) "Yes" else "No")
      ))
}

mbc_section_methylation <- function(r) {
  box(width = NULL, title = "Methylation Evidence (Methylomics)", status = "primary", solidHeader = FALSE,
      tags$table(class = "table table-condensed",
                 mbc_row("Delta beta", mbc_num(r$dbeta, 3)),
                 mbc_row("FDR", mbc_num(r$meth_fdr, 4)),
                 mbc_row("Direction", mbc_val(r$methylation_direction)),
                 mbc_row("Significant", if (isTRUE(r$sig_methylation)) "Yes" else "No")
      ))
}

mbc_section_integrated <- function(r) {
  box(width = NULL, title = "Integrated Evidence", status = "primary", solidHeader = FALSE,
      tags$table(class = "table table-condensed",
                 mbc_row("Direction classification", mbc_val(as.character(r$direction_classification))),
                 mbc_row("Canonical rule", mbc_val(r$canonical_label)),
                 mbc_row("Sample correlation (r)", mbc_num(r$correlation_r, 3)),
                 mbc_row("Correlation p-value", mbc_num(r$correlation_p, 4)),
                 mbc_row("Matched samples (n)", mbc_val(r$correlation_n))
      ),
      p(class = "submodule-desc", MCC_CANONICAL_RULE_TEXT))
}
