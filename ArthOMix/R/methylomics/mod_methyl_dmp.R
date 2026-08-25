## R/methylomics/mod_methyl_dmp.R
## Methylomics sub-module: Differential Methylation (DMPs).
##
## Two sections, both gated on data availability rather than strictly on
## methyl_dataset$preloaded - switching the Dataset tab's source (a new
## "Load preloaded dataset" click vs. an upload) always updates which
## content each shows, but both can be visible together when the preloaded
## dataset's own live matrix is available in this deployment:
##
##   1. "Default analysis (GSE42861)" - rendered only when
##      methyl_dataset$preloaded is TRUE (req() at the top of default_ui).
##      Reproduces the thesis's own published sex-stratified limma +
##      bacon-correction model exactly, at its own documented threshold
##      (FDR<0.05, adjustable here), for BOTH pipeline stages: the plain
##      model (which found zero genome-wide-significant CpGs - a real,
##      documented finding, not a bug) and the SVA-adjusted follow-up that
##      resolves the residual inflation. Nothing here recomputes anything;
##      it filters/plots the existing tables. The stratum/FDR/Δβ/direction
##      controls are always visible, but the value boxes/volcano/table only
##      ever reflect whichever selection was in place the last time "View
##      results" was clicked - via an eventReactive (sva_run), not a plain
##      reactive - so nothing on screen updates just from touching a
##      control, and nothing renders at all before the first click.
##
##   2. "DMP Analysis" - rendered whenever a real beta/M-value matrix is
##      loaded (methyl_dataset$beta not NULL), whether that's an upload, or
##      (when this deployment has the raw ~2.1GB matrix copied in
##      alongside the app) the preloaded dataset's own live matrix. A full,
##      configurable EWAS-style engine (sex selection - including an
##      All-Samples/combined option alongside Female-only/Male-only -
##      dynamic group/reference/comparison, FDR/delta-beta/direction/
##      missingness/variance/SNP filters, optional covariates, limma model,
##      annotation, summary, volcano/Manhattan/top-DMP/beta-distribution
##      plots, reproducibility export) that reads whatever matrix/sample
##      sheet is currently loaded - nothing about it is hardcoded to
##      GSE42861's own group/sex labels or covariates, even when it's
##      GSE42861 supplying the matrix. When only the preloaded dataset's
##      metadata (not its live matrix) is available, this section shows an
##      explanatory message instead, and section 1 above is the only
##      analysis available.

mod_methyl_dmp_config <- list(
  id = "dmp", title = "Differential Methylation (DMPs)", icon = "chart-scatter", group = "Data",
  description = "Finds differentially methylated positions using limma with bacon correction. Uses the bundled whole-blood analysis by default, or a live, configurable run on your own dataset."
)

## Two subtabs, purely a UI reorganization - "SVA" and "DMP" wrap the exact
## same two uiOutputs (default_ui/live_ui) the module already had stacked
## with an <hr> between them; no output ID, reactive, or server logic
## changed. class="tx-menu-wrap" reuses the identical inner-tab look
## already used both by this same Methylomics module's own top-level
## Dataset/Sub-modules tabset (ui.R's mx_menu) and by Transcriptomics'
## Preprocessing tab (mod_preprocessing.R) - no new CSS needed. Shiny's
## tabsetPanel only suspends the hidden tab's *rendering*, not the
## underlying eventReactive results (sva_run()/live_result()), so
## switching subtabs never forces SVA or DMP to recompute.
mod_methyl_dmp_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("dmp_subtabs"), type = "tabs",
      tabPanel("SVA", br(), withSpinner(uiOutput(ns("default_ui")), color = "#2563EB", type = 6)),
      tabPanel("DMP", br(), withSpinner(uiOutput(ns("live_ui")), color = "#2563EB", type = 6))
    )
  )
}

## ---- Shared: filter + summarize + volcano, parameterized by which table --

mod_methyl_dmp_filter <- function(df, fdr_col, effect_col, fdr_max, effect_min, direction) {
  keep <- !is.na(df[[fdr_col]]) & df[[fdr_col]] <= fdr_max & abs(df[[effect_col]]) >= effect_min
  if (identical(direction, "hyper")) keep <- keep & df[[effect_col]] > 0
  if (identical(direction, "hypo"))  keep <- keep & df[[effect_col]] < 0
  df[keep, , drop = FALSE]
}

mod_methyl_dmp_volcano <- function(df, effect_col, p_col, effect_label, fdr_max, effect_min) {
  df$neglog10p <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  df$sig <- !is.na(df[[p_col]]) & df[[p_col]] <= fdr_max & abs(df[[effect_col]]) >= effect_min
  gg <- ggplot(df, aes(x = .data[[effect_col]], y = neglog10p, color = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = ARTHOMIX_STATUS$critical)) +
    labs(x = effect_label, y = expression(-log[10](p)), color = "Significant") +
    theme_arthomix()
  gg
}

## Genomic inflation factor (lambda_gc), standard median-chi-square formula -
## must be fed the FULL raw per-CpG p-value vector from the just-fitted model
## (before any FDR/effect-size/region-seeding filtering), or it is biased.
## The app's own bundled reference pipeline documents severe inflation on
## this exact dataset without SVA/bacon correction (lambda up to 12.1) but
## the live DMP/DMR engines below run plain limma with no such correction and
## previously surfaced no diagnostic - this is that diagnostic, shared by
## both live engines (mod_methyl_dmr.R calls this same function).
mod_methyl_lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  stats::median(stats::qchisq(1 - p, df = 1), na.rm = TRUE) / stats::qchisq(0.5, df = 1)
}

## Observed-vs-expected -log10(p) QQ plot for the same raw p-value vector -
## the standard companion visualization to lambda_gc; systematic departure
## above the y=x line indicates inflation beyond what chance alone predicts.
mod_methyl_qq_plot <- function(p) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  n <- length(p)
  if (n == 0) return(NULL)
  df <- data.frame(observed = -log10(p), expected = -log10(stats::ppoints(n)))
  ggplot(df, aes(x = expected, y = observed)) +
    geom_point(color = ARTHOMIX_COLORS$blue, alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = expression(-log[10] ~ "(expected p)"), y = expression(-log[10] ~ "(observed p)")) +
    theme_arthomix()
}

## ---- Live-engine helpers ---------------------------------------------------

## Row-chunked limma::lmFit(), reused by mod_methyl_dmr.R's live engine too.
##
## limma::lmFit() on a full genome-wide array (e.g. the preloaded dataset's
## real 412,492-probe x 689-sample matrix) can exceed R's vector memory
## limit even on a matrix with no missing values - stats::lm.fit(), which
## lmFit()/lm.series() dispatch to internally, needs the transposed data
## matrix plus its own residuals/effects/fitted-value buffers all at once,
## which for a matrix this size can genuinely need well over 10GB, on top
## of whatever this matrix's own copies already cost - a real crash
## ("vector memory limit of 16.0 Gb reached") reproduced against the
## preloaded dataset's actual matrix on a 16GB machine.
##
## lmFit()'s per-gene fit (coefficients/stdev.unscaled/sigma/df.residual/
## Amean) is independent per row given a shared design matrix, so fitting
## row-chunks separately and concatenating those (much smaller, n_genes x
## n_coef rather than n_genes x n_samples) per-gene summaries back together
## produces a combined MArrayLM object that downstream contrasts.fit()/
## eBayes()/topTable() cannot distinguish from one fit on the whole matrix
## at once - eBayes()'s empirical-Bayes moderation runs once, after
## combining, over every row's real summary statistics, so it sees the
## exact same inputs either way. Verified bit-for-bit identical topTable()
## output (logFC/t/P.Value/adj.P.Val, max abs diff 0) against a
## whole-matrix fit on synthetic data before this was adopted for real use.
methyl_chunked_lmfit <- function(m, design, chunk_size = 20000) {
  n <- nrow(m)
  if (n <= chunk_size) return(limma::lmFit(m, design))
  starts <- seq(1, n, by = chunk_size)
  fits <- lapply(starts, function(s) {
    e <- min(s + chunk_size - 1, n)
    limma::lmFit(m[s:e, , drop = FALSE], design)
  })
  per_row_1d <- c("df.residual", "sigma", "Amean")
  per_row_2d <- c("coefficients", "stdev.unscaled")
  out <- fits[[1]]
  for (nm in per_row_1d) if (!is.null(out[[nm]])) out[[nm]] <- do.call(c, lapply(fits, `[[`, nm))
  for (nm in per_row_2d) if (!is.null(out[[nm]])) out[[nm]] <- do.call(rbind, lapply(fits, `[[`, nm))
  out
}

## Auto-detects a sex/gender column by name (never assumed to exist) - same
## column-name candidates already used elsewhere in this file's folder
## (mod_methyl_qc.R's sex-check panel), so "which column is sex" stays one
## consistent guess across the whole Methylomics module.
mod_methyl_dmp_sex_col <- function(sheet) {
  if (is.null(sheet)) return(NULL)
  cols <- intersect(c("sex", "Sex", "SEX", "gender", "Gender"), colnames(sheet))
  if (length(cols) == 0) return(NULL)
  cols[1]
}

## Builds the "All samples / Female / Male" radio choices from whatever the
## sex column's actual values are - resolves to Female/Male labels only when
## exactly one level looks like each (F/Female, M/Male, case-insensitive);
## otherwise falls back to the raw values themselves rather than mislabeling
## an unfamiliar encoding, so this works for any dataset's own convention
## (spec: never hard-code male/female sample IDs or labels).
mod_methyl_dmp_sex_choices <- function(sheet, sex_col) {
  base <- c("All samples" = "__all__")
  if (is.null(sex_col) || !sex_col %in% colnames(sheet)) return(base)
  lvls <- unique(as.character(stats::na.omit(sheet[[sex_col]])))
  if (length(lvls) == 0) return(base)
  fem <- lvls[grepl("^f(emale)?$", lvls, ignore.case = TRUE)]
  mal <- lvls[grepl("^m(ale)?$", lvls, ignore.case = TRUE)]
  if (length(fem) == 1 && length(mal) == 1) {
    return(c(base, stats::setNames(fem, "Female"), stats::setNames(mal, "Male")))
  }
  c(base, stats::setNames(lvls, lvls))
}

## Candidate covariate columns: any column with at least 2 distinct
## non-missing values, excluding free-text/ID-like character columns (every
## value unique) - continuous covariates (age, cell-type fractions, PCs)
## are numeric and kept even when every value happens to be unique.
mod_methyl_dmp_covariate_cols <- function(sheet, exclude) {
  cols <- setdiff(colnames(sheet), exclude)
  keep <- vapply(cols, function(cl) {
    v <- sheet[[cl]]
    nu <- length(unique(stats::na.omit(v)))
    if (nu < 2) return(FALSE)
    if (!is.numeric(v) && nu == nrow(sheet)) return(FALSE)
    TRUE
  }, logical(1))
  cols[keep]
}

mod_methyl_dmp_manhattan <- function(df, fdr_max = 0.05) {
  d <- df[!is.na(df$chr) & !is.na(df$pos) & !is.na(df$fdr), , drop = FALSE]
  d$chr_num <- suppressWarnings(as.integer(gsub("^chr", "", d$chr, ignore.case = TRUE)))
  d <- d[!is.na(d$chr_num), , drop = FALSE]
  validate(need(nrow(d) > 0, "No tested CpG has both a chromosome and a genomic position from the available annotation - a Manhattan plot needs both."))
  d <- d[order(d$chr_num, d$pos), , drop = FALSE]
  d$idx <- seq_len(nrow(d))
  d$sig <- d$fdr <= fdr_max
  axis_df <- stats::aggregate(idx ~ chr_num, d, mean)
  ggplot(d, aes(x = idx, y = -log10(pmax(fdr, .Machine$double.xmin)))) +
    geom_point(aes(color = factor(chr_num %% 2)), size = 0.6, alpha = 0.5, show.legend = FALSE) +
    { if (any(d$sig)) geom_point(data = d[d$sig, , drop = FALSE], color = ARTHOMIX_STATUS$critical, size = 1.1) } +
    scale_color_manual(values = c("0" = ARTHOMIX_COLORS$blue, "1" = ARTHOMIX_COLORS$ink_muted)) +
    scale_x_continuous(breaks = axis_df$idx, labels = axis_df$chr_num, expand = expansion(mult = 0.01)) +
    geom_hline(yintercept = -log10(fdr_max), linetype = "dashed", color = "grey50") +
    labs(x = "Chromosome", y = expression(-log[10](FDR))) +
    theme_arthomix() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
}

mod_methyl_dmp_topplot <- function(df, rank_by = c("fdr", "dbeta", "combined"), n = 20) {
  rank_by <- match.arg(rank_by)
  d <- df[!is.na(df$fdr) & !is.na(df$dbeta), , drop = FALSE]
  validate(need(nrow(d) > 0, "No tested CpG has both an FDR and a Δβ value to rank."))
  ord <- switch(rank_by,
    fdr      = order(d$fdr, -abs(d$dbeta)),
    dbeta    = order(-abs(d$dbeta), d$fdr),
    combined = order(d$fdr / pmax(abs(d$dbeta), 1e-6))
  )
  d <- utils::head(d[ord, , drop = FALSE], n)
  label <- ifelse(!is.na(d$gene) & nzchar(d$gene), sprintf("%s (%s)", d$cpg, d$gene), d$cpg)
  d$label <- factor(label, levels = rev(label))
  gg <- ggplot(d, aes(x = dbeta, y = label, fill = dbeta > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = ARTHOMIX_STATUS$critical, "FALSE" = ARTHOMIX_COLORS$blue), guide = "none") +
    labs(x = "Δβ", y = NULL) +
    theme_arthomix()
  list(plot = gg, cpgs = as.character(d$cpg))
}

mod_methyl_dmp_betadist <- function(beta_mat, cpgs, grp) {
  cpgs <- intersect(cpgs, rownames(beta_mat))
  validate(need(length(cpgs) > 0, "None of the top-ranked CpGs are present in the analyzed beta matrix."))
  m <- beta_mat[cpgs, , drop = FALSE]
  long <- data.frame(
    cpg = rep(rownames(m), times = ncol(m)),
    beta = as.vector(m),
    group = rep(as.character(grp), each = nrow(m))
  )
  long$cpg <- factor(long$cpg, levels = cpgs)
  ggplot(long, aes(x = cpg, y = beta, fill = group)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = arthomix_pair(levels(grp))) +
    labs(x = NULL, y = "β value", fill = NULL) +
    theme_arthomix() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

mod_methyl_dmp_server <- function(id, methyl_dataset, methyl_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ================= 1. Default analysis (GSE42861) =====================

    default_data <- reactive({
      if (!METH_DATA_AVAILABLE) return(NULL)
      list(
        pheno = load_default_meth_pheno(),
        plain_f = load_default_dmp("plain", "female"), plain_m = load_default_dmp("plain", "male"),
        sva_f   = load_default_dmp("sva", "female"),   sva_m   = load_default_dmp("sva", "male")
      )
    })

    output$default_ui <- renderUI({
      ## Renders nothing at all (not even a placeholder) once the Dataset
      ## tab's upload path is selected - this section is exclusively the
      ## preloaded dataset's own reproduced analysis, not a generic
      ## "here's how to get results" prompt.
      req(isTRUE(methyl_dataset$preloaded))
      d <- default_data()
      req(d)
      n_ra <- sum(d$pheno$group == "RA", na.rm = TRUE); n_ctrl <- sum(d$pheno$group == "Control", na.rm = TRUE)
      n_f <- sum(d$pheno$sex == "F", na.rm = TRUE); n_m <- sum(d$pheno$sex == "M", na.rm = TRUE)
      n_sig_f <- sum(!is.na(d$sva_f$fdr_bacon) & d$sva_f$fdr_bacon < 0.05)
      n_sig_m <- sum(!is.na(d$sva_m$fdr_bacon) & d$sva_m$fdr_bacon < 0.05)
      n_sig_plain_f <- sum(!is.na(d$plain_f$fdr_bacon) & d$plain_f$fdr_bacon < 0.05)
      n_sig_plain_m <- sum(!is.na(d$plain_m$fdr_bacon) & d$plain_m$fdr_bacon < 0.05)
      tagList(
        div(class = "empty-note", icon("flask"),
            sprintf("Default analysis: preloaded whole-blood dataset, %d samples (%d RA / %d Control; %d female / %d male). Sex-stratified limma on M-values (~group + age + smoking + cell-type estimates), bacon-corrected, Benjamini-Hochberg FDR.",
                    nrow(d$pheno), n_ra, n_ctrl, n_f, n_m)),
        div(class = "card",
            div(class = "card-title", icon("check"), "SVA-adjusted model: the usable panel"),
            fluidRow(
              column(4,
                radioButtons(ns("sva_sex"), "Stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female"),
                numericInput(ns("sva_fdr"), "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                numericInput(ns("sva_dbeta"), "Min |Δβ|", value = 0, min = 0, max = 1, step = 0.01),
                radioButtons(ns("sva_direction"), "Direction", inline = TRUE,
                             choices = c("Any" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),
                actionButton(ns("sva_run_btn"), "View results", icon = icon("play"), class = "btn-primary btn-sm")
              ),
              column(8, withSpinner(uiOutput(ns("sva_volcano_ui")), color = "#2563EB", type = 6))
            ),
            uiOutput(ns("sva_valueboxes_ui")),
            uiOutput(ns("sva_table_ui"))
        ),
        div(class = "card",
            div(class = "card-title", icon("circle-info"), "Why an SVA-adjusted stage exists"),
            p(class = "submodule-desc",
              sprintf("The plain (unadjusted) sex-stratified model found %d female-stratum and %d male-stratum genome-wide-significant CpGs - both strata showed residual genomic inflation after bacon-correction, a real documented calibration problem, not a bug or a null biological result. The SVA-adjusted model above resolves it and is the panel actually used downstream (DMR calling, biomarker panels).",
                      n_sig_plain_f, n_sig_plain_m))
        )
      )
    })

    ## Results (value boxes, volcano, table) only ever reflect the stratum/
    ## filter selection at the moment "View results" was last clicked -
    ## changing the Stratum/FDR/Δβ/Direction controls afterward does NOT
    ## silently update anything on screen until the button is clicked again.
    ## eventReactive (not a plain reactive()) is what enforces that: it only
    ## re-evaluates on input$sva_run_btn, so every output below that reads
    ## from it inherits the same click-gating for free.
    sva_run <- eventReactive(input$sva_run_btn, {
      d <- default_data(); req(d)
      df <- if (identical(input$sva_sex, "male")) d$sva_m else d$sva_f
      req(df)
      list(
        df = df, sex = input$sva_sex,
        fdr = input$sva_fdr %||% 0.05, dbeta = input$sva_dbeta %||% 0,
        direction = if (identical(input$sva_direction, "any")) NULL else input$sva_direction
      )
    }, ignoreInit = TRUE)

    sva_has_run <- reactiveVal(FALSE)
    observeEvent(input$sva_run_btn, sva_has_run(TRUE), ignoreInit = TRUE)

    default_filtered <- reactive({
      r <- sva_run()
      mod_methyl_dmp_filter(r$df, "fdr_bacon", "dbeta", r$fdr, r$dbeta, r$direction)
    })

    output$sva_volcano_ui <- renderUI({
      if (!sva_has_run()) {
        return(p(class = "empty-note", icon("circle-info"), "Pick a stratum and filters on the left, then click \"View results\"."))
      }
      withSpinner(plotOutput(ns("default_volcano"), height = 320), color = "#2563EB", type = 6)
    })

    output$default_volcano <- renderPlot({
      r <- sva_run()
      mod_methyl_dmp_volcano(r$df, "dbeta", "fdr_bacon", "Δβ (RA - Control)", r$fdr, r$dbeta)
    })

    output$sva_valueboxes_ui <- renderUI({
      req(sva_has_run())
      r <- sva_run()
      n_sig <- sum(!is.na(r$df$fdr_bacon) & r$df$fdr_bacon < 0.05)
      sex_label <- if (identical(r$sex, "male")) "Male" else "Female"
      div(class = "methyl-stats-row", fluidRow(
        valueBox(format(n_sig), sprintf("%s FDR<0.05", sex_label),
                 icon = icon(if (identical(r$sex, "male")) "mars" else "venus"),
                 color = if (n_sig > 0) "green" else "light-blue", width = 4),
        valueBox(format(nrow(r$df), big.mark = ","), sprintf("CpGs tested (%s)", tolower(sex_label)), icon = icon("dna"), color = "purple", width = 4),
        valueBox(format(nrow(default_filtered()), big.mark = ","), "Passing current filters", icon = icon("filter"), color = "light-blue", width = 4)
      ))
    })

    output$sva_table_ui <- renderUI({
      req(sva_has_run())
      tagList(
        downloadButton(ns("download_default"), "Download filtered CSV", class = "btn-default btn-sm", style = "margin: 8px 0;"),
        DT::dataTableOutput(ns("default_table"))
      )
    })

    output$default_table <- DT::renderDataTable({
      DT::datatable(default_filtered(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("logFC_M", "dbeta", "t", "p_raw", "p_bacon", "fdr_bacon"), digits = 4)
    })
    ## Reproduced/verified bug, pre-existing (not caused by the SVA/DMP
    ## subtab split): this table's own output is nested inside two levels of
    ## dynamically-created renderUI (default_ui -> sva_table_ui). Shiny's
    ## client-side visibility tracking for suspendWhenHidden - which
    ## normally pauses computing/sending an output while its element looks
    ## hidden, and resumes once it's shown - was never correctly detecting
    ## this element as visible in this app's specific layout, so the table's
    ## rendered value was never sent to the browser at all (confirmed: the
    ## data/columns, e.g. "cpg", are correct server-side - DT::renderDataTable()
    ## itself is fine; only the client never receives it). Forcing this output
    ## to always compute/send, regardless of Shiny's own visibility
    ## detection, is what actually fixes it.
    outputOptions(output, "default_table", suspendWhenHidden = FALSE)

    output$download_default <- downloadHandler(
      filename = function() { r <- sva_run(); sprintf("gse42861_dmp_sva_%s_fdr%s.csv", r$sex, r$fdr) },
      content = function(file) utils::write.csv(default_filtered(), file, row.names = FALSE)
    )

    ## ================= 2. DMP Analysis (configurable live engine) =========
    ## Works against whatever methyl_dataset$beta currently holds - the
    ## preloaded dataset's real matrix when available, or an upload. Nothing
    ## in this section runs automatically: the interface (sex/group/filter/
    ## covariate controls) appears as soon as a matrix + sample sheet are
    ## loaded, but the model itself only fits when "Run DMP Analysis" is
    ## clicked (spec: never show results, plots, or empty result cards
    ## before that click).

    sex_col <- reactive({ mod_methyl_dmp_sex_col(methyl_dataset$sample_sheet) })

    id_cols <- reactive({
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet)) return(character(0))
      intersect(c("sample", "Sample", "sample_id", "Sample_ID", "gsm", "GSM"), colnames(sheet))
    })

    anno_result <- reactive({
      req(methyl_dataset$array_type)
      methyl_get_annotation(methyl_dataset$array_type)
    })

    dataset_norm_status <- reactive({
      req(methyl_dataset$beta)
      methyl_norm_status(methyl_dataset$beta, methyl_dataset, anno_result())
    })

    output$live_ui <- renderUI({
      ## Available whenever a real matrix is loaded - either from an
      ## upload, or (when this deployment has the raw ~2.1GB matrix copied
      ## in) the preloaded dataset's own live beta matrix - so an
      ## All-Samples (combined, both sexes) analysis, and single-sex live
      ## reruns, are available here in addition to the Female/Male panels
      ## above, which only ever reproduce the published pipeline's own
      ## precomputed, sex-stratified tables and never recompute anything.
      if (is.null(methyl_dataset$beta)) {
        msg <- if (isTRUE(methyl_dataset$preloaded))
          "The preloaded dataset's live beta matrix isn't available in this deployment, so only the sex-stratified reproduced analysis above is available (no live All-Samples/Female-only/Male-only option)."
        else
          "Upload a beta/M-value matrix or IDAT files on the Methylomics Dataset tab to configure and run a live differential methylation model - including an All-Samples (combined) option."
        return(div(class = "card",
          div(class = "card-title", icon("upload"), "DMP Analysis"),
          p(class = "submodule-desc", msg)
        ))
      }
      if (is.null(methyl_dataset$sample_sheet)) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "No sample sheet"),
          p(class = "submodule-desc", "A DMP model needs a group variable - re-load on the Dataset tab with a sample sheet/phenotype file included.")
        ))
      }
      sheet <- methyl_dataset$sample_sheet
      cols <- colnames(sheet)
      sc <- sex_col()
      ns_ <- dataset_norm_status()
      arr_ok <- isTRUE(anno_result()$ok)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("flask"), "DMP Analysis"),
            p(class = "empty-note", icon("circle-info"),
              sprintf("Dataset: %s. %s probes × %s samples. %s",
                      methyl_dataset$source %||% "(unnamed)",
                      format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta),
                      ns_$message)),
            fluidRow(
              column(4,
                tags$h5("Sex"),
                radioButtons(ns("live_sex"), NULL, inline = TRUE, choices = mod_methyl_dmp_sex_choices(sheet, sc), selected = "__all__"),
                if (is.null(sc)) p(class = "empty-note", icon("triangle-exclamation"), "No sex/gender column detected in the sample sheet - only \"All samples\" is available."),

                tags$h5("Comparison"),
                selectInput(ns("live_group_col"), "Group column", choices = cols,
                            selected = intersect(c("group", "Group", "disease", "Disease"), cols)[1] %||% cols[1]),
                uiOutput(ns("live_level_ui")),

                tags$h5("Filters"),
                numericInput(ns("live_fdr"), "Adjusted p-value (FDR) threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                numericInput(ns("live_dbeta"), "Absolute Δβ threshold", value = 0, min = 0, max = 1, step = 0.01),
                radioButtons(ns("live_direction"), "Direction", inline = TRUE,
                             choices = c("All DMPs" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),
                numericInput(ns("live_min_valid_pct"), "Minimum valid (non-missing) sample %", value = 80, min = 0, max = 100, step = 5),
                numericInput(ns("live_min_variance"), "Minimum methylation variance (optional)", value = 0, min = 0, step = 0.001),
                if (arr_ok) checkboxInput(ns("live_snp_filter"), "Remove SNP-associated probes (manifest Probe_rs/CpG_rs/SBE_rs)", value = FALSE)
                else p(class = "empty-note", icon("circle-info"), anno_result()$reason),

                tags$h5("Covariates (optional)"),
                uiOutput(ns("live_covariate_ui")),

                actionButton(ns("live_run_btn"), "Run DMP Analysis", icon = icon("play"), class = "btn-primary")
              ),
              column(8, withSpinner(uiOutput(ns("live_results_ui")), color = "#2563EB", type = 6))
            )
        )
      )
    })

    output$live_level_ui <- renderUI({
      req(input$live_group_col)
      sheet <- methyl_dataset$sample_sheet
      req(input$live_group_col %in% colnames(sheet))
      levels_available <- sort(unique(as.character(stats::na.omit(sheet[[input$live_group_col]]))))
      validate(need(length(levels_available) >= 2, "This column has fewer than two distinct values - pick a different group column."))
      tagList(
        selectInput(ns("live_ref"), "Reference group", choices = levels_available, selected = levels_available[1]),
        selectInput(ns("live_comp"), "Comparison group", choices = levels_available, selected = levels_available[min(2, length(levels_available))])
      )
    })

    output$live_covariate_ui <- renderUI({
      sheet <- methyl_dataset$sample_sheet
      req(sheet, input$live_group_col)
      exclude <- c(id_cols(), input$live_group_col)
      ## Sex is constant within a Female-only or Male-only run, so it can't
      ## be a covariate there; it's only offered when "All samples" is
      ## selected (spec section 8).
      if (!identical(input$live_sex %||% "__all__", "__all__") && !is.null(sex_col())) exclude <- c(exclude, sex_col())
      cand <- mod_methyl_dmp_covariate_cols(sheet, exclude)
      if (length(cand) == 0) {
        return(p(class = "empty-note", icon("circle-info"), "No additional phenotype columns are available to use as covariates."))
      }
      checkboxGroupInput(ns("live_covariates"), NULL, choices = cand, selected = character(0))
    })

    ## Reset whenever the underlying dataset changes identity (a new Load
    ## click on the Dataset tab, preloaded or uploaded) so a stale result
    ## from a previous data source is never left showing.
    live_has_run <- reactiveVal(FALSE)
    observeEvent(methyl_dataset$beta, live_has_run(FALSE), ignoreNULL = TRUE)

    live_result <- eventReactive(input$live_run_btn, withProgress(message = "Running DMP Analysis - a full genome-wide array (400k+ CpGs) can take a minute or more...", value = 0.2, {
      validate(need(!is.null(methyl_dataset$beta), "Load a dataset first."))
      sheet <- methyl_dataset$sample_sheet
      validate(need(!is.null(sheet), "No sample sheet loaded."))
      validate(need(!is.null(input$live_group_col) && input$live_group_col %in% colnames(sheet), "Pick a group column."))
      validate(need(!is.null(input$live_ref) && !is.null(input$live_comp), "Pick a reference and a comparison group."))
      validate(need(input$live_ref != input$live_comp, "Reference and comparison groups must be different."))

      sample_ids <- methyl_sheet_sample_ids(sheet, colnames(methyl_dataset$beta))
      common <- intersect(colnames(methyl_dataset$beta), sample_ids)
      validate(need(length(common) >= 6, "Fewer than 6 samples match between the matrix and the sample sheet."))
      beta0 <- methyl_dataset$beta[, common, drop = FALSE]
      ## The preloaded dataset's sample sheet is a data.table (pheno.rds) -
      ## single-bracket column selection by a character vector (used below
      ## for the covariate columns) means something different there
      ## (j evaluated as an expression, not a column-name vector) than for
      ## an uploaded sheet's plain data.frame. Coercing once here keeps
      ## every ph0/ph1 access below in ordinary data.frame semantics
      ## regardless of which path the sheet came from.
      ph0 <- as.data.frame(sheet)[match(common, sample_ids), , drop = FALSE]

      ## ---- sex subset --------------------------------------------------
      sex_choice <- input$live_sex %||% "__all__"
      sc <- sex_col()
      sex_label <- "All samples"
      if (!identical(sex_choice, "__all__")) {
        validate(need(!is.null(sc), "No sex column available to subset on."))
        keep_sex <- !is.na(ph0[[sc]]) & as.character(ph0[[sc]]) == sex_choice
        validate(need(sum(keep_sex) >= 6, sprintf("Fewer than 6 samples remain after restricting to sex = \"%s\".", sex_choice)))
        beta0 <- beta0[, keep_sex, drop = FALSE]
        ph0 <- ph0[keep_sex, , drop = FALSE]
        sex_label <- sex_choice
      }
      n_after_sex <- ncol(beta0)

      ## ---- group/reference/comparison subset ----------------------------
      grp_raw <- as.character(ph0[[input$live_group_col]])
      keep_grp <- !is.na(grp_raw) & grp_raw %in% c(input$live_ref, input$live_comp)
      beta1 <- beta0[, keep_grp, drop = FALSE]
      ph1 <- ph0[keep_grp, , drop = FALSE]
      grp <- factor(grp_raw[keep_grp], levels = c(input$live_ref, input$live_comp))
      ## beta0 (still full probe count, up to full sample count) has no
      ## further use once beta1 exists - on a full ~400k-probe array this
      ## is a ~2GB+ matrix left alive for no reason otherwise, one of
      ## several redundant full-size copies that made the live engine
      ## exceed R's vector memory limit on the real preloaded matrix (see
      ## the beta/M-value section below).
      rm(beta0)

      ## ---- optional covariates (complete cases only) --------------------
      covariate_cols <- input$live_covariates %||% character(0)
      cov_df <- NULL
      if (length(covariate_cols) > 0) {
        cc <- ph1[, covariate_cols, drop = FALSE]
        complete <- stats::complete.cases(cc)
        validate(need(sum(complete) >= 6, "Fewer than 6 samples have no missing values in the selected covariates. Deselect a covariate, or pick different ones."))
        beta1 <- beta1[, complete, drop = FALSE]
        ph1 <- ph1[complete, , drop = FALSE]
        grp <- grp[complete]
        cov_df <- ph1[, covariate_cols, drop = FALSE]
        for (cl in covariate_cols) if (!is.numeric(cov_df[[cl]])) cov_df[[cl]] <- factor(as.character(cov_df[[cl]]))
      }

      n_ref <- sum(grp == input$live_ref, na.rm = TRUE)
      n_comp <- sum(grp == input$live_comp, na.rm = TRUE)
      validate(need(n_ref >= 3 && n_comp >= 3,
        sprintf("Each group needs at least 3 samples to fit a model (reference \"%s\": %d, comparison \"%s\": %d).",
                input$live_ref, n_ref, input$live_comp, n_comp)))

      ## ---- probe filters (missingness / variance / SNP), BEFORE the M-value
      ## transform ------------------------------------------------------------
      ## Filtering runs on beta1 itself (already the 0-1 beta scale whenever
      ## input_scale=="beta" - the common case, including the preloaded
      ## dataset's own matrix - so no extra full-size copy is needed there;
      ## only the rarer M-value-upload path needs one temporary beta-scale
      ## conversion of the full matrix first). `m`/`beta_scale` are then
      ## derived only from the already-filtered subset below. On a full
      ## ~400k-probe array, computing the M-value transform AND the beta-
      ## scale conversion on the FULL matrix before filtering (the previous
      ## order) held 3+ full-size copies in memory at once and exceeded R's
      ## vector memory limit on the real preloaded matrix (689 samples) -
      ## this reordering keeps at most one full-size matrix alive at a time.
      is_m_scale <- identical(methyl_dataset$input_scale, "m")
      beta_scale_full <- if (is_m_scale) 2^beta1 / (1 + 2^beta1) else beta1

      max_na_frac <- 1 - (input$live_min_valid_pct %||% 80) / 100
      f_miss <- methyl_filter_missing(beta_scale_full, max_na_frac)
      keep_probe <- f_miss$keep
      f_var <- methyl_filter_variance(beta_scale_full, input$live_min_variance %||% 0)
      keep_probe <- keep_probe & f_var$keep
      snp_note <- NULL
      if (isTRUE(input$live_snp_filter)) {
        ar <- anno_result()
        if (isTRUE(ar$ok)) {
          f_snp <- methyl_filter_snp(beta_scale_full, ar)
          keep_probe <- keep_probe & f_snp$keep
          snp_note <- f_snp$note
        } else snp_note <- ar$reason
      }
      validate(need(sum(keep_probe) >= 10,
        "Fewer than 10 CpGs remain after the missingness/variance/SNP filters. Relax the filters and try again."))

      beta1 <- beta1[keep_probe, , drop = FALSE]
      if (is_m_scale) rm(beta_scale_full)  ## a real extra copy only in this branch; the else-branch aliases beta1, nothing to free
      gc(FALSE)  ## proactively reclaim the freed full-size copies before the memory-heavy limma fit below

      ## ---- beta/M-value matrices (filtered subset only) ---------------------
      m <- if (is_m_scale) beta1 else log2(pmin(pmax(beta1, 1e-6), 1 - 1e-6) / (1 - pmin(pmax(beta1, 1e-6), 1 - 1e-6)))
      beta_scale <- if (is_m_scale) 2^m / (1 + 2^m) else beta1

      ## ---- design + limma fit ---------------------------------------------
      design_grp <- stats::model.matrix(~0 + grp)
      colnames(design_grp) <- levels(grp)
      cov_names <- if (!is.null(cov_df)) colnames(cov_df) else character(0)
      if (length(cov_names) > 0) {
        safe <- sprintf("`%s`", cov_names)
        cov_form <- stats::as.formula(paste("~", paste(safe, collapse = " + ")))
        design_cov <- tryCatch(stats::model.matrix(cov_form, data = cov_df),
                                error = function(e) validate(need(FALSE, paste("Could not build a design matrix for the selected covariates:", conditionMessage(e)))))
        design_cov <- design_cov[, setdiff(colnames(design_cov), "(Intercept)"), drop = FALSE]
        design <- cbind(design_grp, design_cov)
      } else {
        design <- design_grp
      }
      validate(need(qr(design)$rank == ncol(design),
        "The selected group/covariate combination produces a rank-deficient design (e.g. a covariate that perfectly predicts the group, or is constant). Remove a covariate or change the comparison."))

      fit <- methyl_chunked_lmfit(m, design)
      cm <- tryCatch(limma::makeContrasts(contrasts = paste0(input$live_comp, "-", input$live_ref), levels = design),
                      error = function(e) validate(need(FALSE, "Could not build the reference/comparison contrast - group names may contain characters limma can't use directly (try renaming the group levels in your sample sheet).")))
      fit2 <- tryCatch(limma::eBayes(limma::contrasts.fit(fit, cm)),
                        error = function(e) validate(need(FALSE, paste("limma could not fit this model:", conditionMessage(e)))))
      tt <- limma::topTable(fit2, number = Inf, sort.by = "P")

      beta_ref <- rowMeans(beta_scale[, grp == input$live_ref, drop = FALSE], na.rm = TRUE)
      beta_comp <- rowMeans(beta_scale[, grp == input$live_comp, drop = FALSE], na.rm = TRUE)
      dbeta <- (beta_comp - beta_ref)[rownames(tt)]

      df <- data.frame(
        cpg = rownames(tt), t = tt$t, p_raw = tt$P.Value, fdr = tt$adj.P.Val,
        dbeta = dbeta, ref_mean_beta = beta_ref[rownames(tt)], comp_mean_beta = beta_comp[rownames(tt)],
        row.names = NULL, stringsAsFactors = FALSE
      )
      ar <- anno_result()
      df$chr <- NA_character_; df$pos <- NA_real_; df$gene <- NA_character_
      if (isTRUE(ar$ok)) {
        a <- ar$anno
        hit <- df$cpg %in% rownames(a)
        df$chr[hit] <- a[df$cpg[hit], "chr"]
        df$pos[hit] <- a[df$cpg[hit], "pos"]
        df$gene[hit] <- a[df$cpg[hit], "gene"]
      }
      df$direction <- ifelse(df$dbeta > 0, "hyper", "hypo")

      design_formula <- sprintf("Methylation ~ %s%s", input$live_group_col,
                                 if (length(cov_names) > 0) paste0(" + ", paste(cov_names, collapse = " + ")) else "")

      list(
        df = df, ref = input$live_ref, comp = input$live_comp,
        group_col = input$live_group_col, sex_label = sex_label, sex_col = sc,
        covariates = cov_names, design_formula = design_formula,
        n_after_sex = n_after_sex, n_ref = n_ref, n_comp = n_comp,
        n_probes_tested = nrow(df), n_probes_before_filter = nrow(beta1),
        min_valid_pct = input$live_min_valid_pct %||% 80, min_variance = input$live_min_variance %||% 0,
        snp_filter = isTRUE(input$live_snp_filter), snp_note = snp_note,
        missing_note = f_miss$note, variance_note = f_var$note,
        anno_ok = isTRUE(ar$ok), norm_status = dataset_norm_status(),
        dataset_source = methyl_dataset$source %||% "(unnamed)", preloaded = isTRUE(methyl_dataset$preloaded),
        beta_scale = beta_scale, grp = grp,
        lambda_gc = mod_methyl_lambda_gc(df$p_raw),
        run_at = Sys.time()
      )
    }), ignoreInit = TRUE)

    observeEvent(live_result(), {
      r <- live_result()
      methyl_results$dmp <- list(comparison = sprintf("%s vs %s (%s)", r$comp, r$ref, r$sex_label),
                                  n_probes = r$n_probes_tested,
                                  n_sig = sum(!is.na(r$df$fdr) & r$df$fdr < 0.05, na.rm = TRUE))
      live_has_run(TRUE)
    })

    ## FDR/Δβ/direction stay live-adjustable after a run without needing to
    ## refit the model - only the sex/group/covariate/probe-filter controls
    ## (which change what was actually fit) require clicking Run again.
    live_filtered <- reactive({
      r <- live_result()
      direction <- if (identical(input$live_direction, "any")) NULL else input$live_direction
      mod_methyl_dmp_filter(r$df, "fdr", "dbeta", input$live_fdr %||% 0.05, input$live_dbeta %||% 0, direction)
    })

    live_sig_count <- reactive({
      r <- live_result()
      sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05) & abs(r$df$dbeta) >= (input$live_dbeta %||% 0), na.rm = TRUE)
    })

    output$live_results_ui <- renderUI({
      req(live_has_run(), live_result())
      r <- live_result()
      n_sig <- live_sig_count()
      n_hyper <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05) & r$df$dbeta >= (input$live_dbeta %||% 0), na.rm = TRUE)
      n_hypo  <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05) & -r$df$dbeta >= (input$live_dbeta %||% 0), na.rm = TRUE)
      n_pass_fdr <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05), na.rm = TRUE)
      n_pass_dbeta <- sum(!is.na(r$df$dbeta) & abs(r$df$dbeta) >= (input$live_dbeta %||% 0), na.rm = TRUE)

      tagList(
        div(class = "card",
            div(class = "card-title", icon("circle-check"), "Configuration & sample sizes"),
            p(strong("Model: "), tags$code(r$design_formula), " (limma moderated t-test, eBayes)."),
            p(strong("Sex: "), r$sex_label, " | ", strong("Reference: "), r$ref, sprintf(" (n=%d)", r$n_ref),
              " | ", strong("Comparison: "), r$comp, sprintf(" (n=%d)", r$n_comp),
              " | ", strong("Total analyzed: "), r$n_ref + r$n_comp),
            if (length(r$covariates) > 0) p(strong("Covariates: "), paste(r$covariates, collapse = ", ")) else p(class = "submodule-desc", "No covariates selected."),
            p(class = "submodule-desc", r$norm_status$message),
            if ((r$n_ref < 10 || r$n_comp < 10))
              p(class = "empty-note", icon("triangle-exclamation"), "One or both groups have fewer than 10 samples - results may be underpowered."),
            p(class = "submodule-desc", sprintf("%s Probes tested after QC filters: %s of %s.",
                                                  paste(c(r$missing_note, r$variance_note, r$snp_note), collapse = " "),
                                                  format(r$n_probes_tested, big.mark = ","), format(r$n_probes_before_filter, big.mark = ",")))
        ),
        div(class = "card",
            div(class = "card-title", icon("gauge-high"), "Genomic inflation diagnostic"),
            p(strong("Genomic inflation factor (λ): "), if (is.na(r$lambda_gc)) "not available" else sprintf("%.2f", r$lambda_gc),
              " - the ratio of observed to expected median test statistic across all tested CpGs; 1.0 indicates no inflation."),
            if (!is.na(r$lambda_gc) && r$lambda_gc > 1.1)
              p(class = "empty-note", icon("triangle-exclamation"),
                "Genomic inflation is elevated (λ > 1.1). This live engine does not apply SVA/bacon correction - the app's own bundled reference pipeline found severe inflation on unrelated confounding in this dataset without it. Treat the p-values/FDR above as potentially optimistic and cross-check against the default/SVA-corrected panel before drawing conclusions."),
            withSpinner(plotOutput(ns("live_qq"), height = 320), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Summary"),
            div(class = "methyl-stats-row",
              fluidRow(
                valueBox(format(r$n_probes_tested, big.mark = ","), "CpGs tested", icon = icon("dna"), color = "purple", width = 3),
                valueBox(format(n_sig, big.mark = ","), "Significant (FDR + Δβ)", icon = icon("star"), color = if (n_sig > 0) "green" else "light-blue", width = 3),
                valueBox(format(n_hyper, big.mark = ","), "Hypermethylated", icon = icon("arrow-up"), color = "red", width = 3),
                valueBox(format(n_hypo, big.mark = ","), "Hypomethylated", icon = icon("arrow-down"), color = "blue", width = 3)
              ),
              fluidRow(
                valueBox(format(n_pass_fdr, big.mark = ","), sprintf("Pass FDR < %s", input$live_fdr %||% 0.05), icon = icon("filter"), color = "light-blue", width = 4),
                valueBox(format(n_pass_dbeta, big.mark = ","), sprintf("Pass |Δβ| ≥ %s", input$live_dbeta %||% 0), icon = icon("filter"), color = "light-blue", width = 4),
                valueBox(r$sex_label, "Sex subset", icon = icon("venus-mars"), color = "black", width = 4)
              )
            ),
            if (n_sig == 0) p(class = "empty-note", icon("circle-info"),
              "No DMPs passed the selected FDR and Δβ thresholds. Consider relaxing the filtering thresholds or reviewing the sample/group configuration.")
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-scatter"), "Volcano plot"),
            withSpinner(plotOutput(ns("live_volcano"), height = 340), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-column"), "Manhattan plot"),
            if (r$anno_ok) withSpinner(plotOutput(ns("live_manhattan"), height = 320), color = "#2563EB", type = 6)
            else p(class = "empty-note", icon("circle-info"), "No chromosome/position annotation is available for this array type in this deployment.")
        ),
        div(class = "card",
            div(class = "card-title", icon("ranking-star"), "Top DMPs"),
            fluidRow(
              column(6, selectInput(ns("live_rank_by"), "Rank by", choices = c("FDR" = "fdr", "Absolute Δβ" = "dbeta", "Combined (FDR/|Δβ|)" = "combined"), selected = "fdr")),
              column(6, selectInput(ns("live_top_n"), "Show top", choices = c(10, 20, 50, 100), selected = 20))
            ),
            withSpinner(plotOutput(ns("live_topplot"), height = 380), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "β-value distribution (top-ranked CpGs)"),
            withSpinner(plotOutput(ns("live_betadist"), height = 340), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("table"), "Results table"),
            div(class = "table-toolbar",
                downloadButton(ns("download_live"), "Download filtered CSV", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_config"), "Download analysis configuration", class = "btn-default btn-sm")),
            DT::dataTableOutput(ns("live_table"))
        )
      )
    })

    output$live_qq <- renderPlot({
      r <- live_result()
      mod_methyl_qq_plot(r$df$p_raw)
    })

    output$live_volcano <- renderPlot({
      r <- live_result()
      mod_methyl_dmp_volcano(r$df, "dbeta", "fdr", "Δβ", input$live_fdr %||% 0.05, input$live_dbeta %||% 0)
    })

    output$live_manhattan <- renderPlot({
      r <- live_result()
      req(r$anno_ok)
      mod_methyl_dmp_manhattan(r$df, input$live_fdr %||% 0.05)
    })

    live_top <- reactive({
      r <- live_result()
      mod_methyl_dmp_topplot(r$df, rank_by = input$live_rank_by %||% "fdr", n = as.integer(input$live_top_n %||% 20))
    })

    output$live_topplot <- renderPlot({ live_top()$plot })

    output$live_betadist <- renderPlot({
      r <- live_result()
      mod_methyl_dmp_betadist(r$beta_scale, live_top()$cpgs, r$grp)
    })

    output$live_table <- DT::renderDataTable({
      df <- live_filtered()
      show_cols <- c("cpg", "gene", "chr", "pos", "ref_mean_beta", "comp_mean_beta", "dbeta", "p_raw", "fdr", "direction")
      df$significant <- ifelse(!is.na(df$fdr) & df$fdr <= (input$live_fdr %||% 0.05) & abs(df$dbeta) >= (input$live_dbeta %||% 0), "Yes", "No")
      DT::datatable(df[, c(show_cols, "significant")], rownames = FALSE, filter = "top",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("ref_mean_beta", "comp_mean_beta", "dbeta", "p_raw", "fdr"), digits = 4)
    })
    outputOptions(output, "live_table", suspendWhenHidden = FALSE)  ## see default_table's own comment above

    output$download_live <- downloadHandler(
      filename = function() { r <- live_result(); sprintf("dmp_%s_vs_%s_%s.csv", r$comp, r$ref, r$sex_label) },
      content = function(file) utils::write.csv(live_filtered(), file, row.names = FALSE)
    )

    output$download_live_config <- downloadHandler(
      filename = function() "dmp_analysis_configuration.csv",
      content = function(file) {
        r <- live_result()
        cfg <- data.frame(
          parameter = c("dataset", "sex", "reference_group", "comparison_group", "n_reference", "n_comparison",
                        "n_total_analyzed", "covariates", "statistical_method", "design_formula",
                        "fdr_threshold", "dbeta_threshold", "direction_filter", "min_valid_sample_pct",
                        "min_variance", "snp_filter_applied", "n_cpgs_tested", "n_significant", "run_at"),
          value = c(r$dataset_source, r$sex_label, r$ref, r$comp, r$n_ref, r$n_comp, r$n_ref + r$n_comp,
                    if (length(r$covariates) > 0) paste(r$covariates, collapse = ";") else "(none)",
                    "limma (moderated t-test, eBayes)", r$design_formula,
                    input$live_fdr %||% 0.05, input$live_dbeta %||% 0, input$live_direction %||% "any",
                    r$min_valid_pct, r$min_variance, r$snp_filter, r$n_probes_tested, live_sig_count(),
                    format(r$run_at)),
          stringsAsFactors = FALSE
        )
        utils::write.csv(cfg, file, row.names = FALSE)
      }
    )
  })
}
