## R/methylomics/mod_methyl_dmr.R
## Methylomics sub-module: Differentially Methylated Regions (DMRs).
##
## Two sections, gated on data availability (can both be visible together
## when the preloaded dataset's live matrix is available):
##
##   1. "Default analysis (GSE42861)" - shown when methyl_dataset$preloaded
##      is TRUE. Reproduces the published sex-stratified DMRcate region
##      calling (script04_dmr_sexstratified/04_dmr_sexstratified.R:
##      lambda=1000, C=2, on the SVA-adjusted, bacon-corrected per-CpG
##      statistics from the DMP tab's default analysis) at its documented
##      region-level FDR (BH on the Stouffer statistic, <0.05, adjustable).
##      Filters/plots existing tables; nothing is recomputed. Gated behind
##      "Run DMR Analysis" via an eventReactive (mirrors mod_methyl_dmp.R's
##      sva_run).
##
##   2. "DMR Analysis" - shown whenever a real beta/M-value matrix is
##      loaded, whether uploaded or (if this deployment has the raw matrix
##      bundled) the preloaded dataset's own live matrix. A full,
##      configurable region-level engine: the DMP tab's live sex/group/
##      covariate/QC-filter machinery feeds a per-CpG limma fit into
##      DMRcate::dmrcate() (same algorithm the preloaded pipeline uses,
##      with lambda/C/min.cpgs/seeding-p exposed as controls), including an
##      All-Samples option alongside Female-only/Male-only. Shows an
##      explanatory message instead when only metadata (not the live
##      matrix) is available.

mod_methyl_dmr_config <- list(
  id = "dmr", title = "Differentially Methylated Regions (DMRs)", icon = "map-location-dot", group = "Data",
  description = "Performs differentially methylated region analysis."
)

## "SVA" and "DMR" subtabs wrap the same two uiOutputs (default_ui/live_ui).
## Shiny's tabsetPanel only suspends the hidden tab's rendering, not the
## underlying eventReactive results, so switching subtabs never forces the
## default or live DMR analysis to recompute.
mod_methyl_dmr_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("dmr_subtabs"), type = "tabs",
      tabPanel("SVA", br(), withSpinner(uiOutput(ns("default_ui")), color = "#2563EB", type = 6)),
      tabPanel("DMR", br(), withSpinner(uiOutput(ns("live_ui")), color = "#2563EB", type = 6))
    )
  )
}

## ---- Shared helpers --------------------------------------------------------

## Region-table filter: FDR/Δβ/direction via mod_methyl_dmp_filter()
## (reused unchanged from the DMP tab), plus region-specific dimensions
## (CpG count, width) DMPs don't have.
mod_methyl_dmr_filter <- function(df, fdr_col, effect_col, fdr_max, effect_min, direction,
                                   min_cpgs = 0, min_width = 0, max_width = Inf) {
  df <- mod_methyl_dmp_filter(df, fdr_col, effect_col, fdr_max, effect_min, direction)
  if ("no.cpgs" %in% colnames(df)) df <- df[df$no.cpgs >= min_cpgs, , drop = FALSE]
  if ("width" %in% colnames(df)) df <- df[df$width >= min_width & df$width <= max_width, , drop = FALSE]
  df
}

## Region-level counterpart to mod_methyl_dmp_manhattan() - reads
## seqnames/start/dmr_fdr instead of chr/pos/fdr.
mod_methyl_dmr_manhattan <- function(dt, fdr_max = 0.05) {
  d <- dt[!is.na(dt$seqnames) & !is.na(dt$start) & !is.na(dt$dmr_fdr), , drop = FALSE]
  d$chr_num <- suppressWarnings(as.integer(gsub("^chr", "", as.character(d$seqnames), ignore.case = TRUE)))
  d <- d[!is.na(d$chr_num), , drop = FALSE]
  validate(need(nrow(d) > 0, "No candidate region has both an autosomal chromosome and a region-level FDR - a genomic plot needs both."))
  d <- d[order(d$chr_num, d$start), , drop = FALSE]
  d$idx <- seq_len(nrow(d))
  d$sig <- d$dmr_fdr <= fdr_max
  axis_df <- stats::aggregate(idx ~ chr_num, d, mean)
  ggplot(d, aes(x = idx, y = -log10(pmax(dmr_fdr, .Machine$double.xmin)))) +
    geom_point(aes(color = factor(chr_num %% 2)), size = 0.8, alpha = 0.5, show.legend = FALSE) +
    { if (any(d$sig)) geom_point(data = d[d$sig, , drop = FALSE], color = ARTHOMIX_STATUS$critical, size = 1.4) } +
    scale_color_manual(values = c("0" = ARTHOMIX_COLORS$blue, "1" = ARTHOMIX_COLORS$ink_muted)) +
    scale_x_continuous(breaks = axis_df$idx, labels = axis_df$chr_num, expand = expansion(mult = 0.01)) +
    geom_hline(yintercept = -log10(fdr_max), linetype = "dashed", color = "grey50") +
    labs(x = "Chromosome", y = expression(-log[10]~"(region-level FDR)")) +
    theme_arthomix() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
}

## Top-N regions by region-level FDR, as a Δβ bar chart labeled by
## coordinates (+ overlapping gene when annotated).
mod_methyl_dmr_topplot <- function(dt, n = 20) {
  d <- dt[!is.na(dt$dmr_fdr) & !is.na(dt$meandiff), , drop = FALSE]
  validate(need(nrow(d) > 0, "No candidate region has both a region-level FDR and a mean Δβ value to rank."))
  d <- utils::head(d[order(d$dmr_fdr), , drop = FALSE], n)
  has_gene <- !is.na(d$overlapping.genes) & nzchar(d$overlapping.genes)
  label <- ifelse(has_gene,
                   sprintf("%s:%s-%s (%s)", d$seqnames, format(d$start, big.mark = ","), format(d$end, big.mark = ","), d$overlapping.genes),
                   sprintf("%s:%s-%s", d$seqnames, format(d$start, big.mark = ","), format(d$end, big.mark = ",")))
  d$label <- factor(label, levels = rev(label))
  ggplot(d, aes(x = meandiff, y = label, fill = meandiff > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = ARTHOMIX_STATUS$critical, "FALSE" = ARTHOMIX_COLORS$blue), guide = "none") +
    labs(x = "Mean Δβ", y = NULL) +
    theme_arthomix()
}

## Region x sample mean-methylation heatmap for significant DMRs (capped at
## `max_regions` by region-level FDR - a rendering cap only, table/exports
## still carry the full set). Each region's per-sample value is the mean β
## across probes falling inside its coordinates, via GenomicRanges::
## findOverlaps() (same idiom script04_dmr_sexstratified.R uses).
mod_methyl_dmr_heatmap <- function(sig_dt, beta_scale, probe_chr, probe_pos, grp, max_regions = 50) {
  validate(need(nrow(sig_dt) > 0, "No DMR passes the current FDR/Δβ thresholds - there is nothing to show in the heatmap."))
  sig_dt <- sig_dt[order(sig_dt$dmr_fdr), , drop = FALSE]
  sig_dt <- utils::head(sig_dt, max_regions)
  region_ids <- if ("dmr_id" %in% colnames(sig_dt)) sig_dt$dmr_id else sprintf("DMR%03d", seq_len(nrow(sig_dt)))
  region_gr <- GenomicRanges::GRanges(sig_dt$seqnames, IRanges::IRanges(sig_dt$start, sig_dt$end))
  probe_gr <- GenomicRanges::GRanges(probe_chr, IRanges::IRanges(probe_pos, probe_pos))
  hits <- GenomicRanges::findOverlaps(probe_gr, region_gr)
  qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
  validate(need(length(qh) > 0, "None of the significant regions' constituent CpGs could be matched back to the analyzed probe set."))
  mat <- matrix(NA_real_, nrow = nrow(sig_dt), ncol = ncol(beta_scale), dimnames = list(region_ids, colnames(beta_scale)))
  for (i in seq_len(nrow(sig_dt))) {
    probes_i <- qh[sh == i]
    if (length(probes_i) > 0) mat[i, ] <- colMeans(beta_scale[probes_i, , drop = FALSE], na.rm = TRUE)
  }
  sample_order <- colnames(mat)[order(grp)]
  long <- data.frame(
    region = rep(rownames(mat), times = ncol(mat)),
    sample = factor(rep(colnames(mat), each = nrow(mat)), levels = sample_order),
    beta = as.vector(mat),
    group = factor(rep(as.character(grp), each = nrow(mat)), levels = levels(grp)),
    stringsAsFactors = FALSE
  )
  long$region <- factor(long$region, levels = rev(rownames(mat)))
  ggplot(long, aes(x = sample, y = region, fill = beta)) +
    geom_tile() +
    facet_grid(~group, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = ARTHOMIX_COLORS$blue, mid = "white", high = ARTHOMIX_COLORS$red, midpoint = 0.5,
                          na.value = "grey90", limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "Mean β") +
    theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

## Human-readable sex label matching mod_methyl_dmp_sex_choices()'s choice
## names, so the DMR tab's banners say the same thing the DMP tab would.
mod_methyl_dmr_sex_label <- function(sheet, sex_col, value) {
  ch <- mod_methyl_dmp_sex_choices(sheet, sex_col)
  nm <- names(ch)[ch == value]
  if (length(nm) == 0) return(value)
  nm[1]
}

## ChAMPdata::probe.features chromosome/position lookup, used by the
## default section's region-level CpG inspection to resolve which CpGs
## from the DMP tab's per-CpG table fall inside a selected DMR's
## coordinates. Cached for the life of the R process.
.methyl_dmr_probe_pos_cache <- new.env(parent = emptyenv())
methyl_champ_probe_positions <- function() {
  if (!is.null(.methyl_dmr_probe_pos_cache$pos)) return(.methyl_dmr_probe_pos_cache$pos)
  if (!requireNamespace("ChAMPdata", quietly = TRUE)) return(NULL)
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ utils::data("probe.features", package = "ChAMPdata", envir = e); TRUE }, error = function(e) FALSE)
  if (!ok || is.null(e$probe.features)) return(NULL)
  pf <- e$probe.features
  pos <- data.frame(cpg = rownames(pf), chr = paste0("chr", pf$CHR), pos = pf$MAPINFO, stringsAsFactors = FALSE)
  .methyl_dmr_probe_pos_cache$pos <- pos
  pos
}

## Extra packages the live region-calling engine needs beyond limma -
## checked once so live_ui can degrade to a clear message instead of an
## opaque error if a deployment is missing one.
methyl_dmr_engine_pkgs_ok <- function() {
  requireNamespace("DMRcate", quietly = TRUE) &&
    requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("IRanges", quietly = TRUE) &&
    requireNamespace("S4Vectors", quietly = TRUE)
}

mod_methyl_dmr_server <- function(id, methyl_dataset, methyl_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ================= 1. Default analysis (GSE42861) =====================

    default_data <- reactive({
      if (!METH_DATA_AVAILABLE) return(NULL)
      list(
        pheno = load_default_meth_pheno(),
        dmr_f = load_default_dmr("female"), dmr_m = load_default_dmr("male")
      )
    })

    output$default_ui <- renderUI({
      ## Only for the preloaded dataset's own reproduced region calling.
      req(isTRUE(methyl_dataset$preloaded))
      d <- default_data()
      req(d)
      validate(need(!is.null(d$dmr_f) || !is.null(d$dmr_m), "No precomputed DMR results are available in this deployment."))
      n_f <- if (!is.null(d$pheno)) sum(d$pheno$sex == "F", na.rm = TRUE) else NA
      n_m <- if (!is.null(d$pheno)) sum(d$pheno$sex == "M", na.rm = TRUE) else NA
      tagList(
        div(class = "empty-note", icon("map-location-dot"),
            sprintf("Default analysis: DMRcate region calling (lambda=1000, C=2, hg19) on the SVA-adjusted, bacon-corrected per-CpG statistics from the Differentially Methylation (DMPs) tab's own default analysis, sex-stratified%s. Comparison: RA vs Control (the dataset's only case/control contrast).",
                    if (!is.na(n_f)) sprintf(" (%d female / %d male)", n_f, n_m) else "")),
        div(class = "card",
            div(class = "card-title", icon("map-location-dot"), "DMR configuration"),
            fluidRow(
              column(4,
                radioButtons(ns("d_sex"), "Sex / stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female"),
                tags$h5("Statistical significance"),
                numericInput(ns("d_fdr"), "Region-level FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                tags$h5("Effect size"),
                numericInput(ns("d_dbeta"), "Minimum absolute mean Δβ", value = 0, min = 0, max = 1, step = 0.01),
                radioButtons(ns("d_direction"), "Direction", inline = TRUE,
                             choices = c("Any" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),
                tags$h5("Region definition (display filters)"),
                numericInput(ns("d_mincpgs"), "Minimum CpGs per region", value = 3, min = 2, step = 1),
                numericInput(ns("d_minwidth"), "Minimum region size (bp, 0 = no limit)", value = 0, min = 0, step = 50),
                numericInput(ns("d_maxwidth"), "Maximum region size (bp, 0 = no limit)", value = 0, min = 0, step = 500),
                actionButton(ns("d_run_btn"), "Run DMR Analysis", icon = icon("play"), class = "btn-primary btn-sm")
              ),
              column(8, withSpinner(uiOutput(ns("d_volcano_ui")), color = "#2563EB", type = 6))
            ),
            uiOutput(ns("d_valueboxes_ui")),
            uiOutput(ns("d_results_ui"))
        )
      )
    })

    ## Gated behind "Run DMR Analysis" (same idea as mod_methyl_dmp.R's sva_run).
    d_run <- eventReactive(input$d_run_btn, {
      d <- default_data(); req(d)
      df <- if (identical(input$d_sex, "male")) d$dmr_m else d$dmr_f
      validate(need(!is.null(df) && nrow(df) > 0, "No precomputed DMR results are available for this stratum in this deployment."))
      df$direction <- ifelse(df$meandiff > 0, "hyper", "hypo")
      max_w <- if (isTRUE((input$d_maxwidth %||% 0) > 0)) input$d_maxwidth else Inf
      list(
        df = df, sex = input$d_sex,
        fdr = input$d_fdr %||% 0.05, dbeta = input$d_dbeta %||% 0,
        min_cpgs = input$d_mincpgs %||% 2, min_width = input$d_minwidth %||% 0, max_width = max_w,
        direction = if (identical(input$d_direction, "any")) NULL else input$d_direction
      )
    }, ignoreInit = TRUE)

    d_has_run <- reactiveVal(FALSE)
    observeEvent(input$d_run_btn, d_has_run(TRUE), ignoreInit = TRUE)

    d_filtered <- reactive({
      r <- d_run()
      mod_methyl_dmr_filter(r$df, "dmr_fdr", "meandiff", r$fdr, r$dbeta, r$direction,
                             min_cpgs = r$min_cpgs, min_width = r$min_width, max_width = r$max_width)
    })

    output$d_volcano_ui <- renderUI({
      if (!d_has_run()) {
        return(p(class = "empty-note", icon("circle-info"), "Pick a stratum and filters on the left, then click \"Run DMR Analysis\"."))
      }
      withSpinner(plotOutput(ns("d_volcano"), height = 320), color = "#2563EB", type = 6)
    })

    output$d_volcano <- renderPlot({
      r <- d_run()
      mod_methyl_dmp_volcano(r$df, "meandiff", "dmr_fdr", "Mean Δβ (RA - Control)", r$fdr, r$dbeta)
    })

    output$d_valueboxes_ui <- renderUI({
      req(d_has_run())
      r <- d_run(); f <- d_filtered()
      n_sig <- sum(!is.na(r$df$dmr_fdr) & r$df$dmr_fdr < r$fdr & abs(r$df$meandiff) >= r$dbeta, na.rm = TRUE)
      sex_label <- if (identical(r$sex, "male")) "Male" else "Female"
      div(class = "methyl-stats-row",
        fluidRow(
          valueBox(format(nrow(r$df), big.mark = ","), sprintf("Candidate DMRs (%s)", tolower(sex_label)), icon = icon("map-location-dot"), color = "purple", width = 3),
          valueBox(format(n_sig, big.mark = ","), sprintf("Region FDR < %s", r$fdr), icon = icon("star"), color = if (n_sig > 0) "green" else "light-blue", width = 3),
          valueBox(format(nrow(f), big.mark = ","), "Passing all current filters", icon = icon("filter"), color = "light-blue", width = 3),
          valueBox(sex_label, "Sex / stratum", icon = icon(if (identical(r$sex, "male")) "mars" else "venus"), color = "black", width = 3)
        ),
        fluidRow(
          valueBox(format(sum(f$meandiff > 0), big.mark = ","), "Hypermethylated (filtered)", icon = icon("arrow-up"), color = "red", width = 6),
          valueBox(format(sum(f$meandiff < 0), big.mark = ","), "Hypomethylated (filtered)", icon = icon("arrow-down"), color = "blue", width = 6)
        )
      )
    })

    output$d_results_ui <- renderUI({
      req(d_has_run())
      tagList(
        div(class = "card",
            div(class = "card-title", icon("chart-column"), "Genomic (Manhattan) plot"),
            withSpinner(plotOutput(ns("d_manhattan"), height = 320), color = "#2563EB", type = 6)),
        div(class = "card",
            div(class = "card-title", icon("table"), sprintf("%s DMR Analysis: results table", if (identical(d_run()$sex, "male")) "Male" else "Female")),
            p(class = "submodule-desc", "Additional per-column filtering (chromosome, gene, direction, CpG count) is available in the table's own search boxes. Click a row to inspect that region's constituent CpGs below."),
            div(class = "table-toolbar",
                downloadButton(ns("d_download_full"), "Complete results", class = "btn-default btn-sm"),
                downloadButton(ns("d_download_sig"), "Significant DMRs", class = "btn-default btn-sm"),
                downloadButton(ns("d_download_filtered"), "Filtered table", class = "btn-default btn-sm"),
                downloadButton(ns("d_download_config"), "Analysis configuration", class = "btn-default btn-sm")),
            DT::dataTableOutput(ns("d_table"))
        ),
        uiOutput(ns("d_region_ui"))
      )
    })

    output$d_manhattan <- renderPlot({
      r <- d_run()
      mod_methyl_dmr_manhattan(r$df, r$fdr)
    })

    output$d_table <- DT::renderDataTable({
      df <- d_filtered()
      show_cols <- intersect(c("seqnames", "start", "end", "width", "no.cpgs", "meandiff", "maxdiff", "Stouffer", "dmr_fdr", "overlapping.genes", "direction"), colnames(df))
      DT::datatable(df[, show_cols, drop = FALSE], rownames = FALSE, filter = "top", selection = "single",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = intersect(c("meandiff", "maxdiff", "Stouffer", "dmr_fdr"), show_cols), digits = 4)
    })
    ## Workaround: DT tables nested inside two levels of dynamic renderUI
    ## never get sent to the browser under Shiny's default suspendWhenHidden
    ## visibility tracking in this layout (confirmed server-side via
    ## testServer - the data/columns are correct). Forcing this and every
    ## other DT output here to always compute/send fixes it.
    outputOptions(output, "d_table", suspendWhenHidden = FALSE)

    output$d_download_full <- downloadHandler(
      filename = function() sprintf("dmr_default_%s_complete.csv", d_run()$sex),
      content = function(file) utils::write.csv(d_run()$df, file, row.names = FALSE)
    )
    output$d_download_sig <- downloadHandler(
      filename = function() sprintf("dmr_default_%s_significant.csv", d_run()$sex),
      content = function(file) {
        r <- d_run()
        sig <- r$df[!is.na(r$df$dmr_fdr) & r$df$dmr_fdr < r$fdr & abs(r$df$meandiff) >= r$dbeta, , drop = FALSE]
        utils::write.csv(sig, file, row.names = FALSE)
      }
    )
    output$d_download_filtered <- downloadHandler(
      filename = function() sprintf("dmr_default_%s_filtered.csv", d_run()$sex),
      content = function(file) utils::write.csv(d_filtered(), file, row.names = FALSE)
    )
    output$d_download_config <- downloadHandler(
      filename = function() "dmr_default_analysis_configuration.csv",
      content = function(file) {
        r <- d_run()
        cfg <- data.frame(
          parameter = c("dataset", "sex", "comparison", "statistical_method", "region_fdr_threshold", "min_abs_dbeta",
                        "direction_filter", "min_cpgs_per_region", "min_region_width_bp", "max_region_width_bp",
                        "n_candidate_regions", "n_significant_regions"),
          value = c("Preloaded whole-blood dataset (GSE42861)", r$sex, "RA vs Control",
                    "DMRcate (lambda=1000, C=2) on SVA-adjusted, bacon-corrected per-CpG statistics",
                    r$fdr, r$dbeta, input$d_direction %||% "any", r$min_cpgs, r$min_width,
                    if (is.infinite(r$max_width)) 0 else r$max_width,
                    nrow(r$df), sum(!is.na(r$df$dmr_fdr) & r$df$dmr_fdr < r$fdr & abs(r$df$meandiff) >= r$dbeta)),
          stringsAsFactors = FALSE
        )
        utils::write.csv(cfg, file, row.names = FALSE)
      }
    )

    ## ---- Region-level CpG inspection (default section) --------------------
    ## Per-sample β values aren't available for the preloaded pipeline's
    ## reproduced results (only completed result tables are bundled), so
    ## this shows constituent CpGs' per-CpG Δβ/FDR from the DMP tab's
    ## default analysis instead.
    d_region_selected <- reactive({
      sel <- input$d_table_rows_selected
      req(length(sel) == 1)
      row <- d_filtered()[sel, , drop = FALSE]
      r <- d_run()
      dmp_df <- load_default_dmp("sva", r$sex)
      pos <- methyl_champ_probe_positions()
      validate(need(!is.null(dmp_df) && !is.null(pos), "Per-CpG constituent data is not available for region-level inspection in this deployment."))
      m <- merge(dmp_df, pos, by = "cpg")
      m <- m[m$chr == as.character(row$seqnames) & m$pos >= row$start & m$pos <= row$end, , drop = FALSE]
      m <- m[order(m$pos), , drop = FALSE]
      list(row = row, cpgs = m)
    })

    output$d_region_ui <- renderUI({
      div(class = "card",
          div(class = "card-title", icon("magnifying-glass-chart"), "Region-level inspection"),
          if (is.null(input$d_table_rows_selected) || length(input$d_table_rows_selected) == 0)
            p(class = "empty-note", icon("circle-info"),
              "Select a row in the table above to inspect that region's constituent CpGs (per-CpG Δβ and FDR from the DMP tab's default analysis; per-sample values aren't available for the preloaded pipeline's reproduced results).")
          else tagList(
            uiOutput(ns("d_region_summary")),
            DT::dataTableOutput(ns("d_region_table"))
          )
      )
    })

    output$d_region_summary <- renderUI({
      d <- d_region_selected(); row <- d$row
      tagList(
        p(strong("Region: "), sprintf("%s:%s-%s", row$seqnames, format(row$start, big.mark = ","), format(row$end, big.mark = ","))),
        p(strong("Width: "), sprintf("%s bp", format(row$width, big.mark = ",")), " | ", strong("CpGs in region: "), row$no.cpgs,
          " | ", strong("Gene(s): "), if (!is.na(row$overlapping.genes) && nzchar(row$overlapping.genes)) row$overlapping.genes else "(none annotated)"),
        p(strong("Mean Δβ: "), sprintf("%.4f", row$meandiff), sprintf(" (%s)", row$direction),
          " | ", strong("Region FDR: "), signif(row$dmr_fdr, 4), " | ", strong("Constituent CpGs found: "), nrow(d$cpgs))
      )
    })

    output$d_region_table <- DT::renderDataTable({
      d <- d_region_selected()
      validate(need(nrow(d$cpgs) > 0, "No constituent CpG could be matched to this region's coordinates."))
      DT::datatable(d$cpgs[, c("cpg", "pos", "dbeta", "p_bacon", "fdr_bacon")], rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("dbeta", "p_bacon", "fdr_bacon"), digits = 4)
    })
    outputOptions(output, "d_region_table", suspendWhenHidden = FALSE)  ## see d_table's own comment above

    ## ================= 2. DMR Analysis (configurable live engine) =========
    ## Sex/group/covariate/QC-filter machinery reused directly from
    ## mod_methyl_dmp.R's live engine; the per-CpG limma fit below is the
    ## same step, extended with DMRcate region calling on top.

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
      ## Available whenever a real matrix is loaded (upload, or the
      ## preloaded dataset's own live matrix if bundled). Offers an
      ## All-Samples live region call plus single-sex live reruns, in
      ## addition to the Female/Male panels above which only reproduce the
      ## published pipeline's precomputed DMRcate tables.
      if (!methyl_dmr_engine_pkgs_ok()) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "DMR Analysis"),
          p(class = "submodule-desc", "The DMRcate/GenomicRanges packages needed for region-level calling are not installed in this deployment.")
        ))
      }
      if (is.null(methyl_dataset$beta)) {
        msg <- if (isTRUE(methyl_dataset$preloaded))
          "The preloaded dataset's live beta matrix isn't available in this deployment, so only the sex-stratified reproduced analysis above is available (no live All-Samples/Female-only/Male-only option)."
        else
          "Upload a beta/M-value matrix or IDAT files on the Methylomics Dataset tab to configure and run a live region-level differential methylation (DMR) analysis - including an All-Samples (combined) option."
        return(div(class = "card",
          div(class = "card-title", icon("upload"), "DMR Analysis"),
          p(class = "submodule-desc", msg)
        ))
      }
      if (is.null(methyl_dataset$sample_sheet)) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "No sample sheet"),
          p(class = "submodule-desc", "A DMR model needs a group variable - re-load on the Dataset tab with a sample sheet/phenotype file included.")
        ))
      }
      sheet <- methyl_dataset$sample_sheet
      cols <- colnames(sheet)
      sc <- sex_col()
      ns_ <- dataset_norm_status()
      ar <- anno_result()
      arr_ok <- isTRUE(ar$ok)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("map-location-dot"), "DMR Analysis"),
            p(class = "empty-note", icon("circle-info"),
              sprintf("Dataset: %s. %s probes x %s samples. %s",
                      methyl_dataset$source %||% "(unnamed)",
                      format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta),
                      ns_$message)),
            if (!arr_ok) p(class = "empty-note", icon("triangle-exclamation"),
                           sprintf("Region calling needs chromosome/position annotation for every tested CpG, which is unavailable for this array type in this deployment: %s", ar$reason))
            else fluidRow(
              column(4,
                tags$h5("Sex"),
                radioButtons(ns("live_sex"), NULL, inline = TRUE, choices = mod_methyl_dmp_sex_choices(sheet, sc), selected = "__all__"),
                if (length(mod_methyl_dmp_sex_choices(sheet, sc)) <= 1) p(class = "empty-note", icon("circle-info"), "No usable sex information was found for this dataset - showing pooled analysis only."),

                tags$h5("Comparison"),
                selectInput(ns("live_group_col"), "Group column", choices = cols,
                            selected = intersect(c("group", "Group", "disease", "Disease"), cols)[1] %||% cols[1]),
                uiOutput(ns("live_level_ui")),
                uiOutput(ns("live_comparison_label_ui")),

                tags$h5("Statistical significance"),
                numericInput(ns("live_fdr"), "Region-level FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                numericInput(ns("live_seed_p"), "CpG seeding p-value (raw, seeds candidate regions)", value = 0.05, min = 0.0001, max = 1, step = 0.01),

                tags$h5("Effect size"),
                numericInput(ns("live_dbeta"), "Minimum absolute mean Δβ", value = 0.05, min = 0, max = 1, step = 0.01),
                radioButtons(ns("live_direction"), "Direction", inline = TRUE,
                             choices = c("All DMRs" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),

                tags$h5("Region definition"),
                numericInput(ns("live_mincpgs"), "Minimum CpGs per region", value = 3, min = 2, step = 1),
                numericInput(ns("live_minwidth"), "Minimum region size (bp, 0 = no limit)", value = 0, min = 0, step = 50),
                numericInput(ns("live_maxwidth"), "Maximum region size (bp, 0 = no limit)", value = 0, min = 0, step = 500),

                tags$h5("CpG coverage / data quality"),
                numericInput(ns("live_min_valid_pct"), "Minimum valid (non-missing) sample %", value = 80, min = 0, max = 100, step = 5),
                numericInput(ns("live_min_variance"), "Minimum methylation variance (optional)", value = 0, min = 0, step = 0.001),
                checkboxInput(ns("live_snp_filter"), "Remove SNP-associated probes (manifest Probe_rs/CpG_rs/SBE_rs)", value = FALSE),

                tags$h5("Covariates (optional)"),
                uiOutput(ns("live_covariate_ui")),

                checkboxInput(ns("live_show_advanced"), "Show advanced DMR settings", value = FALSE),
                conditionalPanel(condition = sprintf("input['%s']", ns("live_show_advanced")),
                  numericInput(ns("live_lambda"), "Kernel bandwidth / max CpG gap (lambda, nt)", value = 1000, min = 100, step = 100),
                  numericInput(ns("live_c"), "Bandwidth scaling factor (C)", value = 2, min = 0.2, step = 0.1),
                  checkboxInput(ns("live_pcutoff_manual"), "Override candidate-region kernel-smoothed FDR cutoff (not recommended)", value = FALSE),
                  conditionalPanel(condition = sprintf("input['%s']", ns("live_pcutoff_manual")),
                    numericInput(ns("live_pcutoff"), "Candidate-region cutoff", value = 0.05, min = 0, max = 1, step = 0.01))
                ),

                actionButton(ns("live_run_btn"), "Run DMR Analysis", icon = icon("play"), class = "btn-primary")
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
        selectInput(ns("live_ref"), "Group 1 (reference)", choices = levels_available, selected = levels_available[1]),
        selectInput(ns("live_comp"), "Group 2 (comparison)", choices = levels_available, selected = levels_available[min(2, length(levels_available))])
      )
    })

    output$live_comparison_label_ui <- renderUI({
      req(input$live_ref, input$live_comp)
      validate(need(input$live_ref != input$live_comp, "Group 1 and Group 2 must be different."))
      sex_lbl <- mod_methyl_dmr_sex_label(methyl_dataset$sample_sheet, sex_col(), input$live_sex %||% "__all__")
      prefix <- if (identical(sex_lbl, "All samples")) "" else paste0(sex_lbl, " ")
      p(class = "empty-note", icon("code-compare"), strong("Comparison: "),
        sprintf("%s%s vs %s%s", prefix, input$live_comp, prefix, input$live_ref))
    })

    output$live_covariate_ui <- renderUI({
      sheet <- methyl_dataset$sample_sheet
      req(sheet, input$live_group_col)
      exclude <- c(id_cols(), input$live_group_col)
      if (!identical(input$live_sex %||% "__all__", "__all__") && !is.null(sex_col())) exclude <- c(exclude, sex_col())
      cand <- mod_methyl_dmp_covariate_cols(sheet, exclude)
      if (length(cand) == 0) {
        return(p(class = "empty-note", icon("circle-info"), "No additional phenotype columns are available to use as covariates."))
      }
      checkboxGroupInput(ns("live_covariates"), NULL, choices = cand, selected = character(0))
    })

    live_has_run <- reactiveVal(FALSE)
    observeEvent(methyl_dataset$beta, live_has_run(FALSE), ignoreNULL = TRUE)

    live_result <- eventReactive(input$live_run_btn, withProgress(message = "Running DMR Analysis - fitting a full genome-wide array (400k+ CpGs) and calling regions can take a few minutes...", value = 0.2, {
      validate(need(!is.null(methyl_dataset$beta), "Load a dataset first."))
      sheet <- methyl_dataset$sample_sheet
      validate(need(!is.null(sheet), "No sample sheet loaded."))
      validate(need(!is.null(input$live_group_col) && input$live_group_col %in% colnames(sheet), "Pick a group column."))
      validate(need(!is.null(input$live_ref) && !is.null(input$live_comp), "Pick Group 1 and Group 2."))
      validate(need(input$live_ref != input$live_comp, "Group 1 and Group 2 must be different."))
      ar <- anno_result()
      validate(need(isTRUE(ar$ok), sprintf("Region calling needs chromosome/position annotation, unavailable for this array type in this deployment%s.",
                                            if (!is.null(ar$reason)) paste0(": ", ar$reason) else "")))

      sample_ids <- methyl_sheet_sample_ids(sheet, colnames(methyl_dataset$beta))
      common <- intersect(colnames(methyl_dataset$beta), sample_ids)
      validate(need(length(common) >= 6, "Fewer than 6 samples match between the matrix and the sample sheet."))
      beta0 <- methyl_dataset$beta[, common, drop = FALSE]
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
        sex_label <- mod_methyl_dmr_sex_label(sheet, sc, sex_choice)
      }

      ## ---- group subset --------------------------------------------------
      grp_raw <- as.character(ph0[[input$live_group_col]])
      keep_grp <- !is.na(grp_raw) & grp_raw %in% c(input$live_ref, input$live_comp)
      beta1 <- beta0[, keep_grp, drop = FALSE]
      ph1 <- ph0[keep_grp, , drop = FALSE]
      grp <- factor(grp_raw[keep_grp], levels = c(input$live_ref, input$live_comp))
      rm(beta0)  ## no further use once beta1 exists; frees a full-size copy

      ## ---- optional covariates (complete cases only) ----------------------
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
        sprintf("Each group needs at least 3 samples to fit a model (Group 1 \"%s\": %d, Group 2 \"%s\": %d).",
                input$live_ref, n_ref, input$live_comp, n_comp)))

      ## ---- probe filters (missingness / variance / SNP / position), BEFORE
      ## the M-value transform -------------------------------------------------
      ## Filtering runs on beta1 itself (no extra copy needed when
      ## input_scale=="beta", the common case; only the M-value-upload path
      ## needs one temporary beta-scale conversion first). `m`/`beta_scale`
      ## are then derived only from the already-filtered subset - keeps at
      ## most one full-size matrix alive at a time on the ~400k-probe array,
      ## instead of 3+ copies from converting before filtering.
      is_m_scale <- identical(methyl_dataset$input_scale, "m")
      beta_scale_full <- if (is_m_scale) 2^beta1 / (1 + 2^beta1) else beta1

      max_na_frac <- 1 - (input$live_min_valid_pct %||% 80) / 100
      f_miss <- methyl_filter_missing(beta_scale_full, max_na_frac)
      keep_probe <- f_miss$keep
      f_var <- methyl_filter_variance(beta_scale_full, input$live_min_variance %||% 0)
      keep_probe <- keep_probe & f_var$keep
      snp_note <- NULL
      if (isTRUE(input$live_snp_filter)) {
        f_snp <- methyl_filter_snp(beta_scale_full, ar)
        keep_probe <- keep_probe & f_snp$keep
        snp_note <- f_snp$note
      }
      ## Region calling needs a chromosome + position for every tested probe;
      ## probes absent from the manifest can't be placed.
      a <- ar$anno
      hit <- match(rownames(beta_scale_full), rownames(a))
      has_pos <- !is.na(hit) & !is.na(a$chr[hit]) & !is.na(a$pos[hit])
      keep_probe <- keep_probe & has_pos
      validate(need(sum(keep_probe) >= 20,
        "Fewer than 20 CpGs remain after the missingness/variance/SNP/position filters. Relax the filters and try again."))

      beta1 <- beta1[keep_probe, , drop = FALSE]
      if (is_m_scale) rm(beta_scale_full)  ## extra copy only in this branch
      gc(FALSE)  ## reclaim before the memory-heavy limma fit below

      ## ---- beta/M-value matrices (filtered subset only) ---------------------
      m <- if (is_m_scale) beta1 else log2(pmin(pmax(beta1, 1e-6), 1 - 1e-6) / (1 - pmin(pmax(beta1, 1e-6), 1 - 1e-6)))
      beta_scale <- if (is_m_scale) 2^m / (1 + 2^m) else beta1

      ## ---- design + limma fit (per-CpG statistics feeding DMRcate) ---------
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
                      error = function(e) validate(need(FALSE, "Could not build the Group 1/Group 2 contrast - group names may contain characters limma can't use directly (try renaming the group levels in your sample sheet).")))
      fit2 <- tryCatch(limma::eBayes(limma::contrasts.fit(fit, cm)),
                        error = function(e) validate(need(FALSE, paste("limma could not fit this model:", conditionMessage(e)))))
      tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

      beta_ref <- rowMeans(beta_scale[, grp == input$live_ref, drop = FALSE], na.rm = TRUE)
      beta_comp <- rowMeans(beta_scale[, grp == input$live_comp, drop = FALSE], na.rm = TRUE)
      dbeta <- (beta_comp - beta_ref)[rownames(tt)]

      hit2 <- match(rownames(tt), rownames(a))
      chr <- a$chr[hit2]; pos <- a$pos[hit2]

      ## ---- CpGannotated + DMRcate region calling ----------------------------
      ## Same construction as script04_dmr_sexstratified.R's build_annot():
      ## a GRanges of per-CpG stat/rawpval/diff/ind.fdr/is.sig, built
      ## directly (not via cpg.annotate()'s internal refit) so DMRcate
      ## operates on exactly the limma fit above. `is.sig` seeds candidate
      ## regions from the raw CpG-level p-value (a DMRcate tuning parameter,
      ## see ?dmrcate "pcutoff"), separate from `ind.fdr`.
      seed_p <- input$live_seed_p %||% 0.05
      ind_fdr <- stats::p.adjust(tt$P.Value, "BH")
      ## CpGannotated S4 class needs DMRcate's namespace loaded before methods::new() below.
      requireNamespace("DMRcate", quietly = TRUE)
      gr <- GenomicRanges::GRanges(
        seqnames = chr, ranges = IRanges::IRanges(pos, pos),
        stat = tt$t, rawpval = tt$P.Value, diff = dbeta, ind.fdr = ind_fdr,
        is.sig = tt$P.Value < seed_p
      )
      names(gr) <- rownames(tt)
      n_seed <- sum(gr$is.sig)
      validate(need(n_seed >= 2, sprintf("Only %d CpG(s) pass the seeding p-value threshold (raw p < %s) - relax it, or check the comparison/covariates.", n_seed, seed_p)))
      annot <- methods::new("CpGannotated", ranges = gr)

      lambda <- input$live_lambda %||% 1000
      C <- input$live_c %||% 2
      min_cpgs <- max(2, as.integer(input$live_mincpgs %||% 3))
      pcutoff <- if (isTRUE(input$live_pcutoff_manual)) (input$live_pcutoff %||% 0.05) else "fdr"

      dmr_raw <- tryCatch(
        DMRcate::dmrcate(annot, lambda = lambda, C = C, min.cpgs = min_cpgs, pcutoff = pcutoff),
        error = function(e) validate(need(FALSE, paste("DMRcate could not call regions with these settings:", conditionMessage(e))))
      )
      ranges <- tryCatch(
        DMRcate::extractRanges(dmr_raw, genome = "hg19"),
        error = function(e) validate(need(FALSE, paste("Could not extract genomic ranges for the candidate regions:", conditionMessage(e))))
      )
      validate(need(length(ranges) >= 1, "No candidate DMRs were returned with these settings - relax the seeding p-value, lambda, or minimum CpGs."))

      dt <- as.data.frame(ranges)
      dt$seqnames <- as.character(dt$seqnames)
      dt$overlapping.genes <- vapply(dt$overlapping.genes, function(g) {
        if (length(g) == 0 || all(is.na(g))) NA_character_ else paste(unique(as.character(g)), collapse = ";")
      }, character(1))
      dt$dmr_fdr <- stats::p.adjust(dt$Stouffer, method = "BH")
      dt$direction <- ifelse(dt$meandiff > 0, "hyper", "hypo")
      dt <- dt[order(dt$dmr_fdr), , drop = FALSE]
      dt$dmr_id <- sprintf("DMR%03d", seq_len(nrow(dt)))

      ## Region-level per-group mean methylation, computed from the filtered
      ## beta matrix and reported alongside DMRcate's meandiff/maxdiff, via
      ## a vectorized region-CpG join (GenomicRanges::findOverlaps()).
      region_gr <- GenomicRanges::GRanges(dt$seqnames, IRanges::IRanges(dt$start, dt$end))
      probe_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
      hits <- GenomicRanges::findOverlaps(probe_gr, region_gr)
      qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
      dt$ref_mean_beta <- NA_real_; dt$comp_mean_beta <- NA_real_
      if (length(qh) > 0) {
        ref_by_region <- tapply(beta_ref[rownames(tt)[qh]], sh, mean, na.rm = TRUE)
        comp_by_region <- tapply(beta_comp[rownames(tt)[qh]], sh, mean, na.rm = TRUE)
        idx <- as.integer(names(ref_by_region))
        dt$ref_mean_beta[idx] <- as.numeric(ref_by_region)
        dt$comp_mean_beta[idx] <- as.numeric(comp_by_region)
      }

      design_formula <- sprintf("Methylation ~ %s%s", input$live_group_col,
                                 if (length(cov_names) > 0) paste0(" + ", paste(cov_names, collapse = " + ")) else "")

      list(
        dt = dt, ref = input$live_ref, comp = input$live_comp, group_col = input$live_group_col,
        sex_label = sex_label, sex_col = sc, covariates = cov_names, design_formula = design_formula,
        n_ref = n_ref, n_comp = n_comp, n_cpgs_tested = nrow(tt), n_cpgs_before_filter = nrow(beta1),
        n_seed_cpgs = n_seed, seed_p = seed_p, lambda = lambda, C = C, min_cpgs = min_cpgs, pcutoff = pcutoff,
        min_valid_pct = input$live_min_valid_pct %||% 80, min_variance = input$live_min_variance %||% 0,
        snp_filter = isTRUE(input$live_snp_filter), snp_note = snp_note,
        missing_note = f_miss$note, variance_note = f_var$note,
        norm_status = dataset_norm_status(), dataset_source = methyl_dataset$source %||% "(unnamed)",
        preloaded = isTRUE(methyl_dataset$preloaded),
        beta_scale = beta_scale, grp = grp, probe_chr = chr, probe_pos = pos,
        lambda_gc = mod_methyl_lambda_gc(tt$P.Value), cpg_p_raw = tt$P.Value,
        run_at = Sys.time()
      )
    }), ignoreInit = TRUE)

    observeEvent(live_result(), {
      r <- live_result()
      methyl_results$dmr <- list(
        comparison = sprintf("%s vs %s (%s)", r$comp, r$ref, r$sex_label),
        n_regions = nrow(r$dt),
        n_sig = sum(!is.na(r$dt$dmr_fdr) & r$dt$dmr_fdr < 0.05, na.rm = TRUE)
      )
      ## Publishes the full live region table for Candidate CpGs (Module-DMR
      ## Overlap) to pick up automatically, without a re-upload, once this ran
      ## against an uploaded/GEO dataset - skipped for preloaded since
      ## Candidate CpGs' own "Use Preloaded Data" path reproduces the
      ## published DMR results from static files instead.
      if (!isTRUE(methyl_dataset$preloaded)) {
        methyl_results$dmr_table <- r$dt
        methyl_results$dmr_meta <- list(comp = r$comp, ref = r$ref, sex_label = r$sex_label)
      }
      live_has_run(TRUE)
    })

    ## FDR/Δβ/direction/CpG-count/width stay live-adjustable after a run;
    ## only sex/group/covariate/seeding/lambda/C/probe-filter changes (which
    ## change what was actually called) require clicking Run again.
    live_filtered <- reactive({
      r <- live_result()
      max_w <- if (isTRUE((input$live_maxwidth %||% 0) > 0)) input$live_maxwidth else Inf
      direction <- if (identical(input$live_direction, "any")) NULL else input$live_direction
      mod_methyl_dmr_filter(r$dt, "dmr_fdr", "meandiff", input$live_fdr %||% 0.05, input$live_dbeta %||% 0, direction,
                             min_cpgs = input$live_mincpgs %||% 2, min_width = input$live_minwidth %||% 0, max_width = max_w)
    })

    live_sig <- reactive({
      r <- live_result()
      r$dt[!is.na(r$dt$dmr_fdr) & r$dt$dmr_fdr <= (input$live_fdr %||% 0.05) & abs(r$dt$meandiff) >= (input$live_dbeta %||% 0), , drop = FALSE]
    })

    output$live_results_ui <- renderUI({
      req(live_has_run(), live_result())
      r <- live_result()
      sig <- live_sig()
      n_sig <- nrow(sig)
      n_hyper <- sum(sig$meandiff > 0)
      n_hypo <- sum(sig$meandiff < 0)
      sex_heading <- if (identical(r$sex_label, "All samples")) "All-Sample" else r$sex_label

      tagList(
        div(class = "card",
            div(class = "card-title", icon("circle-check"), "Configuration & sample sizes"),
            p(strong("Model: "), tags$code(r$design_formula),
              sprintf(" (limma moderated t-test per CpG; DMRcate kernel region calling, lambda=%s nt, C=%s, min.cpgs=%s).", r$lambda, r$C, r$min_cpgs)),
            p(strong("Sex: "), r$sex_label, " | ", strong("Group 1: "), r$ref, sprintf(" (n=%d)", r$n_ref),
              " | ", strong("Group 2: "), r$comp, sprintf(" (n=%d)", r$n_comp),
              " | ", strong("Total analyzed: "), r$n_ref + r$n_comp),
            if (length(r$covariates) > 0) p(strong("Covariates: "), paste(r$covariates, collapse = ", ")) else p(class = "submodule-desc", "No covariates selected."),
            p(class = "submodule-desc", r$norm_status$message),
            if ((r$n_ref < 10 || r$n_comp < 10))
              p(class = "empty-note", icon("triangle-exclamation"), "One or both groups have fewer than 10 samples - results may be underpowered."),
            p(class = "submodule-desc", sprintf("%s CpGs tested after QC filters: %s of %s. Seeded %s CpG(s) at raw p < %s for region calling.",
                                                  paste(c(r$missing_note, r$variance_note, r$snp_note), collapse = " "),
                                                  format(r$n_cpgs_tested, big.mark = ","), format(r$n_cpgs_before_filter, big.mark = ","),
                                                  format(r$n_seed_cpgs, big.mark = ","), r$seed_p))
        ),
        div(class = "card",
            div(class = "card-title", icon("gauge-high"), "Genomic inflation diagnostic"),
            p(strong("Genomic inflation factor (λ): "), if (is.na(r$lambda_gc)) "not available" else sprintf("%.2f", r$lambda_gc),
              " - the ratio of observed to expected median test statistic across all tested CpGs (before region seeding); 1.0 indicates no inflation."),
            if (!is.na(r$lambda_gc) && r$lambda_gc > 1.1)
              p(class = "empty-note", icon("triangle-exclamation"),
                "Genomic inflation is elevated (λ > 1.1) - this live engine doesn't apply SVA/bacon correction, so treat these region-level p-values/FDR as potentially optimistic."),
            withSpinner(plotOutput(ns("live_qq"), height = 320), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), sprintf("%s DMR Analysis: Summary", sex_heading)),
            div(class = "methyl-stats-row",
              fluidRow(
                valueBox(format(nrow(r$dt), big.mark = ","), "Candidate DMRs detected", icon = icon("map-location-dot"), color = "purple", width = 3),
                valueBox(format(n_sig, big.mark = ","), "Significant (FDR + Δβ)", icon = icon("star"), color = if (n_sig > 0) "green" else "light-blue", width = 3),
                valueBox(format(n_hyper, big.mark = ","), "Hypermethylated", icon = icon("arrow-up"), color = "red", width = 3),
                valueBox(format(n_hypo, big.mark = ","), "Hypomethylated", icon = icon("arrow-down"), color = "blue", width = 3)
              ),
              fluidRow(
                valueBox(format(sum(!is.na(r$dt$dmr_fdr) & r$dt$dmr_fdr <= (input$live_fdr %||% 0.05)), big.mark = ","),
                          sprintf("Pass region FDR < %s", input$live_fdr %||% 0.05), icon = icon("filter"), color = "light-blue", width = 4),
                valueBox(format(sum(!is.na(r$dt$meandiff) & abs(r$dt$meandiff) >= (input$live_dbeta %||% 0)), big.mark = ","),
                          sprintf("Pass |Δβ| ≥ %s", input$live_dbeta %||% 0), icon = icon("filter"), color = "light-blue", width = 4),
                valueBox(r$sex_label, "Sex subset", icon = icon("venus-mars"), color = "black", width = 4)
              )
            ),
            if (n_sig == 0) p(class = "empty-note", icon("circle-info"),
              "No DMRs passed the selected FDR and Δβ thresholds. Consider relaxing the thresholds, the seeding p-value, or reviewing the sample/group configuration.")
        ),
        div(class = "card", div(class = "card-title", icon("chart-scatter"), "Volcano plot"),
            withSpinner(plotOutput(ns("live_volcano"), height = 340), color = "#2563EB", type = 6)),
        div(class = "card", div(class = "card-title", icon("chart-column"), "Genomic (Manhattan) plot"),
            withSpinner(plotOutput(ns("live_manhattan"), height = 320), color = "#2563EB", type = 6)),
        div(class = "card",
            div(class = "card-title", icon("ranking-star"), "Top DMRs by significance"),
            selectInput(ns("live_top_n"), "Show top", choices = c(10, 20, 50, 100), selected = 20),
            withSpinner(plotOutput(ns("live_effectplot"), height = 380), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Methylation heatmap (significant DMRs x samples)"),
            if (n_sig > 50) p(class = "empty-note", icon("circle-info"), sprintf("Showing the top 50 of %d significant DMRs by region-level FDR - all %d appear in the results table and exports below.", n_sig, n_sig)),
            withSpinner(plotOutput(ns("live_heatmap"), height = 420), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("table"), sprintf("%s DMR Analysis: results table", sex_heading)),
            p(class = "submodule-desc", "Additional per-column filtering (chromosome, gene, direction, CpG count) is available in the table's own search boxes. Click a row to inspect that region's constituent CpGs below."),
            div(class = "table-toolbar",
                downloadButton(ns("download_live_full"), "Complete results", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_sig"), "Significant DMRs", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_filtered"), "Filtered table", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_anno"), "Annotation only", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_matrix"), "Region methylation matrix", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_config"), "Analysis configuration", class = "btn-default btn-sm")),
            DT::dataTableOutput(ns("live_table"))
        ),
        div(class = "card",
            div(class = "card-title", icon("magnifying-glass-chart"), "Region-level inspection"),
            uiOutput(ns("live_region_ui"))
        )
      )
    })

    output$live_qq <- renderPlot({
      r <- live_result()
      mod_methyl_qq_plot(r$cpg_p_raw)
    })

    output$live_volcano <- renderPlot({
      r <- live_result()
      mod_methyl_dmp_volcano(r$dt, "meandiff", "dmr_fdr", "Mean Δβ", input$live_fdr %||% 0.05, input$live_dbeta %||% 0)
    })

    output$live_manhattan <- renderPlot({
      r <- live_result()
      mod_methyl_dmr_manhattan(r$dt, input$live_fdr %||% 0.05)
    })

    output$live_effectplot <- renderPlot({
      r <- live_result()
      mod_methyl_dmr_topplot(r$dt, n = as.integer(input$live_top_n %||% 20))
    })

    output$live_heatmap <- renderPlot({
      r <- live_result()
      mod_methyl_dmr_heatmap(live_sig(), r$beta_scale, r$probe_chr, r$probe_pos, r$grp, max_regions = 50)
    })

    output$live_table <- DT::renderDataTable({
      df <- live_filtered()
      show_cols <- intersect(c("dmr_id", "seqnames", "start", "end", "width", "no.cpgs", "ref_mean_beta", "comp_mean_beta",
                                "meandiff", "maxdiff", "Stouffer", "dmr_fdr", "overlapping.genes", "direction"), colnames(df))
      df$significant <- ifelse(!is.na(df$dmr_fdr) & df$dmr_fdr <= (input$live_fdr %||% 0.05) & abs(df$meandiff) >= (input$live_dbeta %||% 0), "Yes", "No")
      DT::datatable(df[, c(show_cols, "significant")], rownames = FALSE, filter = "top", selection = "single",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = intersect(c("ref_mean_beta", "comp_mean_beta", "meandiff", "maxdiff", "Stouffer", "dmr_fdr"), show_cols), digits = 4)
    })
    outputOptions(output, "live_table", suspendWhenHidden = FALSE)  ## see d_table's own comment above

    output$download_live_full <- downloadHandler(
      filename = function() { r <- live_result(); sprintf("dmr_%s_vs_%s_%s_complete.csv", r$comp, r$ref, gsub(" ", "_", r$sex_label)) },
      content = function(file) utils::write.csv(live_result()$dt, file, row.names = FALSE)
    )
    output$download_live_sig <- downloadHandler(
      filename = function() { r <- live_result(); sprintf("dmr_%s_vs_%s_%s_significant.csv", r$comp, r$ref, gsub(" ", "_", r$sex_label)) },
      content = function(file) utils::write.csv(live_sig(), file, row.names = FALSE)
    )
    output$download_live_filtered <- downloadHandler(
      filename = function() { r <- live_result(); sprintf("dmr_%s_vs_%s_%s_filtered.csv", r$comp, r$ref, gsub(" ", "_", r$sex_label)) },
      content = function(file) utils::write.csv(live_filtered(), file, row.names = FALSE)
    )
    output$download_live_anno <- downloadHandler(
      filename = function() "dmr_annotation.csv",
      content = function(file) {
        df <- live_result()$dt
        cols <- intersect(c("dmr_id", "seqnames", "start", "end", "width", "no.cpgs", "overlapping.genes"), colnames(df))
        utils::write.csv(df[, cols, drop = FALSE], file, row.names = FALSE)
      }
    )
    output$download_live_matrix <- downloadHandler(
      filename = function() "dmr_region_methylation_matrix.csv",
      content = function(file) {
        r <- live_result()
        sig <- live_sig()
        validate(need(nrow(sig) > 0, "No significant DMRs to export a region methylation matrix for."))
        region_gr <- GenomicRanges::GRanges(sig$seqnames, IRanges::IRanges(sig$start, sig$end))
        probe_gr <- GenomicRanges::GRanges(r$probe_chr, IRanges::IRanges(r$probe_pos, r$probe_pos))
        hits <- GenomicRanges::findOverlaps(probe_gr, region_gr)
        qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
        mat <- matrix(NA_real_, nrow = nrow(sig), ncol = ncol(r$beta_scale), dimnames = list(sig$dmr_id, colnames(r$beta_scale)))
        for (i in seq_len(nrow(sig))) {
          probes_i <- qh[sh == i]
          if (length(probes_i) > 0) mat[i, ] <- colMeans(r$beta_scale[probes_i, , drop = FALSE], na.rm = TRUE)
        }
        out <- data.frame(dmr_id = rownames(mat), mat, check.names = FALSE)
        utils::write.csv(out, file, row.names = FALSE)
      }
    )
    output$download_live_config <- downloadHandler(
      filename = function() "dmr_analysis_configuration.csv",
      content = function(file) {
        r <- live_result()
        cfg <- data.frame(
          parameter = c("dataset", "sex", "group1", "group2", "n_group1", "n_group2", "n_total_analyzed", "covariates",
                        "statistical_method", "design_formula", "cpg_seeding_pvalue", "n_seed_cpgs", "lambda_bp", "C_scaling",
                        "min_cpgs_per_region", "candidate_region_pcutoff", "region_fdr_threshold", "min_abs_dbeta", "direction_filter",
                        "min_region_width_bp", "max_region_width_bp", "min_valid_sample_pct", "min_variance", "snp_filter_applied",
                        "n_cpgs_tested", "n_candidate_regions", "n_significant_regions", "run_at"),
          value = c(r$dataset_source, r$sex_label, r$ref, r$comp, r$n_ref, r$n_comp, r$n_ref + r$n_comp,
                    if (length(r$covariates) > 0) paste(r$covariates, collapse = ";") else "(none)",
                    "limma (per-CpG) + DMRcate (kernel-smoothed region calling)", r$design_formula,
                    r$seed_p, r$n_seed_cpgs, r$lambda, r$C, r$min_cpgs,
                    if (identical(r$pcutoff, "fdr")) "fdr (DMRcate default)" else r$pcutoff,
                    input$live_fdr %||% 0.05, input$live_dbeta %||% 0, input$live_direction %||% "any",
                    input$live_minwidth %||% 0, input$live_maxwidth %||% 0,
                    r$min_valid_pct, r$min_variance, r$snp_filter,
                    r$n_cpgs_tested, nrow(r$dt), nrow(live_sig()), format(r$run_at)),
          stringsAsFactors = FALSE
        )
        utils::write.csv(cfg, file, row.names = FALSE)
      }
    )

    ## ---- Region-level CpG inspection (live section) ------------------------
    ## Verifies a selected region is a coordinated, multi-CpG change rather
    ## than a single isolated CpG. Unlike the default section, real
    ## per-sample β values are available, so mod_methyl_dmp_betadist() can
    ## be reused directly for the per-CpG-by-group plot.
    live_region_selected <- reactive({
      sel <- input$live_table_rows_selected
      req(length(sel) == 1)
      r <- live_result()
      row <- live_filtered()[sel, , drop = FALSE]
      hit <- !is.na(r$probe_chr) & r$probe_chr == as.character(row$seqnames) & r$probe_pos >= row$start & r$probe_pos <= row$end
      validate(need(any(hit), "No analyzed CpG falls within this region's coordinates."))
      cpg_ids <- rownames(r$beta_scale)[hit]
      ord <- order(r$probe_pos[hit])
      cpg_ids <- cpg_ids[ord]
      list(row = row, cpg = cpg_ids, beta = r$beta_scale[cpg_ids, , drop = FALSE], grp = r$grp, ref = r$ref, comp = r$comp)
    })

    output$live_region_ui <- renderUI({
      if (is.null(input$live_table_rows_selected) || length(input$live_table_rows_selected) == 0) {
        return(p(class = "empty-note", icon("circle-info"),
                  "Select a row in the results table above to inspect that region's constituent CpGs and verify it reflects a coordinated, multi-CpG methylation change rather than a single isolated CpG."))
      }
      tagList(
        uiOutput(ns("live_region_summary")),
        withSpinner(plotOutput(ns("live_region_plot"), height = 320), color = "#2563EB", type = 6),
        DT::dataTableOutput(ns("live_region_table"))
      )
    })

    output$live_region_summary <- renderUI({
      d <- live_region_selected(); row <- d$row
      tagList(
        p(strong("Region: "), sprintf("%s:%s-%s", row$seqnames, format(row$start, big.mark = ","), format(row$end, big.mark = ","))),
        p(strong("Width: "), sprintf("%s bp", format(row$width, big.mark = ",")), " | ", strong("CpGs in region: "), row$no.cpgs,
          " | ", strong("Gene(s): "), if (!is.na(row$overlapping.genes) && nzchar(row$overlapping.genes)) row$overlapping.genes else "(none annotated)"),
        p(strong("Mean Δβ: "), sprintf("%.4f", row$meandiff), sprintf(" (%s)", row$direction),
          " | ", strong(sprintf("%s mean β: ", d$ref)), sprintf("%.4f", row$ref_mean_beta),
          " | ", strong(sprintf("%s mean β: ", d$comp)), sprintf("%.4f", row$comp_mean_beta),
          " | ", strong("Region FDR: "), signif(row$dmr_fdr, 4))
      )
    })

    output$live_region_plot <- renderPlot({
      d <- live_region_selected()
      mod_methyl_dmp_betadist(d$beta, d$cpg, d$grp)
    })

    output$live_region_table <- DT::renderDataTable({
      d <- live_region_selected()
      long <- as.data.frame(t(d$beta))
      long <- cbind(sample = rownames(long), group = as.character(d$grp), long)
      DT::datatable(long, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = d$cpg, digits = 4)
    })
    outputOptions(output, "live_region_table", suspendWhenHidden = FALSE)  ## see d_table's own comment above
  })
}
