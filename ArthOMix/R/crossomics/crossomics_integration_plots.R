## R/crossomics/crossomics_integration_plots.R
## Plotting/rendering helpers for mod_cross_integration.R - kept separate from
## the module file so the Shiny wiring (mod_cross_integration.R) and the
## visualization logic can be read/changed independently. Every plot reuses
## the app's shared theme_arthomix()/ARTHOMIX_COLORS (global.R) so it matches
## the rest of the app, and the quadrant plot is built once as a ggplot
## object (cx_quadrant_ggplot) reused for both the interactive plotly view
## and PNG/PDF/SVG export, so there is exactly one place that draws it.

cx_empty_state <- function() {
  div(class = "empty-note", icon("circle-info"), "Click \"Run Integration\" in the Integration tab to see results here.")
}

cx_fmt_num <- function(x, digits = 3, sci = FALSE) {
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (sci) format(v, digits = 2, scientific = TRUE) else format(round(v, digits), nsmall = digits)
  }, character(1))
}

## ---------------------------------------------------------------------------
## Gene Detail Drawer (modal)
## ---------------------------------------------------------------------------

cx_gene_detail_modal <- function(row, pairing, cpg_rows = NULL) {
  g <- row$gene[1]
  val <- function(field) if (is.null(row[[field]])) NA else row[[field]][1]
  fmt1 <- function(field, digits = 3, sci = FALSE) {
    v <- val(field)
    if (length(v) == 0 || is.na(v)) return("—")
    if (is.numeric(v)) cx_fmt_num(v, digits, sci) else as.character(v)
  }
  cat_label <- if (!is.na(val("category"))) as.character(val("category")) else NA
  cat_desc <- if (!is.na(cat_label)) CX_CATEGORY_LABELS[[cat_label]] else NA
  ev_label <- if (!is.na(val("evidence_level"))) as.character(val("evidence_level")) else NA
  ev_desc <- if (!is.na(ev_label)) CX_EVIDENCE_DESCRIPTIONS[[ev_label]] else NA
  has_cor <- !is.na(val("correlation_r"))
  modalDialog(
    title = paste("Gene:", g), size = "l", easyClose = TRUE, footer = modalButton("Close"),
    tags$h5("Expression"),
    tags$p(sprintf("log2FC: %s (%s)   FDR: %s", fmt1("log2fc"), fmt1("expression_direction"), fmt1("expr_fdr", sci = TRUE))),
    tags$h5("Methylation"),
    tags$p(sprintf("Δβ: %s (%s)   FDR: %s   Primary region: %s   CpGs mapped: %s   Significant CpGs: %s",
                    fmt1("dbeta"), fmt1("methylation_direction"), fmt1("meth_fdr", sci = TRUE),
                    fmt1("primary_region"), fmt1("n_cpg_total"), fmt1("n_cpg_significant"))),
    if (!is.null(cpg_rows) && nrow(cpg_rows) > 0) tagList(
      tags$h5(sprintf("CpG details (%d)", nrow(cpg_rows))),
      DT::datatable(
        cpg_rows[, intersect(c("cpg", "region_fine", "island_context", "dbeta", "fdr", "methylation_direction", "sig_cpg"), colnames(cpg_rows)), drop = FALSE],
        rownames = FALSE, options = list(pageLength = 5, scrollX = TRUE, dom = "tp"), class = "stripe hover compact"
      )
    ),
    if (!is.na(cat_label)) tagList(
      tags$h5("Quadrant classification"),
      tags$p(tags$strong(cat_label), " — ", cat_desc)
    ),
    if (!is.na(ev_label)) tagList(
      tags$h5("Evidence level"),
      tags$p(tags$strong(ev_label), " — ", ev_desc)
    ),
    if (has_cor) tagList(
      tags$h5("Expression–Methylation Association"),
      tags$p(sprintf("%s rho/r = %s, FDR = %s", fmt1("correlation_n"), fmt1("correlation_r"), fmt1("correlation_fdr", sci = TRUE))),
      tags$p(style = "font-style: italic;",
             if (!is.na(val("correlation_r")) && val("correlation_r") < 0) "Interpretation: negative association between methylation and expression." else "Interpretation: positive association between methylation and expression.")
    ) else if (!isTRUE(pairing$paired)) tagList(
      tags$h5("Expression–Methylation Association"),
      tags$p(class = "empty-note", icon("triangle-exclamation"), "Unpaired datasets - sample-level correlation not available for this gene.")
    ),
    tags$hr(),
    tags$p(style = "font-style: italic; color: var(--color-ink-muted, #898781);",
           "This represents an observed statistical association and does not establish causality.")
  )
}

## ---------------------------------------------------------------------------
## Quadrant plot (the hero visualization)
## ---------------------------------------------------------------------------

cx_quadrant_ggplot <- function(df, expr_thresh, meth_thresh, show_lines, highlight_genes = NULL, show_labels = FALSE) {
  df$category <- factor(as.character(df$category), levels = CX_CATEGORY_ORDER)
  df$tooltip <- sprintf(
    "Gene: %s<br>log2FC: %s (FDR %s)<br>Δβ: %s (FDR %s)<br>%s",
    df$gene, cx_fmt_num(df$log2fc), cx_fmt_num(df$expr_fdr, sci = TRUE),
    cx_fmt_num(df$dbeta), cx_fmt_num(df$meth_fdr, sci = TRUE), as.character(df$category)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = dbeta, y = log2fc, color = category, text = tooltip, key = gene)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.8) +
    ggplot2::scale_color_manual(values = CX_CATEGORY_COLORS, drop = FALSE) +
    ggplot2::labs(x = "Δ DNA methylation", y = "Transcriptomic log2 fold change", color = "Category") +
    theme_arthomix()
  if (isTRUE(show_lines)) {
    p <- p +
      ggplot2::geom_vline(xintercept = c(-meth_thresh, meth_thresh), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
      ggplot2::geom_hline(yintercept = c(-expr_thresh, expr_thresh), linetype = "dashed", color = ARTHOMIX_COLORS$axis)
  }
  label_genes <- unique(c(highlight_genes, if (isTRUE(show_labels)) {
    sig <- df[df$sig_expression %in% TRUE & df$sig_methylation %in% TRUE, , drop = FALSE]
    sig <- sig[order(-abs(ifelse(is.na(sig$log2fc), 0, sig$log2fc))), , drop = FALSE]
    utils::head(sig$gene, 15)
  }))
  if (length(label_genes) > 0) {
    lab_df <- df[df$gene %in% label_genes, , drop = FALSE]
    if (nrow(lab_df) > 0) {
      p <- p + ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(label = gene), color = ARTHOMIX_COLORS$ink,
                                          size = 3, show.legend = FALSE, max.overlaps = 30)
    }
  }
  p
}

cx_quadrant_plotly <- function(df, expr_thresh, meth_thresh, show_lines, highlight = NULL, show_labels = FALSE, source_id) {
  gg <- cx_quadrant_ggplot(df, expr_thresh, meth_thresh, show_lines, highlight, show_labels)
  plotly::ggplotly(gg, tooltip = "text", source = source_id) %>%
    plotly::layout(legend = list(orientation = "h", y = -0.2)) %>%
    plotly::event_register("plotly_click")
}

## ---------------------------------------------------------------------------
## Volcano plots
## ---------------------------------------------------------------------------

cx_volcano_expr_plot <- function(df, expr_thresh, expr_fdr_thresh, highlight_hits) {
  d <- df[!is.na(df$log2fc) & !is.na(df$expr_fdr), , drop = FALSE]
  validate(need(nrow(d) > 0, "No expression data to plot."))
  d$neglog10fdr <- -log10(pmax(d$expr_fdr, 1e-300))
  d$hit <- if (isTRUE(highlight_hits)) ifelse(d$sig_methylation %in% TRUE, "Also significant in methylation", "Expression only") else "Gene"
  cols <- if (isTRUE(highlight_hits)) c("Also significant in methylation" = ARTHOMIX_COLORS$red, "Expression only" = ARTHOMIX_COLORS$ink_muted) else c("Gene" = ARTHOMIX_COLORS$blue)
  ggplot2::ggplot(d, ggplot2::aes(x = log2fc, y = neglog10fdr, color = hit)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::geom_vline(xintercept = c(-expr_thresh, expr_thresh), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
    ggplot2::geom_hline(yintercept = -log10(expr_fdr_thresh), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(x = "log2 Fold Change", y = "-log10(FDR)", color = NULL, title = "Transcriptomics volcano") +
    theme_arthomix()
}

cx_volcano_meth_plot <- function(df, meth_thresh, meth_fdr_thresh, highlight_hits) {
  d <- df[!is.na(df$dbeta) & !is.na(df$meth_fdr), , drop = FALSE]
  validate(need(nrow(d) > 0, "No methylation data to plot."))
  d$neglog10fdr <- -log10(pmax(d$meth_fdr, 1e-300))
  d$hit <- if (isTRUE(highlight_hits)) ifelse(d$sig_expression %in% TRUE, "Also significant in expression", "Methylation only") else "Probe/gene"
  cols <- if (isTRUE(highlight_hits)) c("Also significant in expression" = ARTHOMIX_COLORS$red, "Methylation only" = ARTHOMIX_COLORS$ink_muted) else c("Probe/gene" = ARTHOMIX_COLORS$aqua)
  ggplot2::ggplot(d, ggplot2::aes(x = dbeta, y = neglog10fdr, color = hit)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::geom_vline(xintercept = c(-meth_thresh, meth_thresh), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
    ggplot2::geom_hline(yintercept = -log10(meth_fdr_thresh), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(x = "Δβ", y = "-log10(FDR)", color = NULL, title = "Methylomics volcano") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Gene-level correlation plot (spec section 16) - only ever called when
## cx_detect_sample_pairing() found real matched samples.
## ---------------------------------------------------------------------------

cx_gene_correlation_plot <- function(expr_wide, meth_wide, expr_mapping, meth_mapping, gene, common_samples, method) {
  validate(need(!is.null(expr_wide) && !is.null(meth_wide), "Raw per-sample data not available."))
  expr_gene_col <- if (!is.null(expr_mapping) && !is.na(expr_mapping[["gene"]])) expr_mapping[["gene"]] else "gene"
  meth_gene_col <- if (!is.null(meth_mapping) && !is.na(meth_mapping[["gene"]])) meth_mapping[["gene"]] else "gene"
  e_rows <- expr_wide[as.character(expr_wide[[expr_gene_col]]) == gene, common_samples, drop = FALSE]
  m_rows <- meth_wide[as.character(meth_wide[[meth_gene_col]]) == gene, common_samples, drop = FALSE]
  validate(need(nrow(e_rows) > 0 && nrow(m_rows) > 0, "This gene has no matched expression/methylation values."))
  em <- as.matrix(e_rows); storage.mode(em) <- "double"
  mm <- as.matrix(m_rows); storage.mode(mm) <- "double"
  x <- if (nrow(em) == 1) em[1, ] else colMeans(em, na.rm = TRUE)
  y <- if (nrow(mm) == 1) mm[1, ] else colMeans(mm, na.rm = TRUE)
  d <- data.frame(sample = common_samples, expression = as.numeric(x[common_samples]), methylation = as.numeric(y[common_samples]))
  d <- d[is.finite(d$expression) & is.finite(d$methylation), , drop = FALSE]
  validate(need(nrow(d) >= 3, "Not enough matched values to plot a correlation."))
  ct <- stats::cor.test(d$expression, d$methylation, method = method)
  ggplot2::ggplot(d, ggplot2::aes(x = methylation, y = expression, text = sample)) +
    ggplot2::geom_point(color = ARTHOMIX_COLORS$blue, size = 2) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = ARTHOMIX_COLORS$red, linewidth = 0.6, formula = y ~ x) +
    ggplot2::labs(x = "Methylation beta value", y = "Gene expression",
                  title = sprintf("%s: r = %.2f, p = %.3g", gene, unname(ct$estimate), ct$p.value)) +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Static gene-CpG network (spec section 20; best-effort, igraph only).
## ---------------------------------------------------------------------------
## Deliberately does NOT use ggraph::geom_edge_link() - on this deployment's
## installed ggraph/ggplot2 combination it fails on every input, even a
## trivial 3-node graph ("Problem while converting geom to grob ...
## SET_VECTOR_ELT() can only be applied to a 'list', not a 'integer'", a
## known ggraph/grid version-compatibility issue, not specific to this
## module's graphs). Computing the layout with igraph::layout_with_fr()
## directly and drawing it with plain ggplot2::geom_segment()/geom_point()
## sidesteps ggraph's edge-grob code entirely while still producing a real
## static network view.

cx_gene_cpg_network_plot <- function(sub_df) {
  edges <- do.call(rbind, lapply(seq_len(nrow(sub_df)), function(i) {
    cpgs <- strsplit(as.character(sub_df$cpg[i]), ";")[[1]]
    if (length(cpgs) == 0 || !nzchar(cpgs[1])) cpgs <- paste0(sub_df$gene[i], "_probe")
    cpgs <- utils::head(cpgs, 3)
    data.frame(from = sub_df$gene[i], to = cpgs, stringsAsFactors = FALSE)
  }))
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  vtype <- ifelse(names(igraph::V(g)) %in% sub_df$gene, "Gene", "CpG")
  set.seed(1)
  xy <- igraph::layout_with_fr(g)
  nodes <- data.frame(name = names(igraph::V(g)), type = vtype, x = xy[, 1], y = xy[, 2], stringsAsFactors = FALSE)
  el <- igraph::as_edgelist(g, names = TRUE)
  seg <- data.frame(
    x = nodes$x[match(el[, 1], nodes$name)], y = nodes$y[match(el[, 1], nodes$name)],
    xend = nodes$x[match(el[, 2], nodes$name)], yend = nodes$y[match(el[, 2], nodes$name)]
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y, color = type), size = 4) +
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x = x, y = y, label = name), size = 3, max.overlaps = 50) +
    ggplot2::scale_color_manual(values = c(Gene = ARTHOMIX_COLORS$blue, CpG = ARTHOMIX_COLORS$aqua)) +
    ggplot2::labs(color = NULL, title = "Gene - CpG association network (static)") +
    theme_arthomix() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}

## ---------------------------------------------------------------------------
## Genomic (chromosome) view - only rendered when chr is available (spec
## section 21: optional, never mandatory).
## ---------------------------------------------------------------------------

cx_genomic_view_plot <- function(df) {
  d <- df[!is.na(df$chr), , drop = FALSE]
  validate(need(nrow(d) > 0, "No chromosome information available."))
  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
  d$chr <- factor(ifelse(grepl("^chr", d$chr), d$chr, paste0("chr", d$chr)), levels = chr_levels)
  d <- d[!is.na(d$chr), , drop = FALSE]
  validate(need(nrow(d) > 0, "Chromosome values did not match the expected chr1-22/X/Y format."))
  ggplot2::ggplot(d, ggplot2::aes(x = chr, y = dbeta, color = category)) +
    ggplot2::geom_jitter(width = 0.25, alpha = 0.7, size = 1.6) +
    ggplot2::scale_color_manual(values = CX_CATEGORY_COLORS, drop = FALSE) +
    ggplot2::labs(x = "Chromosome", y = "Δβ", color = "Category", title = "Genomic view (by chromosome)") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

## ---------------------------------------------------------------------------
## Pathway enrichment (GO via clusterProfiler::enrichGO, KEGG via enrichKEGG -
## both already-installed dependencies; Reactome needs ReactomePA, which
## isn't installed, so "KEGG"/"BP"/"MF"/"CC" are the only ontology choices
## offered in the UI)
## ---------------------------------------------------------------------------

cx_run_pathway_enrichment <- function(genes, universe, ontology) {
  if (identical(ontology, "KEGG")) {
    ids <- tryCatch(clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db), error = function(e) NULL)
    if (is.null(ids) || nrow(ids) == 0) return(list(ok = FALSE, error = "Could not map any of the selected genes to Entrez IDs for KEGG."))
    uni_ids <- tryCatch(clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db), error = function(e) NULL)
    ek <- tryCatch(
      clusterProfiler::enrichKEGG(gene = ids$ENTREZID, organism = "hsa",
                                   universe = if (!is.null(uni_ids)) uni_ids$ENTREZID else NULL, pAdjustMethod = "BH"),
      error = function(e) e
    )
    if (inherits(ek, "error")) return(list(ok = FALSE, error = paste("KEGG enrichment failed (this needs internet access to the KEGG REST API):", conditionMessage(ek))))
    return(list(ok = TRUE, table = as.data.frame(ek), enrich_obj = ek))
  }
  eg <- tryCatch(
    clusterProfiler::enrichGO(gene = genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL",
                                ont = ontology, universe = universe, pAdjustMethod = "BH"),
    error = function(e) e
  )
  if (inherits(eg, "error")) return(list(ok = FALSE, error = paste("GO enrichment failed:", conditionMessage(eg))))
  list(ok = TRUE, table = as.data.frame(eg), enrich_obj = eg)
}

## ---------------------------------------------------------------------------
## Analysis report (spec section 23/24)
## ---------------------------------------------------------------------------

cx_build_report <- function(df, provenance) {
  counts <- table(df$category)
  top <- df[df$sig_expression %in% TRUE & df$sig_methylation %in% TRUE, , drop = FALSE]
  top <- top[order(-abs(ifelse(is.na(top$log2fc), 0, top$log2fc))), , drop = FALSE]
  top <- utils::head(top, 20)
  c(
    "# Cross-Omics Expression x Methylation Integration Report", "",
    sprintf("Generated: %s", format(Sys.time())), "",
    "## Parameters", "", provenance, "",
    "## Summary", "",
    sprintf("- Genes analyzed: %s", format(nrow(df), big.mark = ",")),
    sprintf("- Significant DEGs: %s", format(sum(df$sig_expression, na.rm = TRUE), big.mark = ",")),
    sprintf("- Significant DMGs: %s", format(sum(df$sig_methylation, na.rm = TRUE), big.mark = ",")),
    sprintf("- Hyper + Down (potential methylation-associated repression): %s", counts[["Hyper + Down"]] %||% 0),
    sprintf("- Hypo + Up (potential methylation-associated activation): %s", counts[["Hypo + Up"]] %||% 0),
    sprintf("- Hyper + Up (concordant-direction / noncanonical): %s", counts[["Hyper + Up"]] %||% 0),
    sprintf("- Hypo + Down (concordant-direction / noncanonical): %s", counts[["Hypo + Down"]] %||% 0),
    "",
    "## Top integrated genes (by |log2FC|, significant in both layers)", "",
    if (nrow(top) > 0) {
      sprintf("- %s: log2FC=%.2f, FDR=%.3g, Δβ=%.3f, FDR=%.3g, category=%s",
              top$gene, top$log2fc, top$expr_fdr, top$dbeta, top$meth_fdr, as.character(top$category))
    } else "(none at current thresholds)",
    "",
    "## Interpretation", "",
    "These results represent statistical associations between differential gene expression and differential DNA methylation.",
    "They do not, on their own, establish a causal regulatory relationship - functional follow-up would be required for that.",
    ""
  )
}

## ---------------------------------------------------------------------------
## Dataset validation checklist (spec section 6)
## ---------------------------------------------------------------------------

cx_validation_checklist_ui <- function(validation) {
  if (is.null(validation)) return(div(class = "empty-note", icon("circle-info"), "Load Transcriptomics and/or Methylomics data to see validation results here."))
  item <- function(chk) tags$li(
    style = sprintf("color:%s; list-style:none; padding:2px 0;", if (isTRUE(chk$ok)) "inherit" else ARTHOMIX_STATUS$critical),
    icon(if (isTRUE(chk$ok)) "circle-check" else "circle-xmark",
         style = sprintf("color:%s; margin-right:6px;", if (isTRUE(chk$ok)) ARTHOMIX_STATUS$good else ARTHOMIX_STATUS$critical)),
    chk$label
  )
  tagList(
    tags$h5("Transcriptomics"),
    tags$ul(style = "padding-left: 4px;", lapply(validation$transcriptomics, item)),
    tags$h5("Methylomics"),
    tags$ul(style = "padding-left: 4px;", lapply(validation$methylomics, item)),
    tags$h5("Cross-omics compatibility"),
    tags$ul(style = "padding-left: 4px;", lapply(validation$compatibility, item)),
    if (!isTRUE(validation$ready)) div(class = "empty-note", icon("triangle-exclamation"),
      "Cannot perform gene-level methylation integration until both a Transcriptomics and a Methylomics dataset are loaded.")
  )
}

## ---------------------------------------------------------------------------
## Sex comparison (spec section 14)
## ---------------------------------------------------------------------------

cx_sex_comparison_summary_ui <- function(cmp) {
  if (!isTRUE(cmp$ok)) return(div(class = "empty-note", icon("circle-info"), cmp$error %||% "Run Female and Male separately above (each already works), then come back here for a combined view across both sexes - this is the recommended way to see \"all samples,\" since no pooled/ALL methylation dataset exists (see Methodology & References)."))
  s <- cmp$summary
  card <- function(label, value) div(class = "card", style = "flex:1 1 140px; text-align:center; padding:10px;",
    div(style = "font-size:1.3em; font-weight:600;", if (is.na(value)) "—" else format(value, big.mark = ",")),
    div(style = "font-size:0.8em; color:var(--color-ink-muted, #898781);", label))
  tagList(
    p(class = "submodule-desc", sprintf("Comparing strata: %s. Differences below are reported as sex-stratified differences in significance, not as confirmed sex-specific biological effects (no matched-sample sex x molecular-effect interaction test was available).", paste(toupper(substring(s$strata, 1, 1)), substring(s$strata, 2), sep = "", collapse = " / "))),
    div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
        card("Genes compared", s$n_genes_compared),
        card("Shared significant", s$n_shared_significant),
        card("Female-only significant", s$n_female_specific),
        card("Male-only significant", s$n_male_specific),
        card("Directionally concordant", s$n_concordant),
        card("Directionally discordant", s$n_discordant))
  )
}

## ---------------------------------------------------------------------------
## Methodology & References (spec section 43) - static content, real
## citations only, no invented DOIs.
## ---------------------------------------------------------------------------

cx_methodology_references_ui <- function() {
  tagList(
    tags$h5("What this module does"),
    tags$p("Gene-level integration of differential expression (Transcriptomics) and differential DNA methylation (Methylomics): standardizes both inputs, harmonizes gene identifiers, maps CpGs to genes with their genomic/CpG-island context, classifies each gene into a quadrant and an Evidence Level, and - only when matched sample-level data are available - computes a methylation-expression correlation. All statements are derived from the data provided; none are generated by an AI model."),
    tags$h5("Inspired by, not a reimplementation of"),
    tags$p("The general approach of combining CpG-to-gene mapping with the observed methylation-expression relationship (rather than assuming every hypermethylated promoter silences its gene) follows the reasoning described by the MethylMix method. This module does not implement MethylMix's beta-mixture-model methylation state detection or its Wilcoxon-based methylation-driven-gene test - it uses a simpler, transparent, threshold-and-correlation-based approach appropriate for the data this application has available."),
    tags$h5("Why there is no combined-sex (ALL) methylation dataset"),
    tags$p("This is a deliberate methodological choice made upstream, not a missing file. The Methylomics pipeline fit fully independent female and male differential-methylation models specifically because a pooled or sex-adjusted analysis can obscure sex-specific methylation associations that a stratified design surfaces (Tesfaye et al. 2024) - so no pooled DMP panel was ever computed, and this module does not fabricate one by pooling the two stratified tables after the fact. Transcriptomics' ALL differential-expression table exists and works normally here; only the methylation side is affected. Run Female and Male separately (both work end-to-end) and use the \"Sex Comparison\" tab for an honest combined view across both sexes - it compares the two independently-run strata rather than pooling their underlying statistics."),
    tags$h5("References"),
    tags$ul(
      tags$li("Gevaert O (2015). \"MethylMix: an R package for identifying DNA methylation-driven genes.\" Bioinformatics, 31(11):1839-1841."),
      tags$li("Cedoz PL, Prunello M, Brennan K, Gevaert O (2018). \"MethylMix 2.0: an R package for identifying DNA methylation genes.\" Bioinformatics, 34(17):3044-3046."),
      tags$li("Bibikova M, Barnes B, Tsan C, et al. (2011). \"High density DNA methylation array with single CpG site resolution.\" Genomics, 98(4):288-295. (Illumina Infinium methylation array design and genomic-context annotation - TSS200/TSS1500/5'UTR/1st Exon/Gene Body/3'UTR and CpG island/shore/shelf definitions.)"),
      tags$li("Krassowski M, Das V, Sahu SK, Misra BB (2020). \"State of the Field in Multi-Omics Research: From Computational Needs to Data Mining and Sharing.\" Frontiers in Genetics, 11:610798."),
      tags$li("Tesfaye M, Spindola LM, Stavrum AK, Shadrin A, Melle I, Andreassen OA, Le Hellard S (2024). \"Sex effects on DNA methylation affect discovery in epigenome-wide association study of schizophrenia.\" Molecular Psychiatry, 29:2467-2477. (Cited by this app's own Methylomics pipeline as the basis for choosing sex-stratified over pooled differential-methylation models.)")
    ),
    tags$h5("Application-specific implementation"),
    tags$p("Gene identifier harmonization (HGNC symbol / Entrez ID / Ensembl Gene ID, with exact-match and alias resolution), the CpG-to-gene aggregation options, the Evidence Level classification (Strong/Moderate/Expression-only/Methylation-only/Discordant/Insufficient), and the ALL/FEMALE/MALE sex-stratified comparison are specific to this application and are not drawn from an external published method.")
  )
}
