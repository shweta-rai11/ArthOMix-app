## R/crossomics/functions/integration/crossomics_integration_plots.R
## Plotting/rendering helpers for mod_cross_integration.R - kept separate from
## the module file so the Shiny wiring (mod_cross_integration.R) and the
## visualization logic can be read/changed independently. Every plot reuses
## the app's shared theme_arthomix()/ARTHOMIX_COLORS (global.R) so it matches
## the rest of the app, and the quadrant plot is built once as a ggplot
## object (cx_quadrant_ggplot) reused for both the interactive plotly view
## and PNG/PDF/SVG export, so there is exactly one place that draws it.

## `message` defaults to the Integration module's own instruction (its
## original hardcoded text) - callers in the OTHER Cross-Omics sub-modules
## (Biomarker Convergence, Cross-Omics MR) must pass their own, since "Run
## Integration" is never the right action there.
cx_empty_state <- function(message = "Click \"Run Integration\" in the Integration tab to see results here.") {
  div(class = "empty-note", icon("circle-info"), message)
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
## Analysis report (spec section 23/24)
## ---------------------------------------------------------------------------

cx_build_report <- function(df, provenance) {
  counts <- table(df$category)
  top <- df[df$sig_expression %in% TRUE & df$sig_methylation %in% TRUE, , drop = FALSE]
  top <- top[order(-abs(ifelse(is.na(top$log2fc), 0, top$log2fc))), , drop = FALSE]
  top <- utils::head(top, 20)
  c(
    "# Cross-Omics Expression and Methylation Integration Report", "",
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

