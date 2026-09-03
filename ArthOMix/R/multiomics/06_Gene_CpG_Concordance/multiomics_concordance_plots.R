## R/multiomics/06_Gene_CpG_Concordance/multiomics_concordance_plots.R
## Plots 1-6 for the "Gene-CpG Concordance" submodule (spec section 20).
## Every function returns NULL (never a placeholder/fake plot) when its

mcc_plot_scatter <- function(df, color_by = "region_fine", meth_value = c("dbeta", "delta_beta")) {
  meth_value <- match.arg(meth_value)
  need <- c("log2fc", meth_value, color_by)
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  d <- df[!is.na(df$log2fc) & !is.na(df[[meth_value]]), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d[[color_by]][is.na(d[[color_by]])] <- "Unknown"
  ylab <- if (meth_value == "dbeta") "Methylation change (delta-M)" else "Methylation change (delta-Beta)"
  ggplot2::ggplot(d, ggplot2::aes(x = log2fc, y = .data[[meth_value]], color = .data[[color_by]])) +
    ggplot2::geom_hline(yintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_vline(xintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_point(alpha = 0.65, size = 1.7) +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(d[[color_by]]))) +
    ggplot2::labs(x = "Expression log2FC", y = ylab, color = color_by) +
    theme_arthomix()
}

mcc_plot_pair_correlation <- function(x, y, gene, cpg, r, p, fdr, n, xlab = "Methylation", ylab = "Expression") {
  if (length(x) == 0 || length(y) == 0 || length(x) != length(y) || n < 3) return(NULL)
  d <- data.frame(x = x, y = y)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = ARTHOMIX_COLORS$blue, alpha = 0.75, size = 2) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = ARTHOMIX_COLORS$red, linewidth = 0.6, formula = y ~ x) +
    ggplot2::labs(x = xlab, y = ylab,
                  title = sprintf("%s x %s: r = %.2f, p = %.3g, FDR = %.3g, n = %d", gene, cpg, r, p, fdr, n)) +
    theme_arthomix()
}

mcc_plot_quadrant <- function(df, meth_value = c("dbeta", "delta_beta"), top_label_n = 15) {
  meth_value <- match.arg(meth_value)
  need <- c("log2fc", meth_value, "gene_symbol")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  d <- df[!is.na(df$log2fc) & !is.na(df[[meth_value]]), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$quadrant <- ifelse(d$log2fc > 0 & d[[meth_value]] < 0, "Up + Hypo",
                 ifelse(d$log2fc < 0 & d[[meth_value]] > 0, "Down + Hyper",
                 ifelse(d$log2fc > 0 & d[[meth_value]] > 0, "Up + Hyper", "Down + Hypo")))
  d$rank_metric <- if ("priority_score" %in% colnames(d)) d$priority_score else abs(d$log2fc) * abs(d[[meth_value]])
  lbl <- d[order(-d$rank_metric), , drop = FALSE]
  lbl <- utils::head(lbl, top_label_n)
  ggplot2::ggplot(d, ggplot2::aes(x = log2fc, y = .data[[meth_value]], color = quadrant)) +
    ggplot2::geom_hline(yintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_vline(xintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_point(alpha = 0.6, size = 1.6) +
    ggrepel::geom_text_repel(data = lbl, ggplot2::aes(label = gene_symbol), size = 2.8, max.overlaps = 30, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(d$quadrant))) +
    ggplot2::labs(x = "Expression log2FC", y = if (meth_value == "dbeta") "Methylation change (delta-M)" else "Methylation change (delta-Beta)", color = "Quadrant") +
    theme_arthomix()
}

mcc_plot_location <- function(df) {
  need <- c("chr", "pos", "gene_symbol", "cpg")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  d <- df[!is.na(df$chr) & !is.na(df$pos), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$chr <- factor(d$chr, levels = paste0("chr", c(1:22, "X", "Y")))
  ggplot2::ggplot(d, ggplot2::aes(x = pos, y = gene_symbol, color = region_fine %||% "Unknown")) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~chr, scales = "free_x") +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(d$region_fine %||% "Unknown"))) +
    ggplot2::labs(x = "Genomic position", y = "Gene", color = "Region") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

mcc_plot_evidence_heatmap <- function(df, top_n = 30) {
  need <- c("gene_symbol", "cpg", "log2fc", "dbeta")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  d <- df
  d$rank_metric <- if ("priority_score" %in% colnames(d)) d$priority_score else abs(d$log2fc) * abs(d$dbeta)
  d <- utils::head(d[order(-d$rank_metric), , drop = FALSE], top_n)
  if (nrow(d) == 0) return(NULL)
  d$label <- paste(d$gene_symbol, d$cpg, sep = " / ")
  ev <- data.frame(
    label = rep(d$label, 5),
    metric = rep(c("Expression |log2FC|", "Methylation |delta-M|", "Correlation |r|", "-log10(FDR)", "Priority score"), each = nrow(d)),
    value = c(abs(d$log2fc), abs(d$dbeta),
              abs(d$correlation_r %||% rep(NA_real_, nrow(d))),
              -log10(pmin(d$expr_fdr %||% NA_real_, d$meth_fdr %||% NA_real_, na.rm = FALSE)),
              d$priority_score %||% rep(NA_real_, nrow(d))),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(ev, ggplot2::aes(x = metric, y = label, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "#f5f4ee", high = ARTHOMIX_COLORS$blue, na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Value") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

mcc_plot_network <- function(df, max_edges = 150) {
  need <- c("gene_symbol", "cpg")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  d <- df[!is.na(df$gene_symbol) & !is.na(df$cpg), , drop = FALSE]
  if (nrow(d) == 0 || nrow(d) > max_edges) return(NULL)
  edges <- data.frame(from = d$gene_symbol, to = d$cpg, stringsAsFactors = FALSE)
  edges <- unique(edges)
  if (nrow(edges) < 2) return(NULL)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  vnames <- names(igraph::V(g))
  vtype <- ifelse(vnames %in% d$gene_symbol, "Gene", "CpG")
  set.seed(1)
  xy <- igraph::layout_with_fr(g)
  nodes <- data.frame(name = vnames, type = vtype, x = xy[, 1], y = xy[, 2], stringsAsFactors = FALSE)
  el <- igraph::as_edgelist(g, names = TRUE)
  seg <- data.frame(
    x = nodes$x[match(el[, 1], nodes$name)], y = nodes$y[match(el[, 1], nodes$name)],
    xend = nodes$x[match(el[, 2], nodes$name)], yend = nodes$y[match(el[, 2], nodes$name)]
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y, color = type), size = 3.5) +
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x = x, y = y, label = name), size = 2.6, max.overlaps = 50) +
    ggplot2::scale_color_manual(values = c(Gene = ARTHOMIX_COLORS$blue, CpG = ARTHOMIX_COLORS$aqua)) +
    ggplot2::labs(color = NULL, title = sprintf("%d gene(s), %d CpG(s), %d edge(s)", sum(vtype == "Gene"), sum(vtype == "CpG"), nrow(edges))) +
    theme_arthomix() +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}
