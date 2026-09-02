## R/multiomics/07_Pathways/multiomics_pathway_plots.R
## The 6 required Pathways-tab plots. Same split/conventions as
## multiomics_plots.R: every plot reuses ARTHOMIX_COLORS/theme_arthomix() and
## is built once as a plot object (or, for the two live pathway-map fetches,
## a saved PNG path) reused for both on-screen rendering and download.

## ---------------------------------------------------------------------------
## Plot 1 - Enrichment dot plot. x = GeneRatio (ORA) or NES (GSEA), y =
## pathway (wrapped label), size = gene count, color = adjusted P/FDR.
## ---------------------------------------------------------------------------

mp_dot_plot <- function(df, method = c("ORA", "GSEA"), top_n = 20) {
  method <- match.arg(method)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df[order(df$p.adjust), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$.label <- vapply(df$Description, function(s) paste(strwrap(s, width = 45), collapse = "\n"), character(1))
  x_var <- if (identical(method, "GSEA")) "NES" else "gene_ratio_numeric"
  x_lab <- if (identical(method, "GSEA")) "Normalized Enrichment Score" else "Gene ratio"
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = stats::reorder(.label, .data[[x_var]]), size = Count, color = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = ARTHOMIX_COLORS$red, high = ARTHOMIX_COLORS$blue) +
    ggplot2::labs(x = x_lab, y = NULL, size = "Gene count", color = "Adj. P / FDR") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Plot 2 - Enrichment bar plot, top pathways, sortable.
## ---------------------------------------------------------------------------

mp_bar_plot <- function(df, top_n = 20, sort_by = c("FDR", "p.adjust", "NES", "Count")) {
  sort_by <- match.arg(sort_by)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  sort_col <- switch(sort_by, "FDR" = "qvalue", "p.adjust" = "p.adjust", "NES" = "NES", "Count" = "Count")
  if (!sort_col %in% colnames(df) || all(is.na(df[[sort_col]]))) sort_col <- "p.adjust"
  decreasing <- sort_col %in% c("NES", "Count")
  df <- df[order(if (decreasing) -abs(df[[sort_col]]) else df[[sort_col]]), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$.label <- vapply(df$Description, function(s) paste(strwrap(s, width = 45), collapse = "\n"), character(1))
  fill_col <- if ("source" %in% colnames(df)) "source" else NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(.label, -log10(pmax(p.adjust, 1e-300))), y = -log10(pmax(p.adjust, 1e-300))))
  p <- if (!is.null(fill_col)) p + ggplot2::geom_col(ggplot2::aes(fill = .data[[fill_col]])) else p + ggplot2::geom_col(fill = ARTHOMIX_COLORS$blue)
  p + ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "-log10(adjusted P)", fill = "Database") + theme_arthomix()
}

## ---------------------------------------------------------------------------
## Plot 3 - Pathway x Omics heatmap. Cell = -log10(FDR) computed for that
## omics layer's own evidence track; genuinely-not-calculated cells are an
## explicit grey tile (NA-aware fill scale), never a fabricated 0.
## ---------------------------------------------------------------------------

mp_omics_heatmap <- function(evidence_df, top_n = 25) {
  if (is.null(evidence_df) || nrow(evidence_df) == 0) return(NULL)
  df <- utils::head(evidence_df[order(evidence_df$p.adjust), , drop = FALSE], top_n)
  long <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    label <- paste(strwrap(df$Description[i], width = 40), collapse = "\n")
    rbind(
      data.frame(pathway = label, track = "Transcriptomics", value = if (!is.na(df$transcript_min_p[i])) -log10(pmax(df$transcript_min_p[i], 1e-300)) else NA_real_),
      data.frame(pathway = label, track = "Methylomics", value = if (!is.na(df$meth_min_p[i])) -log10(pmax(df$meth_min_p[i], 1e-300)) else NA_real_),
      data.frame(pathway = label, track = "Integrated", value = -log10(pmax(df$p.adjust[i], 1e-300)))
    )
  }))
  long$track <- factor(long$track, levels = c("Transcriptomics", "Methylomics", "Integrated"))
  ggplot2::ggplot(long, ggplot2::aes(x = track, y = pathway, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = ARTHOMIX_COLORS$grid, high = ARTHOMIX_COLORS$blue, na.value = "grey85") +
    ggplot2::labs(x = NULL, y = NULL, fill = "-log10(P)") +
    theme_arthomix() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

## ---------------------------------------------------------------------------
## Plot 4 - Gene-Pathway network. Plain igraph::layout_with_fr() + manual
## geom_segment/geom_point (matches cx_gene_cpg_network_plot's own precedent
## in crossomics_integration_plots.R - ggraph is installed but not used
## anywhere real in this app), node color = Transcriptomic/Methylomic/Both.
## ---------------------------------------------------------------------------

mp_gene_pathway_network <- function(enrichment_df, input_df, selected_pathway_ids = NULL, max_genes = 40) {
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)
  df <- if (!is.null(selected_pathway_ids) && length(selected_pathway_ids) > 0) enrichment_df[enrichment_df$ID %in% selected_pathway_ids, , drop = FALSE] else utils::head(enrichment_df[order(enrichment_df$p.adjust), , drop = FALSE], 8)
  if (nrow(df) == 0) return(NULL)

  edges <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    genes <- unique(trimws(unlist(strsplit(df$geneID[i] %||% "", "/"))))
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0) return(NULL)
    data.frame(pathway = df$Description[i], gene = genes, stringsAsFactors = FALSE)
  }))
  if (is.null(edges) || nrow(edges) == 0) return(NULL)
  top_genes <- names(sort(table(edges$gene), decreasing = TRUE))[seq_len(min(max_genes, length(unique(edges$gene))))]
  edges <- edges[edges$gene %in% top_genes, , drop = FALSE]

  has_expr <- toupper(input_df$gene_symbol) %in% toupper(edges$gene) & !is.na(input_df$expr_logFC)
  has_meth <- toupper(input_df$gene_symbol) %in% toupper(edges$gene) & !is.na(input_df$cpg) & nzchar(input_df$cpg)
  expr_genes <- toupper(unique(input_df$gene_symbol[has_expr]))
  meth_genes <- toupper(unique(input_df$gene_symbol[has_meth]))
  gene_evidence <- ifelse(toupper(unique(edges$gene)) %in% expr_genes & toupper(unique(edges$gene)) %in% meth_genes, "Both",
                    ifelse(toupper(unique(edges$gene)) %in% expr_genes, "Transcriptomic",
                    ifelse(toupper(unique(edges$gene)) %in% meth_genes, "Methylomic", "Unmapped")))

  nodes <- data.frame(name = c(unique(df$Description), unique(edges$gene)),
                       type = c(rep("Pathway", length(unique(df$Description))), rep("Gene", length(unique(edges$gene)))),
                       evidence = c(rep("Pathway", length(unique(df$Description))), gene_evidence), stringsAsFactors = FALSE)
  g <- tryCatch(igraph::graph_from_data_frame(edges[, c("pathway", "gene")], vertices = nodes, directed = FALSE), error = function(e) NULL)
  if (is.null(g) || igraph::vcount(g) == 0) return(NULL)
  xy <- igraph::layout_with_fr(g)
  colnames(xy) <- c("x", "y")
  node_df <- cbind(nodes, xy)
  edge_df <- edges
  edge_df$x <- node_df$x[match(edge_df$pathway, node_df$name)]; edge_df$y <- node_df$y[match(edge_df$pathway, node_df$name)]
  edge_df$xend <- node_df$x[match(edge_df$gene, node_df$name)]; edge_df$yend <- node_df$y[match(edge_df$gene, node_df$name)]

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = ARTHOMIX_COLORS$grid, linewidth = 0.3) +
    ggplot2::geom_point(data = node_df, ggplot2::aes(x = x, y = y, color = evidence, shape = type), size = 3) +
    ggplot2::geom_text(data = node_df[node_df$type == "Pathway", ], ggplot2::aes(x = x, y = y, label = name), size = 2.6, vjust = -1, color = ARTHOMIX_COLORS$ink) +
    ggplot2::scale_color_manual(values = c(Pathway = ARTHOMIX_COLORS$ink_secondary, Transcriptomic = ARTHOMIX_COLORS$blue,
                                            Methylomic = ARTHOMIX_COLORS$orange, Both = ARTHOMIX_COLORS$violet, Unmapped = ARTHOMIX_COLORS$ink_muted)) +
    ggplot2::labs(color = "Evidence", shape = NULL) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "bottom")
}

## ---------------------------------------------------------------------------
## Plot 5 - KEGG pathway map. Real pathview render, real KEGG pathway
## structure + gene-level overlay; NULL + reason on failure, never an
## invented diagram.
## ---------------------------------------------------------------------------

mp_kegg_pathway_map <- function(pathway_id, gene_entrez_effect_vec, out_dir = tempdir()) {
  if (!MP_PATHVIEW_AVAILABLE) return(list(ok = FALSE, path = NULL, error = "The pathview package is not installed in this deployment."))
  ## pathview's own internal species.info()/data(bods) lookup only resolves
  ## when the package is on the search path, not merely namespace-loaded via
  ## `::` - confirmed by testing (pathview::pathview() alone throws "object
  ## 'bods' not found"; attaching the namespace first fixes it). This is the
  ## one function in this app that needs its dependency attached rather than
  ## called by `::` throughout, and only for the duration of this call.
  if (!"package:pathview" %in% search()) suppressPackageStartupMessages(attachNamespace("pathview"))
  kegg_id <- sub("^path:", "", pathway_id)
  kegg_id <- sub("^hsa", "", kegg_id)
  suffix <- paste0("mp_", format(Sys.time(), "%H%M%OS3"))
  res <- tryCatch(
    withr_wd_pathview(out_dir, kegg_id, gene_entrez_effect_vec, suffix),
    error = function(e) e
  )
  if (inherits(res, "error")) return(list(ok = FALSE, path = NULL, error = paste("pathview render failed (this needs internet access to the KEGG API):", conditionMessage(res))))
  png_path <- file.path(out_dir, sprintf("hsa%s.%s.png", kegg_id, suffix))
  if (!file.exists(png_path)) return(list(ok = FALSE, path = NULL, error = "pathview did not produce an output image for this pathway."))
  list(ok = TRUE, path = png_path, error = NULL)
}

## pathview::pathview() always writes into the current working directory, so
## this switches into `out_dir` for the call only, restoring the original
## wd even on error - never leaves the Shiny process's cwd changed.
withr_wd_pathview <- function(out_dir, kegg_id, gene_data, suffix) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)
  pathview::pathview(gene.data = gene_data, pathway.id = kegg_id, species = "hsa", kegg.native = TRUE, out.suffix = suffix)
}

## ---------------------------------------------------------------------------
## Plot 6 - Reactome pathway visualization. Live diagram image from the
## official Reactome ContentService exporter; NULL + reason on a pathway
## that has no diagram, never a placeholder image.
## ---------------------------------------------------------------------------

mp_fetch_reactome_diagram_png <- function(stable_id, out_dir = tempdir()) {
  if (!requireNamespace("httr", quietly = TRUE)) return(list(ok = FALSE, path = NULL, error = "The httr package is not installed in this deployment."))
  res <- tryCatch(httr::GET(sprintf("https://reactome.org/ContentService/exporter/diagram/%s.png?quality=7", stable_id), httr::timeout(10)), error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, path = NULL, error = paste("Could not reach the Reactome ContentService:", conditionMessage(res))))
  if (httr::status_code(res) != 200) return(list(ok = FALSE, path = NULL, error = sprintf("No diagram is available for %s on Reactome (HTTP %s).", stable_id, httr::status_code(res))))
  png_path <- file.path(out_dir, sprintf("reactome_%s.png", stable_id))
  writeBin(httr::content(res, as = "raw"), png_path)
  list(ok = TRUE, path = png_path, error = NULL)
}
