## R/multiomics/multiomics_integration_plots.R
## Plot functions for the live DIABLO/SNF/Compare engine
## (multiomics_integration_helpers.R / mod_multi_integration.R). Reuses
## the shared theme_arthomix()/ARTHOMIX_COLORS/arthomix_pair() (global.R),
## multi_empty_state()/multi_plot_or_empty()/multi_png_download()
## (multiomics_plots.R), and - where the data shape already matches -
## multi_diablo_score_plot()/multi_diablo_panel_plot()/
## multi_diablo_variance_plot()/multi_live_pca_plot()/
## multi_live_correlation_heatmap_plot() directly, rather than duplicating
## them. Only genuinely new chart types live here.

## ---------------------------------------------------------------------------
## DIABLO performance: overall + per-class error rate at the fitted model's
## final component (spec section 16's "Performance" card set) - a plain bar
## chart, no fabricated confidence interval (perf() doesn't bootstrap one).
## ---------------------------------------------------------------------------

mi_diablo_error_bar_plot <- function(perf_summary) {
  if (is.null(perf_summary)) return(NULL)
  df <- data.frame(
    metric = c("Overall error", "Overall BER", names(perf_summary$per_class_error)),
    value = c(perf_summary$overall_error, perf_summary$ber, as.numeric(perf_summary$per_class_error)),
    kind = c("Overall", "Overall", rep("Per-class", length(perf_summary$per_class_error)))
  )
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(metric, value), y = value, fill = kind)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
    ggplot2::scale_fill_manual(values = c(Overall = ARTHOMIX_COLORS$blue, `Per-class` = ARTHOMIX_COLORS$aqua)) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = sprintf("Error rate (%s, comp %d, %s)", "cross-validated", perf_summary$ncomp, perf_summary$distance), fill = NULL) +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## SNF fused-network affinity heatmap, samples ordered by cluster - a plain
## visualization of the fused matrix itself, not a claim about which pairs
## are "truly" similar beyond what the matrix already encodes.
## ---------------------------------------------------------------------------

mi_snf_fused_heatmap <- function(W, clusters) {
  if (is.null(W) || is.null(clusters)) return(NULL)
  ord <- names(sort(clusters))
  ord <- intersect(ord, rownames(W))
  if (length(ord) < 2) return(NULL)
  Wo <- W[ord, ord]
  long <- data.frame(
    row = factor(rep(ord, times = length(ord)), levels = ord),
    col = factor(rep(ord, each = length(ord)), levels = ord),
    value = as.numeric(Wo)
  )
  ggplot2::ggplot(long, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = ARTHOMIX_COLORS$blue) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Fused\naffinity") +
    theme_arthomix() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}

## ---------------------------------------------------------------------------
## SNF cluster-number diagnostic - real eigengap/rotation-cost candidates
## from SNFtool::estimateNumberOfClustersGivenGraph(), not an assumed k=2
## (spec section 24).
## ---------------------------------------------------------------------------

mi_snf_cluster_estimate_plot <- function(est) {
  if (is.null(est)) return(NULL)
  df <- data.frame(criterion = names(est), k = as.integer(unlist(est)))
  ggplot2::ggplot(df, ggplot2::aes(x = criterion, y = k)) +
    ggplot2::geom_col(fill = ARTHOMIX_COLORS$aqua) +
    ggplot2::geom_text(ggplot2::aes(label = k), vjust = -0.4, color = ARTHOMIX_COLORS$ink) +
    ggplot2::labs(x = NULL, y = "Candidate number of clusters") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

## ---------------------------------------------------------------------------
## SNF cluster low-dimensional view - PCA over the column-concatenated,
## per-block standardized matrices (real PCA via the existing
## multi_live_pca()/multi_live_pca_plot(), multiomics_dataset_helpers.R /
## multiomics_dataset_plots.R), colored by fused cluster. Not a claim that
## cluster separation in this 2D view equals cluster quality - the
## eigengap/rotation-cost/silhouette diagnostics are the actual evidence.
## ---------------------------------------------------------------------------

mi_snf_pca_cluster_plot <- function(layers, clusters) {
  if (is.null(layers) || is.null(clusters)) return(NULL)
  scaled <- lapply(layers, function(m) scale(m[names(clusters), , drop = FALSE]))
  combined <- do.call(cbind, scaled)
  pca <- multi_live_pca(combined)
  if (!isTRUE(pca$ok)) return(NULL)
  meta <- data.frame(row.names = names(clusters), cluster = factor(clusters))
  multi_live_pca_plot(pca, meta, "cluster")
}

## ---------------------------------------------------------------------------
## Compare subtab: single-omics vs. integrated performance/clustering -
## a plain bar chart, "type" (Single-omics / Integrated / Baseline) sets
## the fill so integration is never visually implied to be the default winner.
## ---------------------------------------------------------------------------

mi_compare_bar_plot <- function(df, value_col = "auroc", ylab = "AUROC") {
  need <- c("model", value_col, "type")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  df$.value <- df[[value_col]]
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(model, .value), y = .value, fill = type)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
    ggplot2::scale_fill_manual(values = c(`Single-omics` = ARTHOMIX_COLORS$blue, Integrated = ARTHOMIX_COLORS$aqua, Baseline = ARTHOMIX_COLORS$ink_muted)) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = ylab, fill = NULL) +
    theme_arthomix()
}
