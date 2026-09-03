## R/multiomics/04_SNF_Clustering/snf_clustering_plots.R
## Plot functions for the live "SNF Clustering" submodule
## (snf_clustering_helpers.R / mod_multi_stratification.R). Reuses

sfc_spectral_embedding <- function(W, clusters) {
  if (is.null(W) || is.null(clusters) || nrow(W) < 3) return(NULL)
  d <- rowSums(W); d[d <= 0] <- 1e-8
  Dinv <- diag(1 / sqrt(d))
  L <- Dinv %*% W %*% Dinv
  eig <- tryCatch(eigen(L, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || ncol(eig$vectors) < 2) return(NULL)
  ids <- rownames(W)
  data.frame(sample = ids, dim1 = eig$vectors[, 1], dim2 = eig$vectors[, 2], cluster = factor(clusters[ids]))
}

sfc_spectral_embedding_plot <- function(W, clusters) {
  df <- sfc_spectral_embedding(W, clusters)
  if (is.null(df)) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, color = cluster)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::scale_color_manual(values = arthomix_pair(levels(df$cluster))) +
    ggplot2::labs(x = "Spectral dimension 1", y = "Spectral dimension 2", color = "Cluster",
                  title = "Projection of the fused patient similarity network") +
    theme_arthomix()
}

sfc_feature_heatmap <- function(layers, clusters, top_n_per_block = 25) {
  if (is.null(layers) || length(layers) == 0 || is.null(clusters)) return(NULL)
  ord <- names(sort(clusters))
  rows <- do.call(rbind, lapply(names(layers), function(b) {
    m <- layers[[b]][ord, , drop = FALSE]
    v <- apply(m, 2, stats::var, na.rm = TRUE)
    top <- names(sort(v, decreasing = TRUE))[seq_len(min(top_n_per_block, ncol(m)))]
    z <- scale(m[, top, drop = FALSE])
    data.frame(
      sample = factor(rep(ord, times = length(top)), levels = ord),
      feature = factor(paste(b, rep(top, each = length(ord)), sep = ": ")),
      block = b, value = as.numeric(z)
    )
  }))
  if (is.null(rows) || nrow(rows) == 0) return(NULL)
  ggplot2::ggplot(rows, ggplot2::aes(x = sample, y = feature, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = ARTHOMIX_COLORS$blue, mid = "white", high = ARTHOMIX_COLORS$red, midpoint = 0, na.value = ARTHOMIX_COLORS$grid) +
    ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Patients (ordered by cluster)", y = NULL, fill = "Z-score") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7))
}

sfc_concordance_bar_plot <- function(conc_df) {
  need <- c("block", "concordance_with_fused")
  if (is.null(conc_df) || nrow(conc_df) == 0 || !all(need %in% colnames(conc_df))) return(NULL)
  ggplot2::ggplot(conc_df, ggplot2::aes(x = stats::reorder(block, concordance_with_fused), y = concordance_with_fused)) +
    ggplot2::geom_col(fill = ARTHOMIX_COLORS$aqua) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Concordance with fused network (NMI)", title = "Modality contribution (SNFtool::concordanceNetworkNMI)") +
    theme_arthomix()
}

sfc_stability_plot <- function(stability) {
  if (is.null(stability) || !isTRUE(stability$ok)) return(NULL)
  df <- data.frame(ari = stability$ari)
  ggplot2::ggplot(df, ggplot2::aes(x = ari)) +
    ggplot2::geom_histogram(bins = min(15, max(5, round(length(stability$ari) / 2))), fill = ARTHOMIX_COLORS$blue, color = "white") +
    ggplot2::geom_vline(xintercept = SFC_STABILITY_THRESHOLDS$stable, linetype = "dashed", color = ARTHOMIX_STATUS$good) +
    ggplot2::geom_vline(xintercept = SFC_STABILITY_THRESHOLDS$moderate, linetype = "dashed", color = ARTHOMIX_STATUS$warning) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Adjusted Rand Index (resample vs. full-cohort clustering)", y = "Resamples",
                  title = sprintf("Cluster stability across %d resamples", stability$n_resamples)) +
    theme_arthomix()
}

sfc_sensitivity_plot <- function(sens) {
  if (is.null(sens) || !isTRUE(sens$ok)) return(NULL)
  df <- sens$detail
  ggplot2::ggplot(df, ggplot2::aes(x = level, y = ari_vs_reference, group = parameter)) +
    ggplot2::geom_col(fill = ARTHOMIX_COLORS$violet) +
    ggplot2::geom_hline(yintercept = SFC_STABILITY_THRESHOLDS$stable, linetype = "dashed", color = ARTHOMIX_STATUS$good) +
    ggplot2::facet_wrap(~parameter, nrow = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "ARI vs. reference clustering", title = "Parameter sensitivity") +
    theme_arthomix()
}

sfc_feature_rank_plot <- function(df, top_n = 20) {
  need <- c("feature", "p_value", "block")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  df <- utils::head(df[order(df$p_value), , drop = FALSE], top_n)
  df$.nlp <- -log10(pmax(df$p_value, 1e-300))
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(feature, .nlp), y = .nlp, fill = block)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(df$block))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "-log10(p), Kruskal-Wallis vs. cluster", fill = "Block") +
    theme_arthomix()
}

sfc_continuous_multi_plot <- function(clusters, sample_meta, vars) {
  if (is.null(sample_meta) || length(vars) == 0) return(NULL)
  rows <- do.call(rbind, lapply(vars, function(v) {
    x <- stats::setNames(sample_meta[[v]], rownames(sample_meta))
    common <- intersect(names(clusters), names(x)); common <- common[!is.na(x[common])]
    if (length(common) == 0) return(NULL)
    data.frame(variable = v, cluster = factor(clusters[common]), value = as.numeric(x[common]))
  }))
  if (is.null(rows) || nrow(rows) == 0) return(NULL)
  ggplot2::ggplot(rows, ggplot2::aes(x = cluster, y = value, fill = cluster)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.4, size = 1) +
    ggplot2::facet_wrap(~variable, scales = "free_y") +
    ggplot2::scale_fill_manual(values = arthomix_pair(levels(rows$cluster))) +
    ggplot2::labs(x = "Cluster", y = NULL, fill = NULL) +
    theme_arthomix() + ggplot2::theme(legend.position = "none")
}

sfc_categorical_multi_plot <- function(clusters, sample_meta, vars) {
  if (is.null(sample_meta) || length(vars) == 0) return(NULL)
  rows <- do.call(rbind, lapply(vars, function(v) {
    x <- stats::setNames(sample_meta[[v]], rownames(sample_meta))
    common <- intersect(names(clusters), names(x)); common <- common[!is.na(x[common])]
    if (length(common) == 0) return(NULL)
    data.frame(variable = v, cluster = factor(clusters[common]), value = as.character(x[common]))
  }))
  if (is.null(rows) || nrow(rows) == 0) return(NULL)
  ggplot2::ggplot(rows, ggplot2::aes(x = cluster, fill = value)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::facet_wrap(~variable, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Cluster", y = "Proportion", fill = NULL) +
    theme_arthomix()
}

sfc_km_plot <- function(surv) {
  if (is.null(surv) || !isTRUE(surv$ok) || !requireNamespace("survminer", quietly = TRUE)) return(NULL)
  tryCatch(
    survminer::ggsurvplot(surv$fit, data = surv$data, pval = FALSE,
                           palette = unname(arthomix_pair(levels(surv$data$cluster))),
                           legend.title = "Cluster", xlab = "Time", ggtheme = theme_arthomix())$plot,
    error = function(e) NULL
  )
}
