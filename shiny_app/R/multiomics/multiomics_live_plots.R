## R/multiomics/multiomics_live_plots.R
## Plotting helpers for the "Live Analysis (Upload & MOFA2)" sub-module -
## same theme_arthomix()/ARTHOMIX_COLORS convention and multi_plot_or_empty()
## empty-state wrapper as the rest of the Multi-Omics module (see
## multiomics_plots.R). Every plot here is drawn from a real, live
## calculation (multiomics_live_helpers.R) over data the user uploaded -
## never a fake/placeholder value.

## ---------------------------------------------------------------------------
## Missingness (spec Plots 1-3)
## ---------------------------------------------------------------------------

multi_live_missingness_by_omics_plot <- function(validations) {
  validations <- Filter(function(v) isTRUE(v$ok), validations)
  if (length(validations) == 0) return(NULL)
  df <- data.frame(omics = vapply(validations, function(v) v$layer, character(1)),
                    pct_missing = vapply(validations, function(v) v$pct_missing, numeric(1)))
  ggplot2::ggplot(df, ggplot2::aes(x = omics, y = pct_missing, fill = omics)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = arthomix_pair(df$omics)) +
    ggplot2::labs(x = NULL, y = "Missing (%)", title = "Missingness by omics layer") +
    theme_arthomix() + ggplot2::theme(legend.position = "none")
}

multi_live_sample_missingness_plot <- function(miss, threshold = NULL) {
  if (is.null(miss) || is.null(miss$per_sample) || nrow(miss$per_sample) == 0) return(NULL)
  df <- miss$per_sample
  df$.over <- if (!is.null(threshold)) df$pct_missing > threshold else FALSE
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(sample, pct_missing), y = pct_missing, fill = .over)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = ARTHOMIX_COLORS$red, `FALSE` = ARTHOMIX_COLORS$blue), guide = "none") +
    (if (!is.null(threshold)) ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) else NULL) +
    ggplot2::labs(x = "Sample", y = "Missing (%)", title = "Sample-level missingness") +
    theme_arthomix() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

multi_live_feature_missingness_plot <- function(miss) {
  if (is.null(miss) || is.null(miss$per_feature) || nrow(miss$per_feature) == 0) return(NULL)
  ggplot2::ggplot(miss$per_feature, ggplot2::aes(x = pct_missing)) +
    ggplot2::geom_histogram(bins = 30, fill = ARTHOMIX_COLORS$blue, color = "white") +
    ggplot2::labs(x = "Missing (%)", y = "Features", title = "Feature-level missingness distribution") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Before/after normalization distributions (spec Plots 4-7)
## ---------------------------------------------------------------------------

multi_live_distribution_plot <- function(mat, kind = c("box", "density"), max_samples = 40) {
  kind <- match.arg(kind)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  m <- mat
  if (nrow(m) > max_samples) m <- m[sample(seq_len(nrow(m)), max_samples), , drop = FALSE]
  long <- data.frame(sample = rep(rownames(m), ncol(m)), value = as.vector(m))
  if (identical(kind, "box")) {
    ggplot2::ggplot(long, ggplot2::aes(x = sample, y = value)) +
      ggplot2::geom_boxplot(fill = ARTHOMIX_COLORS$blue, alpha = 0.6, outlier.size = 0.5) +
      ggplot2::labs(x = "Sample", y = "Value") + theme_arthomix() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  } else {
    ggplot2::ggplot(long, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = ARTHOMIX_COLORS$blue, alpha = 0.5) +
      ggplot2::labs(x = "Value", y = "Density") + theme_arthomix()
  }
}

## ---------------------------------------------------------------------------
## Feature filtering retention (spec section 10)
## ---------------------------------------------------------------------------

multi_live_retention_plot <- function(n_before, n_after) {
  if (is.null(n_before) || is.null(n_after)) return(NULL)
  df <- data.frame(stage = factor(c("Before filtering", "After filtering"), levels = c("Before filtering", "After filtering")), n = c(n_before, n_after))
  ggplot2::ggplot(df, ggplot2::aes(x = stage, y = n, fill = stage)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = format(n, big.mark = ",")), vjust = -0.4) +
    ggplot2::scale_fill_manual(values = c("Before filtering" = ARTHOMIX_COLORS$ink_muted, "After filtering" = ARTHOMIX_COLORS$aqua), guide = "none") +
    ggplot2::labs(x = NULL, y = "Features") + theme_arthomix()
}

## ---------------------------------------------------------------------------
## Cross-omics scale comparison (spec Plots 12-13) - per-layer value-range
## boxplot before vs. after scaling, so a reader can see layers on wildly
## different native scales get pulled onto a common one.
## ---------------------------------------------------------------------------

multi_live_scale_comparison_plot <- function(mat_list, labels) {
  mat_list <- Filter(Negate(is.null), mat_list)
  if (length(mat_list) == 0) return(NULL)
  long <- do.call(rbind, lapply(seq_along(mat_list), function(i) {
    data.frame(layer = labels[i], value = as.vector(mat_list[[i]]))
  }))
  ggplot2::ggplot(long, ggplot2::aes(x = layer, y = value, fill = layer)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.4) +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(long$layer))) +
    ggplot2::labs(x = NULL, y = "Value") + theme_arthomix() + ggplot2::theme(legend.position = "none")
}

## ---------------------------------------------------------------------------
## PCA (spec Plots 8-9, batch diagnostics)
## ---------------------------------------------------------------------------

multi_live_pca_plot <- function(pca, meta = NULL, color_by = NULL) {
  if (is.null(pca) || !isTRUE(pca$ok)) return(NULL)
  df <- pca$scores
  df$.id <- pca$sample_ids
  if (!is.null(meta) && !is.null(color_by) && color_by %in% colnames(meta)) {
    df$.color <- meta[df$.id, color_by]
  } else df$.color <- "all samples"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = factor(.color))) +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(df$.color))) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * pca$var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * pca$var_explained[2]),
      color = color_by %||% NULL
    ) +
    theme_arthomix()
  if (is.null(color_by)) p <- p + ggplot2::theme(legend.position = "none")
  p
}

## ---------------------------------------------------------------------------
## MOFA2 factor results (spec Plots 14-19)
## ---------------------------------------------------------------------------

multi_live_mofa_variance_plot <- function(var_df) {
  need <- c("view", "factor", "variance_explained")
  if (is.null(var_df) || nrow(var_df) == 0 || !all(need %in% colnames(var_df))) return(NULL)
  ggplot2::ggplot(var_df, ggplot2::aes(x = factor, y = variance_explained, fill = view)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(var_df$view))) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Factor", y = "Variance explained", fill = "View (omics)") +
    theme_arthomix()
}

multi_live_factor_score_plot <- function(factors_df, x_factor, y_factor, meta = NULL, color_by = NULL) {
  if (is.null(factors_df) || !all(c(x_factor, y_factor) %in% colnames(factors_df))) return(NULL)
  df <- factors_df
  df$.id <- rownames(df)
  if (!is.null(meta) && !is.null(color_by) && color_by %in% colnames(meta)) {
    df$.color <- meta[df$.id, color_by]
  } else df$.color <- "all samples"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_factor]], y = .data[[y_factor]], color = factor(.color))) +
    ggplot2::geom_hline(yintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_vline(xintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(df$.color))) +
    ggplot2::labs(x = x_factor, y = y_factor, color = color_by %||% NULL) +
    theme_arthomix()
  if (is.null(color_by)) p <- p + ggplot2::theme(legend.position = "none")
  p
}

multi_live_factor_heatmap <- function(factors_df) {
  if (is.null(factors_df) || nrow(factors_df) == 0) return(NULL)
  long <- data.frame(sample = rep(rownames(factors_df), ncol(factors_df)),
                      factor = rep(colnames(factors_df), each = nrow(factors_df)),
                      score = as.vector(as.matrix(factors_df)))
  ggplot2::ggplot(long, ggplot2::aes(x = factor, y = sample, fill = score)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = ARTHOMIX_COLORS$red, mid = "white", high = ARTHOMIX_COLORS$blue) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Factor score") +
    theme_arthomix() + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
}

multi_live_loadings_plot <- function(loadings_df, sign = c("both", "positive", "negative"), top_n = 20) {
  sign <- match.arg(sign)
  need <- c("feature", "value", "view")
  if (is.null(loadings_df) || nrow(loadings_df) == 0 || !all(need %in% colnames(loadings_df))) return(NULL)
  df <- loadings_df
  if (identical(sign, "positive")) df <- df[df$value > 0, , drop = FALSE]
  if (identical(sign, "negative")) df <- df[df$value < 0, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df <- df[order(-abs(df$value)), , drop = FALSE]
  df <- utils::head(df, top_n)
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(feature, value), y = value, fill = view)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(df$view))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Loading", fill = "View (omics)") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Cross-omics correlation (spec Plots 20-21)
## ---------------------------------------------------------------------------

multi_live_correlation_scatter_plot <- function(x, y, xlab, ylab, r, p, n) {
  if (is.null(x) || is.null(y)) return(NULL)
  df <- data.frame(x = x, y = y)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = ARTHOMIX_COLORS$blue, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = ARTHOMIX_COLORS$red, linewidth = 0.6) +
    ggplot2::labs(x = xlab, y = ylab, title = sprintf("r = %.3f, p = %.3g, n = %d", r, p, n)) +
    theme_arthomix()
}

multi_live_correlation_heatmap_plot <- function(cor_df, fdr_threshold = NULL) {
  need <- c("featureA", "featureB", "r")
  if (is.null(cor_df) || nrow(cor_df) == 0 || !all(need %in% colnames(cor_df))) return(NULL)
  df <- cor_df
  if (!is.null(fdr_threshold) && "fdr" %in% colnames(df)) df$r[df$fdr > fdr_threshold] <- NA
  ggplot2::ggplot(df, ggplot2::aes(x = featureB, y = featureA, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = ARTHOMIX_COLORS$red, mid = "white", high = ARTHOMIX_COLORS$blue, na.value = "grey90", limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "r") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7), axis.text.y = ggplot2::element_text(size = 7))
}
