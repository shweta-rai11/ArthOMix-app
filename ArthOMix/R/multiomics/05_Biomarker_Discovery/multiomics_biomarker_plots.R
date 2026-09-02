## R/multiomics/05_Biomarker_Discovery/multiomics_biomarker_plots.R
## Plot functions specific to Biomarker Discovery
## (multiomics_biomarker_helpers.R / mod_multi_biomarker.R). Reuses the
## shared theme_arthomix()/ARTHOMIX_COLORS/arthomix_pair() (global.R),
## multi_empty_state()/multi_plot_or_empty()/multi_png_download()
## (multiomics_plots.R), and - where the data shape already matches -
## multi_diablo_score_plot()/multi_diablo_panel_plot()/
## multi_diablo_variance_plot()/multi_live_correlation_heatmap_plot()
## directly (spec section 13.E's sample/loading/variance/cross-omics-
## relationship plots). Only genuinely new chart types for this submodule
## live here: feature-selection stability, the selected-feature x sample
## heatmap, the between-block component-correlation scatter, and the
## pooled cross-validated ROC curve.

## ---------------------------------------------------------------------------
## Feature-selection stability (spec section 13.C) - one bar per selected
## feature, colored by its evidence-based stability category
## (mb_stability_category(), multiomics_biomarker_helpers.R), faceted by
## omics block so transcriptomic and methylomic stability are never
## conflated on one axis.
## ---------------------------------------------------------------------------

mb_stability_plot <- function(sig_df, top_n = 40) {
  need <- c("feature", "omics", "selection_frequency", "stability_category")
  if (is.null(sig_df) || nrow(sig_df) == 0 || !all(need %in% colnames(sig_df))) return(NULL)
  df <- sig_df[!is.na(sig_df$selection_frequency), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df <- df[!duplicated(df[, c("omics", "feature")]), , drop = FALSE]
  df <- df[order(-df$selection_frequency), , drop = FALSE]
  df <- utils::head(df, top_n)
  cat_colors <- c("Stable" = ARTHOMIX_COLORS$aqua, "Moderately stable" = ARTHOMIX_COLORS$blue, "Low stability" = ARTHOMIX_COLORS$ink_muted)
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(feature, selection_frequency), y = selection_frequency, fill = stability_category)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = c(MB_STABILITY_THRESHOLDS$moderate, MB_STABILITY_THRESHOLDS$stable), linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
    ggplot2::scale_fill_manual(values = cat_colors, breaks = names(cat_colors)) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::facet_wrap(~omics, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Selection frequency across CV repetitions", fill = "Stability") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Selected-feature x sample heatmap, annotated by outcome (spec section
## 13.E.7) - real per-sample values for the top selected features (by
## |loading|), z-scored per feature for a comparable color scale. "Annotated
## by outcome" via facet_grid(. ~ outcome) - samples are grouped by outcome
## class rather than layering a separate annotation track, using only
## ggplot2 (no new plotting dependency).
## ---------------------------------------------------------------------------

mb_heatmap_plot <- function(layers, sig_df, outcome, sample_ids, top_n = 40) {
  need <- c("omics", "feature", "loading")
  if (is.null(sig_df) || nrow(sig_df) == 0 || !all(need %in% colnames(sig_df)) || length(sample_ids) < 2) return(NULL)
  top <- sig_df[!duplicated(sig_df$feature), , drop = FALSE]
  top <- top[order(-abs(top$loading)), , drop = FALSE]
  top <- utils::head(top, top_n)

  vals <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    b <- top$omics[i]; f <- top$feature[i]
    if (!b %in% names(layers) || !f %in% colnames(layers[[b]])) return(NULL)
    x <- layers[[b]][sample_ids, f]
    sdv <- stats::sd(x, na.rm = TRUE)
    z <- if (is.na(sdv) || sdv == 0) rep(0, length(x)) else as.numeric(scale(x))
    data.frame(sample = sample_ids, feature = f, omics = b, value = z, stringsAsFactors = FALSE)
  }))
  if (is.null(vals) || nrow(vals) == 0) return(NULL)

  vals$outcome <- as.character(outcome[vals$sample])
  ord_samples <- sample_ids[order(vals$outcome[match(sample_ids, vals$sample)])]
  vals$sample <- factor(vals$sample, levels = unique(ord_samples))
  vals$feature <- factor(vals$feature, levels = rev(unique(top$feature)))

  ggplot2::ggplot(vals, ggplot2::aes(x = sample, y = feature, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = ARTHOMIX_COLORS$blue, mid = "white", high = ARTHOMIX_COLORS$red, midpoint = 0, name = "Z-score") +
    ggplot2::facet_grid(. ~ outcome, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.spacing = grid::unit(0.4, "lines"))
}

## ---------------------------------------------------------------------------
## Between-block component-correlation scatter (spec section 13.D
## "Component correlations", and the DIABLO component-plot recommended in
## section 13.E.2) - real fit$variates component-1 scores for the two
## blocks, colored by outcome. Restricted to exactly two blocks, this
## submodule's own scope (Transcriptomics + Methylomics).
## ---------------------------------------------------------------------------

mb_component_correlation_plot <- function(fit, outcome) {
  cc <- mb_component_correlation(fit)
  if (is.null(cc)) return(NULL)
  ids <- rownames(fit$variates[[cc$block_a]])
  ids <- intersect(ids, rownames(fit$variates[[cc$block_b]]))
  df <- data.frame(
    x = fit$variates[[cc$block_a]][ids, 1], y = fit$variates[[cc$block_b]][ids, 1],
    outcome = as.character(outcome[ids]), stringsAsFactors = FALSE
  )
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = outcome)) +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(df$outcome))) +
    ggplot2::labs(
      x = sprintf("%s component 1", cc$block_a), y = sprintf("%s component 1", cc$block_b), color = "Outcome",
      title = sprintf("r = %.2f (n = %d)", cc$r, cc$n)
    ) +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Pooled cross-validated ROC curve (spec sections 13.B, 13.E.8) - built
## from mb_cv_roc()'s real out-of-fold pooled predictions
## (multiomics_biomarker_helpers.R), never a resubstitution/training-set
## curve. Orientation matches the app's other ROC plots
## (mod_diagnostic.R::diag_roc_plot_pub) - FPR increasing left to right,
## TPR increasing bottom to top - restyled with theme_arthomix() to match
## this module.
## ---------------------------------------------------------------------------

mb_roc_plot <- function(cv_roc) {
  if (is.null(cv_roc)) return(NULL)
  co <- pROC::coords(cv_roc$roc, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  df <- data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity)
  df <- df[order(df$fpr, df$tpr), ]
  label <- sprintf("AUC = %.3f\n%d-fold pooled out-of-fold, n = %d", cv_roc$auc, cv_roc$folds, cv_roc$n_used)
  ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
    ggplot2::geom_line(color = ARTHOMIX_COLORS$blue, linewidth = 1.1) +
    ggplot2::annotate("text", x = 0.97, y = 0.05, hjust = 1, vjust = 0, size = 3.6, label = label) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "1 - Specificity (false positive rate)", y = "Sensitivity (true positive rate)",
      title = sprintf("Cross-validated ROC: %s vs. %s", cv_roc$pos_class, cv_roc$neg_class)
    ) +
    theme_arthomix()
}
