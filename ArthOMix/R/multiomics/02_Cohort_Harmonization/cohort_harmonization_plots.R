## R/multiomics/02_Cohort_Harmonization/cohort_harmonization_plots.R
## Plotting helpers for the "Cohort Harmonization" sub-module
## (mod_multi_overview.R). PCA and cross-modality correlation reuse

ch_overlap_heatmap_plot <- function(overlap_matrix) {
  if (is.null(overlap_matrix) || nrow(overlap_matrix) == 0) return(NULL)
  df <- as.data.frame(as.table(overlap_matrix))
  colnames(df) <- c("ModalityA", "ModalityB", "n")
  ggplot2::ggplot(df, ggplot2::aes(x = ModalityB, y = ModalityA, fill = n)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = n), size = 3.6) +
    ggplot2::scale_fill_gradient(low = "white", high = ARTHOMIX_COLORS$blue) +
    theme_arthomix() +
    ggplot2::labs(x = NULL, y = NULL, fill = "Shared\nsamples") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

ch_completeness_heatmap_plot <- function(id_sets, max_samples = 150) {
  id_sets <- Filter(function(x) length(x) > 0, id_sets)
  if (length(id_sets) == 0) return(NULL)
  all_ids <- unique(unlist(id_sets))
  truncated <- length(all_ids) > max_samples
  if (truncated) all_ids <- sort(all_ids)[seq_len(max_samples)]
  df <- do.call(rbind, lapply(names(id_sets), function(nm) {
    data.frame(Sample = all_ids, Modality = nm, Present = all_ids %in% id_sets[[nm]])
  }))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Modality, fill = Present)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = ARTHOMIX_COLORS$blue, `FALSE` = "#EEEEEE")) +
    theme_arthomix() +
    ggplot2::labs(x = if (truncated) sprintf("Sample (first %d of %d shown)", max_samples, length(unique(unlist(id_sets)))) else "Sample", y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  p
}

ch_sample_highlight_pca_plot <- function(pca_obj, highlight_id) {
  if (is.null(pca_obj) || !isTRUE(pca_obj$ok)) return(NULL)
  df <- pca_obj$scores
  df$Sample <- pca_obj$sample_ids
  df$Highlight <- ifelse(df$Sample == highlight_id, "Selected sample", "Other samples")
  ve <- pca_obj$var_explained
  ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = Highlight, size = Highlight)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_manual(values = c(`Selected sample` = ARTHOMIX_COLORS$red, `Other samples` = "grey70")) +
    ggplot2::scale_size_manual(values = c(`Selected sample` = 4, `Other samples` = 2), guide = "none") +
    theme_arthomix() +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2]), color = NULL)
}

ch_category_bar_plot <- function(meta, col, title = NULL) {
  if (is.null(meta) || !col %in% colnames(meta)) return(NULL)
  v <- meta[[col]]
  df <- as.data.frame(table(Level = as.character(v)), responseName = "n")
  if (nrow(df) == 0) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Level, -n), y = n)) +
    ggplot2::geom_col(fill = ARTHOMIX_COLORS$blue) +
    theme_arthomix() +
    ggplot2::labs(x = NULL, y = "Samples", title = title %||% col) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
