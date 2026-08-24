## R/multiomics/multiomics_plots.R
## Plotting helpers for the Multi-Omics sub-modules - kept separate from the
## module files themselves, same split crossomics_integration_plots.R uses
## for mod_cross_integration.R. Every plot reuses the app's shared
## theme_arthomix()/ARTHOMIX_COLORS (global.R) so it matches the rest of the
## app; every plot here is built once as a ggplot object and reused for both
## on-screen rendering and PNG export (one function, one place each plot is
## drawn), matching cx_quadrant_ggplot()'s convention.

multi_empty_state <- function(msg = "Load a table (Dataset tab) to see results here.") {
  div(class = "empty-note", icon("circle-info"), msg)
}

## Generic "plot, or an explicit empty state" chooser for a `..._ui`
## renderUI() block. A plot function can legitimately return NULL even after
## the caller's own "is data loaded at all" gate passes - e.g. a user filter
## (region, sex, FDR) zeroing out every remaining row. Without this, that
## case fell through to a bare plotOutput() with nothing ever drawn into it
## (a blank box, not the explicit "no data for this selection" text the
## module aims for everywhere else) - evaluates plot_fn() once, up front, to
## decide which of the two to render; renderPlot() below still re-evaluates
## it for the actual drawing.
multi_plot_or_empty <- function(plot_fn, output_id, msg = "No data for the current selection.", height = "420px") {
  p <- tryCatch(plot_fn(), error = function(e) NULL)
  if (is.null(p)) return(multi_empty_state(msg))
  plotOutput(output_id, height = height)
}

## Generic "ggplot -> PNG download" handler, 7x6in @300dpi, matching
## mod_dge.R's download_volcano_png convention - takes a *function* that
## returns the ggplot object (not the plot itself), so it's re-evaluated at
## download time against whatever is reactive at that moment.
multi_png_download <- function(plot_fn, filename_fn) {
  downloadHandler(
    filename = filename_fn,
    content = function(file) {
      p <- plot_fn()
      if (is.null(p)) { grDevices::png(file, width = 7, height = 6, units = "in", res = 300); grDevices::dev.off(); return() }
      ggplot2::ggsave(file, plot = p, width = 7, height = 6, dpi = 300)
    }
  )
}

## ---------------------------------------------------------------------------
## Integration performance: AUROC + 95% CI, one bar per model/omics, colored
## by whether the CI excludes chance (0.5) - the pipeline's own
## `excludes_chance` column, never recomputed here.
## ---------------------------------------------------------------------------

multi_performance_ci_plot <- function(df, label_col = "omics", title = NULL) {
  need <- c(label_col, "auroc", "ci_lo", "ci_hi")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  df$.label <- df[[label_col]]
  ## Per-row, not per-table: ifelse()'s output length follows its `test`
  ## argument, so wrapping the whole per-row ifelse() in a scalar
  ## isTRUE(all(...)) check (a past bug here) collapsed every row to row 1's
  ## classification. Guard the column's existence once, outside ifelse(),
  ## then classify every row independently.
  df$.excludes <- if ("excludes_chance" %in% colnames(df)) {
    ifelse(df$excludes_chance %in% TRUE, "Excludes chance (CI > 0.5)", "Includes chance")
  } else {
    rep("Includes chance", nrow(df))
  }
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(.label, auroc), y = auroc, color = .excludes)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), linewidth = 0.6, size = 0.6) +
    ggplot2::scale_color_manual(values = c("Excludes chance (CI > 0.5)" = ARTHOMIX_COLORS$aqua, "Includes chance" = ARTHOMIX_COLORS$ink_muted)) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "AUROC (95% CI)", color = NULL, title = title) +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## DIABLO per-patient component score, colored by response - a direct plot
## of the pipeline's own Table29b/34b/40b `score` column (LOOCV-derived, per
## the pipeline's method), not a re-fit.
## ---------------------------------------------------------------------------

multi_diablo_score_plot <- function(df) {
  need <- c("patient_id", "response", "score")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(patient_id, score), y = score, fill = response)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(df$response))) +
    ggplot2::labs(x = "Patient", y = "DIABLO component score", fill = "Response") +
    theme_arthomix() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

## ---------------------------------------------------------------------------
## DIABLO panel loadings, one bar per selected feature, faceted/colored by
## omics view - a direct plot of Table30/35/41's `loading` column.
## ---------------------------------------------------------------------------

multi_diablo_panel_plot <- function(df, top_n = 20) {
  need <- c("feature", "loading", "view")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  df <- df[order(-abs(df$loading)), , drop = FALSE]
  df <- utils::head(df, top_n)
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(feature, loading), y = loading, fill = view)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(df$view))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "DIABLO loading (component 1)", fill = "Omics view") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Gene<->CpG concordance scatter: expression log2FC vs. methylation delta-M,
## colored by the pipeline's own biological_pattern classification - direct
## plot of Table42/Table45, no reclassification performed here.
## ---------------------------------------------------------------------------

multi_concordance_scatter <- function(df, region_filter = NULL) {
  need <- c("expr_logFC", "delta_M", "biological_pattern")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  if (!is.null(region_filter) && length(region_filter) > 0 && "region" %in% colnames(df)) {
    df <- df[df$region %in% region_filter, , drop = FALSE]
  }
  if (nrow(df) == 0) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = expr_logFC, y = delta_M, color = biological_pattern)) +
    ggplot2::geom_hline(yintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_vline(xintercept = 0, color = ARTHOMIX_COLORS$grid) +
    ggplot2::geom_point(alpha = 0.6, size = 1.6) +
    ggplot2::scale_color_manual(values = arthomix_pair(unique(df$biological_pattern))) +
    ggplot2::labs(x = "Expression log2FC (responder vs. non-responder)", y = "Methylation delta-M",
                  color = "Pattern") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## Pathway enrichment dotplot - direct plot of Table43's own
## GeneRatio/p.adjust/Count columns (clusterProfiler/ReactomePA output), no
## enrichment test re-run here.
## ---------------------------------------------------------------------------

multi_pathway_dotplot <- function(df, top_n = 20) {
  need <- c("Description", "p.adjust", "Count", "GeneRatio")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  ratio_to_num <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) if (length(p) == 2 && as.numeric(p[2]) > 0) as.numeric(p[1]) / as.numeric(p[2]) else NA_real_, numeric(1))
  }
  df$.ratio <- ratio_to_num(df$GeneRatio)
  df <- df[order(df$p.adjust), , drop = FALSE]
  df <- utils::head(df, top_n)
  ## Real Table43 Description strings run up to 132 characters - wrapped so
  ## the y-axis labels of a horizontal dotplot don't get clipped or crowd
  ## the plot area off-screen.
  df$.label <- vapply(df$Description, function(s) paste(strwrap(s, width = 45), collapse = "\n"), character(1))
  ggplot2::ggplot(df, ggplot2::aes(x = .ratio, y = stats::reorder(.label, .ratio), size = Count, color = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = ARTHOMIX_COLORS$red, high = ARTHOMIX_COLORS$blue) +
    ggplot2::labs(x = "Gene ratio", y = NULL, size = "Gene count", color = "Adj. P") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## SNF cluster composition: response distribution within each fused cluster
## - direct tabulation of Table_SNFjoint_cluster_assignments_*'s own
## snf_cluster/response columns.
## ---------------------------------------------------------------------------

multi_cluster_composition_plot <- function(df) {
  need <- c("snf_cluster", "response")
  if (is.null(df) || nrow(df) == 0 || !all(need %in% colnames(df))) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = factor(snf_cluster), fill = response)) +
    ggplot2::geom_bar(position = "stack") +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(df$response))) +
    ggplot2::labs(x = "Fused SNF cluster", y = "Patients", fill = "Response") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## DIABLO variance explained per component, per omics block - direct plot of
## the saved block.splsda fit's own `prop_expl_var` (real mixOmics output,
## read via multi_diablo_fit() in multiomics_helpers.R), not re-fit here.
## This is the correctly-scoped, correctly-named ("component", not "factor")
## analog of a MOFA-style variance-explained plot for a method that isn't MOFA.
## ---------------------------------------------------------------------------

multi_diablo_variance_plot <- function(prop_df) {
  need <- c("block", "component", "variance_explained")
  if (is.null(prop_df) || nrow(prop_df) == 0 || !all(need %in% colnames(prop_df))) return(NULL)
  ggplot2::ggplot(prop_df, ggplot2::aes(x = component, y = variance_explained, fill = block)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = arthomix_pair(unique(prop_df$block))) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "DIABLO component", y = "Variance explained (within block)", fill = "Omics block") +
    theme_arthomix()
}

## ---------------------------------------------------------------------------
## SNF per-omics contribution to the fused network - direct plot of the
## pipeline's own Table_SNFjoint_concordanceNMI_* (view_row x view_col NMI
## matrix), the real, already-computed analog of "how much each omics layer
## agrees with the final clustering." Not re-derived from the fused network.
## ---------------------------------------------------------------------------

multi_nmi_heatmap <- function(nmi_df) {
  need <- c("view_row", "view_col", "nmi")
  if (is.null(nmi_df) || nrow(nmi_df) == 0 || !all(need %in% colnames(nmi_df))) return(NULL)
  ggplot2::ggplot(nmi_df, ggplot2::aes(x = view_col, y = view_row, fill = nmi)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", nmi)), color = ARTHOMIX_COLORS$ink, size = 3.5) +
    ggplot2::scale_fill_gradient(low = "white", high = ARTHOMIX_COLORS$blue, limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "NMI") +
    theme_arthomix() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

## ---------------------------------------------------------------------------
## Sample-overlap Venn (RNA-available vs. methylation-available) - reuses
## ggVennDiagram, the same charting library Cross-Omics's own biomarker-
## convergence overlap plot already uses (mod_cross_biomarker_conv.R), not a
## new dependency.
## ---------------------------------------------------------------------------

multi_sample_overlap_venn <- function(matching_df) {
  if (is.null(matching_df) || !all(c("RNA_available_PBMC", "methylation_available", "patient_id") %in% colnames(matching_df))) return(NULL)
  rna_ok <- matching_df$RNA_available_PBMC %in% c(TRUE, "TRUE", "Yes", "yes", 1)
  meth_ok <- matching_df$methylation_available %in% c(TRUE, "TRUE", "Yes", "yes", 1)
  sets <- list(`RNA-seq` = matching_df$patient_id[rna_ok], Methylation = matching_df$patient_id[meth_ok])
  if (all(lengths(sets) == 0)) return(NULL)
  ggVennDiagram::ggVennDiagram(sets, label = "count") +
    ggplot2::scale_fill_gradient(low = "white", high = ARTHOMIX_COLORS$blue) +
    theme_arthomix() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}
