## R/multiomics/functions/multiomics_plots.R
## Plotting helpers for the Multi-Omics sub-modules - kept separate from the
## module files themselves, same split crossomics_integration_plots.R uses

multi_empty_state <- function(msg = "Load a table (Dataset tab) to see results here.") {
  div(class = "empty-note", icon("circle-info"), msg)
}

multi_plot_or_empty <- function(plot_fn, output_id, msg = "No data for the current selection.", height = "420px") {
  p <- tryCatch(plot_fn(), error = function(e) {
    message(sprintf("multi_plot_or_empty(%s): %s", output_id, conditionMessage(e)))
    NULL
  })
  if (is.null(p)) return(multi_empty_state(msg))
  plotly::plotlyOutput(output_id, height = height)
}

multi_render_plotly <- function(plot_fn) {
  plotly::renderPlotly({
    p <- tryCatch(plot_fn(), error = function(e) {
      message(sprintf("multi_render_plotly: %s", conditionMessage(e)))
      NULL
    })
    req(p)
    gp <- plotly::ggplotly(p, tooltip = "all")
    gp <- plotly::layout(gp, hoverlabel = list(bgcolor = "white", font = list(size = 12)))
    plotly::config(gp, displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
  })
}

multi_png_download <- function(plot_fn, filename_fn) {
  downloadHandler(
    filename = filename_fn,
    content = function(file) {
      p <- tryCatch(plot_fn(), error = function(e) {
        message(sprintf("multi_png_download(%s): %s", filename_fn(), conditionMessage(e)))
        NULL
      })
      if (is.null(p)) { grDevices::png(file, width = 7, height = 6, units = "in", res = 300); grDevices::dev.off(); return() }
      ggplot2::ggsave(file, plot = p, width = 7, height = 6, dpi = 300)
    }
  )
}

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

