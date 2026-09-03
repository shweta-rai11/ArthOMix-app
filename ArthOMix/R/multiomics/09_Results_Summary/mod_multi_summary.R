## R/multiomics/09_Results_Summary/mod_multi_summary.R
## Submodule: Results Summary & Reproducibility - rolls up whatever the
## other Multi-Omics sub-modules have published to `multi_results` this

mod_multi_summary_config <- list(
  id = "summary", title = "Results Summary & Reproducibility", icon = "clipboard-list", group = "Interpretation",
  description = "Session summary and downloadable results bundle."
)

mod_multi_summary_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      12,
      box(
        width = NULL, title = "Session dashboard", status = "primary", solidHeader = FALSE,
        uiOutput(ns("dashboard_ui"))
      )
    ),
    column(
      6,
      box(
        width = NULL, title = "Supported vs. not implemented", status = "primary", solidHeader = FALSE,
        tags$h5("Implemented (live analysis, computed on click, on the data source you pick within each tab)"),
        tags$ul(
          tags$li("DIABLO integration (supervised, mixOmics::block.splsda) - performance, per-patient scores, panel loadings; run live on either the Active Multi-Omics Dataset or a preloaded RA anti-TNF cell (rehydrated from its saved fit)"),
          tags$li("DIABLO variance explained per component, per omics block - from that same live run"),
          tags$li("SNF integration and patient clustering (unsupervised, SNFtool::SNF) - likewise run live on either data source"),
          tags$li("SNF per-omics NMI contribution to the fused network - from that same live run"),
          tags$li("Gene <-> CpG concordance with genomic-region awareness - a precomputed browse of the pipeline's own Table42/45 tables when \"Preloaded\" is selected, or a live data-adaptive analysis over the Active Multi-Omics Dataset when \"Active\" is selected; optional BH-FDR recompute over retained raw p-values either way"),
          tags$li("Pathway enrichment (GO/KEGG/Reactome) - run live on click, over the preloaded candidate panel or your own upload"),
          tags$li("Joint biomarker candidate ranking (live DIABLO feature selection) - feature-selection stability is a fixed-threshold evidence label from real cross-validation, never a user-adjustable confidence score")
        ),
        tags$h5("Implemented (Live Analysis tab - your own uploaded data)"),
        tags$ul(
          tags$li("Upload, validation, sample matching, missing-data QC for up to 4 omics layers"),
          tags$li("Omics-appropriate normalization, feature filtering, cross-omics scaling"),
          tags$li("Batch diagnostics (PCA before/after) and correction (ComBat / limma::removeBatchEffect)"),
          tags$li("MOFA2 factor analysis (asynchronous): factor scores/loadings/variance explained/per-view contribution"),
          tags$li("Cross-omics feature correlation (Pearson/Spearman) with FDR-filtered heatmap")
        ),
        tags$h5("Not implemented in this delivery"),
        tags$ul(lapply(MULTI_KNOWN_LIMITATIONS, tags$li))
      )
    ),
    column(
      6,
      box(
        width = NULL, title = "Reproducibility - source scripts", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Precomputed-cohort figures were generated once by these scripts (Research_05_multiomics_sexstratified). Live Analysis results are computed in this app on your own data."),
        tags$ul(lapply(MULTI_REPRODUCIBILITY_SCRIPTS, tags$li)),
        div(class = "table-toolbar", downloadButton(ns("dl_bundle"), "Download session bundle (CSV + report, ZIP)", class = "btn-sm btn-primary"))
      ),
      box(
        width = NULL, title = "Software versions", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Package versions installed in this deployment."),
        DT::dataTableOutput(ns("versions_table"))
      )
    )
  )
}

mod_multi_summary_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$dashboard_ui <- renderUI({
      res <- multi_results %||% reactiveValues()
      overview <- res$overview
      integration <- res$integration
      strat <- res$stratification
      biomarker <- res$biomarker
      concordance <- res$concordance
      pathway <- res$pathway

      cards <- list(
        list(label = "Matched samples", value = if (!is.null(overview)) overview$harmonization$n_matched else NA, color = "blue"),
        list(label = "Loaded integration cell", value = if (!is.null(integration)) integration$cell$label else "None loaded", color = "aqua"),
        list(label = "Candidate biomarker features", value = if (!is.null(biomarker)) length(unique(biomarker$df$feature)) else NA, color = "violet"),
        list(label = "Gene-CpG pairs (concordance)", value = if (!is.null(concordance)) nrow(concordance$df) else NA, color = "orange"),
        list(label = "Enriched pathway terms", value = if (!is.null(pathway)) nrow(pathway$df) else NA, color = "red"),
        list(label = "SNF clusters loaded", value = if (!is.null(strat)) length(unique(strat$assignments$snf_cluster)) else NA, color = "magenta")
      )
      div(style = "display:flex; gap:10px; flex-wrap:wrap;",
          lapply(cards, function(c) {
            val <- c$value
            shown <- if (is.null(val) || (length(val) == 1 && is.na(val))) "Not loaded" else if (is.numeric(val)) format(val, big.mark = ",") else as.character(val)
            div(class = "card", style = "flex:1 1 170px; text-align:center; padding:10px;",
                div(style = sprintf("font-size:1.2em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[c$color]]), shown),
                div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", c$label))
          })
      )
    })

    output$versions_table <- DT::renderDataTable({
      DT::datatable(multi_package_versions(), rownames = FALSE, options = list(dom = "t", pageLength = 10), class = "stripe hover compact")
    })

    output$dl_bundle <- downloadHandler(
      filename = function() paste0("multiomics_session_bundle_", Sys.Date(), ".zip"),
      content = function(file) {
        tmp <- tempfile(); dir.create(tmp)
        res <- if (!is.null(multi_results)) reactiveValuesToList(multi_results) else list()
        write_if_df <- function(x, name) {
          df <- if (is.list(x) && !is.data.frame(x)) x$df %||% x$assignments %||% x$result$panels_wide %||% NULL else x
          if (is.data.frame(df)) utils::write.csv(df, file.path(tmp, paste0(name, ".csv")), row.names = FALSE)
        }
        for (nm in names(res)) write_if_df(res[[nm]], nm)
        writeLines(multi_build_report(res), file.path(tmp, "report.md"))
        old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
        setwd(tmp)
        utils::zip(file, files = list.files(tmp))
      }
    )
  })
}
