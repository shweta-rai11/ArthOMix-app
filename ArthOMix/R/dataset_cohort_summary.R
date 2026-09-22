## R/dataset_cohort_summary.R
## Cohort composition card for a module's Dataset tab: a Cohort x Female/Male/Total count table,
## shown as soon as sample metadata is loaded, plus a note on which sex designs the data supports.
## This is read-only context, not a control - each analysis page already has its own sex-design
## picker (DGE's covariate filter/adjust, DMP's sex radio), and "Start an analysis" (R/workflow_guide.R)
## already walks the user through choosing one with real effect.
##
## Deliberately self-contained (does not call workflow_guide.R's wf_group_col()/wf_sex_col(), which
## would pull that file - and, through it, the methylomics DMP module - into every Dataset tab test).
## The lookup below matches theirs: the same literal column names the Dataset tabs already
## standardise "group"/"sex" to on load.
cohort_col <- function(meta, names) {
  if (!is.data.frame(meta)) return(NULL)
  hit <- intersect(names, colnames(meta))
  if (length(hit)) hit[1] else NULL
}

cohort_summary_ui <- function(id) {
  uiOutput(NS(id, "box"))
}

## meta_reactive: a reactive() (or function) returning the sample metadata data.frame, or NULL/empty
## when nothing is loaded yet. Renders nothing when there is no metadata or no usable group column.
cohort_summary_server <- function(id, meta_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    info <- reactive({
      meta <- meta_reactive()
      if (!is.data.frame(meta) || nrow(meta) == 0) return(NULL)
      gcol <- cohort_col(meta, c("group", "Group", "disease", "Disease"))
      if (is.null(gcol)) return(NULL)

      g <- as.character(meta[[gcol]])
      scol <- cohort_col(meta, c("sex", "Sex", "SEX", "gender", "Gender"))
      has_sex <- !is.null(scol) && scol %in% colnames(meta)

      if (has_sex) {
        raw <- as.character(meta[[scol]])
        s <- ifelse(grepl("^f(emale)?$", raw, ignore.case = TRUE), "Female",
             ifelse(grepl("^m(ale)?$", raw, ignore.case = TRUE), "Male", raw))
        tbl <- as.data.frame.matrix(table(g, s))
        ordered <- intersect(c("Female", "Male"), colnames(tbl))
        tbl <- tbl[, c(ordered, setdiff(colnames(tbl), ordered)), drop = FALSE]
        tbl$Total <- rowSums(tbl)
      } else {
        counts <- table(g)
        tbl <- data.frame(Total = as.integer(counts))
        rownames(tbl) <- names(counts)
      }
      tbl <- cbind(Cohort = rownames(tbl), tbl)
      rownames(tbl) <- NULL

      sex_ok <- has_sex && length(unique(stats::na.omit(as.character(meta[[scol]])))) >= 2
      list(table = tbl, has_sex = has_sex, sex_ok = sex_ok)
    })

    output$table <- DT::renderDataTable({
      inf <- info(); req(inf)
      DT::datatable(inf$table, rownames = FALSE, options = list(dom = "t", scrollX = TRUE),
                    class = "stripe hover compact")
    })

    output$box <- renderUI({
      inf <- info()
      if (is.null(inf)) return(NULL)
      div(
        class = "card",
        div(class = "card-title", icon("venus-mars"), "Cohort composition"),
        DT::dataTableOutput(ns("table")),
        p(class = "submodule-desc",
          if (inf$sex_ok) {
            "This dataset supports pooled (sex-adjusted), sex-stratified and sex-specific analyses. Pick the sex design on the analysis page you open, or use \"Start an analysis\" to be guided through it."
          } else if (inf$has_sex) {
            "The sex column has fewer than two usable values, so only pooled analysis is available for this dataset."
          } else {
            "No usable sex column was found, so only pooled analysis is available for this dataset."
          })
      )
    })
  })
}
