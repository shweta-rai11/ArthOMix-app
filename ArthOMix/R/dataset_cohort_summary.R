## R/dataset_cohort_summary.R
## Cohort composition card for a module's Dataset tab: a Cohort x Female/Male/Total count table,
## an "Analysis strategy" picker, and per-strategy small-strata warnings, all shown as soon as
## sample metadata is loaded (reviewer request, 2026-09-22).
##
## The picker is deliberately a decision-support view, not a run control: each analysis page (DGE's
## covariate filter/adjust, DMP's sex radio) and "Start an analysis" (R/workflow_guide.R) already
## have their own sex-design picker with real effect on a run. Duplicating that here as a second,
## disconnected control that also claims to "run" something would drift out of sync with them. What
## this card adds instead is what neither of those has: the group-by-sex counts and a warning the
## moment a stratum is too small for the strategy you're looking at, before you go pick it there.
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

## Below this many samples in a stratum, flag it as small (matches the wording style already used
## elsewhere in the app for "N or more samples" checks; 10 catches the reviewer's own "only 5
## samples" example while not flagging every reasonably-sized cohort).
COHORT_STRATUM_WARN_N <- 10L
## Below this many, the stratum cannot support the strategy at all (matches wf_check_contrast()'s
## min_per_group elsewhere in the app).
COHORT_STRATUM_MIN_N <- 3L

COHORT_STRATEGIES <- list(
  all     = list(label = "All subjects, sex-adjusted",   needs = "total"),
  female  = list(label = "Females only",                 needs = "female"),
  male    = list(label = "Males only",                   needs = "male"),
  compare = list(label = "Female vs male RA effects",     needs = "both"),
  full    = list(label = "Full disease × sex interaction", needs = "cells")
)

## Eligibility + warning text for one strategy, given the group x {Female, Male, Total} table.
cohort_strategy_note <- function(strategy, tbl) {
  spec <- COHORT_STRATEGIES[[strategy]]
  has_sex <- all(c("Female", "Male") %in% colnames(tbl))
  small <- function(n) n > 0 & n < COHORT_STRATUM_WARN_N
  blocked <- function(n) n < COHORT_STRATUM_MIN_N

  cells <- function(col) stats::setNames(as.integer(tbl[[col]]), tbl$Cohort)

  if (identical(spec$needs, "total")) {
    n <- cells("Total")
    if (any(blocked(n))) return(list(ok = FALSE, warn = FALSE,
      text = sprintf("Blocked: %s has fewer than %d samples.", names(n)[blocked(n)][1], COHORT_STRATUM_MIN_N)))
    if (any(small(n))) return(list(ok = TRUE, warn = TRUE,
      text = sprintf("%s group contains only %d samples. Estimates may be unstable.",
                     names(n)[small(n)][1], n[small(n)][1])))
    return(list(ok = TRUE, warn = FALSE, text = "Every group has enough samples."))
  }
  if (!has_sex) return(list(ok = FALSE, warn = FALSE, text = "No usable sex column - this strategy needs one."))

  if (spec$needs %in% c("female", "male")) {
    col <- if (identical(spec$needs, "female")) "Female" else "Male"
    n <- cells(col)
    if (any(blocked(n))) return(list(ok = FALSE, warn = FALSE,
      text = sprintf("Blocked: %s %s group contains only %d sample%s.",
                     col, names(n)[blocked(n)][1], n[blocked(n)][1], if (n[blocked(n)][1] == 1) "" else "s")))
    if (any(small(n))) return(list(ok = TRUE, warn = TRUE,
      text = sprintf("%s %s group contains only %d samples. Sex-stratified estimates may be unstable.",
                     col, names(n)[small(n)][1], n[small(n)][1])))
    return(list(ok = TRUE, warn = FALSE, text = sprintf("Every group has enough %s samples.", tolower(col))))
  }
  if (identical(spec$needs, "both")) {
    nf <- cells("Female"); nm <- cells("Male")
    if (any(blocked(nf)) || any(blocked(nm))) {
      bad_col <- if (any(blocked(nf))) "Female" else "Male"
      bad_n <- if (any(blocked(nf))) nf else nm
      return(list(ok = FALSE, warn = FALSE,
        text = sprintf("Blocked: %s %s group contains only %d sample%s.",
                       bad_col, names(bad_n)[blocked(bad_n)][1], bad_n[blocked(bad_n)][1],
                       if (bad_n[blocked(bad_n)][1] == 1) "" else "s")))
    }
    if (any(small(nf)) || any(small(nm))) {
      warn_col <- if (any(small(nf))) "Female" else "Male"
      warn_n <- if (any(small(nf))) nf else nm
      return(list(ok = TRUE, warn = TRUE,
        text = sprintf("%s %s group contains only %d samples. Sex-stratified DEG estimates may be unstable.",
                       warn_col, names(warn_n)[small(warn_n)][1], warn_n[small(warn_n)][1])))
    }
    return(list(ok = TRUE, warn = FALSE, text = "Every sex-by-group cell has enough samples to compare."))
  }
  if (identical(spec$needs, "cells")) {
    nf <- cells("Female"); nm <- cells("Male")
    if (any(nf == 0) || any(nm == 0)) return(list(ok = FALSE, warn = FALSE,
      text = "Blocked: every sex-by-group cell needs at least one sample for an interaction model."))
    if (any(small(nf)) || any(small(nm))) {
      warn_col <- if (any(small(nf))) "Female" else "Male"
      warn_n <- if (any(small(nf))) nf else nm
      return(list(ok = TRUE, warn = TRUE,
        text = sprintf("%s %s group contains only %d samples. The interaction estimate may be unstable.",
                       warn_col, names(warn_n)[small(warn_n)][1], warn_n[small(warn_n)][1])))
    }
    return(list(ok = TRUE, warn = FALSE, text = "Every sex-by-group cell has enough samples for an interaction model."))
  }
  list(ok = TRUE, warn = FALSE, text = "")
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
      ## An uploaded dataset's sex column exists but is entirely NA whenever the user leaves the
      ## "Sex column" mapping on "(none)" - the column being present is not enough to prove there
      ## is anything usable in it (table(g, s) on all-NA silently produces zero columns).
      has_sex <- !is.null(scol) && scol %in% colnames(meta) &&
        any(!is.na(meta[[scol]]) & nzchar(as.character(meta[[scol]])))

      ## Total always comes from every row, sex-mapped or not - an upload with only some rows'
      ## sex mapped (or the sex column only partly filled in) must not silently drop the rest of
      ## a group out of its own count just because table(g, s) drops NAs from the cross-tab.
      group_totals <- table(g)
      if (has_sex) {
        raw <- as.character(meta[[scol]])
        s <- ifelse(grepl("^f(emale)?$", raw, ignore.case = TRUE), "Female",
             ifelse(grepl("^m(ale)?$", raw, ignore.case = TRUE), "Male", raw))
        tbl <- as.data.frame.matrix(table(g, s))
        ordered <- intersect(c("Female", "Male"), colnames(tbl))
        tbl <- tbl[, c(ordered, setdiff(colnames(tbl), ordered)), drop = FALSE]
        n_unmapped <- as.integer(group_totals[rownames(tbl)]) - rowSums(tbl)
        if (any(n_unmapped > 0)) tbl$"Sex unknown" <- pmax(n_unmapped, 0L)
        tbl$Total <- as.integer(group_totals[rownames(tbl)])
      } else {
        tbl <- data.frame(Total = as.integer(group_totals))
        rownames(tbl) <- names(group_totals)
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

    ## Every stratum too small for ANY strategy, shown unconditionally - a property of the data,
    ## not of what happens to be selected below.
    output$warnings <- renderUI({
      inf <- info(); req(inf, inf$has_sex)
      notes <- lapply(names(COHORT_STRATEGIES), function(k) cohort_strategy_note(k, inf$table))
      flagged <- Filter(function(n) isTRUE(n$warn) || !isTRUE(n$ok), notes)
      texts <- unique(vapply(flagged, `[[`, character(1), "text"))
      if (!length(texts)) return(NULL)
      tagList(lapply(texts, function(t) div(class = "empty-note wf-note-error",
        icon("triangle-exclamation"), " ", t)))
    })

    output$strategy_ui <- renderUI({
      inf <- info(); req(inf)
      opts <- if (inf$has_sex) COHORT_STRATEGIES else COHORT_STRATEGIES["all"]
      radioButtons(ns("strategy"), "Analysis strategy",
                   choices = stats::setNames(names(opts), vapply(opts, `[[`, character(1), "label")),
                   selected = "all")
    })

    output$strategy_note <- renderUI({
      inf <- info(); req(inf)
      strat <- input$strategy %||% "all"
      if (!strat %in% names(COHORT_STRATEGIES)) return(NULL)
      note <- cohort_strategy_note(strat, inf$table)
      div(class = paste("empty-note", if (!note$ok) "wf-note-error" else if (note$warn) "wf-note-error" else "wf-note-ok"),
          icon(if (!note$ok || note$warn) "triangle-exclamation" else "circle-check"), " ", note$text)
    })

    output$box <- renderUI({
      inf <- info()
      if (is.null(inf)) return(NULL)
      div(
        class = "card",
        div(class = "card-title", icon("venus-mars"), "Cohort composition"),
        DT::dataTableOutput(ns("table")),
        uiOutput(ns("warnings")),
        p(class = "submodule-desc",
          "Pick the strategy you're considering to see whether this dataset supports it. This is context, not a run control - set the sex design for a specific analysis on that analysis's own page, or use \"Start an analysis\" to be guided through it."),
        uiOutput(ns("strategy_ui")),
        uiOutput(ns("strategy_note"))
      )
    })
  })
}
