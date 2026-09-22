## R/workflow_guide.R
## "Start an analysis" guided workflow (external reviewer request, 2026-09-22).
## A thin controller over existing sub-modules: it opens the existing tabs, sets their existing inputs,
## presses their existing run buttons and reads their existing result objects. It computes no statistics.

WF_STEPS <- list(
  list(n = 1L, label = "Choose data",         icon = "database",     blurb = "Pick transcriptomics or methylomics and load a dataset."),
  list(n = 2L, label = "Check samples",       icon = "users",        blurb = "Confirm the sample, group and sex counts."),
  list(n = 3L, label = "Define comparison",   icon = "code-compare", blurb = "Pick the column and the two groups to compare."),
  list(n = 4L, label = "Choose sex analysis", icon = "venus-mars",   blurb = "Pooled, stratified by sex, or sex-specific."),
  list(n = 5L, label = "Run analysis",        icon = "play",         blurb = "Run the existing differential analysis."),
  list(n = 6L, label = "Explore biomarkers",  icon = "star",         blurb = "Open the candidate features and the Biomarker Card."),
  list(n = 7L, label = "Validate",            icon = "flask-vial",   blurb = "Check the signal in an independent cohort.")
)
WF_N_STEPS <- length(WF_STEPS)

## Where each step lands. Ids are existing sub-module ids (TX_MODULES / MX_MODULES).
WF_LAYERS <- list(
  transcriptomics = list(
    label = "Transcriptomics", example_id = "__default_merged__",
    qc = "overview", compare = "dge", interaction = "interaction",
    candidates = "candidates", card = "biomarkercard",
    validate = list(
      list(id = "diagnostic", tabset = "main_tabs", tab = "External Validation"),
      list(id = "crosstissue")
    ),
    feature = "genes"
  ),
  methylomics = list(
    label = "Methylomics", example_id = "gse42861_wholeblood",
    qc = "qc", compare = "dmp", compare_tabset = "dmp_subtabs", compare_tab = "DMP", interaction = "interaction",
    candidates = "candidates", card = "biomarkercard",
    validate = list(list(id = "validation")),
    feature = "CpGs"
  )
)

WF_SEX_MODES <- c(pooled = "Pooled", stratified = "Stratified", sex_specific = "Sex-specific")

WF_CONTROL_PATTERN <- "^(hc|control|controls|ctrl|healthy|healthy control|healthy controls|normal|nc|con)$"

## Group column: the Dataset tabs standardise it to "group" for transcriptomics; methylomics uses the same
## lookup as the DMP tab's own default.
wf_group_col <- function(meta) {
  if (!is.data.frame(meta)) return(NULL)
  hit <- intersect(c("group", "Group", "disease", "Disease"), colnames(meta))
  if (length(hit)) hit[1] else NULL
}

## Sex column: same lookup the DMP and Sex Interaction tabs use.
wf_sex_col <- function(meta) {
  if (!is.data.frame(meta)) return(NULL)
  mod_methyl_dmp_sex_col(meta)
}

## Named vector of usable sex levels (e.g. c(Female = "F", Male = "M")), via the DMP tab's own helper.
wf_sex_levels <- function(meta, sex_col) {
  if (is.null(sex_col) || !is.data.frame(meta) || !sex_col %in% colnames(meta)) return(character(0))
  ch <- mod_methyl_dmp_sex_choices(meta, sex_col)
  ch[ch != "__all__"]
}

## Step 2: does the loaded dataset carry what the rest of the workflow needs?
wf_check_samples <- function(meta, group_col, sex_col, min_per_group = 3L) {
  out <- list(ok = FALSE, n = if (is.data.frame(meta)) nrow(meta) else 0L,
              group_col = group_col, sex_col = sex_col,
              group_counts = NULL, sex_counts = NULL, n_sex_missing = NA_integer_,
              sex_ok = FALSE, problems = character(0), notes = character(0))
  if (!is.data.frame(meta) || nrow(meta) == 0) {
    out$problems <- "Please select a valid dataset before continuing."
    return(out)
  }
  if (is.null(group_col) || !group_col %in% colnames(meta)) {
    out$problems <- "The sample metadata has no group column (looked for group, Group, disease or Disease). Load the data again on the Dataset tab and map the group or diagnosis column."
  } else {
    g <- table(as.character(stats::na.omit(meta[[group_col]])))
    out$group_counts <- g
    if (sum(g >= min_per_group) < 2) {
      out$problems <- c(out$problems, sprintf("At least two comparison groups are required, each with %d or more samples.", min_per_group))
    }
  }
  if (!is.null(sex_col) && sex_col %in% colnames(meta)) {
    s <- table(as.character(stats::na.omit(meta[[sex_col]])))
    out$sex_counts <- s
    out$n_sex_missing <- sum(is.na(meta[[sex_col]]))
    out$sex_ok <- length(s) >= 2
    if (out$sex_ok && out$n_sex_missing > 0) {
      out$notes <- c(out$notes, sprintf("%d sample(s) have no sex value and are left out of sex-based analyses.", out$n_sex_missing))
    }
  }
  if (!out$sex_ok) {
    out$notes <- c(out$notes, "No usable sex column was found (it needs at least two values). Only the pooled analysis can be configured until sex metadata is available.")
  }
  out$ok <- length(out$problems) == 0
  out
}

## Step 3 default: a control-like level as the reference when one exists, otherwise the page's own
## default (first two levels in sorted order).
wf_default_contrast <- function(values) {
  lv <- sort(unique(as.character(stats::na.omit(values))))
  if (length(lv) < 2) return(NULL)
  ctrl <- lv[grepl(WF_CONTROL_PATTERN, lv, ignore.case = TRUE)]
  ref <- if (length(ctrl)) ctrl[1] else lv[1]
  list(ref = ref, comp = setdiff(lv, ref)[1])
}

## Step 3/4: is the chosen comparison runnable, and which sex designs does it support?
wf_check_contrast <- function(meta, col, ref, comp, sex_col, layer, min_per_group = 3L) {
  out <- list(ok = FALSE, problems = character(0), n_ref = 0L, n_comp = 0L,
              strata = NULL, sex_levels = character(0),
              modes = c(pooled = TRUE, stratified = FALSE, sex_specific = FALSE),
              mode_notes = character(0))
  if (!is.data.frame(meta) || is.null(col) || !col %in% colnames(meta)) {
    out$problems <- "Pick a contrast column on the page first."
    return(out)
  }
  if (is.null(ref) || is.null(comp) || !nzchar(ref) || !nzchar(comp)) {
    out$problems <- "Pick the two groups to compare on the page first."
    return(out)
  }
  if (identical(ref, comp)) {
    out$problems <- "The reference and comparison groups must be different."
    return(out)
  }
  v <- as.character(meta[[col]])
  out$n_ref <- sum(v == ref, na.rm = TRUE)
  out$n_comp <- sum(v == comp, na.rm = TRUE)
  if (out$n_ref < min_per_group || out$n_comp < min_per_group) {
    out$problems <- sprintf("At least two comparison groups are required, each with %d or more samples (%s = %d, %s = %d).",
                            min_per_group, ref, out$n_ref, comp, out$n_comp)
    return(out)
  }
  out$ok <- TRUE
  lv <- wf_sex_levels(meta, sex_col)
  out$sex_levels <- lv
  if (length(lv) < 2) {
    out$mode_notes <- c(stratified = "Needs a sex column with at least two values.",
                        sex_specific = "Needs a sex column with at least two values.")
    return(out)
  }
  keep <- v %in% c(ref, comp) & as.character(meta[[sex_col]]) %in% lv
  out$strata <- table(factor(as.character(meta[[sex_col]])[keep], levels = unname(lv)),
                      factor(v[keep], levels = c(ref, comp)))
  small <- out$strata < min_per_group
  out$modes[["stratified"]] <- !any(small)
  if (any(small)) {
    out$mode_notes[["stratified"]] <- sprintf("Each sex needs %d or more samples in both groups.", min_per_group)
  }
  ## The transcriptomics Sex Interaction tab always models the standard "group" column.
  if (identical(layer, "transcriptomics") && !identical(col, "group")) {
    out$mode_notes[["sex_specific"]] <- "The Sex Interaction model uses the \"group\" column; pick that column to use this option."
  } else {
    out$modes[["sex_specific"]] <- !any(out$strata == 0)
    if (any(out$strata == 0)) out$mode_notes[["sex_specific"]] <- "Every sex-by-group cell needs at least one sample."
  }
  out
}

## Step 4 -> 5: the existing analysis runs that implement the chosen sex design.
wf_run_plan <- function(layer, mode, sex_col, sex_levels) {
  kind <- if (identical(layer, "transcriptomics")) "dge" else "dmp"
  if (identical(mode, "sex_specific")) {
    return(list(list(kind = "interaction", label = "Diagnosis-by-sex interaction model (same fit also gives the within-sex effects)")))
  }
  if (identical(mode, "stratified")) {
    nm <- names(sex_levels)
    return(lapply(seq_along(sex_levels), function(i) {
      lab <- if (!is.null(nm) && nzchar(nm[i])) nm[i] else sex_levels[[i]]
      list(kind = kind, sex_level = unname(sex_levels[[i]]), label = sprintf("%s samples only", lab))
    }))
  }
  adjust <- !is.null(sex_col) && length(sex_levels) >= 2
  list(list(kind = kind, sex_level = NULL, adjust_sex = adjust,
            label = if (adjust) "All samples together, adjusted for sex" else "All samples together"))
}

wf_fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

## Step 5 summary: only fields the existing modules actually store.
wf_summary_line <- function(kind, obj, padj = NULL, lfc = NULL) {
  if (is.null(obj)) return(NULL)
  switch(kind,
    dge = sprintf("%s: %s samples, %s genes tested, %s significant (%s up, %s down) at adjusted p < %s and |log2FC| > %s.",
                  obj$contrast, wf_fmt(obj$n_samples), wf_fmt(obj$n_tested), wf_fmt(obj$n_significant),
                  wf_fmt(obj$n_up), wf_fmt(obj$n_down), format(padj %||% "?"), format(lfc %||% "?")),
    dmp = sprintf("%s: %s CpGs tested, %s significant at FDR < 0.05.",
                  obj$comparison, wf_fmt(obj$n_probes), wf_fmt(obj$n_sig)),
    interaction = sprintf("%s: %s features tested, %s with a significant sex-by-diagnosis interaction; %s significant within the reference sex and %s within the comparison sex (adjusted p < %s).",
                          obj$contrast, wf_fmt(obj$n_tested), wf_fmt(obj$n_significant),
                          wf_fmt(obj$n_significant_within_ref_sex), wf_fmt(obj$n_significant_within_comp_sex),
                          format(padj %||% "?")),
    NULL
  )
}

wf_step_state <- function(n, step, done, active) {
  if (!isTRUE(active)) return("future")
  if (n %in% done && n != step) return("done")
  if (n == step) return("current")
  "future"
}

## Home-page tracker: the hero's existing flow-step cards, one per workflow step.
wf_home_tracker_ui <- function(step, done, active) {
  tagList(lapply(WF_STEPS, function(s) {
    st <- wf_step_state(s$n, step, done, active)
    clickable <- isTRUE(active) && (s$n %in% done || s$n == step)
    div(
      class = paste("home-hero-flow-step", paste0("wf-state-", st), if (clickable || (!isTRUE(active) && s$n == 1L)) "wf-clickable"),
      onclick = if (clickable || (!isTRUE(active) && s$n == 1L)) sprintf("Shiny.setInputValue('wf_goto', %d, {priority: 'event'})", s$n),
      title = if (st == "done") "Completed - click to revisit" else NULL,
      div(class = "home-hero-flow-num", if (st == "done") icon("check") else s$n),
      div(
        class = "home-hero-flow-body",
        div(class = "home-hero-flow-title", icon(s$icon), s$label),
        if (st == "current" || (!isTRUE(active) && s$n == 1L)) p(class = "home-hero-flow-blurb", s$blurb)
      )
    )
  }))
}

## Compact tracker in the guide bar: the existing pipeline-summary dots, laid out horizontally.
wf_bar_steps_ui <- function(step, done) {
  tags$ol(
    class = "pipeline-summary-list wf-bar-steps",
    lapply(WF_STEPS, function(s) {
      st <- wf_step_state(s$n, step, done, TRUE)
      clickable <- s$n %in% done || s$n == step
      tags$li(
        class = paste("pipeline-summary-step", paste0("state-", st), if (clickable) "wf-clickable"),
        onclick = if (clickable) sprintf("Shiny.setInputValue('wf_goto', %d, {priority: 'event'})", s$n),
        div(class = "pipeline-summary-dot", if (st == "done") icon("check") else s$n),
        div(class = "pipeline-summary-title", s$label)
      )
    })
  )
}

wf_note <- function(kind, ...) {
  ic <- switch(kind, error = "triangle-exclamation", ok = "circle-check", busy = "spinner", "circle-info")
  div(class = paste("empty-note", paste0("wf-note-", kind)),
      if (identical(kind, "busy")) icon(ic, class = "fa-spin") else icon(ic), " ", ...)
}

wf_counts_text <- function(tbl) {
  if (is.null(tbl) || !length(tbl)) return("none")
  paste(sprintf("%s = %s", names(tbl), wf_fmt(as.integer(tbl))), collapse = ", ")
}

## ---------------------------------------------------------------------------------------------------
## Server. `nav$open(layer, hid)` opens an existing sub-module tab exactly as its "Add" card would,
## `nav$dataset(layer)` opens a layer's Dataset tab.
workflow_guide_server <- function(input, output, session, dataset, results, methyl_dataset, methyl_results, nav) {

  wf <- reactiveValues(
    active = FALSE, step = 1L, done = integer(0), layer = NULL,
    dataset_label = NULL, contrast = NULL, sex_mode = NULL, plan = NULL,
    run_status = "idle", summary = character(0), message = NULL,
    loading = FALSE, complete = FALSE
  )

  L <- function() WF_LAYERS[[wf$layer]]
  pre <- function() if (identical(wf$layer, "transcriptomics")) "tx_" else "mx_"
  mid <- function(hid, id) paste0(pre(), hid, "-", id)
  module_title <- function(hid) {
    reg <- if (identical(wf$layer, "transcriptomics")) TX_MODULES_BY_ID else MX_MODULES_BY_ID
    reg[[hid]]$config$title
  }
  layer_meta <- function() {
    if (identical(wf$layer, "transcriptomics")) dataset$meta else methyl_dataset$sample_sheet
  }
  set_msg <- function(kind, text) wf$message <- list(kind = kind, text = text)

  ## ---- Small action queue: waits for existing inputs to render, sets them, presses buttons. ----
  ## ensure: re-applies `apply()` until `ok()` holds (covers renderUI inputs that appear late).
  ## do:     runs once.  await: polls `check()`, which returns "wait", "ok" or a failure message.
  queue <- reactiveVal(list())
  action_seq <- 0L
  enqueue <- function(...) {
    acts <- lapply(list(...), function(a) { action_seq <<- action_seq + 1L; a$id <- action_seq; a })
    queue(c(queue(), acts))
  }
  fail <- function(msg) {
    queue(list())
    wf$loading <- FALSE
    if (identical(wf$run_status, "running")) wf$run_status <- "failed"
    set_msg("error", msg)
  }
  process_queue <- function() {
    repeat {
      q <- queue()
      if (!length(q)) return(invisible())
      a <- q[[1]]
      if (is.null(a$t0)) { a$t0 <- Sys.time(); q[[1]] <- a; queue(q) }
      res <- switch(a$type,
        do = { a$run(); "ok" },
        ensure = if (isTRUE(a$ok())) "ok" else { a$apply(); "wait" },
        await = a$check()
      )
      ## An action may itself move to another step (which resets the queue): pop only if it is still at the head.
      if (identical(res, "ok")) {
        q2 <- queue()
        if (length(q2) && identical(q2[[1]]$id, a$id)) queue(q2[-1])
        next
      }
      if (identical(res, "wait")) {
        if (as.numeric(difftime(Sys.time(), a$t0, units = "secs")) > (a$timeout %||% 45)) fail(a$fail)
        return(invisible())
      }
      fail(res)
      return(invisible())
    }
  }
  ## Low priority: when a button press and a poll land in the same flush, the module's own observer runs first.
  observe({
    if (!length(queue())) return()
    invalidateLater(300)
    isolate(process_queue())
  }, priority = -100)
  busy <- reactive(length(queue()) > 0)

  ensure_input <- function(id, value, update, fail_msg) {
    list(type = "ensure", ok = function() identical(as.character(input[[id]] %||% ""), as.character(value)),
         apply = function() update(session, id, value), fail = fail_msg)
  }
  upd_select <- function(session, id, value) updateSelectInput(session, id, selected = value)
  upd_radio <- function(session, id, value) updateRadioButtons(session, id, selected = value)
  upd_tabs <- function(session, id, value) updateTabsetPanel(session, id, selected = value)
  keep_open_msg <- function(hid) sprintf("The guide could not set up the %s page. Keep that page open while the guide works, then use the step button again.", module_title(hid))

  ## Presses an existing run button and waits for that module to store a new result.
  click_and_await <- function(btn_id, get_result, expected_same = NULL, on_ok, page) {
    env <- new.env()
    list(
      list(type = "do", run = function() {
        env$clicks <- input[[btn_id]] %||% 0
        env$before <- get_result()
        env$t_click <- Sys.time()
        shinyjs::click(btn_id)
      }),
      list(type = "await", timeout = 3600,
           fail = sprintf("The analysis did not respond. Check the %s page.", page),
           check = function() {
             if ((input[[btn_id]] %||% 0) <= env$clicks) return("wait")
             after <- get_result()
             changed <- !is.null(after) && !identical(after, env$before)
             ## Re-running identical settings stores an identical result; accept it only if it is this run's result.
             same_ok <- !is.null(after) && !is.null(expected_same) && isTRUE(expected_same(after))
             if (changed || same_ok) { on_ok(after); return("ok") }
             errs <- Filter(function(e) identical(e$kind, "error") && isTRUE(e$at >= env$t_click),
                            tryCatch(arthomix_diag_entries(session), error = function(e) list()))
             reason <- if (length(errs)) paste0(" Reported error: ", errs[[length(errs)]]$message) else ""
             sprintf("The analysis did not finish.%s The %s page shows the reason next to its run button. Change the settings there, then press Run analysis again.", reason, page)
           })
    )
  }

  ## ---- Step entry: navigate to the existing page for that step. ----
  enter_step <- function(n) {
    queue(list())
    wf$step <- as.integer(n)
    wf$message <- NULL
    lay <- L()
    switch(as.character(n),
      "1" = nav$dataset(wf$layer),
      "2" = nav$open(wf$layer, lay$qc),
      "3" = open_compare_page(apply_defaults = is.null(wf$contrast)),
      "4" = open_compare_page(apply_defaults = FALSE),
      "5" = open_run_page(),
      "6" = nav$open(wf$layer, lay$candidates),
      "7" = open_validation(lay$validate[[1]])
    )
  }
  complete_step <- function(n) {
    wf$done <- sort(unique(c(wf$done, as.integer(n))))
  }
  clear_after <- function(n) {
    wf$done <- wf$done[wf$done <= n]
    if (n < 5) { wf$run_status <- "idle"; wf$summary <- character(0) }
    wf$complete <- FALSE
  }

  open_validation <- function(v) {
    nav$open(wf$layer, v$id)
    if (!is.null(v$tabset)) enqueue(ensure_input(mid(v$id, v$tabset), v$tab, upd_tabs, keep_open_msg(v$id)))
  }

  ## Contrast inputs on the existing DE / DMP pages.
  contrast_ids <- function() {
    if (identical(wf$layer, "transcriptomics")) {
      list(col = "tx_dge-contrast_col", ref = "tx_dge-ref_group", comp = "tx_dge-comp_group")
    } else {
      list(col = "mx_dmp-live_group_col", ref = "mx_dmp-live_ref", comp = "mx_dmp-live_comp")
    }
  }
  page_contrast <- reactive({
    req(wf$active, wf$layer)
    ids <- contrast_ids()
    list(col = input[[ids$col]], ref = input[[ids$ref]], comp = input[[ids$comp]])
  })

  open_compare_page <- function(apply_defaults) {
    lay <- L()
    nav$open(wf$layer, lay$compare)
    hid <- lay$compare
    if (!is.null(lay$compare_tabset)) enqueue(ensure_input(mid(hid, lay$compare_tabset), lay$compare_tab, upd_tabs, keep_open_msg(hid)))
    if (identical(wf$layer, "transcriptomics")) enqueue(ensure_input("tx_dge-data_source", "pipeline", upd_radio, keep_open_msg(hid)))
    target <- wf$contrast
    if (is.null(target) && apply_defaults) {
      meta <- layer_meta()
      gc <- wf_group_col(meta)
      dc <- if (!is.null(gc)) wf_default_contrast(meta[[gc]]) else NULL
      if (!is.null(dc)) target <- list(col = gc, ref = dc$ref, comp = dc$comp)
    }
    if (!is.null(target)) {
      ids <- contrast_ids()
      enqueue(
        ensure_input(ids$col, target$col, upd_select, keep_open_msg(hid)),
        ensure_input(ids$ref, target$ref, upd_select, keep_open_msg(hid)),
        ensure_input(ids$comp, target$comp, upd_select, keep_open_msg(hid))
      )
    }
  }

  ## Sets the existing sex controls for one planned run (populates the page; does not run it).
  run_settings <- function(run) {
    ct <- wf$contrast
    sc <- wf_sex_col(layer_meta())
    if (identical(run$kind, "interaction")) {
      p <- function(id) mid(L()$interaction, id)
      hid <- L()$interaction
      acts <- list(
        ensure_input(p("ref_group"), ct$ref, upd_select, keep_open_msg(hid)),
        ensure_input(p("comp_group"), ct$comp, upd_select, keep_open_msg(hid))
      )
      if (identical(wf$layer, "methylomics")) acts <- c(list(ensure_input(p("group_col"), ct$col, upd_select, keep_open_msg(hid))), acts)
      return(acts)
    }
    hid <- L()$compare
    if (identical(run$kind, "dge")) {
      base <- list(
        ensure_input("tx_dge-contrast_col", ct$col, upd_select, keep_open_msg(hid)),
        ensure_input("tx_dge-ref_group", ct$ref, upd_select, keep_open_msg(hid)),
        ensure_input("tx_dge-comp_group", ct$comp, upd_select, keep_open_msg(hid))
      )
      sex <- if (!is.null(run$sex_level)) {
        list(ensure_input("tx_dge-covariate_col", sc, upd_select, keep_open_msg(hid)),
             ensure_input("tx_dge-covariate_mode", "filter", upd_radio, keep_open_msg(hid)),
             ensure_input("tx_dge-covariate_level", run$sex_level, upd_select, keep_open_msg(hid)))
      } else if (isTRUE(run$adjust_sex)) {
        list(ensure_input("tx_dge-covariate_col", sc, upd_select, keep_open_msg(hid)),
             ensure_input("tx_dge-covariate_mode", "adjust", upd_radio, keep_open_msg(hid)))
      } else {
        list(ensure_input("tx_dge-covariate_col", "(none)", upd_select, keep_open_msg(hid)))
      }
      return(c(base, sex))
    }
    ## DMP: the Sex radio picks the stratum; pooled adds the sex column as a covariate.
    sex_value <- run$sex_level %||% "__all__"
    acts <- list(
      ensure_input("mx_dmp-live_sex", sex_value, upd_radio, keep_open_msg(hid)),
      ensure_input("mx_dmp-live_group_col", ct$col, upd_select, keep_open_msg(hid)),
      ensure_input("mx_dmp-live_ref", ct$ref, upd_select, keep_open_msg(hid)),
      ensure_input("mx_dmp-live_comp", ct$comp, upd_select, keep_open_msg(hid))
    )
    if (isTRUE(run$adjust_sex)) {
      acts <- c(acts, list(list(
        type = "ensure", fail = keep_open_msg(hid),
        ok = function() sc %in% (input[["mx_dmp-live_covariates"]] %||% character(0)),
        apply = function() updateCheckboxGroupInput(session, "mx_dmp-live_covariates",
                                                    selected = union(input[["mx_dmp-live_covariates"]] %||% character(0), sc))
      )))
    }
    acts
  }

  open_run_page <- function() {
    plan <- wf$plan
    req(length(plan) > 0)
    first <- plan[[1]]
    if (identical(first$kind, "interaction")) {
      nav$open(wf$layer, L()$interaction)
    } else {
      open_compare_page(apply_defaults = FALSE)
    }
    if (!identical(wf$run_status, "done")) do.call(enqueue, run_settings(first))
  }

  run_one <- function(run) {
    lay <- L()
    if (identical(run$kind, "interaction")) {
      hid <- lay$interaction
      btn <- mid(hid, "run_btn")
      getter <- function() if (identical(wf$layer, "transcriptomics")) results$interaction else methyl_results$interaction
      expected <- function(after) identical(after$contrast, sprintf("%s vs %s interaction with sex (%s vs %s)",
        wf$contrast$comp, wf$contrast$ref, input[[mid(hid, "comp_sex")]], input[[mid(hid, "ref_sex")]]))
      padj <- function() input[[mid(hid, "padj_cut")]]
      kind <- "interaction"
    } else if (identical(run$kind, "dge")) {
      hid <- lay$compare
      btn <- "tx_dge-run_btn"
      getter <- function() results$dge_runs
      expected <- NULL
      padj <- function() input[["tx_dge-padj_cut"]]
      kind <- "dge"
    } else {
      hid <- lay$compare
      btn <- "mx_dmp-live_run_btn"
      getter <- function() methyl_results$dmp
      expected <- function(after) identical(after$comparison, sprintf("%s vs %s (%s)", wf$contrast$comp, wf$contrast$ref,
        if (is.null(run$sex_level)) "All samples" else run$sex_level))
      padj <- function() NULL
      kind <- "dmp"
    }
    on_ok <- function(after) {
      obj <- if (identical(kind, "dge")) results$dge else after
      line <- wf_summary_line(kind, obj, padj = padj(), lfc = input[["tx_dge-lfc_cut"]])
      wf$summary <- c(wf$summary, line)
    }
    c(run_settings(run), click_and_await(btn, getter, expected, on_ok, module_title(hid)))
  }

  ## ---- Start: step 1 dialog. ----
  show_start_dialog <- function() {
    tx_src <- dataset$source %||% "none"
    mx_src <- methyl_dataset$source %||% "none loaded yet"
    showModal(modalDialog(
      title = tagList(icon("play"), " Start an analysis"),
      easyClose = TRUE,
      footer = tagList(modalButton("Cancel"), actionButton("wf_start_confirm", "Start", class = "btn-primary")),
      p(class = "submodule-desc", "The guide walks through the existing pages in order: choose data, check samples, define the comparison, choose the sex analysis, run it, explore the biomarkers and validate them. Every page stays usable on its own."),
      radioButtons("wf_layer", "Which omics layer?",
                   choiceNames = list("Transcriptomics (gene expression)", "Methylomics (DNA methylation)"),
                   choiceValues = list("transcriptomics", "methylomics"),
                   selected = wf$layer %||% "transcriptomics"),
      radioButtons("wf_data_choice", "Which data?",
                   choiceNames = list(
                     "Example dataset (loads the bundled rheumatoid arthritis cohort through the Dataset tab)",
                     "Keep the dataset already loaded for that layer (for example your own upload)"
                   ),
                   choiceValues = list("example", "current"), selected = "example"),
      div(class = "empty-note", icon("circle-info"),
          sprintf(" Loaded now - Transcriptomics: %s. Methylomics: %s.", tx_src, mx_src))
    ))
  }

  observeEvent(input$home_start_analysis, show_start_dialog(), ignoreInit = TRUE)

  observeEvent(input$wf_start_confirm, {
    removeModal()
    queue(list())
    wf$active <- TRUE
    wf$layer <- input$wf_layer %||% "transcriptomics"
    wf$done <- integer(0); wf$contrast <- NULL; wf$sex_mode <- NULL; wf$plan <- NULL
    wf$run_status <- "idle"; wf$summary <- character(0); wf$complete <- FALSE; wf$dataset_label <- NULL
    ## With valid data already loaded, go straight to step 2: opening the Dataset tab in the same flush
    ## would win over the QC tab (tab-select messages are applied after tab inserts).
    if (identical(input$wf_data_choice, "current") && dataset_valid()) {
      wf$step <- 1L
      confirm_dataset()
    } else {
      enter_step(1L)
      if (identical(input$wf_data_choice, "current")) confirm_dataset() else load_example()
    }
  }, ignoreInit = TRUE)

  dataset_valid <- function() {
    if (identical(wf$layer, "transcriptomics")) {
      !is.null(dataset$expr) && NCOL(dataset$expr) > 0 && is.data.frame(dataset$meta) && nrow(dataset$meta) > 0
    } else {
      is.data.frame(methyl_dataset$sample_sheet) && nrow(methyl_dataset$sample_sheet) > 0
    }
  }
  confirm_dataset <- function() {
    if (!dataset_valid()) {
      set_msg("error", "Please select a valid dataset before continuing. Load one on this Dataset tab, or restart the guide with the example dataset.")
      return(invisible(FALSE))
    }
    wf$dataset_label <- if (identical(wf$layer, "transcriptomics")) dataset$source else methyl_dataset$source
    complete_step(1L)
    enter_step(2L)
    invisible(TRUE)
  }

  ## Loads the example through the Dataset tab's own "Load this dataset" button.
  load_example <- function() {
    lay <- L()
    p <- if (identical(wf$layer, "transcriptomics")) "tx_dataset-" else "mx_dataset-"
    if (identical(wf$layer, "methylomics")) {
      if (!isTRUE(METH_DATA_AVAILABLE)) {
        set_msg("error", "Please select a valid dataset before continuing. The methylomics example dataset is not installed in this deployment.")
        return(invisible())
      }
      if (isTRUE(methyl_dataset$preloaded) && !is.null(methyl_dataset$beta) && is.data.frame(methyl_dataset$sample_sheet)) {
        ## Already loaded: confirm in a later flush so the QC tab opens after the Dataset tab, not before it.
        enqueue(list(type = "do", run = function() confirm_dataset()))
        return(invisible())
      }
    }
    wf$loading <- TRUE
    set_msg("busy", if (identical(wf$layer, "methylomics"))
      "Loading the example dataset on the Dataset tab. The methylation matrix is about 2 GB, so this can take a few minutes; the rest of the app stays usable."
      else "Loading the example dataset on the Dataset tab.")
    env <- new.env()
    enqueue(
      ensure_input(paste0(p, "preloaded_choice"), lay$example_id, upd_select, "Please select a valid dataset before continuing. The example dataset could not be selected on the Dataset tab."),
      list(type = "do", run = function() {
        env$load_id <- dataset$load_id %||% 0L
        shinyjs::click(paste0(p, "load_preloaded_btn"))
      }),
      list(type = "await", timeout = if (identical(wf$layer, "methylomics")) 1200 else 180,
           fail = "Please select a valid dataset before continuing. The example dataset did not finish loading; the Dataset tab shows the reason.",
           check = function() {
             loaded <- if (identical(wf$layer, "transcriptomics")) {
               (dataset$load_id %||% 0L) > env$load_id
             } else {
               isTRUE(methyl_dataset$preloaded) && is.data.frame(methyl_dataset$sample_sheet)
             }
             if (!loaded) return("wait")
             wf$loading <- FALSE
             confirm_dataset()
             "ok"
           })
    )
  }

  ## A dataset loaded outside the guide replaces the guide's data: go back to checking samples.
  observeEvent(dataset$load_id, {
    if (isTRUE(wf$active) && identical(wf$layer, "transcriptomics") && !isTRUE(wf$loading) && wf$step > 1L) {
      clear_after(1L); wf$contrast <- NULL
      enter_step(2L)
      set_msg("info", "The loaded dataset changed, so the guide went back to Check samples.")
    }
  }, ignoreInit = TRUE)
  observeEvent(methyl_dataset$source, {
    if (isTRUE(wf$active) && identical(wf$layer, "methylomics") && !isTRUE(wf$loading) && wf$step > 1L) {
      clear_after(1L); wf$contrast <- NULL
      enter_step(2L)
      set_msg("info", "The loaded dataset changed, so the guide went back to Check samples.")
    }
  }, ignoreInit = TRUE)

  ## ---- Step 2 check, read from the loaded dataset's own metadata. ----
  sample_check <- reactive({
    req(wf$active, wf$layer)
    meta <- layer_meta()
    wf_check_samples(meta, wf_group_col(meta), wf_sex_col(meta))
  })

  contrast_check <- reactive({
    req(wf$active, wf$layer)
    pc <- page_contrast()
    meta <- layer_meta()
    wf_check_contrast(meta, pc$col, pc$ref, pc$comp, wf_sex_col(meta), wf$layer)
  })

  ## ---- Primary button: complete the current step if its prerequisite holds. ----
  observeEvent(input$wf_next, {
    req(wf$active)
    n <- wf$step
    if (busy()) {
      showNotification("The guide is still setting up this page. Try again in a moment.", type = "message", duration = 4)
      return()
    }
    if (n == 1L) {
      confirm_dataset()
    } else if (n == 2L) {
      chk <- sample_check()
      if (!isTRUE(chk$ok)) { set_msg("error", paste(chk$problems, collapse = " ")); return() }
      complete_step(2L); enter_step(3L)
    } else if (n == 3L) {
      if (identical(wf$layer, "methylomics") && is.null(methyl_dataset$beta)) {
        set_msg("error", "This deployment has only the example dataset's sample sheet, not its live methylation matrix, so the DMP analysis cannot run here.")
        return()
      }
      cc <- contrast_check()
      if (!isTRUE(cc$ok)) { set_msg("error", paste(cc$problems, collapse = " ")); return() }
      pc <- page_contrast()
      new <- list(col = pc$col, ref = pc$ref, comp = pc$comp)
      if (!identical(new, wf$contrast)) { clear_after(3L); wf$sex_mode <- NULL; wf$plan <- NULL }
      wf$contrast <- new
      complete_step(3L); enter_step(4L)
    } else if (n == 4L) {
      mode <- input$wf_sex_mode
      cc <- contrast_check()
      if (is.null(mode) || !isTRUE(cc$modes[[mode]])) { set_msg("error", "Pick one of the available sex analyses first."); return() }
      if (!identical(mode, wf$sex_mode)) clear_after(4L)
      wf$sex_mode <- mode
      meta <- layer_meta()
      wf$plan <- wf_run_plan(wf$layer, mode, wf_sex_col(meta), cc$sex_levels)
      complete_step(4L); enter_step(5L)
    } else if (n == 5L) {
      if (!identical(wf$run_status, "done")) { set_msg("error", "Run the analysis first."); return() }
      complete_step(5L); enter_step(6L)
    } else if (n == 6L) {
      complete_step(6L); enter_step(7L)
    } else if (n == 7L) {
      complete_step(7L); wf$complete <- TRUE
      set_msg("ok", "Guided analysis complete. Every page keeps its results; you can keep working on any of them.")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$wf_run, {
    req(wf$active, identical(wf$step, 5L))
    if (identical(wf$run_status, "running")) return()
    plan <- wf$plan
    req(length(plan) > 0)
    ## Any page set-up still pending stays in the queue and finishes first; the runs are appended after it.
    if (identical(plan[[1]]$kind, "interaction")) nav$open(wf$layer, L()$interaction) else open_compare_page(apply_defaults = FALSE)
    wf$summary <- character(0)
    wf$run_status <- "running"
    set_msg("busy", sprintf("Running %d analysis%s with the page's own run button.", length(plan), if (length(plan) > 1) " runs" else ""))
    acts <- unlist(lapply(plan, run_one), recursive = FALSE)
    acts <- c(acts, list(list(type = "do", run = function() {
      wf$run_status <- "done"
      set_msg("ok", "Analysis finished. The results are on the page below; the summary lists what was stored.")
    })))
    do.call(enqueue, acts)
  }, ignoreInit = TRUE)

  observeEvent(input$wf_back, {
    req(wf$active, wf$step > 1L)
    if (identical(wf$run_status, "running")) return()
    enter_step(wf$step - 1L)
  }, ignoreInit = TRUE)

  observeEvent(input$wf_exit, {
    queue(list())
    wf$active <- FALSE; wf$step <- 1L; wf$done <- integer(0); wf$message <- NULL
    wf$loading <- FALSE; wf$run_status <- "idle"; wf$summary <- character(0); wf$complete <- FALSE
  }, ignoreInit = TRUE)

  observeEvent(input$wf_goto, {
    n <- as.integer(input$wf_goto)
    if (!isTRUE(wf$active)) {
      if (identical(n, 1L)) show_start_dialog()
      return()
    }
    if (identical(wf$run_status, "running") || isTRUE(wf$loading)) return()
    if (n %in% wf$done || n == wf$step) {
      if (n == 1L) show_start_dialog() else enter_step(n)
    } else {
      showNotification("Finish the earlier steps first.", type = "warning", duration = 4)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$wf_open_candidates, { req(wf$active); nav$open(wf$layer, L()$candidates) }, ignoreInit = TRUE)
  observeEvent(input$wf_open_card, { req(wf$active); nav$open(wf$layer, L()$card) }, ignoreInit = TRUE)
  observeEvent(input$wf_open_validation, {
    req(wf$active)
    v <- Find(function(x) identical(x$id, input$wf_open_validation), L()$validate)
    req(v)
    open_validation(v)
  }, ignoreInit = TRUE)

  ## ---- Rendering. ----
  output$wf_home_tracker <- renderUI({
    wf_home_tracker_ui(wf$step, wf$done, wf$active)
  })

  output$wf_bar <- renderUI({
    if (!isTRUE(wf$active)) return(NULL)
    n <- wf$step
    s <- WF_STEPS[[n]]
    lay <- L()
    is_busy <- busy()
    msg <- wf$message
    body <- switch(as.character(n),
      "1" = tagList(
        p(class = "submodule-desc", "The Dataset tab below is the existing one. The guide loads the example dataset with its own Load button; to use other data, load it there and press Continue."),
        if (!is.null(wf$dataset_label)) wf_note("ok", "Loaded: ", wf$dataset_label)
      ),
      "2" = {
        chk <- sample_check()
        tagList(
          p(class = "submodule-desc", sprintf("The %s page below is the existing quality-control page. These counts come from the loaded dataset's sample metadata.", module_title(lay$qc))),
          wf_note(if (chk$ok) "ok" else "error",
                  sprintf("%s samples. Groups (%s): %s. Sex (%s): %s.",
                          wf_fmt(chk$n), chk$group_col %||% "no column", wf_counts_text(chk$group_counts),
                          chk$sex_col %||% "no column", wf_counts_text(chk$sex_counts))),
          lapply(chk$problems, function(t) wf_note("error", t)),
          lapply(chk$notes, function(t) wf_note("info", t))
        )
      },
      "3" = {
        pc <- page_contrast()
        cc <- contrast_check()
        tagList(
          p(class = "submodule-desc", sprintf("Set the comparison on the %s page below. The guide pre-selects the group column and puts a control-like group first as the reference.", module_title(lay$compare))),
          if (!is.null(pc$col)) wf_note(if (cc$ok) "ok" else "error",
            sprintf("Contrast column: %s. Reference: %s (n = %s). Comparison: %s (n = %s).",
                    pc$col, pc$ref %||% "?", wf_fmt(cc$n_ref), pc$comp %||% "?", wf_fmt(cc$n_comp)))
          else wf_note("busy", "Waiting for the page's comparison controls to appear."),
          if (!cc$ok && !is.null(pc$col)) lapply(cc$problems, function(t) wf_note("error", t))
        )
      },
      "4" = {
        cc <- contrast_check()
        avail <- names(cc$modes)[cc$modes]
        desc <- c(
          pooled = "Pooled: all samples in one model, with sex as an adjustment covariate when a sex column exists.",
          stratified = "Stratified: the same comparison run separately in female-only and male-only samples.",
          sex_specific = "Sex-specific: the existing diagnosis-by-sex interaction model, which tests whether the disease effect differs between sexes. A feature significant in one sex but not the other is not by itself sex-specific."
        )
        sel <- intersect(c(isolate(input$wf_sex_mode), wf$sex_mode, "pooled"), avail)[1]
        tagList(
          radioButtons("wf_sex_mode", "Sex analysis", inline = TRUE,
                       choices = stats::setNames(avail, WF_SEX_MODES[avail]), selected = sel),
          tags$ul(class = "wf-mode-list", lapply(names(desc), function(m) tags$li(
            desc[[m]], if (!isTRUE(cc$modes[[m]])) tags$em(paste(" Not available:", if (m %in% names(cc$mode_notes)) cc$mode_notes[[m]] else "check the comparison on the page."))))),
          if (!isTRUE(cc$ok)) lapply(cc$problems, function(t) wf_note("error", t))
        )
      },
      "5" = tagList(
        p(class = "submodule-desc", "Run analysis presses the page's own run button with these settings, so the results are identical to running it by hand."),
        tags$ul(class = "wf-mode-list", lapply(wf$plan, function(r) tags$li(r$label))),
        if (length(wf$summary)) div(class = "wf-summary", lapply(wf$summary, function(t) wf_note("ok", t)))
      ),
      "6" = tagList(
        p(class = "submodule-desc", sprintf("Candidate %s and the Biomarker Card are the existing pages. No candidates is a valid result, not an error; you can still continue to validation.", lay$feature)),
        div(class = "wf-bar-links",
            actionButton("wf_open_candidates", module_title(lay$candidates), icon = icon("star"), class = "btn btn-default btn-sm"),
            actionButton("wf_open_card", module_title(lay$card), icon = icon("id-card"), class = "btn btn-default btn-sm"))
      ),
      "7" = tagList(
        p(class = "submodule-desc", "Validation runs on the existing page below, against an independent cohort."),
        div(class = "wf-bar-links", lapply(lay$validate, function(v) {
          lab <- if (!is.null(v$tab)) sprintf("%s (%s)", v$tab, module_title(v$id)) else module_title(v$id)
          tags$button(type = "button", class = "btn btn-default btn-sm", lab,
                      onclick = sprintf("Shiny.setInputValue('wf_open_validation', '%s', {priority: 'event'})", v$id))
        }))
      )
    )
    primary <- if (n == 5L && !identical(wf$run_status, "done")) {
      actionButton("wf_run", "Run analysis", icon = icon("play"), class = "btn btn-primary btn-sm",
                   disabled = if (identical(wf$run_status, "running")) "disabled" else NULL)
    } else if (n == 7L && isTRUE(wf$complete)) {
      NULL
    } else {
      actionButton("wf_next", if (n == 7L) "Finish" else "Continue", icon = icon(if (n == 7L) "flag-checkered" else "arrow-right"),
                   class = "btn btn-primary btn-sm", disabled = if (is_busy) "disabled" else NULL)
    }
    div(
      class = "card wf-bar",
      div(
        class = "wf-bar-head",
        div(class = "card-title", icon("route"), sprintf(" Start an analysis · %s · Step %d of %d: %s", lay$label, n, WF_N_STEPS, s$label)),
        actionLink("wf_exit", tagList(icon("xmark"), " Exit guide"), class = "wf-bar-exit")
      ),
      wf_bar_steps_ui(n, wf$done),
      body,
      if (!is.null(msg)) wf_note(msg$kind, msg$text),
      div(
        class = "wf-bar-actions",
        if (n > 1L) actionButton("wf_back", "Back", icon = icon("arrow-left"), class = "btn btn-default btn-sm",
                                 disabled = if (identical(wf$run_status, "running")) "disabled" else NULL),
        if (n == 5L && identical(wf$run_status, "done")) actionButton("wf_run", "Run again", icon = icon("play"), class = "btn btn-default btn-sm"),
        primary
      )
    )
  })

  invisible(wf)
}
