## R/workflow_guide.R
## "Start an analysis" guided workflow (external reviewer request, 2026-09-22).
## A thin controller over existing sub-modules: it opens the existing tabs, sets their existing inputs,
## presses their existing run buttons and reads their existing result objects. It computes no statistics.

WF_STEPS <- list(
  list(n = 1L, label = "Choose data",         icon = "database",     blurb = "Pick transcriptomics or methylomics and load a dataset."),
  list(n = 2L, label = "Check samples",       icon = "users",        blurb = "See the cohort composition - group and sex counts - before picking an analysis design."),
  list(n = 3L, label = "Choose analysis",     icon = "list-check",   blurb = "Pick which analyses to run, confirm the comparison (or change it) and the sex design, then lock it in."),
  list(n = 4L, label = "Run analysis",        icon = "play",         blurb = "Run the chosen analyses on their existing pages."),
  list(n = 5L, label = "Explore biomarkers",  icon = "star",         blurb = "Narrow the results to a candidate feature set and see what it is enriched for."),
  list(n = 6L, label = "Biomarker modeling",  icon = "sliders",      blurb = "Turn the candidates into a consensus panel and train a classifier on it."),
  list(n = 7L, label = "Validate",            icon = "flask-vial",   blurb = "Check the panel in an independent cohort, another tissue and another ancestry."),
  list(n = 8L, label = "Interpret biomarkers", icon = "id-card",     blurb = "Bring the dataset evidence and every saved result together on the Biomarker Card.")
)
WF_N_STEPS <- length(WF_STEPS)

WF_STEP_LABEL <- vapply(WF_STEPS, `[[`, character(1), "label")

## Where each step lands. Ids are existing sub-module ids (TX_MODULES / MX_MODULES).
WF_LAYERS <- list(
  transcriptomics = list(
    label = "Transcriptomics", example_id = "__default_merged__",
    qc = "overview", compare = "dge", interaction = "interaction",
    candidates = "candidates", card = "biomarkercard",
    explore = list(list(id = "candidates"), list(id = "enrichment")),
    model = list(list(id = "featureselection"), list(id = "diagnostic")),
    validate = list(
      list(id = "diagnostic", tabset = "main_tabs", tab = "External Validation"),
      list(id = "crosstissue"),
      list(id = "crossancestry")
    ),
    feature = "genes"
  ),
  methylomics = list(
    label = "Methylomics", example_id = "gse42861_wholeblood",
    qc = "qc", compare = "dmp", compare_tabset = "dmp_subtabs", compare_tab = "DMP", interaction = "interaction",
    candidates = "candidates", card = "biomarkercard",
    explore = list(list(id = "candidates")),
    model = list(list(id = "featureselection"), list(id = "diagnostic")),
    validate = list(list(id = "validation")),
    feature = "CpGs"
  )
)

## ---------------------------------------------------------------------------------------------------
## Step map. One entry per navigation slot, so a module that belongs to two steps has two entries with
## the same `id` (one tab) and different `key`s (two sidebar entries, two statuses). Verified against the
## modules' own upstream reads on 2026-09-22 - see the comment on each line that differs from the obvious
## reading of the module's group.
##
##   hard   conjunctive normal form: a list of any-of groups, all of which must hold. A group is a
##          character vector; "dataset", "contrast" and "sex" are settings tokens, anything else is a key.
##   soft   inputs used when present. They never lock a module, but they order the run and they make it
##          stale when they are re-run.
##   slot   the results$<slot> the module writes. `subfield` narrows that to one field inside it.
##   fp     which shared settings this module's result depends on (see wf_fingerprint()).
##   alias  TRUE for the second entry of a two-step module: it gets a sidebar row but no catalogue card,
##          because a second card would duplicate the smcard_<id> element ids.
## ---------------------------------------------------------------------------------------------------

wf_step_entry <- function(key, id, step, hard = list(), soft = character(0), slot = id,
                          subfield = NULL, fp = character(0), tabset = NULL, tab = NULL,
                          alias = FALSE, optional = FALSE, note = NULL) {
  list(key = key, id = id, step = as.integer(step), hard = hard, soft = soft, slot = slot,
       subfield = subfield, fp = fp, tabset = tabset, tab = tab,
       alias = isTRUE(alias), optional = isTRUE(optional), note = note)
}

TX_STEP_MAP <- local({
  e <- list(
    wf_step_entry("overview", "overview", 2, hard = list("dataset"), fp = "dataset"),
    ## Preprocessing writes results$study_meta, not results$preprocessing.
    wf_step_entry("preprocessing", "preprocessing", 2, hard = list("dataset"), slot = "study_meta",
                  fp = "dataset", optional = TRUE,
                  note = "Suggested when Overview shows more than one source dataset or batch."),
    wf_step_entry("dge", "dge", 4, hard = list("dataset", "contrast"), soft = "preprocessing",
                  slot = "dge_runs", fp = c("dataset", "contrast", "sex_mode")),
    wf_step_entry("wgcna", "wgcna", 4, hard = list("dataset"), soft = "preprocessing", fp = "dataset"),
    ## mr and coloc run off bundled eQTL/GWAS summary statistics; they need no dataset.
    wf_step_entry("mr", "mr", 4),
    wf_step_entry("coloc", "coloc", 4),
    wf_step_entry("interaction", "interaction", 4, hard = list("dataset", "sex"),
                  fp = c("dataset", "contrast")),
    wf_step_entry("deconvolution", "deconvolution", 4, hard = list("dataset"), fp = "dataset"),
    ## mod_candidates.R: module_choices() req()s results$wgcna and deg_run_choices() req()s
    ## results$dge_runs; the gene set is the WGCNA module genes intersected with the DEGs.
    ## mr/coloc only narrow it further (causal_refine_ui), so they are soft - but they must run first.
    wf_step_entry("candidates", "candidates", 5, hard = list("dge", "wgcna"), soft = c("mr", "coloc"),
                  fp = c("dataset", "contrast", "sex_mode")),
    ## Enrichment runs off a pasted gene list; the three preset buttons only prefill that textarea from
    ## results$featureselection / results$crossancestry, each with a bundled fallback. Nothing is required.
    wf_step_entry("enrichment", "enrichment", 5, soft = c("featureselection", "crossancestry")),
    wf_step_entry("featureselection", "featureselection", 6, hard = list("candidates"), soft = "wgcna",
                  fp = c("dataset", "contrast", "sex_mode")),
    wf_step_entry("diagnostic", "diagnostic", 6, hard = list(c("featureselection", "candidates")),
                  fp = c("dataset", "contrast", "sex_mode")),
    ## Same tab, opened on its External Validation sub-tab; it needs a trained model and writes into
    ## results$diagnostic[[sex]]$external.
    wf_step_entry("diagnostic_external", "diagnostic", 7, hard = list("diagnostic"),
                  slot = "diagnostic", subfield = "external", fp = c("dataset", "contrast", "sex_mode"),
                  tabset = "main_tabs", tab = "External Validation", alias = TRUE),
    wf_step_entry("crosstissue", "crosstissue", 7, hard = list("featureselection"),
                  soft = c("diagnostic", "dge"), fp = c("dataset", "sex_mode")),
    ## Cross-Ancestry MR replicates off bundled MR35_crossancestry_*.csv (or an upload); it reads no
    ## upstream results at all, so it is never locked.
    wf_step_entry("crossancestry", "crossancestry", 7),
    wf_step_entry("biomarkercard", "biomarkercard", 8, hard = list("dataset"),
                  soft = c("dge", "candidates", "featureselection", "diagnostic", "crosstissue"),
                  fp = c("dataset", "contrast", "sex_mode"))
  )
  stats::setNames(e, vapply(e, `[[`, character(1), "key"))
})

## Registry order, used to break ties inside a step so tabs and cards keep the catalogue's own order.
TX_STEP_MAP_ORDER <- names(TX_STEP_MAP)

WF_STATUS_LABEL <- c(locked = "Locked", ready = "Ready", done = "Done", stale = "Stale")
WF_STATUS_ORDER <- c(locked = 1L, ready = 2L, stale = 3L, done = 4L)

## A step's rolled-up status: stale beats done, ready beats locked, and a step with nothing in it
## (Define comparison, Choose analysis) reports whether the wizard has settled that setting.
wf_step_status <- function(statuses) {
  if (!length(statuses)) return("ready")
  if (any(statuses == "stale")) return("stale")
  if (all(statuses == "done")) return("done")
  if (all(statuses == "locked")) return("locked")
  "ready"
}

WF_SETTING_LABEL <- c(
  dataset  = "a loaded dataset",
  contrast = "a comparison (Step 3)",
  sex      = "a sex column with at least two values"
)

## Human label for a requirement token: a settings token, or a module's own title.
wf_req_label <- function(token, titles = NULL) {
  if (token %in% names(WF_SETTING_LABEL)) return(unname(WF_SETTING_LABEL[[token]]))
  ttl <- titles[[token]]
  if (!is.null(ttl)) ttl else token
}

## Titles for TX_STEP_MAP keys, taken from the modules' own configs (never hard-coded here).
wf_step_titles <- function(map = TX_STEP_MAP, registry = NULL) {
  if (is.null(registry) && exists("TX_MODULES_BY_ID")) registry <- get("TX_MODULES_BY_ID")
  out <- vapply(map, function(e) {
    ttl <- tryCatch(registry[[e$id]]$config$title, error = function(err) NULL)
    base <- ttl %||% e$id
    if (!is.null(e$tab)) sprintf("%s (%s)", e$tab, base) else base
  }, character(1))
  stats::setNames(out, names(map))
}

wf_and_list <- function(x) {
  x <- x[nzchar(x)]
  if (length(x) <= 1) return(paste(x, collapse = ""))
  paste(paste(utils::head(x, -1), collapse = ", "), "and", utils::tail(x, 1))
}
wf_or_list <- function(x) {
  x <- x[nzchar(x)]
  if (length(x) <= 1) return(paste(x, collapse = ""))
  paste(paste(utils::head(x, -1), collapse = ", "), "or", utils::tail(x, 1))
}

## Does results carry a stored result for this entry? `subfield` looks one level into the per-stratum
## lists the Diagnostic module writes (results$diagnostic[[sex]]$external).
wf_has_result <- function(entry, results) {
  v <- tryCatch(results[[entry$slot]], error = function(e) NULL)
  if (is.null(v)) return(FALSE)
  if (is.null(entry$subfield)) return(length(v) > 0)
  if (!is.list(v)) return(FALSE)
  any(vapply(v, function(x) is.list(x) && !is.null(x[[entry$subfield]]), logical(1)))
}

## The fingerprint a stored result is recorded with: the shared settings it used plus the run stamps of
## the upstream results it could have read. Re-running an upstream bumps its stamp, which is what makes
## everything downstream stale.
wf_fingerprint <- function(key, dataset, settings, map = TX_STEP_MAP) {
  e <- map[[key]]
  if (is.null(e)) return(NULL)
  stamps <- settings$stamps %||% list()
  up <- c(unlist(e$hard, use.names = FALSE), e$soft)
  up <- intersect(up, names(map))
  list(
    dataset  = if ("dataset" %in% e$fp) dataset$load_id %||% NA else NULL,
    contrast = if ("contrast" %in% e$fp) {
      ct <- settings$contrast
      if (is.null(ct)) NA_character_ else paste(ct$col %||% "", ct$ref %||% "", ct$comp %||% "", sep = "|")
    } else NULL,
    sex_mode = if ("sex_mode" %in% e$fp) settings$sex_mode %||% NA_character_ else NULL,
    upstream = if (length(up)) {
      vals <- lapply(up, function(k) stamps[[k]] %||% NA_integer_)
      stats::setNames(vals, up)
    } else NULL
  )
}

## "locked" | "ready" | "done" | "stale", with a readable reason.
##   dataset  list(ok =, sex_ok =, load_id =)
##   results  the shared results store (reactiveValues or a plain list)
##   settings list(contrast =, sex_mode =, stamps =, fingerprints =)
wf_module_status <- function(key, dataset = list(), results = list(), settings = list(),
                             map = TX_STEP_MAP, titles = NULL, .seen = character(0)) {
  e <- map[[key]]
  if (is.null(e)) return(list(status = "ready", reason = ""))
  if (key %in% .seen) return(list(status = "ready", reason = ""))   # cycle guard; the map has none
  if (is.null(titles)) titles <- wf_step_titles(map)
  .seen <- c(.seen, key)

  dep_status <- function(k) {
    wf_module_status(k, dataset, results, settings, map, titles, .seen)$status
  }
  satisfied <- function(token) {
    switch(token,
      dataset  = isTRUE(dataset$ok),
      contrast = !is.null(settings$contrast),
      sex      = isTRUE(dataset$sex_ok),
      dep_status(token) %in% c("done", "stale")   # a stale upstream still exists and can be read
    )
  }

  missing <- Filter(Negate(is.null), lapply(e$hard, function(grp) {
    if (any(vapply(grp, satisfied, logical(1)))) NULL
    else wf_or_list(vapply(grp, wf_req_label, character(1), titles = titles))
  }))
  if (length(missing)) {
    return(list(status = "locked",
                reason = paste0("Needs: ", wf_and_list(unlist(missing, use.names = FALSE)), ".")))
  }

  if (!wf_has_result(e, results)) {
    return(list(status = "ready", reason = "Ready to run."))
  }

  recorded <- (settings$fingerprints %||% list())[[key]]
  if (!is.null(recorded)) {
    expected <- wf_fingerprint(key, dataset, settings, map)
    if (!identical(recorded, expected)) {
      return(list(status = "stale", reason = wf_stale_reason(recorded, expected, titles)))
    }
  }

  ups <- intersect(c(unlist(e$hard, use.names = FALSE), e$soft), names(map))
  stale_up <- Filter(function(k) identical(dep_status(k), "stale"), ups)
  if (length(stale_up)) {
    return(list(status = "stale",
                reason = sprintf("Stale - re-run. %s changed since this ran.",
                                 wf_and_list(vapply(stale_up, wf_req_label, character(1), titles = titles)))))
  }
  list(status = "done", reason = "Up to date.")
}

wf_stale_reason <- function(recorded, expected, titles = NULL) {
  bits <- character(0)
  if (!identical(recorded$dataset, expected$dataset)) bits <- c(bits, "the loaded dataset")
  if (!identical(recorded$contrast, expected$contrast)) bits <- c(bits, "the comparison")
  if (!identical(recorded$sex_mode, expected$sex_mode)) bits <- c(bits, "the sex design")
  changed_up <- names(expected$upstream)[
    !vapply(names(expected$upstream), function(k) identical(recorded$upstream[[k]], expected$upstream[[k]]), logical(1))
  ]
  if (length(changed_up)) bits <- c(bits, vapply(changed_up, wf_req_label, character(1), titles = titles))
  if (!length(bits)) bits <- "an upstream result"
  sprintf("Stale - re-run. %s changed since this ran.", wf_and_list(bits))
}

## The selection plus everything it hard-requires, with a reason for each auto-added module.
## For an any-of group the already-selected alternative wins, otherwise the first one listed.
wf_required_closure <- function(selected_ids, map = TX_STEP_MAP, titles = NULL) {
  if (is.null(titles)) titles <- wf_step_titles(map)
  keys <- intersect(TX_STEP_MAP_ORDER, unique(selected_ids))
  added <- list()
  repeat {
    grown <- FALSE
    for (k in keys) {
      for (grp in map[[k]]$hard) {
        grp <- intersect(grp, names(map))
        if (!length(grp) || any(grp %in% keys)) next
        pick <- grp[1]
        keys <- c(keys, pick)
        added[[pick]] <- sprintf("Added automatically: %s needs %s.", titles[[k]], titles[[pick]])
        grown <- TRUE
      }
    }
    keys <- intersect(TX_STEP_MAP_ORDER, unique(keys))
    if (!grown) break
  }
  list(keys = keys, added = added)
}

## Selected modules in dependency order: hard edges always, soft edges when both ends are selected.
## Ties break by step, then by registry order, so the result is deterministic.
wf_run_order <- function(selected_ids, map = TX_STEP_MAP) {
  keys <- intersect(TX_STEP_MAP_ORDER, unique(selected_ids))
  if (!length(keys)) return(character(0))
  rank <- stats::setNames(seq_along(TX_STEP_MAP_ORDER), TX_STEP_MAP_ORDER)
  deps <- lapply(keys, function(k) {
    d <- intersect(c(unlist(map[[k]]$hard, use.names = FALSE), map[[k]]$soft), keys)
    ## Never let a soft edge point backwards into a cycle.
    setdiff(d, k)
  })
  names(deps) <- keys
  out <- character(0)
  remaining <- keys
  while (length(remaining)) {
    free <- Filter(function(k) !length(setdiff(deps[[k]], out)), remaining)
    if (!length(free)) {            # unreachable with this map; keep it total anyway
      free <- remaining[order(vapply(remaining, function(k) map[[k]]$step, integer(1)), rank[remaining])]
      free <- free[1]
    }
    free <- free[order(vapply(free, function(k) map[[k]]$step, integer(1)), rank[free])]
    out <- c(out, free[1])
    remaining <- setdiff(remaining, free[1])
  }
  out
}

## Step 3 presets. Sex-specific adds Sex Interaction and still runs the main DGE alongside it.
TX_PRESETS <- list(
  standard = list(
    label = "Standard biomarker discovery", sex_mode = "pooled",
    blurb = "DGE and WGCNA to candidates, a consensus panel, a classifier, its external cohort and the card.",
    keys = c("overview", "dge", "wgcna", "candidates", "featureselection", "diagnostic",
             "diagnostic_external", "biomarkercard")
  ),
  sex_differences = list(
    label = "Sex differences in RA", sex_mode = "sex_specific",
    blurb = "The standard route plus the diagnosis-by-sex interaction model, per-sex candidates and panels, and the cross-tissue check.",
    keys = c("overview", "preprocessing", "dge", "wgcna", "interaction", "candidates",
             "featureselection", "diagnostic", "diagnostic_external", "crosstissue", "biomarkercard")
  ),
  full = list(
    label = "Full", sex_mode = "pooled",
    blurb = "Every analysis on this page, in pipeline order.",
    keys = TX_STEP_MAP_ORDER
  ),
  custom = list(
    label = "Custom", sex_mode = NULL,
    blurb = "Tick the analyses yourself; anything they need is added for you.",
    keys = character(0)
  )
)

## A preset resolved against the current sex design: the closure of its modules, in run order.
wf_preset_plan <- function(preset, sex_mode = NULL, extra = character(0), map = TX_STEP_MAP, titles = NULL) {
  p <- TX_PRESETS[[preset]] %||% TX_PRESETS$custom
  mode <- sex_mode %||% p$sex_mode
  keys <- unique(c(p$keys, extra))
  if (identical(mode, "sex_specific")) keys <- unique(c(keys, "interaction", "dge"))
  cl <- wf_required_closure(keys, map, titles)
  list(preset = preset, sex_mode = mode, keys = wf_run_order(cl$keys, map), added = cl$added)
}

## The picker's four layers. Only the two with a WF_LAYERS entry can be driven step by step.
WF_PICK_LAYERS <- list(
  list(id = "transcriptomics", label = "Transcriptomics", icon = "dna"),
  list(id = "methylomics",     label = "Methylomics",     icon = "circle-nodes"),
  list(id = "crossomics",      label = "Cross-Omics",     icon = "arrows-left-right"),
  list(id = "multiomics",      label = "Multi-Omics",     icon = "layer-group")
)

## Problems the data can have before anything will run, with what to do about each. The text is the
## hover, so the fix travels with the error instead of needing a separate panel.
WF_FIX_HINTS <- c(
  "no group column" = "Open the Dataset tab, load the data again and map the diagnosis or group column. Without it there is nothing to compare.",
  "At least two comparison groups" = "The loaded data has only one group, or too few samples in one. Load a dataset with at least two groups of three or more samples on the Dataset tab.",
  "valid dataset" = "Nothing usable is loaded. Go to the Dataset tab and load the example cohort, a GEO series, or your own matrix.",
  "sex column" = "Sex-based analyses need a sex column with at least two values. Map it on the Dataset tab, or choose the pooled analysis."
)
wf_fix_hint <- function(text) {
  hit <- names(WF_FIX_HINTS)[vapply(names(WF_FIX_HINTS), function(k) grepl(k, text, fixed = TRUE), logical(1))]
  if (length(hit)) unname(WF_FIX_HINTS[[hit[1]]]) else NULL
}

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

## Step 3: is the chosen comparison runnable, and which sex designs does it support?
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

## Step 3 -> 4: the existing analysis runs that implement the chosen sex design.
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

## Step 4 summary: only fields the existing modules actually store.
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
        ## Only the current step's blurb shows, so all nine rows stay the same height before the
        ## guide starts (i.e. no state is "current" yet) instead of step 1 alone standing taller.
        if (st == "current") p(class = "home-hero-flow-blurb", s$blurb)
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

wf_note <- function(kind, ..., hint = NULL) {
  ic <- switch(kind, error = "triangle-exclamation", ok = "circle-check", busy = "spinner", "circle-info")
  txt <- paste(unlist(list(...)), collapse = "")
  hint <- hint %||% if (identical(kind, "error")) wf_fix_hint(txt) else NULL
  div(class = paste("empty-note", paste0("wf-note-", kind), if (!is.null(hint)) "wf-note-hinted"),
      title = hint,
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
    loading = FALSE, complete = FALSE,
    preset = "standard", selected = TX_PRESETS$standard$keys,
    configured = FALSE, bar_hidden = FALSE
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
      "1" = if (!is.null(nav$design) && identical(wf$layer, "transcriptomics")) nav$design(wf$layer)
            else nav$dataset(wf$layer),
      "2" = if (!is.null(nav$design) && identical(wf$layer, "transcriptomics")) nav$design(wf$layer)
            else nav$open(wf$layer, lay$qc),
      "3" = {
        ## Choose analysis also confirms the comparison and locks the sex design, so the DE/DMP
        ## module is opened in the background here already (its inputs exist as soon as the tab is
        ## inserted) while the guide stays put. The step body links to it for anyone who wants to
        ## change the comparison themselves.
        if (!is.null(nav$design) && identical(wf$layer, "transcriptomics")) {
          if (!is.null(nav$mount)) nav$mount(wf$layer, lay$compare)
          nav$design(wf$layer)
        } else open_compare_page(apply_defaults = FALSE)
      },
      "4" = open_run_page(),
      "5" = nav$open(wf$layer, lay$explore[[1]]$id),
      "6" = nav$open(wf$layer, lay$model[[1]]$id),
      "7" = open_validation(lay$validate[[1]]),
      "8" = nav$open(wf$layer, lay$card)
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

  ## Buttons for the modules a step covers. One generic input so every step uses the same route,
  ## and a module with an inner tab (External Validation) carries it in the payload.
  step_links_ui <- function(entries) {
    if (!length(entries)) return(NULL)
    div(class = "wf-bar-links", lapply(entries, function(v) {
      lab <- if (!is.null(v$tab)) sprintf("%s (%s)", v$tab, module_title(v$id)) else module_title(v$id)
      payload <- paste(v$id, v$tabset %||% "", v$tab %||% "", sep = "|")
      tags$button(
        type = "button", class = "btn btn-default btn-sm", lab,
        onclick = sprintf("Shiny.setInputValue('wf_open_module', '%s', {priority: 'event'})",
                          gsub("'", "\\\\'", payload))
      )
    }))
  }

  ## ---- Step 3: which analyses to run (transcriptomics). ----
  ## Nothing runs here. The preset picks a set of existing sub-modules; wf_required_closure() adds
  ## whatever they cannot run without, and wf_run_order() puts them in dependency order.
  selection_plan <- function() {
    preset <- input$wf_preset %||% "standard"
    extra <- if (identical(preset, "custom")) (input$wf_analyses %||% character(0)) else character(0)
    wf_preset_plan(preset, sex_mode = input$wf_sex_mode %||% wf$sex_mode, extra = extra,
                   titles = wf_step_titles())
  }

  preset_picker_ui <- function() {
    titles <- wf_step_titles()
    preset <- isolate(input$wf_preset) %||% wf$preset %||% "standard"
    catalogue <- Filter(function(e) !e$alias, TX_STEP_MAP)
    tagList(
      radioButtons(
        "wf_preset", "Transcriptomics analysis design",
        choiceNames = lapply(names(TX_PRESETS), function(k) tagList(
          strong(TX_PRESETS[[k]]$label), tags$span(class = "wf-preset-blurb", TX_PRESETS[[k]]$blurb))),
        choiceValues = as.list(names(TX_PRESETS)), selected = preset
      ),
      conditionalPanel(
        condition = "input.wf_preset == 'custom'",
        checkboxGroupInput(
          "wf_analyses", NULL,
          choices = stats::setNames(vapply(catalogue, `[[`, character(1), "key"),
                                    vapply(catalogue, function(e) sprintf("Step %d - %s", e$step, titles[[e$key]]), character(1))),
          selected = isolate(input$wf_analyses) %||% wf$selected
        )
      )
    )
  }

  preset_plan_ui <- function() {
    plan <- selection_plan()
    titles <- wf_step_titles()
    if (!length(plan$keys)) return(wf_note("info", "Tick at least one analysis above."))
    ## One numbered pipeline, starting where the data does - the run order, plain and sequential.
    ## (This used to also print which of the eight wizard steps each row belonged to, but several
    ## rows share a step - e.g. DGE and WGCNA are both "Step 4" - so that label just repeated or
    ## jumped from row to row and read as broken rather than informative. Dropped.)
    items <- c(list(list(n = "1", title = "Load the dataset")),
               lapply(seq_along(plan$keys), function(i) {
                 k <- plan$keys[[i]]
                 list(n = as.character(i + 1L), title = titles[[k]])
               }))
    tagList(
      tags$ul(class = "wf-plan-list", lapply(items, function(it) tags$li(
        tags$span(class = "wf-plan-num", it$n),
        tags$span(class = "wf-plan-title", it$title)))),
      lapply(unname(plan$added), function(t) wf_note("info", t))
    )
  }

  ## The Custom Analysis Design page. Before setup it is the way in; after setup it is where the
  ## analyses can be changed without going through the whole guide again.
  output$wf_design_intro <- renderUI({
    if (isTRUE(wf$active)) return(NULL)
    if (!isTRUE(wf$configured)) {
      return(div(
        class = "card wf-design-intro",
        div(class = "card-title", icon("compass-drafting"), " Custom Analysis Design"),
        p(class = "submodule-desc",
          "Set the analysis up once here - the dataset, the comparison, the sex design and which analyses to run - and every page below follows it. The sidebar then tracks what has run, what is ready and what is waiting on something else."),
        actionButton("wf_design_start", tagList(icon("play"), " Set up an analysis"), class = "btn btn-primary")
      ))
    }
    ct <- wf$contrast
    div(
      class = "card wf-design-intro",
      ## The whole editor folds away on a click, so the page is a summary until it is needed again.
      tags$details(
        class = "wf-design-toggle", open = "open",
        tags$summary(
          class = "card-title",
          icon("compass-drafting"), " Custom Analysis Design",
          tags$span(class = "wf-design-summary", sprintf(
            " - %s, %s",
            if (is.null(ct)) "no comparison set" else sprintf("%s vs %s (%s)", ct$comp, ct$ref, ct$col),
            unname(WF_SEX_MODES[[wf$sex_mode %||% "pooled"]] %||% "pooled")))
        ),
        p(class = "submodule-desc", "Add an analysis you left out, or drop one. The sidebar tracker follows whatever is chosen here."),
        if (identical(wf$layer, "transcriptomics")) tagList(preset_picker_ui(), preset_plan_ui()),
        div(
          class = "wf-bar-actions",
          actionButton("wf_design_restart", tagList(icon("rotate-left"), " Start over"), class = "btn btn-default btn-sm"),
          actionButton("wf_design_update", tagList(icon("check"), " Update analyses"), class = "btn btn-primary btn-sm")
        )
      )
    )
  })
  outputOptions(output, "wf_design_intro", suspendWhenHidden = FALSE)

  ## Re-applies the picker without re-running the guide; new analyses appear in the sidebar at once.
  observeEvent(input$wf_design_update, {
    req(wf$configured, identical(wf$layer, "transcriptomics"), !is.null(nav$apply))
    plan <- selection_plan()
    wf$preset <- plan$preset
    wf$selected <- plan$keys
    nav$apply(wf$layer, plan$keys)
    set_msg("ok", sprintf("Design updated: %d analyses.", length(plan$keys)))
  }, ignoreInit = TRUE)

  observeEvent(input$wf_design_restart, show_start_dialog(), ignoreInit = TRUE)

  ## ---- Start: step 1 dialog. ----
  show_start_dialog <- function() {
    showModal(modalDialog(
      title = tagList(icon("play"), " Start an analysis"),
      easyClose = TRUE,
      ## Radio, not one-click buttons: the choice can be changed, and Cancel leaves everything as it
      ## was. Start is the only thing that commits.
      footer = tagList(
        modalButton("Cancel"),
        actionButton("wf_start_confirm", "Start", class = "btn-primary")
      ),
      p(class = "submodule-desc", "Pick a layer, then Start. Transcriptomics and Methylomics are walked step by step; Cross-Omics and Multi-Omics open their own page. You can come back and change this at any time."),
      radioButtons(
        "wf_layer", "Which omics layer?",
        choiceNames = lapply(WF_PICK_LAYERS, function(L) {
          tagList(icon(L$icon), tags$span(class = "wf-pick-name", L$label))
        }),
        choiceValues = lapply(WF_PICK_LAYERS, `[[`, "id"),
        selected = wf$layer %||% "transcriptomics"
      )
    ))
  }

  observeEvent(input$home_start_analysis, show_start_dialog(), ignoreInit = TRUE)
  observeEvent(input$wf_design_start, show_start_dialog(), ignoreInit = TRUE)

  observeEvent(input$wf_start_confirm, {
    lay <- as.character(input$wf_layer %||% "transcriptomics")
    removeModal()
    ## Cross-Omics and Multi-Omics have no guided pipeline (their pages take DEG/DMP tables and
    ## matched multi-block samples, not one matrix with a group column), so Start opens the page.
    if (!lay %in% names(WF_LAYERS)) {
      if (!is.null(nav$section)) nav$section(lay)
      return()
    }
    start_guide(lay)
  }, ignoreInit = TRUE)

  start_guide <- function(layer) {
    queue(list())
    wf$active <- TRUE
    wf$layer <- layer
    wf$done <- integer(0); wf$contrast <- NULL; wf$sex_mode <- NULL; wf$plan <- NULL
    wf$run_status <- "idle"; wf$summary <- character(0); wf$complete <- FALSE; wf$dataset_label <- NULL
    wf$configured <- FALSE
    ## With valid data already loaded, go straight to step 2: opening the Dataset tab in the same flush
    ## would win over the QC tab (tab-select messages are applied after tab inserts).
    ## Starting opens the layer and stops there. Nothing is chosen, nothing is marked done and no
    ## panel is pushed at the user - the design appears only when they open Custom Analysis Design,
    ## which is where step 4 lives.
    if (dataset_valid()) {
      wf$dataset_label <- if (identical(layer, "transcriptomics")) dataset$source else methyl_dataset$source
      wf$step <- 4L
      if (identical(layer, "transcriptomics") && !is.null(nav$section)) nav$section(layer)
      else enter_step(4L)
    } else {
      enter_step(1L)
      load_example()
    }
  }

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

  ## The comparison the design is about to use. The Differential Expression page wins when it is
  ## open, but its inputs only exist while that tab is shown (Shiny suspends outputs in hidden
  ## tabs), so on the design page this falls back to the dataset's own default - otherwise the sex
  ## options would be stuck on Pooled. It is stated on screen and stays changeable.
  effective_contrast <- reactive({
    req(wf$active, wf$layer)
    pc <- page_contrast()
    if (!is.null(pc$col) && nzchar(pc$ref %||% "") && nzchar(pc$comp %||% "")) {
      return(c(pc, list(source = "page")))
    }
    if (!is.null(wf$contrast)) return(c(wf$contrast, list(source = "guide")))
    meta <- layer_meta()
    gc <- wf_group_col(meta)
    if (is.null(gc)) return(list(col = NULL, ref = NULL, comp = NULL, source = "none"))
    dc <- wf_default_contrast(meta[[gc]])
    if (is.null(dc)) return(list(col = NULL, ref = NULL, comp = NULL, source = "none"))
    list(col = gc, ref = dc$ref, comp = dc$comp, source = "dataset")
  })

  contrast_check <- reactive({
    req(wf$active, wf$layer)
    ec <- effective_contrast()
    meta <- layer_meta()
    wf_check_contrast(meta, ec$col, ec$ref, ec$comp, wf_sex_col(meta), wf$layer)
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
      ## Choose analysis: which analyses to run, the comparison (defaulted from the dataset unless
      ## changed on the page) and the sex design, locked together in one step.
      if (identical(wf$layer, "methylomics") && is.null(methyl_dataset$beta)) {
        set_msg("error", "This deployment has only the example dataset's sample sheet, not its live methylation matrix, so the DMP analysis cannot run here.")
        return()
      }
      cc <- contrast_check()
      if (!isTRUE(cc$ok)) { set_msg("error", paste(cc$problems, collapse = " ")); return() }
      ec <- effective_contrast()
      new <- list(col = ec$col, ref = ec$ref, comp = ec$comp)
      if (!identical(new, wf$contrast)) { clear_after(3L); wf$sex_mode <- NULL; wf$plan <- NULL }
      wf$contrast <- new

      mode <- input$wf_sex_mode
      if (is.null(mode) || !isTRUE(cc$modes[[mode]])) { set_msg("error", "Pick one of the available sex analyses first."); return() }
      if (!identical(mode, wf$sex_mode)) clear_after(3L)
      wf$sex_mode <- mode

      meta <- layer_meta()
      wf$plan <- wf_run_plan(wf$layer, mode, wf_sex_col(meta), cc$sex_levels)
      complete_step(3L)
      if (identical(wf$layer, "transcriptomics") && !is.null(nav$apply)) {
        plan <- selection_plan()
        wf$preset <- plan$preset
        wf$selected <- plan$keys
        ## The design is fixed: open exactly what was chosen, then step out of the way. From here
        ## the sidebar tracker and each module's own Run button carry the analysis.
        wf$configured <- TRUE
        wf$active <- FALSE
        wf$complete <- TRUE
        queue(list())
        nav$apply(wf$layer, plan$keys)
      } else {
        enter_step(4L)
      }
    } else if (n == 4L) {
      if (!identical(wf$run_status, "done")) { set_msg("error", "Run the analysis first."); return() }
      complete_step(4L); enter_step(5L)
    } else if (n %in% c(5L, 6L, 7L)) {
      complete_step(n); enter_step(n + 1L)
    } else if (n == 8L) {
      complete_step(8L); wf$complete <- TRUE
      set_msg("ok", "Guided analysis complete. Every page keeps its results; you can keep working on any of them.")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$wf_run, {
    req(wf$active, identical(wf$step, 4L))
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

  observeEvent(input$wf_bar_close, { wf$bar_hidden <- TRUE }, ignoreInit = TRUE)

  observeEvent(input$wf_exit, {
    queue(list())
    wf$active <- FALSE; wf$step <- 1L; wf$done <- integer(0); wf$message <- NULL
    wf$loading <- FALSE; wf$run_status <- "idle"; wf$summary <- character(0); wf$complete <- FALSE
    wf$configured <- FALSE
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

  ## One route for every step-link button: "<module id>|<inner tabset or empty>|<inner tab or empty>".
  observeEvent(input$wf_open_module, {
    req(wf$active)
    parts <- strsplit(as.character(input$wf_open_module), "|", fixed = TRUE)[[1]]
    parts <- c(parts, "", "")[1:3]
    req(nzchar(parts[1]))
    open_validation(list(id = parts[1],
                         tabset = if (nzchar(parts[2])) parts[2] else NULL,
                         tab = if (nzchar(parts[3])) parts[3] else NULL))
  }, ignoreInit = TRUE)

  ## ---- Rendering. ----
  output$wf_home_tracker <- renderUI({
    wf_home_tracker_ui(wf$step, wf$done, wf$active)
  })

  output$wf_bar <- renderUI({
    if (!isTRUE(wf$active)) return(NULL)
    ## Transcriptomics: the setup box belongs to the Custom Analysis Design page and nowhere else,
    ## and stays closed once the user closes it until that page is opened again.
    if (identical(wf$layer, "transcriptomics") &&
        (isTRUE(wf$bar_hidden) || !identical(input$tx_menu, "Custom Analysis Design"))) return(NULL)
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
          p(class = "submodule-desc", sprintf("The %s page below is the existing quality-control page. These counts come from the loaded dataset's sample metadata - see them here before picking an analysis design.", module_title(lay$qc))),
          wf_note(if (chk$ok) "ok" else "error",
                  sprintf("%s samples. Groups (%s): %s. Sex (%s): %s.",
                          wf_fmt(chk$n), chk$group_col %||% "no column", wf_counts_text(chk$group_counts),
                          chk$sex_col %||% "no column", wf_counts_text(chk$sex_counts))),
          lapply(chk$problems, function(t) wf_note("error", t)),
          lapply(chk$notes, function(t) wf_note("info", t)),
          step_links_ui(list(list(id = lay$qc)))
        )
      },
      "3" = if (identical(wf$layer, "transcriptomics") &&
                !identical(input$tx_menu, "Custom Analysis Design")) {
        ## Off the design page the controls would sit on top of every module, so only the way back
        ## to them is shown here.
        tagList(
          p(class = "submodule-desc", "Pick the analyses, the comparison and the sex design on the Custom Analysis Design page."),
          div(class = "wf-bar-links", tags$button(
            type = "button", class = "btn btn-primary btn-sm",
            onclick = "Shiny.setInputValue('tx_nav_go', '__design__', {priority: 'event'})",
            icon("compass-drafting"), " Open Custom Analysis Design"))
        )
      } else {
        ## The comparison resolves from the dataset's own group column by default (effective_contrast()'s
        ## fallback), so this never has to wait on the DE/DMP page's own inputs to mount before showing
        ## something - it reads whatever is in effect right now and lets you override it on that page.
        ec <- effective_contrast()
        cc <- contrast_check()
        avail <- names(cc$modes)[cc$modes]
        desc <- c(
          pooled = "All samples in one model, adjusted for sex.",
          stratified = "The same comparison run separately in each sex.",
          sex_specific = "Tests whether the disease effect differs between the sexes."
        )
        sel <- intersect(c(isolate(input$wf_sex_mode), wf$sex_mode, "pooled"), avail)[1]
        tagList(
          if (identical(wf$layer, "transcriptomics")) preset_picker_ui(),
          p(class = "submodule-desc", sprintf("The comparison below is the dataset's default; change the column, reference group or comparison group yourself on the %s page if you want a different one. Pick the sex design, then press Lock this design to confirm everything - the guide picks nothing for you.", module_title(lay$compare))),
          if (!is.null(ec$col)) wf_note(if (cc$ok) "ok" else "error",
            sprintf("Contrast column: %s. Reference: %s (n = %s). Comparison: %s (n = %s).",
                    ec$col, ec$ref %||% "?", wf_fmt(cc$n_ref), ec$comp %||% "?", wf_fmt(cc$n_comp)))
          else wf_note("error", "No usable group or diagnosis column was found on the loaded dataset."),
          if (!cc$ok && !is.null(ec$col)) lapply(cc$problems, function(t) wf_note("error", t)),
          step_links_ui(list(list(id = lay$compare))),
          if (isTRUE(cc$ok)) div(
            class = "wf-design-layout",
            div(
              class = "wf-design-main",
              radioButtons("wf_sex_mode", "Sex analysis", inline = TRUE,
                           choices = stats::setNames(avail, WF_SEX_MODES[avail]), selected = sel),
              tags$ul(class = "wf-mode-list", lapply(names(desc), function(m) tags$li(
                tags$strong(unname(WF_SEX_MODES[[m]])), " - ", desc[[m]],
                if (!isTRUE(cc$modes[[m]])) tags$em(paste(" Not available:",
                  if (m %in% names(cc$mode_notes)) cc$mode_notes[[m]] else "check the comparison.")))))
            ),
            ## The resulting pipeline, read top to bottom down the right edge.
            if (identical(wf$layer, "transcriptomics")) div(
              class = "wf-design-aside",
              div(class = "wf-aside-title", "Pipeline"),
              preset_plan_ui()
            )
          )
        )
      },
      "4" = {
        ## Run analysis automates the differential run(s) for the chosen sex design. The other
        ## analyses picked in Step 3 have their own settings, so the guide opens them rather than
        ## pressing their run buttons for you.
        others <- setdiff(intersect(wf$selected, TX_STEP_MAP_ORDER),
                          c("dge", "interaction", "diagnostic_external"))
        others <- Filter(function(k) TX_STEP_MAP[[k]]$step == 4L, others)
        tagList(
          p(class = "submodule-desc", "Run analysis presses the page's own run button with these settings, so the results are identical to running it by hand."),
          tags$ul(class = "wf-mode-list", lapply(wf$plan, function(r) tags$li(r$label))),
          if (identical(wf$layer, "transcriptomics") && length(others)) tagList(
            p(class = "submodule-desc", "Also chosen for this step - open each one and press its own Run button:"),
            step_links_ui(lapply(others, function(k) list(id = TX_STEP_MAP[[k]]$id)))
          ),
          if (length(wf$summary)) div(class = "wf-summary", lapply(wf$summary, function(t) wf_note("ok", t)))
        )
      },
      "5" = tagList(
        p(class = "submodule-desc", sprintf("Candidate %s is the existing page below. No candidates is a valid result, not an error; you can still continue.", lay$feature)),
        step_links_ui(lay$explore)
      ),
      "6" = tagList(
        p(class = "submodule-desc", "Feature Selection turns the candidate list into a consensus panel; the Diagnostic Model trains a classifier on that panel and scores it on its own held-out test samples."),
        step_links_ui(lay$model)
      ),
      "7" = tagList(
        p(class = "submodule-desc", "Validation runs on the existing pages below, against an independent cohort, another tissue and another ancestry."),
        step_links_ui(lay$validate)
      ),
      "8" = tagList(
        p(class = "submodule-desc", "The Biomarker Card brings the dataset evidence and every result saved above together for one feature at a time."),
        step_links_ui(list(list(id = lay$card)))
      )
    )
    ## Whatever currently stops the pipeline, computed once and shown unconditionally.
    blockers <- {
      probs <- character(0)
      if (n >= 2L) {
        chk <- tryCatch(sample_check(), error = function(e) NULL)
        if (!is.null(chk) && !isTRUE(chk$ok)) probs <- c(probs, chk$problems)
      }
      if (n >= 3L) {
        cc <- tryCatch(contrast_check(), error = function(e) NULL)
        if (!is.null(cc) && !isTRUE(cc$ok)) probs <- c(probs, cc$problems)
      }
      probs <- unique(probs[nzchar(probs)])
      if (length(probs)) lapply(probs, function(t) wf_note("error", t)) else NULL
    }

    primary <- if (n == 4L && !identical(wf$run_status, "done")) {
      actionButton("wf_run", "Run analysis", icon = icon("play"), class = "btn btn-primary btn-sm",
                   disabled = if (identical(wf$run_status, "running")) "disabled" else NULL)
    } else if (n == WF_N_STEPS && isTRUE(wf$complete)) {
      NULL
    } else {
      last <- n == WF_N_STEPS
      lbl <- if (last) "Finish" else if (n == 3L) "Lock this design" else "Continue"
      ic <- if (last) "flag-checkered" else if (n == 3L) "lock" else "arrow-right"
      actionButton("wf_next", lbl, icon = icon(ic),
                   class = "btn btn-primary btn-sm", disabled = if (is_busy) "disabled" else NULL)
    }
    div(
      class = "card wf-bar",
      div(
        class = "wf-bar-head",
        div(class = "card-title", icon("route"), sprintf(" Start an analysis · %s · Step %d of %d: %s", lay$label, n, WF_N_STEPS, s$label)),
        div(
          class = "wf-bar-head-actions",
          actionLink("wf_bar_close", tagList(icon("chevron-up"), " Hide"), class = "wf-bar-exit",
                     title = "Hide this box. It comes back when you open Custom Analysis Design."),
          actionLink("wf_exit", tagList(icon("xmark"), " Exit guide"), class = "wf-bar-exit")
        )
      ),
      ## One progress indicator per page by default: on Transcriptomics the sidebar carries the eight
      ## step statuses, so the bar keeps the strip folded away behind a summary the user can open
      ## when they want the whole route at a glance. <details> needs no input id and no round-trip,
      ## so opening it never re-renders the bar. Methylomics has no step-grouped sidebar yet, so it
      ## keeps the strip open as before.
      ## Transcriptomics shows the strip in the page header instead (see tx_settings_banner), so the
      ## bar does not repeat it. Methylomics has no such header slot yet.
      if (!identical(wf$layer, "transcriptomics")) wf_bar_steps_ui(n, wf$done),
      ## Anything blocking is stated in red outside the fold - an error the user cannot see is worse
      ## than no error at all.
      blockers,
      if (identical(wf$layer, "transcriptomics") &&
          !identical(input$tx_menu, "Custom Analysis Design")) {
        ## Away from the design page the guide is a controller, not a panel - the explanation is
        ## there when it is asked for and out of the way otherwise.
        tags$details(class = "wf-steps-toggle wf-body-toggle",
                     tags$summary(icon("circle-info"), " What this step is for"),
                     body)
      } else body,
      if (!is.null(msg)) wf_note(msg$kind, msg$text),
      div(
        class = "wf-bar-actions",
        if (n > 1L) actionButton("wf_back", "Back", icon = icon("arrow-left"), class = "btn btn-default btn-sm",
                                 disabled = if (identical(wf$run_status, "running")) "disabled" else NULL),
        if (n == 4L && identical(wf$run_status, "done")) actionButton("wf_run", "Run again", icon = icon("play"), class = "btn btn-default btn-sm"),
        primary
      )
    )
  })

  invisible(wf)
}
