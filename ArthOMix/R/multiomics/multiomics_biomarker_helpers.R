## R/multiomics/multiomics_biomarker_helpers.R
## Pure data-processing logic for Biomarker Discovery - a live, data-adaptive
## supervised DIABLO (mixOmics::block.splsda) engine scoped specifically to
## this submodule's primary use case: Transcriptomics + Methylomics ->
## supervised multi-omics feature selection -> predictive performance ->
## feature stability -> interpretable biomarker signature. Reuses
## mi_validate_dataset()/mi_diablo_eligibility()/mi_diablo_run()/
## mi_diablo_performance_summary()/mi_diablo_selected_features_df()/
## mi_diablo_stability_df()/mi_outcome_summary()/mi_preloaded_cell_dataset()
## (multiomics_integration_helpers.R) rather than reimplementing the
## verified DIABLO call shapes documented there - Biomarker Discovery differs
## from Multi-omics Integration only in scope (exactly two blocks, always
## labeled Transcriptomics/Methylomics), UI organization, and the
## evidence-based stability labeling / cross-validated ROC / signature
## reporting added below. DIABLO is never fit on anything but matched
## samples with a real categorical outcome (mi_diablo_eligibility()'s own
## rules) - no outcome is ever fabricated.
##
## `mb_*` prefix throughout, matching the `mi_*`/`ch_*` per-engine
## conventions already used elsewhere in this module.

## ---------------------------------------------------------------------------
## 1. Data source - reuses mi_preloaded_cell_dataset() (preloaded) and the
## shared multi_dataset (Active Multi-Omics Dataset, built on the Dataset
## Workspace tab) for uploaded data. No separate upload widget is duplicated
## in this submodule - Dataset Workspace already implements file-type/
## orientation/sample-ID/duplicate/missing-value detection for any number of
## omics layers, exactly what spec section 5 asks for, and this submodule is
## not the place to re-implement it.
## ---------------------------------------------------------------------------

## From whichever N-omics-generic dataset is selected, pick exactly the two
## blocks the user has assigned as the "Transcriptomics"/"Methylomics" roles
## - never assumed from block name or position. Re-labels them to the fixed
## names "Transcriptomics"/"Methylomics" so every downstream helper/plot can
## address them by role rather than by whatever the source dataset called them.
mb_select_blocks <- function(layers, transcript_block, methyl_block) {
  if (is.null(layers) || is.null(transcript_block) || is.null(methyl_block)) return(NULL)
  if (identical(transcript_block, methyl_block)) return(NULL)
  if (!all(c(transcript_block, methyl_block) %in% names(layers))) return(NULL)
  stats::setNames(list(layers[[transcript_block]], layers[[methyl_block]]), c("Transcriptomics", "Methylomics"))
}

## Friendlier DISPLAY label only, for the "Classes to include" dropdown - the
## VALUE submitted/filtered on is always the outcome column's real, unaltered
## string (e.g. the preloaded cohort's own "resp"/"non" codes), never
## replaced or coerced. Recognizes only a short, explicit list of common
## responder/non-responder abbreviations; anything not on this list is shown
## exactly as it appears in the data - this is a label lookup, never a guess,
## and it never invents "Responder"/"Non-responder" as a *selectable* option
## when the data doesn't actually contain a matching class (doing so would
## silently filter every sample out for that empty class).
MB_CLASS_LABEL_SYNONYMS <- c(
  "resp" = "Responder", "responder" = "Responder", "response" = "Responder",
  "non" = "Non-responder", "nonresponder" = "Non-responder", "non_responder" = "Non-responder",
  "non-responder" = "Non-responder", "nonresponse" = "Non-responder", "non_response" = "Non-responder"
)

mb_friendly_class_label <- function(x) {
  hit <- MB_CLASS_LABEL_SYNONYMS[tolower(trimws(as.character(x)))]
  ifelse(is.na(hit), as.character(x), unname(hit))
}

## ---------------------------------------------------------------------------
## 2. Unsupervised feature-count cap (spec section 9's "feature dimensions
## are computationally reasonable" / section 15's "too many features:
## provide an appropriate feature-filtering option rather than crashing").
## Ranks by variance alone - never reads the outcome - so this is not a
## source of outcome leakage; it is applied once, before any CV split, the
## same way a lab would restrict to "the N most variable probes" before
## doing anything supervised.
## ---------------------------------------------------------------------------

mb_variance_prefilter <- function(mat, max_features) {
  if (is.null(mat) || is.null(max_features) || ncol(mat) <= max_features) return(mat)
  v <- apply(mat, 2, stats::var, na.rm = TRUE)
  v[is.na(v)] <- -Inf
  keep <- sort(order(v, decreasing = TRUE)[seq_len(max_features)])
  mat[, keep, drop = FALSE]
}

## Default cap: generous enough to almost never bind for a curated/selected
## panel (like the preloaded cell's own already-selected features), but
## bounded for genome-wide/array-wide uploads so automatic keepX tuning
## below stays within the "may take several minutes" promise.
mb_default_max_features <- function(n_features) as.integer(min(n_features, 5000L))

## ---------------------------------------------------------------------------
## 3. Compact "Data Check" table (spec section 6) - one row per requirement,
## built only from mi_validate_dataset()/mi_outcome_summary()/
## mi_diablo_eligibility()'s own already-computed fields, never a second,
## independently-computed check that could disagree with what DIABLO itself
## will see.
## ---------------------------------------------------------------------------

MB_MISSING_ACCEPTABLE_PCT <- 5

mb_data_check_table <- function(validation, outcome_summary, eligibility) {
  row <- function(check, status) data.frame(Check = check, Status = status, stringsAsFactors = FALSE)
  if (is.null(validation) || !isTRUE(validation$ok)) {
    return(rbind(
      row("Transcriptomics", "Not available"), row("Methylomics", "Not available"),
      row("Matched samples", "-"), row("Outcome", "-"), row("Missing values", "-"),
      row("DIABLO eligibility", if (!is.null(eligibility) && !isTRUE(eligibility$ok)) eligibility$reason else "Not ready")
    ))
  }
  block_status <- function(nm) {
    b <- validation$per_block[[nm]]
    if (is.null(b) || !isTRUE(b$ok)) "Not available" else sprintf("Available (%s samples x %s features)", format(b$n_samples, big.mark = ","), format(b$n_features, big.mark = ","))
  }
  miss <- vapply(validation$per_block, function(b) if (isTRUE(b$ok)) b$pct_missing else NA_real_, numeric(1))
  missing_pct <- suppressWarnings(max(miss, na.rm = TRUE))
  missing_status <- if (!is.finite(missing_pct)) "-" else if (missing_pct <= MB_MISSING_ACCEPTABLE_PCT) sprintf("Acceptable (%.1f%% missing)", missing_pct) else sprintf("High (%.1f%% missing) - resolve before running", missing_pct)
  outcome_status <- if (is.null(outcome_summary)) "No outcome selected" else if (identical(outcome_summary$type, "categorical")) {
    if (is.na(outcome_summary$n_classes) || outcome_summary$n_classes < 2) "Fewer than 2 usable classes" else sprintf("%d classes (n = %d matched)%s", outcome_summary$n_classes, outcome_summary$n, if (isTRUE(outcome_summary$imbalanced)) " - imbalanced" else "")
  } else sprintf("Continuous (%d values) - not usable for DIABLO as-is", outcome_summary$n)
  rbind(
    row("Transcriptomics", block_status("Transcriptomics")),
    row("Methylomics", block_status("Methylomics")),
    row("Matched samples", if (isTRUE(validation$reliable_matching)) format(validation$n_shared, big.mark = ",") else sprintf("%d (not enough to proceed)", validation$n_shared)),
    row("Outcome", outcome_status),
    row("Missing values", missing_status),
    row("DIABLO eligibility", if (isTRUE(eligibility$ok)) "Ready" else sprintf("Not ready - %s", eligibility$reason %||% "requirements not met"))
  )
}

## ---------------------------------------------------------------------------
## 4. Evidence-based feature stability (spec section 11 - replaces "live-
## adjustable confidence relabeling" entirely). Categories are cut from
## perf()'s own real per-repeat selection-frequency table
## (mi_diablo_stability_df(), multiomics_integration_helpers.R) against
## fixed, disclosed thresholds - never a user-adjustable slider presented as
## a statistical result.
## ---------------------------------------------------------------------------

## NOTE on the 0.8 cutoff vs. SNF Clustering's own stability metric: SNF
## Clustering's SFC_STABILITY_THRESHOLDS (snf_clustering_helpers.R) uses a
## different "stable" cutoff (0.75). This is deliberate, not an oversight -
## the two are not the same statistic and are not meant to share one bar.
## MB_STABILITY_THRESHOLDS classifies a per-feature *selection frequency*:
## the raw fraction of CV repeats in which mixOmics::block.splsda's own
## keepX selection kept that feature (mi_diablo_stability_df()), i.e. a
## plain (uncorrected-for-chance) proportion, consistent with the stability-
## selection literature's typical "reliable" cutoffs in the 0.6-0.9 range
## (Meinshausen & Buhlmann 2010). SFC_STABILITY_THRESHOLDS instead classifies
## a mean Adjusted Rand Index between two full sample-partitions (subsampled
## SNF clustering vs. the full-cohort reference clustering) - ARI is
## chance-corrected (0 expected under random labelings, 1 at perfect
## agreement), so a given numeric value reflects a different amount of
## "real" agreement than an uncorrected selection-frequency proportion at
## the same number; conventional ARI cutoffs for "excellent"/"high"
## clustering agreement (e.g. Hubert & Arabie 1985; Ben-Hur et al. 2002) sit
## around 0.75, not 0.8. Forcing these two constants to match would paper
## over that the underlying statistics live on different effective scales.
MB_STABILITY_THRESHOLDS <- list(stable = 0.8, moderate = 0.5)

mb_stability_category <- function(freq) {
  cut(freq, breaks = c(-Inf, MB_STABILITY_THRESHOLDS$moderate, MB_STABILITY_THRESHOLDS$stable, Inf),
      labels = c("Low stability", "Moderately stable", "Stable"), right = FALSE)
}

## Full Signature table (spec sections 12-13.A) - one row per selected
## feature x component, joined with its cross-validated selection frequency
## and evidence-based stability category. Ranked by |loading| within
## omics/component, matching mi_diablo_selected_features_df()'s own order.
mb_signature_table <- function(diablo_res) {
  sel <- mi_diablo_selected_features_df(diablo_res$fit)
  if (is.null(sel)) return(NULL)
  blocks <- unique(sel$block)
  stab_rows <- lapply(blocks, function(b) {
    do.call(rbind, lapply(unique(sel$component[sel$block == b]), function(cmp) {
      s <- mi_diablo_stability_df(diablo_res, block = b, comp = cmp)
      if (is.null(s) || nrow(s) == 0) return(NULL)
      data.frame(block = b, component = cmp, feature = s$feature, selection_frequency = s$stability, stringsAsFactors = FALSE)
    }))
  })
  stab_df <- do.call(rbind, Filter(Negate(is.null), stab_rows))

  out <- sel
  if (!is.null(stab_df) && nrow(stab_df) > 0) {
    out <- merge(out, stab_df, by = c("block", "component", "feature"), all.x = TRUE)
  } else {
    out$selection_frequency <- NA_real_
  }
  out$stability_category <- as.character(mb_stability_category(out$selection_frequency))
  out$stability_category[is.na(out$selection_frequency)] <- "Not available (needs >1 CV repeat)"
  names(out)[names(out) == "block"] <- "omics"
  out <- out[order(out$omics, out$component, -abs(out$loading)), , drop = FALSE]
  out$rank_within_block <- stats::ave(seq_len(nrow(out)), out$omics, out$component, FUN = seq_along)
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------------------
## 5. Cross-validated component-score correlation between the two blocks
## (spec section 13.D "Component correlations") - the real
## fit$variates[[block]][,1] values from the fitted model, not a fabricated
## relationship.
## ---------------------------------------------------------------------------

mb_component_correlation <- function(fit) {
  blocks <- setdiff(names(fit$variates %||% list()), "Y")
  if (length(blocks) != 2) return(NULL)
  ids <- rownames(fit$variates[[blocks[1]]])
  ids <- intersect(ids, rownames(fit$variates[[blocks[2]]]))
  if (length(ids) < 3) return(NULL)
  r <- suppressWarnings(stats::cor(fit$variates[[blocks[1]]][ids, 1], fit$variates[[blocks[2]]][ids, 1]))
  list(block_a = blocks[1], block_b = blocks[2], r = r, n = length(ids))
}

## ---------------------------------------------------------------------------
## 6. Pooled out-of-fold ROC curve (spec sections 13.B, 13.E.8) - a real
## held-out evaluation using the SAME stratified fold count already resolved
## for the main DIABLO run (never a fresh, more favorable split), pooling
## predicted scores across folds exactly like
## cohort_harmonization_helpers.R's own ch_evaluate_binary_outcome()
## out-of-fold idiom. Restricted to binary outcomes ("only when the outcome
## is suitable for ROC analysis") - returns NULL for other class counts
## rather than picking one class arbitrarily. Each fold's own block.splsda()
## is fit on that fold's training samples only, so feature selection
## (keepX) never sees a held-out sample's data - no leakage.
## ---------------------------------------------------------------------------

mb_cv_roc <- function(X, Y, params, seed = 1) {
  if (nlevels(Y) != 2) return(NULL)
  n <- length(Y)
  folds_k <- max(2, min(params$folds %||% 5, min(table(Y))))
  set.seed(seed)
  fold_id <- tryCatch(caret::createFolds(Y, k = folds_k, list = FALSE), error = function(e) NULL)
  if (is.null(fold_id)) return(NULL)
  pos_class <- levels(Y)[2]; neg_class <- levels(Y)[1]
  oof <- stats::setNames(rep(NA_real_, n), names(Y))
  for (f in seq_len(folds_k)) {
    te <- which(fold_id == f); tr <- setdiff(seq_len(n), te)
    ytr <- droplevels(Y[tr])
    if (nlevels(ytr) < 2 || length(te) == 0) next
    Xtr <- lapply(X, function(m) m[tr, , drop = FALSE])
    Xte <- lapply(X, function(m) m[te, , drop = FALSE])
    fit_f <- tryCatch(
      mixOmics::block.splsda(X = Xtr, Y = ytr, ncomp = params$ncomp, keepX = params$keepX, design = params$design, near.zero.var = TRUE, scale = isTRUE(params$scale %||% TRUE)),
      error = function(e) NULL
    )
    if (is.null(fit_f)) next
    pr <- tryCatch(stats::predict(fit_f, newdata = Xte, dist = params$distance), error = function(e) NULL)
    if (is.null(pr) || is.null(pr$WeightedPredict)) next
    wp_arr <- pr$WeightedPredict[, , params$ncomp, drop = FALSE]
    wp <- matrix(as.numeric(wp_arr), nrow = dim(wp_arr)[1], dimnames = list(dimnames(wp_arr)[[1]], dimnames(wp_arr)[[2]]))
    if (!pos_class %in% colnames(wp)) next
    oof[rownames(wp)] <- wp[, pos_class]
  }
  ok <- !is.na(oof)
  if (sum(ok) < 4) return(NULL)
  roc_obj <- tryCatch(pROC::roc(Y[ok], oof[ok], levels = levels(Y), direction = "<", quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NULL)
  list(roc = roc_obj, auc = as.numeric(pROC::auc(roc_obj)), n_used = sum(ok), n_total = n, folds = folds_k, pos_class = pos_class, neg_class = neg_class)
}

## ---------------------------------------------------------------------------
## 7. Reproducibility record (spec section 23) - every value actually used,
## never a placeholder; software versions read from the real installed
## packages (multi_package_versions(), multiomics_helpers.R).
## ---------------------------------------------------------------------------

mb_reproducibility_table <- function(diablo_res, dataset_label, preprocessing_note, seed) {
  p <- diablo_res$params
  keepx_txt <- function(b) paste(p$keepX[[b]], collapse = ",")
  data.frame(
    Parameter = c(
      "Data source", "Samples analyzed (matched)", "Outcome classes", "Number of components (ncomp)",
      "keepX (Transcriptomics)", "keepX (Methylomics)", "Design (Transcriptomics <-> Methylomics)",
      "Prediction distance", "Validation method", "CV folds", "CV repeats", "Random seed", "Preprocessing"
    ),
    Value = c(
      dataset_label %||% "-", p$n_samples, paste(p$classes, collapse = ", "), p$ncomp,
      keepx_txt("Transcriptomics"), keepx_txt("Methylomics"), sprintf("%.2f", p$design["Transcriptomics", "Methylomics"]),
      p$distance, p$validation_method, p$folds, p$nrepeat, seed, preprocessing_note %||% "None beyond what the source dataset already carries."
    ),
    stringsAsFactors = FALSE
  )
}

mb_matched_sample_table <- function(validation, outcome, sample_ids) {
  if (is.null(validation) || length(sample_ids) == 0) return(NULL)
  data.frame(
    sample_id = sample_ids,
    outcome = if (!is.null(outcome)) as.character(outcome[sample_ids]) else NA_character_,
    stringsAsFactors = FALSE
  )
}
