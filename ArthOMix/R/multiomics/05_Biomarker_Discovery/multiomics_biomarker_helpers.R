## R/multiomics/05_Biomarker_Discovery/multiomics_biomarker_helpers.R
## Pure data-processing logic for Biomarker Discovery (live DIABLO block.splsda).

mb_select_blocks <- function(layers, transcript_block, methyl_block) {
  if (is.null(layers) || is.null(transcript_block) || is.null(methyl_block)) return(NULL)
  if (identical(transcript_block, methyl_block)) return(NULL)
  if (!all(c(transcript_block, methyl_block) %in% names(layers))) return(NULL)
  stats::setNames(list(layers[[transcript_block]], layers[[methyl_block]]), c("Transcriptomics", "Methylomics"))
}

MB_CLASS_LABEL_SYNONYMS <- c(
  "resp" = "Responder", "responder" = "Responder", "response" = "Responder",
  "non" = "Non-responder", "nonresponder" = "Non-responder", "non_responder" = "Non-responder",
  "non-responder" = "Non-responder", "nonresponse" = "Non-responder", "non_response" = "Non-responder"
)

mb_friendly_class_label <- function(x) {
  hit <- MB_CLASS_LABEL_SYNONYMS[tolower(trimws(as.character(x)))]
  ifelse(is.na(hit), as.character(x), unname(hit))
}

mb_variance_prefilter <- function(mat, max_features) {
  if (is.null(mat) || is.null(max_features) || ncol(mat) <= max_features) return(mat)
  v <- apply(mat, 2, stats::var, na.rm = TRUE)
  v[is.na(v)] <- -Inf
  keep <- sort(order(v, decreasing = TRUE)[seq_len(max_features)])
  mat[, keep, drop = FALSE]
}

mb_default_max_features <- function(n_features) as.integer(min(n_features, 5000L))

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

MB_STABILITY_THRESHOLDS <- list(stable = 0.8, moderate = 0.5)

mb_stability_category <- function(freq) {
  cut(freq, breaks = c(-Inf, MB_STABILITY_THRESHOLDS$moderate, MB_STABILITY_THRESHOLDS$stable, Inf),
      labels = c("Low stability", "Moderately stable", "Stable"), right = FALSE)
}

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

mb_component_correlation <- function(fit) {
  blocks <- setdiff(names(fit$variates %||% list()), "Y")
  if (length(blocks) != 2) return(NULL)
  ids <- rownames(fit$variates[[blocks[1]]])
  ids <- intersect(ids, rownames(fit$variates[[blocks[2]]]))
  if (length(ids) < 3) return(NULL)
  r <- suppressWarnings(stats::cor(fit$variates[[blocks[1]]][ids, 1], fit$variates[[blocks[2]]][ids, 1]))
  list(block_a = blocks[1], block_b = blocks[2], r = r, n = length(ids))
}

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

## Independent DIABLO evaluation: sealed hold-out split + external-cohort transfer.

MB_HOLDOUT_MIN_PER_CLASS <- 3L

## Stratified hold-out split, sealed before any modelling, tuning or CV.
mb_holdout_split <- function(outcome, sample_ids, frac, seed = 1) {
  frac <- as.numeric(frac %||% 0)
  if (is.na(frac) || frac <= 0) return(list(ok = TRUE, train_ids = sample_ids, test_ids = character(0), frac = 0))
  Y <- droplevels(factor(outcome[sample_ids]))
  set.seed(seed)
  idx <- tryCatch(caret::createDataPartition(Y, p = 1 - frac, list = FALSE)[, 1], error = function(e) NULL)
  if (is.null(idx)) return(list(ok = FALSE, error = "Could not create a stratified hold-out split for this outcome."))
  train_ids <- sample_ids[idx]; test_ids <- setdiff(sample_ids, train_ids)
  tab_tr <- table(droplevels(factor(outcome[train_ids]))); tab_te <- table(factor(outcome[test_ids], levels = levels(Y)))
  if (length(tab_tr) < 2 || min(tab_tr) < MB_HOLDOUT_MIN_PER_CLASS || length(test_ids) < 2 * MB_HOLDOUT_MIN_PER_CLASS || min(tab_te) < 1) {
    return(list(ok = FALSE, error = sprintf(
      "A %.0f%% hold-out leaves too few samples (training %s; hold-out %s). Each class needs >=%d training samples, and the hold-out needs >=%d samples with both classes. Lower the fraction or set it to 0.",
      100 * frac, paste(sprintf("%s = %d", names(tab_tr), as.integer(tab_tr)), collapse = ", "),
      paste(sprintf("%s = %d", names(tab_te), as.integer(tab_te)), collapse = ", "), MB_HOLDOUT_MIN_PER_CLASS, 2 * MB_HOLDOUT_MIN_PER_CLASS)))
  }
  list(ok = TRUE, train_ids = train_ids, test_ids = test_ids, frac = frac)
}

## Scores new samples with a fitted block.splsda: positive-class score, AUROC (DeLong CI), BER, confusion table.
mb_score_new_samples <- function(fit, X_new, Y_new, ncomp, distance = "max.dist") {
  Y_new <- factor(Y_new, levels = levels(fit$Y))
  dist_use <- if (is.null(distance) || identical(distance, "automatic") || identical(distance, "all")) "max.dist" else distance
  pr <- tryCatch(stats::predict(fit, newdata = X_new, dist = dist_use), error = function(e) e)
  if (inherits(pr, "error")) return(list(ok = FALSE, error = paste("Prediction on the new samples failed:", conditionMessage(pr))))
  pos_class <- levels(fit$Y)[2]; neg_class <- levels(fit$Y)[1]
  wp <- tryCatch({
    arr <- pr$WeightedPredict[, , ncomp, drop = FALSE]
    matrix(as.numeric(arr), nrow = dim(arr)[1], dimnames = list(dimnames(arr)[[1]], dimnames(arr)[[2]]))
  }, error = function(e) NULL)
  if (is.null(wp) || !pos_class %in% colnames(wp)) return(list(ok = FALSE, error = "The model did not return a weighted prediction score."))
  score <- stats::setNames(wp[, pos_class], rownames(wp))
  vote <- tryCatch(pr$WeightedVote[[dist_use]][, ncomp], error = function(e) NULL)
  if (is.null(vote)) vote <- ifelse(score >= 0.5, pos_class, neg_class)
  vote <- factor(as.character(vote), levels = levels(fit$Y))
  names(vote) <- rownames(wp)
  y <- Y_new[rownames(wp)]
  ok <- !is.na(y) & !is.na(score)
  if (sum(ok) < 4 || length(unique(y[ok])) < 2) return(list(ok = FALSE, error = "Too few scored samples with both classes to evaluate."))
  roc_obj <- tryCatch(pROC::roc(y[ok], score[ok], levels = levels(fit$Y), direction = "<", quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(list(ok = FALSE, error = "The ROC curve could not be computed."))
  auc <- as.numeric(pROC::auc(roc_obj))
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), error = function(e) c(NA_real_, auc, NA_real_))
  conf <- table(truth = y[ok], predicted = vote[ok])
  per_class_err <- vapply(levels(fit$Y), function(l) { n_l <- sum(y[ok] == l); if (n_l == 0) NA_real_ else mean(vote[ok][y[ok] == l] != l) }, numeric(1))
  list(
    ok = TRUE, n = sum(ok), n_per_class = table(y[ok]), auc = auc, ci_lo = ci[1], ci_hi = ci[3],
    ## excludes_chance is direction-BLIND - kept for backward compatibility only. Anything
    ## that decides pass/fail or color MUST use auroc_call (multi_auroc_call() in
    ## multiomics_helpers.R): a significantly BELOW-chance CI also "excludes chance" but is
    ## the opposite of a validation success.
    excludes_chance = isTRUE(ci[1] > 0.5) || isTRUE(ci[3] < 0.5),
    auroc_call = multi_auroc_call(auc, ci[1], ci[3]),
    ber = mean(per_class_err, na.rm = TRUE), per_class_error = per_class_err, confusion = conf,
    scores = data.frame(sample_id = names(score)[ok], outcome = as.character(y[ok]), score = as.numeric(score[ok]), predicted = as.character(vote[ok]), stringsAsFactors = FALSE),
    pos_class = pos_class, neg_class = neg_class, distance = dist_use, roc = roc_obj
  )
}

## Sealed hold-out evaluation of the model fitted on the training samples.
mb_holdout_evaluate <- function(diablo_res, layers, outcome, test_ids) {
  if (!isTRUE(diablo_res$ok) || length(test_ids) == 0) return(NULL)
  X_new <- lapply(diablo_res$params$blocks, function(b) layers[[b]][test_ids, colnames(diablo_res$fit$X[[b]]), drop = FALSE])
  names(X_new) <- diablo_res$params$blocks
  res <- mb_score_new_samples(diablo_res$fit, X_new, outcome[test_ids], diablo_res$params$ncomp, diablo_res$params$distance)
  if (isTRUE(res$ok)) res$test_ids <- test_ids
  res
}

## Maps Ensembl IDs <-> gene symbols (org.Hs.eg.db) for panel transfer.
mb_gene_id_type <- function(ids) {
  ids <- as.character(ids); ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return("unknown")
  if (mean(grepl("^ENSG[0-9]+", ids)) > 0.5) "ensembl" else if (mean(grepl("^cg[0-9]+$", ids)) > 0.5) "cpg" else "symbol"
}

mb_ensembl_to_symbol <- function(ids) {
  base <- sub("\\..*$", "", ids)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) return(stats::setNames(rep(NA_character_, length(ids)), ids))
  sym <- tryCatch(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(base), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")), error = function(e) NULL)
  if (is.null(sym)) return(stats::setNames(rep(NA_character_, length(ids)), ids))
  stats::setNames(unname(sym[base]), ids)
}

## Renames `mat` columns to match the target ID type (Ensembl<->symbol); drops unmapped/duplicate features.
mb_harmonise_gene_columns <- function(mat, target_ids) {
  src_type <- mb_gene_id_type(colnames(mat)); tgt_type <- mb_gene_id_type(target_ids)
  if (src_type == tgt_type || src_type == "cpg" || tgt_type == "cpg") return(list(mat = mat, mapped = FALSE))
  if (src_type == "ensembl" && tgt_type == "symbol") {
    map <- mb_ensembl_to_symbol(colnames(mat))
    keep <- !is.na(map) & !duplicated(map)
    m <- mat[, keep, drop = FALSE]; colnames(m) <- unname(map[keep])
    return(list(mat = m, mapped = TRUE, note = sprintf("%d of %d Ensembl gene IDs mapped to unique symbols.", sum(keep), ncol(mat))))
  }
  if (src_type == "symbol" && tgt_type == "ensembl") {
    map <- mb_ensembl_to_symbol(target_ids)
    rev <- stats::setNames(names(map)[!is.na(map)], unname(map[!is.na(map)]))
    rev <- rev[!duplicated(names(rev))]
    keep <- colnames(mat) %in% names(rev)
    m <- mat[, keep, drop = FALSE]; colnames(m) <- unname(rev[colnames(m)])
    return(list(mat = m, mapped = TRUE, note = sprintf("%d of %d gene symbols mapped to the panel's Ensembl IDs.", sum(keep), ncol(mat))))
  }
  list(mat = mat, mapped = FALSE)
}

MB_TRANSFER_MIN_COVERAGE <- 0.5

## Transfers the DIABLO panel: refit on discovery training samples (shared panel features only), then score the external cohort (each side standardised separately).
mb_panel_transfer <- function(diablo_res, train_layers, train_ids, train_outcome,
                              external_layers, external_outcome, block_roles = NULL, min_coverage = MB_TRANSFER_MIN_COVERAGE) {
  if (!isTRUE(diablo_res$ok)) return(list(ok = FALSE, error = "Run the discovery analysis first."))
  blocks <- diablo_res$params$blocks
  sel <- mi_diablo_selected_features_df(diablo_res$fit)
  if (is.null(sel) || nrow(sel) == 0) return(list(ok = FALSE, error = "The discovery model selected no features."))
  ext_ids <- Reduce(intersect, lapply(external_layers, rownames))
  ext_ids <- ext_ids[ext_ids %in% names(external_outcome)[!is.na(external_outcome)]]
  if (length(ext_ids) < 6) return(list(ok = FALSE, error = "Fewer than 6 external samples have both layers and an outcome."))

  coverage <- list(); Xtr <- list(); Xext <- list(); notes <- character(0)
  for (b in blocks) {
    ext_b <- external_layers[[b]]
    if (is.null(ext_b)) return(list(ok = FALSE, error = sprintf("The external cohort has no layer assigned to the \"%s\" block.", b)))
    panel_b <- unique(sel$feature[sel$block == b])
    harm <- mb_harmonise_gene_columns(ext_b, panel_b)
    if (isTRUE(harm$mapped) && !is.null(harm$note)) notes <- c(notes, sprintf("%s: %s", b, harm$note))
    ext_b <- harm$mat
    present <- panel_b[panel_b %in% colnames(ext_b)]
    coverage[[b]] <- list(n_panel = length(panel_b), n_present = length(present), missing = setdiff(panel_b, present))
    if (length(present) < 2) return(list(ok = FALSE, error = sprintf("Only %d of the %d panel features of the \"%s\" block are measured in the external cohort.", length(present), length(panel_b), b),
                                         coverage = coverage))
    Xtr[[b]] <- train_layers[[b]][train_ids, present, drop = FALSE]
    m_ext <- ext_b[ext_ids, present, drop = FALSE]
    storage.mode(m_ext) <- "double"
    m_ext <- scale(m_ext)
    m_ext[is.na(m_ext)] <- 0
    Xext[[b]] <- m_ext
  }
  overall_cov <- sum(vapply(coverage, function(c) c$n_present, numeric(1))) / sum(vapply(coverage, function(c) c$n_panel, numeric(1)))
  if (overall_cov < min_coverage) {
    return(list(ok = FALSE, coverage = coverage, overall_coverage = overall_cov, error = sprintf(
      "Only %.0f%% of the panel features are measured in the external cohort (minimum %.0f%%), so the panel can't be evaluated there. Lower the minimum coverage only for an exploratory check.", 100 * overall_cov, 100 * min_coverage)))
  }
  Ytr <- droplevels(factor(train_outcome[train_ids]))
  ncomp <- min(diablo_res$params$ncomp, min(vapply(Xtr, ncol, integer(1))))
  keepX <- lapply(Xtr, function(m) rep(ncol(m), ncomp))
  design <- diablo_res$params$design
  fit <- tryCatch(mixOmics::block.splsda(X = Xtr, Y = Ytr, ncomp = ncomp, keepX = keepX, design = design, near.zero.var = TRUE, scale = TRUE), error = function(e) e)
  if (inherits(fit, "error")) return(list(ok = FALSE, coverage = coverage, error = paste("Panel refit failed:", conditionMessage(fit))))
  Yext <- factor(external_outcome[ext_ids], levels = levels(Ytr))
  res <- mb_score_new_samples(fit, Xext, Yext, ncomp, diablo_res$params$distance)
  if (!isTRUE(res$ok)) { res$coverage <- coverage; return(res) }
  res$coverage <- coverage; res$overall_coverage <- overall_cov; res$notes <- notes
  res$n_train <- length(train_ids); res$ncomp <- ncomp; res$panel_fit <- fit
  res
}

## Reads an external layer file into a samples x features matrix, auto-detecting orientation.
mb_read_external_layer <- function(path, filename = NULL) {
  df <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE, nrows = 200)), error = function(e) NULL)
  total_rows <- if (!is.null(df)) tryCatch(nrow(data.table::fread(path, select = 1L, showProgress = FALSE)), error = function(e) nrow(df)) else NULL
  orient <- multi_live_detect_orientation(df, total_rows = total_rows)$suggested %||% "samples_rows"
  multi_live_read_matrix(path, orientation = orient, filename = filename)
}

## Maps an external outcome column to discovery classes via pos_levels/neg_levels; everything else becomes NA.
mb_map_external_outcome <- function(meta, outcome_col, pos_levels, neg_levels, pos_class, neg_class) {
  if (is.null(meta) || !outcome_col %in% colnames(meta)) return(NULL)
  v <- as.character(meta[[outcome_col]])
  out <- ifelse(v %in% pos_levels, pos_class, ifelse(v %in% neg_levels, neg_class, NA_character_))
  stats::setNames(out, rownames(meta))
}
