## R/multiomics/cohort_harmonization_helpers.R
## Pure data-processing logic for the "Cohort Harmonization" sub-module
## (mod_multi_overview.R) - data-adaptive: every function here inspects
## whatever is actually in the Active Multi-Omics Dataset (multi_dataset,
## built on the Dataset Workspace tab) rather than assuming any fixed
## modality combination. Nothing here silently invents a match, an
## outcome, or a "beats chance" claim - a fact that can't be established
## from the actual data is reported as "Unmatched"/"Unknown"/"Not
## available", never guessed.
##
## Reuses multi_live_validate_matrix()/multi_live_pca()/
## multi_live_correlation_heatmap_data() etc. (multiomics_live_helpers.R)
## wherever they already answer the question - nothing here duplicates
## that logic.

## ---------------------------------------------------------------------------
## 0. Sample Explorer - a per-sample master table (one row per sample seen
## in any selected modality, "Present"/"Missing" per modality, plus every
## detected metadata column) so the user can browse and search individual
## samples rather than only aggregate counts.
## ---------------------------------------------------------------------------

ch_sample_master_table <- function(id_sets, meta = NULL) {
  id_sets <- Filter(function(x) length(x) > 0, id_sets)
  if (length(id_sets) == 0) return(NULL)
  all_ids <- sort(unique(unlist(id_sets)))
  presence <- lapply(id_sets, function(ids) all_ids %in% ids)
  df <- data.frame(`Sample ID` = all_ids, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(id_sets)) df[[nm]] <- ifelse(presence[[nm]], "Present", "Missing")
  df$`Modalities present` <- Reduce(`+`, presence)
  if (!is.null(meta) && nrow(meta) > 0) {
    common <- intersect(all_ids, rownames(meta))
    for (col in colnames(meta)) {
      vals <- stats::setNames(rep(NA_character_, length(all_ids)), all_ids)
      vals[common] <- as.character(meta[common, col])
      df[[col]] <- unname(vals)
    }
  }
  df
}

## ---------------------------------------------------------------------------
## 1. Modality descriptors - one entry per modality in the active dataset,
## built adaptively depending on whether a raw matrix actually exists
## (Upload/GEO) or only a per-patient availability table does (Preloaded -
## no raw expression/methylation matrix is bundled in this deployment).
## ---------------------------------------------------------------------------

## Best-effort, always-labeled-as-a-guess description of what a matrix's
## values likely are - never used to transform data, only to report
## "Preprocessing status: Unknown" honestly when it can't be told (spec
## section 11).
ch_value_scale <- function(mat) {
  if (is.null(mat) || !is.matrix(mat)) return("Unknown")
  v <- as.numeric(mat)
  v <- v[is.finite(v)]
  if (length(v) < 10) return("Unknown")
  mn <- min(v); mx <- max(v); mu <- mean(v)
  frac_int <- mean(abs(v - round(v)) < 1e-6)
  if (mn >= -1e-6 && mx <= 1 + 1e-6) return("Beta-values (0-1 scale, likely)")
  if (mn < -0.01 && abs(mu) < 0.5 && mx < 20) return("Standardized/z-scored (likely)")
  if (mn >= -1e-6 && frac_int > 0.95 && mx > 20) return("Raw counts (likely)")
  if (mn >= -1e-6 && mx < 30) return("Log-transformed (likely)")
  "Unknown"
}

## Upload/GEO: real matrices exist in multi_dataset$layers.
ch_modality_descriptors_live <- function(multi_dataset) {
  layers <- multi_dataset$layers %||% list()
  if (length(layers) == 0) return(list())
  out <- lapply(names(layers), function(nm) {
    mat <- layers[[nm]]
    lm <- multi_dataset$layer_meta[[nm]] %||% list()
    list(
      label = nm, n_samples = nrow(mat), n_features = ncol(mat),
      sample_ids = rownames(mat), has_raw_matrix = TRUE,
      omics_type = lm$omics_type %||% "other", validation = lm$validation,
      value_scale = ch_value_scale(mat), processing = lm$processing %||% "Not processed",
      status = lm$status
    )
  })
  stats::setNames(out, names(layers))
}

## Preloaded: only a per-patient availability table + QC summaries exist -
## no raw sample x feature matrix is bundled in this deployment. Report
## exactly what's determinable (sample counts, feature counts) and nothing
## more; anything needing a raw matrix (ID-format detail, PCA, correlation,
## model evaluation) reports "Insufficient information" downstream.
ch_modality_descriptors_preloaded <- function() {
  matching <- multi_read_registry_table("Patient sample matching (all 80 patients)")
  if (!matching$ok || !"patient_id" %in% colnames(matching$df)) return(list())
  df <- matching$df
  rna <- multi_read_registry_table("RNA-seq QC summary")
  meth <- multi_read_registry_table("Methylation QC summary")

  out <- list()
  if ("RNA_available_PBMC" %in% colnames(df)) {
    rna_ok <- df$RNA_available_PBMC %in% c(TRUE, "TRUE", "Yes", "yes", 1)
    rna_pbmc <- if (rna$ok && all(c("cell_type", "n_genes_retained") %in% colnames(rna$df))) rna$df[rna$df$cell_type == "PBMC", , drop = FALSE] else NULL
    out[["Transcriptomics"]] <- list(
      label = "Transcriptomics", n_samples = sum(rna_ok),
      n_features = if (!is.null(rna_pbmc) && nrow(rna_pbmc) > 0) rna_pbmc$n_genes_retained[1] else NA,
      sample_ids = as.character(df$patient_id[rna_ok]), has_raw_matrix = FALSE,
      omics_type = "rnaseq", validation = NULL,
      value_scale = "Unknown (no raw matrix bundled for the preloaded dataset)",
      processing = "Normalized (pipeline QC)", status = NULL
    )
  }
  if ("methylation_available" %in% colnames(df)) {
    meth_ok <- df$methylation_available %in% c(TRUE, "TRUE", "Yes", "yes", 1)
    out[["Methylomics"]] <- list(
      label = "Methylomics", n_samples = sum(meth_ok),
      n_features = if (meth$ok && "n_probes_retained" %in% colnames(meth$df)) meth$df$n_probes_retained[1] else NA,
      sample_ids = as.character(df$patient_id[meth_ok]), has_raw_matrix = FALSE,
      omics_type = "methylation", validation = NULL,
      value_scale = "Unknown (no raw matrix bundled for the preloaded dataset)",
      processing = "Normalized (pipeline QC)", status = NULL
    )
  }
  out
}

ch_modality_descriptors <- function(multi_dataset) {
  if (is.null(multi_dataset) || !isTRUE(multi_dataset$active)) return(list())
  if (identical(multi_dataset$source, "preloaded")) return(ch_modality_descriptors_preloaded())
  ch_modality_descriptors_live(multi_dataset)
}

## ---------------------------------------------------------------------------
## 2. Pairwise sample overlap - the NxN matrix (spec section 5B), built
## purely from the actual detected sample identifiers per modality, never
## from row order.
## ---------------------------------------------------------------------------

ch_pairwise_overlap_matrix <- function(sample_id_lists) {
  sample_id_lists <- Filter(Negate(is.null), sample_id_lists)
  nm <- names(sample_id_lists)
  n <- length(nm)
  if (n == 0) return(NULL)
  m <- matrix(0L, n, n, dimnames = list(nm, nm))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    m[i, j] <- length(intersect(sample_id_lists[[nm[i]]], sample_id_lists[[nm[j]]]))
  }
  m
}

## ---------------------------------------------------------------------------
## 3. Sample-ID harmonization - safe normalization only (trim + case-fold),
## never applied to the underlying data, never fuzzy-matched. Every ID gets
## an explicit status (spec section 6): Exact match / Normalized match /
## Duplicate / Unmatched / Ambiguous / Invalid.
## ---------------------------------------------------------------------------

ch_id_harmonization_table <- function(sample_id_lists) {
  sample_id_lists <- Filter(Negate(is.null), sample_id_lists)
  if (length(sample_id_lists) == 0) return(NULL)
  norm_id <- function(x) trimws(tolower(as.character(x)))

  rows <- do.call(rbind, lapply(names(sample_id_lists), function(mod) {
    ids <- as.character(sample_id_lists[[mod]])
    dup <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
    data.frame(Modality = mod, Original = ids, Normalized = norm_id(ids), DupInModality = dup, stringsAsFactors = FALSE)
  }))

  status <- character(nrow(rows))
  reason <- character(nrow(rows))
  by_norm <- split(seq_len(nrow(rows)), rows$Normalized)
  for (nrm in names(by_norm)) {
    idx <- by_norm[[nrm]]
    if (!nzchar(nrm)) {
      status[idx] <- "Invalid"; reason[idx] <- "Empty or missing identifier."
      next
    }
    mods_here <- unique(rows$Modality[idx])
    origs_here <- unique(rows$Original[idx])
    for (i in idx) {
      other_mods <- setdiff(mods_here, rows$Modality[i])
      if (rows$DupInModality[i]) {
        if (length(mods_here) >= 2) {
          status[i] <- "Ambiguous"
          reason[i] <- sprintf("Repeated within %s and also appears in %s - correspondence cannot be determined.", rows$Modality[i], paste(other_mods, collapse = ", "))
        } else {
          status[i] <- "Duplicate"
          reason[i] <- sprintf("Repeated within %s.", rows$Modality[i])
        }
      } else if (length(mods_here) < 2) {
        status[i] <- "Unmatched"
        reason[i] <- sprintf("Only present in %s.", rows$Modality[i])
      } else if (length(origs_here) == 1) {
        status[i] <- "Exact match"
        reason[i] <- sprintf("Identical ID in %s.", paste(mods_here, collapse = ", "))
      } else {
        status[i] <- "Normalized match"
        reason[i] <- sprintf("Matched after trimming whitespace/case across %s.", paste(mods_here, collapse = ", "))
      }
    }
  }
  rows$Status <- status
  rows$Reason <- reason
  rows[, c("Modality", "Original", "Normalized", "Status", "Reason")]
}

## ---------------------------------------------------------------------------
## 4. Candidate phenotype/batch columns - name/shape heuristics only, never
## an auto-selected single answer (spec section 16: "do not fabricate an
## outcome"). The user picks from these in a filter.
## ---------------------------------------------------------------------------

## Column names that read as a sample/patient identifier rather than a
## biological or technical grouping variable - shared by ch_detect_candidate_
## columns() and ch_classify_metadata_columns() so both use the same rule.
CH_ID_LIKE_NAME_REGEX <- "^(sample|id|patient|subject)([_.]?id)?$"

## Value-shape classification of one metadata column: "identifier" (near-
## all-unique - not a useful grouping/coloring variable), "continuous"
## (numeric with high cardinality), else "categorical". Report-only, mirrors
## the cardinality rule mi_outcome_summary() (multiomics_integration_live_
## helpers.R) already uses for its own outcome-type check - unified here so
## the two don't quietly disagree on a borderline column.
ch_classify_column <- function(v) {
  vv <- v[!is.na(v) & nzchar(trimws(as.character(v)))]
  n_total <- length(vv)
  if (n_total == 0) return("categorical")
  nu <- length(unique(vv))
  if (nu > 10 && (nu / n_total) > 0.9) return("identifier")
  if (is.numeric(v) && nu > min(10, max(3, floor(n_total / 3)))) return("continuous")
  "categorical"
}

## Classifies every metadata column and suggests a default categorical
## grouping variable (spec: "sensible default, always user-overridable") -
## first a name-keyword match (phenotype/response/group/...), else the first
## reasonably balanced categorical column (no single level over 90% of
## samples - a looser bar than mi_outcome_summary()'s 0.7 "imbalanced" flag
## on purpose, since this only picks a *starting* selection, not an
## eligibility gate).
ch_classify_metadata_columns <- function(meta) {
  if (is.null(meta) || ncol(meta) == 0) return(list(table = NULL, suggested_default = NULL))
  id_like <- grepl(CH_ID_LIKE_NAME_REGEX, colnames(meta), ignore.case = TRUE)
  rows <- lapply(colnames(meta), function(cn) {
    v <- meta[[cn]]
    type <- if (id_like[match(cn, colnames(meta))]) "identifier" else ch_classify_column(v)
    nu <- length(unique(v[!is.na(v) & nzchar(trimws(as.character(v)))]))
    data.frame(column = cn, type = type, n_unique = nu, stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)
  categorical_cols <- tbl$column[tbl$type == "categorical"]
  keyword_hit <- categorical_cols[grepl("phenotype|response|group|condition|treatment|outcome|status", categorical_cols, ignore.case = TRUE)]
  suggested <- if (length(keyword_hit) > 0) keyword_hit[1] else {
    balanced <- Filter(function(cn) {
      tab <- table(as.character(meta[[cn]]))
      length(tab) > 0 && (max(tab) / sum(tab)) <= 0.9
    }, categorical_cols)
    if (length(balanced) > 0) balanced[1] else if (length(categorical_cols) > 0) categorical_cols[1] else NULL
  }
  list(table = tbl, suggested_default = suggested)
}

ch_detect_candidate_columns <- function(meta, kind = c("phenotype", "batch")) {
  kind <- match.arg(kind)
  if (is.null(meta) || ncol(meta) == 0) return(character(0))
  id_like <- grepl(CH_ID_LIKE_NAME_REGEX, colnames(meta), ignore.case = TRUE)
  cols <- colnames(meta)[!id_like]
  if (identical(kind, "batch")) {
    return(cols[grepl("batch|cohort|study|platform|site|processing|run|plate", cols, ignore.case = TRUE)])
  }
  cls <- ch_classify_metadata_columns(meta)
  intersect(cols, cls$table$column[cls$table$type == "categorical"])
}

## ---------------------------------------------------------------------------
## 4b. Matched-sample summary sentence + tri-state status (spec: "42
## transcriptomics / 39 methylomics / 35 matched / 83.3% overlap", and a
## clear Matched/Partially matched/Unmatched call so an unmatched cohort is
## never presented as if it were paired multi-omics data).
## ---------------------------------------------------------------------------

ch_matched_sample_summary <- function(id_sets) {
  id_sets <- Filter(Negate(is.null), id_sets)
  if (length(id_sets) == 0) return(list(per_modality = integer(0), n_matched = 0L, n_union = 0L, pct_overlap = 0, sentence = "No modalities available.", status = "Unmatched"))
  per_modality <- vapply(id_sets, length, integer(1))
  n_matched <- length(Reduce(intersect, id_sets))
  n_union <- length(Reduce(union, id_sets))
  pct_overlap <- if (n_union > 0) 100 * n_matched / n_union else 0
  status <- if (length(id_sets) < 2) "Single modality"
    else if (n_matched == 0) "Unmatched"
    else if (n_matched == min(per_modality)) "Matched"
    else "Partially matched"
  sentence <- paste0(
    paste(sprintf("%s: %s", names(id_sets), format(per_modality, big.mark = ",")), collapse = " / "),
    if (length(id_sets) >= 2) sprintf(" / %s matched / %.1f%% overlap", format(n_matched, big.mark = ","), pct_overlap) else ""
  )
  list(per_modality = per_modality, n_matched = n_matched, n_union = n_union, pct_overlap = pct_overlap, sentence = sentence, status = status)
}

## ---------------------------------------------------------------------------
## 5. Analysis cells - every feasible modality combination (spec sections
## 9-10): never forces complete-case-only integration, never hard-coded to
## a fixed set of modalities. Capped at all-subsets for <= 4 modalities
## (2^N-1 stops being compact beyond that) - singles + pairs + the full
## overlap otherwise, with the omission explicitly stated.
## ---------------------------------------------------------------------------

ch_analysis_cells <- function(id_sets, pheno_available = FALSE, max_full_subsets = 4,
                               min_integration = 3, min_prediction = 6) {
  mods <- names(Filter(Negate(is.null), id_sets))
  n <- length(mods)
  if (n == 0) return(list(cells = list(), omitted_note = NULL))

  all_subsets <- function(x) {
    k <- length(x)
    idx <- unlist(lapply(seq_len(k), function(j) utils::combn(k, j, simplify = FALSE)), recursive = FALSE)
    lapply(idx, function(ix) x[ix])
  }

  omitted_note <- NULL
  if (n <= max_full_subsets) {
    combos <- all_subsets(mods)
  } else {
    singles <- as.list(mods)
    pairs <- utils::combn(mods, 2, simplify = FALSE)
    combos <- c(singles, pairs, list(mods))
    total_possible <- 2^n - 1
    omitted_note <- sprintf("%d modalities selected - showing single modalities, pairs, and the full combination only (%d of %d possible combinations omitted for brevity).",
                             n, total_possible - length(combos), total_possible)
  }

  cells <- lapply(combos, function(cm) {
    ids <- id_sets[cm]
    matched <- Reduce(intersect, ids)
    n_matched <- length(matched)
    methods <- character(0)
    if (length(cm) >= 2 && n_matched >= min_integration) methods <- c(methods, "Unsupervised integration")
    if (n_matched >= min_prediction && isTRUE(pheno_available)) methods <- c(methods, "Supervised prediction")
    if (length(methods) == 0) methods <- "None feasible"
    list(modalities = cm, label = paste(cm, collapse = " + "), n_matched = n_matched, matched_ids = matched, methods = methods)
  })
  list(cells = cells, omitted_note = omitted_note)
}

## ---------------------------------------------------------------------------
## 6. Integration readiness - Ready / Limited / Not suitable (spec section
## 8), applied per analysis cell.
## ---------------------------------------------------------------------------

ch_integration_readiness <- function(cell, min_ready = 10, min_limited = 3) {
  if (length(cell$modalities) < 2) {
    return(list(level = "single", label = "Single modality", reason = "Integration requires at least two modalities."))
  }
  n <- cell$n_matched
  if (n < min_limited) {
    return(list(level = "not_suitable", label = "Not suitable",
                reason = sprintf("Only %d matched sample(s) across %s.", n, cell$label)))
  }
  if (n < min_ready) {
    return(list(level = "limited", label = "Limited",
                reason = sprintf("%d matched samples - integration is possible, but sample size may limit interpretation.", n)))
  }
  list(level = "ready", label = "Ready", reason = sprintf("%d matched samples across %s.", n, cell$label))
}

## ---------------------------------------------------------------------------
## 7. Held-out binary-outcome evaluation (spec sections 16-18) - stratified
## k-fold CV, elastic-net logistic regression (glmnet), feature selection
## and scaling fit inside the training fold only. Per-omics models plus one
## early-fusion (concatenated top-variance features) model; majority-class
## baseline reported alongside. Refuses outright (never silently
## downgrades) when guardrails aren't met.
## ---------------------------------------------------------------------------

## One fold's train/test AUC contribution for a single feature matrix -
## variance-ranked feature selection and z-scoring computed from the
## training rows only, never touching the held-out rows.
ch_fold_predict_view <- function(X, y, train_idx, test_idx, max_features) {
  Xtr <- X[train_idx, , drop = FALSE]; Xte <- X[test_idx, , drop = FALSE]
  ytr <- y[train_idx]
  vars <- apply(Xtr, 2, stats::var, na.rm = TRUE); vars[is.na(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(min(max_features, sum(vars > 0)))]
  if (length(keep) < 2 || length(unique(ytr)) < 2) return(rep(0.5, length(test_idx)))
  Xtr_k <- Xtr[, keep, drop = FALSE]; Xte_k <- Xte[, keep, drop = FALSE]
  mu <- colMeans(Xtr_k, na.rm = TRUE); sdv <- apply(Xtr_k, 2, stats::sd, na.rm = TRUE); sdv[sdv == 0 | is.na(sdv)] <- 1
  Xtr_s <- scale(Xtr_k, center = mu, scale = sdv); Xte_s <- scale(Xte_k, center = mu, scale = sdv)
  Xtr_s[!is.finite(Xtr_s)] <- 0; Xte_s[!is.finite(Xte_s)] <- 0
  fit <- tryCatch(glmnet::cv.glmnet(Xtr_s, ytr, family = "binomial", alpha = 0.5, nfolds = max(3, min(5, min(table(ytr))))), error = function(e) NULL)
  if (is.null(fit)) return(rep(0.5, length(test_idx)))
  pred <- tryCatch(as.numeric(stats::predict(fit, newx = Xte_s, s = "lambda.min", type = "response")), error = function(e) NULL)
  if (is.null(pred)) rep(0.5, length(test_idx)) else pred
}

ch_evaluate_binary_outcome <- function(mat_list, y, k_folds = 5, seed = 1, max_features_per_view = 200) {
  mat_list <- Filter(Negate(is.null), mat_list)
  if (length(mat_list) < 1) return(list(ok = FALSE, error = "No omics layers available for evaluation."))

  common <- Reduce(intersect, lapply(mat_list, rownames))
  common <- intersect(common, names(y))
  if (any(duplicated(common))) return(list(ok = FALSE, error = "Duplicate sample IDs detected across the matched samples - resolve before evaluation."))
  if (length(common) < 10) return(list(ok = FALSE, error = sprintf("Only %d matched samples with a valid outcome - at least 10 are needed for reliable held-out evaluation.", length(common))))

  y <- droplevels(factor(y[common]))
  if (nlevels(y) != 2) return(list(ok = FALSE, error = sprintf("Outcome has %d level(s) - only binary outcomes are supported in this delivery.", nlevels(y))))
  tab <- table(y)
  if (min(tab) < 3) return(list(ok = FALSE, error = sprintf("Smallest class has %d sample(s) - at least 3 per class are needed for reliable held-out evaluation.", min(tab))))
  k_folds <- max(3, min(k_folds, min(tab)))

  mats <- lapply(mat_list, function(m) m[common, , drop = FALSE])

  set.seed(seed)
  folds <- caret::createFolds(y, k = k_folds, list = TRUE, returnTrain = FALSE)

  run_view <- function(X) {
    oof <- rep(NA_real_, length(y))
    for (f in folds) {
      test_idx <- f; train_idx <- setdiff(seq_along(y), f)
      oof[test_idx] <- ch_fold_predict_view(X, y, train_idx, test_idx, max_features_per_view)
    }
    oof
  }
  run_fused <- function() {
    oof <- rep(NA_real_, length(y))
    for (f in folds) {
      test_idx <- f; train_idx <- setdiff(seq_along(y), f)
      parts_tr <- list(); parts_te <- list()
      for (nm in names(mats)) {
        Xtr <- mats[[nm]][train_idx, , drop = FALSE]; Xte <- mats[[nm]][test_idx, , drop = FALSE]
        vars <- apply(Xtr, 2, stats::var, na.rm = TRUE); vars[is.na(vars)] <- 0
        keep <- order(vars, decreasing = TRUE)[seq_len(min(50, sum(vars > 0)))]
        if (length(keep) < 1) next
        Xtr_k <- Xtr[, keep, drop = FALSE]; Xte_k <- Xte[, keep, drop = FALSE]
        mu <- colMeans(Xtr_k, na.rm = TRUE); sdv <- apply(Xtr_k, 2, stats::sd, na.rm = TRUE); sdv[sdv == 0 | is.na(sdv)] <- 1
        parts_tr[[nm]] <- scale(Xtr_k, center = mu, scale = sdv); parts_te[[nm]] <- scale(Xte_k, center = mu, scale = sdv)
      }
      if (length(parts_tr) < 2) { oof[test_idx] <- 0.5; next }
      Xtr_f <- do.call(cbind, parts_tr); Xte_f <- do.call(cbind, parts_te)
      Xtr_f[!is.finite(Xtr_f)] <- 0; Xte_f[!is.finite(Xte_f)] <- 0
      ytr <- y[train_idx]
      fit <- tryCatch(glmnet::cv.glmnet(Xtr_f, ytr, family = "binomial", alpha = 0.5, nfolds = max(3, min(5, min(table(ytr))))), error = function(e) NULL)
      oof[test_idx] <- if (is.null(fit)) 0.5 else {
        p <- tryCatch(as.numeric(stats::predict(fit, newx = Xte_f, s = "lambda.min", type = "response")), error = function(e) NULL)
        if (is.null(p)) 0.5 else p
      }
    }
    oof
  }

  safe_roc <- function(oof) tryCatch(pROC::roc(y, oof, quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)

  per_view_auc <- list(); per_view_roc <- list()
  for (nm in names(mats)) {
    oof <- run_view(mats[[nm]])
    r <- safe_roc(oof)
    per_view_roc[[nm]] <- r
    per_view_auc[[nm]] <- if (!is.null(r)) as.numeric(pROC::auc(r)) else NA_real_
  }

  fused_auc <- NA_real_; fused_ci <- c(NA_real_, NA_real_); fused_roc <- NULL
  if (length(mats) >= 2) {
    fused_roc <- safe_roc(run_fused())
    if (!is.null(fused_roc)) {
      fused_auc <- as.numeric(pROC::auc(fused_roc))
      ci <- tryCatch(as.numeric(pROC::ci.auc(fused_roc, quiet = TRUE)), error = function(e) c(NA, NA, NA))
      fused_ci <- ci[c(1, 3)]
    }
  }

  best_single <- if (length(per_view_auc) > 0 && any(!is.na(unlist(per_view_auc)))) names(which.max(unlist(per_view_auc))) else NA_character_
  best_single_auc <- if (!is.na(best_single)) per_view_auc[[best_single]] else NA_real_

  vs_single_p <- NA_real_
  if (!is.null(fused_roc) && !is.na(best_single) && !identical(best_single, "")) {
    vs_single_p <- tryCatch(pROC::roc.test(fused_roc, per_view_roc[[best_single]], method = "delong", quiet = TRUE)$p.value, error = function(e) NA_real_)
  }

  list(
    ok = TRUE, n = length(y), k_folds = k_folds, levels = levels(y),
    per_view_auc = per_view_auc, best_single = best_single, best_single_auc = best_single_auc,
    fused_auc = fused_auc, fused_ci = fused_ci, vs_single_p = vs_single_p,
    majority_baseline = as.numeric(max(prop.table(tab)))
  )
}
