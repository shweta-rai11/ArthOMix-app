## R/multiomics/snf_clustering_helpers.R
## Pure data-processing logic for the "SNF Clustering" submodule
## (mod_multi_stratification.R) - a live, data-adaptive unsupervised patient
## stratification workflow (Similarity Network Fusion, Wang et al. 2014).
##
## This file deliberately does NOT reimplement SNF itself: the affinity
## construction / fusion / spectral clustering / eigengap estimation /
## network concordance / post-hoc outcome evaluation / ARI primitives
## already exist, tested, in multiomics_integration_live_helpers.R (the
## Multi-omics Integration submodule's own live DIABLO/SNF engine) as
## mi_snf_affinity()/mi_snf_run()/mi_snf_concordance()/mi_snf_posthoc_outcome()/
## mi_ari()/mi_snf_feasible_k_range()/MI_SNF_ALPHA_RANGE, and are called here
## directly (every R/ file is sourced into one shared environment - see
## 0_load_omics_modules.R). This file adds only what that engine does not
## already provide for a dedicated stratification workflow: single-omics
## fallback, preprocessing (transform/missing/filter) ahead of clustering,
## clinical-variable detection + association testing (categorical/continuous/
## survival), resampling-based stability, parameter sensitivity, and
## per-feature cluster-association ranking.
##
## Every function is fail-soft (list(ok, ..., error)); nothing here ever
## silently drops samples/features or fabricates a result.

## ---------------------------------------------------------------------------
## 1. Data-source adapters (spec section 25) - converge on the same shape
## mi_preloaded_cell_dataset()/mi_dataset() already use: list(ok, layers,
## sample_meta, label, provenance).
## ---------------------------------------------------------------------------

## Preloaded RA anti-TNF cohort: reuses mi_preloaded_cell_dataset() (the
## ONLY per-sample matrices bundled for this cohort - see that function's own
## header) and additionally merges the pipeline's own richer per-drug
## clinical table (response/age/RF/anti-CCP/smoking - Table_SNFjoint_cluster_
## assignments_*), matched by patient ID, for the four drug x sex cells that
## have one. The OLD precomputed `snf_cluster` column is dropped explicitly -
## it is a different (precomputed, full-cohort) SNF run and must never be
## confused with this module's own live-computed clusters. Drug-pooled cells
## (drug = NA) have no per-drug clinical table to merge and fall back to the
## base outcome-only metadata.
sfc_preloaded_dataset <- function(cell_key) {
  base <- mi_preloaded_cell_dataset(cell_key)
  if (!isTRUE(base$ok)) return(base)
  cell <- multi_cell_by_key(cell_key)
  if (is.null(cell) || is.na(cell$drug)) return(base)
  clin <- multi_read_registry_table(sprintf("SNF patient clusters - %s", cell$drug))
  if (!isTRUE(clin$ok) || !"patient_id" %in% colnames(clin$df)) return(base)
  df <- multi_filter_cell(clin$df, sex = cell$sex)
  if (nrow(df) == 0) return(base)
  rownames(df) <- df$patient_id
  keep_cols <- setdiff(colnames(df), c("patient_id", "snf_cluster", "sex", "drug"))
  ids <- rownames(base$layers[[1]])
  common <- intersect(ids, rownames(df))
  if (length(common) == 0) return(base)
  meta <- if (!is.null(base$sample_meta)) base$sample_meta else data.frame(row.names = ids)
  for (col in keep_cols) {
    meta[[col]] <- NA
    meta[common, col] <- df[common, col]
  }
  base$sample_meta <- meta
  base$provenance <- paste0(
    base$provenance,
    sprintf(" Clinical variables (%s) merged from the pipeline's own per-patient table for %s, matched by patient ID (%d of %d patients matched).",
            paste(keep_cols, collapse = ", "), cell$drug, length(common), length(ids))
  )
  base
}

## Active Multi-Omics Dataset (Dataset Workspace tab) - identical convergence
## shape to mi_dataset()'s "active" branch in mod_multi_integration.R; kept
## as its own small function (rather than inlined per-server) so the Data
## subtab's reactive stays a one-line dispatch.
sfc_active_dataset <- function(multi_dataset) {
  if (is.null(multi_dataset) || !isTRUE(multi_dataset$active) || length(multi_dataset$layers %||% list()) < 1) {
    return(list(ok = FALSE, error = "No Active Multi-Omics Dataset yet - build one on the Dataset Workspace tab, or switch to \"Preloaded RA anti-TNF cohort\" above."))
  }
  list(
    ok = TRUE, layers = multi_dataset$layers, sample_meta = multi_dataset$sample_meta,
    layer_meta = multi_dataset$layer_meta,
    label = sprintf("Active Multi-Omics Dataset (%s)", paste(names(multi_dataset$layers), collapse = " + ")),
    provenance = sprintf("Active Multi-Omics Dataset from the Dataset Workspace tab (source: %s).", multi_dataset$source %||% "unknown")
  )
}

## ---------------------------------------------------------------------------
## 2. Validation / matching (spec sections 6-8) - delegates to the existing,
## tested mi_validate_dataset() for >=2 blocks; adds the single-block case
## (spec Case C, "Single-Omics Clustering") that mi_validate_dataset() does
## not handle (it requires >=2 blocks to define an "overlap"). For one block
## there is nothing to match across - every one of that block's own samples
## is, trivially, the matched set.
## ---------------------------------------------------------------------------

sfc_validate_dataset <- function(layers, sample_meta = NULL, outcome_col = NULL) {
  layers <- Filter(Negate(is.null), layers %||% list())
  if (length(layers) == 0) return(list(ok = FALSE, error = "No omics blocks are available.", n_blocks = 0))
  if (length(layers) == 1) {
    nm <- names(layers)
    ids <- rownames(layers[[1]])
    per_block <- stats::setNames(list(multi_live_validate_matrix(layers[[1]], layer_label = nm)), nm)
    outcome_summary <- if (!is.null(sample_meta) && !is.null(outcome_col) && outcome_col %in% colnames(sample_meta)) {
      mi_outcome_summary(sample_meta, outcome_col, ids)
    } else NULL
    return(list(
      ok = TRUE, n_blocks = 1L, block_labels = nm, per_block = per_block,
      overlap = list(ok = TRUE, per_layer = stats::setNames(length(ids), nm), n_shared = length(ids), shared_ids = ids, layer_only = stats::setNames(list(character(0)), nm)),
      n_shared = length(ids), shared_ids = ids,
      reliable_matching = length(ids) >= MI_MIN_MATCHED_SAMPLES,
      mismatch_message = if (length(ids) < MI_MIN_MATCHED_SAMPLES) sprintf("Only %d sample(s) in this single omics block - at least %d are required.", length(ids), MI_MIN_MATCHED_SAMPLES) else NULL,
      outcome = outcome_summary
    ))
  }
  mi_validate_dataset(layers, sample_meta, outcome_col)
}

## Eligibility (spec section 26's "honest failure" cases) - single vs.
## multi-block resolves the analysis "mode" (never silently upgraded from
## Single-Omics Clustering to a multi-omics claim, and never run at all with
## unresolved missing values).
sfc_eligibility <- function(validation) {
  if (is.null(validation) || !isTRUE(validation$ok) || validation$n_blocks < 1) {
    return(list(ok = FALSE, reason = "No compatible omics blocks were detected.", mode = NULL))
  }
  if (!isTRUE(validation$reliable_matching)) {
    return(list(ok = FALSE, reason = validation$mismatch_message %||% MI_SAMPLE_MISMATCH_MESSAGE, mode = NULL))
  }
  n_missing <- vapply(validation$per_block, function(v) if (isTRUE(v$ok)) v$n_missing else NA_integer_, numeric(1))
  if (any(!is.na(n_missing) & n_missing > 0)) {
    return(list(ok = FALSE, reason = "Missing values remain in one or more selected blocks. SNF cannot proceed until they are handled (Data tab, Preprocessing).", mode = NULL))
  }
  list(ok = TRUE, reason = NULL, mode = if (validation$n_blocks == 1) "single_omics" else "multi_omics_snf")
}

## ---------------------------------------------------------------------------
## 3. Data type detection + transformation choices (spec sections 9-10) -
## when a layer's own declared omics_type is known (Active Multi-Omics
## Dataset layers carry this in layer_meta from the Dataset Workspace's
## upload flow), the existing, richer MULTI_LIVE_NORM_CHOICES/
## multi_live_normalize() (multiomics_live_helpers.R) is used as-is. Only
## when omics_type is unknown (e.g. a preloaded DIABLO-panel block) is a
## value-shape heuristic used instead, and only to offer transforms that
## cannot be invalid for the shape actually observed - never a fixed list
## regardless of data.
## ---------------------------------------------------------------------------

SFC_TRANSFORM_LABELS <- list(
  count = c("none" = "None (already normalized)", "log2" = "Log2(x + 1)"),
  proportion_0_1 = c("none" = "None (beta values / proportions, 0-1)", "mvalue" = "M-value transform (logit)"),
  binary = c("none" = "None (binary data - no transform is valid)"),
  continuous = c("none" = "None (already normalized / continuous)", "autoscale" = "Autoscale (z-score)"),
  unknown = c("none" = "None")
)

sfc_detect_data_type <- function(mat) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(list(type = "unknown", note = "No finite values to inspect."))
  n_unique <- length(unique(vals))
  if (n_unique <= 2) return(list(type = "binary", note = sprintf("Only %d distinct value(s) detected.", n_unique)))
  all_int <- all(abs(vals - round(vals)) < 1e-8)
  in_unit <- all(vals >= -1e-6 & vals <= 1 + 1e-6)
  if (in_unit) return(list(type = "proportion_0_1", note = "All values lie within [0,1] - consistent with beta values / proportions."))
  if (all_int && all(vals >= 0)) return(list(type = "count", note = "Non-negative integers - consistent with raw counts."))
  list(type = "continuous", note = "Continuous numeric values (already normalized, or an arbitrary scale).")
}

sfc_transform_choices <- function(mat, omics_type = NULL) {
  if (!is.null(omics_type) && omics_type %in% names(MULTI_LIVE_NORM_CHOICES)) {
    label <- names(MULTI_LIVE_OMICS_TYPES)[MULTI_LIVE_OMICS_TYPES == omics_type]
    return(list(type = omics_type, note = sprintf("Declared omics type: %s.", if (length(label) > 0) label[1] else omics_type), choices = MULTI_LIVE_NORM_CHOICES[[omics_type]]))
  }
  det <- sfc_detect_data_type(mat)
  list(type = det$type, note = det$note, choices = SFC_TRANSFORM_LABELS[[det$type]] %||% SFC_TRANSFORM_LABELS$unknown)
}

## Chains the existing missing-value / normalization / feature-filtering
## primitives (multiomics_live_helpers.R) in a fixed, explicit order: resolve
## missingness first (SNF hard-requires complete data), then transform, then
## filter. Never silently replaces missing values - `missing_method = "none"`
## with any remaining NAs is a hard stop, reported plainly (spec section 26).
sfc_preprocess_block <- function(mat, transform = "none", missing_method = "none", missing_threshold = 50,
                                  filter_criterion = "none", filter_top_n = NULL) {
  log <- character(0)
  out <- mat
  if (any(is.na(out))) {
    if (identical(missing_method, "none")) {
      return(list(ok = FALSE, error = "This block has missing values and no missing-value handling was selected. Choose an option under Preprocessing before running SNF.", log = log))
    }
    mh <- multi_live_handle_missing(out, method = missing_method, max_sample_missing_pct = missing_threshold, max_feature_missing_pct = missing_threshold)
    if (!isTRUE(mh$ok)) return(list(ok = FALSE, error = mh$error, log = log))
    out <- mh$mat
    log <- c(log, sprintf("Missing values: %s (threshold %.0f%%) - dropped %d sample(s), %d feature(s); %d value(s) remain missing.",
                           missing_method, missing_threshold, length(mh$dropped_samples), length(mh$dropped_features), mh$n_remaining_na))
  }
  if (!identical(transform, "none") && !is.null(transform)) {
    tr <- multi_live_normalize(out, omics_type = NULL, method = transform)
    if (isTRUE(tr$ok)) { out <- tr$mat; log <- c(log, sprintf("Transformation applied: %s.", transform)) }
  }
  if (!identical(filter_criterion, "none") && !is.null(filter_criterion) && !is.null(filter_top_n) && filter_top_n < ncol(out)) {
    ff <- multi_live_filter_features(out, criterion = filter_criterion, keep_top_n = filter_top_n)
    if (isTRUE(ff$ok)) { out <- ff$mat; log <- c(log, sprintf("Feature filter (%s): kept top %d of %d features.", filter_criterion, ff$n_after, ff$n_before)) }
  }
  if (any(is.na(out))) return(list(ok = FALSE, error = "Missing values remain after preprocessing - SNF cannot proceed until they are fully resolved.", mat = out, log = log))
  if (nrow(out) < MI_MIN_MATCHED_SAMPLES || ncol(out) < 2) return(list(ok = FALSE, error = "Too few samples or features remain after preprocessing.", mat = out, log = log))
  list(ok = TRUE, mat = out, error = NULL, log = log)
}

## ---------------------------------------------------------------------------
## 4. Clustering (spec sections 4, 12-15) - wraps mi_snf_run() for >=2 blocks
## unchanged, and adds a single-omics fallback (spec Case C) that reuses the
## exact same affinity/spectral-clustering/eigengap primitives (never
## reimplemented) with no fusion step. Both paths return the identical shape
## mi_snf_run() already does (ok, Wall, W, clusters, cluster_estimate,
## params) so every downstream plot/concordance/post-hoc/ARI helper works
## unmodified for either mode - only `params$mode` distinguishes them, so the
## UI can label results "Single-Omics Clustering" and never call a
## single-block result "multi-omics integration".
## ---------------------------------------------------------------------------

sfc_snf_run <- function(layers, params = list()) {
  layers <- Filter(Negate(is.null), layers %||% list())
  if (length(layers) == 0) return(list(ok = FALSE, error = "No omics blocks selected."))
  if (length(layers) >= 2) {
    res <- mi_snf_run(layers, params)
    if (isTRUE(res$ok)) res$params$mode <- "multi_omics_snf"
    return(res)
  }

  standardize <- isTRUE(params$standardize %||% TRUE)
  n <- nrow(layers[[1]])
  k_range <- mi_snf_feasible_k_range(n)
  k_manual <- identical(params$k_mode %||% "automatic", "manual")
  alpha_manual <- identical(params$alpha_mode %||% "automatic", "manual")
  k <- if (k_manual) max(2, min(as.integer(params$k %||% k_range$default), n - 1)) else k_range$default
  alpha <- if (alpha_manual) (params$alpha %||% MI_SNF_ALPHA_RANGE$default) else MI_SNF_ALPHA_RANGE$default

  W <- tryCatch(mi_snf_affinity(layers[[1]], k, alpha, standardize), error = function(e) NULL)
  if (is.null(W)) return(list(ok = FALSE, error = "Could not build the similarity network for this block."))
  rownames(W) <- colnames(W) <- rownames(layers[[1]])
  Wall <- stats::setNames(list(W), names(layers))

  max_k_clusters <- min(6, n - 1)
  est <- tryCatch(SNFtool::estimateNumberOfClustersGivenGraph(W, NUMC = 2:max_k_clusters), error = function(e) NULL)
  cluster_mode <- params$cluster_mode %||% "automatic"
  n_clusters <- if (identical(cluster_mode, "manual") && !is.null(params$n_clusters)) {
    max(2, min(as.integer(params$n_clusters), max_k_clusters))
  } else if (!is.null(est)) est[["Eigen-gap best"]] else 2

  clusters <- tryCatch(SNFtool::spectralClustering(W, K = n_clusters), error = function(e) NULL)
  if (is.null(clusters)) return(list(ok = FALSE, error = "Spectral clustering failed on this block's similarity network."))
  names(clusters) <- rownames(layers[[1]])

  list(
    ok = TRUE, Wall = Wall, W = W, clusters = clusters, cluster_estimate = est,
    params = list(
      blocks = names(layers), k = k, alpha = alpha, t = NA_integer_,
      k_mode = if (k_manual) "manual" else "automatic", alpha_mode = if (alpha_manual) "manual" else "automatic",
      t_mode = "not applicable (single-omics, no fusion step)",
      n_clusters = n_clusters, cluster_mode = cluster_mode, n_samples = n, standardize = standardize,
      mode = "single_omics"
    )
  )
}

## ---------------------------------------------------------------------------
## 5. Clinical variable detection (spec sections 20, F-I) - inspects whatever
## sample_meta is actually attached to the selected dataset; never assumes a
## fixed schema. Survival needs a numeric time-like column PLUS a 2-level
## event-like column; short of both, survival is reported unavailable rather
## than guessed at.
## ---------------------------------------------------------------------------

SFC_TIME_COL_PATTERN <- "(?i)(time|duration|_os$|^os$|survival|pfs)"
SFC_EVENT_COL_PATTERN <- "(?i)(event|status|censor|death|relapse|progression)"

sfc_detect_clinical <- function(sample_meta) {
  if (is.null(sample_meta) || ncol(sample_meta) == 0) return(list(categorical = character(0), continuous = character(0), survival = NULL))
  classify <- function(col) {
    x <- sample_meta[[col]]
    x_nn <- x[!is.na(x)]
    if (length(x_nn) < 2) return(NA_character_)
    n_unique <- length(unique(x_nn))
    if (n_unique < 2) return(NA_character_)
    if (is.numeric(x_nn) && n_unique > max(5, floor(length(x_nn) / 4))) "continuous" else "categorical"
  }
  types <- stats::setNames(vapply(colnames(sample_meta), classify, character(1)), colnames(sample_meta))
  cont_cols <- names(types)[!is.na(types) & types == "continuous"]
  cat_cols  <- names(types)[!is.na(types) & types == "categorical"]

  time_candidates <- cont_cols[grepl(SFC_TIME_COL_PATTERN, cont_cols, perl = TRUE)]
  event_candidates <- colnames(sample_meta)[grepl(SFC_EVENT_COL_PATTERN, colnames(sample_meta), perl = TRUE)]
  event_candidates <- Filter(function(col) length(unique(sample_meta[[col]][!is.na(sample_meta[[col]])])) == 2, event_candidates)
  survival <- if (length(time_candidates) > 0 && length(event_candidates) > 0) list(time_col = time_candidates[1], event_col = event_candidates[1]) else NULL

  if (!is.null(survival)) {
    cont_cols <- setdiff(cont_cols, survival$time_col)
    cat_cols <- setdiff(cat_cols, survival$event_col)
  }
  list(categorical = cat_cols, continuous = cont_cols, survival = survival)
}

## ---------------------------------------------------------------------------
## 6. Clinical association tests (spec sections 21-23) - always post-hoc
## (never used to choose K/alpha/T/cluster count), always report which test
## and how many matched samples it ran on. `sfc_clinical_run()` BH-corrects
## across exactly the variables the user selected together, never against a
## silently larger set.
## ---------------------------------------------------------------------------

sfc_test_categorical <- function(clusters, x) {
  common <- intersect(names(clusters), names(x))
  common <- common[!is.na(x[common])]
  if (length(common) < 6) return(list(ok = FALSE, error = "Fewer than 6 matched, non-missing observations for this variable."))
  cl <- factor(clusters[common]); v <- factor(x[common])
  if (nlevels(cl) < 2 || nlevels(v) < 2) return(list(ok = FALSE, error = "Fewer than 2 levels after restricting to matched samples."))
  tab <- table(cluster = cl, value = v)
  p <- tryCatch(stats::fisher.test(tab, simulate.p.value = nrow(tab) * ncol(tab) > 20)$p.value, error = function(e) NA_real_)
  chi <- tryCatch(suppressWarnings(stats::chisq.test(tab)), error = function(e) NULL)
  cramers_v <- if (!is.null(chi)) sqrt(unname(chi$statistic) / (sum(tab) * (min(dim(tab)) - 1))) else NA_real_
  list(ok = TRUE, table = tab, p_value = p, effect = cramers_v, effect_label = "Cramer's V", n = length(common), test = "Fisher's exact test")
}

sfc_test_continuous <- function(clusters, x) {
  common <- intersect(names(clusters), names(x))
  common <- common[!is.na(x[common])]
  if (length(common) < 6) return(list(ok = FALSE, error = "Fewer than 6 matched, non-missing observations for this variable."))
  cl <- factor(clusters[common]); v <- as.numeric(x[common])
  if (nlevels(cl) < 2) return(list(ok = FALSE, error = "Fewer than 2 clusters among matched samples."))
  kw <- tryCatch(stats::kruskal.test(v, cl), error = function(e) NULL)
  if (is.null(kw)) return(list(ok = FALSE, error = "Kruskal-Wallis test failed for this variable."))
  agg <- stats::aggregate(v, by = list(cluster = cl), FUN = function(z) c(median = stats::median(z), n = length(z)))
  list(ok = TRUE, p_value = kw$p.value, statistic = unname(kw$statistic), n = length(common), test = "Kruskal-Wallis test",
       summary = data.frame(cluster = agg$cluster, median = agg$x[, "median"], n = agg$x[, "n"]))
}

sfc_test_survival <- function(clusters, time, event) {
  common <- Reduce(intersect, list(names(clusters), names(time), names(event)))
  common <- common[!is.na(time[common]) & !is.na(event[common])]
  if (length(common) < 6) return(list(ok = FALSE, error = "Fewer than 6 matched samples with complete survival data."))
  df <- data.frame(cluster = factor(clusters[common]), time = as.numeric(time[common]), event = as.numeric(event[common]))
  if (nlevels(df$cluster) < 2) return(list(ok = FALSE, error = "Fewer than 2 clusters among matched samples with survival data."))
  fit <- tryCatch(survival::survfit(survival::Surv(time, event) ~ cluster, data = df), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE, error = "Kaplan-Meier fit failed."))
  diff <- tryCatch(survival::survdiff(survival::Surv(time, event) ~ cluster, data = df), error = function(e) NULL)
  logrank_p <- if (!is.null(diff)) 1 - stats::pchisq(diff$chisq, length(diff$n) - 1) else NA_real_
  hr <- NULL
  if (nlevels(df$cluster) == 2) {
    cox <- tryCatch(survival::coxph(survival::Surv(time, event) ~ cluster, data = df), error = function(e) NULL)
    if (!is.null(cox)) {
      s <- summary(cox)
      hr <- list(hr = unname(s$conf.int[1, "exp(coef)"]), lo = unname(s$conf.int[1, "lower .95"]), hi = unname(s$conf.int[1, "upper .95"]))
    }
  }
  list(ok = TRUE, fit = fit, data = df, logrank_p = logrank_p, hr = hr, n = length(common), n_clusters = nlevels(df$cluster))
}

sfc_km_risk_table <- function(surv, n_points = 6) {
  if (is.null(surv) || !isTRUE(surv$ok)) return(NULL)
  rng <- range(surv$data$time, na.rm = TRUE)
  times <- unique(pmax(0, pretty(rng, n = n_points)))
  times <- times[times <= rng[2]]
  s <- tryCatch(summary(surv$fit, times = times, extend = TRUE), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  data.frame(cluster = sub("^cluster=", "", as.character(s$strata)), time = s$time, n_risk = s$n.risk, n_event = s$n.event, survival = round(s$surv, 3))
}

sfc_clinical_run <- function(clusters, sample_meta, vars, kind = c("categorical", "continuous")) {
  kind <- match.arg(kind)
  results <- lapply(vars, function(v) {
    x <- stats::setNames(sample_meta[[v]], rownames(sample_meta))
    r <- if (identical(kind, "categorical")) sfc_test_categorical(clusters, x) else sfc_test_continuous(clusters, x)
    r$variable <- v
    r
  })
  names(results) <- vars
  pvals <- vapply(results, function(r) if (isTRUE(r$ok)) r$p_value else NA_real_, numeric(1))
  fdr <- stats::p.adjust(pvals, method = "BH")
  for (i in seq_along(results)) results[[i]]$p_fdr <- unname(fdr[i])
  results
}

## ---------------------------------------------------------------------------
## 7. Stability (spec section 18, REQUIRED) - repeated subsampling + ARI
## against the full-cohort ("reference") clustering, using the exact same
## sfc_snf_run()/mi_ari() this module already computed its headline result
## with. Thresholds are fixed constants, stated once, applied identically
## everywhere a verdict is shown (never re-derived per call).
## ---------------------------------------------------------------------------

SFC_STABILITY_THRESHOLDS <- list(stable = 0.75, moderate = 0.5)

sfc_stability_verdict <- function(mean_ari) {
  if (is.na(mean_ari)) return("Not computable")
  if (mean_ari >= SFC_STABILITY_THRESHOLDS$stable) "Stable"
  else if (mean_ari >= SFC_STABILITY_THRESHOLDS$moderate) "Moderately stable"
  else "Unstable"
}

sfc_stability_run <- function(layers, ref_clusters, params, n_resamples = 20, subsample_frac = 0.8, seed = 1) {
  ids <- names(ref_clusters)
  n <- length(ids)
  n_sub <- max(MI_MIN_MATCHED_SAMPLES, round(n * subsample_frac))
  if (n_sub < MI_MIN_MATCHED_SAMPLES + 2 || n_sub >= n) return(list(ok = FALSE, error = "Cohort too small for a meaningful resampling-based stability check at this subsample fraction."))
  set.seed(seed)
  ari_vals <- numeric(0); failures <- 0L
  for (i in seq_len(n_resamples)) {
    sub_ids <- sample(ids, n_sub)
    sub_layers <- lapply(layers, function(m) m[sub_ids, , drop = FALSE])
    res <- tryCatch(sfc_snf_run(sub_layers, params), error = function(e) list(ok = FALSE))
    if (!isTRUE(res$ok) || length(unique(res$clusters)) < 2) { failures <- failures + 1L; next }
    a <- mi_ari(res$clusters, ref_clusters[sub_ids])
    if (!is.na(a)) ari_vals <- c(ari_vals, a) else failures <- failures + 1L
  }
  if (length(ari_vals) < 3) return(list(ok = FALSE, error = "Too few successful resamples to summarize stability - try more resamples or a larger subsample fraction."))
  mean_ari <- mean(ari_vals); sd_ari <- stats::sd(ari_vals)
  list(ok = TRUE, ari = ari_vals, mean_ari = mean_ari, sd_ari = sd_ari, n_resamples = length(ari_vals), n_requested = n_resamples,
       n_failed = failures, subsample_frac = subsample_frac, n_sub = n_sub, seed = seed, verdict = sfc_stability_verdict(mean_ari))
}

## ---------------------------------------------------------------------------
## 8. Parameter sensitivity (spec section 19, best-effort) - reruns clustering
## at the feasible low/high end of each tunable parameter and compares each
## result to the reference clustering by ARI. T is skipped for single-omics
## (no fusion step exists to vary).
## ---------------------------------------------------------------------------

sfc_sensitivity_run <- function(layers, base_params, ref_clusters, seed = 1) {
  set.seed(seed)
  n <- nrow(layers[[1]])
  k_range <- mi_snf_feasible_k_range(n)
  variations <- list(
    K = list(low = max(2, k_range$min), high = k_range$max),
    Alpha = list(low = MI_SNF_ALPHA_RANGE$min, high = MI_SNF_ALPHA_RANGE$max)
  )
  if (length(layers) >= 2) variations$T <- list(low = 10, high = 50)

  rows <- list()
  for (pname in names(variations)) {
    for (level in c("low", "high")) {
      val <- variations[[pname]][[level]]
      p <- base_params
      if (identical(pname, "K")) { p$k_mode <- "manual"; p$k <- val }
      if (identical(pname, "Alpha")) { p$alpha_mode <- "manual"; p$alpha <- val }
      if (identical(pname, "T")) { p$t_mode <- "manual"; p$t <- val }
      res <- tryCatch(sfc_snf_run(layers, p), error = function(e) list(ok = FALSE))
      ari <- if (isTRUE(res$ok) && length(unique(res$clusters)) >= 2) mi_ari(res$clusters, ref_clusters[names(res$clusters)]) else NA_real_
      rows[[length(rows) + 1]] <- data.frame(parameter = pname, level = level, value = val, ari_vs_reference = ari)
    }
  }
  df <- do.call(rbind, rows)
  by_param <- stats::aggregate(ari_vs_reference ~ parameter, data = df, FUN = function(x) mean(x, na.rm = TRUE))
  by_param$sensitivity <- ifelse(is.na(by_param$ari_vs_reference), "Not computable",
                                  ifelse(by_param$ari_vs_reference >= SFC_STABILITY_THRESHOLDS$stable, "Low sensitivity", "High sensitivity"))
  list(ok = TRUE, detail = df, summary = by_param)
}

## ---------------------------------------------------------------------------
## 9. Feature ranking (spec section 24) - per feature, Kruskal-Wallis test
## across the already-computed clusters (never re-derives clusters). Labeled
## as association-with-clusters, distinct from "features used to construct
## the network" (every feature in the selected block was used for that).
## ---------------------------------------------------------------------------

sfc_feature_ranking <- function(mat, clusters, block_label, top_n = 25) {
  common <- intersect(rownames(mat), names(clusters))
  if (length(common) < 6) return(list(ok = FALSE, error = "Too few matched samples for feature ranking."))
  cl <- factor(clusters[common])
  if (nlevels(cl) < 2) return(list(ok = FALSE, error = "Fewer than 2 clusters among matched samples."))
  m <- mat[common, , drop = FALSE]
  pvals <- apply(m, 2, function(col) tryCatch(stats::kruskal.test(col, cl)$p.value, error = function(e) NA_real_))
  fdr <- stats::p.adjust(pvals, method = "BH")
  df <- data.frame(block = block_label, feature = colnames(m), p_value = pvals, p_fdr = fdr, stringsAsFactors = FALSE)
  df <- df[order(df$p_value), , drop = FALSE]
  list(ok = TRUE, table = utils::head(df, top_n), n_features_tested = ncol(m))
}

## ---------------------------------------------------------------------------
## 10. Async run wrapper - bundles the headline clustering with a default
## stability check in one call (mod_multi_stratification.R's Run handler).
## Defined at file scope, not nested inside the module's server closure, so
## future::future_promise()'s automatic global-variable export resolves it
## unambiguously - the same reason mi_diablo_run()/mi_snf_run() are called
## directly (never via a server-local wrapper) inside Multi-omics
## Integration's own ExtendedTask.
## ---------------------------------------------------------------------------

sfc_snf_run_with_stability <- function(layers, params) {
  res <- sfc_snf_run(layers, params)
  if (!isTRUE(res$ok)) return(list(ok = FALSE, error = res$error))
  stab <- tryCatch(sfc_stability_run(layers, res$clusters, params, n_resamples = 20, subsample_frac = 0.8, seed = 1), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  list(ok = TRUE, res = res, stability = stab)
}

## ---------------------------------------------------------------------------
## 11. Analysis summary (spec section 27) - short bullet lines only, no
## fabricated narrative; every value is read back from the actual run.
## ---------------------------------------------------------------------------

sfc_summary_lines <- function(res, stability = NULL, clinical_note = NULL) {
  p <- res$params
  cl_tab <- table(res$clusters)
  c(
    sprintf("Modalities: %s%s", paste(p$blocks, collapse = " + "), if (identical(p$mode, "single_omics")) " (Single-Omics Clustering)" else ""),
    sprintf("Matched patients: %d", p$n_samples),
    sprintf("Clusters: %d (%s; sizes: %s)", p$n_clusters, p$cluster_mode, paste(as.integer(cl_tab), collapse = "/")),
    sprintf("K: %s (%s)", p$k, p$k_mode),
    sprintf("Alpha: %s (%s)", if (is.na(p$alpha)) "-" else sprintf("%.2f", p$alpha), p$alpha_mode),
    sprintf("T: %s", if (is.na(p$t)) "not applicable (single-omics)" else sprintf("%s (%s)", p$t, p$t_mode)),
    sprintf("Cluster stability: %s", if (!is.null(stability) && isTRUE(stability$ok)) sprintf("%s (mean ARI = %.2f across %d resamples)", stability$verdict, stability$mean_ari, stability$n_resamples) else "Not computed"),
    if (!is.null(clinical_note)) clinical_note
  )
}
