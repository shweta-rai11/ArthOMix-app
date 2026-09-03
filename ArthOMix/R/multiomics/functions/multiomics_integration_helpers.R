## R/multiomics/functions/multiomics_integration_helpers.R
## Pure data-processing logic for the live DIABLO/SNF/Compare engine mounted
## in mod_multi_integration.R - the Multi-Omics module's own "run a real

mi_preloaded_cell_dataset <- function(cell_key) {
  cell <- multi_cell_by_key(cell_key)
  if (is.null(cell)) return(list(ok = FALSE, error = "Unknown preloaded analysis cell."))
  fit_res <- multi_diablo_fit(cell)
  if (!fit_res$ok) return(list(ok = FALSE, error = fit_res$error))
  fit <- fit_res$fit
  layers <- fit$X
  if (is.null(layers) || length(layers) < 1) return(list(ok = FALSE, error = "This cell's saved fit has no per-block data."))
  layers <- lapply(layers, function(m) { storage.mode(m) <- "double"; m })
  names(layers) <- vapply(names(layers), function(nm) MULTI_BLOCK_LABELS[[nm]] %||% nm, character(1))
  y <- fit$Y
  meta <- if (!is.null(y)) data.frame(row.names = rownames(layers[[1]]), outcome = as.character(y), stringsAsFactors = FALSE) else NULL
  list(
    ok = TRUE, layers = layers, sample_meta = meta,
    outcome_col = if (!is.null(meta)) "outcome" else NULL,
    label = cell$label,
    provenance = sprintf(
      "Preloaded RA anti-TNF cohort - %s (%s). This re-runs DIABLO/SNF live - a new fit, not the pipeline's saved result. See \"Patient Stratification\" for precomputed SNF clusters, or \"Biomarker Discovery\" for a stability/ROC-reported DIABLO panel.",
      cell$label, paste(sprintf("%s: %d features", names(layers), vapply(layers, ncol, integer(1))), collapse = ", ")
    )
  )
}

MI_MIN_MATCHED_SAMPLES <- 3
MI_SAMPLE_MISMATCH_MESSAGE <- "Integration cannot be performed because the omics datasets do not contain reliably matched samples."

mi_validate_dataset <- function(layers, sample_meta = NULL, outcome_col = NULL) {
  layers <- Filter(Negate(is.null), layers %||% list())
  if (length(layers) == 0) return(list(ok = FALSE, error = "No omics blocks are available.", n_blocks = 0))

  per_block <- stats::setNames(lapply(names(layers), function(nm) multi_live_validate_matrix(layers[[nm]], layer_label = nm)), names(layers))

  overlap <- if (length(layers) >= 2) multi_live_sample_overlap(layers) else list(ok = FALSE, error = "Only one omics block is available - integration needs at least two.")
  n_shared <- if (isTRUE(overlap$ok)) overlap$n_shared else 0
  shared_ids <- if (isTRUE(overlap$ok)) overlap$shared_ids else character(0)
  reliable_matching <- length(layers) >= 2 && isTRUE(overlap$ok) && n_shared >= MI_MIN_MATCHED_SAMPLES

  outcome_summary <- if (!is.null(sample_meta) && !is.null(outcome_col) && outcome_col %in% colnames(sample_meta)) {
    mi_outcome_summary(sample_meta, outcome_col, shared_ids)
  } else NULL

  list(
    ok = TRUE, n_blocks = length(layers), block_labels = names(layers),
    per_block = per_block, overlap = overlap, n_shared = n_shared, shared_ids = shared_ids,
    reliable_matching = reliable_matching, mismatch_message = if (!reliable_matching) MI_SAMPLE_MISMATCH_MESSAGE else NULL,
    outcome = outcome_summary
  )
}

mi_outcome_summary <- function(meta, col, sample_ids = NULL) {
  if (is.null(meta) || !col %in% colnames(meta)) return(NULL)
  vals <- stats::setNames(meta[[col]], rownames(meta))
  if (!is.null(sample_ids) && length(sample_ids) > 0) vals <- vals[intersect(sample_ids, names(vals))]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(list(column = col, n = 0, type = "categorical", n_classes = NA_integer_, class_counts = NULL, imbalanced = FALSE, values = vals))
  n_unique <- length(unique(vals))
  is_numeric_like <- identical(ch_classify_column(vals), "continuous")
  fac <- if (is_numeric_like) NULL else factor(vals)
  tab <- if (!is.null(fac)) table(fac) else NULL
  list(
    column = col, n = length(vals), type = if (is_numeric_like) "continuous" else "categorical",
    n_classes = if (!is.null(tab)) length(tab) else NA_integer_,
    class_counts = tab,
    imbalanced = !is.null(tab) && length(tab) == 2 && (max(tab) / sum(tab)) > 0.7,
    values = vals
  )
}

mi_diablo_eligibility <- function(validation, outcome_summary) {
  if (is.null(validation) || !isTRUE(validation$ok) || validation$n_blocks < 2) {
    return(list(ok = FALSE, reason = "DIABLO needs at least two omics blocks."))
  }
  if (!isTRUE(validation$reliable_matching)) return(list(ok = FALSE, reason = MI_SAMPLE_MISMATCH_MESSAGE))
  if (is.null(outcome_summary) || !identical(outcome_summary$type, "categorical") || is.na(outcome_summary$n_classes) || outcome_summary$n_classes < 2) {
    return(list(ok = FALSE, reason = "DIABLO requires a categorical outcome. No suitable outcome was detected."))
  }
  if (any(outcome_summary$class_counts < 3)) {
    return(list(ok = FALSE, reason = sprintf(
      "The smallest outcome class has %d sample(s) - at least 3 per class are needed for cross-validated DIABLO.",
      min(outcome_summary$class_counts))))
  }
  list(ok = TRUE, reason = NULL)
}

mi_diablo_keepx_grid <- function(n_features) {
  base <- c(5, 10, 20, 50, 100, 150, 200, 300)
  cand <- unique(pmin(base[base <= n_features], n_features))
  if (length(cand) == 0) cand <- n_features
  sort(unique(cand))
}

mi_diablo_feasible_ncomp <- function(n_classes, min_class_n) {
  max_ncomp <- max(1, min(5, n_classes + 1, min_class_n - 1))
  seq_len(max_ncomp)
}

mi_diablo_max_folds <- function(min_class_n) max(2, min(10, min_class_n))

mi_diablo_feasible_folds <- function(min_class_n, requested = NULL) {
  max_folds <- mi_diablo_max_folds(min_class_n)
  if (is.null(requested)) return(min(5, max_folds))
  max(2, min(as.integer(requested), max_folds))
}

MI_DIABLO_LOO_MAX_N <- 60
mi_diablo_loo_feasible <- function(n_samples) n_samples <= MI_DIABLO_LOO_MAX_N

mi_diablo_feasible_repeats <- function(n_samples) {
  if (n_samples < 30) 10L else if (n_samples < 100) 5L else 3L
}

mi_diablo_design <- function(block_names, mode = c("automatic", "custom"), custom = NULL) {
  mode <- match.arg(mode)
  n <- length(block_names)
  design <- matrix(0.1, nrow = n, ncol = n, dimnames = list(block_names, block_names))
  diag(design) <- 0
  if (identical(mode, "custom") && !is.null(custom)) {
    custom <- as.matrix(custom[block_names, block_names, drop = FALSE])
    diag(custom) <- 0
    design <- custom
  }
  design
}

mi_diablo_run <- function(layers, outcome, sample_ids, params = list()) {
  X <- lapply(layers, function(m) m[sample_ids, , drop = FALSE])
  Y <- droplevels(factor(outcome[sample_ids]))
  if (nlevels(Y) < 2) return(list(ok = FALSE, error = "Fewer than two outcome classes remain in the matched samples."))
  block_names <- names(X)

  class_tab <- table(Y)
  min_class_n <- min(class_tab)

  design <- mi_diablo_design(block_names, mode = params$design_mode %||% "automatic", custom = params$design_custom)

  ncomp_choices <- mi_diablo_feasible_ncomp(nlevels(Y), min_class_n)
  ncomp <- if (identical(params$ncomp_mode %||% "automatic", "manual") && !is.null(params$ncomp)) {
    max(1L, min(as.integer(params$ncomp), max(ncomp_choices)))
  } else max(ncomp_choices)

  manual_validation <- identical(params$validation_mode %||% "automatic", "manual")
  requested_method <- if (manual_validation) (params$validation_method %||% "mfold") else "mfold"
  use_loo <- identical(requested_method, "loo") && mi_diablo_loo_feasible(length(Y))
  loo_downgraded <- identical(requested_method, "loo") && !use_loo
  validation_method <- if (use_loo) "loo" else "Mfold"

  folds <- mi_diablo_feasible_folds(min_class_n, if (manual_validation) params$folds else NULL)
  nrepeat <- if (use_loo) 1L else if (manual_validation && !is.null(params$nrepeat)) max(1L, as.integer(params$nrepeat)) else mi_diablo_feasible_repeats(length(Y))

  dist_choice <- params$distance %||% "automatic"
  keepx_mode <- params$keepx_mode %||% "automatic"

  seed <- if (!is.null(params$seed)) as.integer(params$seed) else NULL

  if (identical(keepx_mode, "manual") && !is.null(params$keepx_manual)) {
    keepX <- params$keepx_manual
  } else {
    grid <- stats::setNames(lapply(X, function(m) mi_diablo_keepx_grid(ncol(m))), block_names)
    tuned <- tryCatch(
      mixOmics::tune.block.splsda(
        X = X, Y = Y, ncomp = ncomp, test.keepX = grid, design = design,
        validation = validation_method, folds = folds, nrepeat = nrepeat,
        dist = if (identical(dist_choice, "automatic")) "max.dist" else dist_choice,
        measure = "BER", progressBar = FALSE, near.zero.var = TRUE, seed = seed,
        scale = isTRUE(params$scale %||% TRUE)
      ),
      error = function(e) e
    )
    if (inherits(tuned, "error")) return(list(ok = FALSE, error = paste("DIABLO keepX tuning failed:", conditionMessage(tuned))))
    keepX <- tuned$choice.keepX
  }

  fit <- tryCatch(
    mixOmics::block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design, near.zero.var = TRUE, scale = isTRUE(params$scale %||% TRUE)),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, error = paste("DIABLO fit failed:", conditionMessage(fit))))

  perf_dist_arg <- if (identical(dist_choice, "automatic")) "all" else dist_choice
  perf_res <- tryCatch(
    mixOmics::perf(fit, validation = validation_method, folds = folds, nrepeat = nrepeat, dist = perf_dist_arg, auc = TRUE, progressBar = FALSE, seed = seed),
    error = function(e) e
  )
  if (inherits(perf_res, "error")) return(list(ok = FALSE, error = paste("DIABLO performance assessment failed:", conditionMessage(perf_res))))
  dist_candidates <- names(perf_res$WeightedVote.error.rate)
  ber_by_dist <- vapply(dist_candidates, function(d) {
    tryCatch(perf_res$WeightedVote.error.rate[[d]]["Overall.BER", ncomp], error = function(e) NA_real_)
  }, numeric(1))
  if (all(is.na(ber_by_dist))) return(list(ok = FALSE, error = "Cross-validated performance could not be computed for this configuration."))
  resolved_dist <- names(ber_by_dist)[which.min(ber_by_dist)]

  list(
    ok = TRUE, fit = fit, perf = perf_res,
    params = list(
      blocks = block_names, ncomp = ncomp, keepX = keepX,
      design_mode = params$design_mode %||% "automatic", design = design,
      validation_method = validation_method, folds = folds, nrepeat = nrepeat, distance = resolved_dist,
      distance_mode = dist_choice, n_samples = length(Y),
      classes = names(class_tab), class_counts = as.integer(class_tab),
      loo_downgraded = loo_downgraded, seed = seed,
      scale = isTRUE(params$scale %||% TRUE)
    )
  )
}

mi_diablo_performance_summary <- function(diablo_res) {
  if (!isTRUE(diablo_res$ok)) return(NULL)
  p <- diablo_res$perf; par <- diablo_res$params
  dist <- par$distance; ncomp <- par$ncomp
  mat <- tryCatch(p$WeightedVote.error.rate[[dist]], error = function(e) NULL)
  if (is.null(mat)) return(NULL)
  ber <- mat["Overall.BER", ncomp]; er <- mat["Overall.ER", ncomp]
  per_class <- mat[setdiff(rownames(mat), c("Overall.ER", "Overall.BER")), ncomp, drop = TRUE]
  auc_mat <- tryCatch(p$auc[[paste0("comp", ncomp)]], error = function(e) NULL)
  auc_df <- if (!is.null(auc_mat)) data.frame(comparison = rownames(auc_mat), AUC = auc_mat[, "AUC"], p_value = auc_mat[, "p-value"], row.names = NULL) else NULL
  list(ber = unname(ber), overall_error = unname(er), per_class_error = per_class, auc = auc_df, distance = dist, ncomp = ncomp)
}

mi_diablo_selected_features_df <- function(fit) {
  if (is.null(fit)) return(NULL)
  blocks <- names(fit$X)
  ncomp <- fit$ncomp[1]
  out <- do.call(rbind, lapply(blocks, function(b) {
    do.call(rbind, lapply(seq_len(ncomp), function(c) {
      sv <- tryCatch(mixOmics::selectVar(fit, block = b, comp = c), error = function(e) NULL)
      entry <- sv[[b]]
      if (is.null(entry) || length(entry$name) == 0) return(NULL)
      data.frame(
        block = b, component = c, feature = entry$name, loading = entry$value[, 1],
        selected = TRUE, stringsAsFactors = FALSE
      )
    }))
  }))
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out[order(out$block, out$component, -abs(out$loading)), , drop = FALSE]
}

mi_diablo_sample_scores_df <- function(fit, outcome) {
  blocks <- setdiff(names(fit$variates), "Y")
  if (length(blocks) == 0) return(NULL)
  mats <- lapply(blocks, function(b) fit$variates[[b]][, 1])
  ids <- rownames(fit$variates[[blocks[1]]])
  consensus <- rowMeans(do.call(cbind, mats))
  data.frame(patient_id = ids, response = as.character(outcome[ids]), score = consensus, stringsAsFactors = FALSE)
}

mi_diablo_panel_df_for_plot <- function(sel_df, comp = 1) {
  if (is.null(sel_df)) return(NULL)
  d <- sel_df[sel_df$component == comp, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  data.frame(feature = d$feature, loading = d$loading, view = d$block, stringsAsFactors = FALSE)
}

mi_diablo_stability_df <- function(diablo_res, block, comp) {
  stable <- tryCatch(diablo_res$perf$features$stable, error = function(e) NULL)
  if (is.null(stable)) return(NULL)
  comp_key <- paste0("comp", comp)
  tabs <- lapply(stable, function(rep_entry) rep_entry[[block]][[comp_key]])
  tabs <- Filter(Negate(is.null), tabs)
  if (length(tabs) == 0) return(NULL)
  all_feats <- unique(unlist(lapply(tabs, names)))
  freq <- vapply(all_feats, function(f) mean(vapply(tabs, function(t) if (f %in% names(t)) t[[f]] else 0, numeric(1))), numeric(1))
  df <- data.frame(feature = names(freq), stability = as.numeric(freq), stringsAsFactors = FALSE)
  df[order(-df$stability), , drop = FALSE]
}

mi_snf_eligibility <- function(validation) {
  if (is.null(validation) || !isTRUE(validation$ok) || validation$n_blocks < 2) {
    return(list(ok = FALSE, reason = "SNF unavailable: fewer than two compatible matched omics blocks were detected."))
  }
  if (!isTRUE(validation$reliable_matching)) return(list(ok = FALSE, reason = MI_SAMPLE_MISMATCH_MESSAGE))
  n_missing <- vapply(validation$per_block, function(v) if (isTRUE(v$ok)) v$n_missing else NA_integer_, numeric(1))
  if (any(!is.na(n_missing) & n_missing > 0)) {
    return(list(ok = FALSE, reason = "SNF requires complete data - one or more selected blocks still have missing values. Resolve missing values (impute or remove) before running SNF."))
  }
  list(ok = TRUE, reason = NULL)
}

mi_snf_feasible_k_range <- function(n_samples) {
  lo <- max(3, floor(n_samples / 10))
  hi <- max(lo + 1, min(n_samples - 1, 30))
  list(min = lo, max = hi, default = max(lo, min(hi, round(sqrt(n_samples)) + 5)))
}

MI_SNF_ALPHA_RANGE <- list(min = 0.3, max = 0.8, default = 0.5)
MI_SNF_T_CANDIDATES <- c(10L, 20L, 30L, 50L)

mi_snf_affinity <- function(mat, k, alpha, standardize = TRUE) {
  m <- if (isTRUE(standardize)) SNFtool::standardNormalization(mat) else mat
  d <- SNFtool::dist2(as.matrix(m), as.matrix(m))
  SNFtool::affinityMatrix(d, K = k, sigma = alpha)
}

mi_snf_adaptive_t <- function(Wall, k, candidates = MI_SNF_T_CANDIDATES) {
  fused <- lapply(candidates, function(t) SNFtool::SNF(Wall, K = k, t = t))
  for (i in seq_along(candidates)[-length(candidates)]) {
    r <- suppressWarnings(stats::cor(as.vector(fused[[i]]), as.vector(fused[[i + 1]])))
    if (!is.na(r) && r >= 0.995) return(list(t = candidates[i], W = fused[[i]]))
  }
  list(t = candidates[length(candidates)], W = fused[[length(candidates)]])
}

MI_SNF_CLUSTER_METHODS <- c("Spectral clustering (SNFtool::spectralClustering)" = "spectral",
                            "Hierarchical clustering (average linkage, stats::hclust)" = "hierarchical",
                            "PAM / k-medoids (cluster::pam)" = "pam")

mi_snf_network_dist <- function(W) stats::as.dist(1 - W / max(W))

mi_cluster_from_network <- function(W, k, method = c("spectral", "hierarchical", "pam")) {
  method <- match.arg(method)
  if (identical(method, "spectral")) {
    return(tryCatch(SNFtool::spectralClustering(W, K = k), error = function(e) NULL))
  }
  d <- mi_snf_network_dist(W)
  if (identical(method, "hierarchical")) {
    hc <- tryCatch(stats::hclust(d, method = "average"), error = function(e) NULL)
    if (is.null(hc)) return(NULL)
    grp <- tryCatch(stats::cutree(hc, k = k), error = function(e) NULL)
    return(grp)
  }
  if (!requireNamespace("cluster", quietly = TRUE)) return(NULL)
  pm <- tryCatch(cluster::pam(d, k = k, diss = TRUE), error = function(e) NULL)
  if (is.null(pm)) return(NULL)
  pm$clustering
}

mi_snf_auto_tune <- function(layers, standardize = TRUE, k_fixed = NULL, alpha_fixed = NULL, cluster_method = "spectral") {
  n <- nrow(layers[[1]])
  k_range <- mi_snf_feasible_k_range(n)
  k_grid <- if (!is.null(k_fixed)) max(2, min(as.integer(k_fixed), n - 1)) else unique(round(c(k_range$min, k_range$default, k_range$max)))
  alpha_grid <- if (!is.null(alpha_fixed)) alpha_fixed else c(0.3, 0.5, 0.7)
  best <- NULL
  for (k in k_grid) for (alpha in alpha_grid) {
    Wall <- lapply(layers, mi_snf_affinity, k = k, alpha = alpha, standardize = standardize)
    W <- SNFtool::SNF(Wall, K = k, t = 20)
    est <- tryCatch(SNFtool::estimateNumberOfClustersGivenGraph(W, NUMC = 2:min(6, n - 1)), error = function(e) NULL)
    nc <- if (!is.null(est)) est[["Eigen-gap best"]] else 2
    grp <- mi_cluster_from_network(W, nc, cluster_method)
    if (is.null(grp) || length(unique(grp)) < 2) next
    sil <- tryCatch(mean(cluster::silhouette(grp, mi_snf_network_dist(W))[, "sil_width"]), error = function(e) NA_real_)
    if (is.na(sil)) next
    if (is.null(best) || sil > best$sil) best <- list(k = k, alpha = alpha, sil = sil, Wall = Wall, W = W, nc = nc)
  }
  if (is.null(best)) {
    k <- if (!is.null(k_fixed)) max(2, min(as.integer(k_fixed), n - 1)) else k_range$default
    alpha <- alpha_fixed %||% MI_SNF_ALPHA_RANGE$default
    Wall <- lapply(layers, mi_snf_affinity, k = k, alpha = alpha, standardize = standardize)
    best <- list(k = k, alpha = alpha, sil = NA_real_, Wall = Wall, W = SNFtool::SNF(Wall, K = k, t = 20), nc = NA_integer_)
  }
  best
}

mi_snf_run <- function(layers, params = list()) {
  standardize <- isTRUE(params$standardize %||% TRUE)
  n <- nrow(layers[[1]])
  k_manual <- identical(params$k_mode %||% "automatic", "manual")
  alpha_manual <- identical(params$alpha_mode %||% "automatic", "manual")
  t_manual <- identical(params$t_mode %||% "automatic", "manual")
  cluster_method <- params$cluster_method %||% "spectral"
  seed <- as.integer(params$seed %||% 1)
  set.seed(seed)

  auto <- mi_snf_auto_tune(
    layers, standardize = standardize,
    k_fixed = if (k_manual) params$k else NULL,
    alpha_fixed = if (alpha_manual) params$alpha else NULL,
    cluster_method = cluster_method
  )
  k <- auto$k; alpha <- auto$alpha; Wall <- auto$Wall

  if (t_manual) {
    t_choice <- max(1L, as.integer(params$t %||% 20))
    W <- SNFtool::SNF(Wall, K = k, t = t_choice)
  } else {
    t_res <- mi_snf_adaptive_t(Wall, k)
    t_choice <- t_res$t; W <- t_res$W
  }

  max_k_clusters <- min(6, n - 1)
  est <- tryCatch(SNFtool::estimateNumberOfClustersGivenGraph(W, NUMC = 2:max_k_clusters), error = function(e) NULL)
  cluster_mode <- params$cluster_mode %||% "automatic"
  n_clusters <- if (identical(cluster_mode, "manual") && !is.null(params$n_clusters)) {
    max(2, min(as.integer(params$n_clusters), max_k_clusters))
  } else if (!is.null(est)) est[["Eigen-gap best"]] else 2

  clusters <- mi_cluster_from_network(W, n_clusters, cluster_method)
  if (is.null(clusters)) return(list(ok = FALSE, error = sprintf("%s clustering failed on the fused network.", names(MI_SNF_CLUSTER_METHODS)[MI_SNF_CLUSTER_METHODS == cluster_method])))
  names(clusters) <- rownames(layers[[1]])

  list(
    ok = TRUE, Wall = Wall, W = W, clusters = clusters, cluster_estimate = est,
    params = list(
      blocks = names(layers), k = k, alpha = alpha, t = t_choice,
      k_mode = if (k_manual) "manual" else "automatic", alpha_mode = if (alpha_manual) "manual" else "automatic", t_mode = if (t_manual) "manual" else "automatic",
      n_clusters = n_clusters, cluster_mode = cluster_mode, n_samples = n, standardize = standardize,
      cluster_method = cluster_method, seed = seed
    )
  )
}

mi_snf_concordance <- function(snf_res) {
  if (!isTRUE(snf_res$ok)) return(NULL)
  Wall_plus <- c(snf_res$Wall, list(Fused = snf_res$W))
  nc <- snf_res$params$n_clusters
  mat <- tryCatch(SNFtool::concordanceNetworkNMI(Wall_plus, nc), error = function(e) NULL)
  if (is.null(mat)) return(NULL)
  labels <- names(Wall_plus)
  dimnames(mat) <- list(labels, labels)
  fused_row <- mat["Fused", setdiff(labels, "Fused"), drop = TRUE]
  data.frame(block = names(fused_row), concordance_with_fused = as.numeric(fused_row), row.names = NULL)
}

mi_snf_posthoc_outcome <- function(snf_res, outcome, sample_ids) {
  if (!isTRUE(snf_res$ok) || is.null(outcome)) return(NULL)
  common <- intersect(names(snf_res$clusters), intersect(sample_ids, names(outcome)))
  if (length(common) < 3) return(NULL)
  cl <- snf_res$clusters[common]; y <- factor(outcome[common])
  tab <- table(cluster = cl, outcome = y)
  fisher_p <- tryCatch(stats::fisher.test(tab, simulate.p.value = nrow(tab) * ncol(tab) > 20)$p.value, error = function(e) NA_real_)
  list(
    table = tab, fisher_p = fisher_p,
    nmi = tryCatch(SNFtool::calNMI(as.integer(factor(cl)), as.integer(y)), error = function(e) NA_real_),
    ari = mi_ari(cl, y)
  )
}

mi_ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  choose2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(choose2(tab))
  sum_a <- sum(choose2(rowSums(tab)))
  sum_b <- sum(choose2(colSums(tab)))
  expected <- sum_a * sum_b / choose2(n)
  max_index <- 0.5 * (sum_a + sum_b)
  if (max_index == expected) return(NA_real_)
  (sum_ij - expected) / (max_index - expected)
}

mi_diablo_selected_correlation_data <- function(layers, sel_df, block_a, block_b, sample_ids, method = "pearson") {
  feats_a <- unique(sel_df$feature[sel_df$block == block_a])
  feats_b <- unique(sel_df$feature[sel_df$block == block_b])
  if (length(feats_a) == 0 || length(feats_b) == 0) return(list(ok = FALSE, error = "No selected features for one of these blocks."))
  matA <- layers[[block_a]][sample_ids, feats_a, drop = FALSE]
  matB <- layers[[block_b]][sample_ids, feats_b, drop = FALSE]
  multi_live_correlation_heatmap_data(matA, matB, top_n = max(length(feats_a), length(feats_b)), method = method)
}

mi_compare_supervised <- function(layers, outcome, sample_ids, diablo_res) {
  if (!isTRUE(diablo_res$ok)) return(list(ok = FALSE, error = "Run DIABLO first - Compare needs its cross-validated performance."))
  mat_list <- lapply(layers, function(m) m[sample_ids, , drop = FALSE])
  y <- stats::setNames(outcome[sample_ids], sample_ids)
  single <- ch_evaluate_binary_outcome(mat_list, y, k_folds = diablo_res$params$folds)
  if (!isTRUE(single$ok)) return(list(ok = FALSE, error = single$error))
  perf_summary <- mi_diablo_performance_summary(diablo_res)
  diablo_auroc <- if (!is.null(perf_summary$auc) && nrow(perf_summary$auc) > 0) perf_summary$auc$AUC[1] else NA_real_
  rows <- do.call(rbind, lapply(names(single$per_view_auc), function(nm) {
    data.frame(model = nm, auroc = single$per_view_auc[[nm]], type = "Single-omics", stringsAsFactors = FALSE)
  }))
  rows <- rbind(rows, data.frame(model = "DIABLO (integrated)", auroc = diablo_auroc, type = "Integrated", stringsAsFactors = FALSE))
  rows <- rbind(rows, data.frame(model = "Chance / majority-class baseline", auroc = single$majority_baseline, type = "Baseline", stringsAsFactors = FALSE))
  list(ok = TRUE, table = rows, k_folds = single$k_folds, n = single$n, note = "DIABLO and the single-omics baselines each use their own k-fold splits (same k) - not a shared fold assignment.")
}

mi_compare_unsupervised <- function(snf_res) {
  if (!isTRUE(snf_res$ok)) return(list(ok = FALSE, error = "Run SNF first - Compare needs its fused clustering."))
  nc <- snf_res$params$n_clusters
  per_block <- lapply(names(snf_res$Wall), function(nm) {
    grp <- tryCatch(SNFtool::spectralClustering(snf_res$Wall[[nm]], K = nc), error = function(e) NULL)
    if (is.null(grp)) return(NULL)
    names(grp) <- rownames(snf_res$Wall[[nm]])
    data.frame(
      block = nm,
      nmi_vs_fused = tryCatch(SNFtool::calNMI(as.integer(factor(grp)), as.integer(factor(snf_res$clusters))), error = function(e) NA_real_),
      ari_vs_fused = mi_ari(grp, snf_res$clusters)
    )
  })
  per_block <- do.call(rbind, Filter(Negate(is.null), per_block))
  list(ok = TRUE, table = per_block, n_clusters = nc)
}
