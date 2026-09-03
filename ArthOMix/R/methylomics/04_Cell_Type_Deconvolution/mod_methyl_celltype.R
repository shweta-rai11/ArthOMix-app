## R/methylomics/04_Cell_Type_Deconvolution/mod_methyl_celltype.R
## Cell-Type Deconvolution submodule: estimates cell-type proportions from bulk methylation via EpiDISH
## (Houseman CP, RPC, CBS, hepidish two-stage) against built-in reference panels or a custom upload.

mod_methyl_celltype_config <- list(
  id = "celltype", title = "Cell-Type Deconvolution", icon = "people-group", group = "Data",
  description = "Estimates cell-type proportions from bulk methylation. Works on the loaded dataset, or upload data"
)

methyl_ct_reference_registry <- function() {
  specs <- list(
    list(id = "blood7",       label = "Blood - 7 cell type (Reinius/Houseman)",       object = "centDHSbloodDMC.m", tissue = "blood"),
    list(id = "blood7_compact", label = "Blood - 7 cell type, compact panel",         object = "centBloodSub.m",    tissue = "blood"),
    list(id = "blood12",      label = "Blood - 12 cell type (Salas et al.)",          object = "cent12CT.m",        tissue = "blood"),
    list(id = "blood12_450k", label = "Blood - 12 cell type (450K-restricted panel)", object = "cent12CT450k.m",    tissue = "blood"),
    list(id = "epifib",       label = "Epithelial / Fibroblast / Immune",             object = "centEpiFibIC.m",    tissue = "epithelial"),
    list(id = "epifibfat",    label = "Epithelial / Fibroblast / Fat / Immune",       object = "centEpiFibFatIC.m", tissue = "epithelial"),
    list(id = "blood19",      label = "Blood - extended 19-subtype panel",            object = "centCAB100i.m",     tissue = "blood")
  )
  lapply(specs, function(s) {
    mat <- tryCatch({
      if (!requireNamespace("EpiDISH", quietly = TRUE)) stop("EpiDISH not installed")
      getExportedValue("EpiDISH", s$object)
    }, error = function(e) NULL)
    if (is.null(mat)) return(c(s, list(available = FALSE, ncpg = NA_integer_, celltypes = character(0))))
    c(s, list(available = TRUE, ncpg = nrow(mat), celltypes = colnames(mat)))
  })
}

methyl_ct_get_reference <- function(ref_id) {
  reg <- methyl_ct_reference_registry()
  spec <- Find(function(s) identical(s$id, ref_id), reg)
  if (is.null(spec) || !isTRUE(spec$available)) return(NULL)
  getExportedValue("EpiDISH", spec$object)
}

methyl_ct_unavailable_methods <- function() {
  list(
    list(id = "methylresolver", label = "MethylResolver",
         reason = "The MethylResolver package is not installed in this deployment."),
    list(id = "idol", label = "IDOL-optimized library",
         reason = "IDOLOptimizedCpGs / FlowSorted.Blood.EPIC are not installed in this deployment."),
    list(id = "reffree", label = "Reference-free deconvolution",
         reason = "No reference-free methylation package (e.g. RefFreeEWAS, TOAST) is installed in this deployment.")
  )
}

methyl_ct_detect_scale <- function(mat) {
  v <- mat[is.finite(mat)]
  if (length(v) == 0) return(list(scale = "beta", note = "No finite values found - assuming beta scale."))
  if (length(v) > 20000) v <- v[sample.int(length(v), 20000)]
  q <- stats::quantile(v, c(0.001, 0.5, 0.999), na.rm = TRUE)
  if (q[1] >= -0.05 && q[3] <= 1.05) {
    list(scale = "beta", note = "Values fall within [0,1] - detected as beta values already.")
  } else if (q[1] >= -0.5 && q[3] <= 100.5 && q[3] > 1.5) {
    list(scale = "percent", note = sprintf(
      "Values range up to ~%.1f with a floor near 0 - consistent with percent methylation (0-100), not beta (0-1).", q[3]))
  } else {
    list(scale = "m", note = sprintf(
      "Values range from ~%.2f to ~%.2f (unbounded/negative) - consistent with M-values (log2 methylated/unmethylated ratio).", q[1], q[3]))
  }
}

methyl_ct_pct_to_beta <- function(mat) mat / 100

methyl_ct_m_to_beta <- function(m) 2^m / (1 + 2^m)

methyl_ct_working_summary <- function(mat, array_type = NULL) {
  finite_v <- mat[is.finite(mat)]
  rng <- if (length(finite_v) > 0) range(finite_v) else c(NA_real_, NA_real_)
  n_chr <- NA_integer_
  if (!is.null(array_type)) {
    anno <- methyl_get_annotation(array_type)
    if (isTRUE(anno$ok)) {
      hit <- rownames(mat)[rownames(mat) %in% rownames(anno$anno)]
      if (length(hit) > 0) n_chr <- length(unique(anno$anno[hit, "chr"]))
    }
  }
  list(n_cpg = nrow(mat), n_sample = ncol(mat), missing_pct = mean(is.na(mat)) * 100,
       beta_min = rng[1], beta_max = rng[2], n_duplicated = sum(duplicated(rownames(mat))), n_chr = n_chr)
}

methyl_ct_parse_custom_reference <- function(datapath, filename) {
  p <- methyl_parse_matrix(datapath, filename)
  if (!isTRUE(p$ok)) return(p)
  mat <- p$mat
  if (any(!is.finite(mat))) return(list(ok = FALSE, error = "The custom reference matrix contains non-finite values."))
  if (any(mat < -0.05 | mat > 1.05)) return(list(ok = FALSE, error = "The custom reference matrix must contain beta values (mean methylation per CpG x cell type) in [0,1]."))
  if (ncol(mat) < 2) return(list(ok = FALSE, error = "A reference matrix needs at least 2 cell-type columns."))
  list(ok = TRUE, mat = mat)
}

methyl_ct_marker_rank <- function(ref_mat) {
  ct <- colnames(ref_mat)
  n_ct <- length(ct)
  if (n_ct < 2) stop("A reference matrix needs at least 2 cell types to rank markers.")
  cpg <- rownames(ref_mat)
  n <- length(cpg)
  best_type <- character(n); effect <- rep(-Inf, n)
  direction <- character(n); specificity <- numeric(n); other_mean <- numeric(n)
  for (j in seq_len(n_ct)) {
    target <- ref_mat[, j]
    others <- ref_mat[, -j, drop = FALSE]
    om <- rowMeans(others)
    osd <- apply(others, 1, stats::sd)
    osd[!is.finite(osd) | osd == 0] <- 1e-6
    o_max <- apply(others, 1, max)
    o_min <- apply(others, 1, min)
    eff_hyper <- target - o_max
    eff_hypo <- o_min - target
    eff <- pmax(eff_hyper, eff_hypo)
    dir_j <- ifelse(eff_hyper >= eff_hypo, "hyper", "hypo")
    spec_j <- eff / osd
    better <- eff > effect
    best_type[better] <- ct[j]
    effect[better] <- eff[better]
    direction[better] <- dir_j[better]
    specificity[better] <- spec_j[better]
    other_mean[better] <- om[better]
  }
  btw_var <- apply(ref_mat, 1, stats::var)
  centroid_beta <- ref_mat[cbind(seq_len(n), match(best_type, ct))]
  data.frame(cpg = cpg, cell_type = best_type, effect = effect, direction = direction,
             specificity = specificity, other_mean = other_mean, btw_type_var = btw_var,
             centroid_beta = centroid_beta, row.names = NULL, stringsAsFactors = FALSE)
}

methyl_ct_select_markers <- function(rank_df, dbeta_min = 0, effect_min = 0, direction = "both",
                                      specificity_mode = "all", chr_allowed_ids = NULL) {
  df <- rank_df
  keep <- rep(TRUE, nrow(df))
  if (!is.na(dbeta_min) && dbeta_min > 0) keep <- keep & (df$effect >= dbeta_min)
  if (!is.na(effect_min) && effect_min > 0) keep <- keep & (df$effect >= effect_min)
  if (!identical(direction, "both")) keep <- keep & (df$direction == direction)
  if (!identical(specificity_mode, "all") && any(is.finite(df$specificity))) {
    med <- stats::median(df$specificity[is.finite(df$specificity)], na.rm = TRUE)
    if (identical(specificity_mode, "specific")) keep <- keep & (df$specificity >= med)
    else if (identical(specificity_mode, "shared")) keep <- keep & (df$specificity < med)
  }
  if (!is.null(chr_allowed_ids)) keep <- keep & (df$cpg %in% chr_allowed_ids)
  df[keep, , drop = FALSE]
}

methyl_ct_top_n_balanced <- function(df, sort_col = "effect", top_n = NULL) {
  if (is.null(top_n) || is.na(top_n) || top_n <= 0 || nrow(df) <= top_n) {
    return(df[order(-df[[sort_col]]), , drop = FALSE])
  }
  types <- unique(df$cell_type)
  per_type <- ceiling(top_n / length(types))
  picked <- do.call(rbind, lapply(types, function(t) {
    sub <- df[df$cell_type == t, , drop = FALSE]
    sub <- sub[order(-sub[[sort_col]]), , drop = FALSE]
    sub[seq_len(min(per_type, nrow(sub))), , drop = FALSE]
  }))
  picked <- picked[order(-picked[[sort_col]]), , drop = FALSE]
  if (nrow(picked) > top_n) picked <- picked[seq_len(top_n), , drop = FALSE]
  picked
}

methyl_ct_chr_allowed_ids <- function(cpg_ids, array_type, scope = c("all", "autosomes", "autosomes_x")) {
  scope <- match.arg(scope)
  if (identical(scope, "all")) return(list(ids = cpg_ids, note = "No chromosome restriction applied."))
  anno <- methyl_get_annotation(array_type)
  if (!isTRUE(anno$ok)) return(list(ids = cpg_ids, note = paste("Chromosome filtering unavailable:", anno$reason)))
  a <- anno$anno
  hit <- cpg_ids[cpg_ids %in% rownames(a)]
  chr <- a[hit, "chr"]
  drop_chr <- if (identical(scope, "autosomes")) c("chrX", "chrY", "X", "Y") else c("chrY", "Y")
  keep <- !(chr %in% drop_chr)
  allowed <- hit[keep]
  unresolved <- setdiff(cpg_ids, hit)
  list(ids = union(allowed, unresolved),
       note = sprintf("%d probe(s) removed by chromosome scope; %d unresolved probe(s) kept unfiltered.", sum(!keep), length(unresolved)))
}

methyl_ct_overlap_qc <- function(marker_ids, working_ids) {
  marker_ids <- unique(marker_ids)
  matched <- intersect(marker_ids, working_ids)
  missing <- setdiff(marker_ids, working_ids)
  pct <- if (length(marker_ids) > 0) length(matched) / length(marker_ids) * 100 else 0
  list(n_ref = length(marker_ids), n_matched = length(matched), n_missing = length(missing),
       pct_matched = pct, matched = matched, missing = missing)
}

methyl_ct_overlap_by_type <- function(marker_df, working_ids) {
  matched <- marker_df$cpg %in% working_ids
  types <- sort(unique(marker_df$cell_type))
  rows <- lapply(types, function(t) {
    sel <- marker_df$cell_type == t
    data.frame(cell_type = t, n_markers = sum(sel), n_matched = sum(sel & matched))
  })
  df <- do.call(rbind, rows)
  df$pct_matched <- ifelse(df$n_markers > 0, df$n_matched / df$n_markers * 100, NA_real_)
  df
}

methyl_ct_run_epidish <- function(beta_mat, ref_mat, method = c("RPC", "CBS", "CP"),
                                   maxit = 50, nu.v = c(0.25, 0.5, 0.75), constraint = c("inequality", "equality")) {
  method <- match.arg(method)
  constraint <- match.arg(constraint)
  if (!requireNamespace("EpiDISH", quietly = TRUE)) return(list(ok = FALSE, reason = "EpiDISH is not installed in this deployment."))
  common <- intersect(rownames(beta_mat), rownames(ref_mat))
  if (length(common) < 10) return(list(ok = FALSE, reason = sprintf(
    "Only %d marker CpG(s) overlap between the working matrix and the reference - too few for deconvolution (need at least 10).", length(common))))
  bm <- beta_mat[common, , drop = FALSE]
  rm_ <- ref_mat[common, , drop = FALSE]
  complete <- stats::complete.cases(bm)
  bm <- bm[complete, , drop = FALSE]
  rm_ <- rm_[rownames(bm), , drop = FALSE]
  if (nrow(bm) < 10) return(list(ok = FALSE, reason = "Fewer than 10 complete-case marker CpGs remain after removing missing values."))
  res <- tryCatch(
    EpiDISH::epidish(beta.m = bm, ref.m = rm_, method = method, maxit = maxit, nu.v = nu.v, constraint = constraint),
    error = function(e) e
  )
  if (inherits(res, "error")) return(list(ok = FALSE, reason = paste("epidish() failed:", conditionMessage(res))))
  list(ok = TRUE, fractions = res$estF, method = method, n_markers_used = nrow(bm), ref_used = rm_)
}

methyl_ct_run_hepidish <- function(beta_mat, ref1_mat, ref2_mat, ic_column,
                                    method = c("RPC", "CBS", "CP"), maxit = 50,
                                    nu.v = c(0.25, 0.5, 0.75), constraint = c("inequality", "equality")) {
  method <- match.arg(method)
  constraint <- match.arg(constraint)
  if (!requireNamespace("EpiDISH", quietly = TRUE)) return(list(ok = FALSE, reason = "EpiDISH is not installed in this deployment."))
  h_idx <- match(ic_column, colnames(ref1_mat))
  if (is.na(h_idx)) return(list(ok = FALSE, reason = sprintf("Column \"%s\" not found in the tissue-level reference.", ic_column)))
  common1 <- intersect(rownames(beta_mat), rownames(ref1_mat))
  common2 <- intersect(rownames(beta_mat), rownames(ref2_mat))
  if (length(common1) < 10 || length(common2) < 10) return(list(ok = FALSE, reason = sprintf(
    "Too few overlapping CpGs for the two-stage reference (stage 1: %d, stage 2: %d overlap with the working matrix; need >=10 each).",
    length(common1), length(common2))))
  res <- tryCatch(
    EpiDISH::hepidish(beta.m = beta_mat, ref1.m = ref1_mat, ref2.m = ref2_mat, h.CT.idx = h_idx,
                       method = method, maxit = maxit, nu.v = nu.v, constraint = constraint),
    error = function(e) e
  )
  if (inherits(res, "error")) return(list(ok = FALSE, reason = paste("hepidish() failed:", conditionMessage(res))))
  list(ok = TRUE, fractions = res, method = paste0("hepidish (", method, ")"),
       n_markers_used = length(common1) + length(common2))
}

methyl_ct_reconstruct <- function(ref_mat, fractions) {
  common_ct <- intersect(colnames(ref_mat), colnames(fractions))
  ref_mat[, common_ct, drop = FALSE] %*% t(fractions[, common_ct, drop = FALSE])
}

methyl_ct_validation_metrics <- function(observed, reconstructed) {
  common_cpg <- intersect(rownames(observed), rownames(reconstructed))
  common_sample <- intersect(colnames(observed), colnames(reconstructed))
  if (length(common_cpg) < 2 || length(common_sample) < 1) {
    return(list(ok = FALSE, reason = "Not enough overlapping CpGs/samples between observed and reconstructed data to validate."))
  }
  obs <- observed[common_cpg, common_sample, drop = FALSE]
  rec <- reconstructed[common_cpg, common_sample, drop = FALSE]
  keep <- stats::complete.cases(obs) & stats::complete.cases(rec)
  obs <- obs[keep, , drop = FALSE]
  rec <- rec[keep, , drop = FALSE]
  if (nrow(obs) < 2) return(list(ok = FALSE, reason = "Fewer than 2 complete-case CpGs remain for validation."))
  overall_cor <- stats::cor(as.vector(obs), as.vector(rec))
  overall_rmse <- sqrt(mean((obs - rec)^2))
  overall_mae <- mean(abs(obs - rec))
  overall_r2 <- 1 - sum((obs - rec)^2) / sum((obs - mean(obs))^2)
  per_sample <- data.frame(
    sample = colnames(obs),
    cor = vapply(seq_len(ncol(obs)), function(j) stats::cor(obs[, j], rec[, j]), numeric(1)),
    rmse = vapply(seq_len(ncol(obs)), function(j) sqrt(mean((obs[, j] - rec[, j])^2)), numeric(1)),
    mae = vapply(seq_len(ncol(obs)), function(j) mean(abs(obs[, j] - rec[, j])), numeric(1)),
    row.names = NULL
  )
  list(ok = TRUE, overall = list(cor = overall_cor, rmse = overall_rmse, mae = overall_mae, r2 = overall_r2),
       per_sample = per_sample, observed = obs, reconstructed = rec, n_cpg = nrow(obs), n_sample = ncol(obs))
}

methyl_ct_compare_methods <- function(beta_mat, ref_mat, methods = c("CP", "RPC", "CBS"), ...) {
  results <- list()
  failures <- character(0)
  for (m in methods) {
    r <- methyl_ct_run_epidish(beta_mat, ref_mat, method = m, ...)
    if (isTRUE(r$ok)) results[[m]] <- r$fractions else failures <- c(failures, sprintf("%s: %s", m, r$reason))
  }
  if (length(results) < 2) return(list(ok = FALSE, reason = paste(
    "Fewer than 2 methods produced results to compare.", paste(failures, collapse = "; "))))
  list(ok = TRUE, fractions_by_method = results, failures = failures)
}

methyl_ct_method_correlation <- function(fractions_by_method) {
  methods <- names(fractions_by_method)
  n <- length(methods)
  m <- matrix(1, n, n, dimnames = list(methods, methods))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) next
    a <- fractions_by_method[[i]]; b <- fractions_by_method[[j]]
    common_s <- intersect(rownames(a), rownames(b))
    common_ct <- intersect(colnames(a), colnames(b))
    m[i, j] <- stats::cor(as.vector(a[common_s, common_ct, drop = FALSE]), as.vector(b[common_s, common_ct, drop = FALSE]))
  }
  m
}

methyl_ct_method_agreement_summary <- function(fractions_by_method) {
  methods <- names(fractions_by_method)
  pairs <- utils::combn(methods, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- fractions_by_method[[p[1]]]; b <- fractions_by_method[[p[2]]]
    common_s <- intersect(rownames(a), rownames(b))
    common_ct <- intersect(colnames(a), colnames(b))
    diff <- abs(a[common_s, common_ct, drop = FALSE] - b[common_s, common_ct, drop = FALSE])
    data.frame(method_a = p[1], method_b = p[2], mean_abs_diff = mean(diff), max_abs_diff = max(diff))
  })
  do.call(rbind, rows)
}

methyl_ct_group_stats <- function(fractions, group) {
  group <- as.character(group)
  lv <- sort(unique(stats::na.omit(group)))
  if (length(lv) < 2) return(list(ok = FALSE, reason = "The grouping column needs at least two non-missing levels."))
  cts <- colnames(fractions)
  two_group <- length(lv) == 2
  rows <- lapply(cts, function(ct) {
    x <- fractions[, ct]
    if (two_group) {
      g1 <- x[group == lv[1]]; g2 <- x[group == lv[2]]
      wt <- tryCatch(stats::wilcox.test(g1, g2), error = function(e) NULL)
      n1 <- length(g1); n2 <- length(g2)
      eff <- if (!is.null(wt)) 1 - (2 * unname(wt$statistic)) / (n1 * n2) else NA_real_
      data.frame(cell_type = ct, test = "Wilcoxon rank-sum", statistic = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
                 p = if (!is.null(wt)) wt$p.value else NA_real_, effect_size = eff,
                 mean_diff = mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE))
    } else {
      kt <- tryCatch(stats::kruskal.test(x, factor(group)), error = function(e) NULL)
      n <- sum(!is.na(x))
      eff <- if (!is.null(kt)) unname((kt$statistic - (length(lv) - 1)) / (n - length(lv))) else NA_real_
      data.frame(cell_type = ct, test = "Kruskal-Wallis", statistic = if (!is.null(kt)) unname(kt$statistic) else NA_real_,
                 p = if (!is.null(kt)) kt$p.value else NA_real_, effect_size = eff, mean_diff = NA_real_)
    }
  })
  df <- do.call(rbind, rows)
  df$fdr <- stats::p.adjust(df$p, method = "BH")
  list(ok = TRUE, table = df, test_used = if (two_group) "Wilcoxon rank-sum" else "Kruskal-Wallis", levels = lv)
}

methyl_ct_composition_pca <- function(fractions, n_pcs = 10) {
  m <- stats::na.omit(fractions)
  if (nrow(m) < 3 || ncol(m) < 2) return(list(ok = FALSE, reason = "Not enough samples/cell types for PCA."))
  keep_cols <- apply(m, 2, function(x) stats::sd(x) > 0)
  if (sum(keep_cols) < 2) return(list(ok = FALSE, reason = "Cell-type fractions have no variance to ordinate."))
  pc <- stats::prcomp(m[, keep_cols, drop = FALSE], scale. = TRUE)
  k <- min(n_pcs, ncol(pc$x))
  var_explained <- (pc$sdev^2 / sum(pc$sdev^2))[seq_len(k)]
  list(ok = TRUE, scores = pc$x[, seq_len(k), drop = FALSE], var_explained = var_explained)
}

methyl_ct_composition_mds <- function(fractions, k = 2) {
  m <- stats::na.omit(fractions)
  if (nrow(m) < 3) return(list(ok = FALSE, reason = "Not enough samples for MDS."))
  d <- stats::dist(m)
  k_used <- min(k, nrow(m) - 1)
  fit <- tryCatch(stats::cmdscale(d, k = k_used), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE, reason = "MDS failed to converge on this distance matrix."))
  colnames(fit) <- paste0("Dim", seq_len(ncol(fit)))
  list(ok = TRUE, scores = fit)
}

methyl_ct_plot_marker_bar <- function(marker_df) {
  agg <- as.data.frame(table(cell_type = marker_df$cell_type))
  colnames(agg) <- c("cell_type", "n")
  ggplot(agg, aes(x = cell_type, y = n)) +
    geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.6) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.4, color = ARTHOMIX_COLORS$ink) +
    labs(x = NULL, y = "CpGs retained") + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

methyl_ct_plot_marker_heatmap <- function(ref_mat, marker_ids, max_rows = 200) {
  ids <- intersect(marker_ids, rownames(ref_mat))
  if (length(ids) > max_rows) ids <- ids[seq_len(max_rows)]
  sub <- ref_mat[ids, , drop = FALSE]
  df <- as.data.frame(as.table(sub))
  colnames(df) <- c("cpg", "cell_type", "beta")
  ggplot(df, aes(x = cell_type, y = cpg, fill = beta)) +
    geom_tile() +
    scale_fill_gradient(low = "#EAF3FB", high = ARTHOMIX_COLORS$blue, limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "Centroid beta") + theme_arthomix() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

methyl_ct_plot_marker_scatter <- function(marker_df) {
  pal <- arthomix_pair(marker_df$cell_type)
  ggplot(marker_df, aes(x = effect, y = specificity, color = cell_type)) +
    geom_point(alpha = 0.7, size = 1.6) +
    scale_color_manual(values = pal) +
    labs(x = "Effect size (methylation difference vs. other cell types)", y = "Specificity score", color = "Cell type") +
    theme_arthomix()
}

methyl_ct_plot_stacked_bar <- function(fractions, sample_order = NULL, hide_types = character(0)) {
  df <- as.data.frame(as.table(fractions))
  colnames(df) <- c("sample", "cell_type", "fraction")
  if (length(hide_types) > 0) df <- df[!df$cell_type %in% hide_types, , drop = FALSE]
  if (!is.null(sample_order)) df$sample <- factor(df$sample, levels = sample_order)
  pal <- arthomix_pair(unique(df$cell_type))
  ggplot(df, aes(x = sample, y = fraction, fill = cell_type)) +
    geom_col(position = "stack", width = 0.85) +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = "Estimated fraction", fill = "Cell type") +
    theme_arthomix() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
}

methyl_ct_plot_heatmap <- function(fractions, cluster_rows = TRUE, cluster_cols = TRUE, normalize = FALSE) {
  m <- t(fractions)
  if (isTRUE(normalize)) {
    rng <- apply(m, 1, function(r) diff(range(r)))
    rng[rng == 0] <- 1e-9
    m <- (m - apply(m, 1, min)) / rng
  }
  row_ord <- rownames(m); col_ord <- colnames(m)
  if (isTRUE(cluster_rows) && nrow(m) > 2) row_ord <- rownames(m)[stats::hclust(stats::dist(m))$order]
  if (isTRUE(cluster_cols) && ncol(m) > 2) col_ord <- colnames(m)[stats::hclust(stats::dist(t(m)))$order]
  df <- as.data.frame(as.table(m))
  colnames(df) <- c("cell_type", "sample", "value")
  df$cell_type <- factor(df$cell_type, levels = row_ord)
  df$sample <- factor(df$sample, levels = col_ord)
  ggplot(df, aes(x = sample, y = cell_type, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "#EAF3FB", high = ARTHOMIX_COLORS$blue) +
    labs(x = NULL, y = NULL, fill = if (normalize) "Row-normalized" else "Fraction") +
    theme_arthomix() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
}

methyl_ct_plot_box <- function(fractions, group = NULL, violin = FALSE) {
  df <- as.data.frame(as.table(fractions))
  colnames(df) <- c("sample", "cell_type", "fraction")
  if (!is.null(group)) {
    df$group <- group[as.character(df$sample)]
    p <- ggplot(df, aes(x = cell_type, y = fraction, fill = group))
  } else {
    p <- ggplot(df, aes(x = cell_type, y = fraction, fill = cell_type))
  }
  p <- p + (if (isTRUE(violin)) geom_violin(alpha = 0.8) else geom_boxplot(outlier.size = 0.5, linewidth = 0.3)) +
    labs(x = NULL, y = "Estimated fraction", fill = if (!is.null(group)) "Group" else "Cell type") +
    theme_arthomix() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p
}

methyl_ct_plot_scores <- function(scores, color_by = NULL, color_label = "Group", show_labels = FALSE,
                                   x_lab = "Dim 1", y_lab = "Dim 2") {
  df <- data.frame(sample = rownames(scores), x = scores[, 1], y = scores[, 2])
  if (!is.null(color_by)) {
    df$color <- color_by[df$sample]
    p <- ggplot(df, aes(x = x, y = y, color = color))
  } else {
    p <- ggplot(df, aes(x = x, y = y))
  }
  p <- p + geom_point(size = 2.4, alpha = 0.85, color = if (is.null(color_by)) ARTHOMIX_COLORS$blue else NULL) +
    labs(x = x_lab, y = y_lab, color = color_label) + theme_arthomix()
  if (isTRUE(show_labels)) p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 2.8, color = ARTHOMIX_COLORS$ink_secondary)
  p
}

methyl_ct_plot_corr <- function(cor_mat) {
  df <- as.data.frame(as.table(cor_mat))
  colnames(df) <- c("type1", "type2", "r")
  ggplot(df, aes(x = type1, y = type2, fill = r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3, color = ARTHOMIX_COLORS$ink) +
    scale_fill_gradient2(low = ARTHOMIX_STATUS$critical, mid = "white", high = ARTHOMIX_COLORS$blue, midpoint = 0, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "r") + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

methyl_ct_plot_reconstruction <- function(obs, rec, max_points = 20000) {
  ov <- as.vector(obs); rv <- as.vector(rec)
  if (length(ov) > max_points) {
    idx <- sample.int(length(ov), max_points)
    ov <- ov[idx]; rv <- rv[idx]
  }
  df <- data.frame(observed = ov, reconstructed = rv)
  ggplot(df, aes(x = observed, y = reconstructed)) +
    geom_point(alpha = 0.15, size = 0.6, color = ARTHOMIX_COLORS$blue) +
    geom_abline(slope = 1, intercept = 0, color = ARTHOMIX_STATUS$critical, linetype = "dashed") +
    labs(x = "Observed beta", y = "Reconstructed beta") + theme_arthomix()
}

methyl_ct_plot_bland_altman <- function(frac_a, frac_b, label_a, label_b) {
  common <- intersect(rownames(frac_a), rownames(frac_b))
  cts <- intersect(colnames(frac_a), colnames(frac_b))
  rows <- lapply(cts, function(ct) {
    a <- frac_a[common, ct]; b <- frac_b[common, ct]
    data.frame(sample = common, cell_type = ct, mean = (a + b) / 2, diff = a - b)
  })
  df <- do.call(rbind, rows)
  mean_diff <- mean(df$diff); sd_diff <- stats::sd(df$diff)
  ggplot(df, aes(x = mean, y = diff, color = cell_type)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_hline(yintercept = mean_diff, color = ARTHOMIX_COLORS$ink) +
    geom_hline(yintercept = mean_diff + 1.96 * sd_diff, linetype = "dashed", color = ARTHOMIX_STATUS$warning) +
    geom_hline(yintercept = mean_diff - 1.96 * sd_diff, linetype = "dashed", color = ARTHOMIX_STATUS$warning) +
    scale_color_manual(values = arthomix_pair(cts)) +
    labs(x = "Mean fraction", y = sprintf("%s − %s", label_a, label_b), color = "Cell type") +
    theme_arthomix()
}

methyl_ct_plot_group_diff <- function(fractions, group, stats_df) {
  df <- as.data.frame(as.table(fractions))
  colnames(df) <- c("sample", "cell_type", "fraction")
  df$group <- group[as.character(df$sample)]
  sig_label <- function(p) {
    if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  }
  ann <- stats_df
  ann$label <- vapply(ann$fdr, sig_label, character(1))
  ymax_by_type <- tapply(df$fraction, df$cell_type, max, na.rm = TRUE)
  ann$y <- unname(ymax_by_type[ann$cell_type]) * 1.08
  ggplot(df, aes(x = cell_type, y = fraction, fill = group)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
    geom_text(data = ann, aes(x = cell_type, y = y, label = label), inherit.aes = FALSE, size = 3.6) +
    scale_fill_manual(values = arthomix_pair(unique(stats::na.omit(df$group)))) +
    labs(x = NULL, y = "Estimated fraction", fill = "Group") +
    theme_arthomix() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

mod_methyl_celltype_dataqc_ui <- function(ns) {
  tagList(
    fluidRow(
      column(
        4,
        box(
          width = NULL, title = "1. Data source", status = "primary", solidHeader = FALSE,
          radioButtons(ns("ct_data_source"), NULL, width = "100%",
                       choices = c("Use the dataset loaded on the Dataset tab" = "shared",
                                   "Upload a different matrix here" = "own"),
                       selected = "shared"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'own'", ns("ct_data_source")),
            selectInput(ns("ct_own_array_type"), "Dataset type", choices = METHYL_ARRAY_TYPES, selected = "EPIC", width = "100%"),
            fileInput(ns("ct_own_matrix_file"), "Methylation matrix (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
            fileInput(ns("ct_own_sheet_file"), "Sample sheet / phenotype metadata (optional)", accept = c(".csv", ".tsv", ".txt")),
            uiOutput(ns("ct_own_preview_ui")),
            actionButton(ns("ct_own_load_btn"), "Load dataset", icon = icon("upload"), class = "btn-primary btn-sm")
          )
        )
      ),
      column(
        8,
        box(width = NULL, title = "2. Dataset preview", status = "primary", solidHeader = FALSE,
            uiOutput(ns("ct_source_summary_ui"))),
        uiOutput(ns("ct_scale_ui"))
      )
    ),
    box(
      width = 12, title = "3. QC filters", status = "primary", solidHeader = FALSE,
      fluidRow(
        column(3, selectInput(ns("ct_qc_missing_cpg"), "Missing CpG threshold",
                               choices = c("0%" = 0, "1%" = 0.01, "5%" = 0.05, "10%" = 0.1, "Custom" = "custom"), selected = 0.05),
               conditionalPanel(condition = sprintf("input['%s'] == 'custom'", ns("ct_qc_missing_cpg")),
                                 numericInput(ns("ct_qc_missing_cpg_custom"), NULL, value = 0.05, min = 0, max = 1, step = 0.01))),
        column(3, selectInput(ns("ct_qc_missing_sample"), "Missing sample threshold",
                               choices = c("0%" = 0, "5%" = 0.05, "10%" = 0.1, "20%" = 0.2, "Custom" = "custom"), selected = 0.1),
               conditionalPanel(condition = sprintf("input['%s'] == 'custom'", ns("ct_qc_missing_sample")),
                                 numericInput(ns("ct_qc_missing_sample_custom"), NULL, value = 0.1, min = 0, max = 1, step = 0.01))),
        column(3, selectInput(ns("ct_qc_sexchr"), "Chromosome scope",
                               choices = c("Autosomes only (remove chrX/chrY) - default" = "remove_xy",
                                           "Autosomes + chrX (remove chrY only)" = "remove_y_only",
                                           "Keep all chromosomes" = "keep"), selected = "remove_xy")),
        column(3, checkboxInput(ns("ct_qc_beta_range"), "Validate beta-value range (0-1)", value = TRUE))
      ),
      fluidRow(
        column(3, checkboxInput(ns("ct_qc_snp"), "Remove SNP-associated CpGs", value = FALSE)),
        column(3, checkboxInput(ns("ct_qc_crossreactive"), "Remove cross-reactive probes", value = FALSE)),
        column(3, conditionalPanel(condition = sprintf("input['%s'] == true", ns("ct_qc_crossreactive")),
                                     fileInput(ns("ct_qc_crossreactive_file"), "Exclusion list (probe IDs)", accept = c(".csv", ".tsv", ".txt")))),
        column(3, uiOutput(ns("ct_qc_detp_ui")))
      ),
      uiOutput(ns("ct_qc_cascade_ui"))
    )
  )
}

mod_methyl_celltype_featsel_ui <- function(ns) {
  tagList(
    box(
      width = 12, title = "CpG marker selection", status = "primary", solidHeader = FALSE,
      fluidRow(
        column(3, selectInput(ns("ct_fs_method"), "Method", choices = c(
          "Reference-library markers" = "reference",
          "Variance-based CpGs" = "variance",
          "Custom CpG list" = "custom",
          "Differential methylation markers / DMCs / DMRs (unavailable)" = "dmc_unavailable"
        ), selected = "reference")),
        column(3, uiOutput(ns("ct_fs_dbeta_ui"))),
        column(3, selectInput(ns("ct_fs_topn"), "Number of CpGs",
                               choices = c("50" = 50, "100" = 100, "200" = 200, "333" = 333, "500" = 500, "1000" = 1000, "2000" = 2000, "Custom" = "custom"),
                               selected = 200)),
        column(3, conditionalPanel(condition = sprintf("input['%s'] == 'custom'", ns("ct_fs_topn")),
                                     numericInput(ns("ct_fs_topn_custom"), "Custom count", value = 200, min = 1)))
      ),
      fluidRow(
        column(3, selectInput(ns("ct_fs_direction"), "Marker direction",
                               choices = c("Both" = "both", "Hyper-methylated" = "hyper", "Hypo-methylated" = "hypo"), selected = "both")),
        column(3, selectInput(ns("ct_fs_specificity"), "Specificity",
                               choices = c("All available" = "all", "Cell-type-specific" = "specific", "Shared" = "shared"), selected = "all")),
        column(3, selectInput(ns("ct_fs_chr_scope"), "Chromosome filtering",
                               choices = c("Autosomes only" = "autosomes", "Autosomes + X" = "autosomes_x", "All chromosomes" = "all"), selected = "all")),
        column(3, conditionalPanel(condition = sprintf("input['%s'] == 'custom'", ns("ct_fs_method")),
                                     fileInput(ns("ct_fs_custom_file"), "Custom CpG list (one ID per line, or CSV)", accept = c(".csv", ".tsv", ".txt"))))
      ),
      p(class = "empty-note", icon("circle-info"),
        "Maximum FDR: not applicable - reference centroids have no per-sample replicates, so no p-value/FDR can be computed here."),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'dmc_unavailable'", ns("ct_fs_method")),
        div(class = "empty-note", icon("ban"),
            "Differential methylation markers/DMCs/DMRs need per-sample sorted-cell-type data, which isn't available for the built-in reference panels (mean-beta centroids only, no replicates) - not built in this pass.")
      ),
      div(style = "margin-top:8px;",
          actionButton(ns("ct_fs_run_btn"), "Run CpG Feature Selection", icon = icon("play"), class = "btn-primary"))
    ),
    uiOutput(ns("ct_fs_result_gate"))
  )
}

mod_methyl_celltype_refmethod_ui <- function(ns) {
  reg <- methyl_ct_reference_registry()
  ref_choices <- setNames(vapply(reg, `[[`, character(1), "id"), vapply(reg, `[[`, character(1), "label"))
  blood_refs <- Filter(function(s) identical(s$tissue, "blood"), reg)
  blood_choices <- setNames(vapply(blood_refs, `[[`, character(1), "id"), vapply(blood_refs, `[[`, character(1), "label"))
  tagList(
    fluidRow(
      column(
        5,
        box(
          width = NULL, title = "Reference library", status = "primary", solidHeader = FALSE,
          radioButtons(ns("ct_ref_source"), NULL,
                       choices = c("Built-in reference library" = "registry", "Custom reference matrix (upload)" = "custom"),
                       selected = "registry"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'registry'", ns("ct_ref_source")),
            selectInput(ns("ct_ref_id"), "Reference", choices = ref_choices, selected = "blood7", width = "100%"),
            uiOutput(ns("ct_ref_celltypes_ui"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'custom'", ns("ct_ref_source")),
            p(class = "empty-note", icon("circle-info"),
              "CpG x cell-type mean-beta CSV/TSV - CpGs in rows, one column per cell type, first column the CpG ID."),
            fileInput(ns("ct_custom_ref_file"), "Custom reference matrix", accept = c(".csv", ".tsv", ".txt")),
            uiOutput(ns("ct_custom_ref_preview_ui"))
          )
        )
      ),
      column(
        7,
        box(
          width = NULL, title = "Deconvolution method", status = "primary", solidHeader = FALSE,
          selectInput(ns("ct_method"), NULL, width = "100%",
                      choices = c("CP - Houseman constrained projection" = "CP",
                                  "EpiDISH RPC (robust partial correlations)" = "RPC",
                                  "EpiDISH CBS (CIBERSORT-style)" = "CBS",
                                  "Two-stage (hepidish) - advanced" = "hepidish")),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'hepidish'", ns("ct_method")),
            p(class = "empty-note", icon("circle-info"),
              "Two-stage: the reference above splits the tissue into its top-level components (e.g. Epithelial/Fibroblast/Immune); a second blood reference then further decomposes the chosen component into its own cell subtypes."),
            selectInput(ns("ct_hepidish_ic_col"), "Column to further decompose", choices = NULL, width = "100%"),
            selectInput(ns("ct_hepidish_ref2"), "Second-stage (immune subtype) reference",
                        choices = blood_choices, selected = "blood7", width = "100%")
          ),
          div(class = "empty-note", icon("ban"), "Not available in this deployment:",
              tags$ul(lapply(methyl_ct_unavailable_methods(), function(m) tags$li(strong(m$label), " - ", m$reason))))
        ),
        box(
          width = NULL, title = "Reference-library QC", status = "primary", solidHeader = FALSE,
          selectInput(ns("ct_overlap_threshold"), "Minimum marker overlap to allow deconvolution",
                      choices = c("50%" = 0.5, "60%" = 0.6, "70%" = 0.7, "80%" = 0.8, "90%" = 0.9), selected = 0.5, width = "100%"),
          uiOutput(ns("ct_refqc_ui"))
        ),
        tags$details(
          class = "box box-primary",
          tags$summary(class = "box-header", tags$h3(class = "box-title", icon("sliders"), " Advanced Parameters")),
          div(
            class = "box-body",
            fluidRow(
              column(4, numericInput(ns("ct_adv_maxit"), "Max iterations", value = 50, min = 10, max = 500)),
              column(4, selectInput(ns("ct_adv_constraint"), "Constraint",
                                     choices = c("Inequality (non-negative)" = "inequality", "Equality (sum-to-one)" = "equality"),
                                     selected = "inequality")),
              column(4, numericInput(ns("ct_adv_seed"), "Random seed", value = 1234))
            ),
            p(strong("CBS tuning (nu.v)")),
            fluidRow(
              column(4, numericInput(ns("ct_adv_nu1"), NULL, value = 0.25, min = 0, max = 1, step = 0.05)),
              column(4, numericInput(ns("ct_adv_nu2"), NULL, value = 0.5, min = 0, max = 1, step = 0.05)),
              column(4, numericInput(ns("ct_adv_nu3"), NULL, value = 0.75, min = 0, max = 1, step = 0.05))
            )
          )
        )
      )
    )
  )
}

mod_methyl_celltype_deconv_ui <- function(ns) {
  tagList(
    box(
      width = 12, status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc", "Runs only when clicked - no analysis is triggered automatically by uploading data, selecting a reference, or changing a filter."),
      actionButton(ns("ct_run_decon_btn"), "Run Cell-Type Deconvolution", icon = icon("play"), class = "btn-primary btn-lg")
    ),
    uiOutput(ns("ct_decon_result_gate"))
  )
}

mod_methyl_celltype_composition_ui <- function(ns) {
  uiOutput(ns("ct_composition_gate"))
}

mod_methyl_celltype_validation_ui <- function(ns) {
  tagList(
    box(
      width = 12, status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc", "Compares observed methylation to the methylation reconstructed from the estimated cell fractions and the reference used."),
      actionButton(ns("ct_run_val_btn"), "Run Validation", icon = icon("play"), class = "btn-primary")
    ),
    uiOutput(ns("ct_val_result_gate")),
    box(
      width = 12, title = "Compare Methods", status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc", "Re-runs deconvolution with each selected method on the identical filtered matrix/reference and compares the estimated fractions."),
      checkboxGroupInput(ns("ct_cmpmethods_pick"), "Methods to compare", choices = c("CP", "RPC", "CBS"), selected = c("CP", "RPC", "CBS"), inline = TRUE),
      actionButton(ns("ct_run_cmpmethods_btn"), "Compare Methods", icon = icon("play"), class = "btn-primary")
    ),
    uiOutput(ns("ct_cmpmethods_result_gate"))
  )
}

mod_methyl_celltype_export_ui <- function(ns) {
  tagList(
    box(
      width = 12, title = "Data", status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc", "Each download uses whichever step has already been run - unavailable items are disabled. Figure downloads are available directly beneath each figure in its own tab."),
      fluidRow(
        column(3, downloadButton(ns("ct_export_beta"), "Filtered beta matrix", class = "btn-sm")),
        column(3, downloadButton(ns("ct_export_markers"), "Selected CpG list", class = "btn-sm")),
        column(3, downloadButton(ns("ct_export_ref"), "Reference matrix used", class = "btn-sm")),
        column(3, downloadButton(ns("ct_export_fractions"), "Estimated cell fractions", class = "btn-sm"))
      ),
      fluidRow(
        column(3, downloadButton(ns("ct_export_pheno_fractions"), "Phenotype-linked fractions", class = "btn-sm")),
        column(3, downloadButton(ns("ct_export_comparison"), "Group comparison results", class = "btn-sm"))
      )
    ),
    box(width = 12, title = "Report", status = "primary", solidHeader = FALSE,
        downloadButton(ns("ct_export_report"), "Download analysis summary (.txt)", class = "btn-primary btn-sm"))
  )
}

mod_methyl_celltype_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "submodule-desc", mod_methyl_celltype_config$description),
    tabsetPanel(
      id = ns("ct_subtabs"), type = "tabs",
      tabPanel("1. Data & QC", value = "dataqc", mod_methyl_celltype_dataqc_ui(ns)),
      tabPanel("2. CpG Feature Selection", value = "featsel", mod_methyl_celltype_featsel_ui(ns)),
      tabPanel("3. Reference & Method", value = "refmethod", mod_methyl_celltype_refmethod_ui(ns)),
      tabPanel("4. Deconvolution", value = "deconv", mod_methyl_celltype_deconv_ui(ns)),
      tabPanel("5. Cell Composition", value = "composition", mod_methyl_celltype_composition_ui(ns)),
      tabPanel("6. Validation", value = "validation", mod_methyl_celltype_validation_ui(ns)),
      tabPanel("7. Export", value = "export", mod_methyl_celltype_export_ui(ns))
    )
  )
}

mod_methyl_celltype_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    register_has_run_gate_local <- function(gate_id, has_run_flag_fn, result_output_id, not_run_message) {
      output[[gate_id]] <- renderUI({
        if (isTRUE(has_run_flag_fn())) uiOutput(ns(result_output_id))
        else div(class = "card", p(class = "empty-note", icon("circle-info"), not_run_message))
      })
    }

    make_plot_dl <- function(build_fn, base_name) downloadHandler(
      filename = function() sprintf("%s.png", base_name),
      content = function(file) ggplot2::ggsave(file, plot = build_fn(), width = 9, height = 6, dpi = 300, device = "png")
    )

    plotly_safe <- function(p) {
      plotly::layout(plotly::ggplotly(p),
                      legend = list(orientation = "h", x = 0, y = -0.3, yanchor = "top"),
                      margin = list(b = 110))
    }

    own_raw <- reactiveVal(NULL)
    own_ready <- reactiveVal(NULL)

    own_matrix_parsed <- reactive({
      req(input$ct_own_matrix_file)
      methyl_parse_matrix(input$ct_own_matrix_file$datapath, input$ct_own_matrix_file$name)
    })
    own_sheet_parsed <- reactive({
      req(input$ct_own_sheet_file)
      methyl_parse_sample_sheet(input$ct_own_sheet_file$datapath, input$ct_own_sheet_file$name)
    })
    output$ct_own_preview_ui <- renderUI({
      req(input$ct_own_matrix_file)
      p <- own_matrix_parsed()
      if (!isTRUE(p$ok)) return(div(class = "empty-note", icon("triangle-exclamation"), p$error))
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s CpGs x %s samples.", input$ct_own_matrix_file$name, format(nrow(p$mat), big.mark = ","), ncol(p$mat)))
    })
    observeEvent(input$ct_own_load_btn, {
      p <- own_matrix_parsed()
      if (!isTRUE(p$ok)) { showNotification(p$error, type = "error"); return() }
      sheet_df <- NULL
      if (!is.null(input$ct_own_sheet_file)) {
        s <- own_sheet_parsed()
        if (isTRUE(s$ok)) sheet_df <- s$df
      }
      raw <- list(mat = p$mat, sheet = sheet_df, array_type = input$ct_own_array_type,
                  source_label = sprintf("Uploaded (this module): %s", input$ct_own_matrix_file$name))
      own_raw(raw)
      det <- methyl_ct_detect_scale(raw$mat)
      own_ready(if (identical(det$scale, "beta")) raw else NULL)
    })

    ct_scale_detect <- reactive({
      raw <- own_raw()
      if (!identical(input$ct_data_source, "own") || is.null(raw)) return(NULL)
      methyl_ct_detect_scale(raw$mat)
    })

    output$ct_scale_ui <- renderUI({
      if (!identical(input$ct_data_source, "own")) return(NULL)
      raw <- own_raw()
      req(raw)
      if (!is.null(own_ready())) return(div(class = "empty-note", icon("check"), "Working matrix is on the beta scale."))
      det <- ct_scale_detect()
      box(width = NULL, title = "Methylation scale", status = "warning", solidHeader = FALSE,
          p(class = "empty-note", icon("triangle-exclamation"), det$note),
          actionButton(ns("ct_apply_transform_btn"),
                       sprintf("Apply Transformation (%s -> Beta)", if (identical(det$scale, "percent")) "%" else "M-value"),
                       class = "btn-warning btn-sm"))
    })
    observeEvent(input$ct_apply_transform_btn, {
      raw <- own_raw()
      req(raw)
      det <- ct_scale_detect()
      mat2 <- if (identical(det$scale, "percent")) methyl_ct_pct_to_beta(raw$mat) else methyl_ct_m_to_beta(raw$mat)
      own_ready(list(mat = mat2, sheet = raw$sheet, array_type = raw$array_type,
                      source_label = paste(raw$source_label, sprintf("(%s -> beta applied)", det$scale))))
    })

    ct_source <- reactive({
      if (identical(input$ct_data_source, "own")) {
        r <- own_ready()
        req(r)
        list(mat = r$mat, sheet = r$sheet, array_type = r$array_type, source_label = r$source_label, detp = NULL)
      } else {
        req(dataset$beta)
        mat <- dataset$beta
        if (identical(dataset$input_scale, "m")) mat <- methyl_ct_m_to_beta(mat)
        list(mat = mat, sheet = dataset$sample_sheet, array_type = dataset$array_type %||% "450K",
             source_label = dataset$source %||% "Dataset tab", detp = dataset$detp)
      }
    })

    output$ct_source_summary_ui <- renderUI({
      src <- tryCatch(ct_source(), error = function(e) NULL)
      if (is.null(src)) return(div(class = "empty-note", icon("circle-info"), "No dataset loaded yet - use the Dataset tab, or upload one here."))
      s <- methyl_ct_working_summary(src$mat, src$array_type)
      tagList(
        p(strong(src$source_label)),
        tags$ul(
          tags$li(sprintf("%s CpGs x %s samples", format(s$n_cpg, big.mark = ","), s$n_sample)),
          tags$li(sprintf("Missing values: %.2f%%", s$missing_pct)),
          tags$li(sprintf("Beta range: [%.3f, %.3f]", s$beta_min, s$beta_max)),
          tags$li(sprintf("Duplicated probe IDs: %d", s$n_duplicated)),
          tags$li(sprintf("Chromosomes represented: %s", if (is.na(s$n_chr)) "unknown (no manifest annotation for this array type)" else s$n_chr))
        )
      )
    })

    output$ct_qc_detp_ui <- renderUI({
      src <- tryCatch(ct_source(), error = function(e) NULL)
      if (!is.null(src) && !is.null(src$detp)) {
        checkboxInput(ns("ct_qc_detp"), "Apply detection p-value filtering", value = FALSE)
      } else {
        div(class = "empty-note", icon("circle-info"), "Detection p-values require raw IDAT input - not available for this dataset.")
      }
    })

    ct_filtered <- reactive({
      src <- ct_source()
      mat <- src$mat
      cascade <- list()
      anno <- methyl_get_annotation(src$array_type)

      miss_cpg <- if (identical(input$ct_qc_missing_cpg, "custom")) input$ct_qc_missing_cpg_custom else as.numeric(input$ct_qc_missing_cpg %||% 0.05)
      f1 <- methyl_filter_missing(mat, max_na_frac = miss_cpg); cascade[["Missing-CpG filter"]] <- f1
      mat <- mat[f1$keep, , drop = FALSE]

      f2 <- methyl_filter_sex_chr(mat, anno, mode = input$ct_qc_sexchr %||% "remove_xy"); cascade[["Chromosome scope"]] <- f2
      mat <- mat[f2$keep, , drop = FALSE]

      if (isTRUE(input$ct_qc_snp)) {
        f3 <- methyl_filter_snp(mat, anno); cascade[["SNP-associated"]] <- f3
        mat <- mat[f3$keep, , drop = FALSE]
      }
      if (isTRUE(input$ct_qc_crossreactive)) {
        excl <- NULL
        if (!is.null(input$ct_qc_crossreactive_file)) {
          pl <- methyl_parse_probe_list(input$ct_qc_crossreactive_file$datapath, input$ct_qc_crossreactive_file$name)
          if (isTRUE(pl$ok)) excl <- pl$ids
        }
        f4 <- methyl_filter_cross_reactive(mat, excl); cascade[["Cross-reactive"]] <- f4
        mat <- mat[f4$keep, , drop = FALSE]
      }
      if (isTRUE(input$ct_qc_detp) && !is.null(src$detp)) {
        f5 <- methyl_filter_detection_p(mat, src$detp); cascade[["Detection p-value"]] <- f5
        mat <- mat[f5$keep, , drop = FALSE]
      }

      miss_sample <- if (identical(input$ct_qc_missing_sample, "custom")) input$ct_qc_missing_sample_custom else as.numeric(input$ct_qc_missing_sample %||% 0.1)
      fs <- methyl_fs_sample_missing_ok(mat, max_na_frac = miss_sample)
      mat <- mat[, fs$keep, drop = FALSE]

      cascade_df <- methyl_probe_retention_cascade(nrow(src$mat), cascade)
      list(mat = mat, cascade_df = cascade_df, array_type = src$array_type, sheet = src$sheet, source_label = src$source_label)
    })

    output$ct_qc_cascade_ui <- renderUI({
      f <- tryCatch(ct_filtered(), error = function(e) NULL)
      if (is.null(f)) return(NULL)
      tagList(
        p(class = "empty-note", icon("filter"),
          sprintf("Working matrix after filters: %s CpGs x %s samples.", format(nrow(f$mat), big.mark = ","), ncol(f$mat))),
        DT::dataTableOutput(ns("ct_qc_cascade_table"))
      )
    })
    output$ct_qc_cascade_table <- DT::renderDataTable({
      f <- ct_filtered()
      req(f)
      DT::datatable(f$cascade_df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10, dom = "t"), class = "stripe hover compact")
    })
    outputOptions(output, "ct_qc_cascade_table", suspendWhenHidden = FALSE)

    ct_registry <- methyl_ct_reference_registry()

    output$ct_ref_celltypes_ui <- renderUI({
      ref <- tryCatch(methyl_ct_get_reference(input$ct_ref_id), error = function(e) NULL)
      req(ref)
      checkboxGroupInput(ns("ct_ref_celltypes_keep"), "Cell types to include", choices = colnames(ref), selected = colnames(ref), inline = TRUE)
    })

    active_reference_full <- reactive({
      if (identical(input$ct_ref_source, "custom")) {
        req(ct_custom_ref())
        ct_custom_ref()$mat
      } else {
        req(input$ct_ref_id)
        methyl_ct_get_reference(input$ct_ref_id)
      }
    })

    active_reference <- reactive({
      ref <- active_reference_full()
      req(ref)
      if (identical(input$ct_ref_source, "registry") && !is.null(input$ct_ref_celltypes_keep)) {
        keep <- intersect(input$ct_ref_celltypes_keep, colnames(ref))
        validate(need(length(keep) >= 2, "Select at least 2 cell types."))
        ref <- ref[, keep, drop = FALSE]
      }
      ref
    })

    ct_custom_ref_raw <- reactive({
      req(input$ct_custom_ref_file)
      methyl_ct_parse_custom_reference(input$ct_custom_ref_file$datapath, input$ct_custom_ref_file$name)
    })
    output$ct_custom_ref_preview_ui <- renderUI({
      req(input$ct_custom_ref_file)
      p <- ct_custom_ref_raw()
      if (!isTRUE(p$ok)) return(div(class = "empty-note", icon("triangle-exclamation"), p$error))
      div(class = "empty-note", icon("check"),
          sprintf("%s CpGs x %s cell types: %s", format(nrow(p$mat), big.mark = ","), ncol(p$mat), paste(colnames(p$mat), collapse = ", ")))
    })
    ct_custom_ref <- reactive({
      p <- ct_custom_ref_raw()
      req(isTRUE(p$ok))
      list(mat = p$mat)
    })

    observeEvent(input$ct_method, {
      req(identical(input$ct_method, "hepidish"))
      ref <- tryCatch(active_reference_full(), error = function(e) NULL)
      if (!is.null(ref)) {
        updateSelectInput(session, "ct_hepidish_ic_col", choices = colnames(ref),
                           selected = if ("IC" %in% colnames(ref)) "IC" else colnames(ref)[length(colnames(ref))])
      }
    })

    output$ct_fs_dbeta_ui <- renderUI({
      if (identical(input$ct_fs_method, "reference")) {
        selectInput(ns("ct_fs_dbeta"), "Min |methylation difference|",
                    choices = c("0.05" = 0.05, "0.10" = 0.10, "0.15" = 0.15, "0.20" = 0.20, "Custom" = "custom"), selected = 0.10)
      } else {
        div(class = "empty-note", icon("circle-info"), "Not applicable to this method.")
      }
    })

    ct_fs_custom_ids <- reactive({
      req(input$ct_fs_custom_file)
      pl <- methyl_parse_probe_list(input$ct_fs_custom_file$datapath, input$ct_fs_custom_file$name)
      if (isTRUE(pl$ok)) pl$ids else character(0)
    })

    fs_has_run <- reactiveVal(FALSE)
    register_has_run_gate_local("ct_fs_result_gate", fs_has_run, "ct_fs_result_ui",
                                 "Not run yet - set the filters above and click \"Run CpG Feature Selection\".")

    fs_result <- eventReactive(input$ct_fs_run_btn, {
      ref <- active_reference()
      validate(need(!is.null(ref), "Select a reference library first."))
      rank_df <- methyl_ct_marker_rank(ref)

      chr_ids <- NULL
      if (!identical(input$ct_fs_chr_scope, "all")) {
        src <- ct_source()
        cs <- methyl_ct_chr_allowed_ids(rank_df$cpg, src$array_type, scope = input$ct_fs_chr_scope)
        chr_ids <- cs$ids
      }

      validate(need(!identical(input$ct_fs_method, "dmc_unavailable"),
                    "Differential methylation markers/DMCs/DMRs are not available (see note above) - choose a different method."))
      if (identical(input$ct_fs_method, "custom")) {
        ids <- intersect(ct_fs_custom_ids(), rank_df$cpg)
        validate(need(length(ids) > 0, "None of the uploaded CpG IDs are present in the selected reference."))
        sel <- rank_df[rank_df$cpg %in% ids, , drop = FALSE]
      } else {
        dbeta_min <- if (identical(input$ct_fs_method, "reference")) {
          if (identical(input$ct_fs_dbeta, "custom")) 0 else as.numeric(input$ct_fs_dbeta %||% 0)
        } else 0
        sel <- methyl_ct_select_markers(rank_df, dbeta_min = dbeta_min, direction = input$ct_fs_direction %||% "both",
                                         specificity_mode = input$ct_fs_specificity %||% "all", chr_allowed_ids = chr_ids)
        sort_col <- if (identical(input$ct_fs_method, "variance")) "btw_type_var" else "effect"
        top_n <- if (identical(input$ct_fs_topn, "custom")) input$ct_fs_topn_custom else as.numeric(input$ct_fs_topn %||% 200)
        sel <- methyl_ct_top_n_balanced(sel, sort_col = sort_col, top_n = top_n)
      }
      validate(need(nrow(sel) > 0, "No CpGs passed the current filters - relax the thresholds."))
      list(ranked = rank_df, selected = sel, ref = ref, method = input$ct_fs_method)
    })
    observeEvent(fs_result(), { fs_has_run(TRUE) })

    output$ct_fs_result_ui <- renderUI({
      tagList(
        fluidRow(
          column(6, withSpinner(plotOutput(ns("ct_fs_bar_plot"), height = "320px"), color = "#2563EB", type = 6)),
          column(6, withSpinner(plotOutput(ns("ct_fs_heatmap_plot"), height = "320px"), color = "#2563EB", type = 6))
        ),
        withSpinner(plotly::plotlyOutput(ns("ct_fs_scatter_plot"), height = "440px"), color = "#2563EB", type = 6),
        div(class = "table-toolbar", downloadButton(ns("ct_fs_download"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("ct_fs_table"))
      )
    })
    output$ct_fs_bar_plot <- renderPlot({ req(fs_result()); methyl_ct_plot_marker_bar(fs_result()$selected) })
    output$ct_fs_heatmap_plot <- renderPlot({ req(fs_result()); methyl_ct_plot_marker_heatmap(fs_result()$ref, fs_result()$selected$cpg) })
    output$ct_fs_scatter_plot <- plotly::renderPlotly({ req(fs_result()); plotly_safe(methyl_ct_plot_marker_scatter(fs_result()$selected)) })
    output$ct_fs_table <- DT::renderDataTable({
      req(fs_result())
      df <- fs_result()$selected[, c("cpg", "cell_type", "effect", "direction", "specificity")]
      df$p <- "n/a"; df$fdr <- "n/a (no replicates)"
      colnames(df) <- c("CpG", "Cell Type", "Effect Size", "Direction", "Specificity", "P-value", "FDR")
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("Effect Size", "Specificity"), digits = 4)
    })
    outputOptions(output, "ct_fs_table", suspendWhenHidden = FALSE)
    output$ct_fs_download <- downloadHandler(
      filename = function() "celltype_selected_markers.csv",
      content = function(file) utils::write.csv(fs_result()$selected, file, row.names = FALSE)
    )

    ct_active_markers <- reactive({
      if (isTRUE(fs_has_run())) fs_result()$selected$cpg
      else {
        ref <- tryCatch(active_reference(), error = function(e) NULL)
        if (is.null(ref)) character(0) else rownames(ref)
      }
    })
    ct_active_marker_df <- reactive({
      if (isTRUE(fs_has_run())) fs_result()$selected
      else {
        ref <- tryCatch(active_reference(), error = function(e) NULL)
        if (is.null(ref)) NULL else methyl_ct_marker_rank(ref)
      }
    })

    ct_overlap <- reactive({
      f <- tryCatch(ct_filtered(), error = function(e) NULL)
      ids <- ct_active_markers()
      if (is.null(f) || length(ids) == 0) return(NULL)
      list(overall = methyl_ct_overlap_qc(ids, rownames(f$mat)),
           by_type = { df <- ct_active_marker_df(); if (!is.null(df)) methyl_ct_overlap_by_type(df, rownames(f$mat)) else NULL })
    })

    ct_overlap_ok <- reactive({
      ov <- ct_overlap()
      req(ov)
      thresh <- as.numeric(input$ct_overlap_threshold %||% 0.5)
      (ov$overall$pct_matched / 100) >= thresh
    })

    output$ct_refqc_ui <- renderUI({
      ov <- ct_overlap()
      if (is.null(ov)) return(div(class = "empty-note", icon("circle-info"), "Load a dataset and select a reference to see overlap QC."))
      thresh <- as.numeric(input$ct_overlap_threshold %||% 0.5) * 100
      below <- ov$overall$pct_matched < thresh
      status_icon <- if (!below) icon("circle-check") else icon("triangle-exclamation")
      msg <- sprintf("%.0f%% of reference marker CpGs (%d of %d) were found in the working dataset.",
                      ov$overall$pct_matched, ov$overall$n_matched, ov$overall$n_ref)
      tagList(
        p(class = "empty-note", status_icon, msg,
          if (below) " Deconvolution estimates may be unreliable - Run Deconvolution is disabled until overlap improves." else ""),
        if (!is.null(ov$by_type)) DT::dataTableOutput(ns("ct_refqc_by_type_table"))
      )
    })
    output$ct_refqc_by_type_table <- DT::renderDataTable({
      ov <- ct_overlap()
      req(ov$by_type)
      DT::datatable(ov$by_type, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10, dom = "t"), class = "stripe hover compact") %>%
        DT::formatSignif(columns = "pct_matched", digits = 3)
    })
    outputOptions(output, "ct_refqc_by_type_table", suspendWhenHidden = FALSE)

    observe({
      ok <- tryCatch(isTRUE(ct_overlap_ok()), error = function(e) FALSE)
      if (isTRUE(ok)) shinyjs::enable("ct_run_decon_btn") else shinyjs::disable("ct_run_decon_btn")
    })

    decon_has_run <- reactiveVal(FALSE)
    register_has_run_gate_local("ct_decon_result_gate", decon_has_run, "ct_decon_result_ui",
                                 "Not run yet - click \"Run Cell-Type Deconvolution\" above.")

    decon_result <- eventReactive(input$ct_run_decon_btn, {
      f <- ct_filtered()
      validate(need(!is.null(f) && nrow(f$mat) > 0, "No working matrix available - check the Data & QC tab."))
      ref <- active_reference()
      validate(need(!is.null(ref), "Select a reference library on the Reference & Method tab."))
      markers <- ct_active_markers()
      ref_use <- if (length(markers) > 0) ref[intersect(rownames(ref), markers), , drop = FALSE] else ref
      validate(need(nrow(ref_use) >= 10, "Fewer than 10 marker CpGs available for deconvolution."))

      nu <- c(input$ct_adv_nu1 %||% 0.25, input$ct_adv_nu2 %||% 0.5, input$ct_adv_nu3 %||% 0.75)
      set.seed(input$ct_adv_seed %||% 1234)

      if (identical(input$ct_method, "hepidish")) {
        validate(need(!is.null(input$ct_hepidish_ref2), "Select a second-stage reference for the two-stage method."))
        ref2 <- methyl_ct_get_reference(input$ct_hepidish_ref2)
        validate(need(!is.null(ref2), "Could not load the second-stage reference."))
        res <- methyl_ct_run_hepidish(f$mat, ref_use, ref2, ic_column = input$ct_hepidish_ic_col,
                                       method = "RPC", maxit = input$ct_adv_maxit %||% 50, nu.v = nu,
                                       constraint = input$ct_adv_constraint %||% "inequality")
      } else {
        res <- methyl_ct_run_epidish(f$mat, ref_use, method = input$ct_method,
                                      maxit = input$ct_adv_maxit %||% 50, nu.v = nu,
                                      constraint = input$ct_adv_constraint %||% "inequality")
      }
      validate(need(isTRUE(res$ok), res$reason %||% "Deconvolution failed."))
      c(res, list(ref_used = ref_use, working_mat = f$mat, sheet = f$sheet))
    })

    observeEvent(decon_result(), {
      decon_has_run(TRUE)
      r <- decon_result()
      if (!is.null(results)) {
        results$celltype <- list(
          method = r$method, cell_types = colnames(r$fractions), n_samples = nrow(r$fractions),
          n_markers_used = r$n_markers_used, mean_fraction = round(colMeans(r$fractions), 4)
        )
      }
    })

    output$ct_decon_result_ui <- renderUI({
      tagList(
        div(class = "table-toolbar", downloadButton(ns("ct_decon_download"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("ct_decon_table")),
        h4("Summary statistics per cell type"),
        DT::dataTableOutput(ns("ct_decon_summary_table"))
      )
    })
    output$ct_decon_table <- DT::renderDataTable({
      r <- decon_result()
      req(r)
      df <- as.data.frame(as.table(r$fractions)); colnames(df) <- c("Sample", "Cell Type", "Estimated Fraction")
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = "Estimated Fraction", digits = 4)
    })
    outputOptions(output, "ct_decon_table", suspendWhenHidden = FALSE)
    output$ct_decon_summary_table <- DT::renderDataTable({
      r <- decon_result()
      req(r)
      f <- r$fractions
      df <- data.frame(cell_type = colnames(f), min = apply(f, 2, min), max = apply(f, 2, max), mean = colMeans(f),
                        median = apply(f, 2, stats::median), sd = apply(f, 2, stats::sd))
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10, dom = "t"), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("min", "max", "mean", "median", "sd"), digits = 4)
    })
    outputOptions(output, "ct_decon_summary_table", suspendWhenHidden = FALSE)
    output$ct_decon_download <- downloadHandler(
      filename = function() "celltype_fractions.csv",
      content = function(file) {
        df <- as.data.frame(as.table(decon_result()$fractions)); colnames(df) <- c("sample", "cell_type", "fraction")
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    output$ct_composition_gate <- renderUI({
      if (!isTRUE(decon_has_run())) {
        return(div(class = "card", p(class = "empty-note", icon("circle-info"), "Run deconvolution first (Deconvolution tab) to see cell-composition figures.")))
      }
      uiOutput(ns("ct_composition_ui"))
    })

    output$ct_composition_ui <- renderUI({
      r <- decon_result()
      req(r)
      cts <- colnames(r$fractions)
      sheet_cols <- if (!is.null(r$sheet)) colnames(r$sheet) else character(0)
      tagList(
        box(
          width = 12, title = "A. Stacked bar - cell composition per sample", status = "primary", solidHeader = FALSE,
          fluidRow(
            column(4, selectInput(ns("ct_bar_order"), "Sample order",
                                   choices = c("As in data" = "asis", "Alphabetical" = "alpha", "By dominant cell type" = "dominant"))),
            column(8, checkboxGroupInput(ns("ct_bar_hide"), "Cell types shown", choices = cts, selected = cts, inline = TRUE))
          ),
          withSpinner(plotly::plotlyOutput(ns("ct_bar_plot"), height = "460px"), color = "#2563EB", type = 6),
          downloadButton(ns("ct_bar_download"), "Download PNG", class = "btn-sm")
        ),
        box(
          width = 12, title = "B. Cell-type heatmap", status = "primary", solidHeader = FALSE,
          fluidRow(
            column(4, checkboxInput(ns("ct_heat_cluster_rows"), "Cluster cell types", value = TRUE)),
            column(4, checkboxInput(ns("ct_heat_cluster_cols"), "Cluster samples", value = TRUE)),
            column(4, checkboxInput(ns("ct_heat_normalize"), "Row-normalize", value = FALSE))
          ),
          withSpinner(plotOutput(ns("ct_heat_plot"), height = "360px"), color = "#2563EB", type = 6),
          downloadButton(ns("ct_heat_download"), "Download PNG", class = "btn-sm")
        ),
        box(
          width = 12, title = "C. Cell-type distribution", status = "primary", solidHeader = FALSE,
          fluidRow(
            column(4, radioButtons(ns("ct_box_kind"), NULL, choices = c("Boxplot" = "box", "Violin" = "violin"), inline = TRUE)),
            column(8, selectInput(ns("ct_box_group"), "Group by (optional)", choices = c("None" = "", sheet_cols)))
          ),
          withSpinner(plotOutput(ns("ct_box_plot"), height = "360px"), color = "#2563EB", type = 6),
          downloadButton(ns("ct_box_download"), "Download PNG", class = "btn-sm")
        ),
        box(
          width = 12, title = "D. PCA / MDS of cell composition", status = "primary", solidHeader = FALSE,
          fluidRow(
            column(3, radioButtons(ns("ct_ord_method"), NULL, choices = c("PCA" = "pca", "MDS" = "mds"), inline = TRUE)),
            column(3, selectInput(ns("ct_ord_color"), "Color by",
                                   choices = c("None" = "", stats::setNames(cts, paste0("Cell type: ", cts)),
                                               if (length(sheet_cols) > 0) stats::setNames(sheet_cols, paste0("Phenotype: ", sheet_cols)) else character(0)))),
            column(3, checkboxInput(ns("ct_ord_labels"), "Show sample labels", value = FALSE))
          ),
          withSpinner(plotOutput(ns("ct_ord_plot"), height = "380px"), color = "#2563EB", type = 6),
          downloadButton(ns("ct_ord_download"), "Download PNG", class = "btn-sm")
        ),
        box(
          width = 12, title = "E. Correlation matrix between cell types", status = "primary", solidHeader = FALSE,
          withSpinner(plotOutput(ns("ct_corr_plot"), height = "360px"), color = "#2563EB", type = 6),
          downloadButton(ns("ct_corr_download"), "Download PNG", class = "btn-sm")
        ),
        box(
          width = 12, title = "Group Comparison", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Compare estimated cell fractions between groups from an uploaded phenotype/sample sheet."),
          fluidRow(
            column(6, selectInput(ns("ct_cmp_group_col"), "Grouping variable",
                                   choices = if (length(sheet_cols) > 0) sheet_cols else c("No sample sheet loaded" = ""))),
            column(6, div(style = "margin-top:24px;",
                          actionButton(ns("ct_run_cmp_btn"), "Run Cell-Type Comparison", icon = icon("play"), class = "btn-primary")))
          ),
          uiOutput(ns("ct_cmp_result_gate"))
        )
      )
    })

    ct_sample_order_vec <- reactive({
      r <- decon_result()
      f <- r$fractions
      switch(input$ct_bar_order %||% "asis",
             alpha = sort(rownames(f)),
             dominant = rownames(f)[order(apply(f, 1, which.max), -apply(f, 1, max))],
             rownames(f))
    })
    output$ct_bar_plot <- plotly::renderPlotly({
      r <- decon_result()
      req(r)
      plotly_safe(methyl_ct_plot_stacked_bar(r$fractions, sample_order = ct_sample_order_vec(),
                                              hide_types = setdiff(colnames(r$fractions), input$ct_bar_hide %||% colnames(r$fractions))))
    })
    output$ct_bar_download <- make_plot_dl(function() {
      r <- decon_result()
      methyl_ct_plot_stacked_bar(r$fractions, sample_order = ct_sample_order_vec(),
                                  hide_types = setdiff(colnames(r$fractions), input$ct_bar_hide %||% colnames(r$fractions)))
    }, "celltype_stacked_bar")

    output$ct_heat_plot <- renderPlot({
      r <- decon_result()
      req(r)
      methyl_ct_plot_heatmap(r$fractions, cluster_rows = isTRUE(input$ct_heat_cluster_rows),
                              cluster_cols = isTRUE(input$ct_heat_cluster_cols), normalize = isTRUE(input$ct_heat_normalize))
    })
    output$ct_heat_download <- make_plot_dl(function() {
      r <- decon_result()
      methyl_ct_plot_heatmap(r$fractions, cluster_rows = isTRUE(input$ct_heat_cluster_rows),
                              cluster_cols = isTRUE(input$ct_heat_cluster_cols), normalize = isTRUE(input$ct_heat_normalize))
    }, "celltype_heatmap")

    ct_group_vec <- function(col) {
      r <- decon_result()
      if (is.null(col) || !nzchar(col) || is.null(r$sheet) || !(col %in% colnames(r$sheet))) return(NULL)
      ids <- methyl_sheet_sample_ids(r$sheet, rownames(r$fractions))
      stats::setNames(as.character(r$sheet[[col]]), ids)
    }
    output$ct_box_plot <- renderPlot({
      r <- decon_result()
      req(r)
      methyl_ct_plot_box(r$fractions, group = ct_group_vec(input$ct_box_group), violin = identical(input$ct_box_kind, "violin"))
    })
    output$ct_box_download <- make_plot_dl(function() {
      r <- decon_result()
      methyl_ct_plot_box(r$fractions, group = ct_group_vec(input$ct_box_group), violin = identical(input$ct_box_kind, "violin"))
    }, "celltype_distribution")

    build_ord_plot <- function() {
      r <- decon_result()
      scores_res <- if (identical(input$ct_ord_method, "mds")) methyl_ct_composition_mds(r$fractions) else methyl_ct_composition_pca(r$fractions)
      validate(need(isTRUE(scores_res$ok), scores_res$reason %||% "Not enough samples for this ordination."))
      color_vec <- NULL; color_label <- "Group"
      sel <- input$ct_ord_color
      if (!is.null(sel) && nzchar(sel)) {
        if (sel %in% colnames(r$fractions)) { color_vec <- r$fractions[, sel]; color_label <- sel }
        else { color_vec <- ct_group_vec(sel); color_label <- sel }
      }
      methyl_ct_plot_scores(scores_res$scores, color_by = color_vec, color_label = color_label, show_labels = isTRUE(input$ct_ord_labels),
                             x_lab = if (identical(input$ct_ord_method, "mds")) "Dim 1" else sprintf("PC1 (%.1f%%)", scores_res$var_explained[1] * 100),
                             y_lab = if (identical(input$ct_ord_method, "mds")) "Dim 2" else sprintf("PC2 (%.1f%%)", scores_res$var_explained[2] * 100))
    }
    output$ct_ord_plot <- renderPlot({ req(decon_result()); build_ord_plot() })
    output$ct_ord_download <- make_plot_dl(build_ord_plot, "celltype_pca_mds")

    output$ct_corr_plot <- renderPlot({ r <- decon_result(); req(r); methyl_ct_plot_corr(stats::cor(r$fractions)) })
    output$ct_corr_download <- make_plot_dl(function() methyl_ct_plot_corr(stats::cor(decon_result()$fractions)), "celltype_correlation")

    cmp_has_run <- reactiveVal(FALSE)
    register_has_run_gate_local("ct_cmp_result_gate", cmp_has_run, "ct_cmp_result_ui",
                                 "Not run yet - choose a grouping variable and click \"Run Cell-Type Comparison\".")

    cmp_result <- eventReactive(input$ct_run_cmp_btn, {
      r <- decon_result()
      grp <- ct_group_vec(input$ct_cmp_group_col)
      validate(need(!is.null(grp), "No sample sheet / grouping column available."))
      grp <- grp[rownames(r$fractions)]
      gs <- methyl_ct_group_stats(r$fractions, grp)
      validate(need(isTRUE(gs$ok), gs$reason %||% "Comparison failed."))
      list(stats = gs, group = grp, fractions = r$fractions)
    })
    observeEvent(cmp_result(), { cmp_has_run(TRUE) })

    output$ct_cmp_result_ui <- renderUI({
      cr <- cmp_result()
      tagList(
        p(class = "empty-note", icon("circle-info"),
          sprintf("Test used: %s (%d group(s): %s).", cr$stats$test_used, length(cr$stats$levels), paste(cr$stats$levels, collapse = ", "))),
        withSpinner(plotOutput(ns("ct_cmp_plot"), height = "360px"), color = "#2563EB", type = 6),
        div(class = "table-toolbar", downloadButton(ns("ct_cmp_download"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("ct_cmp_table"))
      )
    })
    output$ct_cmp_plot <- renderPlot({ cr <- cmp_result(); req(cr); methyl_ct_plot_group_diff(cr$fractions, cr$group, cr$stats$table) })
    output$ct_cmp_table <- DT::renderDataTable({
      cr <- cmp_result()
      req(cr)
      DT::datatable(cr$stats$table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("statistic", "p", "effect_size", "mean_diff", "fdr"), digits = 4)
    })
    outputOptions(output, "ct_cmp_table", suspendWhenHidden = FALSE)
    output$ct_cmp_download <- downloadHandler(
      filename = function() "celltype_group_comparison.csv",
      content = function(file) utils::write.csv(cmp_result()$stats$table, file, row.names = FALSE)
    )

    val_has_run <- reactiveVal(FALSE)
    register_has_run_gate_local("ct_val_result_gate", val_has_run, "ct_val_result_ui",
                                 "Not run yet - click \"Run Validation\" above (needs a completed deconvolution run first).")

    val_result <- eventReactive(input$ct_run_val_btn, {
      r <- decon_result()
      validate(need(!is.null(r), "Run Cell-Type Deconvolution first."))
      recon <- methyl_ct_reconstruct(r$ref_used, r$fractions)
      v <- methyl_ct_validation_metrics(r$working_mat, recon)
      validate(need(isTRUE(v$ok), v$reason %||% "Validation failed."))
      v
    })
    observeEvent(val_result(), { val_has_run(TRUE) })

    output$ct_val_result_ui <- renderUI({
      v <- val_result()
      tagList(
        fluidRow(
          column(3, div(class = "card", p(strong("Correlation")), p(sprintf("%.4f", v$overall$cor)))),
          column(3, div(class = "card", p(strong("RMSE")), p(sprintf("%.4f", v$overall$rmse)))),
          column(3, div(class = "card", p(strong("MAE")), p(sprintf("%.4f", v$overall$mae)))),
          column(3, div(class = "card", p(strong("R²")), p(sprintf("%.4f", v$overall$r2))))
        ),
        withSpinner(plotOutput(ns("ct_val_plot"), height = "380px"), color = "#2563EB", type = 6),
        downloadButton(ns("ct_val_download"), "Download PNG", class = "btn-sm"),
        h4("Per-sample reconstruction metrics"),
        DT::dataTableOutput(ns("ct_val_table"))
      )
    })
    output$ct_val_plot <- renderPlot({ v <- val_result(); req(v); methyl_ct_plot_reconstruction(v$observed, v$reconstructed) })
    output$ct_val_download <- make_plot_dl(function() methyl_ct_plot_reconstruction(val_result()$observed, val_result()$reconstructed), "celltype_validation_reconstruction")
    output$ct_val_table <- DT::renderDataTable({
      v <- val_result()
      req(v)
      DT::datatable(v$per_sample, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("cor", "rmse", "mae"), digits = 4)
    })
    outputOptions(output, "ct_val_table", suspendWhenHidden = FALSE)

    cmpm_has_run <- reactiveVal(FALSE)
    register_has_run_gate_local("ct_cmpmethods_result_gate", cmpm_has_run, "ct_cmpmethods_result_ui",
                                 "Not run yet - pick methods and click \"Compare Methods\" above.")

    cmpm_result <- eventReactive(input$ct_run_cmpmethods_btn, {
      r <- decon_result()
      validate(need(!is.null(r), "Run Cell-Type Deconvolution first."))
      validate(need(length(input$ct_cmpmethods_pick) >= 2, "Pick at least 2 methods to compare."))
      cmp <- methyl_ct_compare_methods(r$working_mat, r$ref_used, methods = input$ct_cmpmethods_pick,
                                        maxit = input$ct_adv_maxit %||% 50,
                                        nu.v = c(input$ct_adv_nu1 %||% 0.25, input$ct_adv_nu2 %||% 0.5, input$ct_adv_nu3 %||% 0.75),
                                        constraint = input$ct_adv_constraint %||% "inequality")
      validate(need(isTRUE(cmp$ok), cmp$reason %||% "Method comparison failed."))
      cmp
    })
    observeEvent(cmpm_result(), { cmpm_has_run(TRUE) })

    output$ct_cmpmethods_result_ui <- renderUI({
      cmp <- cmpm_result()
      methods <- names(cmp$fractions_by_method)
      tagList(
        h4("Method-correlation heatmap"),
        withSpinner(plotOutput(ns("ct_cmpmethods_corr_plot"), height = "320px"), color = "#2563EB", type = 6),
        h4("Bland-Altman comparison"),
        fluidRow(
          column(4, selectInput(ns("ct_ba_a"), "Method A", choices = methods, selected = methods[1])),
          column(4, selectInput(ns("ct_ba_b"), "Method B", choices = methods, selected = methods[min(2, length(methods))]))
        ),
        withSpinner(plotOutput(ns("ct_cmpmethods_ba_plot"), height = "360px"), color = "#2563EB", type = 6),
        h4("Agreement summary"),
        DT::dataTableOutput(ns("ct_cmpmethods_summary_table"))
      )
    })
    output$ct_cmpmethods_corr_plot <- renderPlot({ cmp <- cmpm_result(); req(cmp); methyl_ct_plot_corr(methyl_ct_method_correlation(cmp$fractions_by_method)) })
    output$ct_cmpmethods_ba_plot <- renderPlot({
      cmp <- cmpm_result()
      req(cmp, input$ct_ba_a, input$ct_ba_b)
      req(input$ct_ba_a %in% names(cmp$fractions_by_method), input$ct_ba_b %in% names(cmp$fractions_by_method))
      methyl_ct_plot_bland_altman(cmp$fractions_by_method[[input$ct_ba_a]], cmp$fractions_by_method[[input$ct_ba_b]], input$ct_ba_a, input$ct_ba_b)
    })
    output$ct_cmpmethods_summary_table <- DT::renderDataTable({
      cmp <- cmpm_result()
      req(cmp)
      DT::datatable(methyl_ct_method_agreement_summary(cmp$fractions_by_method), rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 10, dom = "t"), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("mean_abs_diff", "max_abs_diff"), digits = 4)
    })
    outputOptions(output, "ct_cmpmethods_summary_table", suspendWhenHidden = FALSE)

    output$ct_export_beta <- downloadHandler(
      filename = function() "celltype_filtered_beta_matrix.csv",
      content = function(file) { f <- ct_filtered(); req(f); utils::write.csv(f$mat, file) }
    )
    output$ct_export_markers <- downloadHandler(
      filename = function() "celltype_selected_markers.csv",
      content = function(file) { req(fs_has_run()); utils::write.csv(fs_result()$selected, file, row.names = FALSE) }
    )
    output$ct_export_ref <- downloadHandler(
      filename = function() "celltype_reference_matrix.csv",
      content = function(file) { req(decon_has_run()); utils::write.csv(decon_result()$ref_used, file) }
    )
    output$ct_export_fractions <- downloadHandler(
      filename = function() "celltype_fractions.csv",
      content = function(file) { req(decon_has_run()); utils::write.csv(decon_result()$fractions, file) }
    )
    output$ct_export_pheno_fractions <- downloadHandler(
      filename = function() "celltype_fractions_with_phenotype.csv",
      content = function(file) {
        req(decon_has_run())
        r <- decon_result()
        df <- as.data.frame(r$fractions); df$sample <- rownames(df)
        if (!is.null(r$sheet)) {
          ids <- methyl_sheet_sample_ids(r$sheet, rownames(r$fractions))
          sheet2 <- r$sheet; rownames(sheet2) <- ids
          extra <- sheet2[df$sample, , drop = FALSE]
          df <- cbind(df, extra)
        }
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    output$ct_export_comparison <- downloadHandler(
      filename = function() "celltype_group_comparison.csv",
      content = function(file) { req(cmp_has_run()); utils::write.csv(cmp_result()$stats$table, file, row.names = FALSE) }
    )
    output$ct_export_report <- downloadHandler(
      filename = function() "celltype_analysis_summary.txt",
      content = function(file) {
        src <- tryCatch(ct_source(), error = function(e) NULL)
        f <- tryCatch(ct_filtered(), error = function(e) NULL)
        lines <- c(
          "ArthOMix - Methylomics Cell-Type Deconvolution - Analysis Summary",
          sprintf("Dataset: %s", if (!is.null(src)) src$source_label else "(not loaded)"),
          if (!is.null(f)) sprintf("Working matrix after QC filters: %s CpGs x %s samples", format(nrow(f$mat), big.mark = ","), ncol(f$mat)),
          if (isTRUE(fs_has_run())) sprintf("CpG feature selection: %s method, %d marker CpG(s) selected", fs_result()$method, nrow(fs_result()$selected)),
          if (isTRUE(decon_has_run())) c(
            sprintf("Reference: %d cell types (%s)", ncol(decon_result()$ref_used), paste(colnames(decon_result()$ref_used), collapse = ", ")),
            sprintf("Method: %s", decon_result()$method),
            sprintf("Marker CpGs used: %d", decon_result()$n_markers_used)
          ),
          if (isTRUE(val_has_run())) sprintf("Validation: correlation=%.4f, RMSE=%.4f, MAE=%.4f, R2=%.4f",
                                              val_result()$overall$cor, val_result()$overall$rmse, val_result()$overall$mae, val_result()$overall$r2)
        )
        writeLines(unlist(lines), file)
      }
    )

    observe({
      shinyjs::toggleState("ct_export_markers", condition = isTRUE(fs_has_run()))
      shinyjs::toggleState("ct_export_ref", condition = isTRUE(decon_has_run()))
      shinyjs::toggleState("ct_export_fractions", condition = isTRUE(decon_has_run()))
      shinyjs::toggleState("ct_export_pheno_fractions", condition = isTRUE(decon_has_run()))
      shinyjs::toggleState("ct_export_comparison", condition = isTRUE(cmp_has_run()))
    })

    NULL
  })
}
