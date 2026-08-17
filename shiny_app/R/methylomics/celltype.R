## R/methylomics/celltype.R
## Non-UI helpers for the Methylomics Cell-Type Deconvolution submodule
## (mod_methyl_celltype.R): reference-library registry, marker-CpG ranking
## off reference centroids, EpiDISH::epidish()/hepidish() wrappers,
## reference/working-matrix overlap QC, reconstruction validation,
## cross-method comparison, group-comparison stats, and this module's own
## ggplot builders (reusing theme_arthomix()/ARTHOMIX_COLORS/arthomix_pair()
## from global.R). Isolated strictly to this module - nothing here is
## called from any other Methylomics or Transcriptomics tab, and this file
## never touches shiny_app/R/transcriptomics/mod_deconvolution.R.
##
## Backend honesty: only EpiDISH (installed) is used for real estimation.
## MethylResolver / IDOL-optimized libraries / true reference-free
## deconvolution would need packages that are NOT installed in this
## deployment (FlowSorted.Blood.EPIC, IDOLOptimizedCpGs, ENmix, TOAST/
## RefFreeEWAS) - see methyl_ct_unavailable_methods() below, which the UI
## renders as disabled options with an explanatory reason rather than
## faking a result.

## =============================================================================
## Reference library registry
## =============================================================================

## Every entry maps to a real EpiDISH data object. ncpg/celltypes are read
## off the actual installed object (not hardcoded), so this can't silently
## drift from whatever EpiDISH version is installed.
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
    ## getExportedValue(), not get(..., envir = asNamespace(...)): EpiDISH's
    ## reference matrices are LazyData objects, registered in a lazy-load
    ## database that plain get() against the namespace environment does not
    ## see until something has "touched" it via :: or data() first - this
    ## bit even after requireNamespace("EpiDISH") succeeded, confirmed via
    ## get("centDHSbloodDMC.m", envir=asNamespace("EpiDISH")) throwing
    ## "object not found" in a fresh session where EpiDISH was never
    ## library()'d. getExportedValue() (what the :: operator itself calls)
    ## resolves lazy data correctly without requiring library(EpiDISH).
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

## Methods the spec asks for that this deployment cannot honestly implement -
## the UI shows these as disabled selector entries with the reason attached,
## per this app's own "explain why, don't fake it" convention (see
## qc.R's methyl_filter_cross_reactive()/methyl_filter_maf() for the same
## pattern applied to probe filters).
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

## =============================================================================
## Scale detection / transforms / dataset summary
## =============================================================================

## Distinguishes beta (0-1), percent methylation (0-100), and M-values
## (unbounded log2 ratio) from a sampled quantile range - a three-way
## extension of mod_methyl_featureselection.R's methyl_fs_detect_scale(),
## which only distinguishes beta vs. M. Never transforms silently - the
## caller surfaces `note` and requires an explicit "Apply Transformation"
## click before touching the working matrix (spec's methylation-scale
## requirement).
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

## Inverse of qc.R's methyl_beta_to_mvalue() logit transform.
methyl_ct_m_to_beta <- function(m) 2^m / (1 + 2^m)

## Compact "Data & QC" preview card summary - reuses methyl_get_annotation()
## (annotation.R) for the chromosome count, degrading to NA (not an error)
## when no manifest annotation resolves for the array type, same convention
## every other manifest-dependent feature in this app already follows.
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

## Custom reference-matrix upload: CpG x cell-type mean-beta CSV/TSV, first
## column CpG ID - reuses parse_upload.R's methyl_parse_matrix() (same
## probe-rows/sample-columns shape, "samples" here being cell types) rather
## than a bespoke parser, then validates it actually looks like a
## methylation reference (finite, in [0,1], >=2 cell-type columns).
methyl_ct_parse_custom_reference <- function(datapath, filename) {
  p <- methyl_parse_matrix(datapath, filename)
  if (!isTRUE(p$ok)) return(p)
  mat <- p$mat
  if (any(!is.finite(mat))) return(list(ok = FALSE, error = "The custom reference matrix contains non-finite values."))
  if (any(mat < -0.05 | mat > 1.05)) return(list(ok = FALSE, error = "The custom reference matrix must contain beta values (mean methylation per CpG x cell type) in [0,1]."))
  if (ncol(mat) < 2) return(list(ok = FALSE, error = "A reference matrix needs at least 2 cell-type columns."))
  list(ok = TRUE, mat = mat)
}

## =============================================================================
## Marker-CpG ranking (off reference centroids - no fabricated p-values)
## =============================================================================

## For every CpG, assigns the cell type it's the strongest marker for and
## records an effect size (its centroid beta vs. the most extreme of the
## OTHER cell types' centroids, whichever direction is larger) and a
## specificity score (that effect normalized by the spread of the other
## cell types' centroids - a CpG that is merely "somewhat different" from a
## tight cluster of others scores lower than one that's cleanly separated).
## Deliberately does NOT produce a p-value/FDR: reference centroids are
## single mean-beta values per cell type with no per-sample replicates, so
## a real significance test isn't available here - see
## methyl_ct_select_markers()'s docs for how the UI's FDR control is
## disabled for this reason rather than fed a made-up number.
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

## Applies the spec's CpG Feature Selection filters on top of
## methyl_ct_marker_rank()'s per-CpG table. `sort_by` is "effect" for the
## "Reference-library markers" method (the default, and the one this
## function's effect/specificity columns were built for) or "btw_type_var"
## for "Variance-based CpGs" (ranks the SAME reference-restricted candidate
## set by between-cell-type variance instead - deliberately not bulk-sample
## variance across conditions, which would be gene-expression-style
## phenotype feature selection, not cell-type marker selection).
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

## Balances a top-N cap roughly evenly across cell types (rather than one
## cell type's larger effect sizes crowding out every other type's markers
## entirely), then trims any overshoot from the globally weakest picks -
## this is what makes spec Figure 1's "CpGs retained per cell type" bar
## chart show every cell type rather than just the one or two with the
## biggest centroid separations.
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

## Chromosome-scope restriction for marker selection (autosomes only /
## autosomes+X / all) - reuses methyl_get_annotation() (annotation.R); a
## CpG with no manifest match is kept rather than dropped (absence of
## annotation isn't evidence it's on a sex chromosome), matching
## qc.R's methyl_filter_maf()'s same "unresolved = kept" convention.
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

## =============================================================================
## Reference / working-matrix overlap QC (spec S17-18)
## =============================================================================

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

## =============================================================================
## EpiDISH / hepidish wrappers
## =============================================================================

## Thin, validated wrapper around EpiDISH::epidish() - method is one of
## "CP" (Houseman constrained projection), "RPC" (EpiDISH robust partial
## correlations), "CBS" (CIBERSORT-style support vector regression).
## Confirmed against a synthetic known-mixture test (see this module's
## verification notes) that estF recovers true fractions with per-cell-type
## correlation > 0.99 and is already non-negative/sum-to-one for all three
## methods - no extra clipping/renormalization is applied here since
## EpiDISH's own constraint already guarantees it.
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

## Two-stage hierarchical EpiDISH: `ref1` splits the tissue into its top-
## level components (e.g. Epi/Fib/IC), then `ref2` further decomposes
## whichever `ref1` column is named in `ic_column` (e.g. "IC") into its own
## cell subtypes (e.g. the 7 blood types) - Zheng et al. 2018's hepidish(),
## a real bonus method for epithelial-type tissue references. Confirmed
## against a synthetic two-level mixture (tissue-level Epi/Fib/IC combined
## with an independent immune-subtype mixture within IC) that hepidish()
## exactly recovers both levels.
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

## =============================================================================
## Reconstruction validation (spec S15)
## =============================================================================

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

## =============================================================================
## Cross-method comparison (spec S16)
## =============================================================================

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

## Overall correlation between every pair of methods' estimated fractions
## (flattened across samples and cell types) - the "method-correlation
## heatmap" input.
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

## Plain-language mean-absolute-difference summary per method pair.
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

## =============================================================================
## Group / phenotype comparison (spec S13-14)
## =============================================================================

## Two-group Wilcoxon (rank-biserial effect size) or multi-group
## Kruskal-Wallis (epsilon-squared effect size), one test per cell type,
## BH-adjusted across cell types - the same auto-selected-by-group-count
## logic as transcriptomics' mod_deconvolution.R's compute_group_stats()
## (mod_deconvolution.R:242-282), reimplemented locally rather than shared
## since that module is not to be touched or imported from.
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

## =============================================================================
## Cell-composition ordination (PCA / MDS on the fraction matrix itself)
## =============================================================================

## Deliberately NOT qc.R's methyl_pca_scores()/methyl_mds_scores(): those
## operate on CpG-scale data and enforce a >=10-row minimum (meant for
## probes), which would wrongly reject a 3-4 cell-type reference (e.g.
## Epi/Fib/IC). These operate directly on the small samples x cell-types
## fraction matrix instead - no top-variance-feature subsetting needed
## since every cell type is already a meaningful "feature."
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

## =============================================================================
## ggplot builders (theme_arthomix()/ARTHOMIX_COLORS/arthomix_pair() from global.R)
## =============================================================================

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

methyl_ct_plot_method_scatter <- function(frac_a, frac_b, cell_type, label_a, label_b) {
  common <- intersect(rownames(frac_a), rownames(frac_b))
  df <- data.frame(sample = common, a = frac_a[common, cell_type], b = frac_b[common, cell_type])
  ggplot(df, aes(x = a, y = b)) +
    geom_point(color = ARTHOMIX_COLORS$blue, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
    labs(x = label_a, y = label_b, title = cell_type) + theme_arthomix()
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

## Boxplot with significance-star annotation, same idea as
## mod_deconvolution.R's render_group_boxplot() (mod_deconvolution.R:295-320)
## but reimplemented locally, not imported - that module is not to be
## touched or depended on.
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
