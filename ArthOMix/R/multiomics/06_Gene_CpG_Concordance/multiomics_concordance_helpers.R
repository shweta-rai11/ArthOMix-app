## R/multiomics/06_Gene_CpG_Concordance/multiomics_concordance_helpers.R
## Data-adaptive engine for the "Gene-CpG Concordance" submodule
## (mod_multi_concordance.R). Pure functions only, no Shiny reactives here.

MCC_REGION_PROMOTER <- c("TSS200", "TSS1500", "5'UTR")
MCC_REGION_BODY <- c("Body", "3'UTR", "1stExon", "ExonBnd")

MCC_DIRECTION_LEVELS <- c(
  "Up expression + Hypomethylation", "Down expression + Hypermethylation",
  "Up expression + Hypermethylation", "Down expression + Hypomethylation",
  "Weak/uncertain", "Not interpretable"
)

mcc_layer_candidates <- function(multi_dataset, omics_type) {
  layers <- multi_dataset$layers %||% list()
  meta <- multi_dataset$layer_meta %||% list()
  if (length(layers) == 0) return(character(0))
  by_meta <- names(layers)[vapply(names(layers), function(nm) identical(meta[[nm]]$omics_type, omics_type), logical(1))]
  if (length(by_meta) > 0) return(by_meta)
  rx <- switch(omics_type, rnaseq = "transcript|rna|express|gene", methylation = "methyl|cpg|beta|meth", NULL)
  if (is.null(rx)) return(character(0))
  names(layers)[grepl(rx, names(layers), ignore.case = TRUE)]
}

mcc_default_layer <- function(candidates, layers) {
  if (length(candidates) == 0) return(NULL)
  candidates[1]
}

mcc_data_status <- function(multi_dataset, multi_results, expr_layer = NULL, meth_layer = NULL, array_type = "450K") {
  layers <- multi_dataset$layers %||% list()
  meta <- multi_dataset$sample_meta
  expr_mat <- if (!is.null(expr_layer) && expr_layer %in% names(layers)) layers[[expr_layer]] else NULL
  meth_mat <- if (!is.null(meth_layer) && meth_layer %in% names(layers)) layers[[meth_layer]] else NULL

  pairing <- if (!is.null(expr_mat) && !is.null(meth_mat)) cx_detect_sample_pairing(rownames(expr_mat), rownames(meth_mat)) else NULL

  has_sex_col <- length(mcc_sex_candidates(meta)) > 0
  anno_pkg_available <- isTRUE(tryCatch(requireNamespace(CX_METH_ANNOTATION_PACKAGES[[array_type]] %||% "", quietly = TRUE), error = function(e) FALSE))
  candidate_pool_available <- !is.null(multi_results$biomarker$df) || !is.null(multi_results$stratification$clusters)

  row <- function(item, available, detail) data.frame(item = item, status = if (isTRUE(available)) "Available" else "Missing", detail = detail, stringsAsFactors = FALSE)
  rbind(
    row("Expression data", !is.null(expr_mat), if (!is.null(expr_mat)) sprintf("%s (%d genes x %d samples)", expr_layer, ncol(expr_mat), nrow(expr_mat)) else "No expression layer selected/available."),
    row("Methylation data", !is.null(meth_mat), if (!is.null(meth_mat)) sprintf("%s (%d CpGs x %d samples)", meth_layer, ncol(meth_mat), nrow(meth_mat)) else "No methylation layer selected/available."),
    row("Matched samples", !is.null(pairing) && isTRUE(pairing$paired), if (!is.null(pairing)) sprintf("%d matched of %d expression / %d methylation samples.", pairing$n_common, pairing$n_expr, pairing$n_meth) else "Needs both expression and methylation data."),
    row("Genes", !is.null(expr_mat), if (!is.null(expr_mat)) format(ncol(expr_mat), big.mark = ",") else "0"),
    row("CpGs", !is.null(meth_mat), if (!is.null(meth_mat)) format(ncol(meth_mat), big.mark = ",") else "0"),
    row("Clinical metadata", !is.null(meta) && ncol(meta) > 0, if (!is.null(meta)) sprintf("%d variable(s): %s", ncol(meta), paste(utils::head(colnames(meta), 6), collapse = ", ")) else "No sample metadata uploaded."),
    row("Sex variable", has_sex_col, if (has_sex_col) paste(mcc_sex_candidates(meta), collapse = ", ") else "No column named sex/gender found in metadata."),
    row("CpG annotation", anno_pkg_available, if (anno_pkg_available) sprintf("%s annotation package installed.", array_type) else sprintf("%s annotation package not installed in this deployment.", array_type)),
    row("Genomic coordinates", anno_pkg_available, "Derived from the same CpG annotation package (chr/pos)."),
    row("Candidate biomarker panel", candidate_pool_available, if (candidate_pool_available) "Biomarker Discovery and/or Patient Stratification results are loaded in this session." else "Run Biomarker Discovery (DIABLO) and/or Patient Stratification (SNF) first, or supply custom genes/CpGs below.")
  )
}

mcc_detect_id_type <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return("unknown")
  frac <- function(rx) mean(grepl(rx, ids))
  if (frac("^ENSG[0-9]+") > 0.5) return("Ensembl Gene ID")
  if (frac("^cg[0-9]+$") > 0.5) return("Illumina CpG probe ID")
  if (frac("^[0-9]+$") > 0.5) return("Entrez ID")
  "Gene symbol"
}

mcc_detect_methylation_value_type <- function(mat) {
  v <- as.numeric(mat)
  v <- v[is.finite(v)]
  if (length(v) == 0) return("unknown")
  if (min(v) >= -0.001 && max(v) <= 1.001) return("beta")
  if (min(v) < 0 && max(v) > 0 && max(abs(v)) < 20) return("M-value")
  "delta/other (already-differenced or non-standard scale)"
}

mcc_match_samples <- function(expr_mat, meth_mat) {
  if (is.null(expr_mat) || is.null(meth_mat)) return(list(ok = FALSE, error = "Both expression and methylation data are required for sample matching."))
  pairing <- cx_detect_sample_pairing(rownames(expr_mat), rownames(meth_mat))
  list(
    ok = isTRUE(pairing$paired), pairing = pairing,
    n_matched = pairing$n_common, n_removed_expr = pairing$n_expr - pairing$n_common,
    n_removed_meth = pairing$n_meth - pairing$n_common, common_samples = pairing$common_samples
  )
}

mcc_diablo_candidates <- function(multi_results) {
  rows <- list()
  df <- multi_results$biomarker$df
  if (!is.null(df) && nrow(df) > 0) {
    rows$freeform <- data.frame(
      feature = df$feature, omics = df$omics, diablo = TRUE, joint = TRUE,
      component = df$component, loading = df$loading,
      selection_frequency = df$selection_frequency, stability_category = df$stability_category,
      stringsAsFactors = FALSE
    )
  }
  strat <- multi_results$integration_stratified$result
  if (!is.null(strat) && isTRUE(strat$ok) && identical(strat$engine, "diablo") && !is.null(strat$panels) && nrow(strat$panels) > 0) {
    sp <- strat$panels
    rows$stratified <- data.frame(
      feature = sp$feature,
      omics = ifelse(sp$view == "expression", "Transcriptomics", ifelse(sp$view == "methylation", "Methylomics", sp$view)),
      diablo = TRUE, joint = TRUE, component = 1L, loading = sp$loading,
      selection_frequency = NA_real_, stability_category = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out[!duplicated(out$feature), , drop = FALSE]
}

mcc_join_preloaded_diablo_panel <- function(df, panel) {
  df$diablo <- FALSE; df$joint <- FALSE
  df$diablo_status <- NA_character_; df$diablo_loading <- NA_real_
  if (is.null(panel) || nrow(panel) == 0) return(df)

  has_sex <- "sex" %in% colnames(df) && "sex" %in% colnames(panel)
  key <- function(sex, feature) toupper(paste(if (has_sex) sex else "", feature))

  gene_panel <- panel[panel$omics == "transcriptomics_gene", , drop = FALSE]
  cpg_panel <- panel[panel$omics == "methylation_CpG", , drop = FALSE]

  idx_g <- match(key(df$sex, df$gene_id), key(gene_panel$sex, gene_panel$feature))
  idx_c <- match(key(df$sex, df$cpg), key(cpg_panel$sex, cpg_panel$feature))
  hit_g <- !is.na(idx_g); hit_c <- !is.na(idx_c)

  df$diablo[hit_g | hit_c] <- TRUE
  df$joint[hit_g | hit_c] <- TRUE
  df$diablo_status[hit_g] <- gene_panel$biomarker_status[idx_g[hit_g]]
  df$diablo_status[hit_c] <- cpg_panel$biomarker_status[idx_c[hit_c]]
  df$diablo_loading[hit_g] <- gene_panel$loading[idx_g[hit_g]]
  df$diablo_loading[hit_c] <- cpg_panel$loading[idx_c[hit_c]]
  df
}

mcc_snf_candidates <- function(multi_results, multi_dataset, expr_layer, meth_layer, fdr_thresh = 0.1, top_n = 500) {
  clusters <- multi_results$stratification$clusters
  if (is.null(clusters) || length(clusters) == 0) return(NULL)
  layers <- multi_dataset$layers %||% list()
  out <- list()
  if (!is.null(expr_layer) && expr_layer %in% names(layers)) {
    r <- sfc_feature_ranking(layers[[expr_layer]], clusters, "Transcriptomics", top_n = top_n)
    if (isTRUE(r$ok)) out$expr <- data.frame(feature = r$table$feature, omics = "Transcriptomics", snf = r$table$p_fdr < fdr_thresh, snf_p_fdr = r$table$p_fdr, stringsAsFactors = FALSE)
  }
  if (!is.null(meth_layer) && meth_layer %in% names(layers)) {
    r <- sfc_feature_ranking(layers[[meth_layer]], clusters, "Methylomics", top_n = top_n)
    if (isTRUE(r$ok)) out$meth <- data.frame(feature = r$table$feature, omics = "Methylomics", snf = r$table$p_fdr < fdr_thresh, snf_p_fdr = r$table$p_fdr, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

mcc_candidate_pool <- function(multi_results, multi_dataset, expr_layer, meth_layer,
                                custom_genes = character(0), custom_cpgs = character(0),
                                snf_fdr_thresh = 0.1) {
  diablo <- mcc_diablo_candidates(multi_results)
  snf <- mcc_snf_candidates(multi_results, multi_dataset, expr_layer, meth_layer, snf_fdr_thresh)
  all_features <- unique(c(diablo$feature, snf$feature, custom_genes, custom_cpgs))
  if (length(all_features) == 0) return(list(ok = FALSE, df = NULL, note = "No candidate biomarkers available."))

  base <- data.frame(feature = all_features, stringsAsFactors = FALSE)
  base$omics <- NA_character_
  base$diablo <- FALSE; base$joint <- FALSE; base$snf <- FALSE; base$custom <- FALSE
  base$selection_frequency <- NA_real_; base$stability_category <- NA_character_

  if (!is.null(diablo)) {
    idx <- match(base$feature, diablo$feature)
    hit <- !is.na(idx)
    base$omics[hit] <- diablo$omics[idx[hit]]
    base$diablo[hit] <- TRUE; base$joint[hit] <- TRUE
    base$selection_frequency[hit] <- diablo$selection_frequency[idx[hit]]
    base$stability_category[hit] <- diablo$stability_category[idx[hit]]
  }
  if (!is.null(snf)) {
    idx <- match(base$feature, snf$feature)
    hit <- !is.na(idx) & snf$snf[idx]
    base$omics[is.na(base$omics) & !is.na(idx)] <- snf$omics[idx[is.na(base$omics) & !is.na(idx)]]
    base$snf[hit] <- TRUE
  }
  base$custom[base$feature %in% custom_genes | base$feature %in% custom_cpgs] <- TRUE
  base$id_type <- vapply(base$feature, function(f) mcc_detect_id_type(f), character(1))
  base$omics[is.na(base$omics) & base$id_type == "Illumina CpG probe ID"] <- "Methylomics"
  base$omics[is.na(base$omics) & base$id_type %in% c("Ensembl Gene ID", "Entrez ID", "Gene symbol")] <- "Transcriptomics"

  list(ok = TRUE, df = base, note = NULL)
}

mcc_filter_source <- function(pool_df, source) {
  if (is.null(pool_df) || nrow(pool_df) == 0) return(pool_df)
  switch(source,
    "All candidates" = pool_df,
    "DIABLO" = pool_df[pool_df$diablo, , drop = FALSE],
    "SNF" = pool_df[pool_df$snf, , drop = FALSE],
    "Joint Biomarker Discovery" = pool_df[pool_df$joint, , drop = FALSE],
    "DIABLO + SNF" = pool_df[pool_df$diablo & pool_df$snf, , drop = FALSE],
    "DIABLO + Joint" = pool_df[pool_df$diablo & pool_df$joint, , drop = FALSE],
    "SNF + Joint" = pool_df[pool_df$snf & pool_df$joint, , drop = FALSE],
    "Shared candidates" = pool_df[(pool_df$diablo + pool_df$snf + pool_df$joint) >= 2, , drop = FALSE],
    "Custom genes" = pool_df[pool_df$custom & pool_df$omics == "Transcriptomics", , drop = FALSE],
    "Custom CpGs" = pool_df[pool_df$custom & pool_df$omics == "Methylomics", , drop = FALSE],
    pool_df
  )
}

mcc_gene_cpg_map <- function(candidate_genes, meth_features = NULL, array_type = "450K") {
  candidate_genes <- unique(candidate_genes[!is.na(candidate_genes) & nzchar(candidate_genes)])
  if (length(candidate_genes) == 0) return(list(ok = FALSE, df = NULL, error = "No candidate genes to map."))
  ar <- cx_get_region_annotation(array_type)
  if (!isTRUE(ar$ok)) return(list(ok = FALSE, df = NULL, error = ar$reason))
  anno <- ar$anno

  harm <- cx_harmonize_gene_ids(candidate_genes)
  symbols <- if (isTRUE(harm$ok)) unique(stats::na.omit(harm$df$canonical_symbol)) else unique(toupper(candidate_genes))
  if (length(symbols) == 0) return(list(ok = FALSE, df = NULL, error = "None of the candidate genes could be resolved to a symbol for annotation matching."))

  hit <- toupper(anno$gene) %in% toupper(symbols)
  sub <- anno[hit, , drop = FALSE]
  if (!is.null(meth_features)) sub <- sub[rownames(sub) %in% meth_features, , drop = FALSE]
  if (nrow(sub) == 0) return(list(ok = FALSE, df = NULL, error = "No CpGs are annotated to these candidate genes in the selected array's annotation (or none of them are present in this dataset's methylation layer)."))

  df <- data.frame(
    gene_symbol = sub$gene, gene_id = NA_character_, transcript_id = NA_character_,
    cpg = rownames(sub), chr = sub$chr, pos = sub$pos, strand = NA_character_,
    region_raw = sub$region_raw, region_fine = cx_region_fine(sub$region_raw),
    island_context = sub$island_context, tss_distance = NA_real_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(harm$ok)) {
    idx <- match(toupper(df$gene_symbol), toupper(harm$df$canonical_symbol))
    df$gene_id <- harm$df$ensembl_id[idx]
  }
  list(ok = TRUE, df = df, error = NULL,
       note = "Transcript ID, strand, and distance-to-TSS are not exposed by this annotation source and are reported as Not available rather than estimated.")
}

mcc_classify_direction <- function(pairs_df, expr_thresh, expr_fdr_thresh, meth_thresh, meth_fdr_thresh) {
  need <- c("log2fc", "expr_fdr", "dbeta", "meth_fdr")
  if (!all(need %in% colnames(pairs_df))) stop("mcc_classify_direction: pairs_df missing required columns.")
  df <- cx_classify(pairs_df, expr_thresh, expr_fdr_thresh, meth_thresh, meth_fdr_thresh)

  region_fine <- df$region_fine %||% rep(NA_character_, nrow(df))
  promoter <- region_fine %in% MCC_REGION_PROMOTER
  body <- region_fine %in% MCC_REGION_BODY
  cat_chr <- as.character(df$category)

  df$direction_classification <- "Not interpretable"
  has_stats <- !is.na(df$log2fc) & !is.na(df$dbeta)
  both_sig <- df$sig_expression & df$sig_methylation
  df$direction_classification[has_stats & !both_sig] <- "Weak/uncertain"
  df$direction_classification[has_stats & both_sig & cat_chr == "Hypo + Up"] <- "Up expression + Hypomethylation"
  df$direction_classification[has_stats & both_sig & cat_chr == "Hyper + Down"] <- "Down expression + Hypermethylation"
  df$direction_classification[has_stats & both_sig & cat_chr == "Hyper + Up"] <- "Up expression + Hypermethylation"
  df$direction_classification[has_stats & both_sig & cat_chr == "Hypo + Down"] <- "Down expression + Hypomethylation"
  df$direction_classification <- factor(df$direction_classification, levels = MCC_DIRECTION_LEVELS)

  sig_relevant <- has_stats & both_sig
  df$canonical <- NA
  df$canonical[promoter & sig_relevant] <- cat_chr[promoter & sig_relevant] %in% c("Hyper + Down", "Hypo + Up")
  df$canonical[body & sig_relevant] <- cat_chr[body & sig_relevant] %in% c("Hyper + Up", "Hypo + Down")
  df$canonical_label <- ifelse(is.na(df$canonical), "Not applicable", ifelse(df$canonical, "Canonical", "Non-canonical"))
  df
}

MCC_CANONICAL_RULE_TEXT <- "Canonical: inverse methylation-expression relationship at promoter/TSS CpGs, concordant at gene-body CpGs. Other regions: Not applicable."

MCC_RAW_DIRECTION_LEVELS <- c(
  "Up expression + Hypomethylation", "Down expression + Hypermethylation",
  "Up expression + Hypermethylation", "Down expression + Hypomethylation"
)

mcc_add_raw_direction <- function(pairs_df) {
  need <- c("methylation_direction", "expression_direction")
  if (!all(need %in% colnames(pairs_df))) stop("mcc_add_raw_direction: pairs_df missing required columns.")
  df <- pairs_df
  flag <- function(col) if (col %in% colnames(df)) df[[col]] %in% TRUE else rep(FALSE, nrow(df))

  meth <- df$methylation_direction; expr <- df$expression_direction
  raw <- rep(NA_character_, nrow(df))
  raw[meth == "Hypo" & expr == "Upregulated"] <- MCC_RAW_DIRECTION_LEVELS[1]
  raw[meth == "Hyper" & expr == "Downregulated"] <- MCC_RAW_DIRECTION_LEVELS[2]
  raw[meth == "Hyper" & expr == "Upregulated"] <- MCC_RAW_DIRECTION_LEVELS[3]
  raw[meth == "Hypo" & expr == "Downregulated"] <- MCC_RAW_DIRECTION_LEVELS[4]
  df$raw_direction <- factor(raw, levels = MCC_RAW_DIRECTION_LEVELS)

  region_fine <- df$region_fine %||% rep(NA_character_, nrow(df))
  promoter <- !is.na(raw) & region_fine %in% MCC_REGION_PROMOTER
  body <- !is.na(raw) & region_fine %in% MCC_REGION_BODY
  raw_canonical <- rep(NA, nrow(df))
  raw_canonical[promoter] <- raw[promoter] %in% c(MCC_RAW_DIRECTION_LEVELS[1], MCC_RAW_DIRECTION_LEVELS[2])
  raw_canonical[body] <- raw[body] %in% c(MCC_RAW_DIRECTION_LEVELS[3], MCC_RAW_DIRECTION_LEVELS[4])
  df$raw_canonical_label <- ifelse(is.na(raw_canonical), "Not applicable", ifelse(raw_canonical, "Canonical", "Non-canonical"))

  src <- data.frame(DIABLO = flag("diablo"), SNF = flag("snf"), Joint = flag("joint"))
  df$biomarker_source <- apply(src, 1, function(r) paste(names(r)[r], collapse = " + "))
  df
}

mcc_join_genome_wide_significance <- function(df, deg, dmp) {
  df$genome_expr_p <- NA_real_; df$genome_expr_fdr <- NA_real_; df$genome_expr_logfc <- NA_real_
  df$genome_meth_p <- NA_real_; df$genome_meth_fdr <- NA_real_; df$genome_meth_dbeta <- NA_real_
  has_sex <- "sex" %in% colnames(df)
  key <- function(sex, id) toupper(paste(if (has_sex) sex else "", id))

  if (!is.null(deg) && nrow(deg) > 0 && "gene_id" %in% colnames(df)) {
    idx <- match(key(df$sex, df$gene_id), key(deg$sex, deg$ensembl_id))
    hit <- !is.na(idx)
    df$genome_expr_p[hit] <- deg$p_response[idx[hit]]
    df$genome_expr_fdr[hit] <- deg$fdr_response[idx[hit]]
    df$genome_expr_logfc[hit] <- deg$logFC_response[idx[hit]]
  }
  if (!is.null(dmp) && nrow(dmp) > 0 && "cpg" %in% colnames(df)) {
    idx <- match(key(df$sex, df$cpg), key(dmp$sex, dmp$CpG))
    hit <- !is.na(idx)
    df$genome_meth_p[hit] <- dmp$p_response[idx[hit]]
    df$genome_meth_fdr[hit] <- dmp$fdr_response[idx[hit]]
    df$genome_meth_dbeta[hit] <- dmp$delta_M_response[idx[hit]]
  }
  df
}

mcc_biomarker_direction_table <- function(pairs_df) {
  need <- c("gene_symbol", "cpg", "methylation_direction", "expression_direction")
  if (!all(need %in% colnames(pairs_df))) stop("mcc_biomarker_direction_table: pairs_df missing required columns.")
  flag <- function(col) if (col %in% colnames(pairs_df)) pairs_df[[col]] %in% TRUE else rep(FALSE, nrow(pairs_df))
  is_biomarker <- flag("diablo") | flag("snf") | flag("joint")
  df <- pairs_df[is_biomarker, , drop = FALSE]
  if (nrow(df) == 0) return(df)
  mcc_add_raw_direction(df)
}

mcc_pair_correlation <- function(expr_mat, meth_mat, pairs_df, common_samples, method = c("pearson", "spearman"), min_n = 3L) {
  method <- match.arg(method)
  if (length(common_samples) < min_n) return(list(ok = FALSE, df = NULL, error = "Not enough matched samples for correlation."))
  rows <- lapply(seq_len(nrow(pairs_df)), function(i) {
    g <- pairs_df$gene_symbol[i]; cpg <- pairs_df$cpg[i]
    if (!(g %in% colnames(expr_mat)) || !(cpg %in% colnames(meth_mat))) return(data.frame(gene_symbol = g, cpg = cpg, r = NA_real_, p = NA_real_, n = 0L))
    x <- as.numeric(expr_mat[common_samples, g]); y <- as.numeric(meth_mat[common_samples, cpg])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_n) return(data.frame(gene_symbol = g, cpg = cpg, r = NA_real_, p = NA_real_, n = sum(ok)))
    ct <- tryCatch(stats::cor.test(x[ok], y[ok], method = method), error = function(e) NULL)
    if (is.null(ct)) return(data.frame(gene_symbol = g, cpg = cpg, r = NA_real_, p = NA_real_, n = sum(ok)))
    data.frame(gene_symbol = g, cpg = cpg, r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
  })
  df <- do.call(rbind, rows)
  df$fdr <- cx_adjust_p(df$p, "BH")
  list(ok = TRUE, df = df, error = NULL)
}

mcc_regression <- function(expr_vec, meth_vec, covariates_df = NULL, covariate_cols = character(0)) {
  d <- data.frame(expression = expr_vec, methylation = meth_vec)
  if (!is.null(covariates_df) && length(covariate_cols) > 0) {
    keep <- intersect(covariate_cols, colnames(covariates_df))
    for (cc in keep) d[[cc]] <- covariates_df[[cc]]
  }
  d <- d[stats::complete.cases(d), , drop = FALSE]
  n_terms <- 2 + length(intersect(covariate_cols, colnames(d)))
  if (nrow(d) < max(10, 5 * n_terms)) return(list(ok = FALSE, error = sprintf("Too few complete observations (N=%d) for a reliable model with %d term(s) - need at least %d.", nrow(d), n_terms, max(10, 5 * n_terms))))
  form <- stats::as.formula(paste("expression ~", paste(c("methylation", intersect(covariate_cols, colnames(d))), collapse = " + ")))
  fit <- tryCatch(stats::lm(form, data = d), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE, error = "Model fit failed."))
  co <- summary(fit)$coefficients
  if (!"methylation" %in% rownames(co)) return(list(ok = FALSE, error = "Methylation term dropped from the model (collinear with a covariate)."))
  list(ok = TRUE, coefficient = co["methylation", "Estimate"], se = co["methylation", "Std. Error"],
       p_value = co["methylation", "Pr(>|t|)"], model_n = nrow(d), covariates_used = intersect(covariate_cols, colnames(d)))
}

mcc_sex_candidates <- function(sample_meta) multi_sex_candidates(sample_meta)

mcc_sex_groups <- function(sample_meta, sex_col, sample_ids) multi_sex_groups(sample_meta, sex_col, sample_ids)

mcc_priority_score <- function(df) {
  norm01 <- function(x) { x <- abs(x); if (all(is.na(x))) return(rep(NA_real_, length(x))); rng <- range(x, na.rm = TRUE); if (diff(rng) == 0) return(ifelse(is.na(x), NA_real_, 0.5)); (x - rng[1]) / diff(rng) }
  expr_component <- norm01(df$log2fc) %||% rep(NA_real_, nrow(df))
  meth_component <- norm01(df$dbeta) %||% rep(NA_real_, nrow(df))
  fdr_component <- 1 - norm01(pmin(df$expr_fdr, df$meth_fdr, na.rm = FALSE))
  cor_component <- norm01(df$correlation_r)
  region_component <- ifelse(is.na(df$canonical), 0.3, ifelse(df$canonical, 1, 0))
  diablo_component <- as.numeric(isTRUE(df$diablo) | df$diablo %in% TRUE)
  snf_component <- as.numeric(df$snf %in% TRUE)
  joint_component <- as.numeric(df$joint %in% TRUE)

  parts <- data.frame(expr_component, meth_component, fdr_component, cor_component, region_component, diablo_component, snf_component, joint_component)
  w <- c(expr_component = 0.15, meth_component = 0.15, fdr_component = 0.15, cor_component = 0.15, region_component = 0.1, diablo_component = 0.1, snf_component = 0.1, joint_component = 0.1)
  score <- rowSums(mapply(function(col, wt) ifelse(is.na(parts[[col]]), 0, parts[[col]]) * wt, names(w), w))
  parts$priority_score <- round(100 * score, 1)
  parts$evidence_label <- ifelse(parts$priority_score >= 60, "Potential Multi-Omics Biomarker", "Candidate Multi-Omics Biomarker")
  parts
}

mcc_design_candidates <- function(sample_meta, sample_ids) {
  if (is.null(sample_meta) || ncol(sample_meta) == 0) return(character(0))
  ok <- vapply(colnames(sample_meta), function(cl) {
    s <- mi_outcome_summary(sample_meta, cl, sample_ids)
    !is.null(s) && identical(s$type, "categorical") && identical(s$n_classes, 2L)
  }, logical(1))
  colnames(sample_meta)[ok]
}

mcc_expression_stats <- function(expr_mat, group) {
  common <- intersect(rownames(expr_mat), names(group))
  if (length(common) < 4) return(list(ok = FALSE, error = "Too few matched samples for a two-group expression comparison."))
  g <- factor(group[common])
  if (nlevels(g) != 2) return(list(ok = FALSE, error = "Design column does not have exactly two classes among the matched samples."))
  m <- expr_mat[common, , drop = FALSE]
  lv <- levels(g)
  rows <- apply(m, 2, function(col) {
    a <- col[g == lv[1]]; b <- col[g == lv[2]]
    lfc <- mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE)
    p <- tryCatch(stats::t.test(b, a)$p.value, error = function(e) NA_real_)
    c(log2fc = lfc, p = p)
  })
  df <- data.frame(feature = colnames(m), log2fc = rows["log2fc", ], p = rows["p", ], stringsAsFactors = FALSE)
  df$fdr <- stats::p.adjust(df$p, method = "BH")
  list(ok = TRUE, df = df, groups = lv, n = length(common))
}

mcc_methylation_stats <- function(meth_mat, group, value_type = c("beta", "M-value", "delta/other (already-differenced or non-standard scale)")) {
  value_type <- match.arg(value_type)
  common <- intersect(rownames(meth_mat), names(group))
  if (length(common) < 4) return(list(ok = FALSE, error = "Too few matched samples for a two-group methylation comparison."))
  g <- factor(group[common])
  if (nlevels(g) != 2) return(list(ok = FALSE, error = "Design column does not have exactly two classes among the matched samples."))
  m <- meth_mat[common, , drop = FALSE]
  beta <- if (identical(value_type, "beta")) m else if (identical(value_type, "M-value")) methyl_ct_m_to_beta(m) else NULL
  mval <- if (identical(value_type, "M-value")) m else if (identical(value_type, "beta")) methyl_beta_to_mvalue(m) else m
  lv <- levels(g)
  rows <- apply(mval, 2, function(col) {
    a <- col[g == lv[1]]; b <- col[g == lv[2]]
    p <- tryCatch(stats::t.test(b, a)$p.value, error = function(e) NA_real_)
    c(dm = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE), p = p)
  })
  dbeta <- if (!is.null(beta)) apply(beta, 2, function(col) mean(col[g == lv[2]], na.rm = TRUE) - mean(col[g == lv[1]], na.rm = TRUE)) else rep(NA_real_, ncol(m))
  df <- data.frame(feature = colnames(m), dbeta = rows["dm", ], delta_beta = dbeta, p = rows["p", ], stringsAsFactors = FALSE)
  df$fdr <- stats::p.adjust(df$p, method = "BH")
  list(ok = TRUE, df = df, groups = lv, n = length(common), value_type = value_type)
}

mcc_summary_counts <- function(pairs_df, sex_col_present = FALSE) {
  if (is.null(pairs_df) || nrow(pairs_df) == 0) return(NULL)
  sig <- !is.na(pairs_df$sig_expression) & !is.na(pairs_df$sig_methylation) & pairs_df$sig_expression & pairs_df$sig_methylation
  list(
    n_genes = length(unique(pairs_df$gene_symbol)),
    n_cpgs = length(unique(pairs_df$cpg)),
    n_pairs = nrow(pairs_df),
    n_significant = sum(sig, na.rm = TRUE),
    n_canonical = sum(pairs_df$canonical %in% TRUE),
    n_noncanonical = sum(pairs_df$canonical %in% FALSE),
    n_potential = sum(pairs_df$evidence_label == "Potential Multi-Omics Biomarker", na.rm = TRUE),
    n_female = if (sex_col_present && "sex" %in% colnames(pairs_df)) sum(tolower(pairs_df$sex) == "female", na.rm = TRUE) else NA_integer_,
    n_male = if (sex_col_present && "sex" %in% colnames(pairs_df)) sum(tolower(pairs_df$sex) == "male", na.rm = TRUE) else NA_integer_,
    n_diablo = sum(pairs_df$diablo %in% TRUE),
    n_snf = sum(pairs_df$snf %in% TRUE),
    n_joint = sum(pairs_df$joint %in% TRUE)
  )
}
