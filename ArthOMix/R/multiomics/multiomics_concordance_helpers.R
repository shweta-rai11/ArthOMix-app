## R/multiomics/multiomics_concordance_helpers.R
## Data-adaptive engine for the "Gene-CpG Concordance" submodule
## (mod_multi_concordance.R). Pure functions only, no Shiny reactives here.
##
## Nothing in this file re-implements logic that already exists elsewhere in
## the app - it composes real, already-verified building blocks:
##   - crossomics_integration_helpers.R (cx_*): sample pairing, gene-ID
##     harmonization, CpG->gene/region/island annotation, Hyper/Hypo x Up/
##     Down classification, evidence tiering, FDR.
##   - multiomics_integration_helpers.R / multiomics_biomarker_helpers.R
##     (mi_*/mb_*): DIABLO selection output, already published live by the
##     Biomarker Discovery submodule to multi_results$biomarker$df.
##   - snf_clustering_helpers.R (sfc_*): per-feature cluster association,
##     reused directly against the real SNF cluster vector published by the
##     Patient Stratification submodule to multi_results$stratification$clusters.
##   - multiomics_dataset_helpers.R (multi_live_*): sample-overlap primitives.
##   - R/methylomics/functions/qc.R::methyl_beta_to_mvalue() for beta->M conversion.
## Every function below fails soft (returns NA/"Not available"/ok=FALSE with
## a reason) rather than fabricating a gene, CpG, coordinate, or statistic -
## per the module's own "never invent" requirement.

## ---------------------------------------------------------------------------
## 0. Small shared vocab
## ---------------------------------------------------------------------------

MCC_REGION_PROMOTER <- c("TSS200", "TSS1500", "5'UTR")
MCC_REGION_BODY <- c("Body", "3'UTR", "1stExon", "ExonBnd")

MCC_DIRECTION_LEVELS <- c(
  "Up expression + Hypomethylation", "Down expression + Hypermethylation",
  "Up expression + Hypermethylation", "Down expression + Hypomethylation",
  "Weak/uncertain", "Not interpretable"
)

## ---------------------------------------------------------------------------
## 1. Layer selection (spec section 5-6) - thin wrapper over the same
## mb_select_blocks() the Biomarker Discovery submodule already uses to turn
## two arbitrarily-named `multi_dataset$layers` entries into fixed
## "Transcriptomics"/"Methylomics" roles. Candidate defaults are guessed the
## same way mod_multi_biomarker.R already guesses them (by omics_type
## metadata first, by name regex as a fallback) - never silently assumed.
## ---------------------------------------------------------------------------

mcc_layer_candidates <- function(multi_dataset, omics_type) {
  layers <- multi_dataset$layers %||% list()
  meta <- multi_dataset$layer_meta %||% list()
  if (length(layers) == 0) return(character(0))
  by_meta <- names(layers)[vapply(names(layers), function(nm) identical(meta[[nm]]$omics_type, omics_type), logical(1))]
  if (length(by_meta) > 0) return(by_meta)
  ## No layer_meta$omics_type recorded (e.g. metadata stripped on import) -
  ## fall back to a name-based guess only, and only for the two types this
  ## module needs; never guessed for other omics types.
  rx <- switch(omics_type, rnaseq = "transcript|rna|express|gene", methylation = "methyl|cpg|beta|meth", NULL)
  if (is.null(rx)) return(character(0))
  names(layers)[grepl(rx, names(layers), ignore.case = TRUE)]
}

mcc_default_layer <- function(candidates, layers) {
  if (length(candidates) == 0) return(NULL)
  candidates[1]
}

## ---------------------------------------------------------------------------
## 2. Data status panel (spec section 4) - inspects the ACTUAL active
## dataset/results, never assumes a field exists. `expr_layer`/`meth_layer`
## are the user's current layer picks (may be NULL before any pick is made).
## ---------------------------------------------------------------------------

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

## ---------------------------------------------------------------------------
## 3. Format detection (spec section 5-6) - report only, never silently coerce.
## ---------------------------------------------------------------------------

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

## ---------------------------------------------------------------------------
## 4. Sample matching (spec section 6) - thin wrapper, never re-implements
## the intersection logic in cx_detect_sample_pairing().
## ---------------------------------------------------------------------------

mcc_match_samples <- function(expr_mat, meth_mat) {
  if (is.null(expr_mat) || is.null(meth_mat)) return(list(ok = FALSE, error = "Both expression and methylation data are required for sample matching."))
  pairing <- cx_detect_sample_pairing(rownames(expr_mat), rownames(meth_mat))
  list(
    ok = isTRUE(pairing$paired), pairing = pairing,
    n_matched = pairing$n_common, n_removed_expr = pairing$n_expr - pairing$n_common,
    n_removed_meth = pairing$n_meth - pairing$n_common, common_samples = pairing$common_samples
  )
}

## ---------------------------------------------------------------------------
## 5. Candidate biomarker pool (spec section 7, 17-19) - unifies whatever the
## app has ACTUALLY already produced. "Joint Biomarker Discovery" and
## "DIABLO" both read multi_results$biomarker$df (the one live DIABLO/
## block.splsda signature this app computes, published by the Biomarker
## Discovery submodule) - see the module's own tooltip for why there is one
## engine behind both facets. "SNF" reads the real cluster vector published
## by Patient Stratification and re-derives per-feature cluster association
## via sfc_feature_ranking() (never re-runs SNF itself).
## ---------------------------------------------------------------------------

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
  ## Sex-Stratified DIABLO panel (Integration's "Sex-Stratified" tab,
  ## multiomics_sexstratified_engine.R::mss_run_stratified()) - folded in
  ## here so a user who ran that tab on an Active/uploaded dataset also gets
  ## DIABLO-panel cross-referencing here, not just from the freeform
  ## Biomarker Discovery signature above. Only the DIABLO engine's panel
  ## (component-1 selectVar() loadings) is comparable to the freeform
  ## engine's own loadings - the Random Forest panel (MeanDecreaseGini
  ## importance) is a different statistic and is never folded in.
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

## Cross-references the Preloaded path's own gene/CpG rows against the
## precomputed, per-sex DIABLO candidate-biomarker panel for the same cohort
## (Table40/44b, read by mod_multi_concordance.R via
## MCC_PRELOADED_DIABLO_PANEL + multi_read_registry_table() - never a live
## refit). `panel` must have columns sex, omics ("transcriptomics_gene" /
## "methylation_CpG"), feature (Ensembl gene ID or CpG probe ID -
## joined against `df$gene_id`/`df$cpg` by exact ID, not gene symbol),
## biomarker_status, loading. Joined per-sex when both sides have a sex
## column, so a female-only panel hit never flags a male row.
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

## Merges DIABLO/Joint + SNF + custom user-entered features into one long
## table, then narrows to `source_filter` (spec section 7's Biomarker Source
## options). Returns list(ok, df, note) - ok = FALSE ("No candidate
## biomarkers available") only when nothing at all is present.
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
  ## Feature identifier type distinguishes a candidate gene from a candidate
  ## CpG when `omics` metadata is missing (e.g. a bare custom entry) -
  ## never guessed beyond the regex already used for id-type detection.
  base$id_type <- vapply(base$feature, function(f) mcc_detect_id_type(f), character(1))
  base$omics[is.na(base$omics) & base$id_type == "Illumina CpG probe ID"] <- "Methylomics"
  base$omics[is.na(base$omics) & base$id_type %in% c("Ensembl Gene ID", "Entrez ID", "Gene symbol")] <- "Transcriptomics"

  list(ok = TRUE, df = base, note = NULL)
}

## Spec section 7's Biomarker Source dropdown, applied to the pool above.
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

## ---------------------------------------------------------------------------
## 6. Gene<->CpG mapping (spec section 9) - one row per (gene, CpG) actually
## annotated to that gene, restricted to CpGs actually present in the
## methylation layer when one is supplied. Never invents a CpG.
## ---------------------------------------------------------------------------

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

## ---------------------------------------------------------------------------
## 7. Direction classification (spec section 10-11) - reuses cx_classify()
## for the base Hyper/Hypo x Up/Down category, then layers a genomic-
## context-aware `canonical` flag on top (never a universal hard-coded
## inverse rule). Promoter/TSS CpGs: canonical = inverse direction
## (Hyper+Down / Hypo+Up). Gene-body CpGs: canonical = concordant direction
## (Hyper+Up / Hypo+Down) - matching the precomputed Table42/45's own
## "concordant (gene-body hypermethylation -> higher expression)" vocabulary,
## so live and preloaded results agree on the rule. Any other/unknown region
## bucket: the rule does not apply (`canonical = NA`), never guessed.
## ---------------------------------------------------------------------------

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

  ## Canonical/non-canonical is only a meaningful verdict once a pair has
  ## actually cleared both significance thresholds - a "Weak/uncertain" or
  ## "Not interpretable" pair has no determined direction to judge against
  ## the region rule, so it stays NA ("Not applicable") rather than being
  ## called "Non-canonical" by default (spec section 11: the rule must not
  ## be applied beyond what the data actually supports).
  sig_relevant <- has_stats & both_sig
  df$canonical <- NA
  df$canonical[promoter & sig_relevant] <- cat_chr[promoter & sig_relevant] %in% c("Hyper + Down", "Hypo + Up")
  df$canonical[body & sig_relevant] <- cat_chr[body & sig_relevant] %in% c("Hyper + Up", "Hypo + Down")
  df$canonical_label <- ifelse(is.na(df$canonical), "Not applicable", ifelse(df$canonical, "Canonical", "Non-canonical"))
  df
}

MCC_CANONICAL_RULE_TEXT <- "Canonical: inverse methylation-expression relationship at promoter/TSS CpGs, concordant at gene-body CpGs. Other regions: Not applicable."

## ---------------------------------------------------------------------------
## 7b. Multi-omics biomarker direction breakdown - `direction_classification`
## (above) only assigns a real Up/Down x Hyper/Hypo label once a pair clears
## BOTH this run's own expression and methylation significance thresholds;
## everything else collapses into "Weak/uncertain", which hides the actual
## direction for a pair that DIABLO/SNF/Joint Biomarker Discovery already
## flagged as a candidate on independent multi-omics evidence. This restricts
## to those already-flagged pairs and reports the raw sign of dbeta/log2fc
## (methylation_direction x expression_direction) regardless of this run's
## significance thresholds - never a new significance verdict, just the
## direction of change for a feature already selected by another method.
## ---------------------------------------------------------------------------

MCC_RAW_DIRECTION_LEVELS <- c(
  "Up expression + Hypomethylation", "Down expression + Hypermethylation",
  "Up expression + Hypermethylation", "Down expression + Hypomethylation"
)

## Adds `raw_direction` (the same 4-way Up/Down x Hyper/Hypo label as
## `direction_classification`, but from the raw sign of dbeta/log2fc alone -
## every row that has both values gets one, never gated by this run's
## significance thresholds), `raw_canonical_label` (the same promoter/
## gene-body region rule as `canonical_label`, applied to `raw_direction`
## instead of the significance-gated category - the un-gated equivalent of
## the precomputed pipeline's own `biological_pattern` column on the
## Preloaded path), and `biomarker_source` (which of DIABLO/SNF/Joint
## flagged this row, "" if none) to every row of `pairs_df`. Never filters
## rows - callers subset first if they only want biomarker-flagged pairs.
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

## ---------------------------------------------------------------------------
## 7c. Genome-wide response-driven significance (Table3/4, script 05) -
## real per-feature nominal p_response/fdr_response/logFC_response
## (genes) and delta_M_response (CpGs), joined by (sex, gene_id) / (sex,
## cpg) against `deg`/`dmp` (multi_read_registry_table() reads of the
## candidate-only lookups above). This is a materially different test family
## than mcc_classify_direction()'s own recompute (which BH-corrects only over
## the small candidate subset already in pairs_df, inflating its apparent
## multiple-testing burden) - use this for a "is this pair genuinely
## significant" verdict, not the FDR columns computed on pairs_df alone.
## ---------------------------------------------------------------------------

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

## Restricted to rows DIABLO/SNF/Joint Biomarker Discovery already flagged as
## a candidate - see mcc_add_raw_direction() above for why raw sign, not the
## significance-gated direction_classification, is used here.
mcc_biomarker_direction_table <- function(pairs_df) {
  need <- c("gene_symbol", "cpg", "methylation_direction", "expression_direction")
  if (!all(need %in% colnames(pairs_df))) stop("mcc_biomarker_direction_table: pairs_df missing required columns.")
  flag <- function(col) if (col %in% colnames(pairs_df)) pairs_df[[col]] %in% TRUE else rep(FALSE, nrow(pairs_df))
  is_biomarker <- flag("diablo") | flag("snf") | flag("joint")
  df <- pairs_df[is_biomarker, , drop = FALSE]
  if (nrow(df) == 0) return(df)
  mcc_add_raw_direction(df)
}

## ---------------------------------------------------------------------------
## 8. Sample-level Gene-CpG correlation (spec section 12) - a genuine
## per-pair primitive: unlike cx_gene_correlation() (which matches
## expr_mat/meth_mat by IDENTICAL rowname, i.e. gene-level aggregated
## methylation vs. gene expression), a Gene-CpG pair needs the individual
## CpG's own value vector against its mapped gene's expression vector - a
## different orientation contract, not a fork of that function. Reuses the
## same stats::cor.test() call and BH-FDR (cx_adjust_p()) it already uses.
## ---------------------------------------------------------------------------

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

## ---------------------------------------------------------------------------
## 9. Optional adjusted association / regression (spec section 13) - only
## uses covariates that actually exist in sample_meta; disabled (returns
## ok=FALSE) rather than silently run when N is too small for the model.
## ---------------------------------------------------------------------------

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

## ---------------------------------------------------------------------------
## 10. Sex-specific analysis (spec section 14) - column detection only;
## values/labels are never invented, only read from real sample_meta.
## Delegates to multi_sex_candidates()/multi_sex_groups()
## (multiomics_helpers.R), the shared implementation every live "stratify by
## sex" path in this module now uses (Integration/Biomarker's Sex-Stratified
## tab, Pathway's upload branch) - kept as thin mcc_-prefixed wrappers so
## this file's own call sites are unchanged.
## ---------------------------------------------------------------------------

mcc_sex_candidates <- function(sample_meta) multi_sex_candidates(sample_meta)

mcc_sex_groups <- function(sample_meta, sex_col, sample_ids) multi_sex_groups(sample_meta, sex_col, sample_ids)

## ---------------------------------------------------------------------------
## 11. Priority score (spec section 16) - interpretable, component-visible,
## never a black-box AI score. Every component is a 0-1 normalized real
## value or a real boolean; final label is always "Potential"/"Candidate",
## never "confirmed".
## ---------------------------------------------------------------------------

mcc_priority_score <- function(df) {
  norm01 <- function(x) { x <- abs(x); if (all(is.na(x))) return(rep(NA_real_, length(x))); rng <- range(x, na.rm = TRUE); if (diff(rng) == 0) return(ifelse(is.na(x), NA_real_, 0.5)); (x - rng[1]) / diff(rng) }
  expr_component <- norm01(df$log2fc) %||% rep(NA_real_, nrow(df))
  meth_component <- norm01(df$dbeta) %||% rep(NA_real_, nrow(df))
  fdr_component <- 1 - norm01(pmin(df$expr_fdr, df$meth_fdr, na.rm = FALSE))
  cor_component <- norm01(df$correlation_r)
  ## isTRUE() is scalar (identical(TRUE, x)) - always FALSE for a vector of
  ## length > 1, which previously collapsed ifelse()'s whole result to a
  ## single recycled constant instead of one value per row. is.na()/a plain
  ## logical vector drive ifelse()'s shape correctly here.
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

## ---------------------------------------------------------------------------
## 12. Two-group differential statistics (spec section 5-6's "raw/normalized
## matrices supplied -> allow the user to select the analysis design and
## calculate the required statistics" path). `multi_dataset$layers` only
## ever holds numeric matrices (never a pre-computed DE-results table - that
## upload shape isn't part of the Dataset Workspace's own ingestion
## contract, which this task does not modify), so this module's live branch
## computes log2FC/dbeta and a per-feature test itself, restricted to the
## already-matched samples and an already-existing 2-class design column
## from real sample_meta - never a fabricated grouping.
## ---------------------------------------------------------------------------

## Candidate 2-class categorical sample_meta columns among the matched
## samples only - reuses mi_outcome_summary()'s own categorical-detection
## logic rather than re-implementing class counting.
mcc_design_candidates <- function(sample_meta, sample_ids) {
  if (is.null(sample_meta) || ncol(sample_meta) == 0) return(character(0))
  ok <- vapply(colnames(sample_meta), function(cl) {
    s <- mi_outcome_summary(sample_meta, cl, sample_ids)
    !is.null(s) && identical(s$type, "categorical") && identical(s$n_classes, 2L)
  }, logical(1))
  colnames(sample_meta)[ok]
}

## Per-feature two-group log2FC (assumes already log2-scale expression, the
## standard convention for the normalized matrices this module receives) +
## Welch t-test + BH-FDR. Returns data.frame(feature, log2fc, p, fdr).
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

## Per-CpG delta-M (Welch t-test on M-values, preferred for the statistical
## test per spec section 6) + delta-Beta (for biological interpretation,
## reported alongside, never in place of, delta-M) + BH-FDR.
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

## ---------------------------------------------------------------------------
## 13. Summary counts (spec section 26)
## ---------------------------------------------------------------------------

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
