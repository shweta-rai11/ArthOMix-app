## R/crossomics/crossomics_biomarkerconv_helpers.R
## Pure data-processing logic for the "Biomarker Convergence" Cross-Omics
## sub-module (mod_cross_biomarker_conv.R) - a live, reconfigurable
## reimplementation of the pipeline's own
## cross_Omics_Sexstratified_COPY/scripts/01_cross_omics_eqtl_mqtl_biomarkers.R,
## which the Dataset tab already lets you browse the ONE-TIME precomputed
## output of (cross_omics_eQTL_mQTL_{female,male,combined}.csv). That script
## does no new statistics - it only JOINS the eQTL-MR gene panel and mQTL-MR
## CpG panel (each already causally screened upstream) against the raw
## DEG/DMP/DMR tables, then labels significance at fixed thresholds. This
## file reproduces that exact join, with every threshold the original script
## hardcoded exposed as a parameter instead (defaulting to the script's own
## value) - see cx_bc_build_join()'s `params` argument.
##
## Deliberately self-contained: reads directly from CX_DATA_ROOT's sibling
## Q2/Q3 result directories rather than reusing cx_load_default_deg()/
## load_default_dmp() (which point at the live Transcriptomics/Methylomics
## apps' own default dataset roots - not necessarily the same physical
## files), matching mod_cross_integration.R's own "independent of whatever
## is loaded elsewhere" principle. Every loader is fail-soft
## (list(ok, df, error)) - never errors, never fabricates.

## ---------------------------------------------------------------------------
## Paths (CX_DATA_ROOT-relative; global.R:439 already defines CX_DATA_ROOT /
## CX_DATA_AVAILABLE)
## ---------------------------------------------------------------------------

CX_BC_ROOT <- if (exists("CX_DATA_ROOT")) dirname(CX_DATA_ROOT) else NA_character_
CX_BC_Q2_DIR <- if (!is.na(CX_BC_ROOT)) file.path(CX_BC_ROOT, "Research_Q2_TRANSCRIPTOMICS_sexstratified_COPY", "results", "tables") else NA_character_
CX_BC_Q3_DIR <- if (!is.na(CX_BC_ROOT)) file.path(CX_BC_ROOT, "Research_Q3_METHYLOMICS_sexstratified_COPY", "methylomics") else NA_character_
CX_BC_DATA_AVAILABLE <- !is.na(CX_BC_Q2_DIR) && dir.exists(CX_BC_Q2_DIR) && !is.na(CX_BC_Q3_DIR) && dir.exists(CX_BC_Q3_DIR)

.cx_bc_read_csv <- function(path) {
  if (is.na(path) || !file.exists(path)) return(list(ok = FALSE, df = NULL, error = sprintf("File not found: %s", path)))
  df <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = sprintf("Could not read %s: %s", basename(path), conditionMessage(df))))
  list(ok = TRUE, df = df, error = NULL)
}

## ---------------------------------------------------------------------------
## Per-layer loaders
## ---------------------------------------------------------------------------

## eQTL-MR gene panel - already FDR<0.05-filtered upstream (this is the same
## MR_causal_FDR_{sex}.csv the report's "eQTL-MR panel" refers to throughout).
cx_bc_load_eqtl_panel <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  res <- .cx_bc_read_csv(file.path(CX_BC_Q2_DIR, sprintf("MR_causal_FDR_%s.csv", sex)))
  if (!res$ok) return(res)
  d <- res$df
  out <- data.frame(
    gene = d$gene, eQTL_MR_beta = d$b, eQTL_MR_pval = d$pval,
    eQTL_MR_OR = d$OR, eQTL_MR_OR_lo = d$OR_lo, eQTL_MR_OR_hi = d$OR_hi,
    eQTL_MR_FDR = d$FDR_pooled, eQTL_MR_direction = d$risk,
    eQTL_MHC_region = isTRUE(as.logical(d$MHC_gene)) | (as.character(d$MHC_gene) %in% c("TRUE", "true", "1")),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene) & nzchar(out$gene), , drop = FALSE]
  list(ok = TRUE, df = out, error = NULL)
}

## mQTL candidate panel: ensemble-vote-selected CpGs (n_votes >= vote_threshold)
## joined with their mQTL-MR estimate. Genes are resolved via the CpG->gene
## annotation table (cx_bc_load_cpg_annotation()) since neither source file
## carries a gene column directly.
cx_bc_load_mqtl_panel <- function(sex = c("female", "male"), vote_threshold = 2L) {
  sex <- match.arg(sex)
  votes <- .cx_bc_read_csv(file.path(CX_BC_Q3_DIR, "script07_ml_feature_selection", "tables", sprintf("ensemble_votes_%s.csv", sex)))
  if (!votes$ok) return(votes)
  est <- .cx_bc_read_csv(file.path(CX_BC_Q3_DIR, "script08_mendelian_randomization", "tables", sprintf("mr_estimates_%s.csv", sex)))
  if (!est$ok) return(est)
  v <- votes$df[votes$df$n_votes >= vote_threshold, , drop = FALSE]
  e <- est$df
  m <- merge(v[, c("cpg", "n_votes"), drop = FALSE], e[, c("exposure", "b", "se", "pval"), drop = FALSE],
             by.x = "cpg", by.y = "exposure", all.x = TRUE)
  colnames(m) <- c("cpg", "n_votes", "mQTL_MR_beta", "mQTL_MR_se", "mQTL_MR_pval")
  list(ok = TRUE, df = m, error = NULL)
}

## The pipeline's own pre-extracted CpG -> gene/chr/pos table (450K array) -
## simpler and more directly matches what 01_cross_omics_eqtl_mqtl_biomarkers.R
## itself used than re-deriving this from the Bioconductor annotation package.
## Cached per R process (same pattern as crossomics_integration_helpers.R's
## .cx_anno_cache).
.cx_bc_anno_cache <- new.env(parent = emptyenv())

cx_bc_load_cpg_annotation <- function() {
  cached <- .cx_bc_anno_cache[["anno"]]
  if (!is.null(cached)) return(list(ok = TRUE, df = cached, error = NULL))
  res <- .cx_bc_read_csv(file.path(CX_DATA_ROOT, "data", "cpg_gene_annotation_450k.csv"))
  if (!res$ok) return(res)
  .cx_bc_anno_cache[["anno"]] <- res$df
  list(ok = TRUE, df = res$df, error = NULL)
}

cx_bc_load_deg <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  res <- .cx_bc_read_csv(file.path(CX_BC_Q2_DIR, sprintf("DEG_%s_full.csv", sex)))
  if (!res$ok) return(res)
  d <- res$df
  out <- data.frame(gene = d$gene, DEG_logFC = d$logFC, DEG_pval = d$P.Value, DEG_adjP = d$adj.P.Val, stringsAsFactors = FALSE)
  out <- cx_bc_dedup_min(out, "gene", "DEG_adjP")
  list(ok = TRUE, df = out, error = NULL)
}

cx_bc_load_dmp <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  res <- .cx_bc_read_csv(file.path(CX_BC_Q3_DIR, "script03_dmp_sva_sexstratified", "tables", sprintf("dmp_%s_full.csv", sex)))
  if (!res$ok) return(res)
  d <- res$df
  ann <- cx_bc_load_cpg_annotation()
  if (!ann$ok) return(ann)
  m <- merge(d[, c("cpg", "dbeta", "p_bacon", "fdr_bacon"), drop = FALSE], ann$df[, c("cpg", "gene", "chr", "pos"), drop = FALSE], by = "cpg", all.x = TRUE)
  m <- m[!is.na(m$gene) & nzchar(m$gene), , drop = FALSE]
  list(ok = TRUE, df = m, error = NULL)
}

cx_bc_load_dmr <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  res <- .cx_bc_read_csv(file.path(CX_BC_Q3_DIR, "script04_dmr_sexstratified", "tables", sprintf("dmr_%s_full.csv", sex)))
  if (!res$ok) return(res)
  d <- res$df
  d$dmr_id <- sprintf("%s:%s-%s", d$seqnames, d$start, d$end)
  genes_split <- strsplit(as.character(d$overlapping.genes), ",\\s*")
  rows <- lapply(seq_len(nrow(d)), function(i) {
    g <- genes_split[[i]]
    g <- g[!is.na(g) & nzchar(g)]
    if (length(g) == 0) return(NULL)
    data.frame(gene = g, dmr_id = d$dmr_id[i], DMR_meandiff = d$meandiff[i], DMR_fdr = d$dmr_fdr[i],
               DMR_ncpgs = d$no.cpgs[i], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) return(list(ok = FALSE, df = NULL, error = "No genes could be parsed from the DMR overlapping-genes column."))
  list(ok = TRUE, df = out, error = NULL)
}

## Generic "keep one row per key, the one with the smallest order_col" dedup -
## used for DEG, DMP-per-gene, DMR-per-gene. NA order_col values sort last.
cx_bc_dedup_min <- function(df, key_col, order_col) {
  if (!any(duplicated(df[[key_col]]))) return(df)
  by_key <- split(seq_len(nrow(df)), df[[key_col]])
  idx <- vapply(by_key, function(ix) {
    ov <- df[[order_col]][ix]
    ix[if (all(is.na(ov))) 1L else which.min(ifelse(is.na(ov), Inf, ov))]
  }, integer(1))
  df[idx, , drop = FALSE]
}

## ---------------------------------------------------------------------------
## The join itself (spec: reproduces build_sex_table() from
## 01_cross_omics_eqtl_mqtl_biomarkers.R, thresholds made configurable)
## ---------------------------------------------------------------------------

CX_BC_DEFAULT_PARAMS <- list(
  vote_threshold = 2L, dmp_genomewide_fdr = 0.05, dmp_suggestive_p = 1e-4,
  mqtl_sig_basis = "nominal_p", mqtl_sig_cutoff = 0.05, dmr_fdr = 0.05
)

## Returns list(ok, df, error). `df` has one row per gene in
## union(eQTL-MR panel genes, mQTL-MR candidate panel genes), left-joined
## with DEG/DMP/DMR evidence and labeled at the given thresholds - never a
## new statistical test, exactly mirroring the original script's own scope.
cx_bc_build_join <- function(sex = c("female", "male"), params = CX_BC_DEFAULT_PARAMS) {
  sex <- match.arg(sex)
  if (!CX_BC_DATA_AVAILABLE) return(list(ok = FALSE, df = NULL, error = "Biomarker Convergence source data is not available in this deployment."))
  p <- utils::modifyList(CX_BC_DEFAULT_PARAMS, params %||% list())

  eqtl <- cx_bc_load_eqtl_panel(sex); if (!eqtl$ok) return(eqtl)
  mqtl <- cx_bc_load_mqtl_panel(sex, p$vote_threshold); if (!mqtl$ok) return(mqtl)
  ann <- cx_bc_load_cpg_annotation(); if (!ann$ok) return(ann)
  deg <- cx_bc_load_deg(sex); if (!deg$ok) return(deg)
  dmp <- cx_bc_load_dmp(sex); if (!dmp$ok) return(dmp)
  dmr <- cx_bc_load_dmr(sex); if (!dmr$ok) return(dmr)

  ## Resolve each mQTL candidate CpG to a gene, keep the best (lowest p) CpG
  ## per gene if a gene has more than one candidate.
  mqtl_g <- merge(mqtl$df, ann$df[, c("cpg", "gene"), drop = FALSE], by = "cpg", all.x = TRUE)
  mqtl_g <- mqtl_g[!is.na(mqtl_g$gene) & nzchar(mqtl_g$gene), , drop = FALSE]
  mqtl_g <- cx_bc_dedup_min(mqtl_g, "gene", "mQTL_MR_pval")

  all_genes <- union(eqtl$df$gene, mqtl_g$gene)
  if (length(all_genes) == 0) return(list(ok = FALSE, df = NULL, error = "No genes in either the eQTL-MR or mQTL-MR candidate panel for this sex."))
  out <- data.frame(sex = sex, gene = all_genes, stringsAsFactors = FALSE)
  out$in_eQTL_MR_panel <- out$gene %in% eqtl$df$gene
  out$in_mQTL_MR_panel <- out$gene %in% mqtl_g$gene
  out <- merge(out, eqtl$df, by = "gene", all.x = TRUE)
  out <- merge(out, mqtl_g[, c("gene", "cpg", "n_votes", "mQTL_MR_beta", "mQTL_MR_pval"), drop = FALSE],
               by = "gene", all.x = TRUE)
  colnames(out)[colnames(out) == "cpg"] <- "mQTL_candidate_cpg"
  out <- merge(out, deg$df, by = "gene", all.x = TRUE)

  ## DMP: the candidate CpG's own value if the gene entered via the mQTL
  ## panel, else the gene's most significant (min p_bacon) annotated CpG -
  ## the asymmetry the original script's header comment calls out explicitly.
  dmp_best_per_gene <- cx_bc_dedup_min(dmp$df, "gene", "p_bacon")
  dmp_by_candidate <- dmp$df[match(out$mQTL_candidate_cpg, dmp$df$cpg), , drop = FALSE]
  dmp_by_gene <- dmp_best_per_gene[match(out$gene, dmp_best_per_gene$gene), , drop = FALSE]
  use_candidate <- !is.na(out$mQTL_candidate_cpg)
  out$DMP_top_cpg <- ifelse(use_candidate, dmp_by_candidate$cpg, dmp_by_gene$cpg)
  out$DMP_dbeta <- ifelse(use_candidate, dmp_by_candidate$dbeta, dmp_by_gene$dbeta)
  out$DMP_fdr_bacon <- ifelse(use_candidate, dmp_by_candidate$fdr_bacon, dmp_by_gene$fdr_bacon)
  out$DMP_p_bacon <- ifelse(use_candidate, dmp_by_candidate$p_bacon, dmp_by_gene$p_bacon)

  ## A gene can overlap more than one DMR - keep only the most significant
  ## (min dmr_fdr) region per gene before the join, matching the original
  ## script's "the gene's most-significant overlapping region" rule; without
  ## this, genes with multiple overlapping DMRs would get duplicated rows.
  dmr_best_per_gene <- cx_bc_dedup_min(dmr$df, "gene", "DMR_fdr")
  out <- merge(out, dmr_best_per_gene[, c("gene", "dmr_id", "DMR_meandiff", "DMR_fdr", "DMR_ncpgs"), drop = FALSE], by = "gene", all.x = TRUE)

  ## Significance labels at the (reconfigurable) thresholds - no new test,
  ## same as the original script's fixed-threshold labeling.
  out$DEG_significant <- !is.na(out$DEG_adjP) & out$DEG_adjP < 0.05
  out$DEG_direction <- ifelse(is.na(out$DEG_logFC), NA_character_, ifelse(out$DEG_logFC > 0, "Up", "Down"))
  out$DMP_genomewide_significant <- !is.na(out$DMP_fdr_bacon) & out$DMP_fdr_bacon < p$dmp_genomewide_fdr
  out$DMP_suggestive <- !is.na(out$DMP_p_bacon) & out$DMP_p_bacon < p$dmp_suggestive_p
  out$DMP_direction <- ifelse(is.na(out$DMP_dbeta), NA_character_, ifelse(out$DMP_dbeta > 0, "Hyper", "Hypo"))
  out$DMR_significant <- !is.na(out$DMR_fdr) & out$DMR_fdr < p$dmr_fdr
  out$DMR_direction <- ifelse(is.na(out$DMR_meandiff), NA_character_, ifelse(out$DMR_meandiff > 0, "Hyper", "Hypo"))
  mqtl_stat <- if (identical(p$mqtl_sig_basis, "fdr")) stats::p.adjust(out$mQTL_MR_pval, method = "BH") else out$mQTL_MR_pval
  out$mQTL_MR_significant <- !is.na(mqtl_stat) & mqtl_stat < p$mqtl_sig_cutoff
  out$eQTL_MR_significant <- !is.na(out$eQTL_MR_FDR)  ## the panel is already FDR<0.05-filtered upstream

  out$methylation_significant <- out$DMP_genomewide_significant %in% TRUE | out$DMR_significant %in% TRUE
  out$n_evidence_layers <- rowSums(cbind(out$eQTL_MR_significant %in% TRUE, out$DEG_significant %in% TRUE,
                                          out$methylation_significant %in% TRUE, out$mQTL_MR_significant %in% TRUE))

  list(ok = TRUE, df = out, error = NULL, params = p)
}

cx_bc_build_provenance <- function(sex, params, run_at) {
  c(
    sprintf("Sex: %s", toupper(sex)),
    sprintf("eQTL-MR gene panel: MR_causal_FDR_%s.csv (already FDR<0.05-filtered upstream)", sex),
    sprintf("mQTL-MR candidate panel: ensemble_votes_%s.csv (n_votes >= %s) joined with mr_estimates_%s.csv", sex, params$vote_threshold, sex),
    sprintf("DEG source: DEG_%s_full.csv", sex),
    sprintf("DMP source: dmp_%s_full.csv (SVA/bacon-adjusted)", sex),
    sprintf("DMR source: dmr_%s_full.csv", sex),
    sprintf("DMP genome-wide-significant threshold: FDR < %s", params$dmp_genomewide_fdr),
    sprintf("DMP suggestive threshold: raw P < %s", params$dmp_suggestive_p),
    sprintf("mQTL-MR significance basis: %s < %s", if (identical(params$mqtl_sig_basis, "fdr")) "BH-FDR (recomputed over this run's test set)" else "nominal P", params$mqtl_sig_cutoff),
    sprintf("DMR significant threshold: FDR < %s", params$dmr_fdr),
    "This stage performs no new statistical test - it joins and labels already-computed per-layer results at the thresholds above.",
    sprintf("Run at: %s", run_at %||% "(not run yet)")
  )
}
