## R/crossomics/functions/integration/crossomics_integration_helpers.R
## Pure data-processing logic for the "Expression and Methylation" Cross-Omics
## sub-module (mod_cross_integration.R) - column auto-detection, gene-level
## methylation aggregation, quadrant classification, sample-pairing/
## correlation, and provenance. No Shiny dependency, so this is unit-testable
## on its own. Every function returns a fail-soft list(ok = TRUE/FALSE, ...)
## rather than throwing, the same sentinel convention as
## R/methylomics/parse_upload.R.

## ---------------------------------------------------------------------------
## Column auto-detection
## ---------------------------------------------------------------------------

CX_FIELD_PATTERNS <- list(
  gene = c("^gene[_ .]?symbol$", "^hgnc[_ .]?symbol$", "^symbol$", "^gene[_ .]?name$",
           "^gene$", "^genes$", "^gene[_ .]?id$", "^geneid$",
           "^ucsc[_ .]?refgene[_ .]?name$", "^overlapping[_ .]?genes?$"),
  cpg = c("^cpg[_ .]?id$", "^cpg$", "^probe[_ .]?id$", "^probeid$", "^probe$", "^illumina[_ .]?id$"),
  log2fc = c("^log2[_ .]?fc$", "^log2fc$", "^logfc$", "^log[_ .]?fold[_ .]?change$",
             "^fold[_ .]?change$", "^log2[_ .]?fold[_ .]?change$", "^lfc$"),
  dbeta = c("^delta[_ .]?beta$", "^dbeta$", "^d[_ .]?beta$", "^meandiff$", "^mean[_ .]?diff$",
            "^beta[_ .]?diff(erence)?$", "^methylation[_ .]?change$", "^methylation[_ .]?difference$"),
  beta = c("^beta$", "^beta[_ .]?value$", "^methylation[_ .]?beta$", "^avg[_ .]?beta$"),
  pvalue = c("^p[_ .]?value$", "^pval$", "^p[_ .]?val$", "^pvalue$", "^p\\.value$",
             "^stouffer$", "^p[_ .]?bacon$"),
  fdr = c("^fdr$", "^adj\\.?p\\.?val$", "^adjusted[_ .]?p[_ .]?value$", "^padj$",
          "^q[_ .]?value$", "^qval$", "^significance$", "^fdr[_ .]?bacon$",
          "^dmr[_ .]?fdr$", "^min[_ .]?smoothed[_ .]?fdr$"),
  chr = c("^chr(omosome)?$", "^seqnames$"),
  pos = c("^pos(ition)?$", "^start$", "^bp$", "^coordinate$"),
  end = c("^end$"),
  n_cpgs = c("^no\\.?cpgs?$", "^n[_ .]?cpgs?$", "^num[_ .]?cpgs?$", "^ncpg$"),
  region_id = c("^dmr[_ .]?id$", "^region[_ .]?id$"),
  region = c("^region$", "^feature$", "^genomic[_ .]?region$",
             "^ucsc[_ .]?refgene[_ .]?group$", "^annotation$"),
  island = c("^cpg[_ .]?island$", "^island$", "^relation[_ .]?to[_ .]?island$",
             "^relation[_ .]?to[_ .]?ucsc[_ .]?cpg[_ .]?island$", "^island[_ .]?context$",
             "^ucsc[_ .]?cpg[_ .]?island$"),
  sample_id = c("^sample[_ .]?id$", "^sample$", "^sample[_ .]?name$", "^id$")
)

CX_EXPRESSION_FIELDS <- c("gene", "log2fc", "pvalue", "fdr", "sample_id")
CX_METHYLATION_FIELDS <- c("cpg", "gene", "dbeta", "beta", "pvalue", "fdr", "chr", "pos",
                            "end", "n_cpgs", "region_id", "region", "island", "sample_id")

cx_match_column <- function(cols, patterns, exclude = character(0)) {
  cols_lower <- tolower(trimws(cols))
  candidates <- setdiff(cols, exclude)
  candidates_lower <- tolower(trimws(candidates))
  for (p in patterns) {
    hit <- grep(p, candidates_lower, perl = TRUE)
    if (length(hit) >= 1) return(candidates[hit[1]])
  }
  NA_character_
}

## Returns a named character vector, one entry per canonical field for `kind`,
## each either a detected column name from `df` or NA_character_. Detection
## is greedy/order-independent (a field already claimed by an earlier match
## is excluded from later candidates so e.g. "fdr" and "pvalue" don't both
## grab the same ambiguous column).
cx_detect_columns <- function(df, kind = c("expression", "methylation")) {
  kind <- match.arg(kind)
  fields <- if (kind == "expression") CX_EXPRESSION_FIELDS else CX_METHYLATION_FIELDS
  cols <- colnames(df)
  detected <- character(0)
  claimed <- character(0)
  for (f in fields) {
    hit <- cx_match_column(cols, CX_FIELD_PATTERNS[[f]], exclude = claimed)
    detected[f] <- hit
    if (!is.na(hit)) claimed <- c(claimed, hit)
  }
  detected
}

## Wide per-sample matrix detection: any numeric column not already claimed
## by a mapped metadata field is treated as a sample. Used to decide whether
## sample-level correlation is even possible (see cx_detect_sample_pairing()).
## Common per-gene/per-probe summary-statistic columns that are numeric but
## are never per-sample measurements - excluded so a plain DEG/DMP summary
## table (gene, logFC, AveExpr, t, P.Value, adj.P.Val, ...) doesn't get
## misread as a wide per-sample matrix.
CX_NON_SAMPLE_NUMERIC_NAMES <- c("aveexpr", "avgexpr", "t", "b", "z", "score", "statistic",
                                  "n", "n_probes", "nprobe", "logfc_m", "meandiff")

cx_detect_sample_columns <- function(df, mapping) {
  claimed <- unname(mapping[!is.na(mapping)])
  candidates <- setdiff(colnames(df), claimed)
  candidates <- candidates[!tolower(trimws(candidates)) %in% CX_NON_SAMPLE_NUMERIC_NAMES]
  is_num <- vapply(df[candidates], function(x) is.numeric(x) || (is.character(x) && suppressWarnings(!any(is.na(as.numeric(x[!is.na(x) & nzchar(x)]))))), logical(1))
  candidates[is_num]
}

## ---------------------------------------------------------------------------
## Standardization
## ---------------------------------------------------------------------------

cx_as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))

## Appends any source columns `mapping` didn't claim (e.g. AveExpr/t/dir on a
## DEG upload, or width/strand/HMFDR on a DMR upload) onto the standardized
## table, so nothing the user uploaded is silently dropped even though only
## the canonical fields drive the integration itself. A leftover column whose
## name collides with one of the standardized table's own reserved names is
## dropped rather than overwriting it. Called before any row-filtering/dedup
## so extras stay aligned with the rows they came from.
cx_bind_extra_columns <- function(out, df, mapping) {
  claimed <- unname(mapping[!is.na(mapping)])
  leftover <- setdiff(colnames(df), c(claimed, colnames(out)))
  if (length(leftover) == 0) return(out)
  cbind(out, df[, leftover, drop = FALSE])
}

## Builds the standardized gene-level expression table: gene, log2fc, pvalue,
## fdr. `mapping` is a named character vector as returned by
## cx_detect_columns()/user overrides (NA entries are simply skipped).
cx_standardize_expression <- function(df, mapping) {
  if (is.na(mapping["gene"])) return(list(ok = FALSE, error = "No Gene ID/Symbol column selected."))
  if (is.na(mapping["log2fc"])) return(list(ok = FALSE, error = "No log2FC/Fold Change column selected."))
  out <- data.frame(
    gene = trimws(as.character(df[[mapping["gene"]]])),
    log2fc = cx_as_numeric_safe(df[[mapping["log2fc"]]]),
    pvalue = if (!is.na(mapping["pvalue"])) cx_as_numeric_safe(df[[mapping["pvalue"]]]) else NA_real_,
    fdr = if (!is.na(mapping["fdr"])) cx_as_numeric_safe(df[[mapping["fdr"]]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- cx_bind_extra_columns(out, df, mapping)
  out <- out[nzchar(out$gene) & !is.na(out$gene), , drop = FALSE]
  if (nrow(out) == 0) return(list(ok = FALSE, error = "No valid gene rows after standardization."))
  out <- cx_dedup_by_gene(out)
  list(ok = TRUE, df = out)
}

## Multiple probes/transcripts/duplicated gene symbols -> one row per gene,
## keeping the row with the smallest FDR (falling back to smallest p-value,
## then the first row) - deterministic and transparent, same idea as the
## methylation side's "min FDR" aggregation option.
cx_dedup_by_gene <- function(df) {
  if (!any(duplicated(df$gene))) return(df)
  rows <- lapply(split(seq_len(nrow(df)), df$gene), function(idx) {
    sub <- df[idx, , drop = FALSE]
    pick <- if (!all(is.na(sub$fdr))) which.min(sub$fdr) else if (!all(is.na(sub$pvalue))) which.min(sub$pvalue) else 1L
    sub[pick, , drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

## Region string -> Promoter / Gene body / Other bucket. Covers both the
## Illumina UCSC_RefGene_Group vocabulary (TSS1500, TSS200, 5'UTR, 1stExon,
## Body, 3'UTR, ExonBnd) and free-text uploads ("promoter", "exon", "intron",
## "gene body", "upstream").
cx_region_bucket <- function(region_raw) {
  r <- tolower(trimws(as.character(region_raw)))
  dplyr::case_when(
    is.na(r) | !nzchar(r) ~ NA_character_,
    grepl("tss|promoter|5'utr|5utr|1stexon|upstream", r) ~ "Promoter",
    grepl("body|exon|intron|3'utr|3utr|exonbnd|gene[_ ]?body", r) ~ "Gene body",
    TRUE ~ "Other"
  )
}

## Fine-grained genomic-context token (spec section 12) - kept distinct from
## cx_region_bucket()'s coarse Promoter/Gene body/Other bucket so callers can
## filter/report/count at either granularity. Illumina's UCSC_RefGene_Group
## values are semicolon-joined (one gene can have several overlapping
## transcript annotations per probe) - this takes the first token, the same
## "first/representative annotation" convention cx_get_region_annotation()
## already uses for the `gene` column itself.
CX_REGION_FINE_VOCAB <- c("TSS200", "TSS1500", "5'UTR", "1stExon", "Body", "3'UTR", "ExonBnd")

cx_region_fine <- function(region_raw) {
  r <- as.character(region_raw)
  ## trimws() has real per-call overhead (argument matching, regex setup) -
  ## calling it once per element inside vapply (as this used to) cost ~16s
  ## on a 319k-row methylation table alone. Extract the first token per
  ## element with vapply (cheap - no trimws inside the loop), then trimws()
  ## ONCE, vectorized, over the whole result.
  first <- vapply(strsplit(r, "[;,]"), function(x) if (length(x) == 0) NA_character_ else x[1], character(1))
  first <- trimws(first)
  ifelse(is.na(first) | !nzchar(first), NA_character_, first)
}

CX_ISLAND_VOCAB <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")

cx_standardize_methylation <- function(df, mapping) {
  has_cpg <- !is.na(mapping["cpg"])
  has_gene <- !is.na(mapping["gene"])
  if (!has_gene) return(list(ok = FALSE, error = "No Gene Symbol column selected (required to join with expression)."))
  dbeta_col <- mapping["dbeta"]
  beta_col <- mapping["beta"]
  if (is.na(dbeta_col) && is.na(beta_col)) {
    return(list(ok = FALSE, error = "No Δβ (methylation change) or beta-value column selected."))
  }
  out <- data.frame(
    cpg = if (has_cpg) as.character(df[[mapping["cpg"]]]) else paste0("row", seq_len(nrow(df))),
    gene = trimws(as.character(df[[mapping["gene"]]])),
    dbeta = if (!is.na(dbeta_col)) cx_as_numeric_safe(df[[dbeta_col]]) else NA_real_,
    pvalue = if (!is.na(mapping["pvalue"])) cx_as_numeric_safe(df[[mapping["pvalue"]]]) else NA_real_,
    fdr = if (!is.na(mapping["fdr"])) cx_as_numeric_safe(df[[mapping["fdr"]]]) else NA_real_,
    chr = if (!is.na(mapping["chr"])) as.character(df[[mapping["chr"]]]) else NA_character_,
    pos = if (!is.na(mapping["pos"])) cx_as_numeric_safe(df[[mapping["pos"]]]) else NA_real_,
    ## end/n_cpgs/region_id are DMR-specific (region span, CpG count, region
    ## identifier) - absent (NA) for a DMP upload, present when auto-detected
    ## on a region-level file (e.g. DMRcate's end/no.cpgs columns).
    end = if (!is.na(mapping["end"])) cx_as_numeric_safe(df[[mapping["end"]]]) else NA_real_,
    n_cpgs = if (!is.na(mapping["n_cpgs"])) cx_as_numeric_safe(df[[mapping["n_cpgs"]]]) else NA_real_,
    region_id = if (!is.na(mapping["region_id"])) as.character(df[[mapping["region_id"]]]) else NA_character_,
    region_raw = if (!is.na(mapping["region"])) as.character(df[[mapping["region"]]]) else NA_character_,
    island_context = if (!is.na(mapping["island"])) as.character(df[[mapping["island"]]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out$region <- cx_region_bucket(out$region_raw)
  out$region_fine <- cx_region_fine(out$region_raw)
  out <- cx_bind_extra_columns(out, df, mapping)
  out <- out[nzchar(out$gene) & !is.na(out$gene), , drop = FALSE]
  if (nrow(out) == 0) return(list(ok = FALSE, error = "No valid gene rows after standardization."))
  list(ok = TRUE, df = out)
}

## ---------------------------------------------------------------------------
## Gene-level methylation aggregation (multiple CpGs/probes -> one row/gene)
## ---------------------------------------------------------------------------

CX_AGGREGATION_METHODS <- c(
  "mean" = "Mean Δβ across CpGs",
  "median" = "Median Δβ across CpGs",
  "min_fdr" = "CpG with the smallest FDR",
  "max_abs_dbeta" = "CpG with the largest |Δβ|",
  "promoter_only" = "Promoter CpGs only (mean)",
  "gene_body_only" = "Gene-body CpGs only (mean)",
  "promoter_and_body" = "Promoter + gene-body CpGs (mean)"
)

## ---------------------------------------------------------------------------
## CpG-level table (Level 2, spec section 9) - every individual CpG, with its
## own region/island context and significance flag, computed on the full
## standardized table BEFORE any gene-level aggregation collapses it. This is
## what "do not collapse CpGs for a gene without retaining their individual
## identities" (spec section 8) actually means in code: this table is never
## discarded, only summarized (see cx_cpg_counts_per_gene() below).
## ---------------------------------------------------------------------------

cx_cpg_level_table <- function(meth_std, meth_thresh, meth_fdr_thresh) {
  df <- meth_std
  df$methylation_direction <- ifelse(is.na(df$dbeta), NA_character_, ifelse(df$dbeta > 0, "Hyper", "Hypo"))
  df$sig_cpg <- !is.na(df$dbeta) & !is.na(df$fdr) & df$fdr < meth_fdr_thresh & abs(df$dbeta) >= meth_thresh
  keep_cols <- c("gene", "cpg", "region_raw", "region", "region_fine", "island_context",
                 "dbeta", "pvalue", "fdr", "methylation_direction", "sig_cpg", "chr", "pos")
  df[, intersect(keep_cols, colnames(df)), drop = FALSE]
}

## Per-gene CpG counting breakdown (spec section 21) - always computed from
## the full CpG-level table above (pre-aggregation), so it reflects every
## CpG mapped to the gene regardless of which single aggregation method the
## user picked for the gene-level dbeta/fdr summary. Vectorized (table()/
## tapply() rather than a per-gene lapply loop) - on this app's own real
## preloaded data (~19k genes, ~320k CpGs) the original per-gene-loop version
## took ~45s per run, too slow to sit behind a "Run Integration" click every
## time; this version does the same aggregation in well under a second.
cx_cpg_counts_per_gene <- function(cpg_level_df) {
  d <- cpg_level_df
  genes <- unique(d$gene)
  gene_f <- factor(d$gene, levels = genes)
  sig <- d$sig_cpg %in% TRUE
  hyper <- !is.na(d$dbeta) & d$dbeta > 0
  hypo <- !is.na(d$dbeta) & d$dbeta < 0

  region_tab <- table(gene_f, factor(d$region_fine, levels = CX_REGION_FINE_VOCAB))
  island_tab <- table(gene_f, factor(d$island_context, levels = CX_ISLAND_VOCAB))
  primary_idx <- max.col(region_tab, ties.method = "first")
  primary_region <- ifelse(rowSums(region_tab) == 0, NA_character_, colnames(region_tab)[primary_idx])

  data.frame(
    gene = genes,
    n_cpg_total = as.integer(table(gene_f)[genes]),
    n_cpg_unique = as.integer(tapply(d$cpg, gene_f, function(x) length(unique(x)))[genes]),
    n_cpg_significant = as.integer(tapply(sig, gene_f, sum)[genes]),
    n_cpg_hyper_sig = as.integer(tapply(sig & hyper, gene_f, sum)[genes]),
    n_cpg_hypo_sig = as.integer(tapply(sig & hypo, gene_f, sum)[genes]),
    n_cpg_hyper_all = as.integer(tapply(hyper, gene_f, sum)[genes]),
    n_cpg_hypo_all = as.integer(tapply(hypo, gene_f, sum)[genes]),
    n_tss200 = as.integer(region_tab[genes, "TSS200"]), n_tss1500 = as.integer(region_tab[genes, "TSS1500"]),
    n_5utr = as.integer(region_tab[genes, "5'UTR"]), n_1stexon = as.integer(region_tab[genes, "1stExon"]),
    n_body = as.integer(region_tab[genes, "Body"]), n_3utr = as.integer(region_tab[genes, "3'UTR"]),
    n_promoter_region = as.integer(tapply(d$region %in% "Promoter", gene_f, sum)[genes]),
    n_gene_body_region = as.integer(tapply(d$region %in% "Gene body", gene_f, sum)[genes]),
    n_island = as.integer(island_tab[genes, "Island"]),
    n_shore = as.integer(island_tab[genes, "N_Shore"]) + as.integer(island_tab[genes, "S_Shore"]),
    n_shelf = as.integer(island_tab[genes, "N_Shelf"]) + as.integer(island_tab[genes, "S_Shelf"]),
    n_open_sea = as.integer(island_tab[genes, "OpenSea"]),
    primary_region = primary_region,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

## `meth_std` is the output of cx_standardize_methylation()$df. Returns
## list(ok, df, note, cpg_level) where `df` has one row per gene (gene,
## dbeta, pvalue, fdr, chr, pos, region, n_probes, cpg - a
## representative/joined CpG list - plus the cx_cpg_counts_per_gene()
## breakdown columns), `cpg_level` is the full un-collapsed CpG-level table
## (spec section 9, Level 2), and `note` describes how the aggregation was
## performed, for on-screen transparency (spec section 4: "clearly show how
## the aggregation was performed").
cx_aggregate_methylation <- function(meth_std, method = "mean", meth_thresh = 0.10, meth_fdr_thresh = 0.05) {
  cpg_level <- cx_cpg_level_table(meth_std, meth_thresh, meth_fdr_thresh)
  df <- meth_std
  if (method %in% c("promoter_only", "gene_body_only", "promoter_and_body")) {
    keep <- switch(method,
      promoter_only = df$region == "Promoter",
      gene_body_only = df$region == "Gene body",
      promoter_and_body = df$region %in% c("Promoter", "Gene body")
    )
    keep[is.na(keep)] <- FALSE
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) {
      return(list(ok = FALSE, error = "No CpGs remain after filtering to the selected methylation region - check that your Region column/mapping is populated."))
    }
  }
  agg_fun <- switch(method,
    mean = , promoter_only = , gene_body_only = , promoter_and_body = "mean",
    median = "median",
    min_fdr = "min_fdr",
    max_abs_dbeta = "max_abs_dbeta"
  )
  by_gene <- split(seq_len(nrow(df)), df$gene)
  rows <- lapply(names(by_gene), function(g) {
    idx <- by_gene[[g]]
    sub <- df[idx, , drop = FALSE]
    pick <- switch(agg_fun,
      mean = ,
      median = NULL,
      min_fdr = which.min(ifelse(is.na(sub$fdr), Inf, sub$fdr)),
      max_abs_dbeta = which.max(ifelse(is.na(sub$dbeta), -Inf, abs(sub$dbeta)))
    )
    if (is.null(pick)) {
      f <- if (agg_fun == "mean") mean else stats::median
      dbeta_v <- f(sub$dbeta, na.rm = TRUE)
      ## dbeta_v aggregates ALL CpGs mapped to this gene, so its significance
      ## must be derived from that same set - not min(fdr)/min(pvalue), which
      ## can (and generically will) be driven by a single different CpG than
      ## the ones producing the mean effect. Combine directional per-CpG
      ## evidence with Stouffer's Z method: each p-value is turned into a
      ## signed z-score (sign = direction of that CpG's own dbeta), the
      ## z-scores are averaged (so CpGs disagreeing in direction partially
      ## cancel, exactly as they do in dbeta_v), and the combined z is
      ## converted back to a two-sided p-value. For a gene with a single CpG
      ## this reduces exactly to that CpG's own p-value. FDR is recomputed
      ## below via BH across all genes' combined p-values, since the gene
      ## (not the CpG) is now the analytical/testing unit.
      p_use <- pmin(pmax(sub$pvalue, .Machine$double.eps), 1 - .Machine$double.eps)
      dir_i <- ifelse(is.na(sub$dbeta) | sub$dbeta == 0, 1, sign(sub$dbeta))
      z_i <- dir_i * stats::qnorm(1 - p_use / 2)
      ok_i <- is.finite(z_i)
      pvalue_v <- if (sum(ok_i) == 0) NA_real_ else 2 * stats::pnorm(-abs(sum(z_i[ok_i]) / sqrt(sum(ok_i))))
      fdr_v <- NA_real_ ## placeholder - replaced by a gene-level BH pass after all genes are aggregated
      cpg_v <- paste(sub$cpg, collapse = ";")
      chr_v <- sub$chr[1]; pos_v <- sub$pos[1]; region_v <- sub$region[1]
    } else {
      dbeta_v <- sub$dbeta[pick]; fdr_v <- sub$fdr[pick]; pvalue_v <- sub$pvalue[pick]
      cpg_v <- sub$cpg[pick]; chr_v <- sub$chr[pick]; pos_v <- sub$pos[pick]; region_v <- sub$region[pick]
    }
    data.frame(gene = g, dbeta = dbeta_v, pvalue = pvalue_v, fdr = fdr_v,
               chr = chr_v, pos = pos_v, region = region_v %||% NA_character_,
               n_probes = length(idx), cpg = cpg_v,
               dbeta_mean = mean(sub$dbeta, na.rm = TRUE), dbeta_median = stats::median(sub$dbeta, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (agg_fun %in% c("mean", "median")) {
    ## Gene-level BH-FDR across the combined per-gene p-values computed above -
    ## the correct multiple-testing universe here is genes, not CpGs, since
    ## each gene's p-value already summarizes every CpG mapped to it.
    out$fdr <- stats::p.adjust(out$pvalue, method = "BH")
  }
  out$dbeta[is.infinite(out$dbeta)] <- NA_real_
  out$fdr[is.infinite(out$fdr)] <- NA_real_
  out$pvalue[is.infinite(out$pvalue)] <- NA_real_
  out$dbeta_mean[is.nan(out$dbeta_mean) | is.infinite(out$dbeta_mean)] <- NA_real_
  out$dbeta_median[is.nan(out$dbeta_median) | is.infinite(out$dbeta_median)] <- NA_real_

  counts <- cx_cpg_counts_per_gene(cpg_level)
  out <- merge(out, counts, by = "gene", all.x = TRUE)

  note <- sprintf(
    "%s (%s CpGs across %s genes; %s gene(s) had more than one CpG mapped).",
    CX_AGGREGATION_METHODS[[method]] %||% method,
    format(nrow(df), big.mark = ","), format(nrow(out), big.mark = ","),
    format(sum(out$n_probes > 1), big.mark = ",")
  )
  list(ok = TRUE, df = out, note = note, cpg_level = cpg_level)
}

## ---------------------------------------------------------------------------
## Classification
## ---------------------------------------------------------------------------

CX_CATEGORY_LABELS <- c(
  "Hyper + Down" = "Potential methylation-associated repression",
  "Hypo + Up" = "Potential methylation-associated activation",
  "Hyper + Up" = "Concordant-direction / potentially noncanonical association",
  "Hypo + Down" = "Concordant-direction / potentially noncanonical association",
  "Not significant" = "No integrated signal at the current thresholds"
)

CX_CATEGORY_ORDER <- c("Hyper + Down", "Hypo + Up", "Hyper + Up", "Hypo + Down", "Not significant")

## Deliberately hand-assigned rather than arthomix_pair(CX_CATEGORY_ORDER) -
## the two "potential regulatory" categories get the strongest (red/blue)
## hues so they visually lead, and "Not significant" is muted rather than a
## saturated series color.
CX_CATEGORY_COLORS <- c(
  "Hyper + Down" = ARTHOMIX_COLORS$red,
  "Hypo + Up" = ARTHOMIX_COLORS$blue,
  "Hyper + Up" = ARTHOMIX_COLORS$orange,
  "Hypo + Down" = ARTHOMIX_COLORS$violet,
  "Not significant" = ARTHOMIX_COLORS$ink_muted
)

CX_FILTER_CHOICES <- c(
  "All genes" = "All",
  "Hyper + Down" = "Hyper + Down", "Hypo + Up" = "Hypo + Up",
  "Hyper + Up" = "Hyper + Up", "Hypo + Down" = "Hypo + Down",
  "Not significant" = "Not significant",
  "Significant in both" = "sig_both",
  "Significant expression only" = "sig_expr_only",
  "Significant methylation only" = "sig_meth_only",
  "Significant correlation only" = "sig_cor_only"
)

## `joined` must have columns log2fc, expr_fdr, dbeta, meth_fdr. Adds
## sig_expression, sig_methylation, category, category_label,
## expression_direction, methylation_direction (spec sections 22/23 - derived
## only from the sign of the already-provided log2fc/dbeta, i.e. whatever
## comparison direction the source contrast used; never silently reversed).
## Never claims causality - category_label is always framed as
## "potential"/"association".
cx_classify <- function(joined, expr_thresh, expr_fdr_thresh, meth_thresh, meth_fdr_thresh) {
  df <- joined
  df$sig_expression <- !is.na(df$log2fc) & !is.na(df$expr_fdr) &
    df$expr_fdr < expr_fdr_thresh & abs(df$log2fc) >= expr_thresh
  df$sig_methylation <- !is.na(df$dbeta) & !is.na(df$meth_fdr) &
    df$meth_fdr < meth_fdr_thresh & abs(df$dbeta) >= meth_thresh
  df$expression_direction <- ifelse(is.na(df$log2fc), NA_character_, ifelse(df$log2fc > 0, "Upregulated", "Downregulated"))
  df$methylation_direction <- ifelse(is.na(df$dbeta), NA_character_, ifelse(df$dbeta > 0, "Hyper", "Hypo"))
  both <- df$sig_expression & df$sig_methylation
  df$category <- "Not significant"
  df$category[both & df$dbeta > 0 & df$log2fc < 0] <- "Hyper + Down"
  df$category[both & df$dbeta < 0 & df$log2fc > 0] <- "Hypo + Up"
  df$category[both & df$dbeta > 0 & df$log2fc > 0] <- "Hyper + Up"
  df$category[both & df$dbeta < 0 & df$log2fc < 0] <- "Hypo + Down"
  df$category <- factor(df$category, levels = CX_CATEGORY_ORDER)
  df$category_label <- CX_CATEGORY_LABELS[as.character(df$category)]
  df
}

## ---------------------------------------------------------------------------
## Evidence-Level classification (spec section 19) - explicit and separate
## from the Hyper/Hypo x Up/Down `category` column above. "Strong candidate"
## additionally requires a significant NEGATIVE sample-level correlation, and
## ONLY when correlation was actually computed (has_correlation = TRUE, i.e.
## cx_detect_sample_pairing() found real matched samples) - it is never
## inferred or assumed when correlation wasn't computed, per spec section 11.
## ---------------------------------------------------------------------------

CX_EVIDENCE_LEVELS <- c("Strong candidate", "Moderate candidate", "Expression-only",
                        "Methylation-only", "Discordant", "Insufficient evidence")

CX_EVIDENCE_DESCRIPTIONS <- c(
  "Strong candidate" = "Significant in both layers, inverse direction, and a significant negative sample-level methylation-expression correlation.",
  "Moderate candidate" = "Significant in both layers with an inverse direction, but sample-level correlation is unavailable or not significant.",
  "Expression-only" = "Significant in expression only - no significant methylation evidence at the current thresholds.",
  "Methylation-only" = "Significant in methylation only - no significant expression evidence at the current thresholds.",
  "Discordant" = "Significant in both layers, but the direction does not fit the simple inverse (methylation-expression) model.",
  "Insufficient evidence" = "Not significant in either layer at the current thresholds, or a required value is missing."
)

cx_classify_evidence <- function(classified_df, has_correlation = FALSE) {
  df <- classified_df
  both_sig <- df$sig_expression & df$sig_methylation
  inverse_dir <- (!is.na(df$dbeta) & !is.na(df$log2fc)) & ((df$dbeta > 0 & df$log2fc < 0) | (df$dbeta < 0 & df$log2fc > 0))
  neg_cor_sig <- if (isTRUE(has_correlation) && "correlation_r" %in% colnames(df)) {
    !is.na(df$correlation_r) & df$correlation_r < 0 & !is.na(df$correlation_fdr) & df$correlation_fdr < 0.05
  } else {
    rep(FALSE, nrow(df))
  }
  ev <- rep("Insufficient evidence", nrow(df))
  ev[df$sig_expression %in% TRUE & !(df$sig_methylation %in% TRUE)] <- "Expression-only"
  ev[!(df$sig_expression %in% TRUE) & df$sig_methylation %in% TRUE] <- "Methylation-only"
  ev[both_sig %in% TRUE & inverse_dir %in% TRUE] <- "Moderate candidate"
  ev[both_sig %in% TRUE & !(inverse_dir %in% TRUE)] <- "Discordant"
  ev[both_sig %in% TRUE & inverse_dir %in% TRUE & neg_cor_sig] <- "Strong candidate"
  factor(ev, levels = CX_EVIDENCE_LEVELS)
}

## Applies one of CX_FILTER_CHOICES' values to an already-classified/
## -correlated integrated table. "sig_cor_only" requires a `correlation_fdr`
## column (added post-correlation) and a `min_abs_r`/`max_fdr` pair; when
## correlation hasn't been computed it simply returns zero rows rather than
## erroring, since the UI already hides that filter option in that case.
cx_filter_by_category <- function(df, filter_key, min_abs_r = 0.3, max_cor_fdr = 0.05) {
  switch(filter_key,
    "All" = df,
    "sig_both" = df[df$sig_expression & df$sig_methylation, , drop = FALSE],
    "sig_expr_only" = df[df$sig_expression & !df$sig_methylation, , drop = FALSE],
    "sig_meth_only" = df[!df$sig_expression & df$sig_methylation, , drop = FALSE],
    "sig_cor_only" = if ("correlation_r" %in% colnames(df)) {
      df[!is.na(df$correlation_fdr) & df$correlation_fdr < max_cor_fdr & abs(df$correlation_r) >= min_abs_r, , drop = FALSE]
    } else df[0, , drop = FALSE],
    df[as.character(df$category) == filter_key, , drop = FALSE]
  )
}

## Advanced filters (spec section 20) - genomic region, CpG island status,
## Evidence Level, and minimum CpG count. Each is a no-op (returns df
## unchanged) when its selector is empty/NULL, or when the required column
## isn't present (e.g. filtering hasn't been run yet) - never errors.
cx_num0 <- function(x) ifelse(is.na(x), 0, x)

cx_filter_by_region <- function(df, regions) {
  if (is.null(regions) || length(regions) == 0 || !"primary_region" %in% colnames(df)) return(df)
  df[df$primary_region %in% regions, , drop = FALSE]
}

cx_filter_by_island <- function(df, islands) {
  if (is.null(islands) || length(islands) == 0 || !"n_island" %in% colnames(df)) return(df)
  keep <- rep(FALSE, nrow(df))
  if ("Island" %in% islands) keep <- keep | cx_num0(df$n_island) > 0
  if (any(c("N_Shore", "S_Shore") %in% islands)) keep <- keep | cx_num0(df$n_shore) > 0
  if (any(c("N_Shelf", "S_Shelf") %in% islands)) keep <- keep | cx_num0(df$n_shelf) > 0
  if ("OpenSea" %in% islands) keep <- keep | cx_num0(df$n_open_sea) > 0
  df[keep, , drop = FALSE]
}

cx_filter_by_evidence <- function(df, levels_selected) {
  if (is.null(levels_selected) || length(levels_selected) == 0 || !"evidence_level" %in% colnames(df)) return(df)
  df[as.character(df$evidence_level) %in% levels_selected, , drop = FALSE]
}

cx_filter_by_min_cpg <- function(df, min_cpg) {
  if (is.null(min_cpg) || min_cpg <= 0 || !"n_cpg_total" %in% colnames(df)) return(df)
  df[!is.na(df$n_cpg_total) & df$n_cpg_total >= min_cpg, , drop = FALSE]
}

cx_filter_by_correlation_direction <- function(df, direction) {
  if (is.null(direction) || identical(direction, "any") || !"correlation_r" %in% colnames(df)) return(df)
  if (identical(direction, "pos")) return(df[!is.na(df$correlation_r) & df$correlation_r > 0, , drop = FALSE])
  if (identical(direction, "neg")) return(df[!is.na(df$correlation_r) & df$correlation_r < 0, , drop = FALSE])
  df
}

## Collapses a probe/CpG-or-gene-level wide matrix (id column + numeric
## sample columns) to one row per unique id, averaging duplicates - used to
## build genes x samples matrices for correlation from either an
## already-gene-level upload or a CpG-level one (post gene mapping).
cx_build_gene_sample_matrix <- function(df, id_col, sample_cols) {
  ids <- as.character(df[[id_col]])
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  if (!any(duplicated(ids))) return(mat)
  out <- rowsum(mat, group = ids, na.rm = TRUE)
  counts <- as.vector(table(ids)[rownames(out)])
  out / counts
}

## ---------------------------------------------------------------------------
## Sample pairing / correlation (spec section 29 - the statistical safeguard)
## ---------------------------------------------------------------------------

## Returns list(paired, common_samples, n_expr, n_meth, n_common) - `paired`
## is TRUE only when both sides expose per-sample columns AND at least
## `min_overlap` sample identifiers are shared between them. Deliberately
## conservative: with no sample-level data at all (e.g. both sides are
## gene-level DEG/DMP summary tables, the normal shape for "Preloaded" mode),
## this always reports paired = FALSE, since there is no per-sample matrix to
## pair in the first place.
cx_detect_sample_pairing <- function(expr_samples, meth_samples, min_overlap = 3L) {
  common <- intersect(expr_samples, meth_samples)
  list(
    paired = length(expr_samples) > 0 && length(meth_samples) > 0 && length(common) >= min_overlap,
    common_samples = common,
    n_expr = length(expr_samples), n_meth = length(meth_samples), n_common = length(common)
  )
}

## Per-gene correlation between an expression value vector and a methylation
## (beta) value vector, both already restricted to the same `common_samples`
## order, for every gene present in both `expr_mat` and `meth_mat` (matrices,
## genes x samples). Returns a data.frame(gene, r, p, n) - FDR is added by the
## caller via cx_adjust_p() once the full gene list is assembled.
cx_gene_correlation <- function(expr_mat, meth_mat, common_samples, method = c("pearson", "spearman"), min_n = 3L) {
  method <- match.arg(method)
  genes <- intersect(rownames(expr_mat), rownames(meth_mat))
  if (length(genes) == 0 || length(common_samples) < min_n) {
    return(list(ok = FALSE, error = "Not enough matched genes/samples to compute correlation."))
  }
  em <- expr_mat[genes, common_samples, drop = FALSE]
  mm <- meth_mat[genes, common_samples, drop = FALSE]
  rows <- lapply(genes, function(g) {
    x <- as.numeric(em[g, ]); y <- as.numeric(mm[g, ])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_n) return(data.frame(gene = g, r = NA_real_, p = NA_real_, n = sum(ok)))
    ct <- tryCatch(stats::cor.test(x[ok], y[ok], method = method), error = function(e) NULL)
    if (is.null(ct)) return(data.frame(gene = g, r = NA_real_, p = NA_real_, n = sum(ok)))
    data.frame(gene = g, r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
  })
  list(ok = TRUE, df = do.call(rbind, rows))
}

cx_adjust_p <- function(p, method = c("BH", "bonferroni")) {
  method <- match.arg(method)
  stats::p.adjust(p, method = if (method == "BH") "BH" else "bonferroni")
}

## ---------------------------------------------------------------------------
## CpG -> gene/region annotation for the Preloaded pathway (450K/EPIC)
## ---------------------------------------------------------------------------
## Deliberately independent of R/methylomics/annotation.R::methyl_get_annotation()
## - that function's `anno` object doesn't expose UCSC_RefGene_Group (only the
## gene symbol), which this module needs for promoter/gene-body aggregation.
## Reads the same underlying Bioconductor data objects directly (Locations +
## Other), the same technique and the same reason (avoiding
## library()-attaching the annotation package - see that file's own comment).
## Cached per array type for the life of the R process.

CX_METH_ANNOTATION_PACKAGES <- list("450K" = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                                     "EPIC" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

.cx_anno_cache <- new.env(parent = emptyenv())

cx_get_region_annotation <- function(array_type = "450K") {
  pkg <- CX_METH_ANNOTATION_PACKAGES[[array_type]]
  if (is.null(pkg)) return(list(ok = FALSE, anno = NULL, reason = sprintf("No annotation package configured for %s.", array_type)))
  cached <- .cx_anno_cache[[pkg]]
  if (!is.null(cached)) return(list(ok = TRUE, anno = cached, reason = NULL))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf("The %s annotation package is not installed in this deployment.", pkg)))
  }
  anno <- tryCatch({
    e <- new.env(parent = emptyenv())
    ## Islands.UCSC (Relation_to_Island: Island/N_Shore/S_Shore/N_Shelf/
    ## S_Shelf/OpenSea) is loaded alongside Locations/Other for CpG-island
    ## context (spec section 12) - same object, same row order (one row per
    ## probe) as the other two, confirmed against this package's actual
    ## column names before writing this.
    utils::data(list = c("Locations", "Other", "Islands.UCSC"), package = pkg, envir = e)
    loc <- as.data.frame(e$Locations)
    oth <- as.data.frame(e$Other)
    isl <- as.data.frame(e$`Islands.UCSC`)
    ids <- rownames(loc)
    gene <- vapply(strsplit(oth[ids, "UCSC_RefGene_Name"], ";"), function(g) if (length(g) == 0 || !nzchar(g[1])) NA_character_ else g[1], character(1))
    group <- vapply(strsplit(oth[ids, "UCSC_RefGene_Group"], ";"), function(g) if (length(g) == 0 || !nzchar(g[1])) NA_character_ else g[1], character(1))
    island_ctx <- as.character(isl[ids, "Relation_to_Island"])
    data.frame(chr = loc[ids, "chr"], pos = loc[ids, "pos"], gene = gene, region_raw = group,
               island_context = island_ctx, row.names = ids, stringsAsFactors = FALSE)
  }, error = function(e) e)
  if (inherits(anno, "error")) return(list(ok = FALSE, anno = NULL, reason = sprintf("Could not extract annotation from %s: %s", pkg, conditionMessage(anno))))
  .cx_anno_cache[[pkg]] <- anno
  list(ok = TRUE, anno = anno, reason = NULL)
}

## ---------------------------------------------------------------------------
## Preloaded transcriptomics DEG loader (the counterpart to global.R's
## load_default_dmp()/load_default_dmr() for methylomics)
## ---------------------------------------------------------------------------

CX_DEG_TABLE_DIR <- file.path(DATA_ROOT, "results", "tables")

## Reads the published sex-stratified (or combined) differential-expression
## table (gene, logFC, AveExpr, t, P.Value, adj.P.Val, sig, dir) - the same
## graceful-degradation contract as load_default_dmp(): NULL, never an error,
## when the file isn't present.
cx_load_default_deg <- function(sex = c("female", "male", "all")) {
  sex <- match.arg(sex)
  path <- file.path(CX_DEG_TABLE_DIR, sprintf("DEG_%s_full.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Preloaded-methylation counterpart: reads global.R's own
## load_default_dmp("sva", sex) (SVA-adjusted, bacon-corrected - the stage
## the project's own comments say resolves the "plain" stage's inflation),
## then annotates cpg -> gene/chr/pos/region via cx_get_region_annotation()
## since that precomputed table (cpg, logFC_M, dbeta, t, p_raw, p_bacon,
## fdr_bacon) carries no gene column at all. Fail-soft: list(ok, df, error).
cx_load_default_methylation <- function(sex = c("female", "male", "all"), array_type = "450K") {
  sex <- match.arg(sex)
  ## sex="all" reads dmp_all_full.csv (script03_dmp_sva_sexstratified/tables/),
  ## a genuine pooled fit - group + sex + age + smoking + cell-type covariates
  ## + surrogate variables from sva::sva(), bacon-corrected - not a merge of
  ## the female/male tables. This mirrors cx_load_default_deg()'s "all", but
  ## note the underlying pipeline's own methods notes
  ## (METHODS_dmp_sexstratified.md, section 2.Z.1) chose sex-stratified models
  ## as the primary analysis specifically because a pooled fit can obscure
  ## sex-specific methylation effects (Tesfaye et al. 2024); this pooled panel
  ## is a genuine, separately-fitted complement to the Female/Male panels for
  ## workflows (like this one) that need a single sex-agnostic table, not a
  ## replacement for the stratified analysis or a claim that it's superior.
  if (!METH_DATA_AVAILABLE) {
    return(list(ok = FALSE, df = NULL, error = "Methylomics preloaded data is not available in this deployment."))
  }
  dmp <- load_default_dmp(stage = "sva", sex = sex)
  if (is.null(dmp)) return(list(ok = FALSE, df = NULL, error = "Could not read the preloaded methylation (DMP) table for this sex stratum."))
  ar <- cx_get_region_annotation(array_type)
  if (!isTRUE(ar$ok)) return(list(ok = FALSE, df = NULL, error = ar$reason))
  a <- ar$anno
  ## One match() against the ~485k-row annotation table's rownames, reused
  ## for every column - rowname-based data.frame indexing (a[ids, "col"])
  ## re-runs its own match() internally on every call, so doing this
  ## separately per column (as this used to) cost ~15-20s on the male DMP
  ## table alone; a single match() + integer indexing is a few seconds.
  idx <- match(dmp$cpg, rownames(a))
  df <- data.frame(
    cpg = dmp$cpg, gene = a$gene[idx], chr = a$chr[idx], pos = a$pos[idx],
    region_raw = a$region_raw[idx], island_context = a$island_context[idx],
    dbeta = dmp$dbeta, pvalue = dmp$p_bacon, fdr = dmp$fdr_bacon, stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "No CpGs in the preloaded methylation table could be annotated to a gene symbol."))
  mapping <- c(cpg = "cpg", gene = "gene", dbeta = "dbeta", beta = NA_character_,
               pvalue = "pvalue", fdr = "fdr", chr = "chr", pos = "pos", region = "region_raw",
               island = "island_context", sample_id = NA_character_)
  std <- cx_standardize_methylation(df, mapping)
  if (!std$ok) return(list(ok = FALSE, df = NULL, error = std$error))
  list(ok = TRUE, df = std$df, error = NULL)
}

## Region-level (DMR) counterpart to cx_load_default_methylation() above -
## same contract, but reads the pipeline's own already-called differentially
## methylated regions (DMRcate output, METH_DMR_DIR/dmr_{sex}_full.csv)
## directly. No cx_get_region_annotation() lookup needed here: this table
## already carries its own gene annotation (`overlapping.genes`), unlike the
## per-CpG DMP table. sex="all" reads dmr_all_full.csv, DMRcate region
## calling on the pooled per-CpG statistics from cx_load_default_methylation()'s
## sex="all" - see that function's own comment. `dbeta` here is a
## region's mean Δβ across its constituent CpGs (`meandiff`), `pvalue` is
## the pre-FDR Stouffer combined-probability statistic, and `fdr` is the
## already-BH-corrected region-level FDR (`dmr_fdr`) - see
## METHODS_dmr_sexstratified.md section 2.BB.4 for how these were computed.
cx_load_default_dmr <- function(sex = c("female", "male", "all")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) {
    return(list(ok = FALSE, df = NULL, error = "Methylomics preloaded data is not available in this deployment."))
  }
  path <- file.path(METH_DMR_DIR, sprintf("dmr_%s_full.csv", sex))
  if (!file.exists(path)) return(list(ok = FALSE, df = NULL, error = "Could not read the preloaded methylation (DMR) table for this sex stratum."))
  dmr <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) e)
  if (inherits(dmr, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the preloaded DMR table:", conditionMessage(dmr))))
  df <- data.frame(
    cpg = sprintf("%s:%s-%s", dmr$seqnames, dmr$start, dmr$end),
    gene = trimws(as.character(dmr$overlapping.genes)),
    dbeta = dmr$meandiff, pvalue = dmr$Stouffer, fdr = dmr$dmr_fdr,
    chr = as.character(dmr$seqnames), pos = dmr$start,
    region_raw = NA_character_, island_context = NA_character_,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "No DMRs in the preloaded region table are annotated to a gene symbol."))
  mapping <- c(cpg = "cpg", gene = "gene", dbeta = "dbeta", beta = NA_character_,
               pvalue = "pvalue", fdr = "fdr", chr = "chr", pos = "pos", region = "region_raw",
               island = "island_context", sample_id = NA_character_)
  std <- cx_standardize_methylation(df, mapping)
  if (!std$ok) return(list(ok = FALSE, df = NULL, error = std$error))
  list(ok = TRUE, df = std$df, error = NULL)
}

## ---------------------------------------------------------------------------
## Provenance
## ---------------------------------------------------------------------------

CX_MODULE_VERSION <- "ArthOMix Cross-Omics Expression and Methylation v1.1"

cx_build_provenance <- function(params) {
  c(
    sprintf("Sex stratum: %s", params$sex_stratum %||% "(not set)"),
    sprintf("Input mode: %s", params$input_mode %||% "(not set)"),
    sprintf("Transcriptomics source: %s", params$expr_source %||% "(not set)"),
    sprintf("Methylomics source: %s", params$meth_source %||% "(not set)"),
    sprintf("Expression threshold: |log2FC| ≥ %s, FDR < %s", params$expr_thresh, params$expr_fdr_thresh),
    sprintf("Methylation threshold: |Δβ| ≥ %s, FDR < %s", params$meth_thresh, params$meth_fdr_thresh),
    sprintf("Methylation region: %s", params$meth_region %||% "All"),
    sprintf("Aggregation method: %s", CX_AGGREGATION_METHODS[[params$agg_method]] %||% params$agg_method %||% "(not set)"),
    sprintf("Correlation method: %s", params$cor_method %||% "(not computed - unpaired data)"),
    sprintf("Multiple-testing adjustment: %s", params$padj_method %||% "(not set)"),
    sprintf("Sample matching: %s", params$sample_matching %||% "Not available"),
    sprintf("Gene annotation source: %s", params$gene_annotation_source %||% "org.Hs.eg.db (Bioconductor) - exact ID/alias lookup only, no fuzzy matching"),
    sprintf("Methylation platform annotation: %s", params$methylation_platform %||% "Not available"),
    if (!is.null(params$methylation_platform) && grepl("Illumina", params$methylation_platform))
      "Note: CpG probes annotated to more than one gene (Illumina UCSC_RefGene_Name) are assigned to their first-listed gene only; other co-annotated genes for those probes are not separately represented in this integration.",
    sprintf("Module version: %s", CX_MODULE_VERSION),
    sprintf("Run at: %s", params$run_at %||% "(not run yet)")
  )
}

## ---------------------------------------------------------------------------
## Gene identifier harmonization (spec section 7) - HGNC symbol / Entrez ID /
## Ensembl Gene ID all resolve to one canonical HGNC symbol where a valid,
## EXACT match exists. Deliberately no fuzzy string matching anywhere in this
## file - an id that isn't an exact symbol/Entrez/Ensembl match is looked up
## in the alias table (org.Hs.eg.db's ALIAS keytype) and, if that resolves to
## exactly one canonical symbol, is reported as "alias_resolved"; if it
## resolves to more than one, it is reported "ambiguous" (never guessed);
## anything left over is "unmatched". Uses the same already-installed
## org.Hs.eg.db/AnnotationDbi dependency the Pathways tab uses
## (cx_run_pathway_enrichment(), crossomics_integration_plots.R) - no new
## package dependency.
## ---------------------------------------------------------------------------

.cx_hgnc_cache <- new.env(parent = emptyenv())

## One-time build of a full SYMBOL -> ENTREZID/ENSEMBL table and a full
## ALIAS -> SYMBOL table from org.Hs.eg.db. Cached for the life of the R
## process (same pattern as .cx_anno_cache above) - building both tables
## takes a few seconds, so this must not be repeated per gene or per run.
cx_build_id_lookup <- function() {
  cached <- .cx_hgnc_cache[["lookup"]]
  if (!is.null(cached)) return(cached)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(NULL)
  }
  db <- org.Hs.eg.db::org.Hs.eg.db
  main <- tryCatch({
    syms <- AnnotationDbi::keys(db, keytype = "SYMBOL")
    m <- AnnotationDbi::select(db, keys = syms, keytype = "SYMBOL", columns = c("ENTREZID", "ENSEMBL"))
    ## SYMBOL -> {ENTREZID, ENSEMBL} is occasionally 1:many (a gene can have
    ## more than one Ensembl Gene ID recorded) - keep one deterministic
    ## representative row per symbol so downstream matching is 1:1.
    m[!duplicated(m$SYMBOL), , drop = FALSE]
  }, error = function(e) NULL)
  aliases <- tryCatch({
    al_keys <- AnnotationDbi::keys(db, keytype = "ALIAS")
    AnnotationDbi::select(db, keys = al_keys, keytype = "ALIAS", columns = "SYMBOL")
  }, error = function(e) NULL)
  if (is.null(main)) return(NULL)
  out <- list(main = main, aliases = aliases)
  .cx_hgnc_cache[["lookup"]] <- out
  out
}

## Given a vector of raw gene identifiers (any mix of HGNC symbols, Entrez
## IDs, or Ensembl Gene IDs - versioned or not), returns
## list(ok, df(input_id, canonical_symbol, entrez_id, ensembl_id,
## match_type), summary(exact, alias_resolved, unmatched, ambiguous)).
## match_type is one of exact_symbol/exact_entrez/exact_ensembl/
## alias_resolved/ambiguous/unmatched. When the annotation package isn't
## available, returns ok = FALSE with an explanatory error rather than
## silently falling back to fuzzy matching.
cx_harmonize_gene_ids <- function(genes) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(trimws(genes))]
  if (length(genes) == 0) return(list(ok = FALSE, error = "No gene identifiers to harmonize.", df = NULL, summary = NULL))
  lut <- cx_build_id_lookup()
  if (is.null(lut)) {
    return(list(ok = FALSE, error = "Gene identifier annotation (org.Hs.eg.db) is not available in this deployment - harmonization skipped; matching falls back to exact text on the Gene column as provided.", df = NULL, summary = NULL))
  }
  main <- lut$main
  clean_id <- sub("\\.[0-9]+$", "", trimws(genes))  ## strip Ensembl version suffix, e.g. ENSG00000141510.5

  out <- data.frame(input_id = genes, clean_id = clean_id, canonical_symbol = NA_character_,
                     entrez_id = NA_character_, ensembl_id = NA_character_, match_type = "unmatched",
                     stringsAsFactors = FALSE)

  .assign_hits <- function(out, idx, hit_row, match_type) {
    out$canonical_symbol[idx] <- main$SYMBOL[hit_row]
    out$entrez_id[idx] <- main$ENTREZID[hit_row]
    out$ensembl_id[idx] <- main$ENSEMBL[hit_row]
    out$match_type[idx] <- match_type
    out
  }

  ## 1. Exact symbol match (case-insensitive on the text only - not fuzzy,
  ## since capitalization variants of the identical symbol are still the
  ## identical gene, per spec section 7's own "handle capitalization" note).
  hit <- match(toupper(out$clean_id), toupper(main$SYMBOL))
  idx <- which(!is.na(hit))
  if (length(idx) > 0) out <- .assign_hits(out, idx, hit[idx], "exact_symbol")

  ## 2. Exact Entrez ID match (remaining unmatched only).
  rem <- which(out$match_type == "unmatched")
  hit <- match(out$clean_id[rem], main$ENTREZID)
  ok <- !is.na(hit)
  if (any(ok)) out <- .assign_hits(out, rem[ok], hit[ok], "exact_entrez")

  ## 3. Exact Ensembl Gene ID match (remaining unmatched only).
  rem <- which(out$match_type == "unmatched")
  hit <- match(out$clean_id[rem], main$ENSEMBL)
  ok <- !is.na(hit)
  if (any(ok)) out <- .assign_hits(out, rem[ok], hit[ok], "exact_ensembl")

  ## 4. Alias resolution (remaining unmatched only) - exactly one distinct
  ## canonical symbol -> alias_resolved; more than one -> ambiguous (every
  ## candidate symbol is listed, never guessed at); zero -> stays unmatched.
  rem <- which(out$match_type == "unmatched")
  if (!is.null(lut$aliases) && length(rem) > 0) {
    al <- lut$aliases[!is.na(lut$aliases$SYMBOL), , drop = FALSE]
    grp <- split(al$SYMBOL, toupper(al$ALIAS))
    rem_key <- toupper(out$clean_id[rem])
    for (k in seq_along(rem)) {
      cand <- unique(grp[[rem_key[k]]])
      if (is.null(cand)) next
      i <- rem[k]
      if (length(cand) == 1) {
        m <- match(cand, main$SYMBOL)
        out$canonical_symbol[i] <- cand
        out$entrez_id[i] <- if (!is.na(m)) main$ENTREZID[m] else NA_character_
        out$ensembl_id[i] <- if (!is.na(m)) main$ENSEMBL[m] else NA_character_
        out$match_type[i] <- "alias_resolved"
      } else {
        out$canonical_symbol[i] <- paste(sort(cand), collapse = ";")
        out$match_type[i] <- "ambiguous"
      }
    }
  }

  summary <- list(
    exact = sum(out$match_type %in% c("exact_symbol", "exact_entrez", "exact_ensembl")),
    alias_resolved = sum(out$match_type == "alias_resolved"),
    unmatched = sum(out$match_type == "unmatched"),
    ambiguous = sum(out$match_type == "ambiguous")
  )
  list(ok = TRUE, df = out[, c("input_id", "canonical_symbol", "entrez_id", "ensembl_id", "match_type")],
       summary = summary, error = NULL)
}

## Rewrites a vector of raw gene ids to their harmonized canonical symbol
## wherever cx_harmonize_gene_ids() found an unambiguous match (exact or
## alias-resolved); leaves ambiguous/unmatched entries as their original
## text unchanged, so the existing exact-text join still applies to them as
## a fallback rather than silently dropping them.
cx_apply_harmonization <- function(genes, harm_df) {
  if (is.null(harm_df)) return(genes)
  usable <- harm_df$match_type %in% c("exact_symbol", "exact_entrez", "exact_ensembl", "alias_resolved")
  lut <- stats::setNames(harm_df$canonical_symbol[usable], harm_df$input_id[usable])
  ifelse(genes %in% names(lut), unname(lut[genes]), genes)
}

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## Dataset validation panel (spec section 6) - a plain-language checklist
## computed directly from what's actually loaded; never inferred or invented.
## ---------------------------------------------------------------------------

cx_validate_dataset <- function(expr_df, meth_df, id_harmonization = NULL) {
  tx_ok <- !is.null(expr_df) && nrow(expr_df) > 0
  mx_ok <- !is.null(meth_df) && nrow(meth_df) > 0
  tx_checks <- list(
    list(ok = tx_ok, label = if (tx_ok) sprintf("Gene identifier detected (%s genes)", format(nrow(expr_df), big.mark = ",")) else "Gene identifier not detected"),
    list(ok = tx_ok && any(!is.na(expr_df$log2fc)), label = if (tx_ok && any(!is.na(expr_df$log2fc))) "log2 Fold Change detected" else "log2 Fold Change not detected"),
    list(ok = tx_ok && any(!is.na(expr_df$fdr)), label = if (tx_ok && any(!is.na(expr_df$fdr))) "Adjusted P value (FDR) detected" else "Adjusted P value not detected (raw P value only, or missing)")
  )
  mx_checks <- list(
    list(ok = mx_ok, label = if (mx_ok) sprintf("CpG ID detected (%s CpG records)", format(nrow(meth_df), big.mark = ",")) else "CpG ID not detected"),
    list(ok = mx_ok && any(!is.na(meth_df$dbeta)), label = if (mx_ok && any(!is.na(meth_df$dbeta))) "Methylation measurement (Δβ) detected" else "Methylation measurement not available"),
    list(ok = mx_ok && any(!is.na(meth_df$gene)), label = if (mx_ok && any(!is.na(meth_df$gene))) "Gene annotation detected" else "Gene annotation not available")
  )
  n_overlap <- if (tx_ok && mx_ok) length(intersect(expr_df$gene, unique(meth_df$gene))) else NA_integer_
  compat <- list(
    list(ok = tx_ok && mx_ok, label = if (tx_ok && mx_ok) "Gene identifiers can be matched" else "Cannot perform gene-level methylation integration - both a Transcriptomics and a Methylomics dataset are required"),
    list(ok = !is.na(n_overlap) && n_overlap > 0, label = if (!is.na(n_overlap)) sprintf("%s transcriptomic genes have methylation annotation", format(n_overlap, big.mark = ",")) else "Overlap not yet calculable")
  )
  if (!is.null(id_harmonization) && isTRUE(id_harmonization$ok)) {
    s <- id_harmonization$summary
    compat <- c(compat, list(list(ok = TRUE, label = sprintf(
      "Gene ID harmonization: %s exact, %s alias-resolved, %s unmatched, %s ambiguous",
      s$exact, s$alias_resolved, s$unmatched, s$ambiguous))))
  }
  list(transcriptomics = tx_checks, methylomics = mx_checks, compatibility = compat, ready = tx_ok && mx_ok)
}
