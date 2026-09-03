## R/multiomics/07_Pathways/multiomics_pathway_helpers.R
## Data-adaptive engine for the "Pathways" sub-module (mod_multi_pathway.R) -
## live GO/KEGG/Reactome/WikiPathways ORA + GSEA, on either the app's

MP_REACTOME_AVAILABLE <- requireNamespace("ReactomePA", quietly = TRUE) && requireNamespace("reactome.db", quietly = TRUE)
MP_FGSEA_AVAILABLE    <- requireNamespace("fgsea", quietly = TRUE)
MP_MSIGDBR_AVAILABLE  <- requireNamespace("msigdbr", quietly = TRUE)
MP_PATHVIEW_AVAILABLE <- requireNamespace("pathview", quietly = TRUE)

MP_DATABASES <- list(
  GO_BP        = list(label = "GO - Biological Process", category = "GO",           topology = FALSE, available = TRUE),
  GO_MF        = list(label = "GO - Molecular Function", category = "GO",           topology = FALSE, available = TRUE),
  GO_CC        = list(label = "GO - Cellular Component", category = "GO",           topology = FALSE, available = TRUE),
  KEGG         = list(label = "KEGG",                     category = "KEGG",        topology = FALSE, available = TRUE),
  Reactome     = list(label = "Reactome",                 category = "Reactome",    topology = TRUE,  available = MP_REACTOME_AVAILABLE),
  WikiPathways = list(label = "WikiPathways",              category = "GMT",         topology = FALSE, available = MP_MSIGDBR_AVAILABLE),
  Hallmark     = list(label = "MSigDB Hallmark", category = "GMT", topology = FALSE, available = MP_MSIGDBR_AVAILABLE)
)

MP_FIELD_PATTERNS <- list(
  id_symbol = c("^gene[_ .]?symbol$", "^hgnc[_ .]?symbol$", "^symbol$", "^gene[_ .]?name$", "^gene$", "^genes$"),
  id_entrez = c("^entrez[_ .]?(id|gene)?$", "^entrezgene$", "^ncbi[_ .]?gene[_ .]?id$"),
  id_ensembl = c("^ensembl[_ .]?(gene)?[_ .]?id$", "^ensg$", "^ensembl$"),
  id_uniprot = c("^uniprot[_ .]?(id|accession)?$", "^swissprot$", "^protein[_ .]?id$"),
  id_cpg = c("^cpg[_ .]?id$", "^cpg$", "^probe[_ .]?id$", "^probeid$", "^probe$", "^illumina[_ .]?id$"),
  effect = c("^log2[_ .]?fc$", "^log2fc$", "^logfc$", "^log[_ .]?fold[_ .]?change$", "^fold[_ .]?change$",
             "^lfc$", "^statistic$", "^t$", "^z$", "^zscore$", "^nes$", "^effect[_ .]?size$", "^delta[_ .]?beta$", "^dbeta$"),
  pvalue = c("^p[_ .]?value$", "^pval$", "^p[_ .]?val$", "^pvalue$", "^p\\.value$"),
  fdr = c("^fdr$", "^adj\\.?p\\.?val$", "^adjusted[_ .]?p[_ .]?value$", "^padj$", "^p\\.adjust$", "^q[_ .]?value$", "^qval$"),
  direction = c("^direction$", "^regulation$", "^expr[_ .]?direction$", "^meth[_ .]?direction$", "^change[_ .]?direction$"),
  omics = c("^omics$", "^omics[_ .]?type$", "^layer$", "^data[_ .]?type$", "^assay$"),
  sex = c("^sex$", "^gender$")
)

MP_PRELOADED_COHORTS <- c(
  "Drug x sex (Etanercept panel)" = "Gene <-> CpG mapping - drug x sex (Etanercept panel)",
  "Response (drug-pooled)" = "Gene <-> CpG mapping - response (drug-pooled)"
)

mp_build_preloaded_input <- function(multi_results, multi_dataset, expr_layer, meth_layer, cohort_label, sex,
                                      custom_genes = character(0), custom_cpgs = character(0)) {
  pool <- mcc_candidate_pool(multi_results, multi_dataset, expr_layer, meth_layer, custom_genes, custom_cpgs)
  if (!isTRUE(pool$ok)) return(list(ok = FALSE, df = NULL, error = pool$note %||% "No candidate biomarkers available - run Biomarker Discovery and/or Patient Stratification first, or supply custom genes/CpGs below."))
  df <- pool$df
  df$evidence_source <- ifelse(df$diablo & df$snf, "DIABLO + SNF", ifelse(df$diablo, "DIABLO", ifelse(df$snf, "SNF", "Custom")))

  reg <- multi_read_registry_table(cohort_label)
  df$gene_symbol <- NA_character_; df$expr_logFC <- NA_real_; df$expr_p <- NA_real_; df$expr_direction <- NA_character_
  df$cpg <- NA_character_; df$meth_delta_M <- NA_real_; df$meth_p <- NA_real_; df$meth_direction <- NA_character_
  df$region <- NA_character_; df$island <- NA_character_; df$biological_pattern <- NA_character_
  source_detail <- "No matching precomputed effect-size table for this cohort/sex - candidate features are still analyzable, but without per-feature effect sizes."

  if (isTRUE(reg$ok)) {
    conc <- reg$df
    if (length(sex) == 1 && sex %in% c("female", "male") && "sex" %in% colnames(conc)) conc <- conc[tolower(conc$sex) == sex, , drop = FALSE]
    gene_idx_by_symbol <- match(toupper(df$feature), toupper(conc$SYMBOL))
    gene_idx_by_ensembl <- if ("ENSEMBL" %in% colnames(conc)) match(df$feature, conc$ENSEMBL) else rep(NA_integer_, nrow(df))
    gene_hit <- !is.na(gene_idx_by_symbol) | !is.na(gene_idx_by_ensembl)
    gi <- ifelse(!is.na(gene_idx_by_symbol), gene_idx_by_symbol, gene_idx_by_ensembl)
    if (any(gene_hit)) {
      df$gene_symbol[gene_hit] <- conc$SYMBOL[gi[gene_hit]]
      df$expr_logFC[gene_hit] <- conc$expr_logFC[gi[gene_hit]]
      df$expr_p[gene_hit] <- conc$expr_p[gi[gene_hit]]
      df$expr_direction[gene_hit] <- conc$expr_direction[gi[gene_hit]]
      df$cpg[gene_hit] <- conc$CpG[gi[gene_hit]]
      df$meth_delta_M[gene_hit] <- conc$delta_M[gi[gene_hit]]
      df$meth_p[gene_hit] <- conc$meth_p[gi[gene_hit]]
      df$meth_direction[gene_hit] <- conc$meth_direction[gi[gene_hit]]
      df$region[gene_hit] <- conc$region[gi[gene_hit]]
      df$island[gene_hit] <- conc$island[gi[gene_hit]]
      df$biological_pattern[gene_hit] <- conc$biological_pattern[gi[gene_hit]]
    }
    cpg_idx <- match(df$feature, conc$CpG)
    cpg_hit <- !is.na(cpg_idx) & is.na(df$cpg)
    if (any(cpg_hit)) {
      df$cpg[cpg_hit] <- conc$CpG[cpg_idx[cpg_hit]]
      df$gene_symbol[cpg_hit] <- conc$SYMBOL[cpg_idx[cpg_hit]]
      df$expr_logFC[cpg_hit] <- conc$expr_logFC[cpg_idx[cpg_hit]]
      df$expr_p[cpg_hit] <- conc$expr_p[cpg_idx[cpg_hit]]
      df$expr_direction[cpg_hit] <- conc$expr_direction[cpg_idx[cpg_hit]]
      df$meth_delta_M[cpg_hit] <- conc$delta_M[cpg_idx[cpg_hit]]
      df$meth_p[cpg_hit] <- conc$meth_p[cpg_idx[cpg_hit]]
      df$meth_direction[cpg_hit] <- conc$meth_direction[cpg_idx[cpg_hit]]
      df$region[cpg_hit] <- conc$region[cpg_idx[cpg_hit]]
      df$island[cpg_hit] <- conc$island[cpg_idx[cpg_hit]]
      df$biological_pattern[cpg_hit] <- conc$biological_pattern[cpg_idx[cpg_hit]]
    }
    source_detail <- sprintf("%s candidate feature(s); %s matched to a precomputed effect size in \"%s\".",
                              nrow(df), sum(gene_hit | cpg_hit), names(MP_PRELOADED_COHORTS)[MP_PRELOADED_COHORTS == cohort_label])
  }
  list(ok = TRUE, df = df, source = "preloaded", source_detail = source_detail)
}

MP_UPLOAD_NA_STRINGS <- c("NA", "", "NaN", "null", "NULL", "#N/A")

mp_read_upload_table <- function(datapath, filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("csv", "tsv", "txt")) {
    df <- tryCatch(as.data.frame(data.table::fread(datapath, showProgress = FALSE, na.strings = MP_UPLOAD_NA_STRINGS)), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0 || ncol(df) < 1) return(list(ok = FALSE, df = NULL, error = "Could not parse this file as a table (CSV/TSV/TXT with at least one column)."))
    return(list(ok = TRUE, df = df, error = NULL))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) return(list(ok = FALSE, df = NULL, error = "The openxlsx package is not installed in this deployment - export your file as CSV/TSV instead."))
    df <- tryCatch(openxlsx::read.xlsx(datapath, sheet = 1), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0 || ncol(df) < 1) return(list(ok = FALSE, df = NULL, error = "Could not parse this XLSX file's first sheet as a table."))
    return(list(ok = TRUE, df = df, error = NULL))
  }
  list(ok = FALSE, df = NULL, error = sprintf("Unsupported file type \".%s\" - upload CSV, TSV, TXT, or XLSX.", ext))
}

mp_match_column <- function(cols, patterns, exclude = character(0)) {
  candidates <- setdiff(cols, exclude)
  candidates_lower <- tolower(trimws(candidates))
  for (p in patterns) {
    hit <- grep(p, candidates_lower, perl = TRUE)
    if (length(hit) >= 1) return(candidates[hit[1]])
  }
  NA_character_
}

mp_detect_upload <- function(df) {
  cols <- colnames(df)
  claimed <- character(0)
  role <- function(f) { hit <- mp_match_column(cols, MP_FIELD_PATTERNS[[f]], exclude = claimed); if (!is.na(hit)) claimed <<- c(claimed, hit); hit }
  id_symbol <- role("id_symbol"); id_entrez <- role("id_entrez"); id_ensembl <- role("id_ensembl")
  id_uniprot <- role("id_uniprot"); id_cpg <- role("id_cpg")
  effect_col <- role("effect"); pvalue_col <- role("pvalue"); fdr_col <- role("fdr")
  direction_col <- role("direction"); omics_col <- role("omics"); sex_col <- role("sex")

  id_col <- if (!is.na(id_cpg)) id_cpg else if (!is.na(id_symbol)) id_symbol else if (!is.na(id_ensembl)) id_ensembl else if (!is.na(id_entrez)) id_entrez else if (!is.na(id_uniprot)) id_uniprot else NA_character_
  id_type <- if (!is.na(id_cpg)) "Illumina CpG probe ID" else if (!is.na(id_ensembl)) "Ensembl Gene ID" else if (!is.na(id_entrez)) "Entrez ID" else if (!is.na(id_uniprot)) "UniProt ID" else if (!is.na(id_symbol)) "Gene symbol" else "unknown"
  if (is.na(id_col) && ncol(df) > 0) {
    id_col <- cols[1]
    id_type <- mcc_detect_id_type(df[[id_col]])
  }

  warnings <- character(0)
  if (is.na(effect_col)) warnings <- c(warnings, "No effect-size/log2FC/statistic column detected - GSEA (which requires a ranked list) will not be available.")
  if (is.na(pvalue_col) && is.na(fdr_col)) warnings <- c(warnings, "No P-value or adjusted-P column detected - ORA will use the full identifier list as-is, with no significance-based selection possible.")
  if (identical(id_type, "unknown")) warnings <- c(warnings, "Could not determine the identifier type of the detected ID column - mapping/enrichment will likely fail until a column is manually selected.")

  species_guess <- if (!is.na(id_ensembl)) {
    vals <- as.character(df[[id_ensembl]])
    if (mean(grepl("^ENSG", vals)) > 0.5) "Homo sapiens" else if (mean(grepl("^ENSMUSG", vals)) > 0.5) "Mus musculus" else if (mean(grepl("^ENSRNOG", vals)) > 0.5) "Rattus norvegicus" else "unknown"
  } else "Homo sapiens (assumed - no Ensembl ID column to confirm against)"

  omics_guess <- if (identical(id_type, "Illumina CpG probe ID")) "Methylomics" else if (id_type %in% c("Gene symbol", "Entrez ID", "Ensembl Gene ID", "UniProt ID")) "Transcriptomics" else "unknown"

  ranked <- !is.na(effect_col) && is.numeric(df[[effect_col]]) && any(df[[effect_col]] < 0, na.rm = TRUE) && any(df[[effect_col]] > 0, na.rm = TRUE)

  list(
    detected = list(id_col = id_col, id_type = id_type, id_symbol_col = id_symbol, id_entrez_col = id_entrez,
                     id_ensembl_col = id_ensembl, id_uniprot_col = id_uniprot, id_cpg_col = id_cpg,
                     effect_col = effect_col, pvalue_col = pvalue_col, fdr_col = fdr_col,
                     direction_col = direction_col, omics_col = omics_col, sex_col = sex_col,
                     ranked = ranked, species_guess = species_guess, omics_guess = omics_guess,
                     n_rows = nrow(df), n_cols = ncol(df)),
    warnings = warnings
  )
}

mp_confirm_upload_mapping <- function(df, mapping) {
  none <- function(x) is.null(x) || is.na(x) || identical(x, "(none)")
  if (none(mapping$id_col) || !mapping$id_col %in% colnames(df)) return(list(ok = FALSE, df = NULL, error = "Select an identifier column before confirming."))
  get_col <- function(name) if (!none(mapping[[name]]) && mapping[[name]] %in% colnames(df)) df[[mapping[[name]]]] else NA

  n <- nrow(df)
  out <- data.frame(
    feature = as.character(df[[mapping$id_col]]),
    id_type = mapping$id_type %||% mcc_detect_id_type(df[[mapping$id_col]]),
    omics = if (!none(mapping$omics_col)) as.character(get_col("omics_col")) else rep(mapping$omics_guess %||% NA_character_, n),
    diablo = FALSE, snf = FALSE, joint = FALSE, custom = TRUE, evidence_source = "Uploaded",
    selection_frequency = NA_real_, stability_category = NA_character_,
    gene_symbol = if (mapping$id_type %in% c("Gene symbol")) as.character(df[[mapping$id_col]]) else NA_character_,
    expr_logFC = if (!none(mapping$effect_col)) suppressWarnings(as.numeric(get_col("effect_col"))) else NA_real_,
    expr_p = if (!none(mapping$pvalue_col)) suppressWarnings(as.numeric(get_col("pvalue_col"))) else NA_real_,
    expr_direction = if (!none(mapping$direction_col)) as.character(get_col("direction_col")) else NA_character_,
    cpg = if (identical(mapping$id_type, "Illumina CpG probe ID")) as.character(df[[mapping$id_col]]) else NA_character_,
    meth_delta_M = NA_real_, meth_p = NA_real_, meth_direction = NA_character_,
    region = NA_character_, island = NA_character_, biological_pattern = NA_character_,
    sex = if (!none(mapping$sex_col)) as.character(get_col("sex_col")) else NA_character_,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$feature) & nzchar(trimws(out$feature)), , drop = FALSE]
  out <- out[!duplicated(out$feature), , drop = FALSE]
  if (nrow(out) == 0) return(list(ok = FALSE, df = NULL, error = "No usable rows after removing empty/duplicate identifiers."))
  attr(out, "fdr_col_present") <- !none(mapping$fdr_col)
  if (!none(mapping$fdr_col)) out$expr_fdr <- suppressWarnings(as.numeric(get_col("fdr_col")))
  list(ok = TRUE, df = out, source = "uploaded", source_detail = sprintf("%d of %d uploaded rows usable after de-duplication.", nrow(out), n))
}

mp_map_uniprot_to_symbol <- function(ids) {
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) == 0 || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
  tryCatch(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ids, keytype = "UNIPROT", columns = c("SYMBOL", "ENTREZID")), error = function(e) NULL)
}

mp_map_candidate_cpgs <- function(cpg_ids, array_type = "450K") {
  cpg_ids <- unique(cpg_ids[!is.na(cpg_ids) & nzchar(cpg_ids)])
  if (length(cpg_ids) == 0) return(list(ok = FALSE, df = NULL, error = "No CpG identifiers to map.", n_mapped = 0L, n_unique_genes = 0L, unmapped_cpgs = character(0)))
  ar <- cx_get_region_annotation(array_type)
  if (!isTRUE(ar$ok)) return(list(ok = FALSE, df = NULL, error = ar$reason, n_mapped = 0L, n_unique_genes = 0L, unmapped_cpgs = cpg_ids))
  anno <- ar$anno
  present <- intersect(cpg_ids, rownames(anno))
  sub <- anno[present, , drop = FALSE]
  has_gene <- !is.na(sub$gene) & nzchar(sub$gene)
  mapped <- sub[has_gene, , drop = FALSE]
  unmapped <- union(setdiff(cpg_ids, rownames(anno)), rownames(sub)[!has_gene])
  df <- data.frame(cpg = rownames(mapped), gene = mapped$gene, region = mapped$region_raw, island = mapped$island_context, stringsAsFactors = FALSE)
  list(ok = nrow(df) > 0, df = df, error = if (nrow(df) == 0) "None of the supplied CpGs are annotated to a gene in the selected array's annotation." else NULL,
       n_mapped = nrow(df), n_unique_genes = length(unique(df$gene)), unmapped_cpgs = unmapped)
}

mp_harmonize_identifiers <- function(input_df, array_type = "450K") {
  df <- input_df
  row_id_type <- df$id_type %||% vapply(df$feature, mcc_detect_id_type, character(1))
  is_cpg <- row_id_type == "Illumina CpG probe ID"
  gene_features <- unique(df$feature[!is_cpg])
  gene_features <- gene_features[!is.na(gene_features) & nzchar(gene_features)]

  harm <- if (length(gene_features) > 0) cx_harmonize_gene_ids(gene_features) else list(ok = FALSE, df = NULL, error = "No gene-like identifiers to harmonize.")

  cpg_map <- NULL
  if (any(is_cpg)) cpg_map <- mp_map_candidate_cpgs(df$feature[is_cpg], array_type)

  out <- data.frame(feature = df$feature, id_type = df$id_type %||% NA_character_, canonical_symbol = NA_character_,
                     entrez_id = NA_character_, ensembl_id = NA_character_, match_type = "unmatched", mapped = FALSE, stringsAsFactors = FALSE)
  if (isTRUE(harm$ok)) {
    idx <- match(toupper(out$feature), toupper(harm$df$input_id))
    hit <- !is.na(idx) & !is_cpg
    out$canonical_symbol[hit] <- harm$df$canonical_symbol[idx[hit]]
    out$entrez_id[hit] <- harm$df$entrez_id[idx[hit]]
    out$ensembl_id[hit] <- harm$df$ensembl_id[idx[hit]]
    out$match_type[hit] <- harm$df$match_type[idx[hit]]
    out$mapped[hit] <- harm$df$match_type[idx[hit]] %in% c("exact_symbol", "exact_entrez", "exact_ensembl", "alias_resolved")
  }
  if (!is.null(cpg_map) && isTRUE(cpg_map$ok)) {
    ci <- match(out$feature, cpg_map$df$cpg)
    chit <- !is.na(ci) & is_cpg
    out$canonical_symbol[chit] <- cpg_map$df$gene[ci[chit]]
    out$match_type[chit] <- "cpg_mapped"
    out$mapped[chit] <- TRUE
    gh <- cx_harmonize_gene_ids(unique(stats::na.omit(out$canonical_symbol[chit])))
    if (isTRUE(gh$ok)) {
      gi <- match(toupper(out$canonical_symbol[chit]), toupper(gh$df$input_id))
      out$entrez_id[chit][!is.na(gi)] <- gh$df$entrez_id[gi[!is.na(gi)]]
      out$ensembl_id[chit][!is.na(gi)] <- gh$df$ensembl_id[gi[!is.na(gi)]]
    }
  }
  list(ok = TRUE, df = out, cpg_map = cpg_map)
}

mp_mapping_summary <- function(mapped_df) {
  n <- nrow(mapped_df)
  n_mapped <- sum(mapped_df$mapped, na.rm = TRUE)
  list(n_input = n, n_mapped = n_mapped, n_unmapped = n - n_mapped,
       mapping_rate = if (n > 0) round(100 * n_mapped / n, 1) else 0,
       unmapped_ids = mapped_df$feature[!mapped_df$mapped])
}

mp_resolve_universe <- function(background_choice, multi_dataset, expr_layer, meth_layer, uploaded_universe_ids = NULL) {
  ids <- switch(background_choice,
    "auto_experimental" = {
      layers <- multi_dataset$layers %||% list()
      feats <- character(0)
      if (!is.null(expr_layer) && expr_layer %in% names(layers)) feats <- c(feats, colnames(layers[[expr_layer]]))
      if (!is.null(meth_layer) && meth_layer %in% names(layers)) {
        cpg_map <- mp_map_candidate_cpgs(colnames(layers[[meth_layer]]))
        if (isTRUE(cpg_map$ok)) feats <- c(feats, unique(cpg_map$df$gene))
      }
      feats
    },
    "uploaded_background" = uploaded_universe_ids,
    "preloaded_universe" = {
      reg <- multi_read_registry_table("Gene <-> CpG mapping - drug x sex (Etanercept panel)")
      if (isTRUE(reg$ok)) unique(reg$df$SYMBOL) else character(0)
    },
    "entire_database" = NULL,
    character(0)
  )
  if (identical(background_choice, "entire_database")) {
    return(list(ok = TRUE, universe_entrez = NULL, universe_label = "Entire selected database - no experimental universe supplied.", n = NA_integer_))
  }
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (identical(background_choice, "auto_experimental") && length(ids) == 0) {
    return(list(ok = TRUE, universe_entrez = NULL,
                universe_label = "No active dataset layers to build an experimental universe from - falling back to the entire selected database.",
                n = NA_integer_))
  }
  if (length(ids) == 0) return(list(ok = FALSE, universe_entrez = NULL, universe_label = NULL, n = 0L, error = "No background/universe identifiers available for the selected option."))
  harm <- cx_harmonize_gene_ids(ids)
  entrez <- if (isTRUE(harm$ok)) unique(stats::na.omit(harm$df$entrez_id)) else character(0)
  if (length(entrez) == 0) return(list(ok = FALSE, universe_entrez = NULL, universe_label = NULL, n = 0L, error = "None of the background identifiers could be mapped to Entrez IDs."))
  label <- switch(background_choice, auto_experimental = "Measured features in the active dataset's layer(s)",
                   uploaded_background = "Uploaded background file",
                   preloaded_universe = "Preloaded cohort's own candidate gene list (already filtered to significant hits, not a genome/panel-wide background)",
                   background_choice)
  list(ok = TRUE, universe_entrez = entrez, universe_label = sprintf("%s (%s genes)", label, format(length(entrez), big.mark = ",")), n = length(entrez))
}

.mp_msigdbr_cache <- new.env(parent = emptyenv())

mp_get_wikipathways_termgene <- function() {
  cached <- .mp_msigdbr_cache[["wikipathways"]]
  if (!is.null(cached)) return(cached)
  if (!MP_MSIGDBR_AVAILABLE) return(NULL)
  wp <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS"), error = function(e) NULL)
  if (is.null(wp) || nrow(wp) == 0) return(NULL)
  term2gene <- unique(wp[, c("gs_exact_source", "ncbi_gene")])
  term2name <- unique(wp[, c("gs_exact_source", "gs_name")])
  out <- list(TERM2GENE = term2gene, TERM2NAME = term2name, db_version = unique(wp$db_version)[1])
  .mp_msigdbr_cache[["wikipathways"]] <- out
  out
}

mp_get_hallmark_termgene <- function() {
  cached <- .mp_msigdbr_cache[["hallmark"]]
  if (!is.null(cached)) return(cached)
  if (!MP_MSIGDBR_AVAILABLE) return(NULL)
  hm <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) NULL)
  if (is.null(hm) || nrow(hm) == 0) return(NULL)
  term2gene <- unique(hm[, c("gs_id", "ncbi_gene")])
  term2name <- unique(hm[, c("gs_id", "gs_name")])
  out <- list(TERM2GENE = term2gene, TERM2NAME = term2name, db_version = unique(hm$db_version)[1])
  .mp_msigdbr_cache[["hallmark"]] <- out
  out
}

mp_normalize_enrich_result <- function(res_df, source_label, method = "ORA") {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  ratio_to_num <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) if (length(p) == 2 && as.numeric(p[2]) > 0) as.numeric(p[1]) / as.numeric(p[2]) else NA_real_, numeric(1))
  }
  data.frame(
    source = source_label, ID = res_df$ID, Description = res_df$Description,
    GeneRatio = res_df$GeneRatio %||% NA_character_, gene_ratio_numeric = if ("GeneRatio" %in% colnames(res_df)) ratio_to_num(res_df$GeneRatio) else NA_real_,
    BgRatio = res_df$BgRatio %||% NA_character_, Count = res_df$Count %||% NA_integer_,
    pvalue = res_df$pvalue, p.adjust = res_df$p.adjust, qvalue = res_df$qvalue %||% NA_real_,
    geneID = res_df$geneID %||% res_df$core_enrichment %||% NA_character_, method = method,
    stringsAsFactors = FALSE
  )
}

mp_finalize_enrich <- function(df, label) {
  if (is.null(df) || nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = sprintf("%s returned no terms - check that enough identifiers mapped to Entrez IDs (see the Mapping Results table).", label)))
  list(ok = TRUE, df = df, error = NULL)
}

mp_run_ora_go <- function(genes_entrez, universe_entrez, ont, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  eg <- tryCatch(clusterProfiler::enrichGO(gene = genes_entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                            ont = ont, universe = universe_entrez, pAdjustMethod = "BH",
                                            pvalueCutoff = pvalueCutoff, qvalueCutoff = 1,
                                            minGSSize = minGSSize, maxGSSize = maxGSSize, readable = TRUE), error = function(e) e)
  if (inherits(eg, "error")) return(list(ok = FALSE, df = NULL, error = paste("GO enrichment failed:", conditionMessage(eg))))
  mp_finalize_enrich(mp_normalize_enrich_result(as.data.frame(eg), paste0("GO_", ont)), paste0("GO ", ont))
}

mp_run_ora_kegg <- function(genes_entrez, universe_entrez, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  ek <- tryCatch(clusterProfiler::enrichKEGG(gene = genes_entrez, organism = "hsa", universe = universe_entrez,
                                              pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = 1,
                                              minGSSize = minGSSize, maxGSSize = maxGSSize), error = function(e) e)
  if (inherits(ek, "error")) return(list(ok = FALSE, df = NULL, error = paste("KEGG enrichment failed (needs internet access to the KEGG REST API):", conditionMessage(ek))))
  mp_finalize_enrich(mp_normalize_enrich_result(as.data.frame(ek), "KEGG"), "KEGG")
}

mp_run_ora_reactome <- function(genes_entrez, universe_entrez, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  if (!MP_REACTOME_AVAILABLE) return(list(ok = FALSE, df = NULL, error = "ReactomePA/reactome.db are not installed in this deployment - Reactome is unavailable."))
  ep <- tryCatch(ReactomePA::enrichPathway(gene = genes_entrez, universe = universe_entrez, organism = "human",
                                            pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = 1,
                                            minGSSize = minGSSize, maxGSSize = maxGSSize, readable = TRUE), error = function(e) e)
  if (inherits(ep, "error")) return(list(ok = FALSE, df = NULL, error = paste("Reactome enrichment failed:", conditionMessage(ep))))
  mp_finalize_enrich(mp_normalize_enrich_result(as.data.frame(ep), "Reactome"), "Reactome")
}

mp_run_ora_wikipathways <- function(genes_entrez, universe_entrez, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  t2g <- mp_get_wikipathways_termgene()
  if (is.null(t2g)) return(list(ok = FALSE, df = NULL, error = "WikiPathways gene sets (msigdbr) are not available in this deployment."))
  ew <- tryCatch(clusterProfiler::enricher(gene = as.character(genes_entrez), universe = if (!is.null(universe_entrez)) as.character(universe_entrez) else NULL,
                                            TERM2GENE = t2g$TERM2GENE, TERM2NAME = t2g$TERM2NAME, pAdjustMethod = "BH",
                                            pvalueCutoff = pvalueCutoff, qvalueCutoff = 1, minGSSize = minGSSize, maxGSSize = maxGSSize), error = function(e) e)
  if (inherits(ew, "error")) return(list(ok = FALSE, df = NULL, error = paste("WikiPathways enrichment failed:", conditionMessage(ew))))
  mp_finalize_enrich(mp_normalize_enrich_result(as.data.frame(ew), "WikiPathways"), "WikiPathways")
}

mp_run_ora_hallmark <- function(genes_entrez, universe_entrez, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  t2g <- mp_get_hallmark_termgene()
  if (is.null(t2g)) return(list(ok = FALSE, df = NULL, error = "MSigDB Hallmark gene sets (msigdbr) are not available in this deployment."))
  eh <- tryCatch(clusterProfiler::enricher(gene = as.character(genes_entrez), universe = if (!is.null(universe_entrez)) as.character(universe_entrez) else NULL,
                                            TERM2GENE = t2g$TERM2GENE, TERM2NAME = t2g$TERM2NAME, pAdjustMethod = "BH",
                                            pvalueCutoff = pvalueCutoff, qvalueCutoff = 1, minGSSize = minGSSize, maxGSSize = maxGSSize), error = function(e) e)
  if (inherits(eh, "error")) return(list(ok = FALSE, df = NULL, error = paste("Hallmark enrichment failed:", conditionMessage(eh))))
  mp_finalize_enrich(mp_normalize_enrich_result(as.data.frame(eh), "Hallmark"), "MSigDB Hallmark")
}

mp_run_ora <- function(database, genes_entrez, universe_entrez, params) {
  switch(database,
    GO_BP = mp_run_ora_go(genes_entrez, universe_entrez, "BP", params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    GO_MF = mp_run_ora_go(genes_entrez, universe_entrez, "MF", params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    GO_CC = mp_run_ora_go(genes_entrez, universe_entrez, "CC", params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    KEGG = mp_run_ora_kegg(genes_entrez, universe_entrez, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    Reactome = mp_run_ora_reactome(genes_entrez, universe_entrez, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    WikiPathways = mp_run_ora_wikipathways(genes_entrez, universe_entrez, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    Hallmark = mp_run_ora_hallmark(genes_entrez, universe_entrez, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    list(ok = FALSE, df = NULL, error = sprintf("Unknown database '%s'.", database))
  )
}

mp_build_ranked_vector <- function(df, entrez_ids, ranking_method = c("log2fc", "statistic", "signed_neglog10p", "custom"), custom_col = NULL, custom_values = NULL) {
  ranking_method <- match.arg(ranking_method)
  vals <- if (identical(ranking_method, "log2fc")) df$expr_logFC
          else if (identical(ranking_method, "signed_neglog10p")) sign(df$expr_logFC %||% 1) * -log10(pmax(df$expr_p, .Machine$double.eps))
          else if (identical(ranking_method, "custom")) custom_values
          else df$expr_logFC
  keep <- !is.na(entrez_ids) & !is.na(vals) & is.finite(vals)
  if (sum(keep) < 10) return(list(ok = FALSE, vec = NULL, error = "Fewer than 10 features have both a mapped Entrez ID and a finite ranking value - not enough for a stable GSEA run."))
  d <- data.frame(entrez = entrez_ids[keep], value = vals[keep])
  d <- d[order(-abs(d$value)), , drop = FALSE]
  d <- d[!duplicated(d$entrez), , drop = FALSE]
  vec <- stats::setNames(d$value, d$entrez)
  list(ok = TRUE, vec = sort(vec, decreasing = TRUE), error = NULL)
}

mp_normalize_gsea_result <- function(res_df, source_label) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  data.frame(
    source = source_label, ID = res_df$ID, Description = res_df$Description,
    GeneRatio = NA_character_, gene_ratio_numeric = NA_real_, BgRatio = NA_character_, Count = res_df$setSize %||% NA_integer_,
    pvalue = res_df$pvalue, p.adjust = res_df$p.adjust, qvalue = res_df$qvalues %||% NA_real_,
    geneID = res_df$core_enrichment %||% NA_character_, method = "GSEA",
    NES = res_df$NES %||% NA_real_, ES = res_df$enrichmentScore %||% NA_real_,
    leading_edge = res_df$leading_edge %||% NA_character_, stringsAsFactors = FALSE
  )
}

mp_run_gsea_go <- function(ranked_vec, ont, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  gg <- tryCatch(clusterProfiler::gseGO(geneList = ranked_vec, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         ont = ont, pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff,
                                         minGSSize = minGSSize, maxGSSize = maxGSSize, eps = 0), error = function(e) e)
  if (inherits(gg, "error")) return(list(ok = FALSE, df = NULL, error = paste("GO GSEA failed:", conditionMessage(gg))))
  mp_finalize_enrich(mp_normalize_gsea_result(as.data.frame(gg), paste0("GO_", ont)), paste0("GO ", ont, " GSEA"))
}

mp_run_gsea_kegg <- function(ranked_vec, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  gk <- tryCatch(clusterProfiler::gseKEGG(geneList = ranked_vec, organism = "hsa", pAdjustMethod = "BH",
                                           pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize, eps = 0), error = function(e) e)
  if (inherits(gk, "error")) return(list(ok = FALSE, df = NULL, error = paste("KEGG GSEA failed (needs internet access):", conditionMessage(gk))))
  mp_finalize_enrich(mp_normalize_gsea_result(as.data.frame(gk), "KEGG"), "KEGG GSEA")
}

mp_run_gsea_reactome <- function(ranked_vec, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  if (!MP_REACTOME_AVAILABLE) return(list(ok = FALSE, df = NULL, error = "ReactomePA/reactome.db are not installed in this deployment - Reactome GSEA is unavailable."))
  gp <- tryCatch(ReactomePA::gsePathway(geneList = ranked_vec, organism = "human", pAdjustMethod = "BH",
                                         pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize, eps = 0), error = function(e) e)
  if (inherits(gp, "error")) return(list(ok = FALSE, df = NULL, error = paste("Reactome GSEA failed:", conditionMessage(gp))))
  mp_finalize_enrich(mp_normalize_gsea_result(as.data.frame(gp), "Reactome"), "Reactome GSEA")
}

mp_run_gsea_wikipathways <- function(ranked_vec, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  t2g <- mp_get_wikipathways_termgene()
  if (is.null(t2g)) return(list(ok = FALSE, df = NULL, error = "WikiPathways gene sets (msigdbr) are not available in this deployment."))
  gw <- tryCatch(clusterProfiler::GSEA(geneList = ranked_vec, TERM2GENE = t2g$TERM2GENE, TERM2NAME = t2g$TERM2NAME,
                                        pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize, eps = 0), error = function(e) e)
  if (inherits(gw, "error")) return(list(ok = FALSE, df = NULL, error = paste("WikiPathways GSEA failed:", conditionMessage(gw))))
  mp_finalize_enrich(mp_normalize_gsea_result(as.data.frame(gw), "WikiPathways"), "WikiPathways GSEA")
}

mp_run_gsea_hallmark <- function(ranked_vec, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500) {
  t2g <- mp_get_hallmark_termgene()
  if (is.null(t2g)) return(list(ok = FALSE, df = NULL, error = "MSigDB Hallmark gene sets (msigdbr) are not available in this deployment."))
  gh <- tryCatch(clusterProfiler::GSEA(geneList = ranked_vec, TERM2GENE = t2g$TERM2GENE, TERM2NAME = t2g$TERM2NAME,
                                        pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize, eps = 0), error = function(e) e)
  if (inherits(gh, "error")) return(list(ok = FALSE, df = NULL, error = paste("Hallmark GSEA failed:", conditionMessage(gh))))
  mp_finalize_enrich(mp_normalize_gsea_result(as.data.frame(gh), "Hallmark"), "MSigDB Hallmark GSEA")
}

mp_run_gsea <- function(database, ranked_vec, params) {
  switch(database,
    GO_BP = mp_run_gsea_go(ranked_vec, "BP", params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    GO_MF = mp_run_gsea_go(ranked_vec, "MF", params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    GO_CC = mp_run_gsea_go(ranked_vec, "CC", params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    KEGG = mp_run_gsea_kegg(ranked_vec, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    Reactome = mp_run_gsea_reactome(ranked_vec, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    WikiPathways = mp_run_gsea_wikipathways(ranked_vec, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    Hallmark = mp_run_gsea_hallmark(ranked_vec, params$pvalueCutoff, params$minGSSize, params$maxGSSize),
    list(ok = FALSE, df = NULL, error = sprintf("Unknown database '%s'.", database))
  )
}

mp_annotate_reactome_hierarchy <- function(reactome_result_df) {
  if (is.null(reactome_result_df) || nrow(reactome_result_df) == 0) return(reactome_result_df)
  df <- reactome_result_df
  df$hierarchy_level <- NA_integer_; df$parent_pathway <- NA_character_; df$child_count <- NA_integer_
  if (!requireNamespace("httr", quietly = TRUE)) return(df)
  for (i in seq_len(nrow(df))) {
    stid <- df$ID[i]
    res <- tryCatch(httr::GET(sprintf("https://reactome.org/ContentService/data/pathway/%s/containedEvents", stid), httr::timeout(5)), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200) next
    body <- tryCatch(httr::content(res, as = "parsed", type = "application/json"), error = function(e) NULL)
    if (is.null(body)) next
    df$child_count[i] <- length(body)
    df$hierarchy_level[i] <- if (length(body) > 0) 1L else 0L
  }
  df
}

mp_run_topology <- function(genes_entrez, universe_entrez, params) {
  res <- mp_run_ora_reactome(genes_entrez, universe_entrez, params$pvalueCutoff, params$minGSSize, params$maxGSSize)
  if (!isTRUE(res$ok) || is.null(res$df)) return(res)
  sig <- res$df[!is.na(res$df$p.adjust) & res$df$p.adjust < (params$fdrCutoff %||% 1), , drop = FALSE]
  if (nrow(sig) == 0) return(res)
  annotated <- mp_annotate_reactome_hierarchy(utils::head(sig[order(sig$p.adjust), , drop = FALSE], 30))
  idx <- match(annotated$ID, res$df$ID)
  res$df$hierarchy_level <- NA_integer_; res$df$parent_pathway <- NA_character_; res$df$child_count <- NA_integer_
  res$df$hierarchy_level[idx] <- annotated$hierarchy_level
  res$df$child_count[idx] <- annotated$child_count
  res$df$method <- "Topology (Reactome ORA + hierarchy)"
  res
}

mp_build_evidence_tracks <- function(enrichment_df, input_df) {
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(enrichment_df)
  df <- enrichment_df
  df$transcript_gene_count <- NA_integer_; df$transcript_direction_summary <- NA_character_
  df$transcript_mean_effect <- NA_real_; df$transcript_min_p <- NA_real_
  df$meth_cpg_count <- NA_integer_; df$meth_mapped_gene_count <- NA_integer_
  df$meth_region_summary <- NA_character_; df$meth_direction_summary <- NA_character_; df$meth_min_p <- NA_real_
  df$integration_label <- NA_character_

  has_expr <- !is.na(input_df$expr_logFC)
  has_meth <- !is.na(input_df$cpg) & nzchar(input_df$cpg)

  for (i in seq_len(nrow(df))) {
    overlap_genes <- unique(trimws(unlist(strsplit(df$geneID[i] %||% "", "/"))))
    if (length(overlap_genes) == 0 || all(!nzchar(overlap_genes))) next
    hit <- toupper(input_df$gene_symbol) %in% toupper(overlap_genes) | toupper(input_df$feature) %in% toupper(overlap_genes)
    sub <- input_df[hit, , drop = FALSE]
    tx <- sub[has_expr[hit], , drop = FALSE]
    mt <- sub[has_meth[hit], , drop = FALSE]

    if (nrow(tx) > 0) {
      df$transcript_gene_count[i] <- nrow(tx)
      dirs <- unique(stats::na.omit(tx$expr_direction))
      df$transcript_direction_summary[i] <- if (length(dirs) == 1) dirs else if (length(dirs) > 1) "Mixed" else NA_character_
      df$transcript_mean_effect[i] <- round(mean(tx$expr_logFC, na.rm = TRUE), 3)
      df$transcript_min_p[i] <- suppressWarnings(min(tx$expr_p, na.rm = TRUE))
    }
    if (nrow(mt) > 0) {
      df$meth_cpg_count[i] <- length(unique(mt$cpg))
      df$meth_mapped_gene_count[i] <- length(unique(mt$gene_symbol[!is.na(mt$gene_symbol)]))
      df$meth_region_summary[i] <- paste(unique(stats::na.omit(mt$region)), collapse = "; ")
      dirs <- unique(stats::na.omit(mt$meth_direction))
      df$meth_direction_summary[i] <- if (length(dirs) == 1) dirs else if (length(dirs) > 1) "Mixed" else NA_character_
      df$meth_min_p[i] <- suppressWarnings(min(mt$meth_p, na.rm = TRUE))
    }
    df$integration_label[i] <- if (nrow(tx) > 0 && nrow(mt) > 0) "RNA + Methylation supported"
                                else if (nrow(tx) > 0) "Transcriptomics-only supported"
                                else if (nrow(mt) > 0) "Methylomics-only supported"
                                else NA_character_
  }
  df[is.infinite(df$transcript_min_p) | is.nan(df$transcript_min_p), "transcript_min_p"] <- NA_real_
  df[is.infinite(df$meth_min_p) | is.nan(df$meth_min_p), "meth_min_p"] <- NA_real_
  df
}

mp_mapping_direction <- function(row) {
  bp <- row$biological_pattern
  if (!is.null(bp) && !is.na(bp) && nzchar(bp)) {
    if (grepl("^concordant", bp, ignore.case = TRUE)) return("Directionally consistent (region-aware)")
    if (grepl("non-canonical", bp, ignore.case = TRUE)) return("Directionally opposite (region-aware)")
    return("Mixed")
  }
  ed <- row$expr_direction; md <- row$meth_direction
  if (is.null(ed) || is.null(md) || is.na(ed) || is.na(md) || !nzchar(ed) || !nzchar(md)) return("Insufficient information")
  up_like <- c("up", "upregulated", "hyper")
  down_like <- c("down", "downregulated", "hypo")
  ed_l <- tolower(ed); md_l <- tolower(md)
  if ((ed_l %in% up_like && md_l %in% up_like) || (ed_l %in% down_like && md_l %in% down_like)) return("Directionally consistent")
  if ((ed_l %in% up_like && md_l %in% down_like) || (ed_l %in% down_like && md_l %in% up_like)) return("Directionally opposite")
  "Mixed"
}

mp_validate_ora_inputs <- function(genes_entrez, universe_entrez) {
  if (is.null(genes_entrez) || length(genes_entrez) < 3) return(list(ok = FALSE, error = "Fewer than 3 identifiers could be mapped to Entrez IDs - not enough for a meaningful over-representation test."))
  if (!is.null(universe_entrez) && length(universe_entrez) < length(genes_entrez)) return(list(ok = FALSE, error = "Background/universe is smaller than the input gene list - switch Background to \"Entire selected database\"."))
  list(ok = TRUE, error = NULL)
}

mp_validate_gsea_inputs <- function(ranked_vec) {
  if (is.null(ranked_vec) || length(ranked_vec) < 10) return(list(ok = FALSE, error = "Fewer than 10 ranked, Entrez-mapped features are available - not enough for a stable GSEA run."))
  list(ok = TRUE, error = NULL)
}

mp_validate_database_choice <- function(database) {
  avail <- MP_DATABASES[[database]]$available %||% FALSE
  if (!isTRUE(avail)) return(list(ok = FALSE, error = sprintf("%s is not available in this deployment.", MP_DATABASES[[database]]$label %||% database)))
  list(ok = TRUE, error = NULL)
}

mp_infer_species <- function(id_type, ids) {
  if (identical(id_type, "Ensembl Gene ID") && length(ids) > 0) {
    if (mean(grepl("^ENSMUSG", ids)) > 0.5) return("Mus musculus")
    if (mean(grepl("^ENSRNOG", ids)) > 0.5) return("Rattus norvegicus")
  }
  "Homo sapiens"
}

mp_build_metadata <- function(database, method, species, background_label, ranking_method, input_count, mapped_count, padj_thresh, min_size, max_size, timestamp = NULL) {
  data.frame(
    Field = c("Database(s)", "Method", "Species", "Background/universe", "Ranking method", "Input features",
              "Mapped features", "Mapping rate (%)", "Adjusted P / FDR threshold", "Min gene-set size", "Max gene-set size", "Analyzed at"),
    Value = c(paste(database, collapse = ", "), method, species, background_label %||% "Not resolved",
              ranking_method %||% "Not applicable (ORA)", input_count, mapped_count,
              if (input_count > 0) round(100 * mapped_count / input_count, 1) else 0,
              padj_thresh, min_size, max_size, format(timestamp %||% Sys.time())),
    stringsAsFactors = FALSE
  )
}
