## R/methylomics/mod_methyl_biomarkercard.R
## Submodule: Biomarker Card - a single-CpG epigenomic profile page.
##
## CpG -> genomic location -> CpG island/regulatory context -> associated
## gene(s) -> evidence from the loaded dataset (incl. sex-specific) ->
## a cautious functional-interpretation chain -> a transparent evidence
## checklist -> sources -> downloadable report.
##
## Everything here is grounded in data already installed/bundled in this
## deployment (Illumina manifest annotation, ChAMPdata, TxDb/org.Hs.eg.db/
## GO.db, a bundled UCSC cytoband table, and this app's own preloaded/
## uploaded methylation results) - no live external-database calls (EWAS
## Catalog / EWAS Atlas / MethBank / KEGG / Reactome / PubMed) are made in
## this build; every place the full spec wants one of those, the UI shows
## an explicit "Not yet connected" placeholder rather than a guess.
##
## Self-contained data intake, same convention mod_methyl_candidates.R
## already established: this file does not read any other module's
## reactive state, and does not modify annotation.R even though it needs a
## richer extraction than methyl_get_annotation() - it reads the Illumina
## annotation packages' Locations/Other/Islands.UCSC objects directly here
## (same wrapper-avoidance rationale documented in annotation.R applies
## regardless of which file calls it).

## ---- Illumina manifest annotation (chr/pos/strand/gene/island/regulatory) ----

.bc_illumina_anno_cache <- new.env(parent = emptyenv())

## 450K's Other object exposes plain Enhancer/DHS flag columns; EPIC's Other
## object has no equivalent single column (it splits enhancer evidence
## across Phantom4/Phantom5/X450k_Enhancer and DHS across
## DNase_Hypersensitivity_NAME/_Evidence_Count) - verified directly against
## both installed packages. Rather than guess an EPIC equivalent, those two
## fields are only ever populated for 450K and are "Not available" for EPIC.
bc_illumina_full_annotation <- function(array_type) {
  pkg <- METHYL_ANNOTATION_PACKAGES[[array_type]]
  if (is.null(pkg)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf(
      "No Bioconductor manifest annotation is available for %s in this deployment.", array_type)))
  }
  cached <- .bc_illumina_anno_cache[[pkg]]
  if (!is.null(cached)) return(list(ok = TRUE, anno = cached, reason = NULL))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf("The %s annotation package is not installed in this deployment.", pkg)))
  }
  anno <- tryCatch({
    e <- new.env(parent = emptyenv())
    utils::data(list = c("Locations", "Other", "Islands.UCSC"), package = pkg, envir = e)
    loc <- as.data.frame(e$Locations)
    oth <- as.data.frame(e$Other)
    isl <- as.data.frame(e$Islands.UCSC)
    ids <- rownames(loc)
    enhancer <- if ("Enhancer" %in% colnames(oth)) ifelse(oth[ids, "Enhancer"] == "TRUE", "TRUE", NA_character_) else rep(NA_character_, length(ids))
    dhs <- if ("DHS" %in% colnames(oth)) ifelse(oth[ids, "DHS"] == "TRUE", "TRUE", NA_character_) else rep(NA_character_, length(ids))
    data.frame(
      chr = loc[ids, "chr"], pos = loc[ids, "pos"], strand = loc[ids, "strand"],
      gene_names = oth[ids, "UCSC_RefGene_Name"], gene_group = oth[ids, "UCSC_RefGene_Group"],
      regulatory_group = oth[ids, "Regulatory_Feature_Group"], regulatory_name = oth[ids, "Regulatory_Feature_Name"],
      enhancer = enhancer, dhs = dhs,
      island_name = isl[ids, "Islands_Name"], island_relation = isl[ids, "Relation_to_Island"],
      row.names = ids, stringsAsFactors = FALSE
    )
  }, error = function(e) e)
  if (inherits(anno, "error")) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf("Could not extract annotation from %s: %s", pkg, conditionMessage(anno))))
  }
  .bc_illumina_anno_cache[[pkg]] <- anno
  list(ok = TRUE, anno = anno, reason = NULL)
}

## Own copy of the ChAMPdata::probe.features extraction - deliberately not
## shared with mod_methyl_candidates.R's mcd_champ_full_annotation(), same
## "each module owns its own data intake" convention that file itself uses.
.bc_champ_anno_cache <- new.env(parent = emptyenv())
bc_champ_annotation <- function() {
  if (!is.null(.bc_champ_anno_cache$anno)) return(.bc_champ_anno_cache$anno)
  if (!requireNamespace("ChAMPdata", quietly = TRUE)) return(NULL)
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ utils::data("probe.features", package = "ChAMPdata", envir = e); TRUE }, error = function(e) FALSE)
  if (!ok || is.null(e$probe.features)) return(NULL)
  pf <- e$probe.features
  out <- data.frame(
    cpg = rownames(pf), chr = paste0("chr", pf$CHR), pos = pf$MAPINFO,
    gene = as.character(pf$gene), feature = as.character(pf$feature), island = as.character(pf$cgi),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$chr) & !is.na(out$pos), , drop = FALSE]
  .bc_champ_anno_cache$anno <- out
  out
}

BC_GROUP_PATTERNS <- c("^group$", "^status$", "^phenotype$", "^disease$", "^condition$", "^diagnosis$", "group")
BC_SEX_PATTERNS <- c("^sex$", "^gender$", "sex")
BC_ID_PATTERNS <- c("^sample$", "^gsm$", "^sample_id$", "^id$")

bc_find_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

## "chr16:53468284-53469209" -> list(chr,start,end,length) - the exact format
## Islands_Name uses in the Illumina annotation packages (verified directly).
bc_island_parse <- function(islands_name) {
  if (is.null(islands_name) || is.na(islands_name) || !nzchar(islands_name)) return(NULL)
  m <- regmatches(islands_name, regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", islands_name))[[1]]
  if (length(m) != 4) return(NULL)
  start <- as.numeric(m[3]); end <- as.numeric(m[4])
  list(chr = m[2], start = start, end = end, length = end - start + 1)
}

BC_FEATURE_MAP <- c(
  "TSS200" = "Promoter/TSS", "TSS1500" = "Promoter/TSS", "5'UTR" = "5' UTR",
  "1stExon" = "1st Exon", "Body" = "Gene body", "3'UTR" = "3' UTR", "ExonBnd" = "Exon boundary"
)
bc_feature_label <- function(gene_group) {
  if (is.null(gene_group) || is.na(gene_group) || !nzchar(gene_group)) return("Intergenic (no annotated RefGene overlap)")
  toks <- unique(trimws(strsplit(gene_group, ";")[[1]]))
  labs <- unique(stats::na.omit(unname(BC_FEATURE_MAP[toks])))
  if (length(labs) == 0) return(gene_group)
  paste(labs, collapse = " / ")
}

bc_island_context_label <- function(relation) {
  if (is.null(relation) || is.na(relation) || !nzchar(relation)) return("Open Sea (no annotated CpG island within range)")
  switch(relation,
    "Island" = "CpG Island", "N_Shore" = "Shore (North)", "S_Shore" = "Shore (South)",
    "N_Shelf" = "Shelf (North)", "S_Shelf" = "Shelf (South)", "OpenSea" = "Open Sea",
    relation)
}

bc_enhancer_flag <- function(r, array_type) {
  if (!identical(array_type, "450K")) return("Not available for this array type in this build")
  if (is.na(r$enhancer)) return("No annotated enhancer overlap")
  "Yes - annotated enhancer region (450K manifest)"
}
bc_dhs_flag <- function(r, array_type) {
  if (!identical(array_type, "450K")) return("Not available for this array type in this build")
  if (is.na(r$dhs)) return("No annotated DHS overlap")
  "Yes - annotated DNase hypersensitive site (450K manifest)"
}

## Merges the Illumina manifest row and the ChAMP row for one probe -
## primary source of truth for coordinates is the Illumina manifest when
## available (it's the array-specific, build-time source); ChAMP fills in
## chr/pos/gene/feature/island when the manifest package isn't available.
bc_resolve_cpg <- function(cpg_id, array_type) {
  a <- bc_illumina_full_annotation(array_type)
  champ <- bc_champ_annotation()
  illumina_row <- if (isTRUE(a$ok) && cpg_id %in% rownames(a$anno)) a$anno[cpg_id, ] else NULL
  champ_row <- if (!is.null(champ) && cpg_id %in% champ$cpg) champ[champ$cpg == cpg_id, , drop = FALSE][1, ] else NULL
  list(
    cpg = cpg_id, found_in_manifest = !is.null(illumina_row), found_in_champ = !is.null(champ_row),
    chr = if (!is.null(illumina_row)) illumina_row$chr else if (!is.null(champ_row)) champ_row$chr else NA_character_,
    pos = if (!is.null(illumina_row)) illumina_row$pos else if (!is.null(champ_row)) champ_row$pos else NA_real_,
    strand = if (!is.null(illumina_row)) illumina_row$strand else NA_character_,
    gene_names = if (!is.null(illumina_row)) illumina_row$gene_names else NA_character_,
    gene_group = if (!is.null(illumina_row)) illumina_row$gene_group else NA_character_,
    regulatory_group = if (!is.null(illumina_row)) illumina_row$regulatory_group else NA_character_,
    regulatory_name = if (!is.null(illumina_row)) illumina_row$regulatory_name else NA_character_,
    enhancer = if (!is.null(illumina_row)) illumina_row$enhancer else NA_character_,
    dhs = if (!is.null(illumina_row)) illumina_row$dhs else NA_character_,
    island_name = if (!is.null(illumina_row)) illumina_row$island_name else NA_character_,
    island_relation = if (!is.null(illumina_row)) illumina_row$island_relation else if (!is.null(champ_row)) champ_row$island else NA_character_,
    champ_gene = if (!is.null(champ_row)) champ_row$gene else NA_character_,
    champ_feature = if (!is.null(champ_row)) champ_row$feature else NA_character_
  )
}

## Case-insensitive whole-token search ("(^|;)SYMBOL(;|$)") across the
## semicolon-joined UCSC_RefGene_Name list - regex over the full ~485k/865k
## row vector is fast; strsplit-per-row would not be.
bc_search_gene <- function(symbol, array_type, max_n = 300) {
  a <- bc_illumina_full_annotation(array_type)
  if (!isTRUE(a$ok)) return(list(ok = FALSE, reason = a$reason, df = NULL))
  anno <- a$anno
  esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", symbol)
  pattern <- sprintf("(^|;)%s(;|$)", esc)
  hit <- !is.na(anno$gene_names) & grepl(pattern, anno$gene_names, ignore.case = TRUE)
  ids <- rownames(anno)[hit]
  if (length(ids) == 0) return(list(ok = TRUE, df = data.frame(cpg = character(0), chr = character(0), pos = numeric(0),
                                                                 gene_names = character(0), gene_group = character(0),
                                                                 island_relation = character(0), stringsAsFactors = FALSE)))
  out <- anno[ids, c("chr", "pos", "gene_names", "gene_group", "island_relation"), drop = FALSE]
  out$cpg <- ids
  out <- out[order(out$pos), , drop = FALSE]
  if (nrow(out) > max_n) out <- out[seq_len(max_n), , drop = FALSE]
  list(ok = TRUE, df = out)
}

## For the multi-gene relationship table: UCSC_RefGene_Name and
## UCSC_RefGene_Group are parallel semicolon-separated lists (one entry per
## overlapping transcript) in the Illumina manifest's Other object.
bc_gene_relationship_table <- function(resolved) {
  names_str <- resolved$gene_names; group_str <- resolved$gene_group
  if (is.null(names_str) || is.na(names_str) || !nzchar(names_str)) return(NULL)
  nm <- trimws(strsplit(names_str, ";")[[1]])
  gp <- if (!is.null(group_str) && !is.na(group_str) && nzchar(group_str)) trimws(strsplit(group_str, ";")[[1]]) else character(0)
  n <- length(nm)
  if (length(gp) < n) gp <- c(gp, rep(NA_character_, n - length(gp)))
  df <- data.frame(
    Gene = nm,
    Relationship = vapply(gp[seq_len(n)], function(x) if (is.na(x) || !nzchar(x)) "Overlapping (unspecified)" else bc_feature_label(x), character(1)),
    stringsAsFactors = FALSE
  )
  df[!duplicated(df$Gene), , drop = FALSE]
}

## ---- Gene structure (real NCBI/Ensembl/exon data, no live network call) ----

.bc_txdb_genes_cache <- new.env(parent = emptyenv())
bc_txdb_genes <- function() {
  if (!is.null(.bc_txdb_genes_cache$g)) return(.bc_txdb_genes_cache$g)
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE) || !requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
  g <- tryCatch(suppressWarnings(suppressMessages(
    GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  )), error = function(e) NULL)
  .bc_txdb_genes_cache$g <- g
  g
}

.bc_txdb_exons_cache <- new.env(parent = emptyenv())
bc_txdb_exons_by_gene <- function() {
  if (!is.null(.bc_txdb_exons_cache$ex)) return(.bc_txdb_exons_cache$ex)
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE) || !requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
  ex <- tryCatch(GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene"), error = function(e) NULL)
  .bc_txdb_exons_cache$ex <- ex
  ex
}

## symbol -> NCBI Gene ID / Ensembl ID / gene name (org.Hs.eg.db) + real
## exon/gene-span/strand structure (TxDb.Hsapiens.UCSC.hg19.knownGene) -
## fail-soft: any missing package/unmapped symbol returns ok=FALSE with a
## reason, never a guessed coordinate.
bc_gene_structure <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No associated gene symbol to resolve."))
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(list(ok = FALSE, reason = "org.Hs.eg.db is not installed in this deployment."))
  res <- tryCatch({
    map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL",
                                                    columns = c("ENTREZID", "ENSEMBL", "GENENAME")))
    map <- map[!is.na(map$ENTREZID), , drop = FALSE]
    if (nrow(map) == 0) return(list(ok = FALSE, reason = sprintf("No NCBI Entrez Gene entry found for symbol \"%s\".", symbol)))
    entrez <- map$ENTREZID[1]
    g <- bc_txdb_genes()
    gene_span <- if (!is.null(g) && entrez %in% names(g)) g[entrez] else NULL
    ex_all <- bc_txdb_exons_by_gene()
    ex <- if (!is.null(ex_all) && entrez %in% names(ex_all)) ex_all[[entrez]] else NULL
    list(ok = TRUE, symbol = symbol, entrez = entrez, ensembl = map$ENSEMBL[1], genename = map$GENENAME[1],
         chr = if (!is.null(gene_span)) as.character(GenomicRanges::seqnames(gene_span))[1] else NA_character_,
         start = if (!is.null(gene_span)) BiocGenerics::start(gene_span)[1] else NA_real_,
         end = if (!is.null(gene_span)) BiocGenerics::end(gene_span)[1] else NA_real_,
         strand = if (!is.null(gene_span)) as.character(BiocGenerics::strand(gene_span))[1] else NA_character_,
         exons = ex)
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, reason = sprintf("Could not resolve gene structure: %s", conditionMessage(res))))
  res
}

bc_go_terms <- function(entrez, n = 5) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(NULL)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
    go <- go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", , drop = FALSE]
    ids <- unique(go$GOALL)
    if (length(ids) == 0) return(NULL)
    terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, keytype = "GOID", columns = "TERM"))
    utils::head(unique(stats::na.omit(terms$TERM)), n)
  }, error = function(e) NULL)
  if (inherits(res, "error")) return(NULL)
  res
}

## ---- UCSC cytoband (bundled, fetched once from hgdownload.soe.ucsc.edu) ----

BC_CYTOBAND_PATH <- get_reference_path("cytoBandIdeo_hg19.txt.gz")
.bc_cytoband_cache <- new.env(parent = emptyenv())
bc_cytoband_hg19 <- function() {
  if (!is.null(.bc_cytoband_cache$df)) return(.bc_cytoband_cache$df)
  if (!file.exists(BC_CYTOBAND_PATH)) return(NULL)
  df <- tryCatch(
    utils::read.delim(gzfile(BC_CYTOBAND_PATH), header = FALSE,
                       col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
                       stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  .bc_cytoband_cache$df <- df
  df
}

BC_GIESTAIN_COLORS <- c(
  gneg = "#FFFFFF", gpos25 = "#C8C8C8", gpos50 = "#969696", gpos75 = "#5A5A5A", gpos100 = "#141414",
  acen = "#B33A3A", gvar = "#A9A9A9", stalk = "#6E8FB0"
)

## ---- Preloaded DMP/DMR lookups (process-wide cache; static preloaded files) ----

.bc_dmp_cache <- new.env(parent = emptyenv())
bc_dmp_table <- function(sex, stage = "sva") {
  key <- paste(stage, sex, sep = "_")
  cached <- .bc_dmp_cache[[key]]
  if (!is.null(cached)) return(cached)
  df <- load_default_dmp(stage, sex)
  if (!is.null(df)) .bc_dmp_cache[[key]] <- df
  df
}
.bc_dmr_cache <- new.env(parent = emptyenv())
bc_dmr_table <- function(sex) {
  cached <- .bc_dmr_cache[[sex]]
  if (!is.null(cached)) return(cached)
  df <- load_default_dmr(sex)
  if (!is.null(df)) .bc_dmr_cache[[sex]] <- df
  df
}
bc_dmp_lookup <- function(cpg, sex, stage = "sva") {
  df <- bc_dmp_table(sex, stage)
  if (is.null(df)) return(NULL)
  row <- df[df$cpg == cpg, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  row[1, ]
}
bc_dmr_lookup <- function(chr, pos, sex) {
  df <- bc_dmr_table(sex)
  if (is.null(df) || is.null(chr) || is.na(chr) || is.null(pos) || is.na(pos)) return(NULL)
  chr_norm <- sub("^chr", "", chr, ignore.case = TRUE)
  seq_norm <- sub("^chr", "", as.character(df$seqnames), ignore.case = TRUE)
  hit <- df[seq_norm == chr_norm & df$start <= pos & df$end >= pos, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit
}

## ---- Live/uploaded-dataset evidence (Welch t-test on the actual beta matrix) ----

bc_pick_case_control <- function(levels) {
  lv <- levels
  case_idx <- which(grepl("case|patient|disease|^ra$|affected|positive", lv, ignore.case = TRUE))
  ctrl_idx <- which(grepl("control|healthy|normal|^hc$|negative", lv, ignore.case = TRUE))
  if (length(case_idx) >= 1 && length(ctrl_idx) >= 1 && case_idx[1] != ctrl_idx[1]) {
    return(list(case = lv[case_idx[1]], control = lv[ctrl_idx[1]]))
  }
  if (length(ctrl_idx) >= 1) {
    other <- setdiff(seq_along(lv), ctrl_idx[1])
    if (length(other) >= 1) return(list(case = lv[other[1]], control = lv[ctrl_idx[1]]))
  }
  list(case = lv[2], control = lv[1])
}

bc_live_stats <- function(beta_row, group_vec, case_label, control_label) {
  keep <- !is.na(beta_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  beta_row <- beta_row[keep]; grp <- group_vec[keep]
  case_vals <- beta_row[grp == case_label]; ctrl_vals <- beta_row[grp == control_label]
  if (length(case_vals) < 2 || length(ctrl_vals) < 2) {
    return(list(ok = FALSE, reason = "Fewer than 2 non-missing samples in one of the two groups - cannot compute a t-test.",
                case_label = case_label, control_label = control_label, n_case = length(case_vals), n_control = length(ctrl_vals)))
  }
  tt <- tryCatch(stats::t.test(case_vals, ctrl_vals), error = function(e) NULL)
  mc <- mean(case_vals); mo <- mean(ctrl_vals)
  list(ok = TRUE, mean_case = mc, mean_control = mo, dbeta = mc - mo,
       p_value = if (!is.null(tt)) tt$p.value else NA_real_,
       ci = if (!is.null(tt)) as.numeric(tt$conf.int) else c(NA_real_, NA_real_),
       n_case = length(case_vals), n_control = length(ctrl_vals),
       direction = if (mc > mo) "Hypermethylated" else if (mc < mo) "Hypomethylated" else "No change",
       case_label = case_label, control_label = control_label)
}

## Works identically whether `dataset` holds the preloaded live matrix or a
## user-uploaded one - methyl_dataset$beta/$sample_sheet is the same shape
## either way (see mod_methyl_dataset.R), so this needs no source-specific
## branching.
bc_dataset_evidence <- function(cpg, dataset) {
  beta <- dataset$beta; sheet <- dataset$sample_sheet
  if (is.null(beta) || is.null(sheet)) return(list(ok = FALSE, reason = "No beta matrix / sample sheet is currently loaded on the Dataset tab."))
  if (!cpg %in% rownames(beta)) return(list(ok = FALSE, reason = "This CpG is not present in the currently loaded beta matrix."))
  samp_ids <- colnames(beta)
  id_col <- bc_find_col(colnames(sheet), BC_ID_PATTERNS)
  if (!is.null(id_col) && all(samp_ids %in% as.character(sheet[[id_col]]))) {
    sheet <- sheet[match(samp_ids, as.character(sheet[[id_col]])), , drop = FALSE]
  } else if (nrow(sheet) != length(samp_ids)) {
    return(list(ok = FALSE, reason = "Could not align the sample sheet to the beta matrix columns (no shared sample-ID column and row counts differ)."))
  }
  beta_row <- beta[cpg, ]
  group_col <- bc_find_col(colnames(sheet), BC_GROUP_PATTERNS)
  if (is.null(group_col)) return(list(ok = FALSE, reason = "Could not auto-detect a case/control grouping column in the loaded sample sheet.", beta_row = beta_row))
  group_vec <- as.character(sheet[[group_col]])
  levels_present <- unique(group_vec[!is.na(group_vec)])
  if (length(levels_present) != 2) {
    return(list(ok = FALSE, reason = sprintf("The detected grouping column (\"%s\") does not have exactly two levels (found: %s).",
                                              group_col, paste(levels_present, collapse = ", ")), beta_row = beta_row))
  }
  cc <- bc_pick_case_control(levels_present)
  overall <- bc_live_stats(beta_row, group_vec, cc$case, cc$control)
  sex_col <- bc_find_col(colnames(sheet), BC_SEX_PATTERNS)
  sex_vec <- if (!is.null(sex_col)) as.character(sheet[[sex_col]]) else NULL
  by_sex <- NULL
  if (!is.null(sex_vec)) {
    sex_levels <- unique(sex_vec[!is.na(sex_vec)])
    if (length(sex_levels) > 0) {
      by_sex <- stats::setNames(lapply(sex_levels, function(s) {
        idx <- which(sex_vec == s)
        bc_live_stats(beta_row[idx], group_vec[idx], cc$case, cc$control)
      }), sex_levels)
    }
  }
  list(ok = TRUE, beta_row = beta_row, group_vec = group_vec, sex_vec = sex_vec,
       group_col = group_col, sex_col = sex_col, case_label = cc$case, control_label = cc$control,
       overall = overall, by_sex = by_sex)
}

## ---- External database lookups (EWAS Catalog, EWAS Atlas, KEGG, Reactome,
## PubMed) - live, keyless REST calls, each independently fail-soft. MethBank
## has no public API (confirmed) and is handled as an outbound link only,
## never a live query. Every one of these is opt-in (its own button, not part
## of card_data()) so a slow/down external service never blocks the fully
## local, already-working core of the card. ----------------------------------

## Column names come from the response's own `fields` array (never hardcoded)
## so a schema change on their end degrades to a parse failure, not silently
## mislabeled data.
bc_ewascatalog_query <- function(type, value) {
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, df = NULL, reason = "httr2 is not installed."))
  url <- sprintf("https://www.ewascatalog.org/api/?%s=%s", type, utils::URLencode(value, reserved = TRUE))
  res <- tryCatch({
    resp <- httr2::request(url) %>% httr2::req_timeout(15) %>% httr2::req_perform()
    httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, df = NULL, reason = sprintf("EWAS Catalog lookup failed: %s", conditionMessage(res))))
  fields <- unlist(res$fields)
  rows <- res$results
  if (is.null(fields)) return(list(ok = FALSE, df = NULL, reason = "EWAS Catalog response did not include a field schema."))
  if (length(rows) == 0) {
    empty <- as.data.frame(matrix(character(0), ncol = length(fields)), stringsAsFactors = FALSE)
    colnames(empty) <- fields
    return(list(ok = TRUE, df = empty, reason = NULL))
  }
  mat <- do.call(rbind, lapply(rows, function(r) {
    r <- lapply(r, function(x) if (is.null(x)) NA_character_ else as.character(x))
    unlist(r)
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- fields
  rownames(df) <- NULL
  list(ok = TRUE, df = df, reason = NULL)
}

bc_ewasatlas_query <- function(endpoint, param, value) {
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, data = NULL, reason = "httr2 is not installed."))
  url <- sprintf("https://ngdc.cncb.ac.cn/ewas/rest/%s?%s=%s", endpoint, param, utils::URLencode(value, reserved = TRUE))
  res <- tryCatch({
    resp <- httr2::request(url) %>% httr2::req_timeout(15) %>% httr2::req_perform()
    httr2::resp_body_json(resp, simplifyVector = TRUE)
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, data = NULL, reason = sprintf("EWAS Atlas lookup failed: %s", conditionMessage(res))))
  code <- res$code
  if (!is.null(code) && !identical(as.integer(code), 0L)) return(list(ok = FALSE, data = NULL, reason = res$msg %||% "EWAS Atlas returned an error."))
  list(ok = TRUE, data = res$data, reason = NULL)
}
bc_ewasatlas_probe <- function(probe_id) bc_ewasatlas_query("probe", "probeId", probe_id)
bc_ewasatlas_gene <- function(gene_symbol) bc_ewasatlas_query("gene", "geneSymbol", gene_symbol)

## keggList("pathway","hsa") is a single ~370-row call, cached process-wide,
## used to resolve pathway names without one keggGet() round trip per pathway.
.bc_kegg_pathway_names_cache <- new.env(parent = emptyenv())
bc_kegg_pathway_names <- function() {
  if (!is.null(.bc_kegg_pathway_names_cache$v)) return(.bc_kegg_pathway_names_cache$v)
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  v <- tryCatch(KEGGREST::keggList("pathway", "hsa"), error = function(e) NULL)
  .bc_kegg_pathway_names_cache$v <- v
  v
}
bc_kegg_pathways_for_gene <- function(entrez) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, pathways = NULL, reason = "No NCBI Gene ID available for KEGG lookup."))
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(list(ok = FALSE, pathways = NULL, reason = "KEGGREST is not installed in this deployment."))
  res <- tryCatch({
    links <- KEGGREST::keggLink("pathway", sprintf("hsa:%s", entrez))
    key <- sub("^path:", "", unname(links))
    names_map <- bc_kegg_pathway_names()
    nm <- if (!is.null(names_map)) unname(names_map[key]) else rep(NA_character_, length(key))
    nm <- sub(" - Homo sapiens.*$", "", nm)
    data.frame(id = key, name = nm, stringsAsFactors = FALSE)
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, pathways = NULL, reason = sprintf("KEGG lookup failed: %s", conditionMessage(res))))
  list(ok = TRUE, pathways = res, reason = NULL)
}

## hsa05323's own gene list, fetched and parsed once (process-wide cache) -
## direct membership check rather than trusting keggLink() (which in testing
## returned empty for some genes that plausibly should show pathway links,
## so the RA-specific spotlight section checks the pathway's own gene list
## directly instead of relying solely on the general per-gene lookup above).
.bc_kegg_ra_pathway_cache <- new.env(parent = emptyenv())
bc_kegg_ra_pathway_genes <- function() {
  if (!is.null(.bc_kegg_ra_pathway_cache$obj)) return(.bc_kegg_ra_pathway_cache$obj)
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    p <- KEGGREST::keggGet("hsa05323")[[1]]
    genes <- p$GENE
    ids <- genes[seq(1, length(genes), by = 2)]
    descs <- genes[seq(2, length(genes), by = 2)]
    syms <- sub(";.*", "", descs)
    list(df = data.frame(entrez = ids, symbol = syms, stringsAsFactors = FALSE),
         name = paste(p$NAME, collapse = " "), description = paste(p$DESCRIPTION, collapse = " "))
  }, error = function(e) e)
  if (inherits(res, "error")) return(NULL)
  .bc_kegg_ra_pathway_cache$obj <- res
  res
}
bc_kegg_ra_pathway_check <- function(entrez, symbol) {
  ra <- bc_kegg_ra_pathway_genes()
  if (is.null(ra)) return(list(ok = FALSE, reason = "KEGG RA pathway lookup failed or is unavailable in this deployment."))
  in_by_id <- !is.null(entrez) && !is.na(entrez) && entrez %in% ra$df$entrez
  in_by_symbol <- !is.null(symbol) && !is.na(symbol) && toupper(symbol) %in% toupper(ra$df$symbol)
  list(ok = TRUE, in_pathway = in_by_id || in_by_symbol, pathway_name = ra$name, pathway_description = ra$description, genes = ra$df)
}

## org.Hs.eg.db's SYMBOL->UNIPROT mapping is 1:many (isoforms/entries); tries
## each until one has Reactome pathway data. NOTE the inner `NULL` (not
## `return(NULL)`) on a non-200 response - `return()` inside a tryCatch()
## expression exits the *enclosing function*, not just that expression, which
## would silently abort the whole loop on the first non-canonical UniProt ID.
bc_reactome_pathways_for_gene <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, pathways = NULL, reason = "No gene symbol available for Reactome lookup."))
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("httr2", quietly = TRUE)) {
    return(list(ok = FALSE, pathways = NULL, reason = "org.Hs.eg.db / httr2 not available in this deployment."))
  }
  uniprot_ids <- tryCatch({
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = "UNIPROT"))
    unique(stats::na.omit(m$UNIPROT))
  }, error = function(e) character(0))
  if (length(uniprot_ids) == 0) return(list(ok = FALSE, pathways = NULL, reason = sprintf("No UniProt ID found for gene symbol \"%s\".", symbol)))
  for (uid in uniprot_ids) {
    res <- tryCatch({
      url <- sprintf("https://reactome.org/ContentService/data/mapping/UniProt/%s/pathways?species=9606", uid)
      resp <- httr2::request(url) %>% httr2::req_timeout(15) %>%
        httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
      if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = TRUE)
    }, error = function(e) NULL)
    if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
      out <- res[, intersect(c("stId", "displayName"), colnames(res)), drop = FALSE]
      return(list(ok = TRUE, pathways = out, reason = NULL, uniprot = uid))
    }
  }
  list(ok = TRUE, pathways = data.frame(stId = character(0), displayName = character(0)), reason = NULL)
}

## Same NCBI E-utilities pubmed_search() already uses (global.R) but keyed by
## exact PMIDs (from EWAS Catalog/Atlas results) instead of a keyword search.
bc_pubmed_summaries <- function(pmids) {
  pmids <- unique(pmids[!is.na(pmids) & nzchar(pmids)])
  empty <- data.frame(pmid = character(0), title = character(0), journal = character(0), year = character(0), stringsAsFactors = FALSE)
  if (length(pmids) == 0) return(list(ok = TRUE, table = empty, reason = NULL))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, table = NULL, reason = "httr2 is not installed."))
  pmids <- utils::head(pmids, 20)
  res <- tryCatch({
    httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "pubmed", id = paste(pmids, collapse = ","), retmode = "json", tool = "arthomix-explorer") %>%
      httr2::req_timeout(15) %>%
      httr2::req_perform() %>%
      httr2::resp_body_json()
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, table = NULL, reason = sprintf("PubMed lookup failed: %s", conditionMessage(res))))
  ids <- unlist(res$result$uids)
  if (is.null(ids)) return(list(ok = TRUE, table = empty, reason = NULL))
  rows <- lapply(ids, function(id) {
    rec <- res$result[[id]]
    data.frame(pmid = id, title = rec$title %||% NA_character_,
               journal = rec$fulljournalname %||% rec$source %||% NA_character_,
               year = substr(rec$pubdate %||% "", 1, 4), stringsAsFactors = FALSE)
  })
  list(ok = TRUE, table = do.call(rbind, rows), reason = NULL)
}

## No public API (confirmed) - outbound link only, never presented as a
## queried result.
bc_methbank_link <- function() "https://ngdc.cncb.ac.cn/methbank/"

BC_DISEASE_CATEGORIES <- list(
  "Rheumatoid Arthritis" = "rheumatoid arthritis|\\bra\\b",
  "Osteoarthritis" = "osteoarthritis|\\boa\\b",
  "SLE / Lupus" = "lupus|\\bsle\\b",
  "Psoriasis" = "psoriasi",
  "Inflammatory Bowel Disease" = "inflammatory bowel|crohn|ulcerative colitis|\\bibd\\b",
  "Type 1 Diabetes" = "type 1 diabetes|\\bt1d\\b",
  "Cancer" = "cancer|carcinoma|tumou?r|neoplasm|leukemia|lymphoma",
  "Cardiovascular Disease" = "cardiovascular|coronary|myocardial|atherosclero|heart failure"
)

## One call per source per CpG/gene; every sub-lookup is independently
## fail-soft so one down/slow service degrades only its own section.
bc_external_evidence <- function(cpg, gene_symbol, entrez) {
  ewc_cpg <- bc_ewascatalog_query("cpg", cpg)
  ewc_gene <- if (!is.na(gene_symbol)) bc_ewascatalog_query("gene", gene_symbol) else list(ok = FALSE, df = NULL, reason = "No associated gene.")
  ewa <- bc_ewasatlas_probe(cpg)
  kegg <- bc_kegg_pathways_for_gene(entrez)
  kegg_ra <- bc_kegg_ra_pathway_check(entrez, gene_symbol)
  reactome <- if (!is.na(gene_symbol)) bc_reactome_pathways_for_gene(gene_symbol) else list(ok = FALSE, pathways = NULL, reason = "No associated gene.")

  ewa_assoc <- NULL
  if (isTRUE(ewa$ok) && !is.null(ewa$data$associationList) && is.data.frame(ewa$data$associationList) && nrow(ewa$data$associationList) > 0) {
    a <- ewa$data$associationList
    ewa_assoc <- data.frame(source = "EWAS Atlas", trait = a$trait, tissue = NA_character_,
                             effect = a$correlation, p = NA_character_, pmid = as.character(a$pmid), stringsAsFactors = FALSE)
  }
  ewc_assoc <- NULL
  if (isTRUE(ewc_cpg$ok) && nrow(ewc_cpg$df) > 0) {
    d <- ewc_cpg$df
    ewc_assoc <- data.frame(source = "EWAS Catalog", trait = d$trait, tissue = d$tissue,
                             effect = d$beta, p = d$p, pmid = as.character(d$pmid), stringsAsFactors = FALSE)
  }
  combined <- do.call(rbind, list(ewc_assoc, ewa_assoc))
  if (!is.null(combined)) combined <- combined[!is.na(combined$trait) & nzchar(combined$trait), , drop = FALSE]

  pmids <- if (!is.null(combined)) unique(combined$pmid[!is.na(combined$pmid) & grepl("^[0-9]+$", combined$pmid)]) else character(0)
  pubs <- bc_pubmed_summaries(pmids)

  ra_rows <- if (!is.null(combined)) combined[grepl("rheumatoid|arthritis", combined$trait, ignore.case = TRUE), , drop = FALSE] else NULL

  disease_counts <- NULL
  if (!is.null(combined) && nrow(combined) > 0) {
    cat_counts <- vapply(BC_DISEASE_CATEGORIES, function(pat) sum(grepl(pat, combined$trait, ignore.case = TRUE)), integer(1))
    other_n <- nrow(combined) - sum(cat_counts)
    disease_counts <- data.frame(category = c(names(BC_DISEASE_CATEGORIES), "Other"), n = c(unname(cat_counts), max(other_n, 0)), stringsAsFactors = FALSE)
    disease_counts <- disease_counts[disease_counts$n > 0, , drop = FALSE]
  }

  tissue_counts <- NULL
  if (isTRUE(ewc_cpg$ok) && nrow(ewc_cpg$df) > 0) {
    tt <- ewc_cpg$df$tissue[!is.na(ewc_cpg$df$tissue) & nzchar(ewc_cpg$df$tissue)]
    if (length(tt) > 0) tissue_counts <- as.data.frame(table(tissue = tt), stringsAsFactors = FALSE)
  }

  list(cpg = cpg, gene_symbol = gene_symbol, ewascatalog_cpg = ewc_cpg, ewascatalog_gene = ewc_gene, ewasatlas = ewa,
       kegg = kegg, kegg_ra = kegg_ra, reactome = reactome, combined_disease = combined, ra_rows = ra_rows,
       disease_counts = disease_counts, tissue_counts = tissue_counts, publications = pubs,
       methbank_link = bc_methbank_link())
}

## ---- Plots ---------------------------------------------------------------

bc_plot_region <- function(chr, pos, island, gene_struct, flank_bp) {
  validate(need(!is.null(chr) && !is.na(chr) && !is.null(pos) && !is.na(pos), "No genomic coordinate is available to plot for this CpG."))
  xmin <- pos - flank_bp; xmax <- pos + flank_bp
  p <- ggplot() +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, 1)) +
    theme_arthomix() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = sprintf("%s position (bp) - hg19", chr), y = NULL)

  if (!is.null(island) && identical(island$chr, chr)) {
    isl_df <- data.frame(xmin = island$start, xmax = island$end, ymin = 0.55, ymax = 0.78)
    p <- p + geom_rect(data = isl_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                        fill = ARTHOMIX_COLORS$aqua, alpha = 0.55, inherit.aes = FALSE) +
      annotate("text", x = (island$start + island$end) / 2, y = 0.88, label = "CpG Island", size = 3.1, color = ARTHOMIX_COLORS$aqua)
  }

  if (isTRUE(gene_struct$ok) && !is.na(gene_struct$chr) && identical(gene_struct$chr, chr)) {
    gy <- 0.22
    gene_df <- data.frame(x = gene_struct$start, xend = gene_struct$end, y = gy, yend = gy)
    p <- p + geom_segment(data = gene_df, aes(x = x, xend = xend, y = y, yend = yend),
                           color = ARTHOMIX_COLORS$ink_muted, linewidth = 1, inherit.aes = FALSE)
    if (!is.null(gene_struct$exons) && length(gene_struct$exons) > 0) {
      ex_df <- as.data.frame(gene_struct$exons)
      ex_df <- ex_df[ex_df$end >= xmin & ex_df$start <= xmax, , drop = FALSE]
      if (nrow(ex_df) > 0) {
        ex_df$ymin <- gy - 0.05; ex_df$ymax <- gy + 0.05
        p <- p + geom_rect(data = ex_df, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
                            fill = ARTHOMIX_COLORS$blue, inherit.aes = FALSE)
      }
    }
    p <- p + annotate("text", x = xmin, y = gy + 0.14, label = gene_struct$symbol, hjust = 0, size = 3.1,
                       fontface = "italic", color = ARTHOMIX_COLORS$ink_muted)
  }

  cpg_df <- data.frame(x = pos, y = 0.42)
  p <- p + geom_point(data = cpg_df, aes(x = x, y = y), color = ARTHOMIX_STATUS$critical, size = 3.2, inherit.aes = FALSE) +
    annotate("text", x = pos, y = 0.32, label = "CpG", size = 3, color = ARTHOMIX_STATUS$critical)
  p
}

bc_plot_ideogram <- function(chr, pos, cytoband_df) {
  validate(need(!is.null(cytoband_df), "The bundled UCSC cytoband reference (hg19) is not available in this deployment."))
  validate(need(!is.null(chr) && !is.na(chr), "No chromosome is available to plot."))
  d <- cytoband_df[cytoband_df$chrom == chr, , drop = FALSE]
  validate(need(nrow(d) > 0, sprintf("No cytoband data available for %s.", chr)))
  d$fillcol <- BC_GIESTAIN_COLORS[d$gieStain]
  d$fillcol[is.na(d$fillcol)] <- "#DDDDDD"
  mk_df <- data.frame(x = pos, y = 0)
  ggplot() +
    geom_rect(data = d, aes(xmin = chromStart, xmax = chromEnd, ymin = -0.35, ymax = 0.35),
              fill = d$fillcol, color = "grey45", linewidth = 0.15) +
    geom_point(data = mk_df, aes(x = x, y = y), color = ARTHOMIX_STATUS$critical, size = 4, shape = 24, fill = ARTHOMIX_STATUS$critical) +
    annotate("text", x = pos, y = 0.62, label = sprintf("%s:%s", chr, format(round(pos), big.mark = ",")),
             size = 3.1, color = ARTHOMIX_STATUS$critical, fontface = "bold") +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = sprintf("%s position (bp) - hg19 cytoBandIdeo", chr), y = NULL, title = sprintf("Chromosome %s", sub("^chr", "", chr))) +
    theme_arthomix() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
}

bc_plot_methylation_dist <- function(df, facet_sex = FALSE) {
  if (isTRUE(facet_sex) && "sex" %in% names(df)) df <- df[!is.na(df$sex), , drop = FALSE]
  df <- df[!is.na(df$beta) & !is.na(df$group), , drop = FALSE]
  validate(need(nrow(df) > 0, "No non-missing beta values are available for this CpG in the current groups."))
  p <- ggplot(df, aes(x = group, y = beta, fill = group)) +
    geom_violin(alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.1) +
    scale_fill_manual(values = arthomix_pair(sort(unique(df$group)))) +
    labs(x = NULL, y = "Beta value (0-1)") + theme_arthomix() + theme(legend.position = "none")
  if (isTRUE(facet_sex) && "sex" %in% names(df)) p <- p + facet_wrap(~sex)
  p
}

bc_ggsave_datauri <- function(plot, width = 8, height = 3.2, dpi = 110) {
  if (is.null(plot) || !requireNamespace("base64enc", quietly = TRUE)) return(NULL)
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggplot2::ggsave(tf, plot = plot, width = width, height = height, dpi = dpi, bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok || !file.exists(tf)) return(NULL)
  uri <- base64enc::dataURI(file = tf, mime = "image/png")
  unlink(tf)
  uri
}

## ---- Display helpers -------------------------------------------------------

bc_fmt_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("Not available")
  if (length(x) == 1 && is.na(x)) return("Not available")
  if (is.character(x) && !nzchar(trimws(x))) return("Not available")
  as.character(x)
}
bc_kv_table <- function(pairs) {
  df <- data.frame(Field = names(pairs), Value = vapply(pairs, bc_fmt_field, character(1)), stringsAsFactors = FALSE)
  DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
}

## ---- Section builders (plain tags - reused identically on-screen and in
## the downloadable report) --------------------------------------------------

bc_section_summary <- function(d) {
  r <- d$resolved
  cand <- list(Female = d$dmp_female, Male = d$dmp_male)
  cand <- cand[!vapply(cand, is.null, logical(1))]
  headline <- NULL; headline_label <- NULL
  if (length(cand) > 0) {
    fdrs <- vapply(cand, function(x) if (is.na(x$fdr_bacon)) Inf else x$fdr_bacon, numeric(1))
    headline_label <- names(cand)[which.min(fdrs)]
    headline <- cand[[headline_label]]
  }
  dbeta <- if (!is.null(headline)) headline$dbeta else NA_real_
  fdr <- if (!is.null(headline)) headline$fdr_bacon else NA_real_
  direction <- if (!is.null(headline) && !is.na(dbeta)) {
    if (dbeta > 0) "Hypermethylated" else if (dbeta < 0) "Hypomethylated" else "No change"
  } else NA_character_
  sig <- if (!is.na(fdr)) (fdr <= 0.05) else NA
  gene_relationship <- if (!is.na(r$gene_names) && nzchar(r$gene_names)) paste(unique(trimws(strsplit(r$gene_names, ";")[[1]])), collapse = ", ") else NA_character_

  pairs <- list(
    "Gene" = d$primary_gene,
    "Chromosome" = r$chr,
    "Genomic position" = if (!is.na(r$pos)) sprintf("%s:%s", r$chr, format(r$pos, big.mark = ",")) else NA_character_,
    "Genome build" = "GRCh37 / hg19",
    "Methylation direction" = direction,
    "Delta Beta" = if (!is.na(dbeta)) sprintf("%.4f (%s stratum, SVA/bacon pipeline)", dbeta, headline_label) else NA_character_,
    "P-value" = if (!is.null(headline) && !is.na(headline$p_bacon)) signif(headline$p_bacon, 4) else NA_real_,
    "Adjusted P-value (FDR)" = if (!is.na(fdr)) signif(fdr, 4) else NA_real_,
    "Effect size (|Delta Beta|)" = if (!is.na(dbeta)) round(abs(dbeta), 4) else NA_real_,
    "CpG context" = bc_island_context_label(r$island_relation),
    "Genomic feature" = bc_feature_label(r$gene_group),
    "Gene relationship" = gene_relationship,
    "Biomarker status" = if (is.na(sig)) NA_character_ else if (sig) "Significant (FDR <= 0.05)" else "Not significant"
  )
  div(class = "card",
      div(class = "card-title", icon("id-card"), sprintf("Biomarker: %s", d$cpg)),
      p(class = "submodule-desc", "Source: SVA/bacon-corrected sex-stratified DMP pipeline (preloaded) + Illumina manifest / ChAMPdata annotation (hg19)."),
      bc_kv_table(pairs)
  )
}

bc_section_genomic_context <- function(d) {
  r <- d$resolved; isl <- d$island
  dist_txt <- if (!is.null(isl) && !is.na(r$pos)) {
    if (r$pos >= isl$start && r$pos <= isl$end) "0 (CpG lies within the island)"
    else if (r$pos < isl$start) sprintf("%s bp upstream of island start", format(isl$start - r$pos, big.mark = ","))
    else sprintf("%s bp downstream of island end", format(r$pos - isl$end, big.mark = ","))
  } else NA_character_

  pairs <- list(
    "CpG ID" = d$cpg, "Chromosome" = r$chr,
    "Exact genomic coordinate" = if (!is.na(r$pos)) format(r$pos, big.mark = ",") else NA_character_,
    "Genome assembly / build" = "GRCh37 / hg19", "Strand" = r$strand,
    "CpG island name / coordinates" = if (!is.null(isl)) sprintf("%s:%s-%s", isl$chr, format(isl$start, big.mark = ","), format(isl$end, big.mark = ",")) else NA_character_,
    "CpG island length (bp)" = if (!is.null(isl)) format(isl$length, big.mark = ",") else NA_character_,
    "Distance from CpG to island boundary" = dist_txt,
    "CpG island relation" = bc_island_context_label(r$island_relation),
    "Promoter / TSS relation" = bc_feature_label(r$gene_group),
    "Gene body relationship" = bc_feature_label(r$gene_group),
    "Enhancer relationship" = bc_enhancer_flag(r, d$array_type),
    "Regulatory feature" = if (!is.na(r$regulatory_group) && nzchar(r$regulatory_group))
      sprintf("%s%s", r$regulatory_group, if (!is.na(r$regulatory_name) && nzchar(r$regulatory_name)) sprintf(" (%s)", r$regulatory_name) else "")
      else "No annotated regulatory feature",
    "DNase hypersensitivity (DHS)" = bc_dhs_flag(r, d$array_type),
    "Nearby / overlapping gene(s)" = if (!is.na(r$gene_names) && nzchar(r$gene_names)) paste(unique(trimws(strsplit(r$gene_names, ";")[[1]])), collapse = ", ") else "None annotated"
  )
  div(class = "card",
      div(class = "card-title", icon("map-location-dot"), "CpG Genomic Context"),
      p(class = "submodule-desc", sprintf("Source: %s manifest annotation (Illumina/UCSC, hg19).", d$array_type)),
      bc_kv_table(pairs)
  )
}

bc_primary_gene_detail <- function(d) {
  gs <- d$gene_struct
  if (!isTRUE(gs$ok)) return(div(class = "empty-note", icon("circle-info"), gs$reason %||% "Gene structure could not be resolved."))
  tss <- if (identical(gs$strand, "-")) gs$end else gs$start
  pos <- d$resolved$pos
  dist <- if (!is.na(pos) && !is.na(gs$start) && !is.na(gs$end)) {
    if (pos >= gs$start && pos <= gs$end) "0 (CpG falls within the gene body)" else format(abs(pos - tss), big.mark = ",")
  } else NA_character_
  pairs <- list(
    "Gene symbol" = gs$symbol, "Gene name" = gs$genename, "NCBI Gene ID (Entrez)" = gs$entrez, "Ensembl Gene ID" = gs$ensembl,
    "Chromosome" = gs$chr,
    "Genomic location" = if (!is.na(gs$start) && !is.na(gs$end)) sprintf("%s:%s-%s", gs$chr, format(gs$start, big.mark = ","), format(gs$end, big.mark = ",")) else NA_character_,
    "Strand / transcription direction" = if (identical(gs$strand, "+")) "+ (forward)" else if (identical(gs$strand, "-")) "- (reverse)" else gs$strand,
    "Distance from CpG to TSS (bp)" = dist,
    "Top GO Biological Process terms" = if (!is.null(d$go_terms) && length(d$go_terms) > 0) paste(d$go_terms, collapse = "; ") else "Not available"
  )
  tagList(
    bc_kv_table(pairs),
    p(style = "font-size:11px; color:var(--color-ink-muted); margin-top:6px;",
      "Source: org.Hs.eg.db (NCBI/Ensembl mapping, GO annotation), TxDb.Hsapiens.UCSC.hg19.knownGene (gene structure), GO.db (GO term names).")
  )
}

bc_section_genes <- function(d) {
  rel_tab <- bc_gene_relationship_table(d$resolved)
  div(class = "card",
      div(class = "card-title", icon("dna"), "Associated Gene(s)"),
      p(class = "submodule-desc", "Source: Illumina manifest UCSC RefGene annotation (hg19) for the relationship table; NCBI/Ensembl/GO detail below is for the primary associated gene only."),
      if (!is.null(rel_tab) && nrow(rel_tab) > 0) DT::datatable(rel_tab, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
      else div(class = "empty-note", icon("circle-info"), "No RefGene-annotated gene overlaps this CpG (intergenic)."),
      tags$hr(), tags$b("Primary gene detail"),
      bc_primary_gene_detail(d)
  )
}

## ---- External evidence sections (rendered only after the "Look Up External
## Evidence" button; `ext` is a bc_external_evidence() result, or NULL before
## it's been fetched) --------------------------------------------------------

bc_ext_not_fetched <- function(title, icon_name) {
  div(class = "card", div(class = "card-title", icon(icon_name), title),
      div(class = "empty-note", icon("circle-info"), "Click \"Look Up External Evidence\" below to query this section."))
}

bc_section_disease_evidence <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Disease Evidence", "notes-medical"))
  cmb <- ext$combined_disease
  body <- if (is.null(cmb) || nrow(cmb) == 0) {
    div(class = "empty-note", icon("circle-info"), "No evidence identified in the connected databases (EWAS Catalog, EWAS Atlas) for this CpG.")
  } else {
    disp <- cmb[order(cmb$source, cmb$trait), , drop = FALSE]
    DT::datatable(disp, colnames = c("Source", "Trait", "Tissue", "Effect (beta/correlation)", "P-value", "PMID"),
                  rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  }
  div(class = "card", div(class = "card-title", icon("notes-medical"), "Disease Evidence"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog + EWAS Atlas (live query, this CpG)."),
      body)
}

bc_section_ra_evidence <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Rheumatoid Arthritis Evidence", "hand-dots"))
  ra <- ext$ra_rows
  found <- !is.null(ra) && nrow(ra) > 0
  banner <- if (found) div(class = "empty-note", style = "border-left-color:#0ca30c;", icon("circle-check"), tags$b("RA-associated biomarker"))
            else div(class = "empty-note", icon("circle-info"), tags$b("No RA-specific evidence identified in the connected databases."))
  body <- if (found) DT::datatable(ra[, c("source", "trait", "effect", "p", "pmid")], colnames = c("Source", "Trait", "Effect", "P-value", "PMID"),
                                    rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") else NULL
  div(class = "card", div(class = "card-title", icon("hand-dots"), "Rheumatoid Arthritis Evidence"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog + EWAS Atlas, filtered to rheumatoid-arthritis-related traits (live query, this CpG)."),
      banner, body)
}

bc_section_disease_comparison <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Disease Comparison", "chart-simple"))
  dc <- ext$disease_counts
  if (is.null(dc) || nrow(dc) == 0) {
    return(div(class = "card", div(class = "card-title", icon("chart-simple"), "Disease Comparison"),
               div(class = "empty-note", icon("circle-info"), "No disease/trait associations were found to compare.")))
  }
  dc <- dc[order(-dc$n), , drop = FALSE]
  dc$category <- factor(dc$category, levels = rev(dc$category))
  p <- ggplot(dc, aes(x = n, y = category)) + geom_col(fill = ARTHOMIX_COLORS$blue) +
    labs(x = "Real associations found (EWAS Catalog + EWAS Atlas)", y = NULL) + theme_arthomix()
  div(class = "card", div(class = "card-title", icon("chart-simple"), "Disease Comparison"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog + EWAS Atlas trait associations for this CpG, bucketed by disease category. Bar length is a real association count, not a synthetic score."),
      renderPlot_static(p))
}

bc_section_tissue_evidence <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Tissue / Cell-Type Evidence", "microscope"))
  tc <- ext$tissue_counts
  body <- if (is.null(tc) || nrow(tc) == 0) div(class = "empty-note", icon("circle-info"), "No tissue information was returned by the connected databases for this CpG.")
          else DT::datatable(tc[order(-tc$Freq), , drop = FALSE], colnames = c("Tissue", "Studies"), rownames = FALSE,
                              options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("microscope"), "Tissue / Cell-Type Evidence"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog (tissue field, live query, this CpG)."),
      body)
}

bc_section_pathways <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Biological Pathways", "diagram-project"))
  kegg <- ext$kegg; reactome <- ext$reactome
  kegg_body <- if (isTRUE(kegg$ok) && nrow(kegg$pathways) > 0) DT::datatable(kegg$pathways, colnames = c("KEGG ID", "Pathway name"), rownames = FALSE,
                                                                              options = list(pageLength = 5, scrollX = TRUE), class = "stripe hover compact")
               else div(class = "empty-note", icon("circle-info"), if (isTRUE(kegg$ok)) "No KEGG pathway annotation found for this gene." else (kegg$reason %||% "KEGG lookup unavailable."))
  reactome_body <- if (isTRUE(reactome$ok) && nrow(reactome$pathways) > 0) DT::datatable(reactome$pathways, colnames = c("Reactome Stable ID", "Pathway name"), rownames = FALSE,
                                                                                          options = list(pageLength = 5, scrollX = TRUE), class = "stripe hover compact")
                   else div(class = "empty-note", icon("circle-info"), if (isTRUE(reactome$ok)) "No Reactome pathway annotation found for this gene." else (reactome$reason %||% "Reactome lookup unavailable."))
  div(class = "card", div(class = "card-title", icon("diagram-project"), "Biological Pathways"),
      p(class = "submodule-desc", "Database annotation (membership lookup for the associated gene) - not a computed enrichment test, which needs a gene set and a background to test against."),
      tags$b("KEGG"), kegg_body, tags$hr(), tags$b("Reactome"), reactome_body,
      div(class = "empty-note", style = "margin-top:8px;", icon("circle-info"), "Gene Ontology terms for this gene are already shown under \"Associated Gene(s)\" above.")
  )
}

bc_section_kegg_ra_pathway <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("KEGG RA Pathway (hsa05323)", "network-wired"))
  kr <- ext$kegg_ra
  if (!isTRUE(kr$ok)) {
    return(div(class = "card", div(class = "card-title", icon("network-wired"), "KEGG RA Pathway (hsa05323)"),
               div(class = "empty-note", icon("triangle-exclamation"), kr$reason %||% "KEGG RA pathway lookup unavailable.")))
  }
  banner <- if (isTRUE(kr$in_pathway)) div(class = "empty-note", style = "border-left-color:#0ca30c;", icon("circle-check"),
                                            tags$b(sprintf("%s is a member gene of the KEGG Rheumatoid Arthritis pathway (hsa05323).", ext$gene_symbol %||% "This gene")))
             else div(class = "empty-note", icon("circle-info"), tags$b(sprintf("%s is not listed as a member gene of the KEGG Rheumatoid Arthritis pathway (hsa05323).", ext$gene_symbol %||% "This gene")))
  div(class = "card", div(class = "card-title", icon("network-wired"), "KEGG RA Pathway (hsa05323)"),
      p(class = "submodule-desc", "Source: KEGG (live query - pathway hsa05323 gene list)."),
      banner,
      p(style = "margin-top:8px;", tags$b(kr$pathway_name)),
      p(style = "font-size:12.5px; color:var(--color-ink-secondary);", kr$pathway_description),
      tags$details(tags$summary(sprintf("All %d genes on this pathway map", nrow(kr$genes))),
                   DT::datatable(kr$genes, colnames = c("NCBI Gene ID", "Symbol"), rownames = FALSE, filter = "top",
                                 options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact"))
  )
}

bc_section_publications <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Supporting Publications", "book"))
  pubs <- ext$publications
  body <- if (!isTRUE(pubs$ok)) div(class = "empty-note", icon("triangle-exclamation"), pubs$reason %||% "PubMed lookup unavailable.")
          else if (is.null(pubs$table) || nrow(pubs$table) == 0) div(class = "empty-note", icon("circle-info"), "No publication PMIDs were returned by the connected databases for this CpG.")
          else tags$ul(lapply(seq_len(nrow(pubs$table)), function(i) {
            r <- pubs$table[i, ]
            tags$li(sprintf("%s (%s). %s. ", r$title %||% "Untitled", r$year %||% "n.d.", r$journal %||% ""),
                    tags$a(href = sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", r$pmid), target = "_blank", sprintf("PMID: %s", r$pmid)))
          }))
  div(class = "card", div(class = "card-title", icon("book"), "Supporting Publications"),
      p(class = "submodule-desc", "Source: PMIDs from EWAS Catalog / EWAS Atlas results, resolved via NCBI PubMed (live query)."),
      body)
}

bc_section_methbank <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("MethBank", "link"))
  div(class = "card", div(class = "card-title", icon("link"), "MethBank"),
      div(class = "empty-note", icon("circle-info"),
          "MethBank has no public API in this deployment - it is referenced here as an outbound link only, not a live query.",
          tags$br(), tags$a(href = ext$methbank_link, target = "_blank", "Search MethBank manually")))
}

## Small ggplot-in-a-card helper for sections built outside a reactive
## render* binding (used by both the on-screen card and the static report).
renderPlot_static <- function(p) {
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggplot2::ggsave(tf, plot = p, width = 7, height = 3, dpi = 110, bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok || !requireNamespace("base64enc", quietly = TRUE)) return(div(class = "empty-note", "Plot unavailable."))
  uri <- base64enc::dataURI(file = tf, mime = "image/png")
  unlink(tf)
  tags$img(src = uri, style = "max-width:100%; height:auto;")
}

bc_section_your_dataset <- function(d, plot_tag) {
  pipeline_rows <- list()
  if (!is.null(d$dmp_female)) pipeline_rows[["Female (pipeline)"]] <- d$dmp_female
  if (!is.null(d$dmp_male)) pipeline_rows[["Male (pipeline)"]] <- d$dmp_male
  pipeline_table <- NULL
  if (length(pipeline_rows) > 0) {
    pipeline_table <- do.call(rbind, lapply(names(pipeline_rows), function(nm) {
      r <- pipeline_rows[[nm]]
      data.frame(Stratum = nm, `Delta Beta` = round(r$dbeta, 4), `P-value` = signif(r$p_bacon, 4), FDR = signif(r$fdr_bacon, 4),
                 Direction = if (is.na(r$dbeta)) NA_character_ else if (r$dbeta > 0) "Hyper" else if (r$dbeta < 0) "Hypo" else "None",
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }

  live_block <- NULL
  if (isTRUE(d$live$ok)) {
    ov <- d$live$overall
    if (isTRUE(ov$ok)) {
      live_pairs <- list(
        "Case group" = sprintf("%s (n=%d)", ov$case_label, ov$n_case), "Control group" = sprintf("%s (n=%d)", ov$control_label, ov$n_control),
        "Case mean beta" = round(ov$mean_case, 4), "Control mean beta" = round(ov$mean_control, 4),
        "Delta Beta" = round(ov$dbeta, 4), "P-value (Welch t-test)" = signif(ov$p_value, 4),
        "95% CI (case - control)" = if (!anyNA(ov$ci)) sprintf("[%.4f, %.4f]", ov$ci[1], ov$ci[2]) else NA_character_,
        "Direction" = ov$direction, "Missing values excluded" = sum(is.na(d$live$beta_row))
      )
      live_block <- tagList(
        p(style = "font-size:12px; color:var(--color-ink-secondary); margin-top:10px;",
          "Computed directly from the loaded beta matrix (Welch t-test, case vs control) - not SVA/bacon-adjusted."),
        bc_kv_table(live_pairs)
      )
    } else {
      live_block <- div(class = "empty-note", icon("triangle-exclamation"), ov$reason %||% "Could not compute live dataset statistics for this CpG.")
    }
  } else {
    live_block <- div(class = "empty-note", icon("circle-info"), d$live$reason %||% "No live/uploaded beta matrix and sample sheet are currently loaded.")
  }

  div(class = "card",
      div(class = "card-title", icon("flask"), "Your Dataset"),
      p(class = "submodule-desc", "Source: this app's own preloaded/uploaded methylation dataset - never an external database."),
      if (!is.null(pipeline_table)) tagList(
        p(style = "font-size:12px; color:var(--color-ink-secondary);", "Pipeline-reported statistics (SVA/bacon-corrected sex-stratified DMP pipeline):"),
        DT::datatable(pipeline_table, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
      ) else div(class = "empty-note", icon("circle-info"), "This CpG was not found in the preloaded sex-stratified DMP tables (not tested on the array, or preloaded pipeline data is unavailable)."),
      live_block,
      tags$div(style = "margin-top:10px;", plot_tag)
  )
}

bc_section_sex_specific <- function(d, plot_tag) {
  mkrow <- function(label, r) {
    if (is.null(r)) return(data.frame(Analysis = label, `Delta Beta` = NA_real_, `P-value` = NA_real_, FDR = NA_real_, Direction = NA_character_, Significant = NA_character_, check.names = FALSE, stringsAsFactors = FALSE))
    data.frame(Analysis = label, `Delta Beta` = round(r$dbeta, 4), `P-value` = signif(r$p_bacon, 4), FDR = signif(r$fdr_bacon, 4),
               Direction = if (is.na(r$dbeta)) NA_character_ else if (r$dbeta > 0) "Hyper" else if (r$dbeta < 0) "Hypo" else "None",
               Significant = if (is.na(r$fdr_bacon)) NA_character_ else if (r$fdr_bacon <= 0.05) "Yes" else "No",
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  pipeline_df <- rbind(mkrow("Female (pipeline)", d$dmp_female), mkrow("Male (pipeline)", d$dmp_male))

  ov <- if (isTRUE(d$live$ok)) d$live$overall else NULL
  overall_row <- if (!is.null(ov) && isTRUE(ov$ok)) {
    data.frame(Analysis = "Overall (your data, computed)", `Delta Beta` = round(ov$dbeta, 4), `P-value` = signif(ov$p_value, 4), FDR = NA_real_,
               Direction = ov$direction, Significant = if (!is.na(ov$p_value) && ov$p_value <= 0.05) "Yes (uncorrected p)" else "No",
               check.names = FALSE, stringsAsFactors = FALSE)
  } else data.frame(Analysis = "Overall (your data, computed)", `Delta Beta` = NA_real_, `P-value` = NA_real_, FDR = NA_real_, Direction = NA_character_, Significant = NA_character_, check.names = FALSE, stringsAsFactors = FALSE)

  live_rows <- NULL
  if (isTRUE(d$live$ok) && !is.null(d$live$by_sex)) {
    live_rows <- do.call(rbind, lapply(names(d$live$by_sex), function(s) {
      r <- d$live$by_sex[[s]]
      if (!isTRUE(r$ok)) return(data.frame(Analysis = sprintf("%s (your data, computed)", s), `Delta Beta` = NA_real_, `P-value` = NA_real_, FDR = NA_real_, Direction = NA_character_, Significant = NA_character_, check.names = FALSE, stringsAsFactors = FALSE))
      data.frame(Analysis = sprintf("%s (your data, computed)", s), `Delta Beta` = round(r$dbeta, 4), `P-value` = signif(r$p_value, 4), FDR = NA_real_,
                 Direction = r$direction, Significant = if (!is.na(r$p_value) && r$p_value <= 0.05) "Yes (uncorrected p)" else "No",
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }
  all_df <- rbind(pipeline_df, overall_row, if (!is.null(live_rows)) live_rows)

  f_sig <- !is.null(d$dmp_female) && !is.na(d$dmp_female$fdr_bacon) && d$dmp_female$fdr_bacon <= 0.05
  m_sig <- !is.null(d$dmp_male) && !is.na(d$dmp_male$fdr_bacon) && d$dmp_male$fdr_bacon <= 0.05
  verdict <- if (f_sig && m_sig) "Significant in both sexes (pipeline, FDR <= 0.05)."
             else if (f_sig && !m_sig) "Female-specific signal (pipeline, FDR <= 0.05 in females only)."
             else if (m_sig && !f_sig) "Male-specific signal (pipeline, FDR <= 0.05 in males only)."
             else if (!is.null(d$dmp_female) && !is.null(d$dmp_male)) "Not significant in either sex-stratified pipeline analysis (FDR <= 0.05)."
             else "Sex-stratified pipeline results are not both available for this CpG."

  div(class = "card",
      div(class = "card-title", icon("venus-mars"), "Sex-Specific Evidence"),
      p(class = "submodule-desc", "Source: sex-stratified SVA/bacon DMP pipeline (Female/Male rows); \"your data, computed\" rows are directly computed from the loaded beta matrix, not statistically adjusted - do not compare their raw p-values directly against the pipeline's FDR."),
      div(class = "empty-note", icon("circle-info"), tags$b(verdict)),
      DT::datatable(all_df, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact"),
      tags$div(style = "margin-top:10px;", plot_tag)
  )
}

bc_section_interpretation <- function(d, ext = NULL) {
  gene_label <- d$primary_gene %||% "Unknown gene"
  region_label <- paste0(bc_feature_label(d$resolved$gene_group), if (!is.null(d$island)) " / CpG island" else "")
  cand <- list(Female = d$dmp_female, Male = d$dmp_male); cand <- cand[!vapply(cand, is.null, logical(1))]
  dir_label <- "Not available"
  if (length(cand) > 0) {
    fdrs <- vapply(cand, function(x) if (is.na(x$fdr_bacon)) Inf else x$fdr_bacon, numeric(1))
    best <- cand[[names(cand)[which.min(fdrs)]]]
    if (!is.na(best$dbeta)) dir_label <- if (best$dbeta > 0) "Hypermethylated" else if (best$dbeta < 0) "Hypomethylated" else "No change"
  }
  pathway_label <- "Not yet connected"; pathway_done <- FALSE
  disease_label <- "Not yet connected"; disease_done <- FALSE
  ra_label <- "Not yet connected"; ra_done <- FALSE
  if (!is.null(ext)) {
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    pathway_label <- if (n_kegg + n_reactome > 0) sprintf("%d pathway(s) found", n_kegg + n_reactome) else "None found"
    pathway_done <- (n_kegg + n_reactome) > 0
    n_disease <- if (!is.null(ext$combined_disease)) nrow(ext$combined_disease) else 0
    disease_label <- if (n_disease > 0) sprintf("%d association(s) found", n_disease) else "None found"
    disease_done <- n_disease > 0
    ra_done <- !is.null(ext$ra_rows) && nrow(ext$ra_rows) > 0
    ra_label <- if (ra_done) "RA-associated" else "No RA evidence"
  }
  nodes <- list(
    list(label = "CpG", value = d$cpg, done = TRUE),
    list(label = "Region", value = region_label, done = TRUE),
    list(label = "Gene", value = gene_label, done = !is.na(d$primary_gene %||% NA)),
    list(label = "Methylation", value = dir_label, done = !identical(dir_label, "Not available")),
    list(label = "Pathway", value = pathway_label, done = pathway_done),
    list(label = "Disease", value = disease_label, done = disease_done),
    list(label = "RA Evidence", value = ra_label, done = ra_done)
  )
  div(class = "card",
      div(class = "card-title", icon("diagram-project"), "Functional Interpretation"),
      p(class = "submodule-desc", "Cautious, evidence-only chain - each arrow means \"associated with\" / \"located within\", never a causal claim."),
      div(style = "display:flex; flex-wrap:wrap; align-items:center; gap:6px;",
          lapply(seq_along(nodes), function(i) {
            n <- nodes[[i]]
            tagList(
              div(style = sprintf("border:1px solid %s; border-radius:8px; padding:8px 12px; min-width:110px; text-align:center; background:%s;",
                                   if (n$done) ARTHOMIX_COLORS$blue else "#cccccc", if (n$done) "#EAF3FB" else "#F5F5F5"),
                  tags$div(style = "font-size:11px; text-transform:uppercase; color:var(--color-ink-secondary);", n$label),
                  tags$div(style = "font-weight:600; font-size:13px;", n$value)
              ),
              if (i < length(nodes)) icon("arrow-right", style = "color:#999;")
            )
          })
      )
  )
}

bc_section_evidence_summary <- function(d, ext = NULL) {
  user_sig <- FALSE
  if (!is.null(d$dmp_female) && !is.na(d$dmp_female$fdr_bacon) && d$dmp_female$fdr_bacon <= 0.05) user_sig <- TRUE
  if (!is.null(d$dmp_male) && !is.na(d$dmp_male$fdr_bacon) && d$dmp_male$fdr_bacon <= 0.05) user_sig <- TRUE
  if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok) && !is.na(d$live$overall$p_value) && d$live$overall$p_value <= 0.05) user_sig <- TRUE
  gene_assoc <- !is.na(d$primary_gene %||% NA)
  island_reg <- !is.null(d$island) || (!is.na(d$resolved$regulatory_group) && nzchar(d$resolved$regulatory_group))
  sex_specific <- !is.null(d$dmp_female) && !is.null(d$dmp_male)

  ewas_repl <- "Not yet checked"; ra_evid <- "Not yet checked"; pathway_assoc <- "Not yet checked"
  if (!is.null(ext)) {
    n_disease <- if (!is.null(ext$combined_disease)) nrow(ext$combined_disease) else 0
    ewas_repl <- if (n_disease > 0) "Yes" else "No"
    ra_evid <- if (!is.null(ext$ra_rows) && nrow(ext$ra_rows) > 0) "Yes" else "No"
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    pathway_assoc <- if ((n_kegg + n_reactome) > 0) "Yes" else "No"
  }

  rows <- list(
    c("User dataset significant (FDR/P <= 0.05)", if (user_sig) "Yes" else "No"),
    c("Gene association", if (gene_assoc) "Yes" else "No"),
    c("CpG island / regulatory annotation", if (island_reg) "Yes" else "No"),
    c("Sex-specific evidence (both sexes tested)", if (sex_specific) "Yes" else "No"),
    c("Independent EWAS replication (EWAS Catalog / EWAS Atlas)", ewas_repl),
    c("RA evidence", ra_evid),
    c("Pathway association (KEGG / Reactome)", pathway_assoc),
    c("Methylation-expression evidence", "Not yet integrated")
  )
  df <- data.frame(Evidence = vapply(rows, `[`, character(1), 1), Status = vapply(rows, `[`, character(1), 2), stringsAsFactors = FALSE)
  div(class = "card",
      div(class = "card-title", icon("clipboard-check"), "Evidence Summary"),
      p(class = "submodule-desc", "Each row is a transparent, individually-computed check - not a combined score. Rows marked \"Not yet integrated\" are planned for a later phase, not evaluated as absent."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
  )
}

bc_section_sources <- function(ext = NULL) {
  used <- c(
    "Illumina manifest annotation (IlluminaHumanMethylation450kanno.ilmn12.hg19 / IlluminaHumanMethylationEPICanno.ilm10b4.hg19) - UCSC hg19-based probe location, CpG island, RefGene, and regulatory-feature tracks",
    "ChAMPdata::probe.features - UCSC hg19-based gene/feature/CpG-island annotation",
    "TxDb.Hsapiens.UCSC.hg19.knownGene - UCSC knownGene transcript/exon structure (hg19)",
    "org.Hs.eg.db - NCBI Gene ID, Ensembl Gene ID, gene name, GO term mapping",
    "GO.db - Gene Ontology term names",
    "UCSC cytoBandIdeo (hg19) - chromosome ideogram/cytoband track, bundled from hgdownload.soe.ucsc.edu",
    "This app's own preloaded sex-stratified SVA/bacon-corrected DMP pipeline results, and/or your uploaded dataset"
  )
  external_available <- c(
    "MRC-IEU EWAS Catalog - live query, https://www.ewascatalog.org/api/",
    "EWAS Atlas - live query, https://ngdc.cncb.ac.cn/ewas/rest/",
    "KEGG - live query via KEGGREST, https://rest.kegg.jp/",
    "Reactome - live query, https://reactome.org/ContentService/",
    "PubMed (NCBI E-utilities) - live query, https://eutils.ncbi.nlm.nih.gov/"
  )
  external_label <- if (is.null(ext)) "Available (click \"Look Up External Evidence\" to query):" else "Used in this Biomarker Card (external, live-queried):"
  div(class = "card",
      div(class = "card-title", icon("database"), "Database Sources"),
      tags$b("Used in this Biomarker Card (local/offline):"), tags$ul(lapply(used, tags$li)),
      tags$b(external_label), tags$ul(lapply(external_available, tags$li)),
      tags$b("Referenced, no public API (outbound link only):"), tags$ul(tags$li("MethBank"))
  )
}

bc_report_css <- function() {
  "body{font-family:-apple-system,Helvetica,Arial,sans-serif; max-width:900px; margin:24px auto; color:#222;}
   .card{border:1px solid #ddd; border-radius:10px; padding:14px 18px; margin-bottom:16px;}
   .card-title{font-weight:700; font-size:15px; margin-bottom:8px;}
   .submodule-desc{color:#666; font-size:12.5px;}
   .empty-note{background:#f6f6f6; border-left:3px solid #999; padding:8px 12px; border-radius:4px; font-size:13px;}
   table{border-collapse:collapse; width:100%;} td,th{border:1px solid #eee; padding:4px 8px; font-size:13px; text-align:left;}"
}

bc_build_report_tags <- function(d, cytoband_df, ext = NULL) {
  region_plot <- tryCatch(bc_plot_region(d$resolved$chr, d$resolved$pos, d$island, d$gene_struct, 5000), error = function(e) NULL)
  ideogram_plot <- tryCatch(bc_plot_ideogram(d$resolved$chr, d$resolved$pos, cytoband_df), error = function(e) NULL)
  dist_plot <- if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok))
    tryCatch(bc_plot_methylation_dist(data.frame(beta = d$live$beta_row, group = d$live$group_vec)), error = function(e) NULL) else NULL
  sex_plot <- if (isTRUE(d$live$ok) && !is.null(d$live$sex_vec))
    tryCatch(bc_plot_methylation_dist(data.frame(beta = d$live$beta_row, group = d$live$group_vec, sex = d$live$sex_vec), facet_sex = TRUE), error = function(e) NULL) else NULL

  img_tag <- function(p) {
    uri <- bc_ggsave_datauri(p)
    if (is.null(uri)) div(class = "empty-note", "Plot unavailable.") else tags$img(src = uri, style = "max-width:100%; height:auto;")
  }

  tagList(
    tags$h2(sprintf("Biomarker Card: %s", d$cpg)),
    tags$p(sprintf("Genome build: GRCh37 / hg19. Array: %s.", d$array_type)),
    bc_section_summary(d), bc_section_genomic_context(d),
    div(class = "card", div(class = "card-title", "Local Region View"), img_tag(region_plot)),
    div(class = "card", div(class = "card-title", "Chromosome View"), img_tag(ideogram_plot)),
    bc_section_genes(d),
    bc_section_disease_evidence(ext), bc_section_ra_evidence(ext), bc_section_disease_comparison(ext),
    bc_section_tissue_evidence(ext), bc_section_pathways(ext), bc_section_kegg_ra_pathway(ext),
    bc_section_publications(ext), bc_section_methbank(ext),
    bc_section_your_dataset(d, img_tag(dist_plot)),
    bc_section_sex_specific(d, img_tag(sex_plot)),
    bc_section_interpretation(d, ext), bc_section_evidence_summary(d, ext), bc_section_sources(ext),
    tags$p(style = "color:#888; font-size:12px; margin-top:16px;",
           if (is.null(ext)) "External evidence (EWAS Catalog, EWAS Atlas, KEGG, Reactome, PubMed) was not looked up before this report was generated - go back to the card and click \"Look Up External Evidence\" first if you want it included."
           else "Not yet available in this build: Methylation<->Expression correlation evidence, and a rendered KEGG pathway diagram image (a real membership/gene-list lookup is shown instead).")
  )
}

## ---- Config / UI -----------------------------------------------------

mod_methyl_biomarkercard_config <- list(
  id = "biomarkercard", title = "Biomarker Card", icon = "id-card", group = "Biomarker modeling",
  description = "A single-CpG epigenomic profile: genomic location, CpG island/gene context, chromosome view, and evidence from your own dataset (including sex-specific results)."
)

mod_methyl_biomarkercard_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("bc_subtabs"), type = "tabs",
      tabPanel("Select Biomarker", br(), uiOutput(ns("select_ui"))),
      tabPanel("Biomarker Card", br(), uiOutput(ns("bc_card_ui")))
    )
  )
}

## ---- Server ------------------------------------------------------------

mod_methyl_biomarkercard_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select_ui <- renderUI({
      tagList(
        div(class = "card",
            div(class = "card-title", icon("magnifying-glass"), "Data source"),
            radioButtons(ns("bc_source"), NULL, inline = TRUE,
                         choiceNames = list(
                           tagList(icon("database"), " Preloaded methylation results (GSE42861, 450K, sex-stratified)"),
                           tagList(icon("upload"), " Your loaded dataset (Dataset tab)")
                         ),
                         choiceValues = list("preloaded", "live"),
                         selected = if (METH_DATA_AVAILABLE) "preloaded" else "live"),
            if (!METH_DATA_AVAILABLE) div(class = "empty-note", icon("triangle-exclamation"), "The preloaded methylomics results folder is not available in this deployment - only your loaded dataset can be used."),
            conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("bc_source")),
                              p(class = "submodule-desc", "Genomic/gene lookups use the 450K array manifest; Delta-Beta/P/FDR lookups use the sex-stratified SVA/bacon-corrected DMP pipeline (Female and Male are always shown together in the generated card)."))
        ),
        div(class = "card",
            div(class = "card-title", icon("list-check"), "Find a biomarker"),
            radioButtons(ns("bc_search_mode"), NULL, inline = TRUE,
                         choices = c("Paste/type a CpG ID" = "cpg", "Search by gene symbol" = "gene", "Browse preloaded DMP results" = "browse",
                                     "Browse diagnostic model panel" = "panel", "Upload a biomarker list" = "upload")),
            conditionalPanel(condition = sprintf("input['%s'] == 'cpg'", ns("bc_search_mode")),
                              textInput(ns("bc_cpg_input"), "CpG ID", placeholder = "cg12277888")),
            conditionalPanel(condition = sprintf("input['%s'] == 'gene'", ns("bc_search_mode")),
                              fluidRow(
                                column(6, textInput(ns("bc_gene_input"), "Gene symbol", placeholder = "TRIM44")),
                                column(6, br(), actionButton(ns("bc_gene_search_btn"), "Search Gene", icon = icon("magnifying-glass"), class = "btn-sm"))
                              ),
                              uiOutput(ns("bc_gene_results_ui"))),
            conditionalPanel(condition = sprintf("input['%s'] == 'browse'", ns("bc_search_mode")),
                              if (!METH_DATA_AVAILABLE) div(class = "empty-note", icon("circle-info"), "Not available - the preloaded results folder is not present in this deployment.")
                              else tagList(
                                fluidRow(
                                  column(6, radioButtons(ns("bc_browse_sex"), "Sex stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female")),
                                  column(6, br(), actionButton(ns("bc_browse_load_btn"), "Load Top Results", icon = icon("play"), class = "btn-sm"))
                                ),
                                uiOutput(ns("bc_browse_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bc_search_mode")),
                              if (!METH_DATA_AVAILABLE) div(class = "empty-note", icon("circle-info"), "Not available - the preloaded results folder is not present in this deployment.")
                              else tagList(
                                p(class = "submodule-desc", "The Diagnostic Classifier's own \"potential biomarkers\": the script07 majority-vote CpG panel (in_LASSO/in_Boruta/in_SVMRFE, n_votes) - the same table that module's own \"Feature Source\" tab reads."),
                                fluidRow(
                                  column(6, radioButtons(ns("bc_panel_sex"), "Sex stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female")),
                                  column(6, br(), actionButton(ns("bc_panel_load_btn"), "Load Panel", icon = icon("play"), class = "btn-sm"))
                                ),
                                uiOutput(ns("bc_panel_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("bc_search_mode")),
                              p(class = "submodule-desc", "Upload a plain CpG ID list (one per line, or the first column of a CSV/TSV) or a Feature Selection RDS export (that module's own \"Save Model as RDS\" download) - auto-detected by file extension."),
                              fileInput(ns("bc_upload_file"), "Biomarker list (.csv, .txt, or .rds)", accept = c(".csv", ".txt", ".rds")),
                              actionButton(ns("bc_upload_load_btn"), "Load Uploaded List", icon = icon("play"), class = "btn-sm"),
                              uiOutput(ns("bc_upload_results_ui")))
        ),
        div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap;",
            uiOutput(ns("bc_selection_status_ui"), inline = TRUE),
            actionButton(ns("bc_generate_btn"), "Generate Biomarker Card", icon = icon("id-card"), class = "btn-primary btn-sm"))
      )
    })

    bc_picked_cpg <- reactiveVal(NULL)
    has_card <- reactiveVal(FALSE)
    observeEvent(input$bc_source, { bc_picked_cpg(NULL); has_card(FALSE) }, ignoreInit = TRUE)
    observeEvent(input$bc_search_mode, bc_picked_cpg(NULL), ignoreInit = TRUE)

    ## Explicit confirmation that a row/CpG is actually selected - without
    ## this, clicking a table row (which highlights it but doesn't otherwise
    ## change the page) is easy to mistake for "nothing happened", especially
    ## once the table has scrolled the Generate button out of view.
    output$bc_selection_status_ui <- renderUI({
      mode <- input$bc_search_mode %||% "cpg"
      cpg <- if (identical(mode, "cpg")) trimws(input$bc_cpg_input %||% "") else bc_picked_cpg()
      if (!is.null(cpg) && nzchar(cpg)) {
        tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("Selected: %s", cpg)), " - click Generate to build the card.")
      } else {
        tagList(icon("circle-info"), "Select a biomarker above (click a table row, or type a CpG ID) before clicking Generate.")
      }
    })

    ## ---- Gene search ----
    gene_search_result <- eventReactive(input$bc_gene_search_btn, {
      req(nzchar(trimws(input$bc_gene_input %||% "")))
      array_type <- if (identical(input$bc_source, "preloaded")) "450K" else (dataset$array_type %||% "450K")
      bc_search_gene(trimws(input$bc_gene_input), array_type)
    }, ignoreInit = TRUE)

    output$bc_gene_results_ui <- renderUI({
      r <- gene_search_result(); req(r)
      if (!isTRUE(r$ok)) return(div(class = "empty-note", icon("triangle-exclamation"), r$reason))
      if (nrow(r$df) == 0) return(div(class = "empty-note", icon("circle-info"), "No CpG probes were found annotated to this gene symbol on this array."))
      tagList(
        p(class = "submodule-desc", sprintf("%d CpG probe(s) found. Click a row to select it, then click \"Generate Biomarker Card\".", nrow(r$df))),
        DT::dataTableOutput(ns("bc_gene_table"))
      )
    })
    output$bc_gene_table <- DT::renderDataTable({
      r <- gene_search_result(); req(r); req(isTRUE(r$ok)); req(nrow(r$df) > 0)
      DT::datatable(r$df[, c("cpg", "chr", "pos", "gene_names", "gene_group", "island_relation")],
                    colnames = c("CpG", "Chr", "Position", "Gene(s)", "Region", "Island relation"),
                    rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bc_gene_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_gene_table_rows_selected, {
      r <- gene_search_result(); req(r)
      idx <- input$bc_gene_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(r$df$cpg[idx])
    })

    ## ---- Browse preloaded DMP results ----
    browse_table <- eventReactive(input$bc_browse_load_btn, {
      validate(need(METH_DATA_AVAILABLE, "The preloaded methylomics results folder is not available in this deployment."))
      df <- bc_dmp_table(input$bc_browse_sex %||% "female")
      validate(need(!is.null(df), "The preloaded DMP results table is not available in this deployment."))
      df <- df[order(df$fdr_bacon), , drop = FALSE]
      utils::head(df, 500)
    }, ignoreInit = TRUE)

    output$bc_browse_results_ui <- renderUI({
      req(input$bc_browse_load_btn)
      tagList(
        p(class = "submodule-desc", "Top 500 CpGs by FDR (SVA/bacon pipeline). Click a row to select it, then click \"Generate Biomarker Card\"."),
        DT::dataTableOutput(ns("bc_browse_table"))
      )
    })
    output$bc_browse_table <- DT::renderDataTable({
      df <- browse_table(); req(df)
      disp <- df[, c("cpg", "dbeta", "p_bacon", "fdr_bacon")]
      dt <- DT::datatable(disp, colnames = c("CpG", "Delta Beta", "P-value", "FDR"), rownames = FALSE, selection = "single",
                           options = list(pageLength = 10, scrollX = TRUE))
      DT::formatSignif(dt, columns = c("dbeta", "p_bacon", "fdr_bacon"), digits = 4)
    })
    outputOptions(output, "bc_browse_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_browse_table_rows_selected, {
      df <- browse_table(); req(df)
      idx <- input$bc_browse_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(df$cpg[idx])
    })

    ## ---- Browse diagnostic model panel ("potential biomarkers") ----
    panel_table <- eventReactive(input$bc_panel_load_btn, {
      validate(need(METH_DATA_AVAILABLE, "The preloaded methylomics results folder is not available in this deployment."))
      df <- load_default_diagnostic_ensemble_votes(input$bc_panel_sex %||% "female")
      validate(need(!is.null(df), "The preloaded diagnostic ensemble-vote panel is not available in this deployment."))
      df[order(-df$n_votes), , drop = FALSE]
    }, ignoreInit = TRUE)

    output$bc_panel_results_ui <- renderUI({
      req(input$bc_panel_load_btn)
      tagList(
        p(class = "submodule-desc", sprintf("%d CpG(s) in the majority-vote panel. Click a row to select it, then click \"Generate Biomarker Card\".", nrow(panel_table()))),
        DT::dataTableOutput(ns("bc_panel_table"))
      )
    })
    output$bc_panel_table <- DT::renderDataTable({
      df <- panel_table(); req(df)
      DT::datatable(df, colnames = c("CpG", "in LASSO", "in Boruta", "in SVM-RFE", "Votes"), rownames = FALSE, selection = "single",
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bc_panel_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_panel_table_rows_selected, {
      df <- panel_table(); req(df)
      idx <- input$bc_panel_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(df$cpg[idx])
    })

    ## ---- Upload a biomarker list (plain ID list, or a Feature Selection RDS
    ## export - same schema check mod_methyl_diagnostic.R's own "Feature
    ## Source: featureselection" upload path uses) ----
    upload_table <- eventReactive(input$bc_upload_load_btn, {
      validate(need(!is.null(input$bc_upload_file), "Upload a .csv, .txt, or .rds file first."))
      path <- input$bc_upload_file$datapath; name <- input$bc_upload_file$name
      if (grepl("\\.rds$", name, ignore.case = TRUE)) {
        obj <- tryCatch(readRDS(path), error = function(e) NULL)
        validate(need(!is.null(obj) && identical(obj$module, "mod_methyl_featureselection"),
                      "Upload a Feature Selection RDS export (from that module's own \"Save Model as RDS\" download), or a plain .csv/.txt CpG ID list."))
        ids <- as.character(obj$final_panel$cpg_ids %||% character(0))
        validate(need(length(ids) > 0, "This RDS export's final_panel$cpg_ids is empty."))
        df <- data.frame(cpg = ids, stringsAsFactors = FALSE)
        ann <- obj$final_panel$annotation
        if (!is.null(ann) && "cpg" %in% colnames(ann)) df <- merge(df, ann, by = "cpg", all.x = TRUE)
      } else if (grepl("\\.txt$", name, ignore.case = TRUE)) {
        pl <- methyl_parse_probe_list(path, name)
        validate(need(isTRUE(pl$ok), pl$error %||% "Could not parse the uploaded file."))
        df <- data.frame(cpg = pl$ids, stringsAsFactors = FALSE)
      } else {
        up <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(up) && nrow(up) > 0, "Could not parse the uploaded file as a delimited table (CSV/TSV)."))
        cpg_col <- intersect(c("cpg", "CpG", "probe", "probe_id", "feature", "ID"), colnames(up))[1]
        ids <- if (!is.na(cpg_col)) as.character(up[[cpg_col]]) else as.character(up[[1]])
        df <- data.frame(cpg = unique(ids[nzchar(ids)]), stringsAsFactors = FALSE)
      }
      validate(need(nrow(df) > 0, "No CpG IDs were found in the uploaded file."))
      df
    }, ignoreInit = TRUE)

    output$bc_upload_results_ui <- renderUI({
      req(input$bc_upload_load_btn)
      tagList(
        p(class = "submodule-desc", sprintf("%d CpG(s) loaded from the uploaded file. Click a row to select it, then click \"Generate Biomarker Card\".", nrow(upload_table()))),
        DT::dataTableOutput(ns("bc_upload_table"))
      )
    })
    output$bc_upload_table <- DT::renderDataTable({
      df <- upload_table(); req(df)
      DT::datatable(df, rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bc_upload_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_upload_table_rows_selected, {
      df <- upload_table(); req(df)
      idx <- input$bc_upload_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(df$cpg[idx])
    })

    ## ---- Generate ----
    ## CpG resolution is inlined directly into card_data() rather than split
    ## into its own eventReactive - calling one eventReactive from inside
    ## another eventReactive bound to the exact same button click event
    ## deadlocks (the inner one's internal "did my event fire" check never
    ## resolves when invoked from inside the outer one's isolate()d handler).
    observeEvent(input$bc_generate_btn, {
      has_card(TRUE)
      updateTabsetPanel(session, "bc_subtabs", selected = "Biomarker Card")
    }, ignoreInit = TRUE)

    card_data <- eventReactive(input$bc_generate_btn, {
      cpg <- switch(input$bc_search_mode,
        cpg = trimws(input$bc_cpg_input %||% ""),
        gene = bc_picked_cpg(),
        browse = bc_picked_cpg(),
        panel = bc_picked_cpg(),
        upload = bc_picked_cpg(),
        NULL
      )
      validate(need(!is.null(cpg) && nzchar(cpg),
                    "Select or enter a CpG ID before generating the Biomarker Card - type a CpG ID, or pick one from a gene search / the preloaded results browser / the diagnostic panel / your uploaded list."))
      source_mode <- isolate(input$bc_source)
      array_type <- if (identical(source_mode, "preloaded")) "450K" else (dataset$array_type %||% "450K")
      resolved <- bc_resolve_cpg(cpg, array_type)
      validate(need(isTRUE(resolved$found_in_manifest) || isTRUE(resolved$found_in_champ),
                    sprintf("\"%s\" could not be found in the %s manifest or the ChAMP reference annotation - check the CpG ID.", cpg, array_type)))

      primary_gene <- if (!is.na(resolved$champ_gene) && nzchar(resolved$champ_gene)) resolved$champ_gene
                       else if (!is.na(resolved$gene_names) && nzchar(resolved$gene_names)) trimws(strsplit(resolved$gene_names, ";")[[1]][1])
                       else NA_character_

      island <- bc_island_parse(resolved$island_name)
      gene_struct <- if (!is.na(primary_gene)) bc_gene_structure(primary_gene) else list(ok = FALSE, reason = "No associated gene to resolve.")
      go_terms <- if (isTRUE(gene_struct$ok)) bc_go_terms(gene_struct$entrez) else NULL

      dmp_female <- bc_dmp_lookup(cpg, "female")
      dmp_male <- bc_dmp_lookup(cpg, "male")
      dmr_female <- bc_dmr_lookup(resolved$chr, resolved$pos, "female")
      dmr_male <- bc_dmr_lookup(resolved$chr, resolved$pos, "male")

      live <- bc_dataset_evidence(cpg, dataset)

      list(cpg = cpg, source_mode = source_mode, array_type = array_type, resolved = resolved,
           primary_gene = primary_gene, island = island, gene_struct = gene_struct, go_terms = go_terms,
           dmp_female = dmp_female, dmp_male = dmp_male, dmr_female = dmr_female, dmr_male = dmr_male, live = live)
    }, ignoreInit = TRUE)

    observeEvent(card_data(), {
      d <- card_data()
      if (!is.null(results)) {
        results$biomarkercard <- list(
          cpg = d$cpg, gene = d$primary_gene %||% NA_character_, chr = d$resolved$chr, pos = d$resolved$pos,
          female_fdr = if (!is.null(d$dmp_female)) d$dmp_female$fdr_bacon else NA_real_,
          male_fdr = if (!is.null(d$dmp_male)) d$dmp_male$fdr_bacon else NA_real_
        )
      }
    }, ignoreInit = TRUE)

    ## External evidence is opt-in (own button/reactive, not part of
    ## card_data()) so a slow/down external service never blocks the fully
    ## local, already-verified core card above. Resets when a new biomarker
    ## is generated so stale evidence from a previous CpG never lingers.
    ## Single observeEvent, not a separate eventReactive+observer pair -
    ## calling an eventReactive bound to input$bc_ext_lookup_btn from inside
    ## an observeEvent bound to that exact same button click deadlocks (same
    ## same-event-nesting hazard as card_data()'s own note above).
    ext_data <- reactiveVal(NULL)
    observeEvent(card_data(), ext_data(NULL), ignoreInit = TRUE)
    observeEvent(input$bc_ext_lookup_btn, {
      d <- card_data(); req(d)
      ext_data(bc_external_evidence(d$cpg, d$primary_gene, d$gene_struct$entrez %||% NA_character_))
    }, ignoreInit = TRUE)

    output$bc_ext_sections_ui <- renderUI({
      tagList(
        bc_section_disease_evidence(ext_data()),
        bc_section_ra_evidence(ext_data()),
        bc_section_disease_comparison(ext_data()),
        bc_section_tissue_evidence(ext_data()),
        bc_section_pathways(ext_data()),
        bc_section_kegg_ra_pathway(ext_data()),
        bc_section_publications(ext_data()),
        bc_section_methbank(ext_data())
      )
    })

    output$bc_card_ui <- renderUI({
      if (!has_card()) return(div(class = "empty-note", icon("circle-info"), "Select a biomarker and click \"Generate Biomarker Card\" on the \"Select Biomarker\" tab."))
      d <- card_data(); req(d)
      tagList(
        bc_section_summary(d),
        bc_section_genomic_context(d),
        div(class = "card",
            div(class = "card-title", icon("magnifying-glass-chart"), "Local Region View"),
            p(class = "submodule-desc", "Source: Illumina manifest annotation (hg19) + TxDb.Hsapiens.UCSC.hg19.knownGene gene structure, ChAMPdata CpG island context."),
            fluidRow(
              column(4, numericInput(ns("bc_flank_kb"), "Window (kb, +/- around the CpG)", value = 5, min = 1, max = 200, step = 1)),
              column(4, br(), actionButton(ns("bc_zoom_tight_btn"), "Zoom to Biomarker", icon = icon("magnifying-glass-plus"), class = "btn-sm")),
              column(4, br(), actionButton(ns("bc_zoom_local_btn"), "View Local Region", icon = icon("expand"), class = "btn-sm"))
            ),
            withSpinner(plotly::plotlyOutput(ns("bc_region_plot"), height = "260px"), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("dna"), "Chromosome View"),
            p(class = "submodule-desc", sprintf("Source: UCSC cytoBandIdeo (hg19), chromosome %s.", d$resolved$chr %||% "?")),
            withSpinner(plotly::plotlyOutput(ns("bc_ideogram_plot"), height = "220px"), color = "#2563EB", type = 6)
        ),
        bc_section_genes(d),
        div(class = "card",
            div(class = "card-title", icon("globe"), "External Evidence"),
            p(class = "submodule-desc", "Live queries to MRC-IEU EWAS Catalog, EWAS Atlas, KEGG, and Reactome (plus PubMed for any publications they cite) - opt-in, since these are external network calls independent of everything above."),
            actionButton(ns("bc_ext_lookup_btn"), "Look Up External Evidence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm")
        ),
        withSpinner(uiOutput(ns("bc_ext_sections_ui")), color = "#2563EB", type = 6),
        bc_section_your_dataset(d, withSpinner(plotly::plotlyOutput(ns("bc_dist_plot"), height = "300px"), color = "#2563EB", type = 6)),
        bc_section_sex_specific(d, withSpinner(plotly::plotlyOutput(ns("bc_sex_dist_plot"), height = "320px"), color = "#2563EB", type = 6)),
        bc_section_interpretation(d, ext_data()),
        bc_section_evidence_summary(d, ext_data()),
        bc_section_sources(ext_data()),
        div(class = "card", div(class = "card-title", icon("download"), "Download"),
            p(class = "submodule-desc", "A self-contained HTML report covering every section above (including external evidence, if you've looked it up), plus a note listing which spec sections are still pending."),
            downloadButton(ns("bc_download_report"), "Download Biomarker Report (HTML)", class = "btn-primary btn-sm"))
      )
    })

    observeEvent(input$bc_zoom_tight_btn, updateNumericInput(session, "bc_flank_kb", value = 1))
    observeEvent(input$bc_zoom_local_btn, updateNumericInput(session, "bc_flank_kb", value = 10))

    output$bc_region_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      flank_bp <- (input$bc_flank_kb %||% 5) * 1000
      p <- bc_plot_region(d$resolved$chr, d$resolved$pos, d$island, d$gene_struct, flank_bp)
      plotly::ggplotly(p) %>% plotly::layout(showlegend = FALSE)
    })

    output$bc_ideogram_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      cb <- bc_cytoband_hg19()
      p <- bc_plot_ideogram(d$resolved$chr, d$resolved$pos, cb)
      plotly::ggplotly(p) %>% plotly::layout(showlegend = FALSE)
    })

    output$bc_dist_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      validate(need(isTRUE(d$live$ok) && isTRUE(d$live$overall$ok), "No live dataset evidence is available to plot for this CpG."))
      df <- data.frame(beta = d$live$beta_row, group = d$live$group_vec, stringsAsFactors = FALSE)
      plotly::ggplotly(bc_plot_methylation_dist(df))
    })

    output$bc_sex_dist_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      validate(need(isTRUE(d$live$ok) && !is.null(d$live$sex_vec), "No sex information is available in the loaded sample sheet to compare female vs male methylation."))
      df <- data.frame(beta = d$live$beta_row, group = d$live$group_vec, sex = d$live$sex_vec, stringsAsFactors = FALSE)
      plotly::ggplotly(bc_plot_methylation_dist(df, facet_sex = TRUE))
    })

    output$bc_download_report <- downloadHandler(
      filename = function() sprintf("biomarker_card_%s.html", card_data()$cpg),
      content = function(file) {
        d <- card_data(); req(d)
        cb <- bc_cytoband_hg19()
        body <- bc_build_report_tags(d, cb, ext_data())
        page <- tags$html(
          tags$head(tags$meta(charset = "utf-8"), tags$title(sprintf("Biomarker Card %s", d$cpg)), tags$style(bc_report_css())),
          tags$body(body)
        )
        htmltools::save_html(page, file = file)
      }
    )
  })
}
