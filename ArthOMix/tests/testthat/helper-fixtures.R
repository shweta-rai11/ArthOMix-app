## Deterministic, in-memory test fixtures shared across transcriptomics test
## files. Kept as generator functions (not static files) so every test that
## needs a "clean" matrix, or one perturbed in one specific way (a duplicate

fx_geo_eset <- function(n_probes = 30, n_samples = 10, seed = 1) {
  set.seed(seed)
  genes <- sprintf("PROBE%d", seq_len(n_probes))
  samples <- sprintf("GSM%04d", 1000 + seq_len(n_samples))
  mat <- matrix(runif(n_probes * n_samples, 3, 12), n_probes, n_samples,
                dimnames = list(genes, samples))
  pdat <- data.frame(
    title = sprintf("sample_%d", seq_len(n_samples)),
    geo_accession = samples,
    status = "Public on Jan 01 2020",
    characteristics_ch1 = "tissue: whole blood",
    `disease state:ch1` = rep(c("RA", "HC"), length.out = n_samples),
    `Sex:ch1` = rep(c("female", "male"), length.out = n_samples),
    check.names = FALSE, row.names = samples, stringsAsFactors = FALSE
  )
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_FIXTURE"
  eset
}

fx_expr_meta <- function(n_genes = 40, n_samples = 16, seed = 1) {
  set.seed(seed)
  genes <- sprintf("GENE%03d", seq_len(n_genes))
  samples <- sprintf("S%02d", seq_len(n_samples))
  expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5), n_genes, n_samples,
                 dimnames = list(genes, samples))
  meta <- data.frame(
    sample = samples,
    group = rep(c("HC", "RA"), length.out = n_samples),
    sex = rep(c("F", "M"), length.out = n_samples),
    batch = rep(c("batch1", "batch2"), each = ceiling(n_samples / 2))[seq_len(n_samples)],
    stringsAsFactors = FALSE
  )
  list(expr = expr, meta = meta)
}

fx_write_expr_csv <- function(expr, path) {
  df <- data.frame(feature = rownames(expr), expr, check.names = FALSE)
  data.table::fwrite(df, path)
  invisible(path)
}

fx_write_meta_csv <- function(meta, path) {
  data.table::fwrite(meta, path)
  invisible(path)
}

fx_expr_with_duplicate_id <- function(expr) {
  expr2 <- rbind(expr, expr[1, , drop = FALSE])
  rownames(expr2)[nrow(expr2)] <- rownames(expr)[1]
  expr2
}

fx_expr_with_all_na_row <- function(expr) {
  expr[1, ] <- NA_real_
  expr
}

fx_expr_with_zero_variance_row <- function(expr) {
  expr[1, ] <- 5
  expr
}

fx_expr_ensembl_ids <- function(expr) {
  set.seed(99)
  rownames(expr) <- sprintf("ENSG%011d", sample.int(9e8, nrow(expr)))
  expr
}

fx_mkfile <- function(path, type = "text/csv") {
  data.frame(name = basename(path), size = file.info(path)$size, type = type,
             datapath = path, stringsAsFactors = FALSE)
}

fx_html_text <- function(output_value) {
  if (is.list(output_value) && !is.null(output_value$html)) return(as.character(output_value$html))
  as.character(output_value)
}

fx_selected_value <- function(html, select_id) {
  html <- fx_html_text(html)
  block <- regmatches(html, regexpr(
    sprintf('(?s)<select[^>]*id="%s"[^>]*>.*?</select>', select_id),
    html, perl = TRUE
  ))
  if (length(block) == 0 || !nzchar(block)) return(NA_character_)
  m <- regmatches(block, regexpr('<option value="([^"]*)" selected', block, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) return(NA_character_)
  sub('<option value="([^"]*)" selected.*', "\\1", m)
}
