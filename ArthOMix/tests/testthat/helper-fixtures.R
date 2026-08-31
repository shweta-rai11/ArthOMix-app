## Deterministic, in-memory test fixtures shared across transcriptomics test
## files. Kept as generator functions (not just static files) so every test
## that needs a "clean" matrix, or one perturbed in one specific way (a
## duplicate ID, an all-NA row, a zero-variance row, an Ensembl-ID rowname
## set), gets it from one seeded source instead of hand-rolling slightly
## different versions per file. tests/fixtures/transcriptomics/ holds the
## on-disk CSV pair generated once from fx_expr_meta() at seed=1, committed
## so a test that specifically wants "a file on disk" (upload/GEO-shaped
## tests) doesn't need to re-materialize it itself.

## n_genes x n_samples numeric matrix + matching sample metadata data.frame
## (sample/group/sex/batch columns), balanced group and sex assignment,
## reproducible for a given seed.
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

## Writes an expr matrix in the exact shape mod_dataset.R's upload path
## expects for a CSV: first column is the gene/feature ID, one column per
## sample thereafter (see mod_dataset.R's expr_raw() reactive).
fx_write_expr_csv <- function(expr, path) {
  df <- data.frame(feature = rownames(expr), expr, check.names = FALSE)
  data.table::fwrite(df, path)
  invisible(path)
}

fx_write_meta_csv <- function(meta, path) {
  data.table::fwrite(meta, path)
  invisible(path)
}

## Perturbations of the base fixture used by boundary/edge-case tests.

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

## Shapes a path into the data.frame a fileInput's `input$x` normally holds
## (name/size/type/datapath), for feeding into session$setInputs() in a
## testServer() block - same convention as
## tests/testthat/test-preprocessing-multi-upload.R's local pp_mkfile().
fx_mkfile <- function(path, type = "text/csv") {
  data.frame(name = basename(path), size = file.info(path)$size, type = type,
             datapath = path, stringsAsFactors = FALSE)
}

## Pulls the currently-selected <option> value out of one non-selectize
## selectInput(id = select_id, ...) inside a block of rendered UI HTML (e.g.
## a testServer() `output$some_renderUI` string) - used to verify a
## guess-the-default-column helper (like mod_dataset.R's local guess_col())
## picked the right column, since testServer() never actually renders a UI
## into live input defaults - session$setInputs() must be explicit, so the
## only observable trace of what a selectInput's `selected =` argument
## computed is in the HTML markup itself.
## Normalizes a testServer() renderUI output into one plain HTML string.
## Shiny represents a rendered UI fragment as a bare character string UNLESS
## it carries an attached html_dependency (e.g. any fontawesome icon()) - in
## that case output$x becomes list(html = <string>, deps = <list>) instead,
## so a plain grepl(pattern, output$x) silently coerces the whole list to
## character (list->character coercion, not the html text) and can compare
## against the wrong thing without erroring. Route every renderUI output
## through this before pattern-matching it.
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
