setwd("/Users/swetarai/ArthOMix/ArthOMix")
suppressMessages(suppressWarnings({
  source("data_paths.R")
  source("global.R")
  for (f in list.files("R", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)) source(f)
}))

OUT_CSV <- "/private/tmp/claude-501/-Users-swetarai-ArthOMix/4d90ea6c-5cd2-4ae4-b26a-b4b077c66b2a/scratchpad/full_205_results.csv"
PROGRESS_LOG <- "/private/tmp/claude-501/-Users-swetarai-ArthOMix/4d90ea6c-5cd2-4ae4-b26a-b4b077c66b2a/scratchpad/full_205_progress.txt"
unlink(PROGRESS_LOG)

tx <- list(
  overview="Overview and Datasets", preprocessing="Preprocessing and Batch Correction", dge="Differential Expression",
  wgcna="WGCNA Co-expression Network", candidates="Candidate Gene Identification", mr="Mendelian Randomization",
  coloc="Colocalization", featureselection="Feature Selection", diagnostic="Diagnostic Model",
  interaction="Sex Interaction Analysis", crosstissue="Cross-Tissue Validation", crossancestry="Cross-Ancestry Validation",
  enrichment="Functional Enrichment", deconvolution="Immune Deconvolution", nomogram="Clinical Utility Nomogram",
  biomarkercard="Biomarker Card"
)
mx <- list(
  qc="Quality Control", normalization="Normalization", celltype="Cell-Type Deconvolution",
  dmp="Differential Methylation (DMPs)", dmr="Differentially Methylated Regions (DMRs)", interaction="Sex Interaction Analysis",
  wgcna="WGCNA (Co-Methylation Network)", candidates="Candidate CpGs (Module-DMR Overlap)", featureselection="ML Feature Selection",
  mr="Mendelian Randomization", coloc="Colocalisation", diagnostic="Diagnostic Classifier",
  validation="Validation", biomarkercard="Biomarker Card"
)
cx <- list(integration="Expression and Methylation", biomarkerconv="Biomarker Convergence", mrstage="Cross-Omics MR")
mo <- list(
  overview="Cohort Harmonization", integration="Multi-omics Integration (DIABLO & SNF)", stratification="SNF Clustering",
  biomarker="Biomarker Discovery", concordance="Gene-CpG Concordance", pathway="Pathways",
  biomarkercard="Biomarker Card", summary="Results Summary & Reproducibility"
)

tx_dataset <- list(expr = matrix(1:20, nrow = 5, dimnames = list(paste0("g",1:5), paste0("s",1:4))),
                    meta = data.frame(sample = paste0("s",1:4), group = c("RA","RA","HC","HC"), sex = c("F","M","F","M")),
                    source = "GSE93272-synthetic-test-cohort")
mx_dataset <- list(source = "GSE42861-synthetic-test-cohort", array_type = "EPIC",
                    beta = matrix(runif(20), nrow = 5, dimnames = list(paste0("cg",1:5), paste0("s",1:4))))
mo_dataset <- list(active = TRUE, source = "Tao2021-synthetic-test-cohort",
                    layers = list(rna = TRUE, methylation = TRUE),
                    sample_meta = data.frame(sample = paste0("s",1:4)))

mk_generic <- function(check_value) list(check_value = check_value, secondary_metric = check_value * 2 + 7,
                                          marker = sprintf("MARK_%d", check_value), n_samples = 24L)
mk_cx_integration <- function(check_value) list(
  params = list(input_mode = "synthetic", sample_matching = sprintf("MARK_%d", check_value)),
  summary = list(n_genes = check_value, n_deg = check_value * 2 + 7, n_dmg = 8, n_integrated = 3,
                 counts = list(`Hyper + Down` = 1, `Hypo + Up` = 1, `Hyper + Up` = 1, `Hypo + Down` = 1))
)
build_results <- function(ids, base) {
  out <- list()
  for (i in seq_along(ids)) {
    cv <- base + i
    out[[ids[i]]] <- if (identical(ids[i], "integration") && base == 5000) mk_cx_integration(cv) else mk_generic(cv)
  }
  out
}
tx_results <- build_results(names(tx), 1000)
mx_results <- build_results(names(mx), 2000)
cx_results <- build_results(names(cx), 5000)
mo_results <- build_results(names(mo), 4000)

ollama_chat <- function(sys_prompt, user_msg) {
  t0 <- Sys.time()
  resp <- tryCatch(
    httr2::request(paste0(ollama_base_url(), "/api/chat")) %>%
      httr2::req_body_json(list(model = ARTHOMIX_OLLAMA_MODEL, stream = FALSE, think = FALSE,
                                 messages = list(list(role = "system", content = sys_prompt),
                                                 list(role = "user", content = user_msg)))) %>%
      httr2::req_timeout(90) %>%
      httr2::req_perform() %>%
      httr2::resp_body_json(),
    error = function(e) list(message = list(content = paste("ERROR:", conditionMessage(e))), total_duration = NA, eval_count = NA)
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(content = resp$message$content %||% "", elapsed = round(elapsed, 2),
       ollama_s = round((resp$total_duration %||% NA) / 1e9, 2), eval_tokens = resp$eval_count %||% NA)
}

kw_hit <- function(answer, title) {
  words <- unique(tolower(unlist(strsplit(gsub("[^A-Za-z ]", " ", title), "\\s+"))))
  words <- words[nchar(words) >= 4]
  any(vapply(words, function(w) grepl(w, tolower(answer), fixed = TRUE), logical(1)))
}

run_one_submodule <- function(vertical, id, title, other_ids_titles, base,
                                dataset_source, ctx_fn) {
  cv <- base + which(names(other_ids_titles) == id)
  sm <- cv * 2 + 7
  ctx <- ctx_fn(id)
  view_label <- sprintf("%s / %s", vertical, title)
  sys_prompt <- paste(ARTHOCHAT_SYSTEM_PROMPT, "", sprintf("## Current view: %s", view_label), "", ctx, sep = "\n")

  other_ids <- setdiff(names(other_ids_titles), id)
  neighbor_id <- other_ids[1]
  neighbor_title <- other_ids_titles[[neighbor_id]]
  neighbor_cv <- base + which(names(other_ids_titles) == neighbor_id)

  qs <- list(
    list(type = "Q1_basic",         q = sprintf("In plain terms, what does the \"%s\" analysis do in this app, and what kind of output does it produce?", title)),
    list(type = "Q2_data_scope",    q = sprintf("Has the \"%s\" analysis been run in this session, and what dataset is currently loaded?", title)),
    list(type = "Q3_numeric_basic", q = "Report the exact check_value number recorded for this sub-module's results, and nothing else's."),
    list(type = "Q4_numeric_advanced", q = "This sub-module's results record a check_value and a secondary_metric. What is the sum of the two?"),
    list(type = "Q5_boundary",      q = sprintf("How does the \"%s\" sub-module's result compare to what's shown here?", neighbor_title))
  )

  rows <- list()
  for (qi in seq_along(qs)) {
    qdef <- qs[[qi]]
    res <- ollama_chat(sys_prompt, qdef$q)
    ans <- res$content

    graded <- NA; correct <- NA; detail <- ""
    if (qdef$type == "Q1_basic") {
      graded <- FALSE; correct <- NA; detail <- sprintf("keyword_hit=%s (informational only, not a correctness grade)", kw_hit(ans, title))
    } else if (qdef$type == "Q2_data_scope") {
      graded <- TRUE; correct <- grepl(dataset_source, ans, fixed = TRUE)
    } else if (qdef$type == "Q3_numeric_basic") {
      graded <- TRUE; correct <- grepl(as.character(cv), ans, fixed = TRUE)
    } else if (qdef$type == "Q4_numeric_advanced") {
      graded <- TRUE; correct <- grepl(as.character(cv + sm), ans, fixed = TRUE)
    } else if (qdef$type == "Q5_boundary") {
      graded <- TRUE; correct <- !grepl(as.character(neighbor_cv), ans, fixed = TRUE)
      detail <- sprintf("neighbor=%s neighbor_check_value=%d", neighbor_id, neighbor_cv)
    }

    rows[[length(rows) + 1]] <- data.frame(
      vertical = vertical, id = id, title = title, question_type = qdef$type,
      question = qdef$q, answer = gsub("\n", " ", ans),
      graded = graded, correct = correct, detail = detail,
      elapsed_s = res$elapsed, ollama_s = res$ollama_s, eval_tokens = res$eval_tokens,
      stringsAsFactors = FALSE
    )
    cat(sprintf("[%s/%s] %-20s graded=%-5s correct=%-5s %.1fs\n", vertical, id, qdef$type, graded, correct, res$elapsed),
        file = PROGRESS_LOG, append = TRUE)
  }
  do.call(rbind, rows)
}

run_vertical <- function(vertical, titles, base, dataset_source, ctx_fn) {
  all_rows <- list()
  for (id in names(titles)) {
    all_rows[[id]] <- run_one_submodule(vertical, id, titles[[id]], titles, base, dataset_source, ctx_fn)
  }
  do.call(rbind, all_rows)
}

tx_out <- run_vertical("transcriptomics", tx, 1000, "GSE93272-synthetic-test-cohort",
                        function(id) build_tx_context(tx_dataset, tx_results, focus_id = id))
mx_out <- run_vertical("methylomics", mx, 2000, "GSE42861-synthetic-test-cohort",
                        function(id) build_mx_context(mx_dataset, mx_results, focus_id = id))
cx_out <- run_vertical("crossomics", cx, 5000, "synthetic",
                        function(id) build_cx_context(cx_results, focus_id = id))
mo_out <- run_vertical("multiomics", mo, 4000, "Tao2021-synthetic-test-cohort",
                        function(id) build_mo_context(mo_dataset, mo_results, focus_id = id))

all_out <- rbind(tx_out, mx_out, cx_out, mo_out)
write.csv(all_out, OUT_CSV, row.names = FALSE)

cat("\n\n==================== SUMMARY ====================\n")
cat(sprintf("Total questions: %d (expected 205)\n", nrow(all_out)))
graded_df <- all_out[all_out$graded == TRUE, ]
cat(sprintf("Graded questions (Q2-Q5): %d\n", nrow(graded_df)))
cat(sprintf("Overall accuracy on graded questions: %d / %d (%.1f%%)\n",
            sum(graded_df$correct), nrow(graded_df), 100*mean(graded_df$correct)))
cat("\nAccuracy by question type:\n")
print(aggregate(correct ~ question_type, graded_df, function(x) sprintf("%d/%d (%.1f%%)", sum(x), length(x), 100*mean(x))))
cat("\nAccuracy by vertical (graded only):\n")
print(aggregate(correct ~ vertical, graded_df, function(x) sprintf("%d/%d (%.1f%%)", sum(x), length(x), 100*mean(x))))
cat(sprintf("\nLatency (s) overall - mean: %.2f, median: %.2f, min: %.2f, max: %.2f\n",
            mean(all_out$elapsed_s, na.rm=TRUE), median(all_out$elapsed_s, na.rm=TRUE),
            min(all_out$elapsed_s, na.rm=TRUE), max(all_out$elapsed_s, na.rm=TRUE)))
cat("\nLatency (s) by question type:\n")
print(aggregate(elapsed_s ~ question_type, all_out, function(x) round(mean(x),2)))
cat("\nAll failures on graded questions:\n")
print(graded_df[!graded_df$correct, c("vertical","id","question_type","answer","detail")])
