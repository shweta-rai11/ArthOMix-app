## ArthOChat sub-module grounding verification (live, manual — NOT part of the
## testthat suite: it requires a running local Ollama server with
## ARTHOMIX_OLLAMA_MODEL pulled, makes real live-model calls, and is
## non-deterministic, so it must never run in CI or the automated test suite.
## Run it by hand from the ArthOMix/ directory:
##   Rscript tests/arthochat_verification/test_submodule_grounding.R
##
## Method — "known answer in, verify it comes back out":
## for every sub-module id in every vertical's own registry, a synthetic
## results object with one unique numeric check_value is injected, the exact
## same build_*_context()/build_scoped_assistant_context() code path
## mod_arthochat.R itself uses is called to build a real system prompt scoped
## to that one sub-module, and a real request is sent to the live Ollama
## model asking it to report that value back. The response is checked for:
##   - correct:      does it contain THIS sub-module's own check_value?
##   - isolation_ok: does it NOT also contain a neighboring sub-module's
##                   check_value (which was never shown to it if the
##                   focus_id scoping is working)?
## Every value used is one this script invented itself, so grading is exact
## string matching against a known ground truth — not a subjective read of
## the model's prose.
##
## Output: results.csv in this same directory (one row per sub-module),
## plus a summary printed to the console.

setwd(normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..")))
if (!file.exists("global.R")) setwd("../..")  ## fallback when sourced via Rscript with a relative path

suppressMessages(suppressWarnings({
  source("data_paths.R")
  source("global.R")
  for (f in list.files("R", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)) source(f)
}))

OUT_CSV <- file.path("tests", "arthochat_verification", "results.csv")

## ---- exact ids, verified against each mod_*_config in the repo (not the registry-file comments) ----
tx_ids <- c("overview", "preprocessing", "dge", "wgcna", "candidates", "mr", "coloc", "featureselection",
            "diagnostic", "interaction", "crosstissue", "crossancestry", "enrichment", "deconvolution",
            "nomogram", "biomarkercard")
mx_ids <- c("qc", "normalization", "celltype", "dmp", "dmr", "interaction", "wgcna", "candidates",
            "featureselection", "mr", "coloc", "diagnostic", "validation", "biomarkercard")
cx_ids <- c("integration", "biomarkerconv", "mrstage")
mo_ids <- c("overview", "integration", "stratification", "biomarker", "concordance", "pathway", "biomarkercard", "summary")

## ---- synthetic dataset stubs, shaped to satisfy each build_*_context()'s own field reads ----
tx_dataset <- list(expr = matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4))),
                    meta = data.frame(sample = paste0("s", 1:4), group = c("RA", "RA", "HC", "HC"), sex = c("F", "M", "F", "M")),
                    source = "synthetic test dataset")
mx_dataset <- list(source = "synthetic methylomics dataset", array_type = "EPIC",
                    beta = matrix(runif(20), nrow = 5, dimnames = list(paste0("cg", 1:5), paste0("s", 1:4))))
mo_dataset <- list(active = TRUE, source = "synthetic multi-omics dataset",
                    layers = list(rna = TRUE, methylation = TRUE),
                    sample_meta = data.frame(sample = paste0("s", 1:4)))

mk_generic <- function(check_value, marker) list(check_value = check_value, marker = marker, n_samples = 24L)

## Cross-Omics' "integration" id is rendered by the bespoke .format_cx_integration()
## formatter (submodules_registry.R), not the generic key/value block, so its
## synthetic result has to match that function's expected shape instead.
mk_cx_integration <- function(check_value, marker) list(
  params = list(input_mode = "synthetic", sample_matching = marker),
  summary = list(n_genes = check_value, n_deg = 10, n_dmg = 8, n_integrated = 3,
                 counts = list(`Hyper + Down` = 1, `Hypo + Up` = 1, `Hyper + Up` = 1, `Hypo + Down` = 1))
)

build_results_list <- function(ids, base) {
  out <- list()
  for (i in seq_along(ids)) {
    id <- ids[i]
    check_value <- base + i
    marker <- sprintf("MARK_%s_%d", id, check_value)
    out[[id]] <- if (identical(id, "integration") && base == 5000) mk_cx_integration(check_value, marker) else mk_generic(check_value, marker)
  }
  out
}

tx_results <- build_results_list(tx_ids, 1000)
mx_results <- build_results_list(mx_ids, 2000)
cx_results <- build_results_list(cx_ids, 5000)   ## base 5000 is what triggers the bespoke integration formatter above
mo_results <- build_results_list(mo_ids, 4000)

## Talks to Ollama's native /api/chat directly (bypassing ellmer) so the only
## thing under test is the system prompt this app actually generates, with no
## extra library behavior in between. Mirrors mod_arthochat.R's own
## api_args = list(think = FALSE) setting.
ollama_chat <- function(sys_prompt, user_msg) {
  t0 <- Sys.time()
  resp <- httr2::request(paste0(ollama_base_url(), "/api/chat")) %>%
    httr2::req_body_json(list(
      model = ARTHOMIX_OLLAMA_MODEL, stream = FALSE, think = FALSE,
      messages = list(list(role = "system", content = sys_prompt), list(role = "user", content = user_msg))
    )) %>%
    httr2::req_timeout(60) %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(content = resp$message$content %||% "", elapsed = elapsed,
       total_duration_s = (resp$total_duration %||% NA) / 1e9,
       eval_count = resp$eval_count %||% NA)
}

run_vertical <- function(vertical, ids, base, build_ctx_call) {
  rows <- list()
  for (i in seq_along(ids)) {
    id <- ids[i]
    check_value <- base + i
    neighbor_idx <- if (i == length(ids)) i - 1 else i + 1
    neighbor_check_value <- base + neighbor_idx

    ctx <- build_ctx_call(id)
    view <- list(module = vertical, view_label = sprintf("%s / test", vertical), submodule_id = id)
    sys_prompt <- paste(ARTHOCHAT_SYSTEM_PROMPT, "", sprintf("## Current view: %s", view$view_label), "", ctx, sep = "\n")
    question <- "Report the exact check_value number for this sub-module's results, and nothing else's."

    res <- tryCatch(
      ollama_chat(sys_prompt, question),
      error = function(e) list(content = paste("ERROR:", conditionMessage(e)), elapsed = NA, total_duration_s = NA, eval_count = NA)
    )

    has_own <- grepl(as.character(check_value), res$content, fixed = TRUE)
    leaked_neighbor <- grepl(as.character(neighbor_check_value), res$content, fixed = TRUE)

    rows[[id]] <- data.frame(
      vertical = vertical, id = id, check_value = check_value,
      correct = has_own, isolation_ok = !leaked_neighbor,
      elapsed_s = round(res$elapsed, 2), ollama_total_s = round(res$total_duration_s, 2),
      eval_tokens = res$eval_count,
      answer = substr(gsub("\n", " ", res$content), 1, 160),
      stringsAsFactors = FALSE
    )
    cat(sprintf("[%s/%s] correct=%s isolation_ok=%s  %.1fs\n", vertical, id, has_own, !leaked_neighbor, res$elapsed))
  }
  do.call(rbind, rows)
}

tx_out <- run_vertical("transcriptomics", tx_ids, 1000, function(id) build_tx_context(tx_dataset, tx_results, focus_id = id))
mx_out <- run_vertical("methylomics", mx_ids, 2000, function(id) build_mx_context(mx_dataset, mx_results, focus_id = id))
cx_out <- run_vertical("crossomics", cx_ids, 5000, function(id) build_cx_context(cx_results, focus_id = id))
mo_out <- run_vertical("multiomics", mo_ids, 4000, function(id) build_mo_context(mo_dataset, mo_results, focus_id = id))

all_out <- rbind(tx_out, mx_out, cx_out, mo_out)
write.csv(all_out, OUT_CSV, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(all_out), OUT_CSV))

cat("\n==================== SUMMARY ====================\n")
cat(sprintf("Total sub-modules tested: %d\n", nrow(all_out)))
cat(sprintf("Correct (own check_value present): %d / %d\n", sum(all_out$correct), nrow(all_out)))
cat(sprintf("Isolation OK (neighbor value absent): %d / %d\n", sum(all_out$isolation_ok), nrow(all_out)))
cat(sprintf("Latency (s) - mean: %.2f, median: %.2f, min: %.2f, max: %.2f\n",
            mean(all_out$elapsed_s, na.rm = TRUE), median(all_out$elapsed_s, na.rm = TRUE),
            min(all_out$elapsed_s, na.rm = TRUE), max(all_out$elapsed_s, na.rm = TRUE)))
cat("\nBy vertical:\n")
print(aggregate(cbind(correct, isolation_ok, elapsed_s) ~ vertical, all_out, function(x) round(mean(as.numeric(x)), 3)))
cat("\nAny failures:\n")
print(all_out[!all_out$correct | !all_out$isolation_ok, c("vertical", "id", "check_value", "correct", "isolation_ok", "answer")])
