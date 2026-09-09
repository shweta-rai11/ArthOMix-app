## Tests the claim (carried in this session's REPRODUCIBILITY.md as "prior
## session, not re-verified") that transcriptomics DEG reproduces from the
## app's own live limma path. Uses the already-bundled combined_expr_
## batchcorrected.rds directly (the same object DEFAULT_EXPR_RDS points at,
## and what a user sees as "the active dataset" by default) - this tests the
## DEG stage itself, not the upstream batch-correction step (which has no
## preloaded-mode pinned hyperparameters in the live UI, unlike DEG).
##
## Live algorithm replicated exactly from mod_dge.R's compute_dge_fit():
## sex-filtered, design = model.matrix(~0+grp), limma::arrayWeights,
## lmFit(weights=aw), makeContrasts("RA-HC"), eBayes (posthoc mode, the
## default - not treat), topTable sort by P.
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/transcriptomics/01_deg.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

expr <- readRDS("data/preloaded/transcriptomics/processed/combined_expr_batchcorrected.rds")
meta <- read.csv("data/preloaded/transcriptomics/processed/combined_meta.csv")
cat(sprintf("Bundled matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("Metadata: n=%d (%s)\n", nrow(meta), paste(names(table(meta$sex)), table(meta$sex), sep = "=", collapse = ", ")))

run_deg_sex <- function(sex_val, precomp_path) {
  cat(sprintf("\n--- %s ---\n", sex_val))
  m <- meta[meta$sex == sex_val & meta$group %in% c("HC", "RA"), ]
  common <- intersect(colnames(expr), m$sample)
  m <- m[match(common, m$sample), ]
  e <- expr[, common, drop = FALSE]
  grp <- factor(m$group, levels = c("HC", "RA"))
  cat(sprintf("n=%d (%s)\n", length(grp), paste(table(grp), collapse = "/")))

  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  aw <- limma::arrayWeights(e, design)
  fit <- limma::lmFit(e, design, weights = aw)
  cm <- limma::makeContrasts(contrasts = "RA-HC", levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  regen <- data.frame(gene = rownames(tt), logFC_regen = tt$logFC, t_regen = tt$t,
                       P_regen = tt$P.Value, adjP_regen = tt$adj.P.Val, row.names = NULL)

  precomp <- read.csv(precomp_path)
  cmp <- merge(precomp, regen, by = "gene")
  cat(sprintf("Matched %d / %d precomputed genes\n", nrow(cmp), nrow(precomp)))

  r_logfc <- cor(cmp$logFC, cmp$logFC_regen, use = "complete.obs")
  r_t <- cor(cmp$t, cmp$t_regen, use = "complete.obs")
  n_sig_precomp <- sum(cmp$adj.P.Val < 0.05, na.rm = TRUE)
  n_sig_regen <- sum(cmp$adjP_regen < 0.05, na.rm = TRUE)
  n_sig_overlap <- sum(cmp$adj.P.Val < 0.05 & cmp$adjP_regen < 0.05, na.rm = TRUE)

  cat(sprintf("logFC correlation: r = %.6f\n", r_logfc))
  cat(sprintf("t-statistic correlation: r = %.6f\n", r_t))
  cat(sprintf("Significant genes (adj.P<0.05): precomputed=%d, regenerated=%d, overlap=%d\n",
              n_sig_precomp, n_sig_regen, n_sig_overlap))
  identical_logfc <- isTRUE(all.equal(cmp$logFC, cmp$logFC_regen, tolerance = 1e-6))
  cat(sprintf("logFC byte/near-identical (tol 1e-6): %s\n", identical_logfc))

  list(sex = sex_val, n_matched = nrow(cmp), r_logfc = r_logfc, r_t = r_t,
       n_sig_precomp = n_sig_precomp, n_sig_regen = n_sig_regen, n_sig_overlap = n_sig_overlap,
       near_identical = identical_logfc)
}

res_f <- run_deg_sex("F", "data/preloaded/transcriptomics/results/tables/DEG_female_full.csv")
res_m <- run_deg_sex("M", "data/preloaded/transcriptomics/results/tables/DEG_male_full.csv")

cat("\n=== DEG reproduction summary ===\n")
print(rbind(as.data.frame(res_f), as.data.frame(res_m)))
saveRDS(list(female = res_f, male = res_m), "reproduce/transcriptomics/output/01_deg_summary.rds")
