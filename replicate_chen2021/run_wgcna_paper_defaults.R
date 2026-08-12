## run_wgcna_paper_defaults.R - standalone replication of exactly what
## shiny_app/R/mod_wgcna.R does when its "uploaded dataset" paper-default
## settings are used (Step1 gene_filter_method = "topmedian" n = 5000,
## Step2 network_type = "signed" / power_mode = "manual" power = 7,
## Step3 min_module_size = 66 / merge_cut_height = 0.3 / deep_split = 2 /
## tom_type = "signed" / pam_respects_dendro = FALSE, Step4 trait = "group"),
## run directly against the merged Chen et al. test dataset from
## merge_and_export.R, so the module-trait result can be checked against the
## paper's Figure 1 without going through the Shiny UI by hand.
suppressMessages(library(WGCNA))
WGCNA::disableWGCNAThreads()
options(stringsAsFactors = FALSE)

expr_raw <- read.csv("upload_csv_merged/chen2021_merged_expression_matrix.csv",
                      row.names = 1, check.names = FALSE)
meta <- read.csv("upload_csv_merged/chen2021_merged_sample_metadata.csv", stringsAsFactors = FALSE)
expr <- as.matrix(expr_raw)
stopifnot(identical(colnames(expr), meta$sample))
expr <- expr[stats::complete.cases(expr), , drop = FALSE]
cat(sprintf("Loaded: %d genes x %d samples\n", nrow(expr), ncol(expr)))

## ---- Step 1: gene filter (paper default: top 5000 by median expression) ----
med <- apply(expr, 1, stats::median)
n_genes <- min(5000, nrow(expr))
keep <- order(med, decreasing = TRUE)[seq_len(n_genes)]
expr_f <- expr[keep, , drop = FALSE]
cat(sprintf("After topmedian filter: %d genes\n", nrow(expr_f)))

## No sample-outlier removal (remove_outliers defaults to FALSE)
texpr <- t(expr_f)

## ---- Step 2: power = 7 manual, signed network, pearson ----
power <- 7
network_type <- "signed"

## ---- Step 3: module detection ----
net <- WGCNA::blockwiseModules(
  texpr,
  power = power, networkType = network_type, TOMType = "signed",
  corType = "pearson",
  deepSplit = 2, minModuleSize = 66,
  reassignThreshold = 0, mergeCutHeight = 0.3,
  pamRespectsDendro = FALSE,
  randomSeed = 1234,
  numericLabels = FALSE, maxBlockSize = ncol(texpr) + 1, verbose = 0
)
module_colors <- net$colors
MEs <- net$MEs
cat(sprintf("Modules detected: %d (excl. grey), sizes:\n", length(unique(module_colors[module_colors != "grey"]))))
print(sort(table(module_colors), decreasing = TRUE))

## ---- Step 4: module-trait correlation, "group" trait, HC/RA as two indicator columns ----
traits <- data.frame(row.names = rownames(texpr))
lv <- sort(unique(meta$group))
for (level in lv) {
  traits[[level]] <- as.numeric(meta$group[match(rownames(texpr), meta$sample)] == level)
}

MEs_use <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
cor_mat <- WGCNA::cor(MEs_use, traits, use = "p")
n_mat <- matrix(nrow(texpr), nrow = nrow(cor_mat), ncol = ncol(cor_mat), dimnames = dimnames(cor_mat))
p_mat <- WGCNA::corPvalueStudent(cor_mat, n_mat)
rownames(cor_mat) <- sub("^ME", "", rownames(cor_mat))
rownames(p_mat) <- sub("^ME", "", rownames(p_mat))

cat("\n--- Module-trait correlation (app's current, unordered) ---\n")
print(round(cor_mat, 2))
cat("\n--- p-values ---\n")
print(signif(p_mat, 2))

## ---- what the paper's figure actually orders modules by: hierarchical
## clustering of the eigengenes (orderMEs), not size-rank/alphabetical ----
ordered <- WGCNA::orderMEs(MEs_use)
ord_names <- sub("^ME", "", colnames(ordered))
cat("\n--- Eigengene-clustering order (top-to-bottom as in a dendrogram-ordered heatmap) ---\n")
print(ord_names)

cat("\n--- Disease-associated modules (|r| >= 0.5) ---\n")
hit <- apply(abs(cor_mat) >= 0.5, 1, any)
print(cor_mat[hit, , drop = FALSE])

saveRDS(list(net = net, cor_mat = cor_mat, p_mat = p_mat, ordered = ord_names, meta = meta),
        "wgcna_paper_defaults_result.rds")
cat("\nSaved wgcna_paper_defaults_result.rds\n")
