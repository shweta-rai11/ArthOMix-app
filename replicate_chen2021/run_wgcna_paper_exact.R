## run_wgcna_paper_exact.R - replication using the paper's ACTUAL stated
## methods text: "WGCNA... Soft-thresholding power (beta) was established
## with the pickSoftThreshold R function... Hierarchical clustering analyses
## were then conducted to detect modules (minmodulesize = 66;
## mergecutheight = 0.3)." Only power/minModuleSize/mergeCutHeight are
## stated as chosen - networkType, TOMType, corType, deepSplit,
## pamRespectsDendro and reassignThreshold are never mentioned, so this run
## uses blockwiseModules()'s own package defaults for every one of those
## (networkType = "unsigned", TOMType = "signed", corType = "pearson",
## deepSplit = 2, pamRespectsDendro = TRUE, reassignThreshold = 1e-6) rather
## than this app's own hardcoded values (signed / reassignThreshold = 0 /
## pamRespectsDendro = FALSE), which were never in the paper's text at all.
suppressMessages(library(WGCNA))
WGCNA::disableWGCNAThreads()
options(stringsAsFactors = FALSE)

expr_raw <- read.csv("upload_csv_merged/chen2021_merged_expression_matrix.csv",
                      row.names = 1, check.names = FALSE)
meta <- read.csv("upload_csv_merged/chen2021_merged_sample_metadata.csv", stringsAsFactors = FALSE)
expr <- as.matrix(expr_raw)
stopifnot(identical(colnames(expr), meta$sample))
expr <- expr[stats::complete.cases(expr), , drop = FALSE]

med <- apply(expr, 1, stats::median)
n_genes <- min(5000, nrow(expr))
keep <- order(med, decreasing = TRUE)[seq_len(n_genes)]
expr_f <- expr[keep, , drop = FALSE]
texpr <- t(expr_f)

power <- 7

net <- WGCNA::blockwiseModules(
  texpr,
  power = power,
  ## everything below this line is left at the WGCNA package's own default -
  ## the paper's text never states any of these, so there is no basis to
  ## override them.
  deepSplit = 2, minModuleSize = 66, mergeCutHeight = 0.3,
  randomSeed = 1234,
  numericLabels = FALSE, verbose = 0
)
module_colors <- net$colors
MEs <- net$MEs
cat(sprintf("Modules detected: %d (excl. grey)\n", length(unique(module_colors[module_colors != "grey"]))))
print(sort(table(module_colors), decreasing = TRUE))

traits <- data.frame(row.names = rownames(texpr))
lv <- sort(unique(meta$group))
for (level in lv) traits[[level]] <- as.numeric(meta$group[match(rownames(texpr), meta$sample)] == level)

MEs_use <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
cor_mat <- WGCNA::cor(MEs_use, traits, use = "p")
n_mat <- matrix(nrow(texpr), nrow = nrow(cor_mat), ncol = ncol(cor_mat), dimnames = dimnames(cor_mat))
p_mat <- WGCNA::corPvalueStudent(cor_mat, n_mat)
rownames(cor_mat) <- sub("^ME", "", rownames(cor_mat))
rownames(p_mat) <- sub("^ME", "", rownames(p_mat))

cat("\n--- Module-trait correlation ---\n")
print(round(cor_mat, 2))
cat("\n--- p-values ---\n")
print(signif(p_mat, 2))

cat("\n--- Module sizes, largest first (for comparison to paper's greenyellow = 87 genes) ---\n")
print(sort(table(module_colors[module_colors != "grey"]), decreasing = TRUE))

saveRDS(list(net = net, cor_mat = cor_mat, p_mat = p_mat, meta = meta),
        "wgcna_paper_exact_result.rds")
