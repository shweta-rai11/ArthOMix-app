## Regenerates the methylomics WGCNA co-methylation network stage from
## beta_raw.rds, using the exact "preloaded"-mode parameters the live app
## itself encodes as its own defaults (mod_methyl_wgcna.R): sex-stratified,
## MAD top-20000 CpGs (M-value scale), no residualization, 5% max missingness,
## pearson correlation, signed network/TOM, custom power vector
## 1,2,...,10,12,14,16,18,20 targeting R^2>=0.85, minModuleSize=20, deepSplit=2,
## mergeCutHeight=0.25, maxBlockSize=5000, pamStage/pamRespectsDendro=TRUE,
## reassignThreshold=1e-6, minKMEtoStay=0.3, minCoreKME=0.5, seed=1234.
## Compares module sizes/assignments against module_assignment_{female,male}.csv.
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/methylomics/04_wgcna.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

beta <- readRDS("data/preloaded/methylomics/matrix/beta_raw.rds")
pheno <- as.data.frame(readRDS("data/preloaded/methylomics/matrix/pheno.rds"))
rownames(pheno) <- pheno$gsm

run_wgcna_sex <- function(sex_letter, sex_label, precomp_path) {
  cat(sprintf("\n--- %s ---\n", sex_label))
  sub <- methyl_qc_subgroup_filter(beta, pheno, "sex", sex_letter, min_n = 15)
  mat <- sub$mat
  cat(sprintf("Samples in stratum: %d\n", ncol(mat)))

  f_miss <- methyl_filter_missing(mat, 0.05)
  mat <- mat[f_miss$keep, , drop = FALSE]
  mat[mat < 1e-6] <- 1e-6; mat[mat > 1 - 1e-6] <- 1 - 1e-6
  m <- log2(mat / (1 - mat))
  rm(mat); gc(FALSE)
  cat(sprintf("CpGs after missingness filter: %d\n", nrow(m)))

  vt <- mx_wgcna_top_variable(m, method = "mad", top_n = 20000)
  m <- vt$mat
  gsg <- WGCNA::goodSamplesGenes(t(m), verbose = 0)
  if (!isTRUE(gsg$allOK)) m <- m[gsg$goodGenes, gsg$goodSamples, drop = FALSE]
  cat(sprintf("CpGs after top-20000 MAD + goodSamplesGenes: %d, samples: %d\n", nrow(m), ncol(m)))

  texpr <- t(m)
  powers <- c(1:10, 12, 14, 16, 18, 20)
  sft <- WGCNA::pickSoftThreshold(texpr, powerVector = powers, networkType = "signed",
                                   corFnc = "cor", RsquaredCut = 0.85, verbose = 0)
  fi <- sft$fitIndices
  max_r_sq <- suppressWarnings(max(fi$SFT.R.sq, na.rm = TRUE))
  reached <- is.finite(max_r_sq) && max_r_sq >= 0.85
  power <- if (reached) min(fi$Power[fi$SFT.R.sq >= 0.85]) else fi$Power[which.max(fi$SFT.R.sq)]
  cat(sprintf("Soft-threshold power selected: %d (R^2=%.3f, cutoff reached=%s)\n", power, max_r_sq, reached))

  net <- WGCNA::blockwiseModules(
    texpr, power = power, networkType = "signed", TOMType = "signed",
    corType = "pearson", deepSplit = 2, minModuleSize = 20,
    mergeCutHeight = 0.25, maxBlockSize = 5000,
    pamStage = TRUE, pamRespectsDendro = TRUE,
    reassignThreshold = 1e-6, minKMEtoStay = 0.3, minCoreKME = 0.5,
    numericLabels = FALSE, saveTOMs = FALSE, randomSeed = 1234, verbose = 0
  )
  module_colors <- net$colors
  names(module_colors) <- colnames(texpr)
  cat(sprintf("Modules found: %d (sizes: %s)\n", length(unique(module_colors)),
              paste(sort(table(module_colors), decreasing = TRUE), collapse = ",")))

  regen <- data.frame(cpg = names(module_colors), module_regen = unname(module_colors), stringsAsFactors = FALSE)
  precomp <- read.csv(precomp_path)
  colnames(precomp)[1] <- "cpg"
  mod_col <- intersect(c("module", "moduleColor", "colors"), colnames(precomp))[1]
  cmp <- merge(precomp, regen, by = "cpg")
  cat(sprintf("Matched %d / %d precomputed CpGs (regen tested %d CpGs total)\n", nrow(cmp), nrow(precomp), nrow(regen)))

  ## Module labels/colors are arbitrary and won't match by name across runs -
  ## measure agreement via Adjusted Rand Index (co-clustering agreement),
  ## the standard way to compare two independent module-assignment/clusterings.
  ari <- tryCatch(mclust::adjustedRandIndex(cmp[[mod_col]], cmp$module_regen), error = function(e) NA_real_)
  cat(sprintf("Adjusted Rand Index (precomputed vs regenerated module assignment): %.4f\n", ari))

  write.csv(cmp, sprintf("reproduce/methylomics/output/04_wgcna_%s_comparison.csv", sex_letter), row.names = FALSE)
  list(sex = sex_label, n_samples = ncol(m), n_cpgs = nrow(m), power = power,
       n_modules_regen = length(unique(module_colors)), n_matched = nrow(cmp), ari = ari)
}

res_f <- run_wgcna_sex("F", "Female", "data/preloaded/methylomics/tables/script05_wgcna_sexstratified/tables/module_assignment_female.csv")
res_m <- run_wgcna_sex("M", "Male", "data/preloaded/methylomics/tables/script05_wgcna_sexstratified/tables/module_assignment_male.csv")

saveRDS(list(female = res_f, male = res_m), "reproduce/methylomics/output/04_wgcna_summary.rds")
cat("\n=== WGCNA regeneration summary ===\n")
print(rbind(as.data.frame(res_f), as.data.frame(res_m)))
