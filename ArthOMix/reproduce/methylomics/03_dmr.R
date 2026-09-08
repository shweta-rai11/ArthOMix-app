## Regenerates the methylomics DMR stage by feeding the already-verified DMP
## regeneration (reproduce/methylomics/output/02_dmp_{F,M}_comparison.csv,
## t_regen/p_bacon_regen/dbeta_regen columns) into DMRcate exactly as the live
## app does (mod_methyl_dmr.R:582-606): lambda=1000, C=2, min.cpgs=3,
## pcutoff="fdr" (DMRcate's own default), seeding p < 0.05 on bacon-corrected
## per-CpG p-values. Compares called regions against dmr_{female,male}_full.csv.
##
## Run from the ArthOMix/ app directory (after reproduce/methylomics/02_dmp.R):
##   Rscript reproduce/methylomics/03_dmr.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))
## methods::new("CpGannotated", ...) needs the DMRcate namespace attached (not
## just loaded via requireNamespace) to resolve the S4 class from top level.
suppressMessages(suppressWarnings({
  library(DMRcate); library(GenomicRanges); library(IRanges); library(S4Vectors)
}))

anno <- methyl_get_annotation("450K")
stopifnot(isTRUE(anno$ok))
a <- anno$anno

run_dmr_sex <- function(sex_letter, sex_label, precomp_path) {
  cat(sprintf("\n--- %s ---\n", sex_label))
  dmp <- read.csv(sprintf("reproduce/methylomics/output/02_dmp_%s_comparison.csv", sex_letter))

  hit <- match(dmp$cpg, rownames(a))
  chr <- a$chr[hit]; pos <- a$pos[hit]
  keep <- !is.na(chr) & !is.na(pos)
  dmp <- dmp[keep, ]; chr <- chr[keep]; pos <- pos[keep]

  ind_fdr <- stats::p.adjust(dmp$p_bacon_regen, "BH")
  gr <- GenomicRanges::GRanges(
    seqnames = chr, ranges = IRanges::IRanges(pos, pos),
    stat = dmp$t_regen, rawpval = dmp$p_bacon_regen, diff = dmp$dbeta_regen,
    ind.fdr = ind_fdr, is.sig = dmp$p_bacon_regen < 0.05
  )
  names(gr) <- dmp$cpg
  cat(sprintf("Seed CpGs (bacon p<0.05): %d\n", sum(gr$is.sig)))
  annot <- methods::new("CpGannotated", ranges = gr)

  dmr_raw <- DMRcate::dmrcate(annot, lambda = 1000, C = 2, min.cpgs = 3, pcutoff = "fdr")
  ranges <- DMRcate::extractRanges(dmr_raw, genome = "hg19")
  dt <- as.data.frame(ranges)
  dt$seqnames <- as.character(dt$seqnames)
  dt$dmr_fdr_regen <- stats::p.adjust(dt$Stouffer, method = "BH")
  cat(sprintf("Candidate DMRs called: %d\n", nrow(dt)))

  precomp <- read.csv(precomp_path)
  cat(sprintf("Precomputed DMRs: %d\n", nrow(precomp)))

  ## Match regions by genomic overlap (region boundaries from a kernel-smoothing
  ## algorithm are not expected to be byte-identical between runs).
  regen_gr <- GenomicRanges::GRanges(dt$seqnames, IRanges::IRanges(dt$start, dt$end))
  precomp_gr <- GenomicRanges::GRanges(precomp$seqnames, IRanges::IRanges(precomp$start, precomp$end))
  ov <- GenomicRanges::findOverlaps(precomp_gr, regen_gr)
  n_matched <- length(unique(S4Vectors::queryHits(ov)))
  cat(sprintf("Precomputed DMRs overlapping a regenerated DMR: %d / %d (%.1f%%)\n",
              n_matched, nrow(precomp), 100 * n_matched / nrow(precomp)))

  sig_precomp <- sum(precomp$dmr_fdr < 0.05, na.rm = TRUE)
  sig_regen <- sum(dt$dmr_fdr_regen < 0.05, na.rm = TRUE)
  sig_precomp_gr <- precomp_gr[which(precomp$dmr_fdr < 0.05)]
  sig_regen_gr <- regen_gr[which(dt$dmr_fdr_regen < 0.05)]
  sig_ov <- if (length(sig_precomp_gr) > 0 && length(sig_regen_gr) > 0) {
    length(unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(sig_precomp_gr, sig_regen_gr))))
  } else 0L
  cat(sprintf("Significant DMRs (FDR<0.05): precomputed=%d, regenerated=%d, overlapping=%d\n",
              sig_precomp, sig_regen, sig_ov))

  write.csv(dt, sprintf("reproduce/methylomics/output/03_dmr_%s_regenerated.csv", sex_letter), row.names = FALSE)
  list(sex = sex_label, n_precomp = nrow(precomp), n_regen = nrow(dt), n_overlap = n_matched,
       sig_precomp = sig_precomp, sig_regen = sig_regen, sig_overlap = sig_ov)
}

res_f <- run_dmr_sex("F", "Female", "data/preloaded/methylomics/tables/script04_dmr_sexstratified/tables/dmr_female_full.csv")
res_m <- run_dmr_sex("M", "Male", "data/preloaded/methylomics/tables/script04_dmr_sexstratified/tables/dmr_male_full.csv")

saveRDS(list(female = res_f, male = res_m), "reproduce/methylomics/output/03_dmr_summary.rds")
cat("\n=== DMR regeneration summary ===\n")
print(rbind(as.data.frame(res_f), as.data.frame(res_m)))
