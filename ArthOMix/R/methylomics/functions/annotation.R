## R/methylomics/functions/annotation.R
## Illumina manifest-based probe annotation (chromosome, SNP overlap) for
## mod_methyl_qc.R's probe filters. Only 450K and EPIC(v1) have an installed

METHYL_ARRAY_TYPES <- c("450K", "EPIC", "EPICv2", "WGBS", "RRBS", "Custom array")

METHYL_ARRAY_TYPES_ILLUMINA <- c("450K", "EPIC", "EPICv2", "Custom array")

METHYL_ANNOTATION_PACKAGES <- list(
  "450K" = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "EPIC" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

.methyl_anno_cache <- new.env(parent = emptyenv())

METHYL_ANNOTATION_SNP_OBJECT <- "SNPs.147CommonSingle"

methyl_get_annotation <- function(array_type) {
  pkg <- METHYL_ANNOTATION_PACKAGES[[array_type]]
  if (is.null(pkg)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf(
      "No Bioconductor manifest annotation is available for %s in this deployment - SNP and sex-chromosome probe filters, and the raw-intensity sex check, are unavailable for this array type.",
      array_type
    )))
  }
  cached <- .methyl_anno_cache[[pkg]]
  if (!is.null(cached)) return(list(ok = TRUE, anno = cached, reason = NULL))

  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf(
      "The %s annotation package is not installed in this deployment.", pkg
    )))
  }
  anno <- tryCatch({
    e <- new.env(parent = emptyenv())
    utils::data(list = c("Locations", "Manifest", "Other", METHYL_ANNOTATION_SNP_OBJECT), package = pkg, envir = e)
    loc <- as.data.frame(e$Locations)
    man <- as.data.frame(e$Manifest)
    snp <- as.data.frame(e[[METHYL_ANNOTATION_SNP_OBJECT]])
    oth <- as.data.frame(e$Other)
    ids <- rownames(loc)
    gene <- vapply(strsplit(oth[ids, "UCSC_RefGene_Name"], ";"), function(g) {
      if (length(g) == 0 || !nzchar(g[1])) NA_character_ else g[1]
    }, character(1))
    data.frame(
      chr = loc[ids, "chr"], pos = loc[ids, "pos"],
      Type = man[ids, "Type"],
      Probe_rs = snp[ids, "Probe_rs"], CpG_rs = snp[ids, "CpG_rs"], SBE_rs = snp[ids, "SBE_rs"],
      gene = gene,
      row.names = ids, stringsAsFactors = FALSE
    )
  }, error = function(e) e)
  if (inherits(anno, "error")) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf("Could not extract probe annotation from %s: %s", pkg, conditionMessage(anno))))
  }
  .methyl_anno_cache[[pkg]] <- anno
  list(ok = TRUE, anno = anno, reason = NULL)
}

methyl_probe_is_cpg <- function(probe_ids) {
  grepl("^cg", probe_ids, ignore.case = TRUE)
}
