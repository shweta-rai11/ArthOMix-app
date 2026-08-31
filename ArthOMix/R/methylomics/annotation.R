## R/methylomics/annotation.R
## Illumina manifest-based probe annotation (chromosome, SNP overlap) for
## mod_methyl_qc.R's probe filters. Only 450K and EPIC(v1) have an installed
## Bioconductor annotation package; other array types degrade to
## list(ok = FALSE, reason = ...) instead of guessing.

METHYL_ARRAY_TYPES <- c("450K", "EPIC", "EPICv2", "WGBS", "RRBS", "Custom array")

## Array types where Illumina probe-ID prefixes (cg.../ch.../rs...) apply - WGBS/RRBS
## probe "IDs" are genomic coordinates, so prefix/manifest filters are hidden for them.
METHYL_ARRAY_TYPES_ILLUMINA <- c("450K", "EPIC", "EPICv2", "Custom array")

## EPICv2 deliberately NOT mapped here even though
## IlluminaHumanMethylationEPICv2anno.20a1.hg38 is genuinely installed in
## some deployments (verified: it does expose the same Locations/Manifest/
## Other/SNPs.147CommonSingle objects methyl_get_annotation() below reads).
## Two things make it unsafe to wire in as a drop-in: (1) it's hg38, while
## every other annotation package here (and DMRcate's downstream calls) is
## hg19 - genomic coordinates would silently mix builds; (2) its probe IDs
## carry a replicate suffix (e.g. "cgXXXXXXXX_BC11"), not the plain
## "cgXXXXXXXX" IDs 450K/EPIC(v1) and most uploaded EPICv2 beta matrices use -
## an ID-based join against this manifest would silently return all-NA
## annotation instead of erroring. mod_methyl_dataset.R's array_type_note
## instead shows an explicit inline warning that EPICv2 has no manifest here.
METHYL_ANNOTATION_PACKAGES <- list(
  "450K" = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "EPIC" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

## Cached per array type for the process lifetime - EPIC's Manifest alone is ~866k rows.
.methyl_anno_cache <- new.env(parent = emptyenv())

## dbSNP build shared by the 450K/EPIC annotation packages - used only as a
## probe-overlaps-known-SNP flag, not an allele-frequency filter, so no need for the newest build.
METHYL_ANNOTATION_SNP_OBJECT <- "SNPs.147CommonSingle"

## Returns list(ok, anno, reason) - `anno` keyed by probe ID with
## chr/pos/Type/Probe_rs/CpG_rs/SBE_rs columns.
##
## Deliberately avoids minfi::getAnnotation(): loading its wrapper object
## requires library()-attaching the package, which pulls in Biostrings and
## masks base::strsplit() app-wide. Reading the underlying Locations/Manifest/
## SNPs.*CommonSingle objects directly avoids that while getting the same data.
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
    ## One representative gene symbol per probe: first token of the semicolon-
    ## separated UCSC_RefGene_Name list (same simplification as ChAMPdata::probe.features$gene).
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

## "cg" prefix = CpG probe. Works from the probe ID alone, no manifest needed - but
## only meaningful for an actual Illumina array (see METHYL_ARRAY_TYPES_ILLUMINA).
methyl_probe_is_cpg <- function(probe_ids) {
  grepl("^cg", probe_ids, ignore.case = TRUE)
}
