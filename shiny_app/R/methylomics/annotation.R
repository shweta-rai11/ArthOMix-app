## R/methylomics/annotation.R
## Illumina manifest-based probe annotation (chromosome, SNP overlap), used
## by the probe-filtering options in mod_methyl_qc.R (shiny_app/). Only
## covers array types with an installed Bioconductor annotation package -
## 450K and EPIC (v1). EPICv2/WGBS/RRBS/Custom array don't have one
## available in this deployment, so every annotation-dependent filter
## degrades to list(ok = FALSE, reason = ...) instead of guessing at probe
## positions or SNP overlap.

METHYL_ARRAY_TYPES <- c("450K", "EPIC", "EPICv2", "WGBS", "RRBS", "Custom array")

## Array types where Illumina probe-ID conventions (cg.../ch.../rs... prefixes,
## a manifest keyed by those IDs) apply at all - WGBS/RRBS probe "IDs" are
## typically genomic coordinates, not Illumina probe IDs, so manifest- and
## prefix-based filters are meaningless for them and are hidden in the UI.
METHYL_ARRAY_TYPES_ILLUMINA <- c("450K", "EPIC", "EPICv2", "Custom array")

METHYL_ANNOTATION_PACKAGES <- list(
  "450K" = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "EPIC" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

## Cached per array type for the life of the R process - not cheap to
## rebuild on every QC/normalization run (EPIC's Manifest data.frame alone
## is ~866k rows).
.methyl_anno_cache <- new.env(parent = emptyenv())

## dbSNP build shared by both the 450K and EPIC (v1) annotation packages -
## see METHYL_ANNOTATION_PACKAGES; there's no need to chase the newest
## build (each package also exports 132/135/.../147CommonSingle variants)
## since this is only used for a probe-overlaps-any-known-SNP flag, not an
## exact allele-frequency filter.
METHYL_ANNOTATION_SNP_OBJECT <- "SNPs.147CommonSingle"

## Returns list(ok, anno, reason) - `anno` a data.frame keyed by probe ID
## (rownames) with chr/pos/Type/Probe_rs/CpG_rs/SBE_rs columns when ok.
##
## Deliberately does NOT go through minfi::getAnnotation() on the
## packages' own IlluminaMethylationAnnotation wrapper object - loading
## that wrapper via get(pkg, envir = asNamespace(pkg)) triggers an internal
## Biobase::updateObject() step that fails unless the package is
## library()-attached (not just namespace-loaded), and attaching it pulls
## in Biostrings, which masks base::strsplit() for the entire R session -
## a real risk to every OTHER module in this app that calls strsplit()
## unqualified. Reading the package's own underlying Locations/Manifest/
## SNPs.*CommonSingle data objects directly (which is genuinely all
## getAnnotation() assembles internally anyway) avoids both problems.
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
    ## One representative gene symbol per probe (first token of the
    ## semicolon-separated UCSC_RefGene_Name list, which repeats once per
    ## overlapping transcript) - "a single nearest/overlapping RefGene
    ## symbol per probe", the same simplification this project's own
    ## methylomics pipeline uses (ChAMPdata::probe.features$gene).
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

## Standard Illumina probe-ID prefixes: "cg" = CpG probe, "ch" = non-CpG
## (CpH) probe, "rs"/"nv" = built-in SNP genotyping control probes. Works
## from the probe ID alone, no manifest needed - but only means something
## for an actual Illumina array (see METHYL_ARRAY_TYPES_ILLUMINA above).
methyl_probe_is_cpg <- function(probe_ids) {
  grepl("^cg", probe_ids, ignore.case = TRUE)
}
