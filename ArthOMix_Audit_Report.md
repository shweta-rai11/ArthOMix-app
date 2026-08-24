# ArthOMix Shiny App — Hostile Scientific & Engineering Audit

**Scope:** `shiny_app/` only, as instructed. Any dependency on the sibling `Research_*`/`methylomics/` pipeline folders is noted as an OUTSIDE-SCOPE DEPENDENCY where relevant, but was checked read-only where necessary to verify a claim the app makes about its own precomputed data.

**Method:** Nine independent deep-read audit passes were run in parallel, each tracing actual UI→server→computation→output code paths with file:line citations, cross-checking numeric claims against real bundled data/pipeline output files where possible. Eight of nine passes completed and are the evidentiary basis for this report. **One pass — covering `mod_mr.R` (1055 ln), `mod_coloc.R` (284 ln), `mod_biomarkercard.R` (1082 ln), and `mod_nomogram.R` (660 ln), all in `R/transcriptomics/`) — was terminated before completion at the user's explicit instruction (15-minute deadline) and never returned findings.** Per the audit's own rule to distinguish verified / inferred / not implemented / cannot verify, these four files are marked **CANNOT VERIFY** throughout this report and are excluded from scored credit rather than guessed at.

---

# 1. Executive Summary

This is a large (~46,400-line), hand-rolled, non-golem R Shiny application implementing four omics modules (Transcriptomics, Methylomics, Cross-Omics, Multi-Omics) plus a local-LLM chat assistant. The dominant finding, repeated across nearly every audited file, is that **this codebase does not fabricate data.** Across roughly 43,000 lines directly audited, no instance was found of a chart, table, or statistic populated by `rnorm()`/`runif()`/hardcoded example values presented as real output. Every method claimed (limma, DESeq2, WGCNA, TwoSampleMR, coloc, DMRcate, EpiDISH, IOBR/CIBERSORT, MOFA2, DIABLO, SNF) is a genuine call into the real package, correctly gated on its real input requirements, and in several cases independently verified bit-for-bit against real precomputed output files on disk. This is unusual and should be stated plainly: **for an app of this scope, the "is any of this fake" question is answered no** for everything actually audited.

That said, this is not a clean bill of health. Three classes of real defect were found:

1. **One CRITICAL, unambiguous scientific-validity gap**: the live (user-run) Methylomics DMP/DMR engines omit the bacon/SVA genomic-inflation correction that the app's own bundled reference pipeline *proves is necessary on this exact dataset* (uncorrected λ up to 12.1), and surface no inflation diagnostic to warn a user. A researcher running their own group comparison on this platform has no way to know their p-values may be severely miscalibrated.
2. **One CRITICAL, reproduced-live engineering bug**: the Multi-Omics "Integration → Performance" chart silently fails for 2 of 6 analysis cells (Female/Male-Adalimumab) because of a data-provenance metadata error (`has_snf=TRUE` asserted for cells the underlying pipeline never actually computed) combined with an R `paste0()` zero-length-recycling quirk, both masked by a blanket `tryCatch(...,error=function(e) NULL)` that converts a real crash into an indistinguishable-from-empty "No data" panel. This is confirmed to be **the exact bug the user reported** ("some figures are not being generated properly").
3. **A cluster of HIGH-severity but narrower issues**: a dead CIBERSORT permutation slider in Deconvolution that silently does nothing; an NA-handling divergence between two copy-pasted preprocessing pipelines that will hard-crash Batch Correction on any preloaded dataset with residual missingness; undisclosed feature-selection leakage in the Methylomics Diagnostic Classifier (internal-test AUC is optimistically biased, though the external-test AUC remains valid); and a total absence of dependency pinning (`renv.lock`/`DESCRIPTION`) or automated tests anywhere in the app, undermining its own stated "reproducible biomarker discovery" value proposition.

**Would I trust this application for exploratory research? YES.** The underlying statistical engines are real and, where checked, correct. A researcher using the preloaded/default panels (which is most of the app) is getting genuine analysis.

**Would I trust it for publication-quality analysis? WITH CONDITIONS.** Only after (a) the Methylomics live-engine inflation gap is fixed or clearly banner-warned, (b) the Multi-Omics Adalimumab Performance bug is fixed, and (c) a user independently verifies any live-computed (non-preloaded) result rather than trusting it blind — the app's own precomputed/default panels are more trustworthy than its live-compute panels precisely because the precomputed ones were built by an offline pipeline the live engines don't fully replicate the safeguards of.

**Would I trust it for clinical research? NO**, not as currently built — there is no dependency pinning, no test suite, no audit trail beyond in-code comments, and at least one confirmed silently-wrong output path. Clinical-grade use requires reproducibility infrastructure this app does not have.

**Would I deploy it publicly? WITH CONDITIONS** — no security vulnerabilities (injection, secrets, unsafe eval) were found, but the 200MB blanket upload cap with no rate limiting, unbounded server-side cache growth, and single-threaded WGCNA/MOFA2 compute blocking the whole process for all concurrent users are real availability risks for a public multi-user deployment, not a single-analyst desktop use case.

---

# 2. Application Architecture

Verified structure (`shiny_app/`):

```
shiny_app/
├── global.R (2096 ln)         — shared state, ~25 package loads, ~60+ helper fns, .arthomix_cache
├── server.R (534 ln)          — pure routing/glue; 4 omics families' reactiveValues + module instantiation
├── ui.R (1668 ln)             — navbarPage shell, home page, static JS canvas animations
├── data/                      — one 8.8KB reference file only (cytoBandIdeo_hg19.txt.gz), not a data root
├── www/                       — custom.css, menuhex.css
└── R/
    ├── 0_load_omics_modules.R — sources every R/{transcriptomics,methylomics,crossomics,multiomics}/*.R
    ├── submodules_registry.R  — TX_MODULES/MX_MODULES/CX_MODULES/MULTI_MODULES config lists
    ├── ui_shell.R              — header/sidebar/pipeline-timeline UI
    ├── transcriptomics/ (19 files, ~18,000 ln)
    ├── methylomics/     (19 files, ~17,700 ln)
    ├── crossomics/       (9 files, ~3,080 ln)
    └── multiomics/       (14 files, ~3,010 ln — UNTRACKED in git as of audit)
```

**Framework**: hand-rolled classic modular Shiny (`ui.R`/`server.R`/`global.R` + sourced `R/` modules) — not golem, not rhino. Loading relies on Shiny's `loadSupport()` non-recursive auto-source plus an explicit recursive loader (`0_load_omics_modules.R`) for the four omics subfolders, with a correctly-understood, documented environment-scoping subtlety (`submodules_registry.R:2-6, 103-108`) around why `build_assistant_context()` can't live in `global.R`. This is genuine, non-cargo-culted understanding of Shiny's loader internals.

**Dependency management**: **none.** No `app.R`, no `DESCRIPTION`, no `renv.lock`, no `.Rprofile`. ~25 analysis packages (`limma`, `WGCNA`, `glmnet`, `caret`, `MendelianRandomization`, `coloc`, `clusterProfiler`, `IOBR`, `MOFA2`, `DMRcate`, `EpiDISH`, etc.) are loaded with zero version pinning.

**Tests**: **none.** No `testthat`/`shinytest`/`shinytest2` files anywhere in `shiny_app/`. Every correctness claim in this audit was established by direct code reading and, in several sub-audits, by independently re-running snippets against real bundled data — never by an existing test suite, because none exists.

**Reactive architecture**: `server.R` itself is a clean, side-effect-only routing layer (31 `observeEvent`/`observe` calls, zero `reactive()`/`eventReactive()`) — all actual computation is correctly delegated to the 39 per-submodule server functions, all of which are instantiated eagerly at session start regardless of tab visibility (a documented, intentional but real performance tradeoff). A ~240-line block (four-times-copy-pasted "sub-module toggle/count/search" logic) is the one significant DRY violation at the shell layer. One latent trap: the `dataset$source` results-reset observer (`server.R:19-21`) fires on any *reassignment*, not just value change, and could silently wipe cached results on a benign reload.

**Global state**: one process-wide cache, `.arthomix_cache` (`global.R:211`), correctly scoped to deterministic read-only reference data (never per-user state) — the right pattern — but with no eviction/TTL (unbounded memory growth risk) and no concurrency guard on disk-cache writes (real but low-probability race given an active `future::multisession` pool).

---

# 3. Transcriptomics Audit

**Coverage**: Data/QC/Preprocessing/Overview/FeatureSelection (fully audited); DGE/Enrichment/WGCNA/Interaction/Diagnostic (fully audited); Cross-Tissue/Cross-Ancestry/Deconvolution/Candidates/ArthOChat (fully audited). **`mod_mr.R`, `mod_coloc.R`, `mod_biomarkercard.R`, `mod_nomogram.R` — CANNOT VERIFY (audit terminated before this slice completed).**

## Data / QC / Preprocessing — key findings
- Normalization is real and correctly type-gated: `limma::normalizeBetweenArrays(method="quantile")` for microarray/log-scale, `edgeR::filterByExpr→calcNormFactors(TMM)→cpm(log=TRUE)` for raw counts, with a `validate()` guard rejecting negative values before TMM — quantile norm is never applied to raw counts nor TMM to already-logged data ([mod_preprocessing.R:1362-1463](shiny_app/R/transcriptomics/mod_preprocessing.R)). ComBat-seq → recomputed TMM factors, matching Zhang/Parmigiani/Johnson 2020 convention.
- **HIGH bug**: `mod_preprocessing.R:1436-1438` — `stats::var()`/`stats::quantile()` in the Batch Correction gene-variance filter lack `na.rm=TRUE`, one line below a correctly-`na.rm=TRUE`'d `rowMeans()` call. `stats::quantile()` hard-errors on any NA. The parallel "Preloaded Data" quick-load box (`mod_preprocessing.R:139-180`) — unlike its sibling per-source pipeline — has **no missing-value handling at all**, so any residual NA in a preloaded dataset reaching Batch Correction crashes with an unhandled R error, not a clean message. Two copy-pasted pipelines that diverged.
- **HIGH gap**: no deduplication of feature IDs anywhere in the real pipeline (only diagnostically counted in the separate EDA tab) — duplicate gene/probe rownames silently keep only the first match on any row-name-keyed merge.
- **MEDIUM**: no transposed-matrix (samples-in-rows) detection anywhere; Feature Selection's SVM/LASSO/RF fits have no pre-fit zero-variance check (inconsistent with the app's own `pca_of()` pattern elsewhere).
- Feature selection math (LASSO via `glmnet::cv.glmnet`, RF via `randomForest`, SVM-RFE via a correctly-derived squared-primal-weight ranking citing Guyon et al. 2002) is all correctly implemented, seeded (1234), and reproducible. The "fast path" precomputed FS results were confirmed present on disk and gated behind a real 4-condition eligibility check, not a blind shortcut.
- **Environment caveat**: raw GEO source files are absent from `data/raw/` in this checkout, blocking per-source preprocessing/live-merge/per-source-QC for all four bundled GSEs (only the default merged/batch-corrected cohort works). May be an environment artifact rather than a shipped defect — needs confirmation against the actual deployment.

## Statistics (DGE / Enrichment / WGCNA / Interaction / Diagnostic) — key findings
- **DGE**: limma (`arrayWeights`→`lmFit`→`makeContrasts`→`eBayes`) and DESeq2 paths both correctly specified, BH-correction never overridden, raw-count guard prevents DESeq2-on-log-data. Volcano plot's dashed threshold line is data-derived from the live fit object, not decorative. |log2FC|>0.1 default threshold **verified to match the project's own published methods chapter exactly** (`Chapter_2_subchapter2_sexstratified.md:102-106`).
- **Enrichment**: universe = the dataset's own mapped genes (correct — avoids universe-inflation), BH explicit in every `enrichGO`/`enrichKEGG`/`enrichPathway` call.
- **WGCNA**: real `WGCNA::` calls throughout; β=12 default, disease-module cutoff, and hub-gene rule **all verified to exactly match the project's own methods doc**; the precomputed hub table's claimed "217 yellow/144 brown" hub genes **independently confirmed via direct `awk` count against the real output CSV** — not fabricated.
- **Diagnostic module (`mod_diagnostic.R`, 2983 ln)** — the leakage check that matters most: train/test standardization uses Train-only mean/SD, CV re-standardizes per-fold, Youden threshold locked from Train never re-optimized on Test, and the "Advanced ML" nested-CV pipeline refits every supervised step inside each outer fold only — **textbook-correct, leakage-free nested CV**, verified line-by-line.
- **MEDIUM**: a fully-implemented Platt-scaling calibration function (`diag_adv_platt_calibrate`) is never wired into the Calibration tab — the displayed reliability diagram is uncalibrated raw probabilities. A gene-level "hub" screen uses raw uncorrected Wilcoxon p-values (explicitly labeled as replicating an external paper's rule, not a corpus-wide FDR claim, but a user could misread it).
- **HIGH performance risk**: WGCNA's "use all genes" default can take tens of minutes and blocks the single-threaded Shiny process for every concurrent user during that time.

## Cross-Tissue / Cross-Ancestry / Deconvolution / Candidates / ArthOChat — key findings
- Cross-Tissue and Cross-Ancestry validation are both genuine, correctly-unpaired, real-cohort comparisons (verified against real bundled `.rds`/`.csv` files — e.g., MR35_crossancestry files' 32/25-gene counts match the module's own header claims). MR reuse (`estimate_mr_set()`) is verified as real shared machinery, bit-for-bit checked against reference numbers, not superficial.
- **MEDIUM**: Cross-Ancestry's "nominal p<threshold" caveat exists only in code comments, never surfaced in the UI slider labels, despite testing ~25-32 genes per sex for "biomarker" status with no FDR correction. An allele-frequency-transferability diagnostic (`eafgap_bbj`) is computed and displayed but never actually used in the classification logic.
- **HIGH — dead UI control**: Deconvolution's CIBERSORT permutation slider (`cib_perm`) has **zero effect on any displayed output** — the only columns it would influence (P-value/Correlation/RMSE fit diagnostics) are explicitly stripped before display, meaning a user has no in-app way to detect a poor tissue/platform match producing plausible-but-unreliable fractions.
- Candidates module: correct hypergeometric enrichment math, no randomization, fully real wiring.
- ArthOChat: genuine live local-Ollama LLM integration (not a stub), no hardcoded secrets found anywhere in a repo-wide scan, responses grounded in real live session state and real document/API sources. **MEDIUM**: the live per-turn chat call lacks a `tryCatch`, inconsistent with the careful degradation used everywhere else in the same file.

**Transcriptomics score: 75/100** (see §19 for rubric breakdown). Reflects strong verified statistics and a clean track record on fabrication, offset by two HIGH-severity code bugs, one dead UI control, and an unverified ~15% of the module's feature surface (MR/coloc/biomarker-card/nomogram).

---

# 4. Methylomics Audit

**Coverage**: all 19 files fully audited across three passes (Data/QC/Normalization/Cell-type; DMP/DMR/FeatureSelection/Candidates; MR/Coloc/Diagnostic/BiomarkerCard/WGCNA).

## Headline: methylation is NOT treated as expression
Beta↔M-value conversion is mathematically verified correct and numerically stable (`M = log2(beta/(1-beta))` with proper `[1e-4, 1-1e-4]` clamping; exact algebraic inverse used everywhere it round-trips — ComBat, RUVm, exports). M-values are used for statistical testing, beta retained separately for interpretable effect size — the textbook Du et al. 2010 approach, applied identically and correctly in both DMP and DMR live engines. Rank-based tests are deliberately run on beta (statistically valid, since rank tests are invariant to the monotonic logit transform — not an oversight). This is the opposite of a copy-pasted, unadapted transcriptomics routine.

## CRITICAL finding
The app's own bundled reference pipeline (`methylomics/script03_dmp_sva_sexstratified/`) documents severe genomic inflation on this exact dataset without SVA/bacon correction — uncorrected λ up to **12.1** (female) / **4.1** (male) — bad enough that **zero CpGs reached genome-wide significance** before correction. The app's default/"SVA" tab exists specifically to demonstrate and correct this. **The live, user-facing DMP and DMR engines run plain `limma`/`DMRcate` with no bacon/SVA correction and surface no inflation diagnostic (λ, QQ plot) anywhere in the UI.** A user running their own group comparison — exactly the scenario this app's own pipeline proved is prone to severe inflation — gets no warning signal. This directly compounds into DMR region-calling, which seeds candidate regions from raw (not bacon-corrected) p-values in the live engine, unlike the reference pipeline it otherwise mirrors exactly.

## Other findings
- **HIGH — undisclosed leakage** (Diagnostic Classifier): the CpG feature panel was selected using label information from the *same* 689-sample GSE42861 cohort later re-split into the module's own internal train/test sets — verified against the pipeline's own methods doc, which is itself candid about this but the live Shiny tool does not carry the caveat forward. Internal-test AUC is optimistically biased for two of three feature sources; the external cohort (GSE111942) remains a valid held-out estimate and is correctly available in the UI, just not visually distinguished from the biased numbers.
- **MEDIUM**: "t-test" univariate mode in Feature Selection is always limma-moderated (eBayes) regardless of user selection — mislabeled, not a classical unmoderated test. All DMP/DMR volcano plots (4 call sites) mislabel the y-axis as `-log10(p)` when the value plotted is always `-log10(FDR)` — cosmetic but present on every such plot; the Candidates module's own volcano gets this right, proving it's a fixable inconsistency within the same codebase. A dead QC checkbox (`ct_qc_beta_range`) is checked-by-default but read nowhere. PCA/MDS/outlier detection run on raw beta rather than M-values, inconsistent with the app's own stated rationale for M-values elsewhere.
- **Strong positives, stated explicitly**: real Illumina Bioconductor manifests joined by CpG ID (never by position); coloc/MR modules genuinely real (`coloc::coloc.abf`, `TwoSampleMR::mr()`, F-stat/heterogeneity/pleiotropy diagnostics all real and correctly gated); cell-type deconvolution uses real EpiDISH reference panels with honest "unavailable" labeling for methods not actually implemented; DMR calling is genuine kernel-smoothed DMRcate verified to match the reference pipeline's parameters exactly; no bundled cross-reactive-probe blacklist is a disclosed limitation, not fabricated data ("fabricating one would violate this project's evidence-based-methods requirement" — an unusually honest in-code stance).

**Methylomics score: 74/100 — capped due to a CRITICAL scientific-validity gap** (undisclosed inflation risk in the live DMP/DMR engine). Without that one finding this module would score in the mid-to-high 80s; the cap is deliberate per the audit's scoring rules.

---

# 5. Cross-Omics Audit

**Coverage**: all 9 files fully audited. This is the strongest-performing module in the app.

- "Cross-omics" is verified to be three distinct sub-features (Biomarker Convergence, Expression×Methylation Integration, Cross-Omics MR), not one integration engine — characterized precisely rather than assumed.
- **Headline verified claim**: the commit message "integrate, not recompute" is **true in code**, not just in the commit message — Biomarker Convergence and Cross-Omics MR both load real pre-existing pipeline CSVs and only relabel already-computed significance flags live; the Integration sub-module genuinely reuses live Transcriptomics DGE results and the Methylomics pipeline's own SVA/bacon-corrected DMP table, with an honestly-disclosed (not silently-broken) limitation that live Methylomics session results can't yet be reused because they only carry summary counts.
- **Multiple-testing correction is present and correctly applied everywhere tested** — a genuinely positive finding relative to the audit's top stated concern for this kind of module.
- Cross-Omics MR is verified genuine single-instrument two-sample MR (Wald ratio, GoDMC cis-mQTL → Ishigaki 2022 RA GWAS), matched exactly against the actual pipeline script, correctly caveated for single-instrument limitations.
- Causal-vs-association language is **unusually disciplined** — every relevant output explicitly states association, not causation, in user-facing text (not just code comments).
- **HIGH**: CpG-to-gene annotation silently keeps only the first gene when a probe maps to multiple genes (`UCSC_RefGene_Name` semicolon-split, first token only) — disclosed only in a source comment, never in the UI/provenance text, so users can't know some genes are systematically undercounted for methylation evidence.
- **MEDIUM**: a user-facing "BH-FDR recomputed" mQTL-MR significance option can operate on as few as 5 p-values (verified: female table has only 5 non-NA `mQTL_MR_pval` of 43 genes) with no small-n warning; a "Strong candidate" correlation threshold is hardcoded inconsistently with the otherwise fully user-configurable FDR method elsewhere in the same module.

**Cross-Omics score: 82/100.**

---

# 6. Multi-Omics Audit

**Coverage**: all 14 files fully audited, including reproducing the reported bug live against real pipeline data (`/Users/swetarai/ArthOMix/Research_05_multiomics_sexstratified`).

## Integration method — real, not concatenation
Two structurally different halves, both verified genuine:
- **Precomputed tabs** (Overview/Integration/Stratification/Biomarker/Concordance/Pathway): real DIABLO (`mixOmics::block.splsda`, confirmed real `prop_expl_var` slot with plausible values) and real SNF (`SNFtool::SNF` outputs read as CSVs, traced to the actual generating pipeline scripts) — but the Shiny app performs **zero** integration computation itself here; it is a read-only browser over an offline pipeline's already-computed fits, and says so honestly in its own module descriptions.
- **Live Analysis tab**: genuine, live-computed **MOFA2** (`create_mofa`→`prepare_mofa`→`run_mofa`→`get_variance_explained`/`get_factors`/`get_weights`, real Python backend confirmed available in this deployment) — real latent-factor decomposition, not stacked-matrix PCA relabeled.
- Factor interpretation is correctly hedged everywhere ("not automatically disease factors, biomarkers, or causal" — repeated at every relevant plot caption).

## CRITICAL — root cause of the reported "figures not generating properly" bug, pinpointed and reproduced
`multiomics_helpers.R:52-59` hardcodes `has_snf=TRUE` for **all four** drug×sex cells, but the real `Table22_snf_integration_performance.csv` on disk contains **only 4 rows, all Etanercept** — the SNF classifier-performance benchmark was simply never computed for Adalimumab by the source pipeline, even though the app's registry metadata claims it was. This causes `multi_filter_cell()` to return a 0-row (not `NULL`) data frame, which then hits an R `paste0()` zero-length-recycling quirk (`paste0("SNF (", character(0), ")")` returns a length-1 string, not length-0) inside a `data.frame()` construction, throwing `"arguments imply differing number of rows: 1, 0"`. This error is **silently swallowed** by a blanket `tryCatch(...,error=function(e) NULL)` wrapping every plot function in the module, producing a generic "No data" panel indistinguishable from a legitimate empty state. **Confirmed live**: this reproduces for exactly "Female-Adalimumab" and "Male-Adalimumab" (2 of 6 cells, 1/3 of the dropdown) and matches the user-reported symptom exactly — the DIABLO performance table below the broken chart shows real data for the same cell, because it doesn't touch the buggy code path, so the user sees "the table has data but the chart above it doesn't."

- **HIGH**: the PNG download handler for this same chart has no `tryCatch` at all — clicking "Download plot" for the broken cells throws an unhandled server error rather than any graceful message.
- **HIGH (systemic)**: the module-wide pattern of swallowing every plot-function error into one generic "No data" message means *any* future bug in ~20 figure functions will look identical to a legitimate empty state, with no logging to tell them apart — this is the actual root cause of why the bug went unnoticed rather than being caught immediately.
- **MEDIUM**: no MOFA2 convergence diagnostic (ELBO trace, iteration count) is ever shown post-training — a barely-converged model looks identical in the UI to a well-converged one.
- Sample alignment across omics layers is genuinely excellent — real order-independent set intersection, row-name-based (never positional) indexing, duplicate IDs surfaced not silently merged.
- No fabricated/randomly-generated placeholder data was found driving any of the ~20 figures in this module.

**Multi-Omics score: 70/100 — capped due to a CRITICAL, reproduced, confirmed-live defect** that matches the exact bug the user reported. The underlying science (real DIABLO/SNF/MOFA2, excellent sample alignment, honest factor-interpretation language) is genuinely strong; the cap reflects the confirmed live-breaking bug and the systemic error-masking pattern that let it hide.

---

# 7. Shared Architecture Audit

See §2 for structural detail. Summary of audit-relevant findings not already covered:

- **Security**: no `system()`/`system2()`/`shell()`, no `eval(parse(...))`, no hardcoded credentials/API keys anywhere in the shared layer or in a repo-wide grep, no unescaped user-input `HTML()` injection at this layer (all 4 `HTML()` calls wrap static JS constants). One hardcoded, machine-specific absolute path fallback (`METH_RAW_DATA_ROOT`, `global.R:392-395`) leaks a specific developer's local filesystem layout into source — fails safe on other machines but is a portability/hygiene smell.
- **Error handling**: generally disciplined — narrow `tryCatch` blocks around genuinely-optional steps, paired with real fallback/messaging. The one clear exception: a weighted-mode MR estimate is silently dropped on error with no user-facing indication (`global.R:1355`), in a helper shared by both `mod_mr.R` and `mod_crossancestry.R`.
- **Code quality**: ~240 duplicated lines (four-times-copy-pasted sub-module toggle/search block); stale documentation in `submodules_registry.R:60-65` describing Cross-Omics sub-modules as "unbuilt placeholders" when they are demonstrably (file size, git log) substantially built — a real documentation-drift risk for future maintainers.
- **Performance**: eager instantiation of all 39 submodule servers on every session regardless of tab use; genuinely well-engineered mitigations already in place elsewhere (two-tier digest-keyed disk+memory cache for the "slowest step in the app," `future::multisession` for the multi-GB methylation matrix load that used to freeze the whole app).

---

# 8. AI-Code / Human-Code Audit

Every one of the nine audit passes independently flagged the same authorship signature, and it is worth stating as a single cross-cutting finding rather than repeating nine times: **this codebase reads as heavily AI-assisted (or AI-written under close human direction — one source methods document literally self-discloses "rewritten 2026-08-03 ... at the user's request"), and the evidence is unusually positive rather than negative.**

Typical negative AI-generated patterns the audit specifically looked for — and did **not** find at meaningful scale:
- No `rnorm`/`runif`/`sample()`-fabricated data driving any chart or table (checked in all nine passes, explicitly grepped for in several).
- No hardcoded example datasets standing in for real results.
- No "coming soon" functionality presented as live (the few genuinely-unfinished features, e.g. some Biomarker Card sections, are explicitly labeled "Not yet connected" and never rendered as if live).
- No generic "# Process data" / "# Generate plot" comments with no scientific content — comments are almost universally specific, technical, and (where independently checked) accurate.

Patterns that **were** found, all flagged per-module above and consolidated here:
- **Recurring "audited on \<date\>" self-report comments** (e.g. `global.R:1524-1586`, headed "audited 2026-08-12," with a structured "gap found and addressed" narrative). These are genuinely substantive, not filler — but they are **self-asserted developer/AI claims, not independent review**, and should be weighted as such by anyone relying on them.
- **Near-verbatim repeated anti-pattern-avoidance phrases** across separate Multi-Omics files ("never render a fake placeholder," "never a fake/placeholder value") — strongly suggests a shared AI-authored spec was followed consistently, and (per §6) the claim largely holds up under direct verification.
- **Duplicated/templated boilerplate** where a pattern was copy-pasted across the four omics families rather than refactored after the second or third repetition (shared-architecture sub-module toggle block; Multi-Omics submodule server signatures carrying an unused `multi_dataset` parameter) — consistent with iterative "add module X following the same pattern as module Y" prompting without a cleanup pass.
- **Dense, uniformly-styled rationale comments** far exceeding typical human commenting discipline in volume and consistency across ~46,000 lines — a stylistic tell, not a defect.
- **The one place this pattern caused real harm**: the blanket `tryCatch(...,error=function(e) NULL)` wrapping every plot function in Multi-Omics (§6) is exactly the kind of "defensive programming that hides the underlying bug" pattern the audit brief warned about, and it is the direct cause of the reported figures bug going undiagnosed.

**Verdict for this category**: mostly "Normal human code" / "Likely AI-assisted, stylistic only" with a small number of genuine "Suspicious pattern" findings, concentrated in defensive-error-masking rather than fabrication. No instance was found across ~43,000 audited lines that would be categorized "Likely AI-assisted" **and** "fabricated/incorrect" simultaneously — the AI-authorship signal and the correctness signal are largely independent in this codebase, which is itself worth noting since they are often conflated.

---

# 9. Statistical Audit

Cross-module synthesis of the statistical findings above:

- **Multiple-testing correction**: applied correctly and consistently in DGE (BH via limma/DESeq2 default), Enrichment (BH explicit), WGCNA module-trait correlation, Cross-Omics (BH/Bonferroni user-selectable, applied at the full-set level not ad hoc), Methylomics DMP/DMR (BH via limma/DMRcate's Stouffer). **Exceptions found**: the transcriptomics "hub gene" screen (`AUC≥0.85, P<0.05`) uses raw Wilcoxon p-values, explicitly labeled as replicating an external paper's specific rule rather than claiming corpus-wide significance (defensible but risks misreading); Cross-Ancestry's "validated biomarker" classification uses nominal, not FDR-corrected, p-values across ~25-32 genes with the caveat only in code comments, not the UI.
- **Data leakage**: rigorously checked in every ML pipeline found (Feature Selection, transcriptomics Diagnostic Classifier's 4-model and nested-CV engines, Cross-Tissue's CV evaluation) — all verified leakage-free with fold-local re-standardization and Train-only threshold locking. **One confirmed leak**: Methylomics Diagnostic Classifier's CpG panel was selected using the same cohort later re-split into train/test, inflating internal-test AUC (external-test AUC remains valid but isn't visually distinguished in the UI).
- **Confounding/batch effects**: genuinely modeled where relevant (ComBat/ComBat-seq/SVA/RUVm real and correctly sequenced relative to normalization; covariates genuinely enter DGE/DMP design matrices, not merely offered and ignored, with rank-deficiency guards).
- **The one systemic statistical gap**: genomic-inflation control (Methylomics live engine, §4) — the single most consequential statistical-validity finding in the whole audit, because it affects the primary "run your own comparison" workflow rather than a secondary feature.
- **Sample size / edge cases**: consistently guarded with real `validate(need(...))` minimums across essentially every module (not decorative — traced to actually gate the reactive and produce a real message).

---

# 10. Bioinformatics Audit

- **Identifier handling**: broadly disciplined — CpG-to-gene, gene-symbol harmonization, and cross-omics ID matching are conservative (ambiguous cases flagged, not guessed), with one real exception (Cross-Omics's undisclosed first-gene-only collapse for multi-gene-annotated CpGs, §5).
- **Platform-appropriate methods**: RNA-seq counts get TMM/CPM, microarray/log data gets quantile normalization, methylation gets platform-appropriate Bioconductor methods (Noob/Funnorm/SWAN/Dasen/BMIQ/PBC), each correctly gated on its real input requirement rather than a one-size-fits-all "normalize" button.
- **Reference data**: real Illumina manifests, real EpiDISH cell-type reference panels, real LM22/MCP-counter deconvolution references, real GO/KEGG/Reactome backgrounds scoped to the dataset's own genes (not the whole database) — no ad hoc or fabricated reference data found anywhere audited.
- **Methylation-specific correctness**: the single strongest bioinformatics finding in the audit is that methylation is *not* treated like expression anywhere audited — beta/M-value separation is textbook-correct and universal (§4).
- **The gap**: genomic inflation control in live methylation analysis (§4, §9) is the one place platform-appropriate best practice (documented and proven necessary by the app's own pipeline) is not carried into the live tool.

---

# 11. Visualization Audit

Across all audited modules, visualizations were checked for real-data provenance, axis/label correctness, and plot/table/download consistency. Findings:

- **No fabricated visualization data found anywhere in ~43,000 audited lines.** Every chart traces to a live reactive built from real computation; downloads consistently match on-screen data (verified explicitly in Transcriptomics, Methylomics, Cross-Omics).
- **Confirmed broken (CRITICAL)**: Multi-Omics Integration → Performance chart for 2 of 6 cells (§6).
- **Mislabeled but not wrong (MEDIUM)**: all four Methylomics DMP/DMR volcano plots label the y-axis `-log10(p)` when the value plotted is always `-log10(FDR)` — data and significance coloring are internally consistent, only the label is wrong, and the Candidates module's own volcano proves this is a fixable, isolated inconsistency.
- **Methodological note, not a bug**: Methylomics QC/Normalization PCA/MDS/outlier detection run on raw beta values rather than M-values, inconsistent with the app's own stated rationale for M-values elsewhere (defensible for outlier detection, but worth exposing as a choice).
- **Legitimate randomness use, correctly scoped**: the only `sample()`/`runif()` calls found across the audit are performance-motivated subsampling of real data for plotting (e.g. EDA tab downsampling to 200k points, WGCNA TOM heatmap gene subsample) — always disclosed, sample size reported alongside, never used to fabricate a value.

---

# 12. Table Audit

Every table checked traces to the same reactive driving its adjacent plot and its download handler — no divergence path was found in any audited module where a table could show different numbers than its chart. Column naming and statistical-value formatting are generally clear (gene/CpG identifiers, FDR/p-value/effect-size columns consistently present). The one table-adjacent finding worth flagging here specifically: in Multi-Omics, the DIABLO performance **table** for the broken Adalimumab cells correctly shows real data even while the **chart** above it is silently blank — meaning the table itself is fine, but its presence right next to a broken chart is precisely what makes the bug read as a rendering glitch rather than a missing-data state (§6).

---

# 13. Security Audit

No CRITICAL or HIGH security findings anywhere in the audited scope.

- No `system()`/`system2()`/`shell()`/`eval(parse(...))` calls found in any audited file.
- No hardcoded credentials, API keys, or tokens found in a repo-wide scan; all external credentials (OpenGWAS JWT, none needed for local Ollama) are sourced via `Sys.getenv()` with graceful degradation when unset.
- No unescaped user-input `HTML()`/XSS injection found at the shared-architecture layer (all static JS constants); ArthOChat's LLM-output/PubMed-title rendering was not independently re-verified for HTML escaping in this audit round and is worth a targeted follow-up given it renders both LLM output and external API text.
- **MEDIUM**: blanket 200MB upload cap applied uniformly app-wide with no per-endpoint differentiation or rate limiting — a generic DoS-surface consideration for a public deployment, not itself a vulnerability.
- **MEDIUM (hygiene, not exploit)**: one hardcoded, machine-specific absolute filesystem path fallback baked into `global.R`.
- ArthOChat's prompt construction correctly passes user input as a distinct message parameter (not string-concatenated), and external tool calls receive only narrow LLM-generated query strings rather than the full session-context object — no structural bulk-exfiltration path, and the data involved is de-identified research data, not PHI.

---

# 14. Performance Audit

| Bottleneck | Location | Severity |
|---|---|---|
| WGCNA "use all genes" (~15,763 genes) live compute can take tens of minutes and blocks the single-threaded Shiny process for all concurrent users | `mod_wgcna.R:975-1002` | **High** |
| `TOMsimilarityFromExpr` full-matrix TOM computation, O(n²) memory, no gene-count safety cap at "all genes" on an uploaded dataset | `mod_wgcna.R:1098-1130` | Medium-High |
| Eager instantiation of all 39 submodule servers on every session regardless of tab use | `server.R:135,95,110-116,127` | Medium (real cost depends on per-module setup work, not fully characterized) |
| Unbounded `.arthomix_cache` growth (no TTL/eviction) over a long-lived deployment | `global.R:211` | Medium |
| No concurrency guard on shared disk-cache writes, combined with an active 2-worker `future::multisession` pool | `global.R:252-267,571-586` | Low-Medium (real but low-probability race) |
| Blanket 200MB upload cap with no per-endpoint tuning or rate limiting | `global.R:9` | Low |
| Process-global WGCNA thread pool contended across concurrent users' analyses | `global.R:85-87` | Low-Medium |

Countervailing positives: a genuinely well-engineered two-tier digest-keyed cache for WGCNA's "single slowest step," and an explicit `future::multisession` fix for what used to be a multi-GB blocking methylation-matrix load — real, hands-on performance engineering exists in this codebase, it's just incomplete.

---

# 15. Reproducibility Audit

Mixed picture, worth stating precisely because the two halves point in opposite directions:

- **Infrastructure: absent.** No `renv.lock`, no `DESCRIPTION`, no `.Rprofile`, no pinned package versions for ~25 analysis packages, no test suite. This directly undermines the app's own "reproducible biomarker discovery" framing.
- **Within-app numerical reproducibility: strong, where checked.** Fixed seeds (1234) are used consistently across every stochastic step found (LASSO CV, RF tuning, SVM-RFE, bootstrap stability selection, MOFA2 factor count logic). Multiple precomputed-panel claims were **independently re-verified against real output files and matched exactly** (WGCNA hub-gene counts, Cross-Omics/Cross-Ancestry gene-count headers, DIABLO variance-explained values) — this is not asserted reproducibility, it's demonstrated reproducibility for the specific things checked.
- **The gap between "this app's numbers are internally consistent" and "this app's environment is reproducible on another machine"** is the real finding here: without a lockfile, the *numbers* verified above are only guaranteed reproducible on a system with the exact same package versions currently installed, which nothing in the repo pins down.

---

# 16. UX Audit

- Error/edge-case messaging is generally strong — real `validate(need(...))` gates with specific, actionable text rather than generic Shiny stack traces, found consistently across every module.
- **Real UX failure modes found**: the Multi-Omics silently-blank chart (§6) is indistinguishable from a legitimate empty state — a researcher has no way to know whether "no data" means "this cell genuinely has none" or "something crashed." The Methylomics live engine's total absence of an inflation diagnostic (§4) is a transparency failure specifically for the class of user most likely to be misled by it (someone running their own comparison, trusting the output at face value). A dead permutation slider (Deconvolution) and a mislabeled "t-test" option (Methylomics FS) both actively mis-set user expectations about what a control does.
- **Terminology discipline is a genuine strength**: causal-vs-associative language is consistently hedged in user-facing text across MR, coloc, and cross-ancestry modules (verified by direct UI-string grep, not just code comments) — a PhD-level researcher would not be misled about association-vs-causation by this app's actual output text, even where the underlying method (MR) is itself a causal-inference framework.

---

# 17. 🚨 Critical Red Flags

### CRITICAL
1. **Methylomics live DMP/DMR engines have no genomic-inflation correction or diagnostic**, despite the app's own bundled pipeline proving this exact dataset shows severe inflation (λ up to 12.1) without it. Affects any user-run group comparison. (`mod_methyl_dmp.R:504-770`, `mod_methyl_dmr.R:618-949`)
2. **Multi-Omics Integration → Performance chart silently fails for the Female/Male-Adalimumab cells** (confirmed live, matches the reported symptom exactly) due to incorrect `has_snf` registry metadata + a `paste0()` recycling quirk, masked by a blanket error-swallowing `tryCatch`. (`mod_multi_integration.R:120-131`, `multiomics_helpers.R:52-59`)

### HIGH
3. Transcriptomics Batch Correction hard-crashes on any preloaded dataset with residual missingness (`stats::var`/`quantile` missing `na.rm=TRUE`, compounded by a sibling preprocessing path with no NA handling at all). (`mod_preprocessing.R:139-180, 1436-1438`)
4. No feature-ID deduplication anywhere in the real transcriptomics pipeline — duplicate gene/probe rows are silently dropped without warning on merge.
5. Deconvolution's CIBERSORT permutation slider is a dead control with zero effect on output; fit-quality diagnostics that would catch a tissue/platform mismatch are computed then discarded. (`mod_deconvolution.R:92-96,145-146,199`)
6. Methylomics Diagnostic Classifier's CpG panel has undisclosed feature-selection leakage against its own internal test set (external test remains valid but isn't visually distinguished in the UI). (`mod_methyl_diagnostic.R:1036-1058`)
7. Multi-Omics PNG-download handler has no error handling — the same crash as #2 is unhandled inside a `downloadHandler`. (`multiomics_plots.R:33-42`)
8. The module-wide blanket `tryCatch`-to-empty-state pattern in Multi-Omics makes any future bug in ~20 figure functions indistinguishable from legitimate empty data — the systemic reason #2 went undiagnosed. (`multiomics_plots.R:23-27`)
9. No dependency pinning (`renv.lock`/`DESCRIPTION`) and no automated test suite anywhere in `shiny_app/`, for an app loading ~25 analysis packages and marketed on reproducibility.
10. Cross-Omics silently collapses multi-gene-annotated CpGs to their first-listed gene, undisclosed outside a source comment. (`crossomics_integration_helpers.R:602-603`)

### MEDIUM
11. Methylomics "t-test" Feature Selection option is always limma-moderated, never the classical test its label implies.
12. All Methylomics DMP/DMR volcano plots mislabel the y-axis as raw p when it's always FDR.
13. Cross-Ancestry's "nominal p" caveat for biomarker classification exists only in code comments, not the UI.
14. No MOFA2 post-training convergence diagnostic shown to users.
15. A weighted-mode MR estimate is silently dropped on error app-wide with no user-facing indication.

### LOW
16. Stale documentation in `submodules_registry.R` describing built Cross-Omics sub-modules as unbuilt placeholders.
17. ~240 duplicated lines of sub-module toggle/search logic across the four omics families.
18. Hardcoded, machine-specific absolute path fallback in `global.R`.
19. Inconsistent boundary-clamp epsilon values (1e-4 vs 1e-6) for the same beta/M-value edge case in two files.
20. Two fully-implemented-but-never-called functions in the Diagnostic module (Platt calibration, threshold-sweep table).

---

# 18. Feature Scorecard

| Module | Feature | UI | Backend | Scientific correctness | Visualization | Edge cases | Reproducibility | Score | Severity |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| Transcriptomics | Data upload/preprocessing | 8 | 6 | 7 | 8 | 6 | 7 | 7/10 | HIGH (NA crash bug) |
| Transcriptomics | DGE (limma/DESeq2) | 9 | 9 | 9 | 9 | 8 | 9 | 9/10 | — |
| Transcriptomics | Enrichment | 9 | 9 | 9 | 8 | 8 | 8 | 9/10 | — |
| Transcriptomics | WGCNA | 8 | 9 | 9 | 8 | 7 | 9 | 8/10 | HIGH (perf/blocking) |
| Transcriptomics | Feature Selection | 8 | 8 | 9 | 8 | 7 | 9 | 8/10 | MEDIUM |
| Transcriptomics | Diagnostic Classifier | 8 | 9 | 9 | 8 | 8 | 8 | 8/10 | MEDIUM (dead calibration) |
| Transcriptomics | Cross-Tissue/Ancestry | 8 | 8 | 8 | 8 | 8 | 8 | 8/10 | MEDIUM |
| Transcriptomics | Deconvolution | 6 | 8 | 7 | 7 | 5 | 8 | 6/10 | HIGH (dead control) |
| Transcriptomics | Candidates | 9 | 9 | 9 | 8 | 8 | 9 | 9/10 | — |
| Transcriptomics | ArthOChat | 8 | 8 | 8 | 7 | 6 | 7 | 7/10 | MEDIUM |
| Transcriptomics | MR / Coloc / Biomarker Card / Nomogram | — | — | — | — | — | — | **N/A** | **CANNOT VERIFY** |
| Methylomics | Data/QC/Normalization | 9 | 9 | 9 | 8 | 8 | 8 | 9/10 | MEDIUM |
| Methylomics | Cell-type Deconvolution | 8 | 9 | 9 | 8 | 8 | 7 | 8/10 | MEDIUM (dead checkbox) |
| Methylomics | DMP (live) | 8 | 6 | 4 | 6 | 8 | 6 | 6/10 | **CRITICAL** (inflation) |
| Methylomics | DMR (live) | 8 | 6 | 4 | 6 | 8 | 6 | 6/10 | **CRITICAL** (inflation) |
| Methylomics | DMP/DMR (default/SVA) | 8 | 9 | 9 | 7 | 8 | 9 | 8/10 | LOW (axis label) |
| Methylomics | Feature Selection | 8 | 8 | 8 | 8 | 8 | 8 | 8/10 | MEDIUM (mislabel) |
| Methylomics | Candidates | 9 | 9 | 9 | 8 | 8 | 9 | 9/10 | — |
| Methylomics | Methylation MR | 8 | 9 | 9 | 8 | 8 | 8 | 8/10 | — |
| Methylomics | Colocalisation | 8 | 9 | 9 | 8 | 8 | 8 | 8/10 | — |
| Methylomics | Diagnostic Classifier | 7 | 8 | 6 | 7 | 8 | 6 | 7/10 | HIGH (leakage undisclosed) |
| Methylomics | Biomarker Card | 8 | 8 | 8 | 8 | 8 | 8 | 8/10 | — |
| Methylomics | WGCNA (methylation) | 8 | 9 | 9 | 8 | 8 | 8 | 8/10 | — |
| Cross-Omics | Biomarker Convergence | 9 | 9 | 9 | 8 | 8 | 9 | 9/10 | — |
| Cross-Omics | Expr×Meth Integration | 8 | 9 | 8 | 9 | 9 | 8 | 8/10 | HIGH (CpG collapse) |
| Cross-Omics | Cross-Omics MR | 8 | 9 | 8 | 8 | 8 | 8 | 8/10 | MEDIUM (small-n) |
| Multi-Omics | Precomputed (DIABLO/SNF) tabs | 7 | 9 | 9 | 6 | 6 | 8 | 7/10 | **CRITICAL** (broken chart) |
| Multi-Omics | Live MOFA2 | 8 | 8 | 8 | 7 | 6 | 6 | 7/10 | MEDIUM (no convergence diag) |
| Multi-Omics | Concordance/Stratification/Pathway | 8 | 9 | 9 | 8 | 8 | 8 | 8/10 | — |

---

# 19. Scores

## Transcriptomics: **75/100**
Data handling 6/10, Preprocessing 7/10, Statistics 17/20, Visualization 13/15, Biological validity 11/15, UI/UX 7/10, Reproducibility 7/10, Code quality 6/10. Reduced from what strong verified statistics alone would suggest by two HIGH-severity bugs, one dead UI control, and 4 of ~24 features (MR/coloc/biomarker card/nomogram) unverified due to the terminated audit pass.

## Methylomics: **74/100** — capped due to a CRITICAL scientific-validity gap
Data handling 9/10, Preprocessing 8/10, Statistics 12/20, Visualization 11/15, Biological validity 13/15, UI/UX 6/10, Reproducibility 7/10, Code quality 8/10. Without the live-engine inflation gap this module would score in the mid-80s.

## Cross-Omics: **82/100**
Data handling 8/10, Preprocessing 8/10, Statistics 17/20, Visualization 13/15, Biological validity 13/15, UI/UX 7/10, Reproducibility 8/10, Code quality 8/10. The strongest module in the app — correctly-applied multiple testing everywhere, verified "integrate not recompute," disciplined causal language.

## Multi-Omics: **70/100** — capped due to a CRITICAL, reproduced-live defect
Data handling 8/10, Preprocessing 7/10, Statistics 16/20, Visualization 10/15, Biological validity 13/15, UI/UX 4/10, Reproducibility 6/10, Code quality 6/10. Genuinely strong underlying science (real DIABLO/SNF/MOFA2, excellent sample alignment) undercut by a confirmed live-breaking bug matching the exact user-reported symptom and the systemic error-masking pattern that hid it.

## Shared Architecture: **78/100**
Clean routing layer, well-engineered caching, no security holes — undercut by the total absence of dependency pinning/tests and one weighted-mode-MR silent-drop finding.

## Statistical Validity: **80/100**
Multiple-testing correction, leakage-avoidance, and design-matrix correctness are strong and independently verified across DGE, WGCNA, Feature Selection, Diagnostic Classifiers, and Cross-Omics. Dragged down by the Methylomics inflation gap and the Diagnostic Classifier leakage disclosure gap.

## Bioinformatics Validity: **83/100**
Platform-appropriate methods throughout, real reference data, correct methylation-specific handling (beta/M separation) — the strongest domain-competence signal in the audit. Dragged down by the same inflation-control gap.

## Visualization: **78/100**
No fabricated visualization data found anywhere; consistent plot/table/download provenance. Dragged down by one confirmed-broken chart and a recurring (cosmetic) volcano-plot axis mislabel.

## Code Quality: **72/100**
Disciplined, well-commented, heavily AI-assisted-looking code with real duplication (diverged preprocessing pipelines, sub-module toggle blocks), dead code (Platt calibration, threshold table, permutation slider), and zero tests.

## Reproducibility: **65/100**
Strong demonstrated within-app numerical reproducibility (independently verified against real output files in several places) undercut by zero dependency pinning and zero automated tests anywhere.

## Security: **88/100**
No injection, no secrets, no unsafe eval found anywhere audited. Minor hygiene issues only (hardcoded path, blanket upload cap).

## Performance: **68/100**
Real engineering exists (caching, async loading) but WGCNA/TOM/MOFA2 blocking risk and unbounded cache growth are genuine availability concerns for multi-user deployment.

## User Experience: **70/100**
Strong, specific error messaging and disciplined causal-language hedging throughout, undercut by one silently-broken chart indistinguishable from empty data and a live analytics engine (Methylomics) that gives no signal when its output may be miscalibrated.

---

# OVERALL SCORE: 74/100

**Calculation**: this is not a simple average. The four core scientific modules average to ~75 (Transcriptomics 75, Methylomics 74, Cross-Omics 82, Multi-Omics 70), and the cross-cutting dimensions average to roughly the mid-70s as well (Statistical 80, Bioinformatics 83, Visualization 78, Code Quality 72, Reproducibility 65, Security 88, Performance 68, UX 70, Shared Architecture 78). A naive average of all these numbers lands around 76.

**The overall score is capped to 74 because of two CRITICAL scientific/engineering validity issues**: the Methylomics live-engine genomic-inflation gap (§4, §9, §17.1) and the confirmed, reproduced-live Multi-Omics chart failure (§6, §17.2). Per the audit's own scoring rule, a critical scientific-validity problem caps the overall score below what dimension-averaging alone would suggest, because either finding — used in the exact live-analysis workflow it affects — could produce a wrong or silently-missing result that a user would have no way to detect from the UI alone. The score also reflects that ~15% of the Transcriptomics feature surface (MR, colocalization, biomarker card, nomogram) could not be verified in this audit round and is scored as unverified rather than assumed-correct.

---

# 20. TOP 20 IMPROVEMENTS

| Priority | Problem | Why it matters | Exact improvement | Expected impact |
|---|---|---|---|---|
| 1 | Live Methylomics DMP/DMR has no inflation control or diagnostic | Could produce miscalibrated FDR on user-run comparisons — the app's own pipeline proves this on this exact data | Add bacon or SVA correction to the live engine, or at minimum compute and display λ + a QQ plot with a warning banner above a threshold (e.g. λ>1.1) | Removes the single most consequential CRITICAL scientific-validity gap |
| 2 | Multi-Omics Adalimumab Performance chart silently broken | Confirmed live, matches reported user symptom exactly | Fix `has_snf` metadata for Adalimumab cells (`multiomics_helpers.R:53-54`) to `FALSE`; change `mod_multi_integration.R:125-127`'s `!is.null()` check to `nrow(...)>0` | Fixes the exact reported bug |
| 3 | Multi-Omics blanket `tryCatch(...,error=function(e) NULL)` masks all plot errors | This exact pattern is why bug #2 went undiagnosed and will hide the next one too | Log the caught error server-side (even just `message()`), and distinguish "no data for this selection" (a real 0-row check before plotting) from "an error occurred" (a distinct, visible UI state) | Prevents future silent failures across ~20 figure functions |
| 4 | No `na.rm=TRUE` in Batch Correction's variance/quantile filter, compounded by a sibling preprocessing path with no NA handling | Hard-crashes the app on any preloaded dataset with residual missingness | Add `na.rm=TRUE` to both calls (`mod_preprocessing.R:1436-1438`) and port the missing-value filter/impute logic from the per-source path into `pp_preloaded_read()` | Fixes a reproducible crash |
| 5 | No dependency pinning (`renv.lock`/`DESCRIPTION`) | The app's entire value proposition is "reproducible biomarker discovery" and currently nothing enforces reproducible package versions | Add `renv::init()`/`renv::snapshot()` and a `DESCRIPTION` file, committed to the repo | Makes the reproducibility claim actually true |
| 6 | No automated tests anywhere in `shiny_app/` | Every correctness claim in this audit required manual code reading; regressions in ~46,000 lines have no safety net | Add `testthat`/`shinytest2` coverage starting with the statistical helper functions (M-value math, `estimate_mr_set()`, normalization dispatch) — these are pure functions and cheap to test | Catches regressions before users do |
| 7 | Methylomics Diagnostic Classifier's internal-test AUC is optimistically biased and not visually distinguished from the valid external-test AUC | Undisclosed leakage risks a user over-trusting internal numbers, especially for models with train AUC≈1.0 | Add an explicit caveat/badge on internal-test AUC rows in the Model Comparison table, and default the headline metric to external-test AUC where available | Prevents a specific, plausible research-integrity mistake |
| 8 | Deconvolution's CIBERSORT permutation slider is a dead control | Users can change it expecting different fit-quality feedback and get nothing — and lose the ability to detect a bad tissue/platform match | Either surface CIBERSORT's P-value/Correlation/RMSE fit diagnostics in the output table, or remove the slider entirely | Restores (or honestly removes) a real quality-control signal |
| 9 | Feature IDs are never deduplicated on ingestion anywhere in Transcriptomics | Silent, undetectable data loss on any dataset with duplicate gene/probe rownames | Add an explicit duplicate-ID check at upload with a clear warning/require-resolution step | Prevents silent wrong results |
| 10 | Methylomics "t-test" Feature Selection option is mislabeled | A user selecting "t-test" is not getting a classical t-test; this is a false expectation about the method actually run | Rename the option to reflect that it's always moderated (e.g. "Moderated t-test (trend off)"), or actually implement a classical `stats::t.test` alternative | Removes a methodological mislabeling |
| 11 | All 4 Methylomics DMP/DMR volcano plots mislabel the y-axis as raw p when it's FDR | Cosmetic but present on plots reviewers/users will screenshot and cite | Parameterize the shared volcano helper to take the actual metric name, mirroring what the Candidates module's own volcano already does correctly | One-line-per-callsite fix, removes a recurring inconsistency |
| 12 | No post-training MOFA2 convergence diagnostic shown | A barely-converged model is indistinguishable from a well-converged one in the UI | Surface the ELBO trace or final delta, with a warning if training hit the iteration cap without converging | Restores a real model-quality signal for the live MOFA2 workflow |
| 13 | Cross-Omics silently collapses multi-gene CpG annotations to the first listed gene | Undercounts methylation evidence for some genes with no user-facing disclosure | Surface this in the provenance/methodology text, or switch to "one row per gene" fan-out with a documented convention | Turns a silent gap into a disclosed, defensible design choice |
| 14 | Weighted-mode MR estimate silently dropped on error app-wide | A user comparing MR methods has no way to know "Weighted mode" is missing because it errored vs. wasn't requested | Surface a visible "not available: \<reason\>" row instead of omitting it | Removes a silent-degradation risk in a causal-inference-facing module |
| 15 | Eager instantiation of all 39 submodule servers on every session | Real, uncharacterized per-session startup cost; likely the single biggest lever on session load time | Profile actual per-module server-setup cost; convert genuinely expensive setup work to lazy (on first tab-open) initialization | Likely session-load-time improvement, currently unquantified |
| 16 | WGCNA "use all genes" default can block the single-threaded app for tens of minutes for all users | Real multi-user availability risk | Move WGCNA computation to `shiny::ExtendedTask`/`future`, matching the pattern already used for the methylation matrix load elsewhere in the app | Prevents one user's analysis from freezing the app for everyone else |
| 17 | `.arthomix_cache` has no eviction/TTL | Unbounded memory growth on a long-lived deployment | Add a simple LRU cap or periodic eviction | Prevents slow memory-growth-driven instability |
| 18 | Stale "unbuilt placeholder" documentation for Cross-Omics sub-modules that are actually built | Misleads future maintainers/auditors about current app state | Update or remove the stale comment in `submodules_registry.R:60-65` | Small, cheap fix, prevents wasted future investigation |
| 19 | ~240 lines of copy-pasted sub-module toggle/search logic across 4 omics families | Maintenance/drift risk; a bugfix to one copy won't propagate to the other three | Refactor into one parameterized helper, following the pattern the codebase already uses correctly for `jump_to_*_submodule` | Reduces future-bug surface area |
| 20 | Two fully-implemented-but-dead functions in the Diagnostic module (Platt calibration, threshold-sweep table) | Implies more functionality than is actually wired to the UI; wastes future-maintainer time figuring out why they're unused | Either wire `diag_adv_platt_calibrate()` into the Calibration tab (it looks ready to use) or delete both dead functions | Either restores a real feature or removes confusing dead code |

---

# 21. WHAT I WOULD DELETE

- **The dead CIBERSORT permutation slider** in Deconvolution — either wire it to something real or remove it; leaving it is actively misleading.
- **The two dead functions in `mod_diagnostic.R`** (`diag_adv_platt_calibrate`, `diag_adv_threshold_table`) — dead code that looks finished invites future confusion about what the Calibration tab actually does.
- **The stale "unbuilt placeholder" comment block** in `submodules_registry.R` — actively wrong documentation is worse than no documentation.
- **The blanket module-wide `tryCatch(...,error=function(e) NULL)` pattern in Multi-Omics plotting** — not the error handling itself (some form is needed), but the specific version that discards the error entirely with no logging. This isn't a feature to delete so much as a pattern to replace everywhere it recurs.

I would **not** delete any analysis feature itself — no module audited was found to be decorative, fake, or scientifically unjustified. This app's problem is narrow, real bugs and one scientific-validity gap, not padding.

---

# 22. WHAT I WOULD KEEP

- **The entire Cross-Omics module** as architecturally exemplary for the rest of the app to follow: genuine "integrate not recompute" reuse, disciplined multiple-testing correction, disciplined causal-vs-associative UI language, honest disclosure of aggregation/collapse decisions (mostly — see the one HIGH finding).
- **The Diagnostic Classifier's nested-CV pipeline** (both transcriptomics and methylomics) — this is textbook-correct leakage-free machine learning, verified line-by-line, better discipline than many published papers achieve.
- **The Methylomics beta/M-value handling** throughout — mathematically verified correct and numerically stable everywhere it's used, and never confused with expression-style analysis.
- **The real DIABLO/SNF/MOFA2 multi-omics integration** — genuinely real joint modeling, not concatenation dressed up in vocabulary, with disciplined non-causal factor-interpretation language.
- **The precomputed-panel-with-verification pattern** (WGCNA hub genes, Cross-Omics/Cross-Ancestry gene counts, DIABLO variance-explained) — repeatedly independently confirmed against real output files rather than trusted on faith, which is exactly the right design for computationally expensive results.
- **ArthOChat's grounding architecture** — genuinely live, genuinely grounded in real session state and real external sources, no hardcoded secrets, sensible choice of a local model to avoid sending research data externally.

---

# 23. WHAT I WOULD REBUILD

- **The Multi-Omics plotting error-handling layer.** Current problem: a single shared `multi_plot_or_empty()` wrapper conflates "no data for this selection" with "an error occurred," with no distinguishing signal anywhere. Patching this one incident (§17.2, §20.2) fixes the Adalimumab bug specifically but leaves the exact same masking pattern in place for the next data-provenance mismatch. Recommended architecture: every plot-producing reactive should return a tagged result (`{status: "ok"|"empty"|"error", data, message}`) rather than `NULL`-or-plot, with the UI layer switching on `status` explicitly and logging every `"error"` case server-side.
- **The app-wide dependency/reproducibility layer.** Current problem: retrofitting `renv` onto an already-large, already-running app with ~25 loaded packages and no `DESCRIPTION` is nontrivial after the fact — package version drift may already have occurred silently. Recommended approach: snapshot the current working environment now (`renv::init()` against the live, working deployment) rather than trying to reconstruct "what versions were actually used" later.
- **The `pp_preloaded_read()` / per-source preprocessing duplication.** Current problem: two independently-maintained copies of essentially the same load→filter→log2 pipeline have already diverged once (missing NA handling in one, §17.3) and will diverge again under future edits. Recommended: collapse both into one parameterized function with the preloaded path as a thin caller.

---

# 24. RECOMMENDED ROADMAP

## Phase 1 — Critical scientific fixes
- Add bacon/SVA correction (or a clearly-surfaced λ/QQ inflation diagnostic) to the live Methylomics DMP/DMR engines.
- Disclose Diagnostic Classifier internal-test leakage in the UI; default to external-test AUC as the headline metric where available.

## Phase 2 — Broken functionality
- Fix the Multi-Omics `has_snf` metadata + `nrow()>0` check to resolve the Adalimumab Performance chart.
- Add error handling to the Multi-Omics PNG download handler.
- Fix the Batch Correction `na.rm=TRUE` crash and unify the two diverged preprocessing pipelines.
- Wire or remove the dead CIBERSORT permutation slider and the two dead Diagnostic-module functions.

## Phase 3 — Statistical improvements
- Surface small-n warnings on Cross-Omics's recomputed-FDR option.
- Add explicit multiple-testing disclosure to Cross-Ancestry's UI (not just code comments).
- Surface a visible "not available" state instead of silently dropping the weighted-mode MR estimate.

## Phase 4 — Visualization improvements
- Fix the 4 mislabeled Methylomics volcano-plot y-axes.
- Add a post-training MOFA2 convergence diagnostic.

## Phase 5 — Performance
- Move WGCNA (and any other similarly expensive live-compute step) onto `shiny::ExtendedTask`, matching the pattern already used elsewhere in the app.
- Add a TOM-computation gene-count safety cap.
- Add eviction/TTL to `.arthomix_cache`.

## Phase 6 — UX
- Replace the app-wide "silent empty state on error" pattern with an explicit, logged, user-visible error state wherever it currently exists.
- Add duplicate-feature-ID detection and warning at upload.

## Phase 7 — Production hardening
- Add `renv.lock`/`DESCRIPTION` and pin dependencies.
- Add a `testthat`/`shinytest2` suite, starting with the pure statistical helper functions.
- Refactor the ~240-line duplicated sub-module toggle block into one parameterized helper.
- Update the stale Cross-Omics documentation.

---

# 25. FINAL SENIOR-LEVEL VERDICT

### Scientific maturity: 76/100
### Software engineering maturity: 70/100
### Bioinformatics maturity: 83/100
### Statistical maturity: 80/100
### Production readiness: 58/100
### Research readiness: 78/100
### Overall: 74/100

**If this application were presented to a senior bioinformatics PI**, the first thing they would criticize is the live Methylomics engine's lack of inflation control — precisely because the app's *own* pipeline proves the problem exists on this exact data, making its absence in the live tool look like an oversight rather than a considered tradeoff. The second thing would be the total absence of dependency pinning for an app that markets itself on reproducibility.

**If this application were reviewed by a computational biology journal reviewer**, they would ask for (a) the internal-vs-external test AUC distinction to be made explicit for the Methylomics Diagnostic Classifier, given the shared-cohort feature-selection leakage, and (b) a demonstration that the live-computed panels (Methylomics DMP/DMR, live MOFA2) produce results consistent with the precomputed/validated panels the app already ships — since right now the precomputed panels are more trustworthy than the "run it yourself" workflows, which inverts the normal expectation for an interactive research tool.

**What would prevent me from trusting its results**, specifically: any result produced by the live Methylomics DMP/DMR engine, until the inflation gap is addressed — I would insist on cross-checking against the default/SVA panel first. Everything else audited — the precomputed panels, the transcriptomics statistics, the cross-omics integration, the MOFA2/DIABLO/SNF machinery — I would trust as a starting point for further validation, which is the appropriate level of trust for exploratory research software.

**The three changes that would increase scientific credibility the most**:
1. Fix or clearly banner-warn the Methylomics live-engine inflation gap.
2. Add `renv`/dependency pinning and a minimal test suite — not because anything found was wrong due to a version mismatch, but because nothing in the repo currently rules that out.
3. Replace the Multi-Omics silent-error-swallowing pattern with a visible, logged error state — not just for the one bug found, but because it's the single architectural choice most likely to hide the *next* one.
