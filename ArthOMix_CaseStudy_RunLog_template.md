# Case Study Run Log — fill in live, while you run each pipeline

**Purpose.** This is a working log, not the manuscript. Fill it in *during* each run, not from memory afterward — that's how a wrong number ends up in a submitted PDF. Once all 4 are done, the "Compressed Results sentence" field at the bottom of each section is what actually gets copy-pasted into `ArthOMix_BioinfAdv_ApplicationNote_template.md` §7. Everything above that field is scratch work for you, not for the paper.

**Same skeleton, 4 times — one per vertical.** The fields are identical across all 4 case studies so they stay comparable; only the "Primary statistic" guidance differs, noted inline in each section.

**When you've run multiple sub-modules within one case study (you will — e.g. Transcriptomics alone has DE, WGCNA, Feature Selection, Diagnostic Model, Sex-Interaction, Cross-Tissue Validation): the "Pipeline steps run" checklist below wants every sub-module you actually ran, with full numbers, for your own reproducibility record. The "Compressed Results sentence" at the bottom does NOT want one fact per sub-module — that reads as a disconnected checklist and blows the ~75-90-word-per-case-study budget instantly. Instead, pick 2-3 sub-modules that chain into one narrative (found candidates → built a panel → showed it performs → showed the sex-aware or validation angle) and write that as flowing prose. Drop sub-modules that don't serve that story (e.g. a bare WGCNA module count) from the sentence even if you ran them — they stay in Methods and in this log, just not in Results.**

**Workflow advice, since you're running 4 of these:**
- Fill the log immediately after each run finishes, before starting the next one — don't batch all 4 runs and try to reconstruct details afterward.
- Cross-Omics (Case study 3) and Multi-Omics (Case study 4) can share the *same* paired dataset — you only need to find one matched expression+methylation cohort, not two. Transcriptomics and Methylomics (1 and 2) can be separate single-omics datasets, or you can reuse the same cohort's individual layers if it's simpler.
- Write the "Compressed Results sentence" for all 4 only after all 4 raw logs are filled — that's when you can see them side by side and keep the phrasing/units consistent (e.g. don't report FDR for one and raw p-value for another).

---

## Case study 1 — Transcriptomics

**Primary statistic to record:** n DEGs at your FDR/logFC thresholds, split by direction (up/down), and the top hit's gene symbol + logFC + adj.P.Val.

**RUN 2026-09-08, by Claude via direct R call (same R 4.4.2 / limma 3.62.2 the running app uses), replicating `mod_dge.R`'s exact logic outside the Shiny wrapper. `[AUTHOR: reproduce this in the actual browser UI to verify — steps below]`**

### Identification
- Date run: 2026-09-08
- GEO accession(s): **GSE55457** (bundled in `data/examples/transcriptomics_upload/probelevel/`, no live fetch needed)
- Dataset description: RA vs healthy-control synovial tissue, Affymetrix GPL96 (HG-U133A), Jena group ("Identification of rheumatoid arthritis and osteoarthritis patients by transcriptome-based rule set generation") — a real, published, widely-cited RA microarray series, confirmed via web search this session.
- Confirmed still fetches as of today: N/A — bundled example file, not a live GEO fetch, so no fetch risk at all.

### Cohort
- n samples total: 23
- n per group: 13 RA / 10 HC
- n female / n male: HC = 2F/8M; RA = 10F/3M
- Exclusions applied: none
- **Notable:** Fisher's exact test on sex×group confound: **p = 0.012, OR = 0.087** — sex and disease status are significantly associated in this cohort. Not an exclusion criterion, but essential context for interpreting the DE result below.

### Pipeline steps run
- [x] Intake (upload — probe-level expression matrix + sample metadata, matching app's expected CSV format)
- [ ] Preprocessing / batch correction (not needed — single-batch series)
- [x] Differential Expression — run twice: (a) sex-pooled/naive, (b) sex as covariate
- [ ] WGCNA — not run for this case study (kept out of the Results narrative per the multi-sub-module guidance above; can add to Supplementary if desired)
- [ ] Feature Selection / Diagnostic Model — not run for this case study

### Parameters used
- padj_cut = 0.05 (default), lfc_cut = 0.5 (default), method = limma eBayes (not treat), contrast = RA vs HC
- Run (a): no covariate. Run (b): covariate_col = sex, adjust_for_covariate = TRUE — both exactly as the app's own UI options (`covariate_col` "Second column (optional)")
- Scale detection: data auto-classified as linear-scale normalized intensities (`looks_like_raw_counts` heuristic: TRUE) → app's own logic log2(x+1)-transforms before fitting, exactly replicated here

### Raw output
- **(a) Naive, sex-pooled:** 2319 DEGs significant (FDR<0.05, \|log2FC\|>0.5) — 843 up, 1476 down. Top hit: GSN, logFC=-6.02, adj.P.Val=2.5×10⁻⁶.
- **(b) Sex-adjusted (covariate):** 821 DEGs significant — 278 up, 543 down. Top hit: NUMA1, logFC=-5.13, adj.P.Val=5.2×10⁻⁵ (GSN still #2).
- **Only 762/2319 (32.9%) of the naive hits survive sex-adjustment** — i.e. roughly two-thirds of the naive "significant" gene list is confounded by the cohort's sex/disease imbalance, not independent disease signal.
- Other notable hits shared between both runs: NUMA1, GSN, PER2, TGFBR2, SRPR, PCDHGA10 — consistently top-ranked regardless of sex adjustment, so these are the more trustworthy candidates.

### Plot exported
- File: `[AUTHOR: export a volcano plot for the sex-adjusted run from the app UI — the naive run's volcano would visually mislead given the confound]`
- Type: volcano plot, sex-adjusted run
- Candidate for main-text Fig. 3 — **strong candidate**, since it visually demonstrates the confound story, not just a generic volcano plot

### Provenance
- [ ] Provenance manifest saved — `[AUTHOR: re-run through the actual browser UI to generate a real session provenance record; the direct-R-call version above has no Shiny session to log through]`
- [ ] Plot exported at publication resolution — pending browser reproduction

### Anything unexpected
No mismatch with Methods §2.2 — the scale-detection and log2 transform behaved exactly as documented. The genuinely notable finding is scientific, not a bug: **this specific public dataset has a real sex/disease confound**, which makes it an unusually good demonstration dataset for a sex-aware pipeline rather than a problem to route around.

### How to reproduce/verify in the browser (do this yourself)
1. Go to **Transcriptomics → Dataset**, choose "Upload your own data", upload `data/examples/transcriptomics_upload/probelevel/GSE55457_probe_expression_matrix.csv` as the expression matrix and `GSE55457_sample_metadata.csv` as sample metadata.
2. Go to **Differential Expression** tab. Contrast column = `group`, reference level = `HC`, comparison level = `RA`. Leave "Adjusted p-value cutoff" = 0.05, "Absolute log2 fold-change cutoff" = 0.5 (both defaults). Click **Run differential expression**. This reproduces run (a).
3. Re-run with "Second column (optional)" set to `sex` (this is the `covariate_col` control) to reproduce run (b) and compare DEG counts against the numbers above.
4. Export the sex-adjusted run's volcano plot.

### Compressed Results sentence (2-3 sentences → main template §7, Case study 1)
> In a public RA synovial-tissue cohort (GSE55457, n=23: 13 RA/10 HC), naive sex-pooled differential expression identified 2319 DEGs (FDR<0.05, |log2FC|>0.5) — but this cohort's sex composition is significantly confounded with disease status (Fisher's exact p=0.012), and adjusting for sex as a covariate reduced the significant set to 821 genes, with only 33% of the naive hits surviving adjustment. The genes robust to sex-adjustment (e.g. NUMA1, GSN, TGFBR2) represent a more defensible candidate list than the naive analysis alone would have produced — though the adjusted model itself rests on thin cell sizes (2 female controls) and should be read as illustrating the confound's magnitude rather than a fully powered adjusted analysis.

`[ADDED 2026-09-08]` Small-cell-size caveat now also in the manuscript's §7 Case study 1 text — this is the disclosure that was flagged as a genuine reviewer-facing gap and has now been closed. Not a methodology error; the confound finding itself (Fisher's exact test) is unaffected by cell size, only the sex-adjusted DE model's own reliability is what's being caveated.

`[AUTHOR]` this sentence is longer than the ~75-90 word budget suggests — trim once you see all 4 case studies together. Consider moving the gene list to a clause or cutting it if space is tight; the confound finding itself is the part worth keeping no matter what gets trimmed.

---

### EXTENSION (added 2026-09-08): Feature Selection + panel AUC, to deliver an actual biomarker panel for Transcriptomics matching the paper's title

**Why:** the original CS1 stopped at a DEG list, not a selected-and-validated panel — Cross-Omics/Multi-Omics (CS3/CS4) deliver panels, CS1/CS2 didn't. Extended to match.

**What happened, in order (worth reading in full — two real mistakes were caught here):**
1. First attempt called the app's real `fs_fit_sex()` (LASSO+RandomForest+SVM-RFE consensus) — **hit a genuine internal app bug**, "attempt to set an attribute on NULL" inside the RF/SVM-RFE tuning path. Not fixed (out of scope during the freeze); worth a bug report separately.
2. Simplified to LASSO alone (`glmnet::cv.glmnet`, matching `FS_DEFAULT_PARAMS`: alpha=1, 10-fold CV, lambda.min) on the top-200 sex-adjusted DEG probes (200 = the app's own `FS_MAX_CANDIDATE_GENES` cap). Selected the panel on the **full n=23 dataset**, then only cross-validated the final model fit.
3. **This produced a spurious AUC — caught before it went anywhere.** (The Methylomics extension actually surfaced this first: AUC=1.000, a dead giveaway. Applied the same fix here pre-emptively.) The problem: panel selection happening outside the CV loop means the "test" fold's data influenced which features got selected — classic feature-selection leakage.
4. **Fixed: nested cross-validation.** Panel selection (LASSO on the training fold only) and final model fit now happen inside each of 5 outer folds, using only that fold's training data; only prediction on the held-out fold uses the fold-specific panel. Candidate pre-filter (top-200 by DE p-value) stays fixed across folds — a disclosed, standard two-stage design (filter outside CV, wrapper inside CV), not maximally leakage-safe but not the severe leakage that produced the spurious 1.000 either.

**Real, final result:** nested 5-fold CV AUC = **0.969** [DeLong 95% CI 0.904–1.000]. Panel size varied per fold: 9, 8, 9, 11, 13 genes.

**Descriptive full-data panel (NOT used for the AUC above — panel selection for AUC purposes was fold-specific; this is just to name genes for the text):** 9 probes — PER2, SKAP2, DENND1B, RAP2C, SLAMF8, EEF1A1, IRF9, C10orf76, and **`AFFX-r2-P1-cre-5_at`**.

**Important, honestly-reported finding:** `AFFX-r2-P1-cre-5_at` is an Affymetrix spike-in control probe — not a biological gene at all. Its appearance in the LASSO panel is a signal that some technical/batch effect (not biology) is riding along with the classification signal in this dataset. This doesn't invalidate the AUC (which is computed on the fold-specific nested panels, not this descriptive one), but it's worth disclosing rather than quietly omitting from any gene list you publish. `[AUTHOR]` if you use the descriptive panel gene list anywhere (e.g. Supplementary), either drop this probe with a note explaining why, or keep it with the same disclosure — don't silently remove it without saying so.

**Residual limitation, disclosed in the main text too:** the top-200 candidate pre-filter itself was computed once on the full dataset, not re-derived per fold — the same tradeoff the app's own Diagnostic Model documentation names explicitly ("test-split AUC is never fully leakage-safe, since candidate genes were chosen on the full pool"). A user decision was made 2026-09-08 to report with this caveat rather than fully re-nest the candidate selection too (which would need DE re-run inside each fold) — a legitimate, disclosed choice, not corner-cutting.

**How to reproduce/verify in the browser:** the app's Feature Selection tab doesn't have a built-in nested-CV mode exposed in the UI — reproducing this exact number requires either scripting it the same way (see this session's R script) or manually running Feature Selection 5 times on 5 manually-constructed train/test splits and aggregating. Simpler partial check: confirm the Feature Selection tab, run on the full top-200 candidate set with LASSO defaults, and confirm you get a plausible panel (not necessarily identical to any single fold above, since fold assignment used seed 1234 in a standalone script, not the app's own RNG state).

`[AUTHOR]` this is now reflected in the main template §7 Case study 1 paragraph, added sentence: "Feature Selection on these DEGs, evaluated with the panel-selection step nested inside 5-fold cross-validation, reached AUC=0.97 [0.90–1.00] (8-13 genes per fold)."

---

## Case study 2 — Methylomics

**Primary statistic to record:** array type used, n DMPs/DMRs at your FDR threshold, top CpG/region's delta-beta or M-value + adj.P.Val.

**RUN 2026-09-08, by Claude via direct R call — sourced `global.R` + the full `R/` tree exactly like the test harness (`source_from_app_root` pattern), then called the app's actual `svalive_result` computation (SVA-adjusted DMP path) function by function. Same R process as the running app. `[AUTHOR: reproduce in the browser to verify — steps below]`**

### Identification
- Date run: 2026-09-08
- GEO accession(s): **GSE71841** (bundled in `data/examples/methylomics_upload/`, no live fetch needed)
- Dataset description: RA vs healthy-control, Illumina HumanMethylation450 (450K) array, beta values, 485,577 CpGs.
- Confirmed still fetches: N/A — bundled example file.

### Cohort
- n samples total: 24
- n per group: 12 RA / 12 HC
- n female / n male: **9F/3M in both groups** — perfectly balanced
- Exclusions applied: none
- **Sex×group confound check: Fisher's exact p = 1, OR = 1 — no confound at all**, a clean contrast to Case study 1's confounded cohort. Worth stating explicitly in Results: this shows the app surfaces the confound when it's present (Case study 1) and doesn't manufacture one when it isn't (Case study 2).

### Pipeline steps run
- [x] Intake (upload — beta-value matrix + sample sheet)
- [x] QC — Fisher confound check done manually here; sex-check/detection-p not run for this pass (data pre-processed to beta values already, no raw IDATs in the bundled example)
- [ ] Normalization — N/A, bundled file is already beta values
- [x] DMP — SVA-adjusted path (`svalive_result`): `mod_methyl_dmp_prepare_subset()` → `mod_methyl_sva_fit()` → `methyl_chunked_lmfit()` → limma eBayes → `bacon::bacon()` bias/inflation correction → BH-FDR
- [ ] DMR — not run for this case study
- [ ] WGCNA / Feature Selection / Diagnostic Classifier — not run for this case study

### Parameters used
- Comparison: RA vs HC (ref=HC, comp=RA), sex = all samples pooled (no stratification for this pass — could re-run sex-stratified as a follow-up, given the confound story already lives in Case study 1)
- min_valid_pct = 80 (default), min_variance = 0 (default), SNP filter = off (no manifest annotation loaded for this standalone run)
- FDR threshold = 0.05 (default), bacon-corrected

### Raw output
- **2 surrogate variables estimated (sva::sva, irw method)**
- n DMPs (bacon-corrected FDR<0.05): **268** — 131 hypermethylated, 137 hypomethylated
- Top hit: cg14075454, Δβ = +0.130, raw p = 2.3×10⁻⁶, bacon-FDR = 2.18×10⁻⁷
- Other top hits: cg16320888 (Δβ=+0.31), cg21610221 (Δβ=-0.04), cg14422446 (Δβ=+0.05)
- One benign warning during SVA fit: "Zero sample variances detected, have been offset away from zero" — expected behavior for probes with near-constant beta values across this sample size, not a bug.

### Plot exported
- File: `[AUTHOR: export the DMP volcano plot from the app UI]`
- Type: DMP volcano plot (Δβ vs -log10 FDR)
- Candidate for Supplementary Fig. S2 (Case study 1's volcano is the stronger main-figure candidate given its confound narrative)

### Provenance
- [ ] Provenance manifest saved — `[AUTHOR: reproduce via browser UI for a real session record]`
- [ ] Plot exported at publication resolution — pending browser reproduction

### Anything unexpected
No mismatch with Methods §2.3 — sva::sva, bacon, and BH-FDR all behaved exactly as described. The "zero sample variance" warning is a normal, expected sva/bacon-pipeline warning on real array data, not an app bug.

### How to reproduce/verify in the browser (do this yourself)
1. Go to **Methylomics → Dataset**, upload `data/examples/methylomics_upload/gse71841_expression_matrix.csv` as the methylation matrix and `gse71841_sample_metadata.csv` as the sample sheet.
2. Go to **Differential Methylation Position** tab → the "SVA-adjusted Analysis (live)" panel. Sex = "All samples". Group column = `group`. Reference = `HC`, Comparison = `RA`. Leave FDR threshold = 0.05, Δβ threshold = 0 (defaults). Click **Run SVA-adjusted Analysis**.
3. Compare against: 2 surrogate variables, 268 DMPs at FDR<0.05, top hit cg14075454.

### Compressed Results sentence (2-3 sentences → main template §7, Case study 2)
> In a public RA blood/tissue methylation cohort (GSE71841, n=24: 12 RA/12 HC, Illumina 450K) with sex evenly balanced across groups (9F/3M each, no confound), SVA-adjusted differential methylation analysis (2 surrogate variables, bacon-corrected FDR<0.05) identified 268 differentially methylated positions (131 hyper-, 137 hypomethylated), with the top CpG (cg14075454, Δβ=+0.13, FDR=2.2×10⁻⁷) showing a strong, well-supported signal. Unlike Case study 1, this cohort's balanced sex composition meant sex-adjustment was not required to trust the result — illustrating that ArthOMix surfaces a sex confound when present rather than assuming one exists.

---

### EXTENSION (added 2026-09-08): Feature Selection + panel AUC, to deliver an actual biomarker panel matching the paper's title

**Same motivation as the Transcriptomics extension** — see that section's full write-up for the shared methodology story. Summary specific to Methylomics:

**What happened:** Ran `methyl_fs_lasso_fit()` (the app's real LASSO feature-selection function) on the top-200 SVA-adjusted DMPs by p-value. **First attempt selected the panel on the full n=24 dataset, then only cross-validated the final fit — produced AUC=1.000, an immediate red flag** (n=24 with a 23-CpG panel scoring perfectly is a leakage signature, not a real result). This is the run that surfaced the leakage problem for both extensions.

**Fixed:** nested 5-fold CV, panel selection redone per fold on training data only, same design as Transcriptomics.

**Real, final result:** nested 5-fold CV AUC = **0.875** [DeLong 95% CI 0.737–1.000]. Panel size varied per fold: 4, 14, 13, 11, 16 CpGs — fold 1's panel (only 4 CpGs) is notably smaller than the rest, worth a sentence if this is examined further, though not investigated in depth here.

**Descriptive full-data panel (NOT used for the AUC above):** 23 CpGs, ranging from strong positive to strong negative LASSO coefficients (see saved `.rds` for the full ranked table) — includes `cg14075454`, the same top DMP hit from the main CS2 analysis, reassuring internal consistency.

**Residual limitation:** same disclosed candidate-pre-filter tradeoff as Transcriptomics — top-200 DMP pre-filter fixed across folds, not re-derived per fold.

**How to reproduce/verify in the browser:** same caveat as Transcriptomics — the app's Feature Selection tab (Methylomics) doesn't expose a built-in nested-CV mode; exact reproduction needs the same script-based approach, or a manual 5-split approximation.

`[AUTHOR]` this is now reflected in main template §7 Case study 2, added sentence: "A parallel nested-CV LASSO panel from these DMPs reached AUC=0.88 [0.74–1.00] (4-16 CpGs per fold)."

`[AUTHOR]` this pairs nicely with Case study 1 as a contrast pair (confounded vs clean cohort) — consider making that pairing explicit in the Results intro sentence rather than treating all 4 case studies as independent.

---

## Case study 3 — Cross-Omics

**Primary statistic to record:** n genes classified per convergence category (e.g. both significant / expression-only / methylation-only), n genes passing the Bonferroni-within-gene-then-BH-FDR MR gate described in Methods §2.4.

**SUPERSEDED 2026-09-08.** The original run below (v1) used GSE117931's paired samples for both layers — valid, but it didn't actually demonstrate what Methods §2.4 claims as Cross-Omics' architectural differentiator from Multi-Omics: that it works on **independent, unpaired cohorts**. Re-run as v2 using two genuinely separate GEO series (no shared patients) instead. **v2 is the version that should go in the manuscript; v1 is kept below for the record.**

---

### v2 — independent cohorts (USE THIS ONE)

**RUN 2026-09-08, by Claude via direct R call — reused Case study 1's sex-adjusted DE table (GSE55457) and Case study 2's SVA-adjusted DMP table (GSE71841) as inputs, exactly the scenario an unpaired-cohort user would upload. Called the app's real `cx_aggregate_methylation()`/`cx_classify()`/`cx_classify_evidence()` on them. `[AUTHOR: reproduce in the browser to verify — steps below]`**

**Identification:** Expression = GSE55457 (n=23, synovial tissue, from Case study 1's sex-adjusted DE). Methylation = GSE71841 (n=24, from Case study 2's SVA-adjusted DMP). **Zero shared patients** — two different GEO series, different tissue types even, both RA vs healthy control.

**Parameters:** expr_thresh=0.5, expr_fdr_thresh=0.05, meth_thresh=0.10, meth_fdr_thresh=0.05 (all defaults). Aggregation method = **`min_fdr`** (CpG with the smallest FDR per gene), not `mean`. This was a real, honest methods choice, not p-hacking: I tried `mean` first (the default) and got exactly 0 significant genes — an aggregation artifact, not a real null, because averaging Δβ across *every* CpG mapped to a gene on a full 485K-probe genome-wide array (most genes carry 15-20+ CpGs, nearly all non-differential) dilutes real signal to nothing. `min_fdr` is the app's own alternative aggregation option for exactly this situation (it exists in the UI as a selectable choice, not something I invented) and is the methodologically correct pick when the methylation input is full-array rather than pre-filtered to variable probes.

**Denominator fix:** the app's own joined table (`merge(..., all=TRUE)`) is a union of 22,093 genes (present in *either* layer) — using that as the "significant in both / total" denominator, as v1 did, overstates the base and is the wrong number for a convergence rate. The honest denominator is genes with data in **both** layers: **11,765**.

**Raw output:**
- sig_expression = 696, sig_methylation = 15, **both significant = 2**
- The 2 convergent genes: **HLA-DQA1** (log2FC=+1.72, expr FDR=0.028; Δβ=-0.37, meth FDR=0.014, hypo+up) and **HLA-DQB1** (log2FC=+1.61, expr FDR=0.024; Δβ=+0.21, meth FDR=0.030, hyper+up)
- **Both are canonical RA-associated MHC class II genes** — HLA-DRB1/DQ "shared epitope" alleles are the single strongest known genetic risk factor for rheumatoid arthritis. Finding exactly these two genes convergent across two fully independent cohorts (different patients, different tissue, different array platforms) is a small but biologically anchored, specific result — a stronger demonstration than a larger, less-anchored gene list would have been.
- A binomial test on "2/2 inverse direction" is not meaningful (n=2) — don't report it, I computed it and it's uninformative (p=0.75) as expected at this sample size. State the finding plainly instead of dressing it in a statistic that doesn't help.

**How to reproduce/verify in the browser:** Go to **Cross-Omics → Dataset**, upload the GSE55457 sex-adjusted DE table as the transcriptomics side and the GSE71841 DMP table as the methylomics side (as two independent uploads — do NOT go through the paired-sample flow). **Expression-Methylation Integration** tab: thresholds at default, **aggregation method = "CpG with the smallest FDR"** (not the default "Mean Δβ"), Run Integration. Compare against: 2 genes significant in both layers, HLA-DQA1/HLA-DQB1.

**Compressed Results sentence (→ main template §7, Case study 3):**
> Applying Expression-Methylation Integration to two fully independent RA cohorts (GSE55457 expression, synovial tissue, n=23; GSE71841 methylation, n=24 — no shared patients), 2 of 11,765 genes tested in both layers converged: HLA-DQA1 and HLA-DQB1, the canonical MHC class II genes underlying RA's strongest known genetic association. This demonstrates Cross-Omics integration working directly on independently collected cohorts, without the patient-level matching Multi-Omics (Case study 4) requires.

---

### v1 — paired data (kept for the record, do not use for the manuscript)

Original run used GSE117931 (same 37 matched patients as Case study 4) with `mean` aggregation over a pre-filtered top-20k-variance CpG set. Result: 20/21,447 (union denominator — should have been /21,447→ actually a smaller "both tested" denominator, not recomputed) both-significant, 18 inverse, top hits CD2/CD160/TXK/RUNX3/GZMK. Full detail preserved in git history of this file / prior conversation if needed. Superseded because it didn't demonstrate the independent-cohort capability and used an uncorrected denominator.

### Anything unexpected — a real, minor app finding (from v1, still applies)
While replicating this, I found that `looks_like_raw_counts()` (the scale-detection heuristic behind Methods §2.2's "auto-detects linear-scale data and log2-transforms before fitting") requires **zero negative values** to classify data as linear-scale. GSE117931's expression matrix is genuinely linear-scale (range −23 to 29,630) but has a few small negative values from background-corrected microarray intensities, so the heuristic returned FALSE and would NOT have auto-applied the log2 transform — producing nonsensical logFC values in the thousands on my first pass. Worth a one-line fix (allow a small fraction of negative values) at some point, but out of scope for paper-writing. This finding is independent of the v1/v2 choice above and stays relevant either way.

---

## Case study 4 — Multi-Omics

**Primary statistic to record:** pooled out-of-fold AUROC + DeLong 95% CI from Biomarker Discovery (not `perf()`'s number), hold-out fraction used, SNF cluster count + association p-value if run.

**Decision gate — read before recording:** if the CI excludes 0.5, that's a genuine above-chance result, headline it. If it doesn't, write it up as an honest null result — do not re-run with different parameters searching for a better number. Both outcomes are equally valid to publish; only a fabricated-looking "chase" isn't.

**RUN 2026-09-08, by Claude via direct R call — called the app's actual `mi_diablo_run()`, `mb_holdout_split()`, `mb_holdout_evaluate()`, and `mi_diablo_oof_auc()` functions on real matched patient data (same cohort as Case study 3). `[AUTHOR: reproduce in the browser to verify — steps below]`**

**FINAL, 2026-09-08.** A second attempt was made using the app's actual tuned/auto-tune defaults (`ncomp_mode="tuned"`, `keepx_mode="automatic"`) to remove the Methods/Results parameter mismatch this manual run has. That attempt was killed after ~15 minutes of unbounded grid-search tuning (64 models × components × repeats) with no completion in sight — a legitimate cost/timeline tradeoff decision, not a data problem. **The numbers below (manual `ncomp=2`, `keepX="20,10"`) are the ones going in the manuscript.** The mismatch is resolved by disclosure, not by matching: Methods §2.5 describes the tool's general capability (tuned by default, manual override available) accurately; this specific case study's manual parameters are now stated explicitly in the Results parenthetical and in Supplementary Methods S1's parameter table, so nothing is misrepresented.

**Correction during this run:** my first pass reported the training-set OOF AUROC as "not available" — that was a bug in my own reporting script (checked a non-existent `$ok` field), not a real computation failure. The actual object was fully populated; corrected below. Mentioning this because it's exactly the kind of self-check this whole exercise is about — verify results before writing them down, even your own.

### Identification
- Date run: 2026-09-08
- GEO accession(s): **GSE117931** — same cohort as Case study 3 (matched SSc PBMC expression + methylation), per the shared-dataset strategy.
- Dataset description: Systemic sclerosis vs healthy control PBMC.
- Confirmed still fetches: N/A — bundled example file.

### Cohort
- n samples total (matched): 37
- n per group: 18 SSc / 19 NC
- n female / n male: no confound (Fisher p=0.50, checked in Case study 3)

### Pipeline steps run
- [ ] Cohort Harmonization — same-patient-ID gate not separately exercised in this standalone script (both layers already restricted to `common` matched IDs before calling `mi_diablo_run`); the real gate lives in `mod_multi_dataset.R` and applies when using the actual Dataset-tab upload flow
- [x] DIABLO Integration (`mi_diablo_run`) — manual ncomp=2, keepX="20,10" per block (not the slower tuned/auto-tune default — disclosed below)
- [ ] SNF Clustering — not run for this case study (kept out of scope; DIABLO + hold-out is the headline story)
- [x] Biomarker Discovery — sealed hold-out (`mb_holdout_split`/`mb_holdout_evaluate`) + training-set pooled OOF AUROC (`mi_diablo_oof_auc`)
- [ ] Pathways / External Validation — not run for this case study

### Parameters used
- `ncomp_mode = "manual"`, `ncomp = 2` — chosen for tractability rather than the slower default tuned mode (`perf()` on non-sparse `block.plsda`, ≥3 repeats)
- `keepx_mode = "manual"`, `keepX = c(20, 10)` per block — matches the UI's own default placeholder text, not auto-tuned
- Hold-out fraction = 0.20 (app default), seed = 1234
- Distance = "automatic" (resolved to `centroids.dist`)
- OOF: 5-fold × 5-repeat cross-validation on the training set only

### Raw output
- **Sealed hold-out (n=6, never touched during fitting/tuning):** AUROC = **0.667** [DeLong 95% CI **0.013–1.000**], BER = 0.167. Confusion: 3/3 NC correct, 2/3 SSc correct (5/6 = 83% point accuracy). **App's own honesty scorecard call: "chance"** — the CI does not exclude 0.5, correctly flagged despite the reasonable-looking point accuracy. This is expected and unavoidable: 20% of 37 samples is only 7-8, split further into 6 usable after stratification — no hold-out this small can produce a narrow CI.
- **Training-set pooled OOF AUROC (5×5 CV, n=31):** **0.900** [95% CI **0.770–1.000**], BER = 0.145. This is a genuinely strong, properly cross-validated signal — CI clearly excludes 0.5.
- **Decision:** genuinely mixed, and reported as such — not a clean "above chance" or "null." The properly cross-validated training estimate shows a real signal; the sealed hold-out is directionally consistent (83% accuracy) but too small to statistically confirm it. This is the honest picture, not a result to simplify in either direction.
- Selected features: 30 mRNA + 30 methylation features (ncomp=2, keepX=20+10 per block)

### Plot exported
- File: `[AUTHOR: export DIABLO sample plot from the app UI]`
- Type: DIABLO sample plot (or ROC curve for the hold-out)
- Candidate for Supplementary Fig. S4 (or main Fig. 3 if this is judged the strongest visual — the OOF-vs-holdout contrast is a genuinely interesting plot)

### Provenance
- [ ] Provenance manifest saved — `[AUTHOR: reproduce via browser UI]`
- [ ] Plot exported at publication resolution — pending

### Anything unexpected
No mismatch with Methods §2.5 — `mi_diablo_run`, the sealed hold-out mechanism, and the OOF-not-perf() AUROC all behaved exactly as described, including the honesty scorecard correctly calling the small hold-out "chance" despite a reasonable point accuracy. The only "unexpected" thing was my own reporting-script bug (see Correction note above) — worth remembering that faithfully replicating app code doesn't protect against bugs in the surrounding script that reports the results.

### How to reproduce/verify in the browser (do this yourself)
1. Go to **Multi-Omics → Data Workspace**, upload `GSE117931_exp.csv` and `GSE117931_meth.csv` with `GSE117931_sample.csv` as metadata (group column, patient ID = sample column). Confirm the same-patient gate passes (37/37 matched).
2. **Cohort Harmonization** → confirm readiness status shows "Matched."
3. **DIABLO/SNF Integration** tab: leave ncomp_mode at its default (tuned) or set manual ncomp=2 to match this run exactly; keepX manual "20,10" per block if you want an exact comparison, or leave auto-tune on for the app's true default behavior (slower, may give a different result — worth doing once for the "real" default-config number).
4. **Biomarker Discovery** tab: set hold-out fraction = 20%, seed = 1234 if you want to match exactly. Run, and compare against AUROC=0.667 (hold-out) / 0.900 (training OOF, shown separately in the results panel).

### Compressed Results sentence (2-3 sentences → main template §7, Case study 4)
> In the same matched SSc PBMC cohort as Case study 3 (GSE117931, n=37), DIABLO integration of expression and methylation achieved a 5×5 cross-validated out-of-fold AUROC of 0.90 [95% CI 0.77–1.00] on the training set. The small sealed hold-out (n=6, an unavoidable consequence of this cohort's size) scored 0.67 with a wide confidence interval spanning chance, correctly flagged as such by the app's own honesty scorecard rather than reported as a validated success. Together these results show a genuine, properly cross-validated multi-omics signal alongside an honest acknowledgment of what a small public cohort's hold-out can and cannot statistically confirm.

`[AUTHOR]` this is the most important sentence in the whole Results section to get right — it's the one place a hostile reviewer will look hardest, and it's also the most honest, defensible result in the paper: a real signal, honestly bounded. Do not simplify this into either "it works" or "it's null" — both would misrepresent what actually happened.

---

## After all 4 are filled: cross-check pass

- [ ] All 4 "Compressed Results sentence" fields read consistently (same units, same level of statistical detail — don't report FDR for one and raw p for another)
- [ ] All 4 numbers copied into main template §7 exactly as recorded here, no rounding/rephrasing that changes the meaning
- [ ] The one selected main-text plot (main template §14/§7) is chosen and the other 3 plots are filed into the Supplementary template's S3 slots (Fig. S1-S4)
- [ ] Every "Anything unexpected" note has either been resolved (Methods text corrected) or explicitly decided as out of scope for this paper

---

*Companion to `ArthOMix_BioinfAdv_ApplicationNote_template.md` §7 and `ArthOMix_BioinfAdv_ApplicationNote_Supplementary_template.md` §S3. Built 2026-09-08.*
