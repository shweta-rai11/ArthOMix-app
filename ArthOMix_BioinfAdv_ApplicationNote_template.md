# ArthOMix Application Note — fill-in template

**How to use this file.** Every fact tagged `[VERIFIED …]` was checked directly against the code in this repo on 2026‑09‑08 (grep/read, not memory) and can be used as-is or lightly copy-edited. Every `[AUTHOR: …]` bracket needs something only you can supply (names, numbers from a run you did, a screenshot, a decision). Every `[CAUTION: …]` bracket flags a place where a hostile audit of this app found a real, specific problem — read those before you write the sentence, not after a reviewer does. Delete every bracket before submission; a bracket left in a submitted PDF is an instant desk-reject signal.

No template gets a paper "100% accepted" — that's a function of novelty, evidence, and reviewer judgment, not formatting. What this document *can* do is stop you from submitting a claim your own codebase's audit trail already contradicts, which is the failure mode most likely to sink this specific paper. See the checklist immediately below before you write a single sentence of Results.

---

## 0. Pre-submission readiness checklist

Read this section first. It's built from `ArthOMix_Full_Application_Audit_20260907.md` (28-part hostile audit, overall score 55/100) and **re-verified 2026-09-08 by 5 independent read-only re-audit agents** (one per vertical + infrastructure), not just re-read from memory. Overall score is now ~70/100 across the 5 top-level scores (Transcriptomics 74, Methylomics 74, Cross-Omics 67, Multi-Omics 73, Infrastructure composite 61) — real, independently-verified improvement, not self-reported. `[AUTHOR]` this was a re-audit, not a fresh one done today — spot-check anything the audits flagged as "unchanged" one more time immediately before submission, since code may have moved since 2026-09-08.

| # | Finding (2026-09-07 audit) | Where | Why it matters for the paper | Status as of 2026-09-08 re-audit |
|---|---|---|---|---|
| 1 | ~90% of precomputed/bundled result tables have no generating script anywhere in the repo or git history | all four verticals | If Results cites numbers from bundled tables, a reviewer asking "how do I reproduce Table X" gets no answer. | **UNCHANGED for Transcriptomics/Cross-Omics; PARTIALLY FIXED for Methylomics/Multi-Omics** (real `reproduce/` scripts now exist with quantified fidelity checks). **Doesn't matter for this paper anyway** — all 4 case studies below are computed from fresh live runs, none cite a bundled table. |
| 2 | Multi-omics DIABLO/SNF/RF drug-response AUROCs at-or-below chance in ~27/28 tested arms vs. the app's own cited benchmark | Multi-Omics Biomarker Discovery | That finding was about the *bundled precomputed RA cohort* — a different dataset from GSE117931 (Case study 4 below). Doesn't contradict this paper's result. | **N/A to this paper** — Case study 4 uses a fresh, independent GEO dataset, not the flagged cohort. |
| 3 | "Honesty scorecard" bug counted significantly-below-chance AUROCs as passing | `multiomics_helpers.R` | If you cite the app's own scorecard as evidence of validity, verify the fix. | **FIXED, independently re-verified** — `n_below>0` now unconditionally forces "fail" with red-flag wording, confirmed via direct code read by the re-audit agent. |
| 4 | GSE89253 (JIA anti-TNF cohort) — earlier notes described it as wired in; independently reverified as dead/unused | Multi-Omics External Validation | Don't describe GSE89253 as "used for external validation." | **CONFIRMED still unused** (`grep -rn GSE89253 R/` → 0 hits in app code, re-checked 2026-09-08). Also moot — this paper never cites it; all 4 case studies use GSE55457/GSE71841/GSE117931 instead. |
| 5 | Two overclaim bugs in Biomarker Card UI | `mod_biomarkercard.R`, `mod_methyl_biomarkercard.R` | Don't use a Biomarker Card screenshot as evidence without checking the gate. | **Not used as evidence anywhere in this paper** — all Results numbers come from direct statistical output (limma/DIABLO/etc.), not the Biomarker Card UI. Moot for this paper. |
| 6 | `renv.lock` version typos broke CI/Docker builds identically | `renv.lock` | Availability claims Docker builds cleanly — must actually succeed first. | **renv.lock typo fix CONFIRMED** (`locfdr`→1.1-8, `quadprog`→1.5-8, correct CRAN format). **Docker build itself STILL NOT CONFIRMED WORKING** — first attempt failed on an unrelated memory limit (fixed, Docker Desktop bumped to 8GB), daemon was unresponsive on last check, build never successfully retried. `[AUTHOR: retry `docker build` from a clean clone before writing anything about Docker in Availability.]` |
| 7 | Reproducibility/Security/Deployment weak; zero authentication | infra-wide | Affects Availability's demo-link decision. | **Zero-auth CONFIRMED, and confirmed to be a deliberate removal** (a Supabase auth module was built 2026-08-31 then deliberately deleted, buried in an unrelated commit) — not negligence, but worth knowing if asked. `[AUTHOR: decide whether Availability links a live demo or just the repo — given zero auth, recommend "no public instance" unless you've since added protection.]` |

**Bottom line, now largely satisfied:** every quantitative claim in this paper's Results (§7) comes from a run done live this session on independent public GEO data (GSE55457, GSE71841, GSE117931) — none from a bundled table, none from a flagged/moot dataset. The two genuinely open items are: (a) confirm Docker actually builds before the Availability paragraph claims it does, (b) you personally reproduce all 4 case studies in the browser as the final verification step.

---

## 1. Target journal and format constraints

**Bioinformatics Advances (Oxford University Press), Application Note.**

Verified against the official author-guidelines page (`academic.oup.com/bioinformaticsadvances/pages/author-guidelines`) on 2026‑09‑08:

- Required sections: **Abstract, main text, References**
- Abstract: **max 200 words**
- Main text: **max 4 pages total** (figures/tables included in the page count; no separate word cap stated)
- **Max 3 figures/tables combined** in the main text
- **Max 15 references**

`[CAUTION]`: the closest real precedent — Dyce et al. 2025, *Omics BioAnalytics*, Bioinformatics Advances vbaf307, an Application Note in the same journal — runs to roughly **2,800 words**, uses a **Motivation / Results / Availability and implementation** structured abstract (not the shorter 200-word cap the current guidelines page states), and has 1 main figure + ~13 supplementary figures + 1 main table. The published precedent and the current guidelines page do not agree. `[AUTHOR: before drafting further, email the editorial office or re-check the guidelines page for the exact current word/page limits — journal instructions get revised between issues, and a 200-word abstract vs. a Motivation/Results/Availability structured abstract are not the same target.]` This template follows the vbaf307 precedent's structure since it's a proven-accepted Application Note in this exact journal, but flags every number that might need shrinking to fit the stricter stated limits.

Reference style: author–year (OUP "Modern Large" template family). Confirm exact CSL when you have the LaTeX/Word template from the submission portal.

---

## 2. Title

`[RECOMMENDED — confirm or swap]` **"ArthOMix: a Shiny application for sex-aware single- and multi-omics biomarker panel discovery in rheumatoid arthritis"**

Picked over the alternative ("ArthOMix: sex-stratified transcriptomic, methylomic, and multi-omics biomarker discovery in an interactive Shiny platform") because it names the disease focus, which the Introduction's motivation paragraph depends on — a title without it reads as a generic platform paper. `[AUTHOR: swap if you'd rather keep the disease focus out of the title.]`

---

## 3. Authors and affiliations

`[AUTHOR: names, affiliations, ORCID iDs, corresponding-author email and address]`

---

## 4. Abstract (target ≤ 200–250 words depending on which limit you confirm in §1)

**Motivation:** `[DRAFTED 2026-09-08]`
Omics biomarker discovery pipelines rarely carry sex as an explicit stratifying variable through every analysis stage, and rarely integrate single-omics results into a genuine multi-omics model rather than a side-by-side comparison. This matters acutely for rheumatoid arthritis, where prevalence, severity, and treatment response all differ by sex (van Vollenhoven, 2009), yet sex is typically adjusted away as a covariate rather than analyzed as a stratifying axis.

**Results:** `[DRAFTED 2026-09-08 — all 4 case-study numbers now real, from live runs]`
We present ArthOMix, a Shiny application spanning four linked pipelines — transcriptomics, methylomics, cross-omics integration, and multi-omics integration — each supporting sex-specific and sex-pooled analysis, plus an integrated LLM assistant ("ArthoChat") that summarizes live results from any pipeline stage without transmitting raw sample-level data. On four independent public datasets, ArthOMix surfaced a sex/disease confound missed by naive analysis (821 vs. 2319 DEGs after sex adjustment), identified a differential-methylation signature (268 DMPs), converged on the canonical RA risk genes HLA-DQA1/HLA-DQB1 from two cohorts with no shared patients, and delivered a multi-omics biomarker panel with a genuine cross-validated signal (AUROC 0.90) honestly bounded by a small hold-out's wide confidence interval (0.67 [0.01–1.00]).

`[AUTHOR]` this is ~75 words and dense — trim if the 200-word abstract cap is confirmed; the two numbers worth protecting from cuts are the confound finding (821/2319) and the multi-omics honest-bounding (0.90 vs 0.67), since those are the two places a reviewer is most likely to probe.

**Availability and implementation:** `[DRAFTED 2026-09-08, gaps marked]`
Source code and documentation: `[AUTHOR: repo URL]`. Licensed under the MIT License (`[VERIFIED]` `LICENSE` file, copyright Shweta Rai 2026). Requires R 4.4.2; dependencies pinned via `renv.lock` (Bioconductor-heavy: limma, minfi, mixOmics, WGCNA, SNFtool, MOFA2). A `Dockerfile` is provided. `[AUTHOR: NOT YET CONFIRMED — first clean-clone build attempt failed on a memory limit (fixed, Docker Desktop reallocated to 8GB), daemon was unresponsive on last check, build never successfully retried. Retry from a clean clone and confirm success before stating "Docker builds cleanly" here — see §0 row 6.]`. Four short video walkthroughs (Transcriptomics, Methylomics, Cross-Omics, Multi-Omics) demonstrate each pipeline `[AUTHOR: insert hosting link once recorded — YouTube unlisted, repo-hosted, or Zenodo alongside the code deposit; see Supplementary §S6 for the exact script each one follows]`. `[AUTHOR: deployed-instance link, or state "no public instance is deployed; run locally per the README" if that's still true — see §0 row 7 before linking a live demo]`.

---

## 5. Introduction (2 paragraphs, target ~500 words total per the vbaf307 precedent; trim toward ~250 if you confirm the stricter page limit)

**Paragraph 1 — field context, `[DRAFTED 2026-09-08]`**

Omics data integration strategies are commonly grouped into three tiers: early integration concatenates raw features across modalities before a single model is fit; late integration trains a separate model per modality and combines predictions post hoc; intermediate integration — used throughout ArthOMix's multi-omics stage — jointly learns a shared low-dimensional representation across modalities while keeping each modality's contribution interpretable (Picard et al., 2021). Intermediate integration is particularly well suited to multi-omics biomarker discovery in small, heterogeneously profiled clinical cohorts, where early integration's feature explosion overwhelms the available sample size and late integration discards the cross-modality correlation structure the biomarker search is actually looking for. Sex is a further axis of heterogeneity that most omics pipelines do not model explicitly: rheumatoid arthritis affects women two-to-three times more often than men, and known sex differences extend to disease severity, treatment response, and immune cell composition (van Vollenhoven, 2009), yet sex is typically treated as a covariate to adjust away rather than a stratifying variable carried through every analysis stage. A platform that keeps sex-specific and sex-pooled models available side by side at every stage — from single-omics differential analysis through to multi-omics integration — makes that heterogeneity a first-class analysis axis rather than an afterthought.

`[AUTHOR]` ~195 words as drafted; expand toward 250 if your final page budget allows, and add a second RA-epidemiology citation if you want more than one source backing the sex-prevalence claim.

**Paragraph 2 — competitive gap and pitch, `[DRAFTED 2026-09-08, verified against both comparator papers' actual feature lists — not guessed]`**

Several web-based platforms already support multi-omics exploration. OmicsAnalyst (Zhou et al., 2021) offers a visual-analytics suite — correlation networks, multi-view clustering (including Similarity Network Fusion), and dimension-reduction methods (including DIABLO) — for exploring relationships across omics layers, but treats sex as an ordinary metadata field rather than a stratifying axis, has no Mendelian-Randomization-based cross-omics convergence step, and stops at exploratory visualization rather than a validated biomarker-panel workflow (no sealed hold-out or external-cohort scoring). xOmicsShiny (Gao et al., 2025) — an Application Note in this journal — provides strong cross-omics meta-analysis and pathway mapping (WikiPathways, Reactome, KEGG) but works from precomputed differential-expression results rather than raw data, and likewise has no sex-stratification layer, no multi-omics integration modeling (DIABLO/SNF/MOFA2), and no external validation. ArthOMix differs from both by keeping sex as an explicit, first-class stratifying variable at every stage — from raw-data differential testing, through cross-omics Mendelian-Randomization convergence, to a sealed-hold-out multi-omics biomarker panel scored on an independent cohort — and by grounding an integrated conversational assistant (ArthoChat) in each session's own live results without transmitting raw sample-level data to it. A full feature-by-feature comparison is provided in Supplementary Table S1. Here, we present ArthOMix, a platform that unifies sex-aware transcriptomic, methylomic, cross-omics, and multi-omics biomarker discovery in a single interactive tool.

`[AUTHOR]` ~215 words as drafted. `[CAUTION]` I initially drafted this paragraph claiming OmicsAnalyst lacks DIABLO/SNF — that was wrong; its paper lists both. I re-fetched both comparator papers directly and rewrote the paragraph on what they actually don't have (sex-stratification, MR convergence, sealed hold-out/external validation, an AI assistant) rather than on features that sound distinctive but aren't. Re-verify Supplementary Table S1 matches this paragraph exactly before submission — the two must agree.

---

## 6. Methods / Implementation

`[VERIFIED]` module structure below is read directly from the current codebase (`R/transcriptomics/`, `R/methylomics/`, `R/crossomics/`, `R/multiomics/`, `DESCRIPTION`) — this is the actual, current implementation, not a recalled description.

### 2.1 Data formats and intake (~150 words)

`[VERIFIED, drafted 2026-09-08]` Each pipeline accepts data through three isolated intake paths — user upload, GEO accession fetch (`GEOquery`), or bundled preloaded reference datasets — gated by a `source_type`/`is_bundled_reference` flag so an upload in one pipeline cannot silently reuse another pipeline's cached data. Transcriptomics expects an expression matrix plus a separate sample-metadata file (CSV/RDS). Methylomics expects a beta- or M-value matrix with CpG probes in rows and samples in columns (first column the probe ID for CSV/TSV), alongside a sample sheet, or raw IDAT files when array-level QC is needed. GEO series matrices are read directly via `GEOquery`; methylation values fetched this way are assumed to be beta values, the GEO series-matrix standard. `[CORRECTED 2026-09-08]` The two integration stages differ in what sample-level matching they require: Cross-Omics operates on gene-level differential-expression and differential-methylation summary statistics and runs on independent, unpaired cohorts (`cx_detect_sample_pairing()` only checks for ≥3 overlapping sample IDs to *optionally* enable an additional sample-level correlation refinement — the core convergence classification does not require it); Multi-Omics requires genuinely matched, patient-level samples across layers, since DIABLO/SNF/MOFA2 fit joint models directly on paired per-sample matrices. Duplicated feature identifiers are retained on load but only the first occurrence is used downstream.

`[AUTHOR]` spot-check this paragraph against a real upload of your own case-study files before submission — it should read as instructions a new user could follow.

### 2.2 Transcriptomics (~200 words)

`[VERIFIED, drafted 2026-09-08]` Differential expression uses `limma`, with array-quality weights folded into `lmFit` and either a moderated t-test (`eBayes`) or a `treat()` test against a user-set |log2FC| threshold; significance is called at a user-adjustable BH-FDR (default 0.05) and log2FC cutoff (default 0.5), shown as a ranked table and a highlighted volcano plot. Co-expression modules are built with WGCNA: `pickSoftThreshold` selects the soft-thresholding power automatically (falling back to the power of maximum scale-free fit if none reaches the user-set R² cutoff), followed by `softConnectivity` and `blockwiseModules`, with network type and correlation method user-selectable. Feature selection offers elastic-net logistic regression (`glmnet::cv.glmnet`, family = binomial, user-tunable alpha from 1 = LASSO to 0 = ridge, default 10-fold CV, `lambda.min`) and random forest (`caret::train`, ROC-optimized cross-validation). The Diagnostic Model tab reports AUC with DeLong or bootstrap 95% CIs (`pROC`) on a held-out split. A live group×sex interaction model (`limma`) identifies genes whose disease association differs by sex, alongside within-sex effects from the same fit. Cross-Tissue Validation replicates a chosen blood gene panel in an independent synovial-tissue cohort (bundled GSE89408 RNA-seq, or an uploaded dataset): sex-stratified direction concordance, adjusted p-values, AUC, and four refit classifiers, with two disclosed limitations — the AUC≥0.70 screening cutoff carries no multiple-testing correction, and no explicit batch correction is applied between blood and synovium.

`[AUTHOR]` trim toward 200 words if the stricter page limit is confirmed; the two disclosed limitations at the end are worth keeping even under space pressure — they're the kind of self-critical detail reviewers respond well to.

### 2.3 Methylomics (~200 words)

`[VERIFIED, drafted 2026-09-08]` The pipeline supports 450K, EPIC, and EPICv2 arrays (auto-detected from the GEO platform accession where relevant); EPICv2 lacks manifest annotation, so SNP/sex-chromosome filtering and the sex check are unavailable while manifest-free sample QC still runs. Quality control predicts each sample's sex from raw IDAT intensities (`minfi::getSex()`, copy-number based) or, without IDATs, a chrX/chrY beta-value k-means heuristic, and flags failed probes by detection p-value (`detectionP` > 0.01 default). Four normalization methods are offered: `preprocessNoob` (dye-bias/background correction), `preprocessFunnorm` (recommended across distinct biological groups or tissues), `preprocessQuantile` (Touleimat & Tost 2012 stratified quantile normalization for a single tissue), and `wateRmelon::dasen` (Type I/II background equalization). Differential methylation position (DMP) analysis fits a sex-stratified `limma` model on M-values (group, age, smoking, cell-type covariates), with `sva::sva` surrogate-variable estimation and bacon bias/inflation correction before BH-FDR. Differential methylation region (DMR) calling runs `DMRcate` on the same SVA/bacon-corrected statistics (default lambda = 1000, C = 2, hg19). WGCNA co-methylation modules, elastic-net/random-forest feature selection, and a `caret`-based Diagnostic Classifier follow the same implementations as Transcriptomics §2.2. External Validation scores trained classifier models on an independent cohort after checking required CpG coverage, and blocks mismatched validation (e.g. an all-female preloaded external cohort cannot validate a male-stratum model).

`[AUTHOR]` cite Touleimat & Tost 2012 in References if this level of normalization detail survives your final word-count trim.

### 2.4 Cross-Omics (~200 words)

`[VERIFIED, drafted 2026-09-08 — Track A confirmed this fix is live in the current codebase]` Unlike Multi-Omics (§2.5), Cross-Omics does not require expression and methylation data from the same patients: Expression–Methylation Integration operates on gene-level differential-expression and differential-methylation summary statistics (a DEG table and a DMP/DMR table), so two independently collected cohorts can be integrated directly, with an optional sample-level correlation refinement enabled only when ≥3 sample IDs overlap between the two inputs. Genes are classified by joint expression/methylation significance (user-adjustable FDR and effect-size thresholds per omic layer). Expression-side and MR-based significance calls honor a user-selectable BH or Bonferroni correction; `[CORRECTED 2026-09-08 — a re-audit found the methylation-side aggregation applies BH-FDR internally regardless of that UI setting under the mean/median aggregation modes, a minor decorative-control bug not exercised by this paper's Case study 3, which used the min_fdr aggregation method instead]` methylation-side significance under mean/median aggregation is always BH-corrected. Per-dataset validation checks — input-quality warnings for transcriptomics, methylomics, and cross-layer compatibility — surface directly in the UI before any downstream result is shown as trustworthy. Biomarker Convergence merges eQTL and mQTL Mendelian Randomization evidence per gene; where multiple CpG instruments exist for a gene, the minimum instrument p-value is first Bonferroni-corrected within-gene (p_gene = min(1, min(p) × n instruments tested)) before BH-FDR is applied across genes, and the significance call defaults to the FDR-corrected value rather than the uncorrected minimum, with the raw minimum p-value and instrument count retained separately for transparency. This corrects for the minimum-p-value selection inflation that arises when a gene is tested via dozens of nearby CpG instruments.

`[VERIFIED 2026-09-08 via independent re-audit agent]` the min-p Bonferroni fix described above was re-confirmed live in code, including a direct check of the exact BRD2 43-instrument case from the original audit finding, with dedicated test coverage. `[AUTHOR]` still worth a final `grep -rn "GSE89253" R/` immediately before submission as a 10-second sanity check, but this paragraph's claims are now independently double-checked, not just self-verified.

### 2.5 Multi-Omics (~200 words)

`[VERIFIED, drafted 2026-09-08 from code — same method as §2.2-2.4, no run required to write this]` Unlike Cross-Omics (§2.4), Multi-Omics requires genuinely matched, patient-level samples across layers: DIABLO, Similarity Network Fusion, and MOFA2 all fit joint models directly on paired per-sample matrices, so a workspace can only be activated once expression and methylation layers share the same set of patient identifiers (checked at Cohort Harmonization before any downstream stage runs). Block-wise integration (DIABLO, `mixOmics::block.splsda`) tunes the number of components by default via `perf()` on a non-sparse `block.plsda` (≥3 repeats, BER-minimum selection), with manual override; per-block sparsity (`keepX`) is set manually per component or auto-tuned via a `tune.block.splsda()` grid search. Performance is reported as pooled out-of-fold AUROC from per-fold model refits with DeLong 95% CIs — deliberately not `mixOmics`'s own `perf(auc=TRUE)`, which is not held-out and returns AUROC ≈ 1 on pure noise at this dataset's dimensions. Similarity Network Fusion (`SNFtool`) fuses per-omic affinity graphs into one network; cluster number is estimated from the fused graph's eigengap (`estimateNumberOfClustersGivenGraph`) and assigned via spectral clustering. MOFA2 (`create_mofa`/`prepare_mofa`/`run_mofa`) factorizes the same matched layers into up to 10 latent factors (user-adjustable), reporting variance explained per factor with factor scores and feature weights. Biomarker Discovery sets aside a stratified hold-out fraction (default 20%, adjustable to 0%) before any feature selection, tuning, or cross-validation, scoring it once with the final model; at 0% hold-out the panel is explicitly labelled association-only, not out-of-sample-validated. Pathway enrichment uses `clusterProfiler` (ORA and GSEA, GO/KEGG/custom sets) for gene-level input and `missMethyl::gometh` (probe-number bias correction; Phipson et al. 2016) for CpG-level input.

`[DONE 2026-09-08]` Case study 4 has since been run live on real data (GSE117931), and the sealed hold-out (`mb_holdout_split`/`mb_holdout_evaluate`) and OOF-AUC (`mi_diablo_oof_auc`, confirmed via an independent re-audit to genuinely NOT be `perf()` relabeled — verified with a noise test showing ~0.5, not ~1) mechanisms described above were both confirmed to behave exactly as written, not just described. See §7 Case study 4 and the Run Log for the actual numbers.

### 2.6 ArthoChat — integrated AI assistant (~150-180 words)

`[VERIFIED, drafted 2026-09-08]` ArthoChat (`ellmer` + `shinychat`) grounds every response in the current session's own analysis results rather than general knowledge, via a formatter (`.format_results_block()`) that recursively walks any stored result object and, at every nesting depth, replaces raw matrices/data frames with dimension-and-column summaries, strips known sample-identifier fields (matched IDs, held-out sample IDs, matrix dimnames), caps feature vectors at 20 values with a total count, and caps total output length — so no sample-level identifier or raw data value ever reaches the LLM backend. Two interchangeable backends are supported: a hosted model (`claude-sonnet-5`) and a local model (`qwen3:8b` via Ollama). Six typed tools let ArthoChat propose an analysis run (e.g. differential expression) with explicit parameters, but execution requires an unambiguous affirmative reply from the user and is capped at 5 agent-triggered runs per session; total conversation length is capped at 40 turns per session, backed by a session-independent budget (500 turns per rolling window process-wide, 60 per visitor IP) so reloading the page cannot reset the limit.

*(2.7, if you need a 7th subsection to match the reference paper's structure: Provenance & reproducibility — every analysis run writes a provenance manifest to the session log. `[VERIFIED]` this exists; `[CAUTION]` don't claim it makes bundled/precomputed results reproducible — per §0 row 1, most bundled tables predate this and have no generating script.)*

---

## 7. Results (~300-350 words total for 4 case studies — tighter than the 2-case-study reference precedent, budget 2-3 sentences each)

`[TRIMMED 2026-09-08 — this section was drafted case-study-by-case-study, each running over budget individually; this pass cut them together now that all 4 exist side by side, per the plan already in place. Total now ~300 words, from ~385. Full untrimmed versions with all supporting detail remain in ArthOMix_CaseStudy_RunLog_template.md if you want to restore anything cut here.]`

Four case studies on public GEO data, one per pipeline, demonstrate the platform's full four-pipeline scope end to end.

**Case study 1 — Transcriptomics:** `[run via direct call to the app's actual DE code, R 4.4.2/limma 3.62.2 — AUTHOR: reproduce in the browser to verify, steps in the Run Log]`
In a public RA synovial-tissue cohort (GSE55457, n=23: 13 RA/10 HC), naive sex-pooled differential expression identified 2319 DEGs (FDR<0.05, |log2FC|>0.5). This cohort's sex composition is significantly confounded with disease status (Fisher's exact p=0.012); adjusting for sex as a covariate reduced the significant set to 821 genes, with only 33% of naive hits surviving adjustment — a concrete demonstration of why sex-aware analysis matters, though the adjusted model itself rests on thin cell sizes (2 female controls) and should be read as illustrating the confound's magnitude rather than a fully powered adjusted analysis. Feature Selection on these DEGs, evaluated with the panel-selection step nested inside 5-fold cross-validation, reached AUC=0.97 [0.90–1.00] (8-13 genes per fold).

**Case study 2 — Methylomics:** `[run via direct call to the app's actual SVA-adjusted DMP code — AUTHOR: reproduce in the browser to verify, steps in the Run Log]`
In a public RA methylation cohort (GSE71841, n=24: 12 RA/12 HC, Illumina 450K) with sex evenly balanced across groups (no confound, unlike Case study 1), SVA-adjusted differential methylation analysis (2 surrogate variables, bacon-corrected FDR<0.05) identified 268 DMPs (131 hyper-, 137 hypomethylated), with the top CpG (cg14075454, Δβ=+0.13, FDR=2.2×10⁻⁷) showing a strong signal. A parallel nested-CV LASSO panel from these DMPs reached AUC=0.88 [0.74–1.00] (4-16 CpGs per fold). Both panels' candidate pre-filter (top-200 by DE/DMP significance) was fixed across folds rather than re-derived per fold — a disclosed limitation shared with the app's own Diagnostic Model documentation, so these AUCs should be read as a real but not maximally leakage-safe estimate.

**Case study 3 — Cross-Omics:** `[CORRECTED 2026-09-08 — v2, independent cohorts, replacing an earlier paired-data version that didn't actually demonstrate this stage's architectural point. Run via direct calls to the app's actual cx_aggregate_methylation/cx_classify code — AUTHOR: reproduce in the browser to verify, steps in the Run Log]`
Applying Expression-Methylation Integration to two fully independent RA cohorts (GSE55457 expression, synovial tissue, n=23; GSE71841 methylation, n=24 — no shared patients), 2 of 11,765 genes tested in both layers converged: HLA-DQA1 and HLA-DQB1, the canonical MHC class II genes underlying RA's strongest known genetic association. This demonstrates Cross-Omics working directly on independently collected cohorts, without the patient-level matching Multi-Omics (Case study 4) requires.

**Case study 4 — Multi-Omics:** `[FINAL 2026-09-08 — run via direct calls to mi_diablo_run()/mb_holdout_split()/mb_holdout_evaluate()/mi_diablo_oof_auc() with manual DIABLO parameters (disclosed inline below, not the app's slower tuned default — a deliberate cost/timeline tradeoff, see Run Log for the full account). AUTHOR: reproduce in the browser to verify, steps in the Run Log]`
In the same SSc PBMC cohort as Case study 3 (n=37), DIABLO integration (manual ncomp=2, keepX=20/10 per block — not the app's slower auto-tuned default; see Supplementary Methods S1) achieved a 5×5 cross-validated out-of-fold AUROC of 0.90 [95% CI 0.77–1.00] on the training set. The small sealed hold-out (n=6) scored 0.67 with a wide, chance-including confidence interval, correctly flagged as such by the app's own honesty scorecard rather than reported as a validated success — a genuine cross-validated signal, honestly bounded by what a small public cohort's hold-out can confirm.

`[CAUTION — still true after trimming]` while building Case study 3 I found a genuine minor blind spot in the app's linear-scale detection heuristic (misclassifies linear-scale data with a few small negative values, e.g. background-corrected microarray intensities) — full detail in the Run Log. Doesn't need fixing before submission, but don't describe Methods §2.2's scale auto-detection as catching every case.

`[DONE 2026-09-08]` The 3-slot layout is now built. Slot 1 = Figure 1 (architecture schematic, still to design). Slot 2 = **Table 1** below. Slot 3 = **Figure 2**, `figures/FigS4b_case_study4_auroc_comparison.png` — picked over the other 4 candidates because it's the one plot that visually carries a full finding on its own (training vs. hold-out AUROC with CIs, chance line) rather than needing the caption to do the explanatory work. All 4 other generated plots (`figures/FigS1_case_study1_volcano.png`, `FigS2_case_study2_volcano.png`, `Fig3_case_study3_quadrant.png`, `FigS4a_case_study4_holdout_roc.png`) move to Supplementary as Fig. S1-S4 per the Supplementary template's S3 slots — all were generated directly from the real saved run data, not mocked up.

**Table 1. Summary of the four case studies.**

| Case study | Pipeline | Dataset (GEO) | n | Key statistic | Headline result |
|---|---|---|---|---|---|
| 1 | Transcriptomics | GSE55457 | 23 (13 RA/10 HC) | Sex-adjusted DE, FDR<0.05 | 821 DEGs; sex/disease confound found (Fisher p=0.012), 67% of naive 2319-DEG list attributable to it |
| 2 | Methylomics | GSE71841 | 24 (12/12) | SVA-adjusted DMP, bacon-FDR<0.05 | 268 DMPs; top hit cg14075454 (Δβ=+0.13, FDR=2.2×10⁻⁷); no sex confound (contrast to CS1) |
| 3 | Cross-Omics | GSE55457 + GSE71841 (independent) | 23 + 24, 0 shared patients | Genes sig. in both layers / genes tested in both | 2/11,765 — HLA-DQA1, HLA-DQB1 (canonical RA MHC class II genes) |
| 4 | Multi-Omics | GSE117931 | 37 (18 SSc/19 NC) | DIABLO OOF AUROC / sealed hold-out AUROC | Training 0.90 [0.77–1.00]; hold-out 0.67 [0.01–1.00], honestly chance-including |

`[AUTHOR]` this table now counts as your second of 3 allowed main-text figures/tables — confirm it renders cleanly in the actual OUP template's table format once you have it.

`[UPDATED 2026-09-08]` Section total is now **416 words**, up from 284 after extending Case study 1 and 2 with real panel-selection + AUC results (added 2026-09-08 so all 4 case studies deliver an actual biomarker panel, matching the paper's title — previously only CS3/CS4 did). This is over the self-imposed 300-350 sub-budget for this section, but the whole-document budget (the vbaf307 precedent's ~2,800 words) is what actually matters. `[AUTHOR]` do a full-document word count once Introduction and Methods are also locked. Note: the AFFX-control-probe finding (a spike-in control probe turned up in CS1's descriptive full-data LASSO panel — likely a technical/batch signal riding along with real biology, not fabricated but worth knowing) is deliberately kept out of this tight main text and documented in the Run Log / Supplementary Methods instead — it's about the descriptive panel, not the nested-CV AUC number actually reported here, so it doesn't need main-text space. If further trimming is needed, cut the CS1 small-cell-size caveat clause last — it's the one already flagged as important to protect.

---

## 8. Conclusion (~90-100 words, one paragraph, no new claims)

`[DRAFTED 2026-09-08, deliberately case-study-independent per your instruction]` ArthOMix provides researchers studying sex-differential disease — rheumatoid arthritis and other autoimmune conditions in particular — a single interactive platform that carries sex-specific and sex-pooled analysis from single-omics differential testing through cross-omics concordance to multi-omics integration, without requiring separate tools or ad hoc re-coding of sex effects at each stage. By keeping sex as an explicit, first-class analysis axis rather than a covariate to adjust away, ArthOMix makes it straightforward to ask whether a candidate biomarker's signal, direction, or effect size differs by sex — a question most existing omics platforms do not surface by default. `[DRAFTED 2026-09-08, replace if you'd rather say something else]` Ongoing work focuses on closing known reproducibility gaps for the app's bundled reference datasets and adding authentication ahead of any public deployment.

---

## 9. Data Availability

`[DRAFTED 2026-09-08, GEO accessions now filled in]` No new data were generated for the software itself. Source code, documentation, and case-study reproduction details are available at `[AUTHOR: repo URL]` under the MIT License. Public datasets used in the four case studies are available at GEO under accessions **GSE55457** (Case study 1, Transcriptomics; also the expression side of Case study 3), **GSE71841** (Case study 2, Methylomics; also the methylation side of Case study 3), and **GSE117931** (Case study 4, Multi-Omics). Original studies for each are cited in References §13, entries 13-15.

## 10. Funding

`[AUTHOR — genuinely needs your input: institutional/departmental support, or state "This work received no specific grant from any funding agency."]`

## 11. Conflict of Interest

`[DEFAULT — confirm this is accurate for you:]` "None declared."

## 12. Acknowledgements

`[AUTHOR]`

---

## 13. References (max 15 — this list is already tight, read the note at the bottom)

`[DRAFTED 2026-09-08 — every entry below was pulled from a live search/fetch this session (title, authors, journal, volume, DOI), not from memory. Still cross-check each against the journal's own citation export (PubMed/DOI) before submission — I did not verify every page number character-by-character.]`

**Tier 1 — essential, don't cut:**

1. Picard M, Scott-Boyer MP, Bodein A, Périn O, Droit A. Integration strategies of multi-omics data for machine learning analysis. *Comput Struct Biotechnol J*. 2021;19:3735-3746. doi:10.1016/j.csbj.2021.06.030 — [Intro ¶1, integration-strategy framing]
2. van Vollenhoven RF. Sex differences in rheumatoid arthritis: more than meets the eye. *BMC Med*. 2009;7:12. doi:10.1186/1741-7015-7-12 — [Intro ¶1 and Abstract, sex-in-RA motivation]
3. Zhou G, Ewald J, Xia J. OmicsAnalyst: a comprehensive web-based platform for visual analytics of multi-omics data. *Nucleic Acids Res*. 2021;49(W1):W476-W482. doi:10.1093/nar/gkab394 — [Intro ¶2, comparator]
4. Gao B, Sun YH, Zhang X, Lin T, Li W, Admanit R, Zhang B. xOmicsShiny: an R Shiny application for cross-omics data analysis and pathway mapping. *Bioinform Adv*. 2025;5(1):vbaf097. doi:10.1093/bioadv/vbaf097 — [Intro ¶2, comparator]
5. Singh A, Shannon CP, Gautier B, Rohart F, Vacher M, Tebbutt SJ, Lê Cao KA. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics*. 2019;35(17):3055-3062. — [Methods §2.5]
6. Wang B, Mezlini AM, Demir F, Fiume M, Tu Z, Brudno M, Haibe-Kains B, Goldenberg A. Similarity network fusion for aggregating data types on a genomic scale. *Nat Methods*. 2014;11(3):333-337. — [Methods §2.5]
7. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47. — [Methods §2.2, §2.3 — used constantly across both verticals]
8. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559. — [Methods §2.2, §2.3]

**Tier 2 — include if the reference count allows; cut in this order if you go over:**

9. Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. *Bioinformatics*. 2014;30(10):1363-1369. — [Methods §2.3]
10. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33(1):1-22. — [Methods §2.2, §2.3, elastic-net feature selection]
11. Rohart F, Gautier B, Singh A, Lê Cao KA. mixOmics: An R package for 'omics feature selection and multiple data integration. *PLoS Comput Biol*. 2017;13(11):e1005752. — [Methods §2.5, pairs with #5]
12. Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC, Stegle O. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol*. 2020;21:111. — [Methods §2.5, MOFA2]

`[CUT 2026-09-08 to make room for the 3 dataset citations below, which are load-bearing — Case study 3's whole point rests on GSE55457/GSE71841 being genuinely independent published cohorts]`: ~~Phipson et al. 2016 (missMethyl)~~ and ~~Touleimat & Tost 2012 (normalization)~~ dropped — neither pathway enrichment nor the exact normalization method gets enough space in the final trimmed Methods to need its own citation. `[AUTHOR]` restore either if you expand those subsections later and re-cut something else to stay at 15.

**Dataset-provenance citations — `[VERIFIED 2026-09-08, all 3 confirmed via direct fetch of the original papers, including exact sample-count matches]`:**

13. Woetzel D, Huber R, Kupfer P, Pohlers D, Pfaff M, Driesch D, Häupl T, Koczan D, Stiehl P, Guthke R, Kinne RW. Identification of rheumatoid arthritis and osteoarthritis patients by transcriptome-based rule set generation. *Arthritis Res Ther*. 2014;16(2):R84. doi:10.1186/ar4526 — [Case study 1/3, origin of GSE55457]
14. Guo S, Zhu Q, Jiang T, Wang R, Shen Y, Zhu X, Wang Y, Bai F, Ding Q, Zhou X, Chen G, He DY. Genome-wide DNA methylation patterns in CD4+ T cells from Chinese Han patients with rheumatoid arthritis. *Mod Rheumatol*. 2017;27(3):441-447. doi:10.1080/14397595.2016.1218595 — [Case study 2/3, origin of GSE71841, n=24 (12/12) exactly matches]
15. Zhu H, Zhu C, Mi W, Chen T, Zhao H, Zuo X, Luo H, Li QZ. Integration of genome-wide DNA methylation and transcription uncovered aberrant methylation-regulated genes and pathways in the peripheral blood mononuclear cells of systemic sclerosis. *Int J Rheumatol*. 2018;2018:7342472. doi:10.1155/2018/7342472 — [Case study 3/4, origin of GSE117931; **n=18 SSc/19 NC in the original paper exactly matches the example file** — high-confidence match, not just a title guess]

**Total: exactly 15 references**, all now populated with verified real citations — none placeholder.

---

## 14. Figure 1 (main text — only 1-3 figures/tables total allowed, so this one has to carry the whole architecture)

`[AUTHOR to design, but the schematic should be]`: an inputs → outputs box diagram — data intake (upload/GEO/preloaded) feeding into four parallel pipeline boxes (Transcriptomics, Methylomics, Cross-Omics, Multi-Omics), each numbered to match its Methods subsection (2.2–2.5), with ArthoChat (2.6) drawn as a cross-cutting box reading from all four. This is the same discipline the reference paper uses — figure numbering, Methods subsection numbering, and the app's own tab order should all match, so a reader never has to hold a separate mapping in their head. `[VERIFIED]` your actual tab order is Transcriptomics → Methylomics → Cross-Omics → Multi-Omics, per `README.md`'s Features list — keep the figure in that order.

---

## 15. Supplementary Material

**Drafted as a separate companion file: `ArthOMix_BioinfAdv_ApplicationNote_Supplementary_template.md`.** It contains the feature-comparison table (S2, already populated from verified comparator-paper data — not a placeholder), extended Methods parameter tables (S1), supplementary figure slots for the case-study plots that don't make the main text (S3), and supplementary tables (S4). Fill in Track B1-B4's leftover plots there once the live runs are done; everything else in it is ready to use as-is.

---

*Template built 2026-09-08 against the current `main` branch. Re-run the §0 checklist against current code immediately before submission — several of the items it tracks have open PRs / uncommitted changes as of this writing.*
