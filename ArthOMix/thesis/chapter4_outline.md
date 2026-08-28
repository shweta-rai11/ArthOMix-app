# Thesis Chapter Outline: An Integrated, Reproducibility-Audited Multi-Omics Analytical Framework for Rheumatoid Arthritis Anti-TNF Response

**Status:** Outline / structural draft — Phase 1 of 2 (per your instruction, no full prose yet).
**Grounding:** Every module/submodule/method named below is verified against the actual ArthOMix codebase (`ArthOMix/R/`) and `Research_05_multiomics_sexstratified/AUDIT.md`, not invented. Every biological/statistical *result* is marked `[INSERT RESULT]` unless it is quoted directly from `AUDIT.md`, in which case it is cited as such — those are real, already-computed, independently-audited numbers from your own pipeline, not fabrications.

---

## 0. Working title options

1. "A Reproducibility-Audited Web Framework for Multi-Omics Integration in Rheumatoid Arthritis Anti-TNF Response"
2. "From Single-Omics Discovery to Cross-Platform Convergence: An Integrated Analytical Architecture for RA Multi-Omics"
3. "Design and Validation of a Modular Multi-Omics Framework: Transcriptomic, Epigenomic, and Cross-Omics Integration in Anti-TNF-Treated RA"

**Recommendation:** Option 1 — it foregrounds the thing that actually differentiates this work from Quickomics/OmicsAnalyst/etc.: the audited, leakage-checked validation layer (§9), not just "another Shiny app."

---

## 1. How "Module 1–4" maps onto the real codebase

This is the backbone decision, confirmed with you: the four chapter modules are ArthOMix's four top-level analytical areas, each a real folder under `ArthOMix/R/`, each with its own dataset ingestion, config/ui/server-registered submodules, and shared reactive results object read by the app's AI-context builder (`build_assistant_context()`, `ArthOMix/R/submodules_registry.R:117`).

| Chapter section | Codebase area | Registry list | # submodules | Maturity |
|---|---|---|---|---|
| Module 1 | `R/transcriptomics/` | `TX_MODULES` | 16 | Fully built |
| Module 2 | `R/methylomics/` | `MX_MODULES` | 12 | QC/Normalization/DMP built; Cell-type→Diagnostic partially scaffolded — **flag in chapter as a limitation, not hidden** |
| Module 3 | `R/crossomics/` | `CX_MODULES` | 3 (Integration built; Biomarker Convergence, MR Stage scaffolded) |
| Module 4 | `R/multiomics/` | `MULTI_MODULES` | 6 registry submodules + Dataset Workspace's embedded live MOFA2 | Fully built, most methodologically dense — gets the extra depth the template calls for |

This ordering is also the one biological logic already imposes: single-omics discovery (1–2) → pairwise cross-omics linkage (3) → joint multi-omics integration (4). §3's workflow diagram should make this a strict pipeline, not four parallel silos.

---

## 2. Section-by-section outline

### 1. Introduction
**Purpose:** Motivate why RA anti-TNF non-response is a multi-omics problem, and why the *analytical infrastructure* (not just the algorithms) is itself a research gap.
- Background: RA as a heterogeneous, partially sex-dimorphic autoimmune disease; anti-TNF non-response rate [INSERT CITATION — typically cited 30–40%]; heterogeneity as the reason single-omics predictors underperform.
- Challenges in multi-omics: dimensionality, platform heterogeneity, sample non-overlap, batch effects, the "curse of integration" (block methods overfitting on small n).
- The reproducibility problem specifically: most multi-omics web tools present a result with no leakage audit trail — the reader cannot tell whether reported performance survived proper train/test separation. This is the concrete gap your framework closes (§9 is the receipt).
- Existing web platforms: named here, detailed in §2.
- Research gap (see §2's Positioning subsection for the fully argued version — this is the compressed intro version).
- Overall objective: design, implement, and *validate the validation of* a modular framework spanning single-omics discovery through joint multi-omics integration for one concrete disease cohort (Tao et al. 2021 RA anti-TNF), with an explicit, auditable leakage-safety layer.
- Specific objectives (map 1:1 onto Modules 1–4 + cross-module + validation):
  1. Build reusable single-omics transcriptomic discovery/validation pipeline (Module 1).
  2. Build parallel methylomic pipeline with array-specific QC/normalization (Module 2).
  3. Quantify gene–CpG regulatory concordance directly (Module 3).
  4. Perform supervised (DIABLO) and unsupervised (SNF, MOFA2) multi-omics integration, patient stratification, and pathway-level convergence (Module 4).
  5. Independently audit every stage for train/test leakage and report leakage-corrected performance (§9) — the methodological contribution most existing platforms don't attempt.
- Contributions: bullet list, cross-referenced to §12's Discussion "Methodological contribution" for the full argument. Draft list:
  - A four-tier (single→single→pairwise→joint) modular architecture with a shared reactive results contract (`multi_results`) enabling any later module to consume any earlier module's typed output.
  - A dual-track design separating **precomputed, independently audited** pipeline results from a **live, data-adaptive** computational track (DIABLO/SNF/MOFA2) operating on user-uploaded data — see [[project-arthomix-overview]].
  - A documented leakage audit (AUDIT.md) that *found and fixed* two critical CV leakage bugs and reported the corrected, chance-inclusive performance rather than the inflated pre-fix numbers — a reproducibility case study in its own right (§9.1).

### 2. Existing Web-Based Omics Analysis Landscape
**Purpose:** Position the framework against the field using the five named platforms as *methodological* reference points, not a feature checklist.
- Single-omics platforms (brief: GEO2R, iDEP, DEBrowser class of tools) — establishes the baseline your Module 1/2 exceed.
- Multi-omics platforms — framing paragraph on why integration tools trail single-omics tools in maturity.
- Per-platform subsections (Omics BioAnalytics, Quickomics, BRIDGE, XOmicsShiny, OmicsAnalyst): for each — intended scope, integration method(s) used, what it validates, what it doesn't. **[INSERT CITATION]** for each platform's primary paper — do not paraphrase from memory without the source in hand; I have not verified these five platforms' current feature sets against their papers, so flag every specific claim about them for citation-check.
- Comparison table (see §"Tables" below) across: data input/preprocessing, QC, statistical analysis, multi-omics integration method(s), biomarker discovery, patient stratification, pathway analysis, visualization, validation/reproducibility.
- Limitations/gaps: (a) integration method is typically singular and fixed (no supervised/unsupervised choice); (b) validation is rarely more than a single train/test split, almost never nested CV or an explicit leakage audit; (c) precomputed-only or live-only, rarely both; (d) sex/covariate-stratified analysis is uncommon.
- Positioning: state plainly — the gap is not "no tool does DIABLO+SNF" (several do one or the other) but "no tool ships an audited leakage trail alongside its integration results, and none pair a precomputed, previously-validated pipeline with a live re-analysis track on the same architecture."

### 3. Overall Framework Architecture
**Purpose:** Establish the technical/data architecture once so Modules 1–4 don't each re-explain plumbing.
- Conceptual framework: single-omics discovery → cross-omics linkage → joint integration → cross-module convergence → validation, as a directed pipeline (Figure 1).
- System architecture: R/Shiny; per-area folder → one `mod_<area>_<name>.R` per submodule exporting `config`/`ui`/`server`; `submodules_registry.R` assembles per-area module lists (`TX_MODULES`/`MX_MODULES`/`CX_MODULES`/`MULTI_MODULES`) consumed by `ui_shell.R`; `0_load_omics_modules.R` auto-sources every area folder alphabetically.
- Data architecture:
  - Acquisition: bundled/preloaded cohort data (`data/preloaded/`) plus live GEO fetch or user upload (`mod_dataset.R`, `mod_multi_dataset.R`).
  - Preprocessing: area-specific (RNA-seq QC/TMM or limma-ready normalization for Module 1; array-specific Noob/Functional normalization/SWAN/Dasen/BMIQ/PBC/quantile for Module 2 — `mod_methyl_normalization.R`).
  - QC: per-area QC submodules plus the Overview tab's cross-module `multi_qc_scorecard()`.
  - Metadata handling / sample harmonization: `cohort_harmonization_helpers.R` (Module 4) — this is the piece that actually enforces sample-ID matching across omics layers before any integration runs; worth a subsection of its own since it's the thing that prevents a whole class of silent misalignment bugs.
  - Multi-omics compatibility: the `multi_dataset` shared reactiveValues (`$layers`, `$sample_meta`, `$active`, `$source`) is the single contract every Module 4 submodule and Module 3's cross-omics integration reads from — this is the architectural device that makes "Cross-Module Integration" (§8) more than a narrative claim.
- Statistical analysis / Integration / Biological interpretation / Validation / Reproducibility / Computational environment: one short paragraph each, cross-referencing where the full treatment lives (Modules 1–4, §8, §9).
- **Explicitly explain data flow between the four modules** (this is a named requirement in the brief): Module 1 output (DEGs, WGCNA modules, deconvolution fractions) and Module 2 output (DMPs/DMRs, cell-type-adjusted methylation) are the two inputs Module 3's gene–CpG concordance consumes; Module 3's concordant gene set is one of the candidate-feature inputs offered to Module 4's DIABLO/SNF; Module 4's DIABLO loadings and SNF cluster-associated features feed back into Module 3's "Biomarker Convergence" panel and forward into §8's cross-module prioritization. This loop (not a one-way pipeline) is the actual architecture — say so explicitly, it's the strongest structural claim in the chapter.

### 4. Module 1 — Transcriptomics
**Objective:** Genome-wide discovery of expression-level RA anti-TNF response signal, with mechanistic, causal, and cell-composition corroboration.
**Scientific questions:** Which genes differ by response/group status? Are they organized into coherent co-expression programs? Is any candidate causally implicated (not just associated) via genetic instruments? Does apparent expression difference reflect a true regulatory effect or immune cell-composition shift?
**Data/input:** bundled RA anti-TNF cohort (Tao et al. 2021) or user GEO/upload; counts or normalized matrix + sample metadata.
**Preprocessing:** `mod_preprocessing.R` / `mod_preprocessing_explore.R` — auto-detection of raw-counts vs. normalized input, driving the DGE method choice below.

Submodules (each gets full Research Question → ... → Take-Home structure in the final chapter):
1. **Overview** (`mod_overview.R`) — cohort/dataset characterization, not itself a hypothesis test; establishes the sample sizes and group balance every downstream test's power depends on.
2. **Differential Gene Expression** (`mod_dge.R`) — moderated t-test (limma, Ritchie et al. 2015) for log-scale/normalized data or negative-binomial GLM (DESeq2, Love et al. 2014) for raw counts, auto-selected by input type; `arrayWeights` (Ritchie et al. 2006) applied in the project's own pipeline. State design formula, covariate handling, BH-FDR, effect-size threshold — all **[INSERT PARAMETER]** where the live-run defaults need to be pulled from `mod_dge.R`'s config.
3. **WGCNA co-expression network** (`mod_wgcna.R`) — soft-threshold power selection (`pickSoftThreshold`), `blockwiseModules`, correlation function choice (Pearson vs. `bicor`), module–trait association. Addresses systems-level organization DGE alone cannot.
4. **Candidate gene prioritization** (`mod_candidates.R`) — how DGE + WGCNA hits are combined into a ranked candidate list; this is the hand-off artifact to Module 3/4.
5. **Mendelian Randomization** (`mod_mr.R`) — TwoSampleMR (`format_data`/`harmonise_data`) + `MendelianRandomization::mr_ivw` (random-effects IVW) — addresses causality, not just association, for candidate genes with eQTL instruments.
6. **Colocalization** (`mod_coloc.R`) — `coloc.abf`, testing whether an eQTL and RA-GWAS signal at a candidate locus share one causal variant vs. two independent ones — the causal-inference complement to MR.
7. **Feature selection** (`mod_featureselection.R`) — ML-facing dimensionality reduction distinct from DGE's hypothesis-testing framing.
8. **Diagnostic / Interaction / Cross-tissue / Cross-ancestry** (`mod_diagnostic.R`, `mod_interaction.R`, `mod_crosstissue.R`, `mod_crossancestry.R`) — robustness/generalizability submodules: does the signature hold across tissue, ancestry, or in the presence of covariate interactions?
9. **Enrichment** (`mod_enrichment.R`) — `clusterProfiler::enrichGO`/`enrichKEGG` ORA on the DGE/WGCNA candidate set.
10. **Deconvolution** (`mod_deconvolution.R`) — CIBERSORT/LM22 (via IOBR) as primary estimator, MCP-counter as an independent corroborating check — directly addresses the confound "is this a cell-composition shift, not a per-cell expression change?"
11. **Nomogram / Biomarker Card** (`mod_nomogram.R`, `mod_biomarkercard.R`) — synthesis/reporting layer, not new statistics; the module's terminal output artifact.

**Module-level integration:** how DGE, WGCNA, MR, coloc, and deconvolution jointly triangulate a candidate gene (association + network centrality + causal support + cell-composition-independence) before it's called a "candidate," not any single test alone.
**Findings / Biological interpretation:** `[INSERT RESULT]` throughout — no DGE/MR/coloc output has been reviewed in this session.
**Contribution to subsequent modules:** candidate gene list → Module 3 (gene–CpG concordance) and Module 4 (DIABLO/SNF feature space).

### 5. Module 2 — Methylomics
**Objective:** Parallel epigenomic discovery track, explicitly modeling cell-type composition and array-processing artifacts that transcriptomics does not face.
**Scientific questions:** Which CpGs/regions differ by response/group, correcting for the technical and biological confounds specific to methylation arrays? Does DMP signal organize into co-methylation modules paralleling Module 1's WGCNA? Is any DMP causally upstream of a candidate gene (MR/coloc on methylation QTLs)?
**Data/input:** raw IDAT (unlocks the widest normalization choice) or pre-processed beta/M-value matrix + metadata.
**Preprocessing:** `mod_methyl_normalization.R` — normalization method availability is itself data-adaptive: Noob/Functional normalization/SWAN/Dasen/Stratified-quantile require an `RGChannelSet` (raw IDAT only); BMIQ/PBC/Noob+BMIQ need a Type I/II manifest; plain quantile/no-normalization always available. This adaptivity is a real methodological point — state explicitly that the *tool itself enforces* which normalization is statistically valid for the input actually supplied, rather than letting a user apply an IDAT-only method to an already-processed matrix.

Submodules:
1. **QC** (`mod_methyl_qc.R`, `qc.R`, `idat_metrics.R`) — array-specific QC (detection p-values, bisulfite conversion, sex-check) preceding any statistical test.
2. **Normalization** (above).
3. **Cell-type deconvolution** (`mod_methyl_celltype.R`, `celltype.R`) — EpiDISH (Houseman/Reinius reference, blood 7-cell-type), explicitly noted in-code as the only *real* backend ("Backend honesty: only EpiDISH ... is used"), with unavailable reference methods surfaced rather than silently omitted. This submodule is the methylation-side analogue of Module 1's CIBERSORT deconvolution and its output is a required covariate for DMP, not optional.
4. **DMP** (`mod_methyl_dmp.R`) — limma on M-values with **bacon correction** (genomic-inflation correction for EWAS), row-chunked `lmFit` for genome-wide arrays; the project's own preloaded default: sex-stratified `~ group + age + smoking + cell-type estimates`, BH-FDR. Live/configurable run available for user data; explicitly note the live engine does **not** apply bacon correction — a documented methodological difference between the precomputed and live tracks worth stating plainly (not softening).
5. **DMR** (`mod_methyl_dmr.R`) — region-level aggregation beyond single-CpG DMP calls.
6. **Methylation WGCNA** (`mod_methyl_wgcna.R`) — co-methylation network, module 2's parallel to Module 1's network submodule.
7. **Candidates / Feature selection / MR / Coloc / Diagnostic / Biomarker Card** (`mod_methyl_candidates.R`, `mod_methyl_featureselection.R`, `mod_methyl_mr.R`, `mod_methyl_coloc.R`, `mod_methyl_diagnostic.R`, `mod_methyl_biomarkercard.R`) — mirrors Module 1's structure; **flag in the chapter which of these are fully built vs. registry-scaffolded placeholders** (per `submodules_registry.R`'s own comments, Cell-type through Diagnostic were queued for incremental build-out) — a PhD committee will read a claimed-but-scaffolded submodule as an integrity problem if not disclosed.

**Module-level integration / findings / interpretation:** `[INSERT RESULT]`.
**Contribution to subsequent modules:** DMP/DMR list + cell-type-adjusted matrix → Module 3 concordance, Module 4 integration.

### 6. Module 3 — Cross-Omics
**Objective:** Test direct regulatory concordance between transcriptomic and methylomic signal at the gene level — the pairwise linkage step between single-omics discovery and joint multi-omics integration.
**Scientific questions:** Do hyper/hypomethylation and up/down-regulation co-occur at the same genes more than chance, and in the biologically expected direction (hypermethylation-associated repression, hypomethylation-associated activation)?
**Data/input:** Module 1's DEG output + Module 2's DMP output, joined at the gene level.
**Preprocessing:** gene-level aggregation of CpG-level `dbeta` signs (`crossomics_integration_helpers.R`) — hyper defined as `dbeta > 0`, hypo as `dbeta < 0`, tabulated per gene against DEG direction.

Submodules:
1. **Integration** (`mod_cross_integration.R`) — the four-way concordance categories actually implemented: **Hyper+Down** ("potential methylation-associated repression"), **Hypo+Up** ("potential methylation-associated activation") — the two *canonical* expected directions — versus **Hyper+Up**/**Hypo+Down** ("concordant-direction/potentially noncanonical") — explicitly labeled as association, not mechanism (`build_assistant_context()` appends a standing caveat: "these are statistical associations, not established causal relationships — never state that methylation 'causes' an expression change" — this is a real, code-enforced epistemic guardrail worth quoting directly in the chapter as evidence of methodological discipline).
2. **Biomarker Convergence** (`mod_cross_biomarker_conv.R`) — registry-scaffolded per `submodules_registry.R`; disclose build status.
3. **MR Stage** (`mod_cross_mr_stage.R`) — registry-scaffolded; disclose build status.

**Module-level integration:** counts/proportions in each of the 4 concordance categories, tested against a null (chance co-occurrence) — **[INSERT PARAMETER]** for the actual test used (hypergeometric/Fisher's exact, most likely, but verify in code before stating).
**Findings:** `[INSERT RESULT]` — n genes analyzed, n significant DEGs/DMGs/both, per-category counts (the live app's own AI-context builder already computes exactly these fields at runtime — `sprintf("Significant DEGs: %s, Significant DMGs: %s, Significant in both: %s", ...)` — so this is a case where pulling the real number from a live app run, once available, is trivial and should replace the placeholder rather than staying abstract).
**Contribution to subsequent modules:** concordant gene set is one of the candidate feature pools offered to Module 4's DIABLO/SNF panels.

### 7. Module 4 — Multi-Omics Integration
**Objective:** Joint (not pairwise) integration of transcriptomic and methylomic layers for supervised biomarker discovery, unsupervised patient stratification, and pathway-level convergence — the module the template calls out for extra depth.
**Scientific questions:** Can a supervised multi-block model discriminate response/outcome using jointly-selected features across omics layers, and does its performance survive leakage-safe validation? Do unsupervised, omics-fused patient similarity structures reveal clinically meaningful subgroups independent of any outcome label? Do these two independent lines of evidence (supervised DIABLO panel, unsupervised SNF clusters) converge on the same genes/CpGs?
**Data/input:** `multi_dataset` shared reactiveValues — sample-matched, harmonized transcriptomic + methylomic layers (`cohort_harmonization_helpers.R` enforces the matching before either method runs).
**Preprocessing:** sample harmonization/cohort matching; DIABLO and SNF maintain **entirely separate parameter panels**, a deliberate design choice worth stating explicitly (their optimal preprocessing genuinely differs — supervised sparse selection vs. unsupervised affinity-network construction — so the tool does not force one shared preprocessing path onto both).

Submodules — this is the section to give the template's requested "additional detail":
1. **Overview** (`mod_multi_overview.R`) — cross-layer QC scorecard (`multi_qc_scorecard()`), cohort characterization at the point of entry to joint integration.
2. **Integration: DIABLO + SNF** (`mod_multi_integration.R`) —
   - **DIABLO** (`mixOmics::block.splsda`, supervised): a multi-block sparse discriminant model, `X_m` blocks per omics layer jointly regularized against a categorical outcome; requires a genuinely categorical (not continuous) outcome column, enforced in-app before the run is allowed. State the design (block correlation weighting, `keepX` per block, component count) as **[INSERT PARAMETER]** pending a code read of the live-run defaults.
   - **SNF** (`SNFtool::SNF`, unsupervised): per-omics-layer similarity/affinity network construction, fused via the SNF algorithm (Wang et al. 2014) into one joint patient-similarity network; **K/Alpha/T** are the three tunable parameters — neighbor count, fusion diffusion iterations, hyperparameter — surfaced directly in the UI (`mod_multi_integration.R`), so document their live-run values as **[INSERT PARAMETER]** rather than default library values.
   - Single-omics comparison mode and sex-stratified analysis are both explicit, code-level options — cite this as the mechanism by which the module tests whether integration *adds* signal over either omics layer alone (directly answerable via §9's already-existing benchmark, see below).
3. **Patient Stratification / SNF Clustering** (`mod_multi_stratification.R`) — downstream of SNF's fused network: spectral clustering into patient subgroups, with (per in-code workflow comment) an explicit Data→Validate→Configure→Run→Clusters→**Stability**→Clinical→Features pipeline — clusters are unsupervised molecular groups first, clinical/outcome association tested only afterward, never used to *choose* the cluster count (a leakage-avoidance discipline worth naming explicitly — this is the same discipline the AUDIT.md leakage fixes exist to enforce, and stating it here strengthens the cross-reference to §9).
4. **Biomarker Discovery** (`mod_multi_biomarker.R`) — the panel of features (genes + CpGs) DIABLO/SNF jointly implicate; this is the artifact §8's cross-module prioritization consumes.
5. **Gene–CpG Concordance** (`mod_multi_concordance.R`) — a joint-integration-aware version of Module 3's pairwise concordance: cross-references DIABLO's and SNF's implicated features against the CpG↔gene annotation map, restricted to features actually present in the dataset (explicit code-level refusal to fabricate a gene/CpG/coordinate/statistic when absent — another quotable methodological-discipline point).
6. **Pathway** (`mod_multi_pathway.R`) — live GO/KEGG/Reactome/WikiPathways/MSigDB Hallmark ORA and GSEA on the joint candidate panel, with per-database availability gated on installed annotation packages (`ReactomePA`/`reactome.db`/`fgsea`/`msigdbr`) — document actual availability in the deployed environment rather than assuming all four ran.
7. **(Embedded in Dataset Workspace) MOFA2** (`mod_multi_live_mofa.R`) — unsupervised latent-factor model, `X_m ≈ Z W_m^T + E_m` (Argelaguet et al. 2018, 2020); explicitly documented in-app as capturing *variation*, not disease status or causality — a third, fully independent unsupervised method (alongside SNF) whose factor structure can be cross-checked against SNF's clusters as a convergence/robustness check.

**Module-level integration:** the module's own internal convergence check — do DIABLO's supervised panel, SNF's cluster-associated features, and MOFA2's top-loading features overlap more than chance on the same gene/CpG set?
**Findings:** here is where §9's real, audited numbers belong, not a placeholder — quote directly from `AUDIT.md`:
> Leakage-safe, sex-stratified nested-CV AUROCs, Models A (transcriptomics-only) / B (methylomics-only) / C (integrated): female (n=56) 0.338 / 0.427 / 0.430; male (n=23) 0.439 / 0.523 / 0.409 — none exclude chance in either sex, and the integrated model does not outperform either single-omics model (DeLong female A-vs-C p=0.41, B-vs-C p=0.92).
> Drug-type DIABLO (sex-stratified, leakage-corrected): female AUROC 0.579 [0.428, 0.729], male 0.504 [0.263, 0.745].
This is a **scientifically important negative/null result**, and the chapter should present it exactly as AUDIT.md does — plainly, with the caveat that these numbers come from the underlying research pipeline's own nested-CV benchmark (not the live Shiny DIABLO/SNF run itself, which is a separate, data-adaptive computational path — state that distinction explicitly so the two are never conflated).
**Biological interpretation:** `[INSERT RESULT]` for anything beyond the AUROC/leakage story above (e.g., which genes/CpGs the DIABLO loadings/SNF clusters actually implicate).
**Contribution to subsequent modules:** feeds §8 (cross-module prioritization) directly; is itself the module §9's validation section is built around.

### 8. Cross-Module Integration
**Purpose:** Make explicit, with a diagram, the actual convergent-evidence logic already implied by §3's data-flow description — this section should read as the payoff of the loop, not a summary.
- Candidate genes/features: Module 1 (DGE+WGCNA+MR/coloc-supported) ∩ Module 2 (DMP+cell-type-adjusted) ∩ Module 3 (concordant direction) ∩ Module 4 (DIABLO/SNF/MOFA2-implicated).
- Evidence convergence table: rows = candidate genes, columns = which of the 4 modules' independent lines of evidence support each one — this is the figure/table that operationalizes "multi-omics convergence" as a countable thing rather than a claim.
- Patient stratification ↔ pathway convergence: do SNF clusters differ in pathway enrichment (Module 4 Pathway submodule) in a way that's consistent with the DIABLO panel's biology?
- Clinical associations: tested only after unsupervised structure is fixed (per §7's stated discipline) — never used to tune cluster count.
- Final candidate prioritization: explicit ranking rule (e.g., number of independent modules supporting a feature, or a weighted score) — **[INSERT PARAMETER]** for whatever rule is actually implemented/chosen; do not invent a scoring formula not present in the code.
- Integrated biological model: a synthesis figure/diagram (see Figures list) — the "final biological/conceptual model" the template asks for, built only from confirmed findings once §4–§7's `[INSERT RESULT]` placeholders are filled.

### 9. Validation and Benchmarking
**This section is unusually strong already** because real audit work exists — build it around `AUDIT.md`, cited directly, not around a hypothetical validation plan.
- Discovery vs. validation datasets: state directly what exists (Tao et al. 2021 discovery cohort; note whether any independent validation cohort was used — `[INSERT DATA]` if not, since inventing one would violate the no-fabrication rule).
- Cross-dataset / cross-platform validation: Module 1's Cross-Tissue/Cross-Ancestry submodules are the actual generalizability tests available — describe what they test, `[INSERT RESULT]` for outcomes.
- **Statistical validation / nested CV — the chapter's strongest subsection.** Structure it as a case study:
  1. What the pipeline's own methodology claims (leakage-safe, train-only feature selection, per-fold imputation).
  2. What the independent audit found: two **critical** data-leakage bugs (global instead of per-fold methylation imputation) in two independently-written scripts (script 14 and script 06), plus one **moderate**, documented-but-unfixed leak (`filterByExpr` using full-cohort outcome labels for gene-universe selection, script 02).
  3. What changed: pre-fix vs. post-fix AUROCs, quoted exactly as in AUDIT.md (§7 above has the numbers) — frame this explicitly as *the leakage correction moved the numbers modestly and did not reverse the qualitative conclusion*, which is itself evidence the surrounding pipeline discipline (nested CV, per-fold imputation elsewhere) was already largely sound.
  4. Benchmark vs. published: this project's own integrated-model accuracy (female 0.411, male 0.565) against Tao et al.'s reported 0.79–0.89 — state the audit's own honest reading: a sex-stratified, genome-wide, in-fold-selected approach at this sample size does not replicate the published accuracy, and say why that's a legitimate, disclosable finding rather than a failure to hide (differently-validated pipeline, smaller per-sex n, stricter leakage control).
- Feature stability / sensitivity / parameter robustness: **[INSERT RESULT]** — check whether AUDIT.md's later sections (script 07b SNF, script 13 Random Forest) or other project docs contain stability analyses before writing this as pure placeholder; worth a follow-up code/doc read.
- Overall audit score for transparency: quote the AUDIT.md summary table (methodological rigor 22/25, scientific honesty 20/20, code structure 16/20, documentation 10/15, reproducibility 10/12, robustness/guardrails 0/8 — total 78/100) as a **validation artifact in its own right**: few multi-omics platforms publish an external audit of their own analysis code, and doing so is itself a reproducibility contribution worth naming in §12's Discussion.
- Comparison with existing methods/platforms: cross-reference §10, don't duplicate.
- Usability evaluation: `[INSERT DATA]` — no user study exists unless you tell me otherwise.

### 10. Comparison With Existing Platforms
- Systematic table (see Tables list) across: data handling, QC, statistical analysis, multi-omics integration, biomarker discovery, patient stratification, pathway analysis, visualization, reproducibility, validation, extensibility.
- Three explicit columns per capability: **exists elsewhere** / **ArthOMix has it** / **novel contribution** — the brief's own required distinction. Novel-contribution candidates to argue for (pending literature check against the five platforms, `[INSERT CITATION]`):
  - Dual precomputed-audited + live-adaptive computational tracks on one shared data contract.
  - Both supervised (DIABLO) *and* unsupervised (SNF, MOFA2) integration with cross-method convergence checking, rather than one fixed integration method.
  - A published, independently-scored leakage audit with before/after numbers.
  - Sex-stratified analysis as a first-class option across all four modules, not a post-hoc subgroup slice.
  - Explicit epistemic guardrails enforced in code (association-not-causation caveats, refusal to fabricate absent features) rather than left to documentation alone.

### 11. Biological and Scientific Interpretation
All subsections `[INSERT RESULT]`-gated except the two already-established negative findings from §9 (integration not outperforming single-omics; project's own accuracy not replicating Tao et al.'s published number) — lead with those as the chapter's most defensible biological/statistical claims, then layer in `[INSERT RESULT]` for gene/pathway-level specifics once Modules 1–4 are run and reviewed.

### 12. Discussion
- Principal findings: anchor on the audited, real findings from §7/§9, not on invented biology.
- Scientific contribution vs. Methodological contribution — keep these genuinely separate subsections; the methodological one (dual-track architecture, leakage audit, sex-stratification) is the stronger, better-evidenced of the two right now.
- Comparison with previous research: Tao et al. 2021 benchmark comparison belongs here too (cross-ref §9).
- Strengths: audited leakage discipline; modular reusable architecture; dual precomputed/live design; explicit epistemic guardrails.
- Limitations: honestly state which submodules are registry-scaffolded, not built (Module 2's Cell-type→Diagnostic tail, Module 3's Biomarker Convergence/MR Stage); the one unfixed moderate leak (script 02); small per-sex n (23–56) limiting statistical power — this is exactly the kind of self-aware disclosure AUDIT.md itself models, and a thesis committee will read consistent self-disclosure across the whole chapter as a strength, not a weakness.
- Biases/confounding: sex-stratification itself, cell-type composition, batch effects.
- Generalizability: single cohort, `[INSERT DATA]` on any second cohort.
- Future directions: complete the scaffolded submodules; extend the leakage audit to the live/data-adaptive track (currently the audit only covers the precomputed pipeline, per your project memory — the live DIABLO/SNF/MOFA2 track has its own separate verification history, see [[project-arthomix-integration-rewrite]], but not the same nested-CV leakage audit AUDIT.md performed).

### 13. Conclusion
Concise synthesis — no new claims, only restated from §4–§12, each with its section back-reference.

---

## 3. Recommended Figures

| # | Figure | Scientific information it must communicate |
|---|---|---|
| 1 | Overall framework/workflow diagram | The four-module pipeline (§3) as a directed graph with the actual feedback loop (Module 3/4 outputs re-entering earlier modules' candidate lists), not a one-way flowchart |
| 2 | System architecture | `config`/`ui`/`server` registry pattern, per-area module lists, shared `multi_dataset`/`multi_results`/`results`/`cross_results` reactive objects — what actually enables cross-module data flow, not a generic Shiny diagram |
| 3 | Data flow between modules | Concretely: DEG list + DMP list → Module 3 concordance categories → candidate pool → Module 4 DIABLO/SNF/MOFA2 → biomarker panel → back into Module 3's convergence panel |
| 4 | Module 1 workflow | DGE → WGCNA → MR/coloc → deconvolution triangulation, ending at the candidate gene artifact |
| 5 | Module 2 workflow | IDAT/matrix input-adaptive normalization branch → cell-type deconvolution as a required DMP covariate → DMP/DMR → candidates |
| 6 | Module 3 workflow | Gene-level aggregation of CpG-level `dbeta` sign vs. DEG direction into the 4 concordance categories, with the canonical-vs-noncanonical distinction visually explicit |
| 7 | Module 4 multi-omics workflow | Parallel DIABLO (supervised) / SNF (unsupervised) / MOFA2 (unsupervised, embedded) branches converging on one biomarker panel + stratification + pathway layer |
| 8 | Cross-module integration diagram | The evidence-convergence table (§8) rendered as a Venn/UpSet-style diagram over the 4 modules' independent candidate sets |
| 9 | Final biological/conceptual model | Only buildable once §11's `[INSERT RESULT]`s are filled — placeholder until real findings exist |
| 10 | Leakage audit before/after | AUROC forest plot, pre-fix vs. post-fix, both models (A/B/C) both sexes — this is a genuinely novel figure type for a multi-omics platform chapter and worth its own callout |

## 4. Recommended Tables

| # | Table | Content |
|---|---|---|
| 1 | Dataset characteristics | Tao et al. 2021 cohort: n per sex/drug/response category (from AUDIT.md: good=22/moderate=17/no=41 pooled — confirm sex-stratified breakdown, e.g. n=56 female / n=23 male, against the actual DIABLO/nested-CV tables before citing) |
| 2 | Input/output requirements per module | What each of the 4 modules needs as input and emits as output (the §3 data-flow table, formalized) |
| 3 | Analytical methods by module | limma/DESeq2, WGCNA, TwoSampleMR/mr_ivw, coloc.abf, CIBERSORT/MCP-counter (Module 1); Noob/Funnorm/SWAN/Dasen/BMIQ/PBC, EpiDISH, limma+bacon (Module 2); concordance categorization test (Module 3); DIABLO/SNF/MOFA2, GO/KEGG/Reactome/Hallmark ORA+GSEA (Module 4) |
| 4 | Parameters by method | `[INSERT PARAMETER]` for every live-run default — needs a dedicated code read of each module's config before filling in |
| 5 | Submodule objectives and outputs | One row per submodule across all 4 modules — essentially §4–§7's tables consolidated |
| 6 | Comparison with existing platforms | §10's exists-elsewhere / ArthOMix-has-it / novel-contribution table |
| 7 | Validation strategy | Nested CV design (5-fold × 5/10-repeat), leakage bugs found/fixed/documented, pre/post-fix AUROCs — directly from AUDIT.md |
| 8 | Final prioritized biomarkers/features | `[INSERT RESULT]` — populated only after §8's convergence analysis is actually run |

## 5. Logical connections between modules (for §3's explicit treatment)

Single-omics discovery (1, 2) run independently but on harmonized, sample-matched data → their outputs are paired at the gene level for direct regulatory concordance (3) → the concordant/candidate pool becomes one of several feature inputs to joint multi-omics integration (4), which itself runs three independent methods (supervised DIABLO, unsupervised SNF, unsupervised MOFA2) whose agreement is itself evidence → Module 4's outputs feed back into Module 3's convergence panel and forward into cross-module prioritization (8) → every quantitative claim anywhere in this loop is subject to the leakage audit discipline demonstrated concretely in Module 4 and formalized in §9.

## 6. Methodological contributions (draft list, to sharpen in §12)

1. A shared, typed reactive-results contract (`multi_results`, `cross_results`, per-area `results`) that makes cross-module consumption of prior outputs a first-class architectural feature rather than manual re-wiring.
2. Dual-track design: precomputed/audited pipeline results alongside a live, data-adaptive computational path on the same UI, explicitly distinguished rather than blended.
3. A published, independently-scored leakage audit (AUDIT.md) with concrete before/after performance numbers — the concrete instantiation of "reproducibility" the field's other tools mostly gesture at.
4. Input-adaptive method gating (methylation normalization choices constrained by whether raw IDAT or a processed matrix was supplied; DIABLO's categorical-outcome requirement enforced before the run is allowed) as a mechanism for preventing statistically invalid analyses, not just documenting the constraint.
5. Code-enforced epistemic guardrails (association-not-causation caveats auto-appended to cross-omics output; refusal to fabricate absent gene/CpG identifiers) as a design pattern, not a documentation promise.

## 7. Expected scientific contributions (draft list, `[INSERT RESULT]`-gated where real findings are needed)

1. A concrete, sex-stratified negative/null finding on integrated-model predictive value for anti-TNF response in this cohort — scientifically useful precisely because it's disclosed rather than suppressed.
2. `[INSERT RESULT]`: specific candidate genes/CpGs surviving 4-module convergence (§8), once §4–§7 are populated with real outputs.
3. `[INSERT RESULT]`: patient stratification structure (SNF/MOFA2) and its relationship (if any) to clinical/drug-response covariates, established only after the unsupervised structure is fixed.
4. A methodological/tooling contribution independent of any specific biological finding: the architecture and audit discipline itself, defensible even if every `[INSERT RESULT]` above ultimately reads "no significant signal" — worth stating explicitly in §12 so the chapter's contribution doesn't rise or fall on any single biological result.

---

## 8. Open items before full-prose drafting

- **[INSERT PARAMETER]** items above need a code read of each module's live-run defaults (DGE design formula/threshold, DIABLO `keepX`/component count, SNF K/Alpha/T, concordance significance test) before §4/§5/§7 can state them precisely.
- **[INSERT RESULT]** items need either (a) an actual run of the live app on the bundled dataset, or (b) a review of `Research_05_multiomics_sexstratified/analyses/*/results/tables/` for already-computed values beyond what AUDIT.md quotes — likely a rich additional source I haven't read yet.
- **[INSERT CITATION]** items (the five comparator platforms' primary papers, RA anti-TNF non-response rate, DIABLO/SNF/MOFA2/limma/DESeq2/EpiDISH/CIBERSORT/MCP-counter original method papers) need a literature pass — the `claude.ai PubMed`/`bioRxiv` connectors in this environment aren't authorized yet, so that pass needs to happen either via your own reference manager or once those connectors are enabled.
- Confirm the sex-stratified cohort n's (56 female / 23 male, from AUDIT.md) against Table 1's intended source before it goes in the chapter as a hard number.

---

**Next step, when you're ready:** tell me which module (or the cross-module/validation sections) to draft into full prose first, and whether to pull additional real result values from `Research_05_multiomics_sexstratified/analyses/*/results/` before writing rather than leaving them as `[INSERT RESULT]`.
