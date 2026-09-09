# Supplementary Material for: [Title] — fill-in template

**How this fits together.** Bioinformatics Advances submissions are two files: the main text (`ArthOMix_BioinfAdv_ApplicationNote_template.md`, ≤4 pages, ≤3 figures/tables, ≤15 references) and a separate Supplementary Material PDF with no length cap. Everything trimmed from the main text for space — extra case-study plots, the full feature-comparison table, extended parameter detail — lives here instead. The main text references this file's contents as "Fig. S1", "Table S1", etc.; every one of those references in the main text must resolve to something that actually exists below, and vice versa — don't leave an orphaned Fig. S4 with nothing citing it, and don't cite Fig. S4 from the main text without it existing here.

`[AUTHOR]` confirm Bioinformatics Advances' exact supplementary-file format requirement (single combined PDF vs. separate files per item) from the submission portal — this template assumes a single combined PDF, which is the more common convention, but verify.

---

## S1. Supplementary Methods

Extended versions of anything compressed in main-text Methods §2.1-2.6 to hit the word budget. Concretely, once you've trimmed the main text, whatever got cut goes here in full. Candidates already flagged as compressible in the main template:

- **Full parameter tables** — every user-adjustable control per sub-module (not just the ones named in main-text Methods), with its default value and range. One table per pipeline is cleanest:

  `[AUTHOR — build from the actual UI, one row per control]`

  | Pipeline | Sub-module | Control | Default | Range/options |
  |---|---|---|---|---|
  | Transcriptomics | Differential Expression | BH-FDR cutoff | 0.05 | 0-1 |
  | Transcriptomics | Differential Expression | \|log2FC\| cutoff | 0.5 | any positive |
  | Transcriptomics | WGCNA | R² cutoff (soft-power selection) | `[AUTHOR: check UI default]` | — |
  | Transcriptomics | Feature Selection | Elastic-net alpha | 1 (LASSO) | 0-1 |
  | Transcriptomics | Feature Selection | CV folds | 10 | — |
  | Methylomics | QC | Detection p-value threshold | 0.01 | — |
  | Methylomics | DMP | Covariates | group, age, smoking, cell-type | — |
  | Methylomics | DMR (DMRcate) | lambda / C | 1000 / 2 | — |
  | Multi-Omics | DIABLO | ncomp mode | tuned (perf, ≥3 repeats) | manual override available |
  | Multi-Omics | DIABLO | keepX | manual (e.g. "20,10") | or auto-tuned grid search |
  | Multi-Omics | DIABLO (**Case study 4, specifically**) | ncomp / keepX | **manual: ncomp=2, keepX="20,10" per block** | Not the tuned/auto-tune default — a deliberate cost/timeline tradeoff (the tuned default's grid search ran >15 min with no completion in sight and was stopped). Disclosed here per the Methods §2.5/Results §7 cross-reference. |
  | Multi-Omics | Biomarker Discovery | Hold-out fraction | 0.20 | 0-0.5 |
  | `[AUTHOR]` | … | … | … | … |

- **Full array-platform / normalization detail** for Methylomics if trimmed from main text (§2.3 already flags Touleimat & Tost 2012 as the first thing to cut there).
- **Full ArthoChat tool list** — the 6 typed tools main-text §2.6 counts but doesn't name individually: `[AUTHOR: list each tool's name and one-line purpose, e.g. "propose_run_dge — proposes a differential expression run with explicit parameters, requires user confirmation before executing"]`.
- **Provenance manifest schema** — what fields get logged per analysis run, if you want to back up the "reproducible run" claims made in Methods.

---

## S2. Table S1 — Feature comparison

`[VERIFIED 2026-09-08 — every "No" below is sourced from actually fetching and reading the OmicsAnalyst and xOmicsShiny papers this session, not guessed. This table must stay in sync with Introduction §2's competitive-gap paragraph in the main template — they were drafted from the same verification pass; if you edit one, check the other.]`

| Feature | ArthOMix | OmicsAnalyst (Zhou et al. 2021) | xOmicsShiny (Gao et al. 2025) |
|---|---|---|---|
| Raw-data differential analysis (DE/DMP), not precomputed input only | Yes (limma, sex-stratified) | Yes (ANOVA/t-test) | **No** (accepts precomputed DE results only) |
| Sex-stratified / sex-pooled analysis at every stage | **Yes, unique among the three** | No | No |
| WGCNA co-expression/co-methylation | Yes | No | Yes |
| Multi-omics integration: DIABLO | Yes | Yes | No |
| Multi-omics integration: Similarity Network Fusion | Yes | Yes | No |
| Multi-omics integration: MOFA2 | Yes | No | No |
| Cross-omics Mendelian Randomization convergence | **Yes, unique among the three** | No | No |
| Sealed hold-out biomarker evaluation | **Yes, unique among the three** | No | No |
| External independent-cohort validation | **Yes, unique among the three** | No | No |
| Pathway enrichment (ORA/GSEA) | Yes | Yes (Reactome shown in case study) | Yes (WikiPathways/Reactome/KEGG) |
| GEO direct accession fetch | Yes | No (manual upload) | No (manual upload) |
| Docker deployment | Yes | No | No |
| Integrated conversational AI assistant grounded in live results | **Yes, unique among the three** | No | No |

`[AUTHOR]` this table is the entire evidentiary backing for the "uniquely supports" language in Introduction §2 and the Abstract — if a reviewer challenges that language, this table (and the two source papers) is your answer. Do not add a row claiming a fourth capability without re-verifying it the same way (fetch the actual comparator paper, don't infer from its abstract). Add a 4th comparator column (e.g. the reference paper this template is modeled on, Dyce et al. 2026 Omics BioAnalytics) if you want a stronger table — its case studies (SARS-CoV-2, heart failure) suggest it also lacks sex-stratification and cross-omics MR, but `[AUTHOR: verify by reading that paper directly before asserting it in a table, same rule as above]`.

---

## S3. Supplementary Figures

One entry per figure not selected for the main text's single representative plot (main template §7/§14). Candidates, one per case study not chosen as the main-text figure:

`[DONE 2026-09-08 — all 4 generated directly from the real saved run data (not mocked), at repo root `figures/`. Figure 2 in the main text is FigS4b (the training-vs-hold-out AUROC comparison) — everything below is the rest.]`

- **Fig. S1.** `figures/FigS1_case_study1_volcano.png` — Case study 1 (Transcriptomics) volcano plot, sex-adjusted DE on GSE55457 (821 DEGs highlighted, FDR<0.05/|log2FC|>0.5).
- **Fig. S2.** `figures/FigS2_case_study2_volcano.png` — Case study 2 (Methylomics) volcano plot, SVA-adjusted DMP on GSE71841 (268 DMPs highlighted, bacon-FDR<0.05; non-significant points subsampled to 40k for file size, all significant points shown).
- **Fig. S3.** `figures/Fig3_case_study3_quadrant.png` — Case study 3 (Cross-Omics) quadrant plot, expression log2FC (GSE55457) vs. methylation Δβ (GSE71841), HLA-DQA1/HLA-DQB1 labeled as the convergent genes.
- **Fig. S4.** `figures/FigS4a_case_study4_holdout_roc.png` — Case study 4 (Multi-Omics) sealed hold-out ROC curve (n=6, AUROC=0.67).

`[AUTHOR]` all 4 need final captions written in full supplementary style (dataset, n, method, what the reader is looking at, per the S3 instructions above) before submission — the labels above are working titles, not final captions.
- **Fig. S5+.** `[AUTHOR: architecture/deployment diagram, if you want one beyond main Figure 1 — e.g. showing the Docker/renv dependency stack]`

Each entry needs, at minimum: a figure number, the plot itself (publication resolution), and a full caption (unlike main-text figures, supplementary captions can be as long as needed — state the dataset, n, method, and what the reader is looking at without needing the main text open).

`[AUTHOR]` once Track B1-B4 (the four live runs) are done, you'll have more raw plots than you need — this section is where the "leftover" ones go rather than being discarded. Number them in the order the case studies appear in main-text §7, not the order you happened to generate them.

---

## S4. Supplementary Tables (beyond Table S1)

`[DONE 2026-09-08]` **Table S2 — Transcriptomics pipeline: methodology decision guide.** Full content in `ArthOMix_Transcriptomics_DecisionGuide_Supplementary.md` (too long to inline here — one table per submodule, 17 submodules, every user-facing control that changes which statistical method/model runs). Built by reading the entire `R/transcriptomics/` codebase and extracting, for every such control: what the options mean, what the app actually does differently for each, and when to pick which and why — quoting the app's own validation/tooltip text rather than inventing rationale, and stating plainly wherever no such guidance exists in the code. Motivated directly by a live incident this session: uploading GSE77298 (a microarray dataset) and running it through DESeq2 (meant for raw sequencing counts) produced an opaque `"some values in assay are negative"` error with no indication of the real cause (a Data Type / Method mismatch). Table S2's Differential Expression section documents exactly this decision point. Its closing section lists 7 documentation/UI issues surfaced as a byproduct of the audit: 6 were real and have since been **fixed in code** (each independently re-verified, plus the full existing test suite re-run with zero failures) — most importantly the exact root cause of the GSE77298 incident (`looks_like_raw_counts()` now requires near-integer values, not just a wide non-negative range, so linear-scale microarray intensities can no longer be misclassified as raw sequencing counts) and a stale statistical-criterion description in Candidate Gene Identification that no longer matched the code. One initial finding (a Cross-Tissue AUC-orientation claim) was checked algebraically and **retracted** as not actually a bug rather than "fixed" incorrectly — worth citing as evidence of the verification standard applied throughout.

`[AUTHOR — optional, add only if you have material that doesn't fit a figure]`:
- A per-case-study full results table (e.g. every significant DEG/DMP, not just the top hits shown in main-text Table 1) — common practice, keeps main text lean while giving reviewers full reproducibility.
- The full parameter-default table from S1 above, if you'd rather present it as a table than inline in Methods.

---

## S5. Supplementary References

`[AUTHOR — check requirement]` Some OUP journals want supplementary-only citations numbered separately from the main reference list; others fold everything into one shared list. If you cite something in Supplementary Methods (S1) or a figure caption that isn't already in the main text's References (§13), it needs to appear either here (if separate numbering is required) or added to §13 (if shared). Confirm which convention applies from the submission portal before finalizing either list.

---

## S6. Video Walkthroughs

`[ADDED 2026-09-08]` Four short videos, one per pipeline, mirroring the reference paper's own Availability precedent ("video walkthroughs" alongside code and demo). **Each video's script is already written** — it's the "How to reproduce/verify in the browser" steps already in `ArthOMix_CaseStudy_RunLog_template.md` for that case study. This isn't a coincidence: walking through the exact steps that reproduce the paper's own Results is more valuable to a reader than a generic feature tour, and it doubles as your own final verification pass while you record.

| # | Pipeline | Script source | Target length |
|---|---|---|---|
| 1 | Transcriptomics | Run Log, Case study 1 "How to reproduce/verify in the browser" — upload GSE55457 files, run DE with defaults, re-run with sex as covariate, compare 2319 vs 821 DEGs | 3-5 min |
| 2 | Methylomics | Run Log, Case study 2 — upload GSE71841 files, SVA-adjusted Analysis panel, compare 268 DMPs, top hit cg14075454 | 3-5 min |
| 3 | Cross-Omics | Run Log, Case study 3 v2 — upload the two independent DEG/DMP tables (not the paired v1 data), Integration tab, min_fdr aggregation, HLA-DQA1/HLA-DQB1 result | 4-6 min (the aggregation-method choice is worth explaining on camera, it's a real methods decision not a default click-through) |
| 4 | Multi-Omics | Run Log, Case study 4 — GSE117931 upload, Data Workspace same-patient gate, DIABLO Integration (mention both the tuned-default and the manual-parameter path used in the paper, and why), Biomarker Discovery hold-out, compare 0.90/0.67 AUROC | 5-7 min (the longest one — has the most to honestly explain) |

`[AUTHOR]`:
- Record against the exact frozen commit (`paper-case-studies-freeze-20260908`, `2627687`) so what's on screen matches what the paper describes — don't record after resuming normal development post-freeze without re-checking the numbers still match.
- Host externally (YouTube unlisted, or alongside the code deposit on Zenodo/OSF if you want a citable, permanent link rather than a platform-dependent one) — video files don't belong in the Supplementary PDF itself.
- Add the resulting link(s) to the main manuscript's Availability paragraph (§4) — already has a placeholder waiting for this.
- Captions/narration should say what's being clicked *and* why (e.g. "we choose sex as covariate here because the app's own confound check flagged it" for video 1) — a silent screen recording is a feature tour, not a walkthrough that helps someone trust the paper's numbers.

---

*Companion to `ArthOMix_BioinfAdv_ApplicationNote_template.md`. Built 2026-09-08. Table S1 and Introduction §2 of the main template were verified together in the same session — keep them in sync if either changes.*
