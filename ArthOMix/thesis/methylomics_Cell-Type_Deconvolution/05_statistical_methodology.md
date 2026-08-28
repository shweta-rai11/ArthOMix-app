# 05. Statistical Methodology — Methylomics Cell-Type Deconvolution

## 5.1 The reference-based deconvolution model (general theory)

Every method this module offers shares the same underlying linear mixture model. For a bulk sample with observed methylation vector `β_obs` over a set of `m` marker CpGs, a reference matrix `R` (`m` CpGs × `K` cell types, each column the mean methylation of one purified/sorted cell type at those CpGs), and an unknown proportion vector `w` (length `K`):

```
β_obs ≈ R %*% w,     subject to  w_k ≥ 0  for all k   (and often  Σ w_k = 1)
```

The methods differ in how they fit `w`:

- **Houseman constrained projection (CP).** The original 2012 method: minimize `‖β_obs − Rw‖²` subject to `w ≥ 0` (and optionally `Σw = 1`), solved as a quadratic program. This is the textbook description; see §5.2 for exactly how this codebase invokes it.
- **RPC (robust partial correlations).** Fits `w` via a robust linear regression of `β_obs` on `R` (down-weighting outlier CpGs rather than a plain least-squares fit), then clips negative coefficients to zero and renormalizes to sum to one.
- **CBS (CIBERSORT-style).** Fits `w` via nu-support-vector regression (linear kernel) of `β_obs` on `R`, again clipping and renormalizing. Originally a gene-expression deconvolution idea (Newman et al. 2015's CIBERSORT), adapted to methylation data by `EpiDISH`.
- **hepidish (two-stage/hierarchical).** Runs one of the above three methods twice: once against a top-level reference (e.g. Epithelial/Fibroblast/Immune), and independently against a second, finer-grained reference for one specific top-level component (e.g. blood-cell subtypes within "Immune"), then multiplies the top-level component's fraction into the second-stage sub-fractions to produce one combined, fully-resolved proportion vector.

## 5.2 Exactly how this codebase invokes each method (verified against the installed `EpiDISH` package's own source)

All three single-stage methods are reached through one call, `EpiDISH::epidish(beta.m, ref.m, method, maxit, nu.v, constraint)`, wrapped by `methyl_ct_run_epidish()` (`celltype.R:299-319`). This documentation set independently read the installed `EpiDISH` package's exported and internal source (`EpiDISH::epidish`, `EpiDISH:::DoCP`, `EpiDISH:::DoRPC`, `EpiDISH:::DoCBS`, `EpiDISH::hepidish`) in this environment to confirm the following, rather than taking the code's own comments on faith:

**CP (`method = "CP"`)** — `EpiDISH:::DoCP()` builds a quadratic program and calls `quadprog::solve.QP()` per sample:
- `constraint = "inequality"` (this module's UI default, `mod_methyl_celltype.R:193-195`): the QP constraint matrix enforces `w_k ≥ 0` for every cell type and `Σw_k ≤ 1` (an *inequality* on the sum — not a strict equality), passed as `meq = 0` (zero leading equality constraints) to `solve.QP()`.
- `constraint = "equality"`: enforces `w_k ≥ 0` and `Σw_k = 1` exactly, passed as `meq = 1`.
- No random-number generation occurs anywhere in `DoCP()` — it is a deterministic numerical optimization for a fixed input.

**RPC (`method = "RPC"`)** — `EpiDISH:::DoRPC()` calls `MASS::rlm(β_obs ~ R, maxit = maxit)` once per sample (a robust IRLS linear fit), then sets any negative fitted coefficient to zero and divides by their sum so the result always sums to exactly 1, **regardless of the `constraint` argument** — `constraint` is validated by `epidish()`'s top-level `match.arg()` call but is never read inside the `RPC` branch of `epidish()`'s own dispatch. Deterministic; `maxit` (this module's `ct_adv_maxit`, default 50) is `MASS::rlm()`'s own IRLS iteration cap.

**CBS (`method = "CBS"`)** — `EpiDISH:::DoCBS()` calls `e1071::svm(x = R, y = β_obs, type = "nu-regression", kernel = "linear", nu = each value in nu.v)` once per sample per candidate `nu`, then, for each sample, keeps whichever `nu` minimized the reconstruction RMSE (`ref %*% t(estF)` vs. the observed data) on that sample — the `nu.v` grid this module exposes as `ct_adv_nu1`/`nu2`/`nu3` (default `0.25/0.5/0.75`) is therefore a per-sample-selected tuning grid, not a single fixed hyperparameter. Coefficients are clipped negative-to-zero and renormalized to sum to 1, same as RPC, and — again — `constraint` is not consulted for this method. Deterministic (linear-kernel nu-SVR via `e1071`'s LIBSVM backend has no random initialization for this formulation).

**hepidish (`ct_method = "hepidish"`)** — `EpiDISH::hepidish(beta.m, ref1.m, ref2.m, h.CT.idx, method, maxit, nu.v, constraint)`, wrapped by `methyl_ct_run_hepidish()` (`celltype.R:329-350`). Internally it runs the *same* method (CP/RPC/CBS) independently against `ref1` and `ref2`, then replaces column `h.CT.idx` of the `ref1` result with `frac1[, h.CT.idx] * frac2` (element-wise product of the top-level component's fraction and the second-stage sub-fractions), concatenated with the other `ref1` columns unchanged. **This codebase's own call site hard-codes `method = "RPC"` for hepidish** (`mod_methyl_celltype.R:718-720`) — there is no UI control letting a user pick CP or CBS as hepidish's internal method; see `07_code_audit_findings.md`.

**Reproducibility / the "Random seed" control.** `decon_result()` calls `set.seed(input$ct_adv_seed %||% 1234)` immediately before every estimation call (`mod_methyl_celltype.R:711-712`). Having independently read `DoCP`/`DoRPC`/`DoCBS`'s source, none of the three call any R random-number-generating function (`runif`/`rnorm`/`sample`/similar) — `quadprog::solve.QP()`, `MASS::rlm()`, and `e1071::svm()`'s linear-kernel nu-SVR solve are all deterministic for a fixed input matrix. **The seed control therefore currently has no observable effect on the estimated fractions for any of CP, RPC, CBS, or hepidish** — re-running with a different seed value and identical data/reference/method will reproduce bit-identical results. This is not a defect in the seed-setting code itself (it is a harmless, forward-compatible no-op — if a future EpiDISH version or a reference-free method introduced stochastic elements, the seed would then start mattering), but the UI presents it as an "Advanced Parameter" without clarifying that none of the currently available methods are stochastic.

## 5.3 Reconstruction-validation math

`methyl_ct_reconstruct(ref_mat, fractions)` (`celltype.R:356-359`) computes the forward mixture model directly: `reconstructed = R[, ct] %*% t(fractions[, ct])`, restricted to the cell types common to both matrices. This is the same linear model §5.1 describes, run forward with the *estimated* `w` rather than solved backward for it — a self-consistency check, not an independent validation against ground truth (there is no ground-truth cell count in this pipeline; "validation" here means "how well does the fitted mixture explain the observed bulk signal," not "how accurate are the fractions against a gold standard").

`methyl_ct_validation_metrics(observed, reconstructed)` (`celltype.R:361-386`) aligns the two matrices on common CpGs/samples, drops any row incomplete in either, and computes:
- **Overall Pearson correlation** (`stats::cor()`) between the flattened observed and reconstructed values.
- **RMSE** = `sqrt(mean((obs − rec)²))`.
- **MAE** = `mean(|obs − rec|)`.
- **R²** = `1 − Σ(obs − rec)² / Σ(obs − mean(obs))²`, the standard coefficient-of-determination definition, computed against the observed data's own mean (not a null model's).
- The identical four metrics computed per sample (one row per sample, `vapply()` over columns).

## 5.4 Marker-ranking math (`methyl_ct_marker_rank()`, `celltype.R:160-193`)

For reference matrix `R` (CpG × cell type) and a target cell type `j`:
- `o_max`/`o_min` = the row-wise max/min of every *other* cell type's centroid at that CpG.
- `eff_hyper = R[, j] − o_max` (how much higher `j`'s methylation is than the highest of the others).
- `eff_hypo = o_min − R[, j]` (how much lower `j`'s methylation is than the lowest of the others).
- `effect = max(eff_hyper, eff_hypo)` — the larger of the two, i.e. however `j` is most extreme relative to every other cell type, in whichever direction.
- `direction = "hyper"` if `eff_hyper ≥ eff_hypo`, else `"hypo"`.
- `specificity = effect / SD(other cell types' centroids at that CpG)` (SD floored at `1e-6` to avoid division by zero) — normalizes the raw effect size by how tightly clustered the *other* cell types already are, so a CpG that is merely "somewhat different" from a tightly-clustered group of other types scores lower than one that is separated by the same raw magnitude from a group that is already spread out.
- Each CpG is assigned to whichever cell type `j` gives it the largest `effect` across the loop over all `j`.
- `btw_type_var` = `stats::var()` of that CpG's centroid values across *all* cell types (a simple between-group variance, independent of which type "won" the assignment above) — this is what "Variance-based CpGs" sorts by, as an alternative to sorting by `effect`.

**Why no p-value/FDR is computed, and why that is the scientifically honest choice.** A p-value for "is this CpG differentially methylated between cell type `j` and the others" requires a null distribution derived from per-sample (or per-replicate) variability within each group — the standard deviation *of the estimate*, not just the spread of the *other groups' point estimates*. The reference centroids this function operates on are single mean-beta values per cell type (the shape every EpiDISH built-in reference panel ships in, and the shape this module's custom-reference upload also requires — verified: `methyl_ct_parse_custom_reference()` accepts exactly one row per CpG per cell-type *column*, with no mechanism for multiple replicate columns per cell type). With no per-sample replicate data, there is no valid way to estimate sampling variability for a real hypothesis test, and any p-value computed anyway (e.g. by inventing a placeholder sample size) would be fabricated. The module's own UI is verified consistent with this: the marker-selection results table hardcodes the `P-value`/`FDR` columns to `"n/a"`/`"n/a (no replicates)"` literal strings (`mod_methyl_celltype.R:626`) rather than computing anything, and an explicit note states "Maximum FDR: not applicable — reference centroids have no per-sample replicates" (`mod_methyl_celltype.R:118-119`). This is judged scientifically sound and consistently enforced.

## 5.5 Group-comparison statistics (`methyl_ct_group_stats()`, `celltype.R:445-472`)

Unlike marker ranking (which has no per-sample data to test), the *estimated fraction matrix* genuinely does have one value per sample per cell type, so a real hypothesis test against an uploaded phenotype grouping column is valid here.

- **2 groups:** `stats::wilcox.test(g1, g2)` (Mann-Whitney/Wilcoxon rank-sum, non-parametric — appropriate since cell-type fractions are bounded in `[0,1]` and often non-normally distributed, especially near 0). Effect size is the **rank-biserial correlation**, computed as `1 − 2·W / (n1·n2)` where `W` is the Wilcoxon statistic — a standard, distribution-free effect-size measure ranging from −1 to 1.
- **>2 groups:** `stats::kruskal.test(x, group)` (the non-parametric generalization of one-way ANOVA). Effect size is **epsilon-squared**, computed as `(H − (k−1)) / (n − k)` where `H` is the Kruskal-Wallis statistic, `k` the number of groups, and `n` the number of non-missing observations — a standard, bounded-below-by-0 effect-size measure for this test.
- One test is run per cell type, and `stats::p.adjust(method = "BH")` (Benjamini-Hochberg) is applied **across cell types** within one run — controlling the false discovery rate over the family of cell-type comparisons for that single grouping variable, the correct scope given the analysis structure (one grouping variable tested against every cell type at once).

## 5.6 What this module does *not* implement, restated precisely

- **Houseman's original 2012 R implementation** (as published, using `minfi`'s `estimateCellCounts()` machinery) — this codebase uses `EpiDISH`'s own from-scratch CP reimplementation of the same constrained-projection idea, not the original Houseman code.
- **MethylResolver, IDOL-optimized libraries, and true reference-free deconvolution** — surfaced as disabled options with an installation reason, never approximated or faked (`methyl_ct_unavailable_methods()`, `celltype.R:68-77`).
- **A significance test for marker CpGs** — deliberately omitted for the reasons in §5.4, not an oversight.
- **A user-selectable internal method for hepidish** — always RPC in this codebase's call site, not user-configurable (§5.2).
