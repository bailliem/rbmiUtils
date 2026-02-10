# Project Research Summary

**Project:** rbmiUtils v3 milestone — ARD Enrichment & Publication Polish
**Domain:** Clinical trial R package enhancement — MI-specific metadata and regulatory reporting
**Researched:** 2026-02-10
**Confidence:** HIGH

## Executive Summary

rbmiUtils v3 adds publication-quality polish and MI-specific diagnostic metadata to an established clinical trial analysis package. The research reveals that this milestone can be delivered entirely with the existing dependency set — no new packages needed. The three feature categories (MI diagnostics in ARD, imputation diagnostic helpers, publication styling) touch different architectural layers with minimal coupling, making parallel development feasible.

The core technical challenge is computing MI diagnostics (FMI, lambda, relative increase in variance) that rbmi's pool() function calculates internally but discards. These must be recomputed from the analysis object's per-imputation estimates using Rubin's rules formulas — a straightforward 40-60 lines of R code. The diagnostic helpers (describe_draws, describe_imputation) are standalone functions that inspect rbmi objects but don't modify existing workflows. Publication styling refinements are pure visual polish with no data flow changes.

The primary risk is coupling to rbmi internal object structures. While the per-imputation data access pattern for MI diagnostics is stable (used by rbmi's own pool() function), the draws and imputation object structures accessed by diagnostic helpers are implementation details. Defensive programming with null-checks and method validation guards against fragility. The second risk is ARD validation — new MI diagnostic rows must preserve list-column types to pass cards::check_ard_structure(). Following the existing pool_to_ard() construction pattern prevents this.

## Key Findings

### Recommended Stack

**No new dependencies required.** All v3 features are implementable with the current stack. The existing dependencies (rbmi, cards, gt, ggplot2, patchwork, cli) provide all needed capabilities. Rubin's rules formulas are 6 lines of base R arithmetic — adding mice or mitml for pool.scalar() would introduce 30+ transitive dependencies for functionality we can compute directly.

**Core technologies:**
- **rbmi (>= 1.4):** Provides analysis object with per-imputation estimates/SEs — stable access pattern
- **cards (>= 0.4.0):** ARD validation via check_ard_structure() — accepts custom stat_name values
- **gt (>= 0.10.0):** All publication styling features available via tab_options() and tab_style()
- **ggplot2 (>= 3.4.0):** Standard API unaffected by 4.0 S7 migration — backward compatible
- **cli (>= 3.6.0):** Formatted console output for diagnostic helpers

**Version compatibility:** All packages current on CRAN as of 2026-02-10. ggplot2 4.0+ S7 migration does not affect rbmiUtils (uses only public geom/scale/theme API). No code changes needed.

### Expected Features

**Must have (table stakes):**
- **FMI in ARD** — The single most important MI diagnostic. Regulatory reviewers expect it. mice::pool() reports it. Users cannot compute it from rbmi pool objects.
- **Relative increase in variance (RIV) in ARD** — Standard MI diagnostic reported alongside FMI.
- **Lambda (proportion of variance due to missingness) in ARD** — Third standard diagnostic from Rubin's rules.
- **Pooling method identification in ARD** — Already present as method row; verify standardization.
- **describe_draws()** — Imputation model summary (method, samples, failures, formula, MCMC diagnostics for Bayes).
- **describe_imputation()** — Imputation statistics (method, M, subjects imputed, visits with missing data).

**Should have (competitive):**
- **Within/between/total variance components in ARD** — Full Rubin's rules decomposition (V_w, V_b, V_t). mice reports ubar/b/t.
- **Adjusted degrees of freedom in ARD** — Barnard-Rubin df for each parameter. Important for small-sample MI.
- **Relative efficiency (RE) in ARD** — RE = 1/(1 + FMI/M). Assesses imputation adequacy.
- **Efficacy table font/spacing refinements** — Publication defaults via theme_regulatory() helper.
- **Forest plot font/spacing/alignment refinements** — Expose font_family and panel_widths parameters.
- **README with realistic clinical workflow** — Full ADEFF → draws → impute → analyse → pool → table/plot with clinically meaningful data.
- **Function documentation examples with real usage** — Replace minimal examples with ADMI/ADEFF data.

**Defer (v2+):**
- Sensitivity analysis comparison features (side-by-side tables, forest overlays)
- Generic imputation diagnostic plots (trace plots, density overlays — specialized visualization)
- BMLMI-specific MI diagnostics (requires BMLMI-modified Rubin's rules)
- gtsummary tbl_ard_*() integration (ARD format doesn't match gtsummary's assumptions for regulatory tables)

### Architecture Approach

The v3 features integrate into the existing 7-layer architecture at three distinct points with clear boundaries. No changes to Layers 1-4 (data preparation, analysis execution, results processing). All features are additive.

**Major components:**
1. **MI diagnostics computation (new Layer 5 enhancement)** — compute_mi_diagnostics() internal helper recomputes FMI/lambda/RIV from analysis object per-imputation estimates. Uses Rubin's rules formulas (40-60 LOC). Integrates into pool_to_ard() via optional analysis_obj parameter.
2. **Diagnostic helpers (new cross-cutting layer)** — describe_draws() and describe_imputation() standalone functions. Accept rbmi objects, return rbmiUtils-owned summary classes with S3 print methods. No coupling to analysis pipeline.
3. **Publication styling (Layer 6 refinements)** — theme_regulatory() internal helper for gt tables, theme_forest() extension for ggplot2 plots. Add optional parameters, preserve defaults for backward compatibility.

**Data flow:**
- MI diagnostics: analysis_obj → compute_mi_diagnostics() → pool_to_ard(pool_obj, analysis_obj) → enriched ARD
- Diagnostic helpers: draws_obj → describe_draws() → summary (print), impute_obj → describe_imputation() → summary (print)
- Styling: pool_obj → efficacy_table(theme = "regulatory") → themed table, pool_obj → plot_forest(font_family, widths) → polished plot

### Critical Pitfalls

1. **Adding FMI rows breaks check_ard_structure() validation** — New ARD rows must use I(list(...)) pattern for stat/fmt_fn/warning/error columns. A single non-list value corrupts the entire data frame. Prevention: Follow existing pool_to_ard() construction pattern, validate with check_ard_structure() after adding rows.

2. **FMI not available in rbmi pool object** — rbmi computes var_w, var_b, lambda internally but discards them. Must recompute from analysis object's per-imputation estimates/SEs. Prevention: Accept analysis_obj parameter in pool_to_ard(), compute via Rubin's rules formulas, validate against rbmi::rubin_rules() output.

3. **describe_draws() couples to rbmi internal object structure** — draws object fields (samples, fit, method) lack stability guarantee. Prevention: Minimize internal field access, wrap in tryCatch, use rstan public API for MCMC diagnostics, test with multiple rbmi versions.

4. **describe_imputation() memory exhaustion with full extraction** — Extracting all M imputed datasets for diagnostics is wasteful. Prevention: Accept stacked ADMI data (already in memory), sample subset of imputations if extracting from impute object.

5. **Enriched pool_to_ard() API signature change breaks existing callers** — Adding required analysis_obj parameter is breaking. Prevention: Make analysis_obj optional (default NULL), return base ARD when omitted, enriched ARD when provided. Backward compatible.

## Implications for Roadmap

Based on research, suggested phase structure mirrors architectural boundaries and risk profile:

### Phase 1: MI Diagnostic Statistics (Core Value)

**Rationale:** Highest-value work filling a real gap. MI diagnostics (FMI, lambda, RIV) are not available from rbmi but are essential for regulatory reporting. Must come first because the diagnostic computation is the foundation for ARD enrichment. Critical risk area (Pitfalls 1, 2, 5, 9, 12) requiring careful API design.

**Delivers:**
- compute_mi_diagnostics() internal helper
- Enriched pool_to_ard() with optional analysis_obj parameter
- ARD with FMI, lambda, RIV, V_w, V_b, V_t, df_adjusted, RE as new stat_name rows

**Addresses:** TS-1/2/3 (FMI/RIV/Lambda), D-1/2/3 (variance components/df/RE), TS-4 (method standardization)

**Avoids:** Pitfall 1 (ARD validation), Pitfall 2 (recomputation strategy), Pitfall 5 (API signature)

**Research flags:** Standard patterns (Rubin's rules formulas well-documented). Skip research-phase.

### Phase 2: Describe Helpers (New Functions)

**Rationale:** Self-contained diagnostic functions with no downstream consumers in package. Moderate complexity but high documentation value. Independent of ARD enrichment — can be built in parallel with Phase 1 but sequenced after for focus. Critical risk: coupling to rbmi internals (Pitfall 3).

**Delivers:**
- describe_draws(draws_obj) — method, samples, failures, formula, MCMC diagnostics (ESS, Rhat)
- describe_imputation(imputed_data, original_data, vars) — method, M, missingness by visit/arm, reference mapping
- S3 print methods for rbmiUtils-owned summary classes

**Addresses:** TS-5 (describe_draws), TS-6 (describe_imputation)

**Avoids:** Pitfall 3 (rbmi internals coupling via defensive access), Pitfall 4 (memory via stacked data input)

**Research flags:** Skip research-phase (rbmi object structures verified, API design complete).

### Phase 3: Publication Polish (Visual Refinements)

**Rationale:** Lower complexity, moderate user-facing impact. Must come after functional changes to avoid regenerating images twice. Independent of Phases 1-2 — pure styling with no data flow changes. Moderate risk (Pitfalls 6, 7) but visually impactful.

**Delivers:**
- theme_regulatory() internal helper for gt styling
- efficacy_table() with theme/font_size/row_striping parameters
- theme_forest() extension with font_family control
- plot_forest() with font_family and panel_widths parameters

**Addresses:** D-4 (efficacy table styling), D-5 (forest plot refinements)

**Avoids:** Pitfall 6 (styling changes) by adding parameters not changing defaults, Pitfall 7 (patchwork alignment) by exposing widths parameter

**Research flags:** Skip research-phase (gt and ggplot2 APIs well-documented).

### Phase 4: Documentation Overhaul (Last)

**Rationale:** Must come after all functional changes are complete. Regenerate images, update examples, rebuild pkgdown site. Lower risk (Pitfalls 8, 10, 13) but process-heavy.

**Delivers:**
- README overhaul with realistic clinical workflow (ANCOVA + binary endpoints)
- Function documentation examples with ADMI/ADEFF data
- Vignette updates for MI diagnostics
- Pre-rendered images regenerated
- pkgdown site rebuild
- NEWS.md updates

**Addresses:** D-6 (README workflow), D-7 (documentation examples)

**Avoids:** Pitfall 8 (pkgdown build) by using pre-rendered images and testing both GitHub and pkgdown rendering

**Research flags:** Skip research-phase (documentation patterns, not domain research).

### Phase Ordering Rationale

- **Functional correctness before visual polish before documentation** — Each layer depends on previous being stable.
- **MI diagnostics first** — Foundation for ARD enrichment. Critical risk area requiring careful design. Highest user value.
- **Describe helpers second** — Self-contained, no downstream coupling. Can validate independently.
- **Styling third** — Avoid regenerating images twice. Add parameters, preserve defaults.
- **Documentation last** — Documents finalized features. Regenerate all images once.

**Phase dependencies:** Phases 2-3 can technically run in parallel after Phase 1 (describe helpers and styling are independent). Phase 4 depends on all.

### Research Flags

**Phases with standard patterns (skip research-phase):**
- **Phase 1:** Rubin's rules formulas well-documented (Rubin 1987, van Buuren 2018, mice package). cards ARD format verified.
- **Phase 2:** rbmi object structures verified via source inspection. cli formatting patterns established.
- **Phase 3:** gt and ggplot2 APIs well-documented. Styling patterns verified.
- **Phase 4:** Documentation patterns, not domain research.

**No phases need /gsd:research-phase.** All research complete. Ready for requirements definition.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All features verified against existing package APIs; formulas are base R arithmetic |
| Features | HIGH | Feature landscape mapped from codebase inspection, rbmi internals, mice benchmarking |
| Architecture | HIGH | Integration points identified; no changes to core layers; clear boundaries |
| Pitfalls | MEDIUM-HIGH | Critical pitfalls grounded in source code analysis; memory/platform concerns inferred |

**Overall confidence:** HIGH

### Gaps to Address

**During Phase 1 implementation:**
- Validate computed FMI values against mice::pool() output for cross-validation
- Test with all pooling methods (rubin, bootstrap, jackknife, bmlmi) — only rubin should return MI diagnostics
- Verify var_t from manual computation matches pool_obj$pars[[x]]$se^2

**During Phase 2 implementation:**
- Confirm rbmi draws object structure across rbmi 1.4, 1.5, 1.6 (multiple version testing)
- Test MCMC diagnostics extraction with actual Stan fits (requires rstan installation)
- Validate memory behavior of describe_imputation() with M=100, M=1000 imputation scenarios

**During Phase 3 implementation:**
- Test styled tables/plots at multiple output sizes and DPI settings (300 DPI for publication)
- Verify cross-platform font rendering (macOS vs Windows)
- Check patchwork alignment with various text_size and font_family values

**During Phase 4 implementation:**
- Validate pkgdown image paths render correctly on both GitHub and built site
- Test README rendering on GitHub (markdown) and pkgdown (HTML)

## Sources

### Primary (HIGH confidence)

**Codebase inspection:**
- rbmiUtils R/ source files (ard_conversion.R, efficacy_table.R, plot_forest.R, pool_methods.R, tidiers.R, analyse_mi_data.R, utils.R)
- rbmiUtils tests/testthat/ (test-ard_conversion.R, test-pool_methods.R)
- rbmiUtils .planning/PROJECT.md, .planning/MILESTONES.md

**rbmi internals (source inspection):**
- rbmi:::rubin_rules — computes var_w, var_b, var_t, df; returns est_point, var_t, df
- rbmi:::rubin_df — computes lambda internally; returns df only
- rbmi:::pool_internal.rubin — calls rubin_rules then parametric_ci
- rbmi:::print.draws — shows formula, method, samples, failures
- rbmi:::print.imputation — shows datasets count, missingness %, references
- [rbmi pool.R source](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R)
- [rbmi draws.R source](https://github.com/insightsengineering/rbmi/blob/main/R/draws.R)

**Package documentation:**
- [cards 0.7.1 CRAN](https://cran.r-project.org/web/packages/cards/index.html) — check_ard_structure() validation
- [cards check_ard_structure source](https://raw.githubusercontent.com/insightsengineering/cards/main/R/check_ard_structure.R)
- [gt 1.3.0 CRAN](https://cran.r-project.org/web/packages/gt/index.html) — tab_options() reference
- [ggplot2 4.0.2 CRAN](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [rbmi 1.6.0 CRAN](https://cran.r-project.org/web/packages/rbmi/index.html)

### Secondary (MEDIUM confidence)

**Statistical formulas:**
- [Rubin's Rules — Book of MI (Heymans)](https://bookdown.org/mwheymans/bookmi/rubins-rules.html)
- [Measures of Missing Data Information (Heymans)](https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html)
- [mice::pool documentation](https://amices.org/mice/reference/pool.html)
- [van Buuren: MI in a Nutshell](https://stefvanbuuren.name/fimd/sec-nutshell.html)
- Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
- Barnard, J. and Rubin, D.B. (1999). Small-sample degrees of freedom with multiple imputation. Biometrika, 86(4), 948-955.

**Package ecosystem:**
- [ggplot2 4.0.0 announcement](https://tidyverse.org/blog/2025/09/ggplot2-4-0-0/) — S7 migration compatibility
- [patchwork layout guide](https://patchwork.data-imaginist.com/articles/guides/layout.html)
- [gt tab_options() reference](https://gt.rstudio.com/reference/tab_options.html)
- [pkgdown reference documentation](https://pkgdown.r-lib.org/articles/pkgdown.html)
- [R Packages (2e)](https://r-pkgs.org/) — dependency management, README/pkgdown

### Tertiary (LOW confidence)

**Web search results:**
- [R Consortium: cards package and ARD standard](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/)
- [Introduction to clinical tables with gt](https://www.r-bloggers.com/2024/02/introduction-to-clinical-tables-with-the-gt-package/)
- [PharmaSUG 2023 QT-263: R Tables via GT for Regulatory Submissions](https://pharmasug.org/proceedings/2023/QT/PharmaSUG-2023-QT-263.pdf)

---
*Research completed: 2026-02-10*
*Ready for roadmap: yes*
