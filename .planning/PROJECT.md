# rbmiUtils — Reporting & Robustness

## What This Is

An R package extending rbmi (reference-based multiple imputation) for clinical trial analysis workflows. rbmiUtils provides utilities for data preparation, analysis execution across imputations, results tidying, formatting, and a reporting layer that bridges rbmi results into publication-ready regulatory tables and forest plots via the pharmaverse ARD ecosystem. Includes MI diagnostic metadata enrichment and pipeline introspection helpers for regulatory review workflows.

## Core Value

Clinical trial results from rbmi flow seamlessly into publication-ready regulatory tables and figures — no manual wrangling between pooled results and final output.

## Requirements

### Validated

- ✓ Data validation before imputation (validate_data) — existing
- ✓ ICE preparation for intercurrent events (prepare_data_ice) — existing
- ✓ Missingness summarization (summarise_missingness) — existing
- ✓ Imputed data extraction with ID mapping (get_imputed_data) — existing
- ✓ Analysis execution across imputations (analyse_mi_data) — existing
- ✓ Binary outcome g-computation (gcomp_responder, gcomp_binary) — existing
- ✓ Multi-visit g-computation (gcomp_responder_multi) — existing
- ✓ Pooled results tidying (tidy_pool_obj) — existing
- ✓ P-value and estimate formatting (format_pvalue, format_estimate) — existing
- ✓ Results table formatting (format_results_table, format_results) — existing
- ✓ Result combination and extraction (combine_results, extract_trt_effects, extract_lsm) — existing
- ✓ Imputation storage optimization (reduce_imputed_data, expand_imputed_data) — existing
- ✓ IMPID creation for non-rbmi imputations (create_impid) — existing
- ✓ Convert tidy pool results to ARD format (pool_to_ard) — v1
- ✓ Generate efficacy summary tables using gt (efficacy_table) — v1
- ✓ Create forest plot function (plot_forest) — v1
- ✓ Improve print/summary methods for pool objects — v1
- ✓ Improve print/summary methods for analysis class — v1
- ✓ Refactor analyse_mi_data() to use inherits() and deprecate internal helpers — v1
- ✓ Harden gcomp functions (input validation, beeca output pinning) — v1
- ✓ Harden storage functions (round-trip digest verification) — v1
- ✓ Harden tidier functions (regex-based parameter parsing) — v1
- ✓ Harden data preparation functions (cli messaging, interaction term validation, factor handling, edge cases) — v2
- ✓ End-to-end clinical trial analysis vignette (rbmi → rbmiUtils → table + forest plot) — v2
- ✓ README enhancement with visual teaser and pipeline overview — v2
- ✓ Full pkgdown site polish (navbar, grouped references, hex logo, social cards, footer) — v2
- ✓ Rendered examples for plot_forest() and efficacy_table() — v2
- ✓ NEWS.md version history — v2
- ✓ Inline cross-references to rbmi and beeca in vignettes — v2
- ✓ FMI, lambda, RIV, Barnard-Rubin df, and relative efficiency in ARD output (pool_to_ard with analysis_obj) — v3
- ✓ Pooling method as stat_name row in ARD output — v3
- ✓ MI diagnostics omitted for non-Rubin pooling methods with cli message — v3
- ✓ describe_draws() for method, samples, failures, formula, MCMC diagnostics — v3
- ✓ describe_imputation() for method, M, missingness by visit/arm — v3
- ✓ Structured S3 objects with cli-formatted print methods for describe helpers — v3
- ✓ Font family and size controls for efficacy_table() — v3
- ✓ Row padding control for efficacy_table() — v3
- ✓ Font family and panel width controls for plot_forest() — v3
- ✓ Forest plot left-panel alignment fix — v3
- ✓ README with complete ADEFF pipeline workflow — v3
- ✓ Function examples using ADMI/ADEFF package datasets — v3
- ✓ MI diagnostics and pipeline inspection vignette — v3
- ✓ Pre-rendered images reflecting v3 styling — v3
- ✓ NEWS.md v0.3.0 entries — v3

### Active

(None — planning next milestone)

### Deferred (v4+)

- Responder bar chart function (proportion responding by arm and visit)
- Forest plot with sensitivity analysis overlay
- Responder chart with treatment difference annotations
- Column formatting controls for gt theming
- Sensitivity analysis comparison table
- as_gt() / as_gtsummary() S3 methods for pool objects
- BMLMI-specific MI diagnostics (modified Rubin's rules)
- Imputation diagnostic plots (trace, density, stripplot)

### Out of Scope

- Interactive/Shiny dashboards — focus is on static report generation
- Word/PowerPoint output via flextable — gt + HTML/PDF is sufficient
- Custom pooling methods — use rbmi's built-in pooling
- Spaghetti/trajectory plots — not requested
- Safety tables (AE listings) — outside rbmi efficacy workflow
- rtables/tern integration — different ecosystem, incompatible with ARD paradigm
- Re-implementing rbmi::pool() to store FMI — utility layer, not a fork
- Automatic quality thresholds for MI diagnostics — context-dependent, let user interpret
- MI diagnostics for non-Rubin methods — FMI/lambda/RIV are Rubin-specific

## Context

- **Shipped:** v1 Reporting & Robustness (2026-02-08), v2 Documentation & Hardening (2026-02-10), v3 ARD Enrichment & Polish (2026-02-11)
- **Codebase:** 32 exported functions across 8 layers (data prep, analysis, utilities, tidying, formatting, storage, reporting, introspection)
- **Source code:** 4,980 lines R, 6,668 lines tests
- **Test coverage:** 14+ test files, 400+ new v3 expectations, comprehensive coverage of all functions
- **Dependencies:** cli, lifecycle as Imports; cards, gt, ggplot2, patchwork as Suggests with dependency guards
- **cards/cardx ecosystem:** pool_to_ard() produces valid ARD with MI diagnostics passing cards::check_ard_structure()
- **Documentation:** 5 vignettes (pipeline, analysis, data prep, storage, MI diagnostics); visual README with ADEFF pipeline Quick Start; versioned NEWS.md (0.3.0/0.2.0/0.1.0)
- **Site:** pkgdown at openpharma.github.io/rbmiUtils with hex logo, 10-group reference index (including Introspection), navbar, open graph cards, pharmaverse footer
- **Known gaps:** No validation of delta data subject-visit uniqueness; table images not regenerated (Chromium unavailable, existing images still valid)
- **Tech debt:** Table pre-rendered images from v2 (minor; functional but don't reflect font_size/row_padding defaults)

## Constraints

- **Tech stack**: R package (CRAN-compatible), must pass R CMD check on macOS/Windows/Ubuntu
- **Dependencies**: cards, gt, ggplot2, patchwork as Suggests (guarded with requireNamespace)
- **Compatibility**: rbmi >= 1.4, R >= 4.1
- **Testing**: testthat edition 3, must maintain existing test coverage and add new tests for all additions
- **CI/CD**: GitHub Actions (R-CMD-check, test-coverage, pkgdown)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| ARD via cards (not custom format) | Pharmaverse standard, flows into gtsummary ecosystem | ✓ Good — passes check_ard_structure() |
| gt direct for tables (not gtsummary layer) | Custom table layout needs gt directly; gtsummary is for regression tables | ✓ Good — cleaner API |
| ggplot2 + patchwork for figures | Three-panel composition needs patchwork; ggplot2 standard for R publication figures | ✓ Good — customizable via & theme() |
| inherits() + deprecation (not full rbmi::analyse wrapper) | Achieves stability without reimplementation risk | ✓ Good — reduced maintenance burden |
| All reporting deps as Suggests | Keeps core package lightweight; dependency guards prevent installation failures | ✓ Good — clean install path |
| Two-pass regex parsing for parameters | Single regex can't handle both ANCOVA and gcomp formats cleanly | ✓ Good — handles underscores correctly |
| Okabe-Ito palette for LSM colors | Maximally distinguishable for color vision deficiency | ✓ Good — accessible default |
| Filled vs open circles for significance | Clear visual distinction without relying on color alone | ✓ Good — accessible |
| All-NA covariates warn+skip, not error | validate_data returns TRUE/error, cannot modify vars | ✓ Good — non-breaking |
| All-NA outcome is hard error | Nothing to model, analysis cannot proceed | ✓ Good — fail-fast |
| Pre-generated static images for README | Avoids slow MCMC during README rendering | ✓ Good — fast builds |
| Static images via \figure{} for help pages | gt tables cannot render in base R help | ✓ Good — works everywhere |
| Tutorial tone for pipeline vignette | analyse2 is the reference; pipeline.Rmd is getting-started | ✓ Good — clear roles |
| NEWS.md organized by version (0.3.0, 0.2.0, 0.1.0) | Tidyverse convention for pkgdown::build_news() | ✓ Good — standard |
| Inline hyperlinks, not callout boxes | Natural reading flow, avoids visual clutter | ✓ Good — clean prose |
| Keep original rbmiUtils.png alongside logo.png | Avoid breaking external links | ✓ Good — backward compat |
| Added show_pvalues parameter to plot_forest() | User feedback: cleaner without p-value panel | ✓ Good — flexible |
| Restructured README: code before images | Better narrative flow showing pipeline then output | ✓ Good — clearer |
| Curated MI diagnostics only (no V_w/V_b/V_t) | Keep ARD concise for regulatory review; variance components are internal | ✓ Good — clean ARD |
| mice-convention stat naming (dot-separated lowercase) | Consistent with established MI ecosystem conventions | ✓ Good — familiar to users |
| Non-Rubin pooling omits diagnostics entirely | FMI/lambda/RIV undefined for non-Rubin methods; NA rows would mislead | ✓ Good — principled |
| Adjusted FMI follows mice formula | (riv + 2/(df+3))/(1+riv) distinct from lambda; matches established R packages | ✓ Good — reproducible |
| base R for missingness tables (not dplyr) | Minimizes dependencies; expand.grid + loop is sufficient | ✓ Good — lightweight |
| cli output via stderr (type="message") | cli writes to stderr connection; tests capture via type="message" | ✓ Good — testable |
| NULL defaults for styling params | Backward compatible; existing output unchanged unless user opts in | ✓ Good — non-breaking |
| panel_widths c(3,4,1.5) default | Balanced three-panel layout preserving existing visual proportions | ✓ Good — backward compat |
| \donttest{} for ADMI, \dontrun{} for ADEFF examples | ADMI runs fast but slow; ADEFF requires MCMC unavailable in check | ✓ Good — pragmatic |

## Milestone History

- **v1 Reporting & Robustness** (2026-02-08) — Phases 1-4: Foundation hardening, print/summary, ARD, efficacy tables, forest plots
- **v2 Documentation & Hardening** (2026-02-10) — Phases 5-7: Data prep hardening, end-to-end vignette, pkgdown site polish
- **v3 ARD Enrichment & Polish** (2026-02-11) — Phases 8-11: MI diagnostic metadata, describe helpers, publication styling, documentation overhaul

See `.planning/MILESTONES.md` for full details.

---
*Last updated: 2026-02-11 after v3 milestone*
