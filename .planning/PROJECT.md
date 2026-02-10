# rbmiUtils — Reporting & Robustness

## What This Is

An R package extending rbmi (reference-based multiple imputation) for clinical trial analysis workflows. rbmiUtils provides utilities for data preparation, analysis execution across imputations, results tidying, formatting, and a reporting layer that bridges rbmi results into publication-ready regulatory tables and forest plots via the pharmaverse ARD ecosystem.

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

### Active

## Current Milestone: v3 ARD Enrichment & Polish

**Goal:** Enrich ARD output with MI-specific metadata, add imputation diagnostics helpers, and polish tables, plots, and documentation to publication quality.

**Target features:**
- MI-specific metadata in ARD (FMI, relative increase in variance, pooling method)
- describe_draws() summary diagnostics for draws objects
- describe_imputation() summary diagnostics for impute objects
- Efficacy table styling refinements for publication quality
- Forest plot visual refinements (fonts, spacing, alignment)
- README overhaul with realistic clinical trial workflow examples
- Documentation examples that reflect real usage patterns

### Deferred (v4+)

- Responder bar chart function (proportion responding by arm and visit)
- Forest plot with sensitivity analysis overlay
- Responder chart with treatment difference annotations
- Column formatting controls for gt theming
- Sensitivity analysis comparison table
- as_gt() / as_gtsummary() S3 methods for pool objects

### Out of Scope

- Interactive/Shiny dashboards — focus is on static report generation
- Word/PowerPoint output via flextable — gt + HTML/PDF is sufficient
- Custom pooling methods — use rbmi's built-in pooling
- Spaghetti/trajectory plots — not requested
- Safety tables (AE listings) — outside rbmi efficacy workflow
- rtables/tern integration — different ecosystem, incompatible with ARD paradigm

## Context

- **Shipped:** v1 Reporting & Robustness (2026-02-08), v2 Documentation & Hardening (2026-02-10)
- **Codebase:** 28 exported functions across 7 layers (data prep, analysis, utilities, tidying, formatting, storage, reporting)
- **Source code:** 4,146 lines R, 5,285 lines tests
- **Test coverage:** 14 test files, 95+ data prep tests, comprehensive coverage of all functions
- **Dependencies:** cli, lifecycle as Imports; cards, gt, ggplot2, patchwork as Suggests with dependency guards
- **cards/cardx ecosystem:** pool_to_ard() produces valid ARD passing cards::check_ard_structure()
- **Documentation:** 4 vignettes (pipeline, analysis, data prep, storage); visual README with rendered table/plot teasers; versioned NEWS.md (0.2.0/0.1.0)
- **Site:** pkgdown at openpharma.github.io/rbmiUtils with hex logo, 9-group reference index, navbar, open graph cards, pharmaverse footer
- **Known gaps:** No validation of delta data subject-visit uniqueness; show_pvalues parameter added late (not in original spec)
- **Tech debt:** Some pre-generated images may need regeneration if plot functions change

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
| NEWS.md organized by version (0.2.0, 0.1.0) | Tidyverse convention for pkgdown::build_news() | ✓ Good — standard |
| Inline hyperlinks, not callout boxes | Natural reading flow, avoids visual clutter | ✓ Good — clean prose |
| Keep original rbmiUtils.png alongside logo.png | Avoid breaking external links | ✓ Good — backward compat |
| Added show_pvalues parameter to plot_forest() | User feedback: cleaner without p-value panel | ✓ Good — flexible |
| Restructured README: code before images | Better narrative flow showing pipeline then output | ✓ Good — clearer |

## Milestone History

- **v1 Reporting & Robustness** (2026-02-08) — Phases 1-4: Foundation hardening, print/summary, ARD, efficacy tables, forest plots
- **v2 Documentation & Hardening** (2026-02-10) — Phases 5-7: Data prep hardening, end-to-end vignette, pkgdown site polish
- **v3 ARD Enrichment & Polish** (in progress) — ARD metadata, describe helpers, table/plot polish, documentation overhaul

See `.planning/MILESTONES.md` for full details.

---
*Last updated: 2026-02-10 after v3 milestone started*
