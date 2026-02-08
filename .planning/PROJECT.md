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

### Active

- [ ] Harden data preparation functions (validation gaps, factor handling, edge cases)
- [ ] End-to-end clinical trial analysis vignette (rbmi → rbmiUtils → table + forest plot)
- [ ] README enhancement with visual teaser and pipeline overview
- [ ] Full pkgdown site polish (navbar, grouped references, hex logo, social cards, footer)
- [ ] Rendered examples for plot_forest() and efficacy_table()
- [ ] NEWS.md version history
- [ ] Inline cross-references to rbmi and beeca in vignettes

### Deferred (v3+)

- Responder bar chart function (proportion responding by arm and visit)
- Forest plot with sensitivity analysis overlay
- Responder chart with treatment difference annotations
- Column formatting controls for gt theming
- Sensitivity analysis comparison table
- MI-specific metadata in ARD (FMI, pooling method)
- as_gt() / as_gtsummary() S3 methods for pool objects
- describe_draws() helper for draws objects
- describe_imputation() helper for impute objects

### Out of Scope

- Interactive/Shiny dashboards — focus is on static report generation
- Word/PowerPoint output via flextable — gt + HTML/PDF is sufficient
- Custom pooling methods — use rbmi's built-in pooling
- Spaghetti/trajectory plots — not requested
- Safety tables (AE listings) — outside rbmi efficacy workflow
- rtables/tern integration — different ecosystem, incompatible with ARD paradigm

## Context

- **Shipped:** v1 Reporting & Robustness (2026-02-08)
- **Codebase:** 28 exported functions across 7 layers (data prep, analysis, utilities, tidying, formatting, storage, reporting)
- **Source code:** 3,976 lines R, 4,916 lines tests
- **Test coverage:** 14 test files, comprehensive coverage of all functions
- **Dependencies:** cli, lifecycle as Imports; cards, gt, ggplot2, patchwork as Suggests with dependency guards
- **cards/cardx ecosystem:** pool_to_ard() produces valid ARD passing cards::check_ard_structure()
- **Known gaps:** No validation of delta data subject-visit uniqueness; silent defaults for vars$strategy; no vignettes for efficacy_table/plot_forest
- **Documentation:** 3 vignettes (analysis, data prep, storage); README with basic example; minimal pkgdown config; no hex logo; no NEWS.md
- **Site:** pkgdown deployed at openpharma.github.io/rbmiUtils but no function grouping, navbar, or visual polish

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

## Current Milestone: v2 Documentation & Hardening

**Goal:** Make rbmiUtils discoverable and trustworthy — polished site, end-to-end examples, and hardened data prep functions that give clear errors on bad input.

**Target features:**
- Harden data preparation functions (interaction term validation, factor handling, edge cases, batched warnings)
- End-to-end clinical trial vignette showing full pipeline with rbmi/beeca integration
- README visual teaser with table and forest plot output
- Full pkgdown site polish (hex logo, custom navbar, grouped references, social cards)
- NEWS.md, rendered function examples, inline cross-references to rbmi/beeca

---
*Last updated: 2026-02-08 after v2 milestone start*
