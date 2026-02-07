# rbmiUtils — Reporting & Robustness Milestone

## What This Is

An R package extending rbmi (reference-based multiple imputation) for clinical trial analysis workflows. rbmiUtils provides utilities for data preparation, analysis execution across imputations, results tidying, and formatting. This milestone adds a reporting layer that bridges rbmi results into the cards/ARD ecosystem for clinical trial table and figure generation, improves print/summary methods for key objects, refactors the analyse_mi_data() wrapper, and hardens all recent additions for production use.

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

### Active

- [ ] Convert tidy pool results to ARD format (cards package)
- [ ] Generate efficacy summary tables using gtsummary + gt (regulatory style)
- [ ] Create forest plot function (treatment effects across visits/subgroups with CIs)
- [ ] Create responder bar chart function (proportion responding by arm and visit)
- [ ] Improve print/summary methods for rbmi pool objects
- [ ] Improve print/summary methods for rbmi draws/impute objects
- [ ] Improve print/summary methods for the analysis class
- [ ] Refactor analyse_mi_data() to wrap rbmi::analyse() instead of reimplementing internals
- [ ] Harden gcomp functions (edge cases, input validation, error messages)
- [ ] Harden storage functions (reduce/expand — type coercion, attribute preservation)
- [ ] Harden data preparation functions (validation gaps, factor handling)
- [ ] Harden tidier/formatter functions (parameter parsing fragility, edge cases)

### Out of Scope

- Interactive/Shiny dashboards — focus is on static report generation
- Word/PowerPoint output via flextable — gt + HTML/PDF is sufficient for this milestone
- Custom pooling methods — use rbmi's built-in pooling
- Spaghetti/trajectory plots — not requested for this milestone
- Safety tables (AE listings) — outside rbmi efficacy workflow

## Context

- **Existing codebase:** 23 exported functions across 6 layers (data prep, analysis, utilities, tidying, formatting, storage)
- **Test coverage:** 9 test files with ~159 test blocks, good coverage of core functions
- **Known fragility:** tidy_pool_obj() parameter parsing relies on "_" separator; beeca output format coupling; rbmi class hierarchy dependency not version-pinned
- **Known gaps:** No validation of delta data subject-visit uniqueness; silent defaults for vars$strategy; formula construction without input sanitization
- **cards/cardx ecosystem:** Analysis Results Datasets (ARD) are the pharmaverse standard for representing analysis results as structured data frames. cardx provides extensions for common statistical analyses. gtsummary consumes ARDs to produce tables.
- **rbmi analyse() wrapper:** Current analyse_mi_data() contains internal helper functions (extract_covariates2, as_simple_formula2, as_analysis2) copied from rbmi. Wrapping rbmi::analyse() directly would reduce maintenance burden and drift risk.

## Constraints

- **Tech stack**: R package (CRAN-compatible), must pass R CMD check on macOS/Windows/Ubuntu
- **Dependencies**: cards, cardx, gtsummary, gt as Suggests vs Imports — to be decided during planning
- **Compatibility**: rbmi >= 1.4, R >= 4.1
- **Testing**: testthat edition 3, must maintain existing test coverage and add new tests for all additions
- **CI/CD**: GitHub Actions (R-CMD-check, test-coverage, pkgdown)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| ARD via cards (not custom format) | Pharmaverse standard, flows into gtsummary ecosystem | — Pending |
| gtsummary + gt for tables (not flextable) | Natural fit with ARD pipeline, HTML/PDF output sufficient | — Pending |
| ggplot2 for figures | Already suggested dependency, standard for R publication figures | — Pending |
| Wrap rbmi::analyse() (not reimplement) | Reduces drift risk, maintenance burden | — Pending |
| Dependency classification (Suggests vs Imports) | Affects installation footprint — defer to planning phase | — Pending |

---
*Last updated: 2026-02-07 after initialization*
