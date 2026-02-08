# Milestone v1: Reporting & Robustness

**Status:** SHIPPED 2026-02-08
**Phases:** 1-4
**Total Plans:** 9

## Overview

This roadmap delivered a reporting layer for rbmiUtils that bridges rbmi pooled results into publication-ready regulatory tables and figures via the pharmaverse ARD ecosystem. The work proceeded in four phases: first hardening fragile foundations (parameter parsing, class indexing, input validation), then adding print/summary methods and ARD conversion as infrastructure, then building efficacy tables consuming those ARDs, and finally adding forest plot visualization. Each phase delivered a coherent, testable capability.

## Phases

### Phase 1: Foundation Hardening

**Goal**: Existing functions are robust against edge cases, version drift, and malformed inputs so new layers build on stable foundations
**Depends on**: Nothing (first phase)
**Requirements**: HRD-01, HRD-02, HRD-03, HRD-04
**Success Criteria:**
  1. tidy_pool_obj() correctly parses parameter names containing underscores without data corruption
  2. analyse_mi_data() uses inherits() instead of fragile class()[[2]] indexing
  3. gcomp functions validate inputs and produce clear error messages for edge cases
  4. reduce/expand_imputed_data preserve column types through round-trip
  5. R CMD check passes with zero errors or warnings
**Plans**: 3 plans

Plans:
- [x] 01-01-PLAN.md — Add cli/lifecycle deps and fix tidy_pool_obj parameter parsing (HRD-01)
- [x] 01-02-PLAN.md — Refactor analyse_mi_data to use inherits() and deprecate internal helpers (HRD-02)
- [x] 01-03-PLAN.md — Harden gcomp input validation and storage round-trip digest (HRD-03, HRD-04)

**Completed:** 2026-02-07

### Phase 2: Print/Summary & ARD Conversion

**Goal**: Users get informative console output from key rbmi objects, and tidy pool results convert cleanly to the pharmaverse ARD standard for downstream table generation
**Depends on**: Phase 1 (stable tidy tibble output required for ARD conversion)
**Requirements**: PRT-01, PRT-02, PRT-03, PRT-04, ARD-01, ARD-02, ARD-03
**Success Criteria:**
  1. Printing a pool object shows rounded estimates, formatted CIs, and parameter labels
  2. summary() on a pool object produces a visit-level breakdown with significance flags
  3. Printing an analysis object shows parameter count, visits covered, and function name
  4. pool_to_ard() converts a pool object to a valid cards ARD data frame
  5. ARD output preserves visit, parameter type, and arm as grouping columns
**Plans**: 3 plans

Plans:
- [x] 02-01-PLAN.md — Create print.pool and summary.pool S3 methods with cli formatting (PRT-01, PRT-02)
- [x] 02-02-PLAN.md — Modernize print.analysis and summary.analysis with cli formatting (PRT-03, PRT-04)
- [x] 02-03-PLAN.md — Create pool_to_ard() ARD conversion function (ARD-01, ARD-02, ARD-03)

**Completed:** 2026-02-08

### Phase 3: Efficacy Tables

**Goal**: Users produce regulatory-style efficacy summary tables directly from rbmi pool objects with a single function call
**Depends on**: Phase 2 (ARD conversion layer required)
**Requirements**: TBL-01, TBL-02, TBL-03, TBL-04
**Success Criteria:**
  1. efficacy_table(pool_obj) produces a formatted table with LS means, treatment differences, CIs, and p-values by visit
  2. The table renders as gt output suitable for HTML and PDF
  3. The table includes footnotes documenting imputations, pooling method, and confidence level
  4. Users can override default formatting via function arguments
**Plans**: 2 plans

Plans:
- [x] 03-01-PLAN.md — Create efficacy_table() function with gt rendering and tests (TBL-01, TBL-02, TBL-03, TBL-04)
- [x] 03-02-PLAN.md — Add edge case handling, visit ordering, and integration tests

**Completed:** 2026-02-08

### Phase 4: Visualization

**Goal**: Users produce publication-quality forest plots of treatment effects across visits from rbmi pool objects
**Depends on**: Phase 1 (stable tidy tibble output); independent of Phases 2-3
**Requirements**: VIZ-01, VIZ-02, VIZ-03
**Success Criteria:**
  1. A forest plot function produces a figure showing treatment effect point estimates and CIs across visits
  2. The function returns a ggplot2 object that users can further customize
  3. The plot includes a reference line at zero by default, configurable to any user-specified value
**Plans**: 1 plan

Plans:
- [x] 04-01-PLAN.md — Implement plot_forest() with three-panel patchwork composition, tests, and visual verification (VIZ-01, VIZ-02, VIZ-03)

**Completed:** 2026-02-08

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Foundation Hardening | 3/3 | Complete | 2026-02-07 |
| 2. Print/Summary & ARD Conversion | 3/3 | Complete | 2026-02-08 |
| 3. Efficacy Tables | 2/2 | Complete | 2026-02-08 |
| 4. Visualization | 1/1 | Complete | 2026-02-08 |

---

## Milestone Summary

**Key Decisions:**
- 01-01-D1: Two-pass parsing (separate_wider_regex then mutate) for parameter names with underscores
- 01-02-D1: Suppress lifecycle deprecation warning for internal as_analysis2() calls
- 01-03-D1: Integrity check scoped to shared columns (intersect approach)
- 02-01-D1: cli output captured via withCallingHandlers for test verification
- 02-03-D1: as_card() direct construction over ard_identity() for batch ARD building
- 02-03-D2: Extracted is_cards_available() helper for testable dependency guards
- 03-01-D1: Letter-digit boundary spacing in visit label cleaning
- 04-01-D1: Horizontal orientation with visits on y-axis (clinical convention)
- 04-01-D2: Filled vs open circles for significance highlighting
- 04-01-D3: Okabe-Ito blue/vermilion for LSM arm colors

**Issues Resolved:**
- Fragile parameter parsing with underscore-containing names (tidy_pool_obj)
- Brittle class()[[2]] positional indexing (analyse_mi_data)
- Missing input validation on gcomp functions
- No integrity verification on storage round-trips
- No console output for pool/analysis objects
- No ARD interchange format for pharmaverse integration
- No regulatory-style table generation
- No forest plot visualization

**Issues Deferred:**
- Vignettes for efficacy_table() and plot_forest() (documentation gap)
- pool_to_ard() → gtsummary integration examples
- Responder bar charts (VIZ-04, v2)
- Sensitivity analysis overlays (VIZ-05, v2)
- MI-specific ARD metadata (ARD-04, v2)

**Technical Debt Incurred:**
- None critical. Package builds cleanly with 0 errors, 0 warnings.

---

_For current project status, see .planning/PROJECT.md_
