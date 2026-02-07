# Roadmap: rbmiUtils Reporting & Robustness

## Overview

This roadmap delivers a reporting layer for rbmiUtils that bridges rbmi pooled results into publication-ready regulatory tables and figures via the pharmaverse ARD ecosystem. The work proceeds in four phases: first hardening fragile foundations (parameter parsing, class indexing, input validation), then adding print/summary methods and ARD conversion as infrastructure, then building efficacy tables consuming those ARDs, and finally adding forest plot visualization. Each phase delivers a coherent, testable capability.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3, 4): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 1: Foundation Hardening** - Fix fragile parameter parsing, wrap rbmi::analyse(), harden gcomp and storage functions
- [ ] **Phase 2: Print/Summary & ARD Conversion** - Improve console output for key objects and build ARD interchange layer
- [ ] **Phase 3: Efficacy Tables** - Generate regulatory-style efficacy summary tables via gtsummary + gt
- [ ] **Phase 4: Visualization** - Forest plot function for treatment effects across visits

## Phase Details

### Phase 1: Foundation Hardening
**Goal**: Existing functions are robust against edge cases, version drift, and malformed inputs so new layers build on stable foundations
**Depends on**: Nothing (first phase)
**Requirements**: HRD-01, HRD-02, HRD-03, HRD-04
**Success Criteria** (what must be TRUE):
  1. `tidy_pool_obj()` correctly parses parameter names containing underscores (e.g., "Week_24", "Follow_Up") without data corruption
  2. `analyse_mi_data()` delegates to `rbmi::analyse()` directly, eliminating internal copies of rbmi helper functions
  3. `gcomp_responder()` and `gcomp_binary()` validate inputs and produce clear error messages for edge cases (single-arm data, missing columns, zero-variance outcomes)
  4. `reduce_imputed_data()` and `expand_imputed_data()` preserve column types and attributes through round-trip compression
  5. R CMD check passes with zero errors or warnings after all hardening changes
**Plans**: 3 plans

Plans:
- [ ] 01-01-PLAN.md — Add cli/lifecycle deps and fix tidy_pool_obj parameter parsing (HRD-01)
- [ ] 01-02-PLAN.md — Refactor analyse_mi_data to use inherits() and deprecate internal helpers (HRD-02)
- [ ] 01-03-PLAN.md — Harden gcomp input validation and storage round-trip digest (HRD-03, HRD-04)

### Phase 2: Print/Summary & ARD Conversion
**Goal**: Users get informative console output from key rbmi objects, and tidy pool results convert cleanly to the pharmaverse ARD standard for downstream table generation
**Depends on**: Phase 1 (stable tidy tibble output required for ARD conversion)
**Requirements**: PRT-01, PRT-02, PRT-03, PRT-04, ARD-01, ARD-02, ARD-03
**Success Criteria** (what must be TRUE):
  1. Printing a pool object at the console shows rounded estimates, formatted CIs, and parameter labels without requiring manual formatting
  2. `summary()` on a pool object produces a visit-level breakdown with significance flags
  3. Printing an analysis object shows parameter count, visits covered, and analysis function name
  4. `pool_to_ard()` converts a pool object to a valid cards ARD data frame containing estimate, SE, CI bounds, and p-value per parameter
  5. ARD output preserves visit, parameter type (trt/lsm), and arm as grouping columns, and passes `cards::check_ard_structure()` validation
**Plans**: TBD

Plans:
- [ ] 02-01: TBD
- [ ] 02-02: TBD
- [ ] 02-03: TBD

### Phase 3: Efficacy Tables
**Goal**: Users produce regulatory-style efficacy summary tables directly from rbmi pool objects with a single function call
**Depends on**: Phase 2 (ARD conversion layer required)
**Requirements**: TBL-01, TBL-02, TBL-03, TBL-04
**Success Criteria** (what must be TRUE):
  1. `efficacy_table(pool_obj)` produces a formatted table showing LS means by arm, treatment difference, CIs, and p-values organized by visit
  2. The table renders as gt output suitable for HTML and PDF inclusion in clinical study reports
  3. The table includes footnotes documenting number of imputations, pooling method, and model description
  4. Users can override default formatting (decimal precision, CI bracket style) via function arguments
**Plans**: TBD

Plans:
- [ ] 03-01: TBD
- [ ] 03-02: TBD

### Phase 4: Visualization
**Goal**: Users produce publication-quality forest plots of treatment effects across visits from rbmi pool objects
**Depends on**: Phase 1 (stable tidy tibble output); independent of Phases 2-3
**Requirements**: VIZ-01, VIZ-02, VIZ-03
**Success Criteria** (what must be TRUE):
  1. A forest plot function produces a figure showing treatment effect point estimates and confidence intervals across visits
  2. The function returns a ggplot2 object that users can further customize with standard ggplot2 layers
  3. The plot includes a reference line at zero by default, configurable to any user-specified value
**Plans**: TBD

Plans:
- [ ] 04-01: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4
(Phase 4 is independent of Phases 2-3 and could execute after Phase 1, but is sequenced last for simplicity)

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Foundation Hardening | 0/3 | Planned | - |
| 2. Print/Summary & ARD Conversion | 0/3 | Not started | - |
| 3. Efficacy Tables | 0/2 | Not started | - |
| 4. Visualization | 0/1 | Not started | - |
