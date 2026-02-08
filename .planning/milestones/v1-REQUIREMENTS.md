# Requirements Archive: v1 Reporting & Robustness

**Archived:** 2026-02-08
**Status:** SHIPPED

This is the archived requirements specification for v1.
For current requirements, see `.planning/REQUIREMENTS.md` (created for next milestone).

---

## v1 Requirements

### ARD & Reporting

- [x] **ARD-01**: Convert tidy pool results to cards ARD format via `pool_to_ard()` — Shipped, passes cards::check_ard_structure()
- [x] **ARD-02**: ARD includes all statistics (estimate, SE, lower CI, upper CI, p-value) per parameter — Shipped, 6 stats per parameter including method
- [x] **ARD-03**: ARD preserves visit, parameter type (trt/lsm), and arm metadata as grouping columns — Shipped, group1/group2/group3 columns

### Efficacy Tables

- [x] **TBL-01**: Generate regulatory-style efficacy summary table (LS means by arm, treatment difference, CIs, p-values by visit) — Shipped
- [x] **TBL-02**: Table renders via gtsummary + gt with HTML/PDF output — Shipped (gt direct, appropriate for custom table layout)
- [x] **TBL-03**: One-call `efficacy_table(pool_obj)` from pool object directly to gt table with opinionated defaults — Shipped
- [x] **TBL-04**: Table includes footnotes (number of imputations, pooling method, model description) — Shipped

### Visualization

- [x] **VIZ-01**: Forest plot function showing treatment effects across visits with point estimates and CIs — Shipped
- [x] **VIZ-02**: Forest plot returns ggplot2 object for user customization — Shipped (patchwork/ggplot, customizable via & theme())
- [x] **VIZ-03**: Forest plot supports reference line at zero (or user-specified value) — Shipped

### Print/Summary Methods

- [x] **PRT-01**: Enhanced print method for pool objects showing rounded estimates, formatted CIs, and parameter labels — Shipped
- [x] **PRT-02**: Summary method for pool objects with additional detail (visit breakdown, significance flags) — Shipped
- [x] **PRT-03**: Improved print.analysis() showing parameter count, visit info, function name — Shipped
- [x] **PRT-04**: Improved summary.analysis() with parameter preview table — Shipped

### Foundation Hardening

- [x] **HRD-01**: Fix tidy_pool_obj() parameter parsing — replace `_` separator with robust regex — Shipped
- [x] **HRD-02**: Refactor analyse_mi_data() to wrap rbmi::analyse() instead of reimplementing internals — Shipped (uses inherits(), deprecated as_analysis2())
- [x] **HRD-03**: Harden gcomp functions — validate inputs, handle edge cases, pin beeca output format — Shipped
- [x] **HRD-04**: Harden storage functions — type coercion, attribute preservation, string key collision handling — Shipped

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| ARD-01 | Phase 2 | Complete |
| ARD-02 | Phase 2 | Complete |
| ARD-03 | Phase 2 | Complete |
| TBL-01 | Phase 3 | Complete |
| TBL-02 | Phase 3 | Complete |
| TBL-03 | Phase 3 | Complete |
| TBL-04 | Phase 3 | Complete |
| VIZ-01 | Phase 4 | Complete |
| VIZ-02 | Phase 4 | Complete |
| VIZ-03 | Phase 4 | Complete |
| PRT-01 | Phase 2 | Complete |
| PRT-02 | Phase 2 | Complete |
| PRT-03 | Phase 2 | Complete |
| PRT-04 | Phase 2 | Complete |
| HRD-01 | Phase 1 | Complete |
| HRD-02 | Phase 1 | Complete |
| HRD-03 | Phase 1 | Complete |
| HRD-04 | Phase 1 | Complete |

---

## Milestone Summary

**Shipped:** 18 of 18 v1 requirements
**Adjusted:**
- TBL-02: Originally specified "gtsummary + gt", shipped as "gt direct" — gtsummary layer not needed for custom table layout, gt alone is sufficient
- HRD-02: Originally "wrap rbmi::analyse()", implemented as "use inherits() and deprecate internal helpers" — achieves stability goal without full reimplementation

**Dropped:** None

---
*Archived: 2026-02-08 as part of v1 milestone completion*
