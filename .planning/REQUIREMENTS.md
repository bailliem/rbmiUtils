# Requirements: rbmiUtils Reporting & Robustness

**Defined:** 2026-02-07
**Core Value:** Clinical trial results from rbmi flow seamlessly into publication-ready regulatory tables and figures

## v1 Requirements

### ARD & Reporting

- [ ] **ARD-01**: Convert tidy pool results to cards ARD format via `pool_to_ard()`
- [ ] **ARD-02**: ARD includes all statistics (estimate, SE, lower CI, upper CI, p-value) per parameter
- [ ] **ARD-03**: ARD preserves visit, parameter type (trt/lsm), and arm metadata as grouping columns

### Efficacy Tables

- [ ] **TBL-01**: Generate regulatory-style efficacy summary table (LS means by arm, treatment difference, CIs, p-values by visit)
- [ ] **TBL-02**: Table renders via gtsummary + gt with HTML/PDF output
- [ ] **TBL-03**: One-call `efficacy_table(pool_obj)` from pool object directly to gt table with opinionated defaults
- [ ] **TBL-04**: Table includes footnotes (number of imputations, pooling method, model description)

### Visualization

- [ ] **VIZ-01**: Forest plot function showing treatment effects across visits with point estimates and CIs
- [ ] **VIZ-02**: Forest plot returns ggplot2 object for user customization
- [ ] **VIZ-03**: Forest plot supports reference line at zero (or user-specified value)

### Print/Summary Methods

- [ ] **PRT-01**: Enhanced print method for pool objects showing rounded estimates, formatted CIs, and parameter labels
- [ ] **PRT-02**: Summary method for pool objects with additional detail (visit breakdown, significance flags)
- [ ] **PRT-03**: Improved print.analysis() showing parameter count, visit info, function name
- [ ] **PRT-04**: Improved summary.analysis() with parameter preview table

### Foundation Hardening

- [ ] **HRD-01**: Fix tidy_pool_obj() parameter parsing — replace `_` separator with robust regex
- [ ] **HRD-02**: Refactor analyse_mi_data() to wrap rbmi::analyse() instead of reimplementing internals
- [ ] **HRD-03**: Harden gcomp functions — validate inputs, handle edge cases, pin beeca output format
- [ ] **HRD-04**: Harden storage functions — type coercion, attribute preservation, string key collision handling

## v2 Requirements

### Visualization

- **VIZ-04**: Responder analysis bar chart (proportion responding by arm and visit)
- **VIZ-05**: Forest plot with sensitivity analysis overlay (primary + delta-adjusted side-by-side)
- **VIZ-06**: Responder chart with treatment difference annotations on plot

### Tables

- **TBL-05**: Column formatting controls for gt/gtsummary theming (decimal precision, CI style)
- **TBL-06**: Sensitivity analysis comparison table (primary vs. delta-adjusted side-by-side)

### ARD Extensions

- **ARD-04**: MI-specific metadata in ARD (number of imputations, pooling method, FMI)
- **ARD-05**: as_gt() / as_gtsummary() S3 methods for pool objects

### Print/Summary Extensions

- **PRT-05**: describe_draws() helper for draws objects
- **PRT-06**: describe_imputation() helper for impute objects

## Out of Scope

| Feature | Reason |
|---------|--------|
| rtables/tern integration | Different ecosystem (cell-based), incompatible with ARD paradigm. Users can extract tidy tibble for rtables. |
| Interactive Shiny dashboards | teal/teal.modules.clinical covers this space. Static figures are the deliverable. |
| Word/PowerPoint direct export | gt handles HTML/PDF. gt::as_rtf() available for Word needs. No flextable dependency. |
| Safety tables (AE, labs, vitals) | Outside rbmi efficacy domain. Handled by tern/gtreg. |
| Custom pooling methods | Use rbmi's built-in pooling. Reimplementing pool() is out of scope. |
| Generic table layout engine | gt/gtsummary are mature table frameworks. rbmiUtils provides content, not layout. |
| Spaghetti/trajectory plots | Different visualization category from pooled results. get_imputed_data() available for users. |
| Reimplementing cardx functions | Use cardx where applicable. Only build custom ARD for rbmi-specific output. |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| ARD-01 | — | Pending |
| ARD-02 | — | Pending |
| ARD-03 | — | Pending |
| TBL-01 | — | Pending |
| TBL-02 | — | Pending |
| TBL-03 | — | Pending |
| TBL-04 | — | Pending |
| VIZ-01 | — | Pending |
| VIZ-02 | — | Pending |
| VIZ-03 | — | Pending |
| PRT-01 | — | Pending |
| PRT-02 | — | Pending |
| PRT-03 | — | Pending |
| PRT-04 | — | Pending |
| HRD-01 | — | Pending |
| HRD-02 | — | Pending |
| HRD-03 | — | Pending |
| HRD-04 | — | Pending |

**Coverage:**
- v1 requirements: 18 total
- Mapped to phases: 0
- Unmapped: 18

---
*Requirements defined: 2026-02-07*
*Last updated: 2026-02-07 after initial definition*
