# Requirements: rbmiUtils v3 ARD Enrichment & Polish

**Defined:** 2026-02-10
**Core Value:** Clinical trial results from rbmi flow seamlessly into publication-ready regulatory tables and figures

## v3 Requirements

Requirements for v3 milestone. Each maps to roadmap phases.

### MI Diagnostics

- [ ] **MIDIAG-01**: User can obtain FMI (fraction of missing information) per parameter in ARD output from pool_to_ard()
- [ ] **MIDIAG-02**: User can obtain relative increase in variance (RIV) per parameter in ARD output
- [ ] **MIDIAG-03**: User can obtain lambda (proportion of variance due to missingness) per parameter in ARD output
- [ ] **MIDIAG-04**: ~~User can obtain within/between/total variance components per parameter in ARD output~~ REMOVED: User decided to exclude V_w/V_b/V_t from ARD output (curated essentials only). Variance components are computed internally but not exposed.
- [ ] **MIDIAG-05**: User can obtain Barnard-Rubin adjusted degrees of freedom per parameter in ARD output
- [ ] **MIDIAG-06**: User can obtain relative efficiency (RE) per parameter in ARD output
- [ ] **MIDIAG-07**: Pooling method is standardized as a proper stat_name row in ARD output
- [ ] **MIDIAG-08**: MI diagnostic rows are omitted entirely (not returned as NA rows) with informative cli message for non-Rubin pooling methods

### Describe Helpers

- [ ] **DESC-01**: User can call describe_draws() on a draws object to get method, sample count, failures, formula summary
- [ ] **DESC-02**: User can see MCMC convergence diagnostics (ESS, Rhat) from describe_draws() for Bayesian draws
- [ ] **DESC-03**: User can call describe_imputation() on imputed + original data to get method, M, missingness by visit/arm
- [ ] **DESC-04**: Both describe functions return structured objects with informative print methods

### Publication Styling

- [ ] **STYLE-01**: User can specify font family for efficacy_table() output
- [ ] **STYLE-02**: User can control font size and row padding for compact regulatory table layout
- [ ] **STYLE-03**: User can specify font_family parameter for plot_forest()
- [ ] **STYLE-04**: User can control panel width ratios in plot_forest() via panel_widths parameter
- [ ] **STYLE-05**: Forest plot left panel alignment is improved for consistent visit label/estimate text positioning

### Documentation

- [ ] **DOCS-01**: README shows realistic clinical trial workflow from ADEFF through rbmi pipeline to table and forest plot
- [ ] **DOCS-02**: Function documentation examples use ADMI/ADEFF sample data showing real usage patterns
- [ ] **DOCS-03**: Vignettes updated to cover MI diagnostics and describe helpers
- [ ] **DOCS-04**: Pre-rendered images regenerated reflecting styling improvements
- [ ] **DOCS-05**: NEWS.md updated with v0.3.0 entries

## Future Requirements

Deferred to v4+. Tracked but not in current roadmap.

### Visualization Extensions

- **VIZ-01**: Responder bar chart function (proportion responding by arm and visit)
- **VIZ-02**: Forest plot with sensitivity analysis overlay
- **VIZ-03**: Responder chart with treatment difference annotations

### Advanced Reporting

- **RPT-01**: Sensitivity analysis comparison table
- **RPT-02**: Column formatting controls for gt theming
- **RPT-03**: as_gt() / as_gtsummary() S3 methods for pool objects

### Advanced Diagnostics

- **DIAG-01**: BMLMI-specific MI diagnostics (modified Rubin's rules)
- **DIAG-02**: Imputation diagnostic plots (trace, density, stripplot)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Re-implementing rbmi::pool() to store FMI | rbmiUtils is a utility layer, not a fork; recomputation from analysis object is cleaner |
| Generic imputation diagnostic plots | Specialized visualization; users can use bayesplot/rstan directly |
| Automatic quality thresholds for MI diagnostics | Context-dependent; let the user interpret |
| gtsummary tbl_ard_*() integration | ARD format doesn't match gtsummary's assumptions for regulatory tables |
| MI diagnostics for non-Rubin methods | FMI/lambda/RIV are Rubin-specific; bootstrap/jackknife have different decompositions |
| Sensitivity analysis features | Deferred to v4+; independent of v3 MI metadata work |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| MIDIAG-01 | Phase 8 | Pending |
| MIDIAG-02 | Phase 8 | Pending |
| MIDIAG-03 | Phase 8 | Pending |
| MIDIAG-04 | Phase 8 | Pending |
| MIDIAG-05 | Phase 8 | Pending |
| MIDIAG-06 | Phase 8 | Pending |
| MIDIAG-07 | Phase 8 | Pending |
| MIDIAG-08 | Phase 8 | Pending |
| DESC-01 | Phase 9 | Pending |
| DESC-02 | Phase 9 | Pending |
| DESC-03 | Phase 9 | Pending |
| DESC-04 | Phase 9 | Pending |
| STYLE-01 | Phase 10 | Pending |
| STYLE-02 | Phase 10 | Pending |
| STYLE-03 | Phase 10 | Pending |
| STYLE-04 | Phase 10 | Pending |
| STYLE-05 | Phase 10 | Pending |
| DOCS-01 | Phase 11 | Pending |
| DOCS-02 | Phase 11 | Pending |
| DOCS-03 | Phase 11 | Pending |
| DOCS-04 | Phase 11 | Pending |
| DOCS-05 | Phase 11 | Pending |

**Coverage:**
- v3 requirements: 22 total
- Mapped to phases: 22
- Unmapped: 0

---
*Requirements defined: 2026-02-10*
*Last updated: 2026-02-10 (traceability updated with phase mappings)*
