# Requirements: rbmiUtils

**Defined:** 2026-02-14
**Core Value:** Clinical trial results from rbmi flow seamlessly into publication-ready regulatory tables and figures

## v4 Requirements

Requirements for CRAN release readiness. Each maps to roadmap phases.

### Documentation

- [ ] **DOC-01**: Binary responder analysis exists as standalone vignette showcasing imputed data storage, modification, and reanalysis workflow
- [ ] **DOC-02**: NEWS.md older version entries cleaned up to CRAN standards with consistent formatting

### Visual Polish

- [ ] **VIZ-01**: Forest plot has refined typography, spacing, and styling for publication quality

### CRAN Compliance

- [ ] **CRAN-01**: DESCRIPTION passes CRAN policy checks (Title, Description, Authors, License, URL formatting)
- [ ] **CRAN-02**: All vignettes build without Stan compilation warnings or other noise
- [ ] **CRAN-03**: Tests and examples complete within CRAN time limits (no timeouts)
- [ ] **CRAN-04**: R CMD check --as-cran passes with 0 errors, 0 warnings, 0 notes

## v5+ Requirements

Deferred to future release. Not in current roadmap.

### Visualization

- Responder bar chart function (proportion responding by arm and visit)
- Forest plot with sensitivity analysis overlay
- Responder chart with treatment difference annotations

### Theming

- Column formatting controls for gt theming

### Analysis

- Sensitivity analysis comparison table
- as_gt() / as_gtsummary() S3 methods for pool objects

### Diagnostics

- BMLMI-specific MI diagnostics (modified Rubin's rules)
- Imputation diagnostic plots (trace, density, stripplot)

## Out of Scope

| Feature | Reason |
|---------|--------|
| New exported functions | v4 is pure polish and compliance |
| Interactive/Shiny dashboards | Focus on static report generation |
| Word/PowerPoint output | gt + HTML/PDF is sufficient |
| rtables/tern integration | Different ecosystem |
| New features | CRAN readiness only |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| DOC-01 | — | Pending |
| DOC-02 | — | Pending |
| VIZ-01 | — | Pending |
| CRAN-01 | — | Pending |
| CRAN-02 | — | Pending |
| CRAN-03 | — | Pending |
| CRAN-04 | — | Pending |

**Coverage:**
- v4 requirements: 7 total
- Mapped to phases: 0
- Unmapped: 7

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-14 after initial definition*
