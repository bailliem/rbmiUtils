# Requirements: rbmiUtils v2 Documentation & Hardening

**Defined:** 2026-02-08
**Core Value:** Clinical trial results from rbmi flow seamlessly into publication-ready regulatory tables and figures — no manual wrangling between pooled results and final output.

## v2 Requirements

### Hardening (HRD)

- [x] **HRD-01**: validate_data() rejects malformed interaction terms (e.g., "A*", ":B", "") with clear error message
- [x] **HRD-02**: prepare_data_ice() errors when vars$strategy is NULL instead of silently defaulting to "strategy"
- [x] **HRD-03**: prepare_data_ice() warns when visit column is character (not factor) with ordering guidance
- [x] **HRD-04**: validate_data() rejects empty (zero-row) data frames with informative minimum data requirements error
- [x] **HRD-05**: validate_data() warns on all-NA covariate columns with actionable guidance
- [x] **HRD-06**: validate_data() batch-displays type coercion warnings instead of emitting them one at a time
- [x] **HRD-07**: Edge case test coverage for single subject, single visit, all-NA outcome, all-complete data across data prep functions

### Documentation (DOC)

- [ ] **DOC-01**: End-to-end clinical trial analysis vignette covering raw data → rbmi draws/impute/analyse/pool → rbmiUtils tidy/format → efficacy table + forest plot with inline rbmi/beeca links
- [ ] **DOC-02**: README enhanced with visual teaser showing rendered efficacy table and forest plot output, linking to end-to-end vignette
- [ ] **DOC-03**: Rendered @examples for plot_forest() showing treatment effect forest plot
- [ ] **DOC-04**: Rendered @examples for efficacy_table() showing regulatory summary table
- [ ] **DOC-05**: NEWS.md tracking v1 and v2 changes with version history
- [ ] **DOC-06**: Inline cross-references to rbmi and beeca documentation in all existing vignettes where relevant

### Site (SITE)

- [ ] **SITE-01**: Hex logo for rbmiUtils package
- [ ] **SITE-02**: Custom pkgdown navbar with articles menu, getting started link, and news
- [ ] **SITE-03**: Grouped function reference organized by layer (data prep, analysis, utilities, tidying, formatting, storage, reporting)
- [ ] **SITE-04**: Social/open graph cards for link sharing
- [ ] **SITE-05**: Custom footer with openpharma/pharmaverse links

## Future Requirements (v3+)

Deferred from v1 Active list. Not in current roadmap.

- Responder bar chart function (proportion responding by arm and visit)
- Forest plot with sensitivity analysis overlay
- Responder chart with treatment difference annotations
- Column formatting controls for gt theming
- Sensitivity analysis comparison table
- MI-specific metadata in ARD (FMI, pooling method)
- as_gt() / as_gtsummary() S3 methods for pool objects
- describe_draws() helper for draws objects
- describe_imputation() helper for impute objects

## Out of Scope

| Feature | Reason |
|---------|--------|
| Interactive/Shiny dashboards | Focus is on static report generation |
| Word/PowerPoint output via flextable | gt + HTML/PDF is sufficient |
| Custom pooling methods | Use rbmi's built-in pooling |
| Spaghetti/trajectory plots | Not requested |
| Safety tables (AE listings) | Outside rbmi efficacy workflow |
| rtables/tern integration | Different ecosystem, incompatible with ARD paradigm |
| Ecosystem overview page | Inline cross-references in vignettes sufficient for v2 |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| HRD-01 | Phase 5 | Complete |
| HRD-02 | Phase 5 | Complete |
| HRD-03 | Phase 5 | Complete |
| HRD-04 | Phase 5 | Complete |
| HRD-05 | Phase 5 | Complete |
| HRD-06 | Phase 5 | Complete |
| HRD-07 | Phase 5 | Complete |
| DOC-01 | Phase 6 | Pending |
| DOC-02 | Phase 6 | Pending |
| DOC-03 | Phase 6 | Pending |
| DOC-04 | Phase 6 | Pending |
| DOC-05 | Phase 6 | Pending |
| DOC-06 | Phase 6 | Pending |
| SITE-01 | Phase 7 | Pending |
| SITE-02 | Phase 7 | Pending |
| SITE-03 | Phase 7 | Pending |
| SITE-04 | Phase 7 | Pending |
| SITE-05 | Phase 7 | Pending |

**Coverage:**
- v2 requirements: 18 total
- Mapped to phases: 18
- Unmapped: 0

---
*Requirements defined: 2026-02-08*
*Last updated: 2026-02-08 (HRD-01 through HRD-07 complete)*
