# Roadmap: rbmiUtils

## Milestones

- [x] **v1 Reporting & Robustness** - Phases 1-4 (shipped 2026-02-08)
- [ ] **v2 Documentation & Hardening** - Phases 5-7 (in progress)

## Phases

<details>
<summary>v1 Reporting & Robustness (Phases 1-4) - SHIPPED 2026-02-08</summary>

### Phase 1: Foundation Hardening
**Goal**: Existing foundation functions are robust against edge cases and internal API changes
**Plans**: 3 plans

Plans:
- [x] 01-01: Harden parameter parsing and class detection
- [x] 01-02: Harden gcomp validation and beeca output pinning
- [x] 01-03: Harden storage round-trip integrity

### Phase 2: Print/Summary & ARD Conversion
**Goal**: Users get informative object inspection and pharmaverse-standard ARD output
**Plans**: 2 plans

Plans:
- [x] 02-01: Print/summary S3 methods for pool and analysis objects
- [x] 02-02: pool_to_ard() ARD conversion via cards package

### Phase 3: Efficacy Tables
**Goal**: Users can produce regulatory-style efficacy summary tables from pooled results
**Plans**: 2 plans

Plans:
- [x] 03-01: efficacy_table() gt-based regulatory table
- [x] 03-02: Format helpers and table customization

### Phase 4: Visualization
**Goal**: Users can create publication-quality forest plots from pooled results
**Plans**: 2 plans

Plans:
- [x] 04-01: plot_forest() three-panel patchwork composition
- [x] 04-02: Forest plot modes and accessibility features

</details>

### v2 Documentation & Hardening (In Progress)

**Milestone Goal:** Make rbmiUtils discoverable and trustworthy -- polished site, end-to-end examples, and hardened data prep functions that give clear errors on bad input.

- [x] **Phase 5: Data Prep Hardening** - Validate_data and prepare_data_ice give clear, actionable errors on bad input
- [ ] **Phase 6: Documentation** - End-to-end vignette, rendered examples, and reference docs make the package learnable
- [ ] **Phase 7: Site Polish** - pkgdown site is professional, navigable, and shareable

## Phase Details

### Phase 5: Data Prep Hardening
**Goal**: Users get clear, actionable error messages when data preparation functions receive bad input, and edge cases are handled gracefully
**Depends on**: Phase 4 (v1 complete)
**Requirements**: HRD-01, HRD-02, HRD-03, HRD-04, HRD-05, HRD-06, HRD-07
**Success Criteria** (what must be TRUE):
  1. validate_data() rejects malformed interaction terms (e.g., "A*", ":B", "") and empty data frames with specific, informative error messages
  2. prepare_data_ice() errors immediately when vars$strategy is NULL rather than silently using a default column name
  3. prepare_data_ice() warns when visit column is character with guidance to convert to factor for correct ordering
  4. validate_data() displays all type coercion warnings in a single batched message and warns on all-NA covariate columns
  5. Data prep functions handle edge cases (single subject, single visit, all-NA outcome, all-complete data) without crashing
**Plans**: 2 plans

Plans:
- [x] 05-01: Migrate validate_data() and prepare_data_ice() to cli messaging + add HRD-01 through HRD-06 validation checks
- [x] 05-02: Edge case test coverage (HRD-07) for single subject, single visit, all-NA outcome, all-complete data

### Phase 6: Documentation
**Goal**: Users can learn the full rbmiUtils pipeline from a single end-to-end vignette, see rendered output in function docs, and find version history in NEWS.md
**Depends on**: Phase 5 (hardened functions demonstrated in vignette)
**Requirements**: DOC-01, DOC-02, DOC-03, DOC-04, DOC-05, DOC-06
**Success Criteria** (what must be TRUE):
  1. An end-to-end vignette walks through raw data to rbmi draws/impute/analyse/pool to rbmiUtils tidy/format to efficacy table and forest plot, with inline links to rbmi and beeca documentation
  2. README shows rendered efficacy table and forest plot output as visual teasers and links to the end-to-end vignette
  3. plot_forest() and efficacy_table() help pages display rendered example output
  4. NEWS.md documents v1 and v2 changes with version numbers
  5. Existing vignettes contain inline cross-references to rbmi and beeca documentation where relevant
**Plans**: TBD

Plans:
- [ ] 06-01: TBD
- [ ] 06-02: TBD

### Phase 7: Site Polish
**Goal**: The pkgdown site is professional, well-organized, and produces good social media previews when shared
**Depends on**: Phase 6 (documentation content exists for navbar/reference grouping)
**Requirements**: SITE-01, SITE-02, SITE-03, SITE-04, SITE-05
**Success Criteria** (what must be TRUE):
  1. Package has a hex logo displayed on the pkgdown site header and README
  2. pkgdown navbar includes articles menu, getting started link, news section, and GitHub link
  3. Function reference page groups functions by layer (data prep, analysis, utilities, tidying, formatting, storage, reporting)
  4. Sharing a link to the site on social media/Slack shows a branded open graph card with package name and description
  5. Site footer includes openpharma and pharmaverse links
**Plans**: TBD

Plans:
- [ ] 07-01: TBD

## Progress

**Execution Order:** Phase 5 -> Phase 6 -> Phase 7

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Foundation Hardening | v1 | 3/3 | Complete | 2026-02-07 |
| 2. Print/Summary & ARD | v1 | 2/2 | Complete | 2026-02-07 |
| 3. Efficacy Tables | v1 | 2/2 | Complete | 2026-02-08 |
| 4. Visualization | v1 | 2/2 | Complete | 2026-02-08 |
| 5. Data Prep Hardening | v2 | 2/2 | Complete | 2026-02-08 |
| 6. Documentation | v2 | 0/? | Not started | - |
| 7. Site Polish | v2 | 0/? | Not started | - |

---
*Roadmap created: 2026-02-08*
*Last updated: 2026-02-08 (Phase 5 complete)*
