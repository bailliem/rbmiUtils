# Roadmap: rbmiUtils

## Milestones

- ✅ **v1 Reporting & Robustness** -- Phases 1-4 (shipped 2026-02-08)
- ✅ **v2 Documentation & Hardening** -- Phases 5-7 (shipped 2026-02-10)
- **v3 ARD Enrichment & Polish** -- Phases 8-11 (in progress)

## Phases

<details>
<summary>v1 Reporting & Robustness (Phases 1-4) -- SHIPPED 2026-02-08</summary>

- [x] Phase 1: Foundation Hardening (3/3 plans) -- completed 2026-02-07
- [x] Phase 2: Print/Summary & ARD (2/2 plans) -- completed 2026-02-07
- [x] Phase 3: Efficacy Tables (2/2 plans) -- completed 2026-02-08
- [x] Phase 4: Visualization (2/2 plans) -- completed 2026-02-08

</details>

<details>
<summary>v2 Documentation & Hardening (Phases 5-7) -- SHIPPED 2026-02-10</summary>

- [x] Phase 5: Data Prep Hardening (2/2 plans) -- completed 2026-02-08
- [x] Phase 6: Documentation (3/3 plans) -- completed 2026-02-08
- [x] Phase 7: Site Polish (2/2 plans) -- completed 2026-02-08

</details>

### v3 ARD Enrichment & Polish (In Progress)

**Milestone Goal:** Enrich ARD output with MI-specific diagnostic metadata, add imputation diagnostic helpers, and polish tables, plots, and documentation to publication quality.

- [x] **Phase 8: MI Diagnostic Statistics** -- Enriched ARD with FMI, lambda, variance decomposition, and Rubin's rules metadata -- completed 2026-02-10
- [x] **Phase 9: Describe Helpers** -- describe_draws() and describe_imputation() with cli-formatted print methods and 98 tests -- completed 2026-02-11
- [ ] **Phase 10: Publication Styling** -- Font, spacing, and layout controls for tables and forest plots
- [ ] **Phase 11: Documentation Overhaul** -- README, examples, vignettes, and images reflecting all v3 features

## Phase Details

### Phase 8: MI Diagnostic Statistics
**Goal**: Users can access MI-specific diagnostic metadata (FMI, lambda, degrees of freedom, relative efficiency) directly from ARD output, enabling regulatory reviewers to assess imputation quality without manual recomputation
**Depends on**: Phase 7 (v2 complete -- existing pool_to_ard() is the integration point)
**Requirements**: MIDIAG-01, MIDIAG-02, MIDIAG-03, MIDIAG-05, MIDIAG-06, MIDIAG-07, MIDIAG-08 (MIDIAG-04 removed: V_w/V_b/V_t excluded per user decision)
**Success Criteria** (what must be TRUE):
  1. User can call pool_to_ard(pool_obj, analysis_obj) and the resulting ARD contains FMI, lambda, RIV, Barnard-Rubin adjusted df, and relative efficiency as stat_name rows for each parameter
  2. The enriched ARD passes cards::check_ard_structure() validation without errors
  3. User calling pool_to_ard(pool_obj) without analysis_obj gets the same base ARD as before (backward compatible)
  4. When pooling method is not Rubin's rules, MI diagnostic rows are omitted entirely with an informative cli message
**Plans**: 2 plans

Plans:
- [x] 08-01-PLAN.md -- TDD: Rubin's rules diagnostic computation (compute_rubin_diagnostics)
- [x] 08-02-PLAN.md -- Integrate MI diagnostics into pool_to_ard() with tests and docs

### Phase 9: Describe Helpers
**Goal**: Users can inspect draws and imputation objects to understand what happened during the MI pipeline -- method used, sample counts, convergence, missingness patterns -- without reading raw object internals
**Depends on**: Nothing (independent of Phase 8; both depend on v2 baseline)
**Requirements**: DESC-01, DESC-02, DESC-03, DESC-04
**Success Criteria** (what must be TRUE):
  1. User can call describe_draws(draws_obj) and see a formatted summary showing method, number of samples, number of failures, and the model formula
  2. User can see MCMC convergence diagnostics (ESS, Rhat) from describe_draws() when the draws used Bayesian methods
  3. User can call describe_imputation(imputed_data, original_data) and see method, number of imputations (M), and missingness breakdown by visit and treatment arm
  4. Both describe functions return structured S3 objects with informative print() output using cli formatting
**Plans**: 2 plans

Plans:
- [x] 09-01-PLAN.md -- TDD: describe_draws() with MCMC diagnostics and print method
- [x] 09-02-PLAN.md -- TDD: describe_imputation() with missingness breakdown and print method

### Phase 10: Publication Styling
**Goal**: Users can produce publication-quality tables and forest plots with controlled typography, spacing, and layout without post-hoc manual adjustments
**Depends on**: Nothing (independent of Phases 8-9; pure visual refinement of existing functions)
**Requirements**: STYLE-01, STYLE-02, STYLE-03, STYLE-04, STYLE-05
**Success Criteria** (what must be TRUE):
  1. User can specify font_family and font_size parameters in efficacy_table() and get a gt table rendered in that font with controlled row padding
  2. User can specify font_family parameter in plot_forest() and get consistent typography across all three panels
  3. User can control panel width ratios in plot_forest() via a panel_widths parameter to adjust relative sizes of table, forest, and p-value panels
  4. Visit labels and estimate text in the forest plot left panel align consistently regardless of label length or font choice
**Plans**: 2 plans

Plans:
- [ ] 10-01-PLAN.md -- Add font_family, font_size, row_padding parameters to efficacy_table()
- [ ] 10-02-PLAN.md -- Add font_family, panel_widths to plot_forest() and fix left-panel alignment

### Phase 11: Documentation Overhaul
**Goal**: All documentation reflects the finalized v3 features with realistic clinical trial examples, updated vignettes, regenerated images, and a versioned changelog
**Depends on**: Phases 8, 9, 10 (documents finalized features; regenerates images after styling changes)
**Requirements**: DOCS-01, DOCS-02, DOCS-03, DOCS-04, DOCS-05
**Success Criteria** (what must be TRUE):
  1. README demonstrates a realistic clinical trial workflow from ADEFF data through the full rbmi pipeline to efficacy_table() and plot_forest() output
  2. Exported function documentation examples use ADMI/ADEFF sample datasets showing real clinical analysis patterns (not minimal synthetic data)
  3. Vignettes cover MI diagnostics from pool_to_ard() and the describe_draws()/describe_imputation() helpers with worked examples
  4. Pre-rendered images for README and help pages reflect the v3 styling improvements
  5. NEWS.md contains v0.3.0 entries documenting all v3 additions
**Plans**: TBD

Plans:
- [ ] 11-01: TBD
- [ ] 11-02: TBD

## Progress

**Execution Order:**
Phases 8-11 execute sequentially. Phases 8, 9, and 10 have no mutual dependencies but Phase 11 depends on all three being complete. Recommended order: 8 -> 9 -> 10 -> 11.

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Foundation Hardening | v1 | 3/3 | Complete | 2026-02-07 |
| 2. Print/Summary & ARD | v1 | 2/2 | Complete | 2026-02-07 |
| 3. Efficacy Tables | v1 | 2/2 | Complete | 2026-02-08 |
| 4. Visualization | v1 | 2/2 | Complete | 2026-02-08 |
| 5. Data Prep Hardening | v2 | 2/2 | Complete | 2026-02-08 |
| 6. Documentation | v2 | 3/3 | Complete | 2026-02-08 |
| 7. Site Polish | v2 | 2/2 | Complete | 2026-02-08 |
| 8. MI Diagnostic Statistics | v3 | 2/2 | Complete | 2026-02-10 |
| 9. Describe Helpers | v3 | 2/2 | Complete | 2026-02-11 |
| 10. Publication Styling | v3 | 0/2 | Not started | - |
| 11. Documentation Overhaul | v3 | 0/TBD | Not started | - |

---
*Roadmap created: 2026-02-08*
*Last updated: 2026-02-11 (Phase 10 planned: 2 plans in 1 wave)*
