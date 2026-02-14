# Roadmap: rbmiUtils

## Milestones

- [x] **v1 Reporting & Robustness** - Phases 1-4 (shipped 2026-02-08)
- [x] **v2 Documentation & Hardening** - Phases 5-7 (shipped 2026-02-10)
- [x] **v3 ARD Enrichment & Polish** - Phases 8-11 (shipped 2026-02-11)
- [ ] **v4 CRAN Release Readiness** - Phases 12-14 (in progress)

## Phases

<details>
<summary>v1 Reporting & Robustness (Phases 1-4) - SHIPPED 2026-02-08</summary>

Phases 1-4: Foundation hardening, print/summary, ARD, efficacy tables, forest plots.
See .planning/MILESTONES.md for details.

</details>

<details>
<summary>v2 Documentation & Hardening (Phases 5-7) - SHIPPED 2026-02-10</summary>

Phases 5-7: Data prep hardening, end-to-end vignette, pkgdown site polish.
See .planning/MILESTONES.md for details.

</details>

<details>
<summary>v3 ARD Enrichment & Polish (Phases 8-11) - SHIPPED 2026-02-11</summary>

Phases 8-11: MI diagnostic metadata, describe helpers, publication styling, documentation overhaul.
See .planning/MILESTONES.md for details.

</details>

### v4 CRAN Release Readiness (In Progress)

**Milestone Goal:** Polish, harden, and prepare the package for initial CRAN submission -- zero new features, pure quality and compliance.

- [x] **Phase 12: Content & Visual Polish** - Binary responder vignette and forest plot refinements (completed 2026-02-14)
- [ ] **Phase 13: CRAN Compliance** - Policy fixes, warning suppression, timing audit, NEWS cleanup
- [ ] **Phase 14: Final Validation** - R CMD check --as-cran passes clean

## Phase Details

### Phase 12: Content & Visual Polish
**Goal**: Users have a standalone binary responder vignette demonstrating the imputed data storage workflow, and forest plots meet publication-quality visual standards
**Depends on**: Phase 11 (v3 complete)
**Requirements**: DOC-01, VIZ-01
**Success Criteria** (what must be TRUE):
  1. A standalone vignette for binary responder analysis exists, separate from the pipeline vignette, demonstrating imputed data storage, modification, and reanalysis
  2. The binary responder vignette builds without errors and appears in the pkgdown articles listing
  3. Forest plot output shows refined typography, spacing, and styling suitable for regulatory submissions
  4. Existing forest plot tests still pass after visual refinements
**Plans:** 2 plans

Plans:
- [x] 12-01-PLAN.md -- Standalone binary responder vignette (deriving-endpoints.Rmd) and pkgdown nav update
- [x] 12-02-PLAN.md -- Forest plot visual polish (typography, spacing, dashed ref line, panel border)

### Phase 13: CRAN Compliance
**Goal**: All CRAN policy requirements are satisfied -- DESCRIPTION metadata, vignette build cleanliness, timing limits, and NEWS formatting
**Depends on**: Phase 12 (content must exist before compliance audit)
**Requirements**: CRAN-01, CRAN-02, CRAN-03, DOC-02
**Success Criteria** (what must be TRUE):
  1. DESCRIPTION Title, Description, Authors@R, License, URL, and BugReports fields pass CRAN policy checks without notes
  2. All vignettes build without Stan compilation warnings or other extraneous console output
  3. All tests and examples complete within CRAN time limits (no single test file exceeds 60 seconds, total check under 10 minutes)
  4. NEWS.md has consistent formatting across all version entries (0.1.0, 0.2.0, 0.3.0, 0.4.0) suitable for CRAN
**Plans**: TBD

Plans:
- [ ] 13-01: DESCRIPTION audit and NEWS cleanup
- [ ] 13-02: Vignette warning suppression and timing audit

### Phase 14: Final Validation
**Goal**: R CMD check --as-cran passes with zero errors, zero warnings, and zero notes -- the package is ready for CRAN submission
**Depends on**: Phase 13 (all compliance fixes applied)
**Requirements**: CRAN-04
**Success Criteria** (what must be TRUE):
  1. `R CMD check --as-cran` returns 0 errors, 0 warnings, 0 notes on the development machine
  2. All vignettes, tests, and examples execute successfully during check
  3. Package tarball builds cleanly with no unexpected files or oversized components
**Plans**: TBD

Plans:
- [ ] 14-01: R CMD check clean pass

## Progress

**Execution Order:**
Phases execute in numeric order: 12 -> 13 -> 14

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 12. Content & Visual Polish | v4 | 2/2 | ✓ Complete | 2026-02-14 |
| 13. CRAN Compliance | v4 | 0/3 | Not started | - |
| 14. Final Validation | v4 | 0/1 | Not started | - |

---
*Roadmap created: 2026-02-14*
*Last updated: 2026-02-14 (Phase 12 complete)*
