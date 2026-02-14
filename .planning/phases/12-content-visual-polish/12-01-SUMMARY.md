---
phase: 12-content-visual-polish
plan: 01
subsystem: docs
tags: [vignette, binary-responder, gcomp, pkgdown, rmarkdown]

requires:
  - phase: none
    provides: "Pre-built ADMI dataset with imputed data"
provides:
  - "Standalone binary responder vignette (deriving-endpoints.Rmd)"
  - "pkgdown navigation entry for new article"
affects: [12-02-content-visual-polish]

tech-stack:
  added: []
  patterns:
    - "Post-imputation binary endpoint derivation from continuous CHG"
    - "gcomp_responder_multi() for multi-visit responder analysis"

key-files:
  created:
    - vignettes/deriving-endpoints.Rmd
  modified:
    - _pkgdown.yml

key-decisions:
  - "Used CHG > 5 as clinical cutoff threshold to demonstrate flexibility of imputed data reuse"
  - "Included cards/ARD section with eval guard for optional dependency"

patterns-established:
  - "Binary responder vignette pattern: load ADMI, factor prep, analyse_mi_data + gcomp_responder_multi, pool, efficacy_table"

duration: 2min
completed: 2026-02-14
---

# Phase 12 Plan 01: Deriving Endpoints Vignette Summary

**Standalone binary responder vignette demonstrating threshold-based (CHG > 3) and clinical cutoff (CHG > 5) endpoint derivation from pre-imputed data using gcomp_responder_multi()**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-14T21:21:22Z
- **Completed:** 2026-02-14T21:23:11Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Created standalone vignette covering two binary responder workflows from imputed continuous data
- Demonstrated imputed data reuse flexibility by deriving a new endpoint (CHG > 5) without re-running draws/impute
- Added ARD conversion section showing pool_to_ard() for pharmaverse integration
- Included caveats on imputation model assumptions, threshold pre-specification, and reference-based conditioning

## Task Commits

Each task was committed atomically:

1. **Task 1: Create deriving-endpoints vignette** - `6c6aa58` (feat)
2. **Task 2: Add vignette to pkgdown navigation** - `771113b` (chore)

## Files Created/Modified

- `vignettes/deriving-endpoints.Rmd` - Standalone binary responder vignette (273 lines)
- `_pkgdown.yml` - Added article entry under Workflow Guides

## Decisions Made

- Used CHG > 5 (not a clinical scale cutoff) as the second threshold -- keeps the example simple and clearly distinct from the pre-derived CHG > 3
- Guarded efficacy_table() chunks with `eval = requireNamespace("gt", quietly = TRUE)` and ARD chunk with cards guard, matching pipeline.Rmd pattern

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

- Pre-existing test failures in test-analysis_utils.R due to namespace loading when using testthat::test_dir() directly (gcomp_responder / gcomp_responder_multi not found). These are not caused by this plan's changes and exist on the prior commit.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Vignette builds cleanly and is listed in pkgdown navigation
- Ready for Phase 12 Plan 02 (content/visual polish)

## Self-Check: PASSED

- [x] vignettes/deriving-endpoints.Rmd exists
- [x] _pkgdown.yml exists
- [x] Commit 6c6aa58 found
- [x] Commit 771113b found

---
*Phase: 12-content-visual-polish*
*Completed: 2026-02-14*
