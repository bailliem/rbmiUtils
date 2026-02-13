---
phase: quick-1
plan: 01
subsystem: analysis
tags: [vignette, stan, deprecation, gcomp, warnings]

# Dependency graph
requires: []
provides:
  - "Warning-free vignette build for pipeline.Rmd"
  - "Deprecation-free internal calls in gcomp_responder()"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Local helper functions to avoid calling deprecated exports internally"

key-files:
  created: []
  modified:
    - "vignettes/pipeline.Rmd"
    - "R/analysis_utils.R"

key-decisions:
  - "Used local .extract_covars() helper inside gcomp_responder rather than removing deprecation from exported functions"
  - "Kept deprecated functions exported unchanged for backwards compatibility"

patterns-established:
  - "Internal code should not call its own deprecated exports; inline the logic instead"

# Metrics
duration: 2min
completed: 2026-02-13
---

# Quick Task 1: Fix Vignette Warnings Summary

**Suppressed Stan compilation warnings in pipeline vignette and inlined deprecated helper logic in gcomp_responder to eliminate deprecation warnings**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-13T21:55:26Z
- **Completed:** 2026-02-13T21:57:12Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Stan compilation noise no longer leaks through as warnings in the pipeline vignette draws chunk
- gcomp_responder() no longer triggers deprecation warnings from extract_covariates2() and as_simple_formula2()
- All 1006 existing tests pass with zero failures
- Deprecated functions remain exported for external callers' backwards compatibility

## Task Commits

Each task was committed atomically:

1. **Task 1: Suppress Stan compilation warnings in vignette** - `6690aba` (fix)
2. **Task 2: Inline deprecated helper logic in gcomp_responder** - `5785499` (fix)

## Files Created/Modified
- `vignettes/pipeline.Rmd` - Added `warning = FALSE` to the draws chunk header
- `R/analysis_utils.R` - Added local `.extract_covars()` helper and inlined formula construction in `gcomp_responder()`

## Decisions Made
- Used a local `.extract_covars()` helper defined inside `gcomp_responder()` to keep the inlined logic clean and avoid duplicating the strsplit/trimws pattern inline at two call sites
- Kept deprecated `extract_covariates2()` and `as_simple_formula2()` functions completely unchanged for backwards compatibility

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Vignette and R CMD check should now be free of Stan compilation and deprecation warnings
- No blockers or concerns

## Self-Check: PASSED

- FOUND: vignettes/pipeline.Rmd
- FOUND: R/analysis_utils.R
- FOUND: commit 6690aba
- FOUND: commit 5785499

---
*Quick Task: 1-fix-vignette-warnings-suppress-stan-comp*
*Completed: 2026-02-13*
