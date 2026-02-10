---
phase: 08-mi-diagnostic-statistics
plan: 01
subsystem: statistics
tags: [rubin-rules, mi-diagnostics, fmi, lambda, barnard-rubin, variance-decomposition]

# Dependency graph
requires:
  - phase: 04-ard-conversion
    provides: pool_to_ard() function and R/ard_conversion.R file structure
provides:
  - compute_rubin_diagnostics() internal function for MI variance decomposition
  - Comprehensive test suite for Rubin's rules diagnostic computation
affects:
  - 08-02 (pool_to_ard enrichment will call compute_rubin_diagnostics)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Pure computational helper pattern: internal function with no side effects"
    - "stats::var() prefix convention for base R stats functions"
    - "Edge case handling for MI diagnostics: lambda=0, Inf df, NA inputs"

key-files:
  created: []
  modified:
    - R/ard_conversion.R
    - tests/testthat/test-ard_conversion.R

key-decisions:
  - "Use stats::var() prefix to match codebase convention and avoid R CMD check NOTE"
  - "Adjusted FMI formula follows mice package convention: (riv + 2/(df+3)) / (1 + riv)"
  - "var_b computed in all-NA-SE case (estimates still available), all other diagnostics set NA"

patterns-established:
  - "Rubin's rules edge case logic: check NA v_com, then Inf+var_b==0, then lambda==0 branching"

# Metrics
duration: 10min
completed: 2026-02-10
---

# Phase 8 Plan 01: Rubin's Rules Diagnostic Computation Summary

**Pure Rubin's rules variance decomposition function (FMI, lambda, RIV, Barnard-Rubin df, RE) with 7 TDD test cases covering known values and edge cases**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-10T22:07:10Z
- **Completed:** 2026-02-10T22:17:40Z
- **Tasks:** 3 (TDD: RED, GREEN, REFACTOR)
- **Files modified:** 2

## Accomplishments
- Implemented `compute_rubin_diagnostics()` internal function computing var_w, var_b, var_t, lambda, riv, df_adj, dfcom, fmi, and re
- Wrote 7 comprehensive test cases: known values with hand-calculated expectations, lambda=0, v_com=Inf, v_com=Inf+lambda=0, all-NA SEs, fmi-vs-lambda distinction, NA v_com
- All formulas verified against rbmi:::rubin_rules and rbmi:::rubin_df source code
- R CMD check passes with 0 errors, 0 warnings

## Task Commits

Each task was committed atomically (TDD cycle):

1. **RED: Write failing tests** - `ff5fbda` (test)
2. **GREEN: Implement function** - `4b6a806` (feat)
3. **REFACTOR: stats:: prefix fix** - `c3581c9` (refactor)

_TDD cycle: 3 commits (test -> feat -> refactor)_

## Files Created/Modified
- `R/ard_conversion.R` - Added `compute_rubin_diagnostics()` internal function (non-exported, @noRd)
- `tests/testthat/test-ard_conversion.R` - Added 7 test_that blocks with 49+ expectations for diagnostic computation

## Decisions Made
- Used `stats::var()` prefix instead of bare `var()` to match existing codebase convention (all stats functions use stats:: prefix)
- Adjusted FMI formula `(riv + 2/(df+3)) / (1+riv)` follows mice package convention, distinct from lambda
- In the all-NA-SEs case, `var_b` is still computed from estimates (available), while all other diagnostics return NA

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Used stats::var() to resolve R CMD check NOTE**
- **Found during:** REFACTOR phase (R CMD check verification)
- **Issue:** Using bare `var()` in internal function caused R CMD check NOTE about undefined global
- **Fix:** Prefixed with `stats::` to match codebase convention
- **Files modified:** R/ard_conversion.R
- **Verification:** R CMD check 0 errors, 0 warnings, NOTE resolved
- **Committed in:** c3581c9

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Minor namespace convention fix. No scope change.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `compute_rubin_diagnostics()` is ready to be called by Plan 02's `pool_to_ard()` enrichment
- Function signature matches the API designed in RESEARCH.md Pattern 2: `compute_rubin_diagnostics(ests, ses, v_com, M)`
- All edge cases handled and tested, matching rbmi:::rubin_df() behavior

---
*Phase: 08-mi-diagnostic-statistics*
*Completed: 2026-02-10*
