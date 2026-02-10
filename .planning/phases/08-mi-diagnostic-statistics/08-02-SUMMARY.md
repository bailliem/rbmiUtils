---
phase: 08-mi-diagnostic-statistics
plan: 02
subsystem: statistics
tags: [ard, mi-diagnostics, fmi, lambda, pool-to-ard, cards, enriched-ard]

# Dependency graph
requires:
  - phase: 04-ard-conversion
    provides: pool_to_ard() function and base ARD construction
  - phase: 08-01
    provides: compute_rubin_diagnostics() internal function for MI variance decomposition
provides:
  - Enriched pool_to_ard(pool_obj, analysis_obj) producing ARD with MI diagnostic stat rows
  - compute_mi_diagnostics() orchestrator for per-parameter diagnostic computation
  - build_diagnostic_ard_rows() ARD row constructor for diagnostic statistics
  - Comprehensive test suite (7 new test blocks, 176 expectations) for enriched ARD
affects:
  - 09 (efficacy_table may encounter enriched ARD)
  - 10 (forest plot may encounter enriched ARD)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "ARD enrichment pattern: optional analysis_obj parameter enables diagnostic rows"
    - "Internal orchestrator pattern: compute_mi_diagnostics dispatches to compute_rubin_diagnostics per parameter"
    - "Non-Rubin graceful degradation: cli_inform message, no NA rows, clean omission"

key-files:
  created: []
  modified:
    - R/ard_conversion.R
    - tests/testthat/test-ard_conversion.R
    - man/pool_to_ard.Rd

key-decisions:
  - "Curated essentials only in ARD: fmi, lambda, riv, df.adjusted, df.complete, re, m.imputations (no var.within/var.between/var.total)"
  - "Non-Rubin pooling returns base ARD with cli_inform message, no diagnostic rows or NA rows"
  - "Diagnostic stat names follow mice convention: lowercase dot-separated (df.adjusted, df.complete, m.imputations)"
  - "analysis_obj validation checks parameter name match between pool_obj$pars and analysis_obj$results[[1]]"

patterns-established:
  - "Optional enrichment via second parameter: pool_to_ard(pool_obj, analysis_obj) backward-compatible API"
  - "Per-parameter diagnostics extracted via vapply over analysis_obj$results"
  - "Diagnostic ARD rows share same grouping columns as base ARD rows for consistent filtering"

# Metrics
duration: 9min
completed: 2026-02-10
---

# Phase 8 Plan 02: ARD Diagnostic Enrichment Summary

**pool_to_ard() enriched with MI diagnostics via analysis_obj parameter: 7 diagnostic stat rows per parameter (FMI, lambda, RIV, Barnard-Rubin df, RE) with backward-compatible API and cards validation**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-10T22:19:29Z
- **Completed:** 2026-02-10T22:27:58Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments
- Extended pool_to_ard() with optional analysis_obj parameter for MI diagnostic enrichment
- Created compute_mi_diagnostics() and build_diagnostic_ard_rows() internal helpers
- Added 7 comprehensive test blocks (176 new expectations) covering diagnostic presence, backward compatibility, ARD validation, non-Rubin handling, method row, numeric reasonableness, and parameter name validation
- Enriched ARD passes cards::check_ard_structure() validation
- R CMD check: 0 errors, 0 warnings

## Task Commits

Each task was committed atomically:

1. **Task 1: Modify pool_to_ard() signature and add parameter validation** - `0b46bef` (feat)
2. **Task 2: Create compute_mi_diagnostics() and build_diagnostic_ard_rows() helpers** - `50d27d2` (feat)
3. **Task 3: Add comprehensive tests for enriched ARD diagnostics** - `737efad` (test)

## Files Created/Modified
- `R/ard_conversion.R` - Added analysis_obj parameter, validation, compute_mi_diagnostics(), build_diagnostic_ard_rows(), and diagnostic enrichment integration
- `tests/testthat/test-ard_conversion.R` - Added make_mock_analysis() helper and 7 test blocks for enriched ARD
- `man/pool_to_ard.Rd` - Updated roxygen documentation for analysis_obj parameter

## Decisions Made
- Curated essentials only: 7 diagnostic stat rows (fmi, lambda, riv, df.adjusted, df.complete, re, m.imputations) -- no variance components (var.within, var.between, var.total) per locked context decision
- Non-Rubin pooling omits diagnostics entirely (returns NULL from compute_mi_diagnostics) with cli_inform message -- no NA placeholder rows
- Stat naming convention follows mice package: lowercase dot-separated names (df.adjusted, df.complete, m.imputations)
- Parameter validation uses sorted name comparison between pool_obj$pars and analysis_obj$results[[1]]
- Complete-data df extracted as unique(dfs)[1] (constant across imputations for a given parameter)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 8 (MI Diagnostic Statistics) is complete: both plans executed successfully
- Users can now call pool_to_ard(pool_obj, analysis_obj) to get enriched ARD with MI diagnostics
- Existing pool_to_ard(pool_obj) API unchanged (backward compatible)
- Ready for phase 9+ work (enriched ARD may flow through efficacy_table and forest plots)

---
*Phase: 08-mi-diagnostic-statistics*
*Completed: 2026-02-10*
