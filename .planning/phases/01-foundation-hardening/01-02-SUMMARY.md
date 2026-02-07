---
phase: 01-foundation-hardening
plan: 02
subsystem: analyse_mi_data
tags: [inherits, cli, lifecycle, deprecation, error-handling]

dependency_graph:
  requires:
    - "01-01 (cli and lifecycle available as Imports)"
  provides:
    - "analyse_mi_data() uses inherits() for robust method type detection"
    - "All errors use cli::cli_abort() with custom error classes"
    - "as_analysis2() deprecated with lifecycle warning"
  affects:
    - "Phase 2 (print/summary methods can rely on inherits-based detection)"

tech_stack:
  added: []
  patterns:
    - "inherits() for S3 class detection instead of positional class()[[2]] indexing"
    - "lifecycle::deprecate_warn() for internal helper deprecation"
    - "suppressWarnings(..., classes = 'lifecycle_warning_deprecated') for internal calls"
    - "cli::cli_warn() for non-fatal warnings"

file_tracking:
  created: []
  modified:
    - "R/analyse_mi_data.R"
    - "tests/testthat/test-analyse_mi_data.R"

decisions:
  - id: "01-02-D1"
    decision: "Suppress lifecycle deprecation warning when as_analysis2() is called internally by analyse_mi_data()"
    reason: "The deprecation is meant for direct callers of as_analysis2(). Without suppression, every analyse_mi_data() call triggers a noisy deprecation warning that confuses users."

metrics:
  duration: "~5 minutes (recovery from interrupted agent)"
  completed: "2026-02-07"
---

# Phase 01 Plan 02: Refactor analyse_mi_data to use inherits() Summary

**Replace fragile class()[[2]] positional indexing with inherits()-based method detection, migrate all errors to cli::cli_abort(), and deprecate internal helpers**

## Performance

- **Completed:** 2026-02-07
- **Duration:** ~5 minutes (recovery from interrupted agent)

## Accomplishments

1. **Replaced all class()[[2]] usage with inherits()** -- Four locations in analyse_mi_data.R converted: `analyse_mi_data()` n_expected extraction, `as_analysis2()` next_class determination, `print.analysis()` method display, and `summary.analysis()` method display. This is robust against rbmi class vector reordering.

2. **Migrated all error messages to cli::cli_abort()** -- All `stop()` and `assertthat::assert_that()` calls replaced with `cli::cli_abort()` using custom error classes: `rbmiUtils_error_validation`, `rbmiUtils_error_type`, `rbmiUtils_error_dependency`, `rbmiUtils_error_internal`. The too-many-imputations warning upgraded to `cli::cli_warn()`.

3. **Added lifecycle deprecation to as_analysis2()** -- `lifecycle::deprecate_warn("0.2.0", "as_analysis2()")` added with roxygen badge. Internal call from `analyse_mi_data()` wrapped in `suppressWarnings(..., classes = "lifecycle_warning_deprecated")` to avoid user-facing noise.

4. **Updated all tests to class-based error matching** -- Existing tests updated from message-string matching to `class = "rbmiUtils_error_validation"` and `class = "rbmiUtils_error_type"`. Added new tests for: inherits()-based method detection (bayes + approxbayes), deprecation warning detection, imputation count warning, and rbmi::pool() compatibility.

## Task Commits

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Refactor analyse_mi_data method detection and error messaging | `8419b54` | R/analyse_mi_data.R |
| 2 | Update tests for cli error classes and add inherits/deprecation tests | `06c4130` | tests/testthat/test-analyse_mi_data.R, R/analyse_mi_data.R |

## Files Modified

- **R/analyse_mi_data.R** -- Replaced class()[[2]] with inherits(), stop()/assertthat with cli::cli_abort(), added lifecycle deprecation to as_analysis2(), suppressed internal deprecation warning
- **tests/testthat/test-analyse_mi_data.R** -- Updated 7 existing error tests to use class matching, added 4 new tests (inherits detection, deprecation, imputation warning, pool compatibility)

## Decisions Made

1. **Suppress internal deprecation warning** (01-02-D1): `as_analysis2()` deprecation warning suppressed when called from `analyse_mi_data()` since the deprecation targets direct callers, not the primary API.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Noisy deprecation warnings in all tests**
- **Found during:** Task 2 (test execution)
- **Issue:** Every `analyse_mi_data()` call triggered a lifecycle deprecation warning from the internal `as_analysis2()` call
- **Fix:** Wrapped the internal `as_analysis2()` call in `suppressWarnings(..., classes = "lifecycle_warning_deprecated")`
- **Files modified:** R/analyse_mi_data.R
- **Commit:** 06c4130

**2. [Rule 1 - Bug] Pool compatibility test required factor columns**
- **Found during:** Task 2 (test execution)
- **Issue:** `rbmi::ancova` requires group variable to be a factor with 2 levels
- **Fix:** Added factor conversion for TRT, USUBJID, AVISIT in the pool compatibility test
- **Commit:** 06c4130

## Issues Encountered

None blocking.

## Next Phase Readiness

- **Blockers:** None
- **Ready for:** Phase 2 (print/summary methods already use inherits(), ARD conversion can rely on stable analysis objects)
- **Concerns:** None
