---
phase: 05
plan: 01
subsystem: data-validation
tags: [cli, error-handling, validation, data-prep]

dependency_graph:
  requires: []
  provides:
    - "Hardened validate_data() with cli messaging and 6 new validation checks"
    - "Hardened prepare_data_ice() with cli messaging and HRD-02/HRD-03 checks"
    - "75 passing tests covering all new and existing validation behavior"
  affects:
    - "05-02 (edge case test coverage builds on this hardened foundation)"
    - "06-xx (documentation vignette will reference improved error messages)"

tech_stack:
  added: []
  patterns:
    - "cli::cli_abort with custom error classes (rbmiUtils_error_validation)"
    - "cli::cli_warn with custom warning classes (rbmiUtils_warning_coercion)"
    - "cli::cli_inform with informational class (rbmiUtils_info)"
    - "Batched warnings via cli_warn (collect then emit once)"
    - "Fail-fast for fatal errors, collect-then-report for non-fatal"

key_files:
  created: []
  modified:
    - R/data_helpers.R
    - tests/testthat/test-data_helpers.R

decisions:
  - id: D-05-01-01
    choice: "All-NA covariate columns warn and skip validation (not error)"
    rationale: "validate_data() returns invisible(TRUE), not modified vars; warn guides user to remove column"
  - id: D-05-01-02
    choice: "All-NA outcome column is a hard error (collected issue)"
    rationale: "Nothing to impute or model -- analysis cannot proceed"
  - id: D-05-01-03
    choice: "Used stats::setNames() to avoid R CMD check NOTE"
    rationale: "setNames is from stats, not base; needs explicit namespace qualification"

metrics:
  duration: "8 minutes"
  completed: "2026-02-08"
---

# Phase 5 Plan 01: Data Prep Hardening -- cli Messaging & Validation Checks

**One-liner:** Migrated validate_data() and prepare_data_ice() from stop()/assertthat to cli messaging with 6 new validation checks (HRD-01 through HRD-06) and 75 passing tests.

## Performance

- **Duration:** ~8 minutes
- **Tests:** 75 pass (was 63, added 12 new)
- **R CMD check:** 0 errors, 0 warnings, 1 note (unrelated timestamp)

## Accomplishments

1. **validate_data() fully migrated to cli messaging** -- all `stop()`, `warning()`, and `assertthat` calls replaced with `cli::cli_abort()`, `cli::cli_warn()`, and `cli::cli_inform()` using custom error/warning classes.

2. **prepare_data_ice() fully migrated to cli messaging** -- all `assertthat::assert_that()` and `stop()` calls replaced with cli equivalents using typed error classes (`rbmiUtils_error_type`, `rbmiUtils_error_validation`).

3. **HRD-01: Malformed interaction term validation** -- empty terms and terms with leading/trailing/consecutive operators (`A*`, `:B`, `A**B`) rejected with clear cli-formatted errors before parsing.

4. **HRD-02: NULL strategy error** -- `prepare_data_ice()` now errors immediately when `vars$strategy` is NULL instead of silently defaulting to `"strategy"`, with guidance to use `rbmi::set_vars()`.

5. **HRD-03: Character visit warning** -- `prepare_data_ice()` warns when visit column is character with guidance on factor conversion and explicit level ordering.

6. **HRD-04: Empty data frame check** -- `validate_data()` errors immediately on 0-row data frames with guidance to provide at least one subject-visit observation.

7. **HRD-05: All-NA covariate warning** -- all-NA covariate columns now warn (not error) and are excluded from validation, while partially-NA covariates still error as before.

8. **HRD-06: Batched character column warnings** -- multiple character columns are reported in a single `cli::cli_warn()` call instead of individual `warning()` calls per column.

9. **All-NA outcome and all-complete informational messages** -- all-NA outcome is a hard error; all-complete data emits `cli::cli_inform()` confirming no missing data.

10. **12 new tests added** covering all HRD checks plus all-NA outcome; all existing tests updated to use class-based matching.

## Task Commits

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Harden validate_data() | `778e226` | R/data_helpers.R |
| 2 | Harden prepare_data_ice() | `9902194` | R/data_helpers.R |
| 3 | Update and add tests | `a4b691d` | tests/testthat/test-data_helpers.R, R/data_helpers.R |

## Files Modified

| File | Changes |
|------|---------|
| `R/data_helpers.R` | Rewrote validate_data() and prepare_data_ice() with cli messaging, added 6 new validation checks |
| `tests/testthat/test-data_helpers.R` | Updated 12 existing tests for class-based matching, added 12 new tests for HRD-01 through HRD-06 |

## Decisions Made

| ID | Decision | Rationale |
|----|----------|-----------|
| D-05-01-01 | All-NA covariates warn+skip, not error | validate_data returns TRUE/error, cannot modify vars; user guidance via warning |
| D-05-01-02 | All-NA outcome is hard error | Nothing to model, analysis cannot proceed |
| D-05-01-03 | stats::setNames() for namespace compliance | Avoids R CMD check NOTE for undefined global |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed cli pluralization error in all-NA covariate warning**
- **Found during:** Task 3 (test execution)
- **Issue:** cli `{?s}` pluralization in multi-bullet warning messages crashed with "Cannot pluralize without a quantity" when the quantity was in a different bullet than the pluralization marker
- **Fix:** Restructured the warning message to keep quantity and pluralization in the same bullet line
- **Files modified:** R/data_helpers.R
- **Commit:** `a4b691d`

**2. [Rule 3 - Blocking] Fixed R CMD check NOTE for setNames**
- **Found during:** Task 3 (R CMD check)
- **Issue:** `setNames()` used without namespace qualification triggers NOTE: "no visible global function definition for 'setNames'"
- **Fix:** Changed to `stats::setNames()`
- **Files modified:** R/data_helpers.R
- **Commit:** `a4b691d`

## Issues Identified

None.

## Next Phase Readiness

- **05-02 (Edge Case Tests):** Ready. The hardened functions provide the foundation for edge case testing (single subject, single visit, all-complete data scenarios).
- **Blockers:** None.
- **Concerns:** None.
