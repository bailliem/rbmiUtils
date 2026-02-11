---
phase: 09-describe-helpers
plan: 02
subsystem: introspection
tags: [rbmi, imputation, S3, cli, missingness, TDD]

# Dependency graph
requires:
  - phase: 09-describe-helpers
    plan: 01
    provides: "describe_draws() pattern, R/describe.R file, test infrastructure"
provides:
  - "describe_imputation() exported function for inspecting rbmi imputation objects"
  - "print.describe_imputation() S3 method with cli formatting"
  - "Missingness breakdown by visit and treatment arm"
affects: [future reporting/CLI phases, vignette documentation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "describe_imputation pattern: extract metadata + missingness from rbmi imputation objects"
    - "cli::cli_verbatim() for routing base print output through cli message connection"

key-files:
  created:
    - man/describe_imputation.Rd
    - man/print.describe_imputation.Rd
  modified:
    - R/describe.R
    - tests/testthat/test-describe.R
    - NAMESPACE
    - inst/WORDLIST

key-decisions:
  - "09-02-D1: Missingness table uses base expand.grid + loop aggregation rather than dplyr, keeping zero external dependencies"
  - "09-02-D2: cli::cli_verbatim() routes data.frame print output through stderr message connection for consistent cli capture"
  - "09-02-D3: Mock impute helper derives IDs from groups when custom groups provided, avoiding ID mismatch"

patterns-established:
  - "describe_* functions: complete set for draws + imputation pipeline objects"
  - "Mock impute helper: make_mock_impute() with environment-based longdata mock"

# Metrics
duration: 7min
completed: 2026-02-11
---

# Phase 9 Plan 02: describe_imputation Summary

**describe_imputation() extracts method, M, references, and missingness breakdown by visit/arm from rbmi imputation objects with cli-formatted print method**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-11T05:52:21Z
- **Completed:** 2026-02-11T05:59:18Z
- **Tasks:** 3 (TDD: RED, GREEN, REFACTOR)
- **Files modified:** 6

## Accomplishments
- describe_imputation() exported function extracting method, n_imputations, references, subjects, visits, and missingness
- Missingness data.frame with cross-tabulation of is_missing by visit and treatment group (n_total, n_miss, pct_miss)
- print.describe_imputation() with cli formatting: method, M, subjects, references table, missingness breakdown
- Graceful NULL references handling with "No explicit references" message
- 98 total passing tests (53 describe_draws + 45 describe_imputation)
- R CMD check: 0 errors, 0 warnings

## Task Commits

Each task was committed atomically (TDD cycle):

1. **RED: Failing tests** - `38e6b4c` (test)
2. **GREEN: Implementation** - `63a90fb` (feat)
3. **REFACTOR: Cleanup** - `deb6581` (refactor)

_TDD cycle: test -> feat -> refactor_

## Files Created/Modified
- `R/describe.R` - describe_imputation() and print.describe_imputation() appended to existing file
- `tests/testthat/test-describe.R` - 45 new tests with make_mock_impute() helper
- `man/describe_imputation.Rd` - Generated roxygen2 documentation
- `man/print.describe_imputation.Rd` - Generated roxygen2 documentation
- `NAMESPACE` - Added export(describe_imputation) and S3method(print,describe_imputation)
- `inst/WORDLIST` - Updated spelling wordlist

## Decisions Made
- **09-02-D1:** Missingness aggregation uses base R expand.grid + for loop rather than dplyr, keeping the package dependency-light
- **09-02-D2:** Used cli::cli_verbatim() to route data.frame print output through the stderr/message connection, ensuring consistent output capture in tests
- **09-02-D3:** Mock impute helper derives subject IDs from groups parameter names when custom groups are provided, preventing ID mismatch errors

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Mock impute helper ID mismatch**
- **Found during:** Task 2 (GREEN phase)
- **Issue:** make_mock_impute() generated SUBJ1..SUBJN IDs but custom groups/missing_pattern used S1..SN keys, causing vapply to fail on lookup
- **Fix:** Derive IDs from names(groups) when custom groups are provided
- **Files modified:** tests/testthat/test-describe.R
- **Verification:** All 98 tests pass
- **Committed in:** 63a90fb (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug in test helper)
**Impact on plan:** Minor test helper fix, no scope creep.

## Issues Encountered
- Missingness table printed via base print() outputs to stdout while cli outputs to stderr. Refactored to use cli::cli_verbatim() in REFACTOR phase to keep all output on consistent connection.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Both describe_* functions complete and coexisting in R/describe.R
- Phase 9 (describe-helpers) fully complete
- No blockers for future phases

---
*Phase: 09-describe-helpers*
*Completed: 2026-02-11*
