---
phase: 09-describe-helpers
plan: 01
subsystem: introspection
tags: [rbmi, draws, S3, cli, MCMC, rstan, TDD]

# Dependency graph
requires:
  - phase: 08-mi-diagnostic-statistics
    provides: "MI diagnostic patterns and ARD enrichment"
provides:
  - "describe_draws() exported function for inspecting rbmi draws objects"
  - "print.describe_draws() S3 method with cli formatting"
  - "Mock draws object helpers for testing"
affects: [09-02-PLAN (describe_imputation), future CLI/reporting phases]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "describe_* pattern: extract metadata from rbmi pipeline objects into S3 classes with print methods"
    - "Mock rbmi object construction for unit testing draws/impute objects"
    - "cli output capture via type='message' in tests"

key-files:
  created:
    - R/describe.R
    - tests/testthat/test-describe.R
    - man/describe_draws.Rd
    - man/print.describe_draws.Rd
  modified:
    - NAMESPACE
    - inst/WORDLIST

key-decisions:
  - "09-01-D1: Method name includes subtype for condmean: 'Conditional Mean (jackknife)' not just 'Conditional Mean'"
  - "09-01-D2: Condmean sample display uses n_primary/n_resampled fields with '1 + N' print format"
  - "09-01-D3: All-NA Rhat converged = NA (not TRUE), handled explicitly to avoid misleading result"
  - "09-01-D4: approxbayes has no bayes_control (only bayes has it) since approxbayes uses different internals"
  - "09-01-D5: cli output captured via type='message' in tests since cli writes to stderr connection"

patterns-established:
  - "describe_* functions: take rbmi pipeline object, return S3 class with print method"
  - "Mock draws helpers: make_mock_draws_{condmean,bayes,approxbayes} for testing"

# Metrics
duration: 7min
completed: 2026-02-11
---

# Phase 9 Plan 01: describe_draws Summary

**describe_draws() extracts method, formula, samples, failures, covariance, and MCMC diagnostics from rbmi draws objects with cli-formatted print method**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-11T05:43:31Z
- **Completed:** 2026-02-11T05:50:26Z
- **Tasks:** 3 (TDD: RED, GREEN, REFACTOR)
- **Files modified:** 6

## Accomplishments
- describe_draws() exported function covering all 3 rbmi method types (condmean, bayes, approxbayes)
- print.describe_draws() with cli formatting showing "1 + N" sample format for condmean
- MCMC diagnostics extraction from stanfit when rstan available, with graceful degradation
- 53 passing tests covering input validation, all method types, print output, and edge cases
- R CMD check: 0 errors, 0 warnings

## Task Commits

Each task was committed atomically (TDD cycle):

1. **RED: Failing tests** - `616bf7f` (test)
2. **GREEN: Implementation** - `c37a380` (feat)
3. **REFACTOR: Cleanup** - `0b300a8` (refactor)

_TDD cycle: test -> feat -> refactor_

## Files Created/Modified
- `R/describe.R` - describe_draws() and print.describe_draws() with roxygen2 docs
- `tests/testthat/test-describe.R` - 53 tests with mock draws object helpers
- `man/describe_draws.Rd` - Generated roxygen2 documentation
- `man/print.describe_draws.Rd` - Generated roxygen2 documentation
- `NAMESPACE` - Added export(describe_draws) and S3method(print,describe_draws)
- `inst/WORDLIST` - Added technical terms (ESS, Rhat, MCMC, stanfit, etc.)

## Decisions Made
- **09-01-D1:** Method name includes subtype for condmean: "Conditional Mean (jackknife)" rather than just "Conditional Mean", matching rbmi's own display convention
- **09-01-D2:** Condmean sample display stores n_primary and n_resampled fields; print method formats as "1 + N (1 primary + N jackknife/bootstrap)"
- **09-01-D3:** When all Rhat values are NA, converged is set to NA (not TRUE), avoiding misleading convergence claims
- **09-01-D4:** approxbayes does not expose bayes_control since it uses different internals than full Bayesian; only method_bayes has control$warmup etc.
- **09-01-D5:** cli output writes to stderr, so test capture uses capture.output(type = "message")

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- cli output capture: `capture.output()` defaults to stdout, but cli writes to stderr/message connection. Fixed by using `capture.output(..., type = "message")` in tests. Not a deviation -- discovered during normal RED-GREEN iteration.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- describe_draws() pattern established for 09-02 (describe_imputation)
- Mock object helpers in test file can be referenced/reused
- NAMESPACE and documentation infrastructure updated
- No blockers for describe_imputation implementation

---
*Phase: 09-describe-helpers*
*Completed: 2026-02-11*
