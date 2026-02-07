---
phase: 01-foundation-hardening
plan: 03
subsystem: gcomp, imputation_storage
tags: [input-validation, cli, digest, round-trip, beeca, error-handling]

dependency_graph:
  requires:
    - "01-01 (cli available as Import)"
  provides:
    - "gcomp_responder() and gcomp_binary() validate inputs with clear error messages"
    - "beeca output schema validated before column access"
    - "reduce/expand_imputed_data() perform round-trip digest verification"
    - "Custom error classes: rbmiUtils_error_validation, rbmiUtils_error_type, rbmiUtils_error_dependency, rbmiUtils_error_integrity"
  affects:
    - "Phase 2 (stable storage functions for ARD workflows)"
    - "Phase 3 (gcomp results can be reliably produced for table generation)"

tech_stack:
  added: []
  patterns:
    - "rlang::hash() for round-trip digest verification (XXH128)"
    - "Attribute-based metadata storage (rbmiUtils_digest, rbmiUtils_col_metadata, rbmiUtils_col_names)"
    - "Fail-early validation with cli::cli_abort() and custom error classes"
    - "beeca output schema validation (marginal_results, STAT, STATVAL, TRTVAL)"

file_tracking:
  created: []
  modified:
    - "R/analysis_utils.R"
    - "R/utils.R"
    - "R/imputation_storage.R"
    - "tests/testthat/test-analysis_utils.R"
    - "tests/testthat/test-imputation_storage.R"

decisions:
  - id: "01-03-D1"
    decision: "Integrity check only verifies columns present in both stored metadata AND original_data"
    reason: "imputed_data may contain derived columns (e.g., CRIT1FLN) not present in original_data. These columns won't appear in expanded result since expansion builds from original_data. Checking them would cause false positives."

metrics:
  duration: "~8 minutes (recovery from interrupted agent)"
  completed: "2026-02-07"
---

# Phase 01 Plan 03: Harden gcomp Validation and Storage Digest Summary

**Add fail-early input validation to gcomp functions, pin beeca output schema, and implement round-trip digest verification for imputation storage**

## Performance

- **Completed:** 2026-02-07
- **Duration:** ~8 minutes (recovery from interrupted agent)

## Accomplishments

1. **Added input validation to gcomp_responder()** -- Validates: data.frame type, required columns exist, >= 2 treatment arms, numeric outcome, non-zero variance outcome. All errors use `cli::cli_abort()` with custom error classes.

2. **Added input validation to gcomp_binary()** -- Same validation pattern as gcomp_responder with adapted parameter names (explicit `outcome`, `treatment`, `covariates` instead of `vars`).

3. **Added beeca output schema validation** -- After `beeca::get_marginal_effect()` call, validates that `marginal_results` exists and contains expected columns (STAT, STATVAL, TRTVAL). Uses `rbmiUtils_error_dependency` class for beeca version incompatibility.

4. **Added minimal validation to gcomp_responder_multi()** -- Checks visit_var exists in data and at least one visit present. Per-visit validation delegated to gcomp_responder().

5. **Implemented round-trip digest verification** -- `reduce_imputed_data()` computes `rlang::hash()` digest of column metadata, first imputation hash, column names, and imputation count. Stored as attributes. `expand_imputed_data()` verifies column names and types match on expansion.

6. **Migrated storage function errors to cli::cli_abort()** -- All `stop()` and `assertthat::assert_that()` calls replaced with `cli::cli_abort()` using `rbmiUtils_error_type`, `rbmiUtils_error_validation`, and `rbmiUtils_error_integrity` classes.

7. **Added comprehensive tests** -- Edge case tests for gcomp functions (single-arm, missing columns, zero-variance, non-numeric). Round-trip tests for digest storage, exact data preservation, factor level preservation, type mismatch detection, and error class validation.

## Task Commits

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Add input validation to gcomp functions and beeca output pinning | `06e4075` | R/analysis_utils.R, R/utils.R, tests/testthat/test-analysis_utils.R |
| 2 | Add round-trip digest verification to storage functions | `55a97ad` | R/imputation_storage.R, tests/testthat/test-imputation_storage.R |
| 3 | Fix integrity check for mismatched column sets | `a526f8e` | R/imputation_storage.R, tests/testthat/test-imputation_storage.R |

## Files Modified

- **R/analysis_utils.R** -- Added fail-early validation to gcomp_responder() and gcomp_responder_multi(), beeca output schema check
- **R/utils.R** -- Added fail-early validation to gcomp_binary(), beeca output schema check
- **R/imputation_storage.R** -- Added rlang::hash() digest computation/verification, replaced stop()/assertthat with cli::cli_abort(), relaxed column check to intersect with original_data
- **tests/testthat/test-analysis_utils.R** -- Added 5 edge case tests (single-arm, missing cols, zero-variance, non-numeric, gcomp_binary)
- **tests/testthat/test-imputation_storage.R** -- Added 7 new tests (digest storage, col_metadata, round-trip exact, factor levels, extra columns tolerance, type mismatch, cli error classes), updated 4 existing tests to class matching

## Decisions Made

1. **Integrity check scoped to shared columns** (01-03-D1): Column name verification in `expand_imputed_data()` only checks columns that appear in both stored metadata and `original_data`. This prevents false positives when `imputed_data` had derived columns (e.g., CRIT1FLN, CRIT1FL) not present in `original_data`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Column name integrity check too strict for ADMI/ADEFF data**
- **Found during:** Task 2 verification (test execution)
- **Issue:** ADMI contains columns (CRIT1FLN, CRIT1FL, CRIT) not present in ADEFF (original_data). The strict column name check rejected valid round-trips.
- **Fix:** Changed check to `intersect(stored_col_names, orig_cols)` — only verify columns that exist in both
- **Files modified:** R/imputation_storage.R
- **Commit:** a526f8e

**2. [Rule 1 - Bug] cli pluralization error in type mismatch message**
- **Found during:** Task 2 verification
- **Issue:** `"column type mismatch{?es}"` caused cli pluralization error ("Cannot pluralize without a quantity")
- **Fix:** Changed to static string and used named vector for bullet points
- **Commit:** a526f8e

**3. [Rule 1 - Bug] Existing tests used string matching incompatible with cli formatting**
- **Found during:** Task 2 verification
- **Issue:** Old tests used `expect_error(..., "must be a data.frame")` but cli formats as `must be a <data.frame>`
- **Fix:** Updated to `expect_error(..., class = "rbmiUtils_error_type")`
- **Commit:** a526f8e

## Issues Encountered

None blocking.

## Next Phase Readiness

- **Blockers:** None
- **Ready for:** Phase 2 (stable, validated gcomp + storage functions)
- **Concerns:** None
