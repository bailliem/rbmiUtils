---
phase: 02-print-summary-ard-conversion
plan: 03
subsystem: api
tags: [cards, ard, pharmaverse, pool, conversion]

# Dependency graph
requires:
  - phase: 01-foundation-hardening
    provides: tidy_pool_obj() for tidying pool objects
provides:
  - pool_to_ard() function converting pool objects to pharmaverse ARD format
  - cards package added to Suggests
  - Testable dependency guard pattern (is_cards_available helper)
affects:
  - 03-gtsummary-table-generation (ARD input for gtsummary table building)

# Tech tracking
tech-stack:
  added: [cards (Suggests)]
  patterns: [as_card() direct ARD construction, requireNamespace guard with testable helper, list columns for ARD stat/level/warning/error]

key-files:
  created:
    - R/ard_conversion.R
    - tests/testthat/test-ard_conversion.R
    - man/pool_to_ard.Rd
  modified:
    - DESCRIPTION
    - NAMESPACE

key-decisions:
  - "02-03-D1: Used as_card() direct construction over ard_identity() for batch efficiency"
  - "02-03-D2: Extracted is_cards_available() helper for testable dependency guard"
  - "02-03-D3: Method rows included per parameter to satisfy check_ard_structure() fully"

patterns-established:
  - "ARD construction: build data.frame with I() list columns then as_card() + tidy_ard_column_order()"
  - "Suggests dependency: requireNamespace wrapped in internal helper for testability"

# Metrics
duration: 6min
completed: 2026-02-08
---

# Phase 02 Plan 03: ARD Conversion Summary

**pool_to_ard() converts rbmi pool objects to pharmaverse ARD standard using cards::as_card() with visit/parameter_type/lsm_type grouping columns**

## Performance

- **Duration:** ~6 min
- **Started:** 2026-02-08T06:58:58Z
- **Completed:** 2026-02-08T07:04:37Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- pool_to_ard() converts pool objects to valid ARD data frames passing cards::check_ard_structure()
- ARD contains 6 statistics per parameter (estimate, std.error, conf.low, conf.high, p.value, method)
- Grouping columns preserve visit (group1), parameter_type (group2), lsm_type (group3) per ARD-03 spec
- 8 test blocks with 75 assertions, R CMD check clean (0 errors, 0 warnings)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add cards to Suggests and create pool_to_ard() function** - `40651c1` (feat)
2. **Task 2: Add comprehensive tests for pool_to_ard()** - `0ec075c` (test)

## Files Created/Modified
- `R/ard_conversion.R` - pool_to_ard() function with ARD construction, is_cards_available() helper
- `tests/testthat/test-ard_conversion.R` - 8 test blocks covering structure, grouping, list columns, NA handling, dependency guard, validation, numeric accuracy
- `man/pool_to_ard.Rd` - Generated documentation
- `DESCRIPTION` - cards added to Suggests (alphabetical order)
- `NAMESPACE` - export(pool_to_ard) added

## Decisions Made

- **02-03-D1: as_card() direct construction over ard_identity()** - Direct data.frame construction with I() list columns then as_card() is more efficient for batch conversion of multiple parameters than per-row ard_identity() calls
- **02-03-D2: Extracted is_cards_available() helper** - requireNamespace("cards") wrapped in internal helper function so tests can mock it via local_mocked_bindings without infinite recursion from mocking base::requireNamespace
- **02-03-D3: Method rows per parameter** - Each parameter gets its own method row (stat_name = "method", stat = pool_obj$method) to fully satisfy check_ard_structure() without notes

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Extracted is_cards_available() internal helper for testable dependency guard**
- **Found during:** Task 2 (test for dependency error)
- **Issue:** Mocking base::requireNamespace via local_mocked_bindings caused infinite recursion (stack overflow) because the mock called base::requireNamespace for non-cards packages
- **Fix:** Extracted requireNamespace("cards") into internal helper `is_cards_available()` in R/ard_conversion.R, making it mockable via local_mocked_bindings without touching base namespace
- **Files modified:** R/ard_conversion.R, tests/testthat/test-ard_conversion.R
- **Verification:** Test passes, mock correctly triggers dependency error class
- **Committed in:** 0ec075c (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Minimal structural change to support testability. Function behavior unchanged. No scope creep.

## Issues Encountered
- Spelling check diff showed new technical terms (ARD, pharmaverse, lsm, etc.) but wordlist was already current -- no action needed

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- ARD conversion is complete and ready for Phase 3 (gtsummary table generation)
- pool_to_ard() output passes cards::check_ard_structure() validation
- ARD structure includes grouping columns that gtsummary will use for table layout
- cards is in Suggests, consistent with planned approach for cardx and gtsummary in Phase 3

---
*Phase: 02-print-summary-ard-conversion*
*Completed: 2026-02-08*
