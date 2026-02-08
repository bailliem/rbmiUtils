---
phase: 03-efficacy-tables
plan: 01
subsystem: tables
tags: [gt, efficacy-table, regulatory, clinical-trial]

dependency-graph:
  requires: [01-01, 01-02]
  provides: [efficacy_table-function, gt-integration]
  affects: [03-02]

tech-stack:
  added: [gt]
  patterns: [dependency-guard, tidy-then-render]

file-tracking:
  key-files:
    created:
      - R/efficacy_table.R
      - tests/testthat/test-efficacy_table.R
      - man/efficacy_table.Rd
    modified:
      - DESCRIPTION
      - NAMESPACE
      - inst/WORDLIST

decisions:
  - id: 03-01-D1
    decision: "Added letter-digit boundary spacing in visit label cleaning (regex: ([a-zA-Z])(\\d) -> \\1 \\2)"
    rationale: "Pool object visit names like 'Week4' need space insertion, not just underscore replacement"

metrics:
  duration: ~5 min
  completed: 2026-02-08
---

# Phase 3 Plan 1: Efficacy Table Core Function Summary

**One-liner:** Regulatory-style efficacy table from rbmi pool object using gt with visit row groups, arm labels, and footnotes.

## What Was Built

Created `efficacy_table()` -- a single-call function that takes an rbmi pool object and produces a publication-ready gt table matching CDISC/ICH Table 14.2.x format. The function bridges rbmi's analysis pipeline output directly into formatted regulatory tables.

### Key Components

1. **`efficacy_table()` function** (exported): Takes a pool object, tidies it via `tidy_pool_obj()`, constructs visit-grouped gt table with LS Means and Treatment Differences per visit.

2. **`is_gt_available()` helper** (internal): Dependency guard mirroring the `is_cards_available()` pattern from `R/ard_conversion.R`.

3. **Comprehensive test suite**: 22 tests covering all function behaviors.

## Implementation Details

### Data Flow
```
pool_obj -> tidy_pool_obj() -> clean visit labels -> create row labels -> build gt table
```

### Table Structure
- **Row groups:** Visit labels (cleaned: underscore to space, letter-digit boundary spacing, title case)
- **Rows per visit:** LS Mean (Reference), LS Mean (Treatment), Treatment Difference
- **Columns:** Estimate, Std. Error, {ci_level}% CI, P-value
- **P-values:** Em dash for LSM rows (NA), formatted value for treatment rows
- **Footnotes:** Pooling method, number of imputations, confidence level

### Parameters
- `pool_obj`: Required pool object
- `title`, `subtitle`: Optional table headers
- `digits`: Decimal places for estimates (default: 2)
- `ci_level`: Confidence level for CI label (auto-extracted from pool_obj)
- `arm_labels`: Named vector `c(ref = "Placebo", alt = "Drug A")` for custom arm names
- `pval_digits`, `pval_threshold`: P-value formatting controls

## Decisions Made

### 03-01-D1: Letter-digit boundary spacing in visit labels
- **Context:** Pool object visit names like "Week4" come from parameter parsing without separators
- **Decision:** Added regex `gsub("([a-zA-Z])(\\d)", "\\1 \\2", ...)` to insert space between letters and digits
- **Rationale:** `gsub("_", " ", "Week4")` alone produces "Week4", not the desired "Week 4"

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Visit label cleaning missed letter-digit boundaries**
- **Found during:** Task 2 (test execution)
- **Issue:** Plan specified `gsub("_", " ", visit)` + `tools::toTitleCase()` but mock pool visits "Week4"/"Week8" have no underscores -- space insertion between letters and digits was needed
- **Fix:** Added `gsub("([a-zA-Z])(\\d)", "\\1 \\2", visit_clean)` before title case
- **Files modified:** R/efficacy_table.R
- **Commit:** 051bff2

**2. [Rule 3 - Blocking] Spelling wordlist needed CDISC, ICH, MMRM terms**
- **Found during:** Task 2 (R CMD check)
- **Issue:** Roxygen documentation references CDISC, ICH, MMRM acronyms flagged by spelling checker
- **Fix:** Ran `spelling::update_wordlist()` to add these terms
- **Files modified:** inst/WORDLIST
- **Commit:** 051bff2

## Verification Results

| Check | Result |
|-------|--------|
| `devtools::document()` | Pass |
| `devtools::test(filter = 'efficacy_table')` | 22/22 pass |
| `devtools::check(args = c('--no-vignettes'))` | 0 errors, 0 warnings, 2 notes (pre-existing) |
| gt in DESCRIPTION Suggests | Confirmed |
| efficacy_table in NAMESPACE | Confirmed |
| efficacy_table(mock_pool) returns gt_tbl | Confirmed |

## Test Coverage

| Test Area | Tests | Status |
|-----------|-------|--------|
| Dependency guard | 1 | Pass |
| Input validation | 2 | Pass |
| Returns gt_tbl | 1 | Pass |
| Visit labels | 2 | Pass |
| Row labels | 2 | Pass |
| Custom arm labels | 2 | Pass |
| Default arm labels | 2 | Pass |
| Footnotes | 3 | Pass |
| Title/subtitle | 2 | Pass |
| P-value formatting | 2 | Pass |
| Digits parameter | 1 | Pass |
| CI level from pool | 1 | Pass |
| CI level override | 1 | Pass |

## Next Phase Readiness

Plan 03-02 can proceed. The `efficacy_table()` function is fully operational and tested. Phase 3 Plan 2 (if it covers customization/theming or additional table variants) can build on this foundation.
