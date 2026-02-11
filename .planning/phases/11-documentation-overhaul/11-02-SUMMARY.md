---
phase: 11-documentation-overhaul
plan: 02
subsystem: documentation
tags: [roxygen2, examples, ADMI, ADEFF, rbmi-pipeline, donttest, dontrun]

# Dependency graph
requires:
  - phase: 08-mi-diagnostic-statistics
    provides: analyse_mi_data, pool_to_ard with MI diagnostics
  - phase: 09-describe-functions
    provides: describe_draws, describe_imputation functions
  - phase: 10-publication-styling
    provides: efficacy_table and plot_forest styling parameters
provides:
  - Executable @examples for pool_to_ard, efficacy_table, plot_forest using ADMI
  - Realistic @examples for describe_draws, describe_imputation using ADEFF
affects: [11-03-PLAN, vignettes, pkgdown site]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "ADMI-based examples use \\donttest{} (executable but slow)"
    - "ADEFF-based examples use \\dontrun{} (require full MCMC pipeline)"
    - "requireNamespace guards for Suggests dependencies"
    - "n_samples=20 with warmup=20, thin=1 for minimal MCMC in examples"

key-files:
  created: []
  modified:
    - R/ard_conversion.R
    - R/efficacy_table.R
    - R/plot_forest.R
    - R/describe.R
    - man/pool_to_ard.Rd
    - man/efficacy_table.Rd
    - man/plot_forest.Rd
    - man/describe_draws.Rd
    - man/print.describe_draws.Rd
    - man/describe_imputation.Rd
    - man/print.describe_imputation.Rd

key-decisions:
  - "11-02-D1: ADMI examples use \\donttest{} since they CAN execute with pre-imputed data"
  - "11-02-D2: ADEFF examples use \\dontrun{} since they require full MCMC computation"
  - "11-02-D3: print.* examples use short references to parent function instead of repeating full pipeline"

patterns-established:
  - "Example pattern: data() -> factor prep -> set_vars -> method -> analyse_mi_data/draws -> pool/impute -> function"
  - "requireNamespace guards for Suggests packages in examples"

# Metrics
duration: 2min
completed: 2026-02-11
---

# Phase 11 Plan 02: Function Examples Upgrade Summary

**Replaced commented-out pseudocode in 5 exported functions with executable/realistic examples using ADMI and ADEFF package datasets**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-11T19:43:44Z
- **Completed:** 2026-02-11T19:45:55Z
- **Tasks:** 2
- **Files modified:** 11 (4 R source + 7 Rd man pages)

## Accomplishments
- pool_to_ard(), efficacy_table(), plot_forest() examples now use executable code with ADMI data inside \donttest{}
- describe_draws() and describe_imputation() examples show full ADEFF pipeline from data prep to programmatic metadata access
- print.describe_draws() and print.describe_imputation() examples are self-contained with no undefined variables
- All roxygen2 processing completes without errors

## Task Commits

Each task was committed atomically:

1. **Task 1: Upgrade pool_to_ard, efficacy_table, plot_forest examples** - `76a0471` (docs)
2. **Task 2: Upgrade describe_draws and describe_imputation examples** - `957a6d8` (docs)

## Files Created/Modified
- `R/ard_conversion.R` - Executable pool_to_ard() example with ADMI pipeline
- `R/efficacy_table.R` - Executable efficacy_table() example with ADMI pipeline and publication styling demo
- `R/plot_forest.R` - Executable plot_forest() example with ADMI pipeline showing trt and lsm modes
- `R/describe.R` - Realistic describe_draws/describe_imputation examples with ADEFF pipeline
- `man/pool_to_ard.Rd` - Regenerated from roxygen2
- `man/efficacy_table.Rd` - Regenerated from roxygen2
- `man/plot_forest.Rd` - Regenerated from roxygen2
- `man/describe_draws.Rd` - Regenerated from roxygen2
- `man/print.describe_draws.Rd` - Regenerated from roxygen2
- `man/describe_imputation.Rd` - Regenerated from roxygen2
- `man/print.describe_imputation.Rd` - Regenerated from roxygen2

## Decisions Made
- **11-02-D1:** ADMI-based examples (pool_to_ard, efficacy_table, plot_forest) use `\donttest{}` since they can execute with pre-imputed data -- just slow due to MCMC
- **11-02-D2:** ADEFF-based examples (describe_draws, describe_imputation) use `\dontrun{}` since they require full rbmi draws/impute pipeline with MCMC computation
- **11-02-D3:** print.describe_draws() and print.describe_imputation() examples use short references to parent function examples rather than repeating the full pipeline, avoiding undefined variables

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All 5 key exported function examples are now upgraded
- Function help pages demonstrate real clinical analysis patterns
- Ready for Plan 03 (remaining documentation overhaul tasks)

---
*Phase: 11-documentation-overhaul*
*Completed: 2026-02-11*
