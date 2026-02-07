---
phase: 01-foundation-hardening
plan: 01
subsystem: tidiers
tags: [parameter-parsing, regex, cli, lifecycle, error-handling]

dependency_graph:
  requires: []
  provides:
    - "cli and lifecycle available as Imports for all subsequent hardening plans"
    - "Robust regex-based parameter parsing in tidy_pool_obj()"
    - "Custom error class hierarchy (rbmiUtils_error_validation)"
  affects:
    - "01-02 (analyse_mi_data refactor can use cli::cli_abort)"
    - "01-03 (gcomp hardening can use cli::cli_abort)"
    - "02-xx (ARD conversion depends on stable tidy_pool_obj output)"

tech_stack:
  added:
    - "cli >= 3.6.0 (error messaging with inline formatting)"
    - "lifecycle >= 1.0.4 (deprecation support for future plans)"
  patterns:
    - "tidyr::separate_wider_regex() with named capture groups for column extraction"
    - "cli::cli_abort() with custom error class hierarchy"
    - "Two-pass parsing: regex extraction then conditional logic"

file_tracking:
  created: []
  modified:
    - "DESCRIPTION"
    - "R/tidiers.R"
    - "R/rbmiUtils-package.R"
    - "tests/testthat/test-tidiers.R"
    - "NEWS.md"

decisions:
  - id: "01-01-D1"
    decision: "Use two-pass parsing (separate_wider_regex then dplyr::mutate) rather than single regex"
    reason: "Parameter name formats vary between ANCOVA (trt_Visit) and gcomp (lsm_ArmName_Visit). A single regex cannot cleanly capture both without backtracking."
  - id: "01-01-D2"
    decision: "Treatment comparison description does not include visit name"
    reason: "Maintains backward compatibility with existing test expectations. Treatment rows already have visit in the visit column."
  - id: "01-01-D3"
    decision: "Mock pool objects must use rbmi internal structure ($pars list) not flat data.frames"
    reason: "rbmi's as.data.frame.pool() method requires $pars list with $est, $se, $ci, $pvalue elements. Flat data.frame mock would fail as_tibble conversion."

metrics:
  duration: "~9 minutes"
  completed: "2026-02-07"
---

# Phase 01 Plan 01: Fix tidy_pool_obj Parameter Parsing Summary

**Add cli/lifecycle dependencies and replace fragile underscore-based parameter splitting with regex-based parsing using tidyr::separate_wider_regex()**

## Performance

- **Started:** 2026-02-07T21:21:42Z
- **Completed:** 2026-02-07T21:30:55Z
- **Duration:** ~9 minutes

## Accomplishments

1. **Added cli and lifecycle to DESCRIPTION Imports** -- cli (>= 3.6.0) and lifecycle (>= 1.0.4) added in alphabetical order, enabling all subsequent hardening plans to use `cli::cli_abort()` and `lifecycle::deprecate_warn()`.

2. **Replaced fragile parameter parsing** -- `tidyr::separate(sep = "_")` replaced with `tidyr::separate_wider_regex()` using a two-pass approach: first extract `parameter_type` (trt|lsm) and remainder, then parse remainder conditionally based on parameter type. This correctly handles parameter names containing underscores (e.g., "Week_24", "Follow_Up") and gcomp-style arm names (e.g., "lsm_Placebo_Week 24").

3. **Upgraded error handling** -- `stop()` replaced with `cli::cli_abort()` using custom error class `c("rbmiUtils_error_validation", "rbmiUtils_error")` for programmatic error catching.

4. **Added comprehensive tests** -- New tests verify underscore-containing visit names parse correctly, gcomp-style parameter names with arm names as lsm_type work, and error class hierarchy is catchable. Updated existing error tests to use class matching instead of string matching.

5. **Documented breaking change in NEWS.md** -- Breaking change section prepended to development version documenting the parameter parsing change and new dependencies.

## Task Commits

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Add cli/lifecycle deps and update tidy_pool_obj parsing | `96986a9` | DESCRIPTION, R/tidiers.R |
| 2 | Update tests and add NEWS.md entry | `68615eb` | tests/testthat/test-tidiers.R, NEWS.md, R/rbmiUtils-package.R |

## Files Modified

- **DESCRIPTION** -- Added cli (>= 3.6.0) and lifecycle (>= 1.0.4) to Imports, reordered alphabetically
- **R/tidiers.R** -- Replaced separate() with separate_wider_regex(), replaced stop() with cli::cli_abort()
- **R/rbmiUtils-package.R** -- Added 'remainder' to globalVariables to silence R CMD check NOTE
- **tests/testthat/test-tidiers.R** -- Added 3 new test blocks (underscore params, gcomp params, error class), updated 3 existing error tests
- **NEWS.md** -- Added Breaking Changes and Dependencies sections

## Decisions Made

1. **Two-pass parsing approach** (01-01-D1): Used `separate_wider_regex()` for initial type extraction, then `dplyr::mutate()` with `case_when()` for context-dependent remainder parsing. A single regex cannot handle both ANCOVA and gcomp formats cleanly.

2. **Treatment description backward compatibility** (01-01-D2): Kept `description = "Treatment Comparison"` without appending visit name, matching existing test expectations. Visit information is available in the `visit` column.

3. **Mock pool object structure** (01-01-D3): Discovered that `dplyr::as_tibble()` on pool objects dispatches through `rbmi::as.data.frame.pool()` which requires internal `$pars` list structure. Tests must use proper mock construction, not flat data.frames with class attribute.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed R CMD check NOTE for unbound global variable**
- **Found during:** Task 2 (R CMD check verification)
- **Issue:** `remainder` column created by `separate_wider_regex()` triggered "no visible binding for global variable" NOTE
- **Fix:** Added `"remainder"` to `utils::globalVariables()` in `R/rbmiUtils-package.R`
- **Files modified:** R/rbmiUtils-package.R
- **Commit:** 68615eb

**2. [Rule 1 - Bug] Mock pool objects require rbmi internal structure**
- **Found during:** Task 1 verification
- **Issue:** Simple `data.frame` with `class(x) <- "pool"` fails because `dplyr::as_tibble()` dispatches through `rbmi::as.data.frame.pool()` which requires `$pars` list
- **Fix:** Created `make_mock_pool()` helper in tests that builds proper internal structure
- **Files modified:** tests/testthat/test-tidiers.R
- **Commit:** 68615eb

## Issues Encountered

None blocking. The lifecycle import generates a NOTE ("Namespace in Imports field not imported from: 'lifecycle'") because it is not yet used in code -- this is expected and will resolve when Plan 01-02 adds `lifecycle::deprecate_warn()` calls.

## Next Phase Readiness

- **Blockers:** None
- **Ready for:** 01-02 (analyse_mi_data refactor) -- cli::cli_abort() now available for error messaging
- **Concerns:** None
