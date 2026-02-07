# Coding Conventions

**Analysis Date:** 2026-02-07

## Naming Patterns

**Files:**
- snake_case with `.R` extension
- Example: `data_helpers.R`, `analyse_mi_data.R`, `test-formatting.R`
- Test files prefixed with `test-` in `tests/testthat/` directory

**Functions:**
- snake_case for exported and internal functions
- Examples: `validate_data()`, `prepare_data_ice()`, `format_pvalue()`, `get_imputed_data()`, `gcomp_responder()`
- Internal helper functions also snake_case: `extract_covariates2()`, `as_simple_formula2()`

**Variables:**
- snake_case for local variables and parameters
- Examples: `imputed_dfs`, `test_data`, `missing_cols`, `id_map`
- Column names follow clinical data conventions (UPPERCASE for CDISC-like vars): `USUBJID`, `AVISIT`, `TRT`, `CHG`, `BASE`
- Temporary working variables use snake_case: `visit_order`, `is_ice`, `n_unique`

**Types/Classes:**
- No explicit S3/S4 class definitions; functions return standard data.frame/tibble/list structures
- Output objects documented in `@return` tags with class hints (e.g., "tibble", "data.frame", "list")

## Code Style

**Formatting:**
- Roxygen2 documentation system for all exported functions (roxygen2 7.3.2)
- 2-space indentation (standard in R packages)
- Function bodies prefer tidyverse style with pipe operator `|>` (modern R 4.1+)

**Linting:**
- No explicit linting configuration detected (.lintr, eslintrc absent)
- Implicit standard: follow tidyverse style guide
- Code adheres to clean, readable patterns without aggressive abbreviations

**Indentation Example:**
```r
result <- data |>
  dplyr::group_by(visit = .data[[visit_var]]) |>
  dplyr::summarise(
    n = dplyr::n(),
    n_miss = sum(is.na(.data[[outcome_var]])),
    .groups = "drop"
  )
```

## Import Organization

**Namespace Use:**
- Explicit namespace references: `dplyr::`, `purrr::`, `rbmi::`, `beeca::`, `assertthat::`
- Example from `data_helpers.R`:
```r
result <- data[is_ice, , drop = FALSE]
ice_data <- ice_data |>
  dplyr::arrange(.data[[subjid_var]], .data$.visit_order) |>
  dplyr::distinct(.data[[subjid_var]], .keep_all = TRUE)
```
- Dependency order in DESCRIPTION: assertthat, dplyr, purrr, rbmi, beeca, rlang, tidyr
- No `library()` calls in function code; imports declared in DESCRIPTION
- Examples in documentation use explicit `library()` calls for clarity

**Path Aliases:**
- No custom path aliases; direct relative imports via namespace prefixes

## Error Handling

**Primary Strategy:**
- Defensive validation using `assertthat::assert_that()` for preconditions
- `stop()` with `call. = FALSE` for user-friendly error messages
- Multi-issue error collection before reporting (seen in `validate_data()`)

**Pattern:**
```r
issues <- character(0)

# Accumulate all issues
if (!is.data.frame(data)) {
  issues <- c(issues, "`data` must be a data.frame")
}
if (!is.numeric(outcome_var)) {
  issues <- c(issues, "`outcome` must be numeric")
}

# Report all together
if (length(issues) > 0) {
  stop(paste(c("Validation failed:", issues), collapse = "\n  - "),
       call. = FALSE)
}
```

**assertthat Usage:**
```r
assertthat::assert_that(
  is.data.frame(data),
  msg = "`data` must be a data.frame"
)
assertthat::assert_that(
  is.numeric(digits) && length(digits) == 1 && digits >= 1,
  msg = "`digits` must be a positive integer"
)
```

**Warnings:**
- Used for non-fatal issues: character-to-factor conversions, invalid p-values
- Pattern: `warning(sprintf(...), call. = FALSE)`

## Logging

**Framework:** Base R `message()` and `warning()`

**Patterns:**
- Validation functions report all issues in a single error message
- Status messages: `message("Processing imputation ", i)` (implicit, rarely used)
- Diagnostic information: provided via return values and attributes, not logging
- No explicit logging framework (no logger, glue, or stringr logging)

## Comments

**When to Comment:**
- Function-level documentation: Roxygen2 `#'` tags (comprehensive)
- Complex algorithmic steps: inline comments for non-obvious logic
- Example from `data_helpers.R`:
```r
# For each subject, take the first visit by factor level order
if (is.factor(ice_data[[visit_var]])) {
  visit_order <- as.integer(ice_data[[visit_var]])
} else {
  visit_order <- seq_len(nrow(ice_data))
}
```

**Roxygen/Documentation:**
- Mandatory for all exported functions (`@export` tag)
- Includes: `@param`, `@return`, `@details`, `@examples`, `@seealso`, `@export`
- Example from `format_pvalue()`:
```r
#' Format P-values for Publication
#'
#' Formats p-values according to common publication standards, with configurable
#' thresholds and decimal places.
#'
#' @param x A numeric vector of p-values.
#' @param digits Integer. Number of decimal places for rounding. Default is 3.
#' @param threshold Numeric. P-values below this threshold are displayed as
#'   "< threshold". Default is 0.001.
#'
#' @return A character vector of formatted p-values.
#'
#' @details
#' The function applies the following rules:
#' \itemize{...}
#' @examples
#' \donttest{...}
#'
#' @export
```

## Function Design

**Size:**
- Prefer functions under 50 lines for simple tasks
- Longer functions (100+ lines) for complex multi-step processes like `validate_data()`, `summarise_missingness()`
- Helper functions extracted for reusability

**Parameters:**
- Named parameters with defaults where sensible
- Example: `format_pvalue(x, digits = 3, threshold = 0.001, html = FALSE)`
- Validation of all parameters at function start using `assertthat::assert_that()`

**Return Values:**
- Single return per function (either implicit last value or explicit `return()`)
- Multiple values returned as list or data.frame
- Invisible returns for side-effect functions: `invisible(TRUE)` in `validate_data()`
- Examples: returns tibble from `summarise_missingness()`, named list from `gcomp_responder()`

**Example from `analyse_mi_data.R`:**
```r
analyse_mi_data <- function(
  data = NULL,
  vars,
  method,
  fun = rbmi::ancova,
  delta = NULL,
  ...
) {
  # Validation
  assertthat::assert_that(is.data.frame(data), msg = "...")
  assertthat::assert_that(is.list(vars), msg = "...")

  # Processing
  # ...

  # Return analysis object
  out
}
```

## Module Design

**Exports:**
- Single function per file (rare) or grouped related functions
- Example: `formatting.R` exports `format_pvalue()`, `format_estimate()`, `format_results_table()`
- All exported functions marked with `@export` in Roxygen

**Barrel Files:**
- Not used; no index files like `index.R` or `__init__.R`
- NAMESPACE file auto-generated by Roxygen2 from `@export` tags
- Users import via `library(rbmiUtils)` and access all exported functions

## Code Example: Typical Pattern

From `data_helpers.R`, function `validate_data()`:

```r
#' Validate Data Before Imputation
#'
#' Pre-flight validation of data, variable specification, and intercurrent event
#' data before calling [rbmi::draws()]...
#'
#' @param data A data.frame containing the analysis dataset.
#' @param vars A `vars` object as created by [rbmi::set_vars()].
#' @param data_ice An optional data.frame...
#'
#' @return Invisibly returns `TRUE` if all checks pass...
#' @export
validate_data <- function(data, vars, data_ice = NULL) {

  issues <- character(0)

  # --- Basic type checks ---
  if (!is.data.frame(data)) {
    issues <- c(issues, "`data` must be a data.frame")
  }

  assertthat::assert_that(
    is.list(vars),
    msg = "`vars` must be a list as returned by `rbmi::set_vars()`"
  )

  # ... more validation ...

  # --- Report results ---
  if (length(issues) > 0) {
    stop(paste(c("Data validation failed:", issues), collapse = "\n  - "),
         call. = FALSE)
  }

  invisible(TRUE)
}
```

---

*Convention analysis: 2026-02-07*
