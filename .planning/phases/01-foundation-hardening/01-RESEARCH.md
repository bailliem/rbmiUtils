# Phase 1: Foundation Hardening - Research

**Researched:** 2026-02-07
**Domain:** R package hardening -- input validation, error messaging, deprecation, round-trip data integrity
**Confidence:** HIGH

## Summary

This phase hardens four function groups in rbmiUtils: (1) `tidy_pool_obj()` parameter parsing, (2) `analyse_mi_data()` refactor to delegate to `rbmi::analyse()`, (3) `gcomp_responder()`/`gcomp_binary()` input validation and beeca output pinning, and (4) `reduce_imputed_data()`/`expand_imputed_data()` round-trip fidelity with digest verification. The locked decisions specify using `cli` for error formatting with custom error classes, `lifecycle` for deprecation warnings, fail-early validation (stop at first error), and exact `identical()` round-trip reproduction with a built-in digest integrity check.

The standard approach uses three new dependencies (`cli`, `lifecycle`, and `digest` -- or alternatively `rlang::hash()` which is already available) combined with structured regex parsing via `tidyr::separate_wider_regex()` to replace the fragile `_`-separator parsing. The `analyse_mi_data()` refactor cannot simply delegate to `rbmi::analyse()` because that function requires an `imputations` object from `rbmi::impute()`, while `analyse_mi_data()` operates on data frames with an IMPID column. The refactor must instead use `rbmi::analyse()`'s internal `as_analysis()` constructor or replicate its class-building logic using `inherits()` checks rather than positional `class(method)[[2]]` indexing.

**Primary recommendation:** Use `cli::cli_abort()` with custom classes for all errors, `tidyr::separate_wider_regex()` for parameter parsing, `inherits()` for method type detection, beeca output schema validation, and `rlang::hash()` for round-trip digest verification.

## Standard Stack

The established libraries/tools for this phase:

### Core (New Dependencies)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| cli | >= 3.6.0 | `cli_abort()`, `cli_warn()` for formatted error messages with custom classes | Tidyverse/pharmaverse standard for user-facing messages. Supports inline markup, pluralization, bullet formatting. `class` parameter passed through to `rlang::abort()`. |
| lifecycle | >= 1.0.4 | `deprecate_warn()` for analyse_mi_data internals; `badge()` for roxygen2 | Standard R package deprecation workflow. Badge renders in HTML docs. No need to add as Import for badge-only use, but required as Import for `deprecate_warn()` at runtime. |

### Core (Already Imported)

| Library | Version | Purpose | Why Relevant |
|---------|---------|---------|--------------|
| rlang | >= 1.1.0 | `rlang::hash()` for round-trip digest, `rlang::abort()` underlies `cli_abort()`, `.data` pronoun | Already an Import. `hash()` uses XXH128 algorithm, reproducible within same R version and platform. No need for separate `digest` dependency. |
| tidyr | >= 1.3.0 | `separate_wider_regex()` replaces superseded `separate()` | Already an Import. Version >= 1.3.0 required for `separate_wider_regex()`. Verify current version constraint in DESCRIPTION. |
| assertthat | current | Existing validation (will be gradually replaced by `cli_abort()` in new code) | Already an Import. Existing code uses it. New validation code should use `cli_abort()` instead. |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| beeca | >= 0.1.3 | G-computation via `get_marginal_effect()` | Already imported. Add version constraint `>= 0.1.3` to pin output format. |
| rbmi | >= 1.4 | `as_class()`, `pool()`, `ancova()`, method constructors | Already imported. `as_analysis()` is internal to rbmi (not exported). Must replicate class-building logic. |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| rlang::hash() | digest::digest() | digest supports more algorithms (SHA-256, MD5) but adds a dependency. rlang::hash() uses XXH128 which is fast and sufficient for integrity checks. Prefer rlang since already imported. |
| cli::cli_abort() | rlang::abort() directly | cli_abort adds formatted output (bullets, markup). rlang::abort() is lower-level. Decision locks cli formatting, so use cli_abort(). |
| tidyr::separate_wider_regex() | stringr::str_match() with manual extraction | Both work. separate_wider_regex() integrates into tidyverse pipeline and produces tidy columns directly. str_match requires manual column assignment. Prefer tidyr. |

**Installation (additions to DESCRIPTION Imports):**
```
cli (>= 3.6.0),
lifecycle (>= 1.0.4)
```

**Verify existing version constraints:**
```
tidyr (>= 1.3.0)   # Currently unversioned, needs >= 1.3.0 for separate_wider_regex
beeca (>= 0.1.3)   # Currently unversioned, needs version pin
```

## Architecture Patterns

### Current vs Target Error Handling Pattern

```
Current (assertthat + stop):
  assertthat::assert_that(is.data.frame(data), msg = "`data` must be a data.frame")
  stop("`data` cannot be NULL.", call. = FALSE)

Target (cli):
  cli::cli_abort(
    c("x" = "{.arg data} must be a data frame, not {.cls {class(data)}}.",
      "i" = "Pass a data frame with columns {.field IMPID}, {.field USUBJID}, etc."),
    class = "rbmiUtils_error_validation",
    call = rlang::caller_env()
  )
```

### Custom Error Class Hierarchy

```
rbmiUtils_error                    # Base class for all rbmiUtils errors
  rbmiUtils_error_validation       # Input validation failures (wrong type, missing columns)
  rbmiUtils_error_type             # Type mismatch errors (wrong class of object)
  rbmiUtils_error_integrity        # Round-trip integrity check failures
  rbmiUtils_error_dependency       # Upstream package (beeca, rbmi) output format violations
```

Each error is created via `cli::cli_abort()` with the `class` parameter:
```r
cli::cli_abort(
  message = c("Error headline", "x" = "Detail", "i" = "Hint"),
  class = c("rbmiUtils_error_validation", "rbmiUtils_error"),
  call = rlang::caller_env()
)
```

### Pattern 1: Fail-Early Validation Block

**What:** All input validation at function start, stopping at first error
**When to use:** Every public function entry point

```r
# Source: Locked decision -- stop at first validation error
my_function <- function(data, vars, method) {
  # --- Input validation ---
  if (!is.data.frame(data)) {
    cli::cli_abort(
      c("{.arg data} must be a data frame.",
        "x" = "You supplied a {.cls {class(data)}} object."),
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  if (!"IMPID" %in% names(data)) {
    cli::cli_abort(
      c("{.arg data} must contain an {.field IMPID} column.",
        "i" = "Use {.fn get_imputed_data} or {.fn expand_imputed_data} to create data with IMPID."),
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  # --- Business logic starts after all validation ---
  ...
}
```

### Pattern 2: Deprecation Cycle for analyse_mi_data Internals

**What:** Old internal helpers (`extract_covariates2`, `as_simple_formula2`, `as_analysis2`) deprecated with one-version warning
**When to use:** analyse_mi_data() refactor

```r
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This internal function has been deprecated in favour of direct
#' delegation to rbmi's class construction. It will be removed in
#' the next minor version.
#'
#' @keywords internal
as_analysis2 <- function(results, method, delta = NULL, fun = NULL, fun_name = NULL) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "as_analysis2()",
    details = "analyse_mi_data() now constructs the analysis object directly."
  )
  # ... existing implementation for backward compat ...
}
```

### Pattern 3: Regex-Based Parameter Parsing

**What:** Replace `tidyr::separate()` with `tidyr::separate_wider_regex()` using named capture groups
**When to use:** tidy_pool_obj() fix

```r
# Source: tidyr docs -- https://tidyr.tidyverse.org/reference/separate_wider_delim.html
# Current parameter format: "trt_Week 24", "lsm_ref_Week 24", "lsm_alt_Week 48"
df <- df |>
  tidyr::separate_wider_regex(
    parameter,
    patterns = c(
      parameter_type = "trt|lsm",  # Known prefix
      "_",                          # Literal separator
      lsm_type = "(?:ref|alt)",    # Optional ref/alt (only for lsm)
      "_?",                         # Optional separator
      visit = ".*"                  # Everything else is the visit
    ),
    too_few = "align_start",
    cols_remove = FALSE
  )
```

**Note:** The exact regex pattern must be validated against actual parameter names produced by `rbmi::ancova()` and `gcomp_responder()`. The parameter names follow this structure:
- Treatment comparison: `"trt_<visit>"` (e.g., `"trt_Week 24"`)
- LSM reference: `"lsm_ref_<visit>"` (e.g., `"lsm_ref_Week 24"`)
- LSM alternative: `"lsm_alt_<visit>"` (e.g., `"lsm_alt_Week 24"`)

A simpler approach using `stringr::str_match()`:
```r
parsed <- stringr::str_match(
  df$parameter,
  "^(trt|lsm)_(?:(ref|alt)_)?(.+)$"
)
df$parameter_type <- parsed[, 2]
df$lsm_type <- parsed[, 3]
df$visit <- parsed[, 4]
```

### Pattern 4: beeca Output Schema Validation

**What:** Validate beeca output structure before accessing columns
**When to use:** gcomp_responder(), gcomp_binary()

```r
marginal_fit <- beeca::get_marginal_effect(model, ...)
res <- marginal_fit$marginal_results

# Validate beeca output schema
required_cols <- c("STAT", "STATVAL", "TRTVAL")
missing_cols <- setdiff(required_cols, names(res))
if (length(missing_cols) > 0) {
  cli::cli_abort(
    c("Unexpected output from {.fn beeca::get_marginal_effect}.",
      "x" = "Missing column{?s}: {.field {missing_cols}}.",
      "i" = "This may indicate an incompatible version of {.pkg beeca}.",
      "i" = "rbmiUtils requires {.pkg beeca} >= 0.1.3."),
    class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
  )
}
```

### Pattern 5: Round-Trip Digest Verification

**What:** Store hash on reduce, verify on expand
**When to use:** reduce_imputed_data() / expand_imputed_data()

```r
# In reduce_imputed_data():
reduced <- imputed_data[keep_rows, , drop = FALSE]
# Store digest as attribute
attr(reduced, "rbmiUtils_digest") <- rlang::hash(original_data)
attr(reduced, "rbmiUtils_original_ncol") <- ncol(original_data)
attr(reduced, "rbmiUtils_original_colnames") <- names(original_data)

# In expand_imputed_data():
stored_digest <- attr(reduced_data, "rbmiUtils_digest")
if (!is.null(stored_digest)) {
  current_digest <- rlang::hash(original_data)
  if (stored_digest != current_digest) {
    cli::cli_abort(
      c("Integrity check failed: {.arg original_data} has changed since reduction.",
        "x" = "Stored digest: {.val {stored_digest}}",
        "x" = "Current digest: {.val {current_digest}}",
        "i" = "The original data passed to {.fn expand_imputed_data} must be identical to what was passed to {.fn reduce_imputed_data}."),
      class = c("rbmiUtils_error_integrity", "rbmiUtils_error")
    )
  }
}
```

### Pattern 6: Method Type Detection via inherits()

**What:** Replace `class(method)[[2]]` with `inherits()` checks
**When to use:** analyse_mi_data() and as_analysis2() refactor

```r
# Current (fragile):
next_class <- switch(class(method)[[2]], ...)

# Target (robust):
determine_pooling_class <- function(method) {
  if (inherits(method, "bayes") || inherits(method, "approxbayes")) {
    return("rubin")
  }
  if (inherits(method, "condmean")) {
    return(if (method$type == "jackknife") "jackknife" else "bootstrap")
  }
  if (inherits(method, "bmlmi")) {
    return("bmlmi")
  }
  cli::cli_abort(
    c("Unrecognized method type: {.cls {class(method)}}.",
      "i" = "Expected a method created by {.fn rbmi::method_bayes}, {.fn rbmi::method_approxbayes}, {.fn rbmi::method_condmean}, or {.fn rbmi::method_bmlmi}."),
    class = c("rbmiUtils_error_type", "rbmiUtils_error")
  )
}
```

### Anti-Patterns to Avoid

- **Positional class indexing:** Never use `class(x)[[2]]` or `class(x)[[1]]` for method dispatch. Use `inherits()`.
- **Silent type coercion:** Never silently convert types (e.g., factor to character) without validation. Check types match, error if not.
- **Batching errors when decision says stop-at-first:** The locked decision is stop at first validation error. Do not collect multiple errors (except in `validate_data()` which has different existing behavior that is out of scope for this phase).
- **Using `stop()` with `call. = FALSE` in new code:** Always use `cli::cli_abort()` for new error paths. Existing `stop()` calls in functions not being modified stay as-is.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Object hashing for integrity | Custom serialization + MD5 | `rlang::hash()` | Uses XXH128, handles arbitrary R objects including data frames with attributes. Already a dependency. |
| Deprecation warnings | Manual `warning("Deprecated...")` | `lifecycle::deprecate_warn()` | Tracks first-use, integrates with roxygen2 badges, follows R ecosystem conventions. Users expect this pattern. |
| Formatted error messages | `sprintf()` + `stop()` | `cli::cli_abort()` with named bullets | Handles pluralization, inline styling ({.arg}, {.cls}, {.fn}), glue interpolation. Consistent with tidyverse/pharmaverse ecosystem. |
| Regex column extraction | `strsplit()` + manual assignment | `tidyr::separate_wider_regex()` or `stringr::str_match()` | Named capture groups, handles edge cases (too few matches), integrates with data frame pipeline. |
| Custom error classes | `structure(list(), class = "error")` | `cli_abort(..., class = "pkg_error_type")` | Inherits from rlang condition classes. Catchable via `tryCatch(expr, pkg_error_type = handler)`. |

**Key insight:** Every problem in this phase has a well-tested existing solution in the tidyverse/rlang ecosystem. The primary risk is not the implementation complexity but rather the careful mapping of existing fragile patterns to their robust replacements without changing external behavior.

## Common Pitfalls

### Pitfall 1: analyse_mi_data Cannot Simply Delegate to rbmi::analyse()

**What goes wrong:** Attempting to refactor `analyse_mi_data()` to call `rbmi::analyse()` directly fails because `rbmi::analyse()` requires an `imputations` object from `rbmi::impute()`, not a data frame with an IMPID column. The function signatures are fundamentally different: rbmi::analyse(imputations, fun, delta) vs analyse_mi_data(data, vars, method, fun, delta).

**Why it happens:** The requirement says "wrap rbmi::analyse()" but the actual API contract does not permit this. `rbmi::analyse()` internally calls `extract_imputed_df()` on an imputations object, applies delta, then runs the analysis function. `analyse_mi_data()` already has the extracted data frames (from `get_imputed_data()` or `expand_imputed_data()`).

**How to avoid:** The refactor should:
1. Replace `as_analysis2()` (which copies rbmi's `as_analysis()` logic) with `inherits()`-based class detection
2. Eliminate the copied helper functions (`extract_covariates2()`, `as_simple_formula2()`)
3. Use `rbmi::as_class()` for result subclass assignment (this IS exported by rbmi)
4. Keep the data-frame-with-IMPID interface since it serves users who store/retrieve imputed data

**Warning signs:** Tests that try to pass a data frame to `rbmi::analyse()` will fail with "Input must be an imputations object."

### Pitfall 2: separate_wider_regex() Requires Exact Pattern Match

**What goes wrong:** `separate_wider_regex()` errors when the regex does not match all rows. Unlike `separate()` which has `fill = "right"`, `separate_wider_regex()` uses `too_few = "error"` by default.

**Why it happens:** The parameter names from `gcomp_responder()` follow a different format than `rbmi::ancova()`. ANCOVA produces `"trt_Week 24"`, `"lsm_ref_Week 24"`. G-computation produces `"trt_Drug A vs Placebo_Week 24"`, `"lsm_Placebo_Week 24"`. A single regex must handle both formats.

**How to avoid:** Use `too_few = "align_start"` to handle partial matches gracefully. Test the regex against actual parameter names from both ANCOVA and g-computation outputs before finalizing. Consider using `stringr::str_match()` instead if the two-format complexity is too high for a single regex.

**Warning signs:** `separate_wider_regex()` errors with "Expected X pieces, got Y" or produces NA columns.

### Pitfall 3: rlang::hash() Is R-Version-Dependent

**What goes wrong:** A digest stored with R 4.3 may not match when verified with R 4.4, because `rlang::hash()` uses R's internal serialization (version 3), which can change between R versions.

**Why it happens:** The hash depends on the exact binary serialization of the R object. Different R versions may serialize the same object differently.

**How to avoid:** Document that the digest is valid only within the same R session or R version. If cross-version compatibility is required, use `digest::digest(algo = "xxhash128", serialize = TRUE, serializeVersion = 3)` which explicitly pins the serialization version. For this phase, `rlang::hash()` is sufficient since reduce/expand typically happens within the same R session.

**Warning signs:** Integrity check failures when users upgrade R between reduce and expand operations.

### Pitfall 4: cli::cli_abort() Class Must Include Parent Classes

**What goes wrong:** Creating an error with `class = "rbmiUtils_error_validation"` but not including `"rbmiUtils_error"` in the class vector means users cannot catch all rbmiUtils errors with a single handler.

**Why it happens:** R condition handling dispatches on the first matching class. If a user writes `tryCatch(expr, rbmiUtils_error = handler)`, it only matches if `"rbmiUtils_error"` is in the condition's class vector.

**How to avoid:** Always pass the full class hierarchy: `class = c("rbmiUtils_error_validation", "rbmiUtils_error")`. The more specific class comes first.

**Warning signs:** `tryCatch(expr, rbmiUtils_error = handler)` does not catch validation errors.

### Pitfall 5: Breaking Existing Test Expectations

**What goes wrong:** Existing tests match specific error messages (e.g., `expect_error(expr, "must be of class 'pool'")`) using regex. Changing to `cli_abort()` changes the exact message text, breaking tests.

**Why it happens:** cli formatting adds markup characters and changes message structure (bullets, inline formatting).

**How to avoid:** Update tests to match on error class instead of message text: `expect_error(expr, class = "rbmiUtils_error_validation")`. For backward compatibility during transition, also check that the message contains key words. Use `conditionMessage()` for message text assertions.

**Warning signs:** Existing tests fail with "Error message does not match" after switching to cli_abort.

### Pitfall 6: Attribute Preservation in Subsetting Operations

**What goes wrong:** R drops custom attributes when subsetting data frames with `[`. After `reduced <- imputed_data[keep_rows, , drop = FALSE]`, custom column attributes (labels, formats) may be lost.

**Why it happens:** Base R's `[.data.frame` preserves some attributes but not all custom ones. The current code manually restores attributes in a loop, but the exclude list (`c("names", "dim", "dimnames")`) may not cover all framework-specific attributes.

**How to avoid:** The current attribute restoration loop is correct in principle. Ensure the exclude list is minimal (only `names`, `dim`, `dimnames`). Add tests that verify specific attributes (e.g., `label` attribute on columns) survive the round-trip.

**Warning signs:** `identical(original, roundtripped)` fails due to missing attributes.

## Code Examples

Verified patterns from official sources:

### cli_abort with Custom Class and Formatted Message

```r
# Source: https://cli.r-lib.org/reference/cli_abort.html + https://rlang.r-lib.org/reference/abort.html
validate_pool_input <- function(pool_obj) {
  if (!inherits(pool_obj, "pool")) {
    cli::cli_abort(
      c("{.arg pool_obj} must be a {.cls pool} object.",
        "x" = "You supplied a {.cls {class(pool_obj)}} object.",
        "i" = "Create a pool object with {.code rbmi::pool(analysis_obj)}."),
      class = c("rbmiUtils_error_type", "rbmiUtils_error")
    )
  }
}

# Catching:
tryCatch(
  validate_pool_input(mtcars),
  rbmiUtils_error_type = function(e) message("Type error: ", conditionMessage(e)),
  rbmiUtils_error = function(e) message("Any rbmiUtils error: ", conditionMessage(e))
)
```

### lifecycle::deprecate_warn with Roxygen Badge

```r
# Source: https://lifecycle.r-lib.org/articles/communicate.html

#' Construct an rbmi analysis object
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `as_analysis2()` was deprecated in rbmiUtils 0.2.0.
#' `analyse_mi_data()` now constructs the analysis object directly
#' using `inherits()`-based method detection.
#'
#' @keywords internal
as_analysis2 <- function(results, method, delta = NULL, fun = NULL, fun_name = NULL) {
  lifecycle::deprecate_warn("0.2.0", "as_analysis2()")
  # ... implementation ...
}
```

### tidyr::separate_wider_regex for Parameter Parsing

```r
# Source: https://tidyr.tidyverse.org/reference/separate_wider_delim.html
# Handles: "trt_Week 24", "lsm_ref_Week 24", "lsm_alt_Week 48"
df |>
  tidyr::separate_wider_regex(
    parameter,
    patterns = c(
      parameter_type = "trt|lsm",
      "_",
      lsm_type = "(?:ref|alt)?",
      "_?",
      visit = ".+"
    ),
    too_few = "align_start",
    cols_remove = FALSE
  )
```

### rlang::hash for Data Frame Digest

```r
# Source: https://rlang.r-lib.org/reference/hash.html
# Store digest on reduce
digest_value <- rlang::hash(original_data)
attr(reduced, "rbmiUtils_digest") <- digest_value

# Verify on expand
stored <- attr(reduced_data, "rbmiUtils_digest")
current <- rlang::hash(original_data)
if (!is.null(stored) && stored != current) {
  cli::cli_abort(
    c("Round-trip integrity check failed.",
      "x" = "The {.arg original_data} has changed since {.fn reduce_imputed_data} was called.",
      "i" = "Stored digest: {.val {stored}}",
      "i" = "Current digest: {.val {current}}"),
    class = c("rbmiUtils_error_integrity", "rbmiUtils_error")
  )
}
```

### beeca Output Validation

```r
# Source: https://search.r-project.org/CRAN/refmans/beeca/html/get_marginal_effect.html
# beeca marginal_results columns: TRTVAR, TRTVAL, PARAM, ANALTYP1, STAT, STATVAL, ANALMETH, ANALDESC
# STAT values: "N", "n", "%", "risk", "risk_se", "diff", "diff_se"

BEECA_REQUIRED_COLS <- c("STAT", "STATVAL", "TRTVAL")
BEECA_STAT_VALUES <- c("diff", "diff_se", "risk", "risk_se")

validate_beeca_output <- function(marginal_results) {
  missing <- setdiff(BEECA_REQUIRED_COLS, names(marginal_results))
  if (length(missing) > 0) {
    cli::cli_abort(
      c("Unexpected output from {.pkg beeca}.",
        "x" = "Missing column{?s}: {.field {missing}}.",
        "i" = "Expected columns from {.fn beeca::get_marginal_effect}: {.field {BEECA_REQUIRED_COLS}}.",
        "i" = "Check that {.pkg beeca} version >= 0.1.3 is installed."),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }

  available_stats <- unique(marginal_results$STAT)
  missing_stats <- setdiff(BEECA_STAT_VALUES, available_stats)
  if (length(missing_stats) > 0) {
    cli::cli_warn(
      c("Some expected statistics missing from {.pkg beeca} output.",
        "i" = "Missing: {.val {missing_stats}}.",
        "i" = "Available: {.val {available_stats}}.")
    )
  }
}
```

### inherits()-Based Method Detection

```r
# Source: https://adv-r.hadley.nz/s3.html (S3 class checking best practices)
determine_pooling_class <- function(method) {
  if (inherits(method, "bayes")) return("rubin")
  if (inherits(method, "approxbayes")) return("rubin")
  if (inherits(method, "condmean")) {
    if (identical(method$type, "jackknife")) return("jackknife")
    return("bootstrap")
  }
  if (inherits(method, "bmlmi")) return("bmlmi")

  cli::cli_abort(
    c("Unrecognized imputation method: {.cls {class(method)}}.",
      "i" = "Use {.fn rbmi::method_bayes}, {.fn rbmi::method_approxbayes}, {.fn rbmi::method_condmean}, or {.fn rbmi::method_bmlmi}."),
    class = c("rbmiUtils_error_type", "rbmiUtils_error")
  )
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `tidyr::separate()` | `tidyr::separate_wider_regex()` | tidyr 1.3.0 (Jan 2023) | `separate()` is superseded, still works but no new features. New code should use `separate_wider_regex()`. |
| `stop(..., call. = FALSE)` | `cli::cli_abort(..., class = ...)` | cli 3.4.0+ / tidyverse convention | Standard across tidyverse. Provides inline formatting, pluralization, custom classes. |
| `assertthat::assert_that()` | `cli::cli_abort()` | Ongoing migration in tidyverse | assertthat is no longer actively developed. cli_abort is the replacement for new code. Existing assertthat calls can remain until touched. |
| `class(x)[[2]]` for S3 dispatch | `inherits(x, "classname")` | Always was best practice | Positional indexing breaks when class vector changes. inherits() is stable. |
| Manual deprecation warnings | `lifecycle::deprecate_warn()` | lifecycle 1.0.0 (2021) | Standard R package deprecation workflow. Integrates with roxygen2 badges. |

**Deprecated/outdated:**
- `tidyr::separate()`: Superseded by `separate_wider_delim()`, `separate_wider_position()`, `separate_wider_regex()`. Still works but new code should not use it.
- `assertthat` package: No longer actively developed. Still functional, but `cli` + `rlang` is the current standard.
- `stringsAsFactors = FALSE` in `data.frame()`: Default since R 4.0. Present in current code but unnecessary.

## Open Questions

Things that could not be fully resolved:

1. **Exact parameter name format from gcomp_responder()**
   - What we know: ANCOVA produces `"trt_Week 24"`, `"lsm_ref_Week 24"`. G-computation produces `"trt_<treatment_comparison>_<visit>"` and `"lsm_<treatment_level>_<visit>"` based on source code at lines 91-117 of `analysis_utils.R`.
   - What is unclear: Whether the treatment comparison string in gcomp output (e.g., `"Drug A vs Placebo"`) contains underscores that would break even the regex approach.
   - Recommendation: During implementation, run the actual gcomp functions on test data and inspect the parameter names. The regex should anchor on `"trt_"` and `"lsm_"` prefixes and treat everything after the known prefix as a visit-or-context string. If gcomp parameter names differ structurally, consider separate parsing logic for ANCOVA vs gcomp outputs.

2. **rbmi::as_analysis() exportability**
   - What we know: `rbmi::as_analysis()` is an internal function (not exported). `rbmi::as_class()` IS exported and used by the current `as_analysis2()` code.
   - What is unclear: Whether rbmi will export `as_analysis()` in a future version, which would simplify the refactor.
   - Recommendation: Do not depend on unexported functions. Replicate the class construction using exported `rbmi::as_class()` and `inherits()` checks. This is what the current code does; just fix the positional indexing.

3. **Cross-R-version digest stability**
   - What we know: `rlang::hash()` uses R serialization v3, which may differ across R versions.
   - What is unclear: Whether users realistically reduce in one R version and expand in another.
   - Recommendation: Use `rlang::hash()` for now. Document that the integrity check is valid within the same R major version. If cross-version support is later needed, switch to `digest::digest()` with explicit serialization version.

4. **beeca version 0.2.0 output format**
   - What we know: beeca 0.1.3 is the current CRAN version. Version 0.2.0 exists on r-universe. The output columns (STAT, STATVAL, TRTVAL) are confirmed for 0.1.3.
   - What is unclear: Whether 0.2.0 changes the `marginal_results` structure.
   - Recommendation: Pin `beeca (>= 0.1.3)` in DESCRIPTION. The output schema validation code will catch any breaking changes from future versions.

## Sources

### Primary (HIGH confidence)
- [cli::cli_abort() reference](https://cli.r-lib.org/reference/cli_abort.html) -- message formatting, .envir parameter
- [rlang::abort() reference](https://rlang.r-lib.org/reference/abort.html) -- custom class parameter, metadata attachment
- [rlang::hash() reference](https://rlang.r-lib.org/reference/hash.html) -- XXH128 algorithm, data frame hashing
- [lifecycle::deprecate_warn() reference](https://lifecycle.r-lib.org/articles/communicate.html) -- deprecation workflow, badge integration
- [lifecycle::badge() reference](https://rdrr.io/cran/lifecycle/man/badge.html) -- roxygen2 badge rendering
- [tidyr::separate_wider_regex() reference](https://tidyr.tidyverse.org/reference/separate_wider_delim.html) -- named patterns, too_few handling
- [tidyr::separate() supersession notice](https://tidyr.tidyverse.org/reference/separate.html) -- confirms superseded status
- [beeca::get_marginal_effect() reference](https://search.r-project.org/CRAN/refmans/beeca/html/get_marginal_effect.html) -- return value structure
- [beeca vignette](https://cran.r-project.org/web/packages/beeca/vignettes/estimand_and_implementations.html) -- marginal_results columns: TRTVAR, TRTVAL, PARAM, ANALTYP1, STAT, STATVAL, ANALMETH, ANALDESC
- [rbmi analyse() source](https://github.com/insightsengineering/rbmi/blob/main/R/analyse.R) -- requires imputations object, not data frame
- [rbmi quickstart vignette](https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html) -- analyse() workflow

### Secondary (MEDIUM confidence)
- [R Packages (2e) Lifecycle chapter](https://r-pkgs.org/lifecycle.html) -- deprecation best practices
- [Advanced R: S3 chapter](https://adv-r.hadley.nz/s3.html) -- inherits() vs class() indexing
- [R-hub blog: cli package cliff notes](https://blog.r-hub.io/2023/11/30/cliff-notes-about-cli/) -- practical cli usage patterns

### Tertiary (LOW confidence)
- beeca 0.2.0 output format -- only r-universe, not yet on CRAN, output format unverified

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries verified via official documentation and CRAN
- Architecture patterns: HIGH -- patterns follow established tidyverse/rlang conventions with code examples from official docs
- Pitfalls: HIGH -- pitfall #1 (analyse cannot delegate) verified by reading rbmi source code; all others grounded in codebase analysis

**Research date:** 2026-02-07
**Valid until:** 2026-03-07 (30 days -- all dependencies are stable mature packages)
