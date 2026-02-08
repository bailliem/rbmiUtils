# Phase 5: Data Prep Hardening - Research

**Researched:** 2026-02-08
**Domain:** R package input validation hardening -- cli error messaging, edge case handling, type coercion policy for data preparation functions
**Confidence:** HIGH

## Summary

This phase hardens two existing data preparation functions (`validate_data()` and `prepare_data_ice()`) in `R/data_helpers.R`. The current implementations use a mix of `assertthat::assert_that()`, `stop()` with `call. = FALSE`, and `warning()` with `call. = FALSE` for error reporting. The target is to migrate all error/warning/informational messaging to `cli::cli_abort()`, `cli::cli_warn()`, and `cli::cli_inform()` with custom error classes following the pattern already established in newer files (`analyse_mi_data.R`, `imputation_storage.R`, `efficacy_table.R`).

The seven requirements (HRD-01 through HRD-07) map directly to specific, bounded code changes:
- **HRD-01** (malformed interaction terms): Add regex validation of covariate term syntax before the existing `strsplit()` parsing in `validate_data()`.
- **HRD-02** (NULL strategy): Change `prepare_data_ice()` line 309 from silently defaulting to erroring with `cli::cli_abort()`.
- **HRD-03** (character visit warning): Add a `cli::cli_warn()` check early in `prepare_data_ice()` when the visit column is character.
- **HRD-04** (empty data frame): Add a `nrow(data) == 0` check near the top of `validate_data()`.
- **HRD-05** (all-NA covariates): Modify the existing covariate NA loop to detect all-NA columns and emit `cli::cli_warn()` instead of erroring.
- **HRD-06** (batched type warnings): Collect character-column warnings into a vector and emit one `cli::cli_warn()` call instead of individual `warning()` calls in a loop.
- **HRD-07** (edge case tests): Add test cases for single subject, single visit, all-NA outcome, and all-complete data scenarios.

No new functions are created. No new dependencies are added (`cli` is already in Imports at `>= 3.6.0`). The work is entirely within `R/data_helpers.R` and `tests/testthat/test-data_helpers.R`.

**Primary recommendation:** Migrate all `stop()`/`warning()` calls in `validate_data()` and `prepare_data_ice()` to `cli::cli_abort()`/`cli::cli_warn()` with custom error classes (`rbmiUtils_error_validation`, `rbmiUtils_error`), following the exact pattern established in `R/analyse_mi_data.R`. Add the seven specific validation checks, and update tests to assert on error class where possible.

## Standard Stack

### Core (Already Imported -- No Changes)

| Library | Version | Purpose | Why Relevant |
|---------|---------|---------|--------------|
| cli | >= 3.6.0 | `cli_abort()`, `cli_warn()`, `cli_inform()` for formatted error messages | Already in DESCRIPTION Imports. Extensively used in `analyse_mi_data.R`, `imputation_storage.R`, `efficacy_table.R`, `pool_methods.R`. No version change needed. |
| rlang | current | Underlies `cli_abort()` (passes `class` to `rlang::abort()`), `.data` pronoun | Already in DESCRIPTION Imports. |
| assertthat | current | Existing validation calls remain for now | Already in DESCRIPTION Imports. Only the targeted functions get migrated; `assertthat` calls in `summarise_missingness()` stay as-is (out of scope). |
| testthat | >= 3.0.0 | Test framework with `expect_error(class = ...)` support | Already in DESCRIPTION Suggests. Edition 3 config confirmed. |

### Not Needed

| Library | Why Not |
|---------|---------|
| lifecycle | No deprecations in this phase -- we are hardening existing functions, not removing/renaming anything. |
| digest/rlang::hash | No integrity checks needed for data prep functions. |

### Alternatives Considered

| Instead of | Could Use | Decision |
|------------|-----------|----------|
| `cli::cli_abort()` | `stop()` with manual formatting | Locked: use cli. Project already uses cli extensively in newer files. |
| `cli::cli_warn()` | `warning()` | Locked: use cli for consistency. |
| `cli::cli_inform()` | `message()` | Use cli_inform for informational messages (e.g., all-complete data notification). |
| Collecting all errors | Fail-fast | Locked: fail fast on first fatal error. Exception: `validate_data()` already collects non-fatal issues -- the existing multi-issue pattern stays for warnings/type checks within validate_data, but fatal errors stop immediately. |

## Architecture Patterns

### Current File Structure (No Changes)

```
R/
├── data_helpers.R          # Contains validate_data(), prepare_data_ice(), summarise_missingness()
tests/testthat/
├── test-data_helpers.R     # Tests for all three functions
```

Both target functions live in `data_helpers.R`. No file moves or splits needed.

### Pattern 1: cli_abort with Custom Error Classes (Established Pattern)

**What:** Replace `stop(..., call. = FALSE)` and `assertthat::assert_that()` with `cli::cli_abort()` including named bullets and custom classes.
**When to use:** Every error path in `validate_data()` and `prepare_data_ice()`.
**Source:** Established in `R/analyse_mi_data.R` lines 89-170, `R/imputation_storage.R` lines 67-113.

```r
# Before (current in data_helpers.R):
stop(paste(c("Data validation failed:", issues), collapse = "\n  - "),
     call. = FALSE)

# After (target):
cli::cli_abort(
  c(
    "Data validation failed.",
    "x" = issues  # cli_abort handles character vectors as separate bullets
  ),
  class = c("rbmiUtils_error_validation", "rbmiUtils_error")
)
```

For single-issue errors:
```r
# Before:
assertthat::assert_that(
  is.data.frame(data),
  msg = "`data` must be a data.frame"
)

# After:
if (!is.data.frame(data)) {
  cli::cli_abort(
    c(
      "{.arg data} must be a {.cls data.frame}, not {.cls {class(data)}}."
    ),
    class = c("rbmiUtils_error_validation", "rbmiUtils_error")
  )
}
```

### Pattern 2: Batched Warnings via cli_warn (New for This Phase)

**What:** Collect multiple non-fatal type issues into a single `cli::cli_warn()` call with bullet points.
**When to use:** HRD-06 -- character column warnings in `validate_data()`.

```r
# Before (current -- one warning() per character column):
} else if (is.character(data[[col]])) {
  warning(sprintf(
    "Column `%s` is character and will be converted to factor by `rbmi::draws()`",
    col
  ), call. = FALSE)
}

# After (collect and batch):
char_cols <- character(0)
for (col in factor_cols) {
  if (is.character(data[[col]])) {
    char_cols <- c(char_cols, col)
  }
}
if (length(char_cols) > 0) {
  cli::cli_warn(
    c(
      "{length(char_cols)} column{?s} {?is/are} character and will be coerced
       to factor by {.fn rbmi::draws}.",
      "i" = "Column{?s}: {.field {char_cols}}.",
      "i" = "Convert to factor explicitly for control over level ordering:
             {.code data${char_cols[1]} <- factor(data${char_cols[1]})}"
    ),
    class = c("rbmiUtils_warning_coercion", "rbmiUtils_warning")
  )
}
```

### Pattern 3: cli_inform for Informational Messages

**What:** Use `cli::cli_inform()` for non-warning, non-error informational output.
**When to use:** When all data is complete (no missing outcomes), confirming no ICE imputations needed.

```r
cli::cli_inform(
  c(
    "v" = "All outcome values are complete -- no ICE imputations needed.",
    "i" = "You can proceed directly to {.fn rbmi::draws} without {.arg data_ice}."
  ),
  class = "rbmiUtils_info"
)
```

### Pattern 4: Fail-Fast on First Fatal Error

**What:** Stop immediately when a fatal condition is detected, do not collect.
**When to use:** All new validation checks (empty data frame, NULL strategy, etc.).
**Important nuance:** The existing `validate_data()` collects multiple issues before stopping. The user decision says "fail fast on first fatal error." The reconciliation is:
- **Fatal errors** (can't proceed at all): fail immediately. E.g., empty data frame, not a data.frame.
- **Non-fatal issues** (can detect more): the existing collection pattern stays. E.g., missing columns + type mismatches can all be collected and reported together.
- New fatal checks (HRD-04 empty data, HRD-01 malformed interaction terms) should fail immediately since they prevent further validation.

```r
# Fatal: stop immediately (before any other checks)
if (nrow(data) == 0) {
  cli::cli_abort(
    c(
      "{.arg data} has 0 rows.",
      "i" = "Provide a data frame with at least one subject-visit observation."
    ),
    class = c("rbmiUtils_error_validation", "rbmiUtils_error")
  )
}
```

### Pattern 5: Interaction Term Validation (HRD-01)

**What:** Validate covariate term syntax before parsing.
**When to use:** In `validate_data()`, after extracting `covariate_cols` from `vars$covariates`.

A well-formed interaction term matches the pattern: `variable_name`, `var1*var2`, `var1:var2`, or `var1*var2*var3` (arbitrarily deep). Each variable name must be a non-empty R identifier-like string. Malformed patterns include: `"A*"`, `":B"`, `""`, `"*"`, `"A*:B"`, `"A**B"`.

```r
# Validate each covariate term
if (!is.null(covariate_cols) && length(covariate_cols) > 0) {
  for (term in covariate_cols) {
    if (!nzchar(trimws(term))) {
      cli::cli_abort(
        c(
          "Empty covariate term found in {.arg vars$covariates}.",
          "i" = "Remove empty strings from the covariates vector."
        ),
        class = c("rbmiUtils_error_validation", "rbmiUtils_error")
      )
    }
    # Check for malformed interaction: leading/trailing operators, consecutive operators
    if (grepl("^[*:]|[*:]$|[*:]{2}", trimws(term))) {
      cli::cli_abort(
        c(
          "Malformed interaction term in {.arg vars$covariates}: {.val {term}}.",
          "x" = "Interaction terms must have variable names on both sides of {.code *} or {.code :}.",
          "i" = "Example: {.code c(\"BASE\", \"TRT*STRATA\", \"TRT:BASE\")}"
        ),
        class = c("rbmiUtils_error_validation", "rbmiUtils_error")
      )
    }
  }
}
```

### Custom Error/Warning Class Hierarchy

Extends the hierarchy established in Phase 1:

```
rbmiUtils_error                    # Base class for all rbmiUtils errors
  rbmiUtils_error_validation       # Input validation failures
  rbmiUtils_error_type             # Type mismatch errors

rbmiUtils_warning                  # Base class for all rbmiUtils warnings
  rbmiUtils_warning_coercion       # Type coercion warnings

rbmiUtils_info                     # Informational messages
```

The `rbmiUtils_warning` and `rbmiUtils_info` classes are new to this phase. The error classes already exist in the codebase.

### Anti-Patterns to Avoid

- **Collecting errors then stopping:** For fatal errors (empty data, malformed terms), stop immediately. Do not add them to the `issues` vector.
- **Using `assertthat::assert_that()` in new code:** All new validation paths use `cli::cli_abort()`. Existing `assertthat` calls in `prepare_data_ice()` lines 281-293 and `summarise_missingness()` lines 440-444 should be migrated while touching those functions, but `summarise_missingness()` is out of scope.
- **Auto-coercing types:** Locked decision -- character visit columns warn but do NOT coerce. Character-to-numeric columns are hard errors.
- **Checking unused columns:** Locked decision -- only validate what the model references.
- **Emitting individual warnings in a loop:** Use the batched pattern (Pattern 2) instead.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Formatted error messages with bullets | `paste0()` + `stop()` | `cli::cli_abort()` with named vector | Already established pattern in 7+ files. Handles pluralization `{?s}`, inline markup `{.arg}`, `{.cls}`, `{.field}`. |
| Pluralized messages | `if (length > 1) "s" else ""` | cli pluralization `{?s}`, `{?is/are}` | Built into cli. Zero-effort pluralization. |
| Warning batching | Manual accumulation + `paste(collapse)` | `cli::cli_warn()` with `{.field {vec}}` | cli automatically formats vectors with Oxford commas. |
| Regex validation of R identifiers | Character-by-character parsing | `grepl()` with anchor patterns | Simple regex patterns suffice for interaction term validation. No need for a full parser. |
| Custom condition classes | `structure(list(), class = ...)` | `cli_abort(..., class = c("specific", "general"))` | cli passes through to rlang::abort which creates proper condition objects. |

**Key insight:** This phase is primarily a migration of existing validation logic from `stop()`/`warning()`/`assertthat` to `cli`, plus adding ~7 targeted new checks. The cli package handles all formatting complexity. The work is bounded and well-defined.

## Common Pitfalls

### Pitfall 1: validate_data() Currently Collects Issues -- Fail-Fast Needs Care

**What goes wrong:** Blindly changing `validate_data()` to fail on first error breaks existing behavior where multiple issues are reported together. The existing test `"validate_data collects multiple issues in one error"` verifies multi-issue collection.

**Why it happens:** The user decision says "fail fast on first fatal error." But the existing function is designed to collect non-fatal issues.

**How to avoid:** Distinguish between:
- **Pre-checks** (fatal, fail immediately): empty data frame (HRD-04), data not a data.frame, malformed interaction terms (HRD-01). These go BEFORE the issues collection loop.
- **Validation issues** (collectible): missing columns, type mismatches, duplicate rows, etc. These continue to be collected into the `issues` vector and reported together at the end with `cli::cli_abort()`.

**Warning signs:** The test `"validate_data collects multiple issues in one error"` fails if you accidentally make type checks fail-fast instead of collectible.

### Pitfall 2: cli_abort Message Vector Semantics Differ from paste()

**What goes wrong:** The current `stop(paste(c("header:", issues), collapse = "\n  - "))` creates a single string. `cli::cli_abort(c("header", issues))` treats each element as a separate bullet. Unnamed elements get no bullet prefix; named elements (x, i, !) get specific icons.

**Why it happens:** Different formatting semantics between `paste()` and `cli_bullets()`.

**How to avoid:** When converting the multi-issue error, use named bullets:
```r
# Convert issues vector to named bullets
bullets <- setNames(issues, rep("x", length(issues)))
cli::cli_abort(
  c("Data validation failed:", bullets),
  class = c("rbmiUtils_error_validation", "rbmiUtils_error")
)
```

**Warning signs:** Error messages appear as a flat list without bullet formatting, or show raw vector names instead of symbols.

### Pitfall 3: Existing Tests Match Error Messages by Regex

**What goes wrong:** Tests like `expect_error(validate_data(dat, vars), "must be a data.frame")` match the literal text of `stop()` output. After migration to `cli::cli_abort()`, the message text changes (cli adds formatting, changes wording to use `{.cls data.frame}` markup).

**Why it happens:** cli renders `{.cls data.frame}` as `<data.frame>` in non-ANSI contexts. The test regex `"must be a data.frame"` might still match, or it might not depending on exact wording changes.

**How to avoid:** Update tests in two ways:
1. **Prefer class-based matching:** `expect_error(expr, class = "rbmiUtils_error_validation")` -- this is the established pattern in newer tests (see `test-analyse_mi_data.R`, `test-tidiers.R`).
2. **Keep keyword matching as secondary:** `expect_error(expr, "data.frame")` still works if the cli message contains that keyword.
3. Update existing regex patterns where they would break.

**Warning signs:** Tests that use `expect_error(expr, "exact old message text")` fail after migration.

### Pitfall 4: prepare_data_ice() assertthat Calls Need Migration Too

**What goes wrong:** `prepare_data_ice()` uses both `assertthat::assert_that()` (lines 281-293, 326-332) and `stop()` (lines 302-305, 314-317, 319-322, 344-347). If only the `stop()` calls are migrated, the function has inconsistent error styles.

**Why it happens:** The function was originally written with `assertthat` and `stop()` before the project adopted `cli`.

**How to avoid:** Migrate ALL validation code in `prepare_data_ice()` to `cli::cli_abort()` in a single pass. The `assertthat` calls are straightforward 1:1 conversions.

**Warning signs:** Some errors in the function show cli formatting while others show bare `assertthat` messages.

### Pitfall 5: Character Visit Warning Must Not Auto-Coerce

**What goes wrong:** After warning about a character visit column, developer reflexively adds `data[[vars$visit]] <- factor(data[[vars$visit]])` to "fix" the issue.

**Why it happens:** Natural instinct to resolve the warning. But the locked decision says "warn but do NOT auto-coerce" because auto-coercion would use alphabetical ordering for factor levels, which is wrong for visit sequences like "Week 4", "Week 8", "Week 12" (alphabetical: Week 12, Week 4, Week 8).

**How to avoid:** The warning message must explicitly guide the user to convert with explicit levels:
```
! Visit column `AVISIT` is character, not factor.
i Factor ordering controls visit sequence in the imputation model.
i Convert with explicit levels: `data$AVISIT <- factor(data$AVISIT, levels = c("Week 4", "Week 8", "Week 12"))`
```

**Warning signs:** Tests pass but the imputation model silently uses wrong visit ordering because auto-coercion applied alphabetical levels.

### Pitfall 6: All-NA Covariate Must Warn+Drop, Not Error

**What goes wrong:** The current code errors when covariates have ANY missing values. HRD-05 requires that all-NA columns specifically should warn and be dropped, while partially-NA columns remain errors.

**Why it happens:** The current loop treats all NA covariates the same way.

**How to avoid:** Split the covariate check into two passes:
1. **First pass:** Detect all-NA covariate columns. Warn and note them for dropping.
2. **Second pass:** For remaining covariates, check for partial NAs (existing error behavior).

The "drop from the model" means removing the all-NA column names from the covariate list that validate_data reports/uses -- but validate_data does not modify the data. The warn+drop means: warn the user, and do not report the all-NA column as an error (effectively "dropping it from consideration" in validation). The user must still handle it in their actual data.

**Warning signs:** validate_data errors on all-NA covariates instead of warning.

## Code Examples

### Example 1: Malformed Interaction Term Validation (HRD-01)

```r
# Regex to detect malformed terms:
# - leading operator: "^[*:]"
# - trailing operator: "[*:]$"
# - consecutive operators: "[*:]{2}"
# - empty after trimming: nzchar check

validate_interaction_terms <- function(covariate_cols) {
  if (is.null(covariate_cols) || length(covariate_cols) == 0) return(invisible(NULL))

  for (term in covariate_cols) {
    trimmed <- trimws(term)
    if (!nzchar(trimmed)) {
      cli::cli_abort(
        c(
          "Empty covariate term found in {.arg vars$covariates}.",
          "i" = "Remove empty strings from the covariates vector."
        ),
        class = c("rbmiUtils_error_validation", "rbmiUtils_error")
      )
    }
    if (grepl("^[*:]|[*:]$|[*:]{2}", trimmed)) {
      cli::cli_abort(
        c(
          "Malformed interaction term: {.val {term}}.",
          "x" = "Interaction terms must have variable names on both sides of {.code *} or {.code :}.",
          "i" = "Example valid terms: {.code c(\"BASE\", \"TRT*STRATA\", \"TRT:BASE\")}"
        ),
        class = c("rbmiUtils_error_validation", "rbmiUtils_error")
      )
    }
  }
}
```

### Example 2: NULL Strategy Error in prepare_data_ice (HRD-02)

```r
# Before (silently defaults):
if (is.null(vars$strategy) || !nzchar(vars$strategy)) {
  vars$strategy <- "strategy"
}

# After (error with guidance):
if (is.null(vars$strategy) || !nzchar(vars$strategy)) {
  cli::cli_abort(
    c(
      "{.field vars$strategy} must be defined when preparing ICE data.",
      "i" = "Set it via: {.code rbmi::set_vars(strategy = \"strategy_column_name\", ...)}",
      "i" = "The strategy column in your data should contain values like
             {.val MAR}, {.val JR}, {.val CR}, {.val CIR}, or {.val LMCF}."
    ),
    class = c("rbmiUtils_error_validation", "rbmiUtils_error")
  )
}
```

Note: `rbmi::set_vars()` defaults strategy to `"strategy"` (the string literal). So `vars$strategy` is only NULL if the user manually tampered with the vars object or constructed their own list without a strategy field. The check protects against hand-built lists.

### Example 3: Character Visit Warning in prepare_data_ice (HRD-03)

```r
# After validating visit column exists:
if (is.character(data[[vars$visit]])) {
  cli::cli_warn(
    c(
      "Visit column {.field {vars$visit}} is character, not factor.",
      "!" = "Character visits use alphabetical ordering, which is often wrong
             for clinical visit sequences.",
      "i" = "Convert to factor with explicit level ordering:",
      " " = "{.code data${vars$visit} <- factor(data${vars$visit},
             levels = c(\"Week 4\", \"Week 8\", \"Week 12\"))}"
    ),
    class = c("rbmiUtils_warning_coercion", "rbmiUtils_warning")
  )
}
```

### Example 4: Empty Data Frame Error (HRD-04)

```r
# Add immediately after the is.data.frame check passes:
if (nrow(data) == 0) {
  cli::cli_abort(
    c(
      "{.arg data} has 0 rows.",
      "i" = "Provide a data frame with at least one subject-visit observation."
    ),
    class = c("rbmiUtils_error_validation", "rbmiUtils_error")
  )
}
```

### Example 5: All-NA Covariate Warning + Drop (HRD-05)

```r
# Replace the existing covariate NA check loop:
existing_covars <- intersect(all_covar_names, existing_cols)
all_na_covars <- character(0)
partial_na_covars <- character(0)

for (col in existing_covars) {
  n_na <- sum(is.na(data[[col]]))
  if (n_na == nrow(data)) {
    all_na_covars <- c(all_na_covars, col)
  } else if (n_na > 0) {
    partial_na_covars <- c(partial_na_covars, col)
    issues <- c(issues, sprintf(
      "Covariate `%s` has %d missing value(s) out of %d rows", col, n_na, nrow(data)
    ))
  }
}

if (length(all_na_covars) > 0) {
  cli::cli_warn(
    c(
      "{length(all_na_covars)} covariate{?s} {?is/are} entirely NA and will
       be excluded from validation.",
      "i" = "Column{?s}: {.field {all_na_covars}}.",
      "i" = "Consider removing {?this/these} column{?s} from {.arg vars$covariates}
             or investigating why all values are missing."
    ),
    class = c("rbmiUtils_warning_coercion", "rbmiUtils_warning")
  )
}
```

### Example 6: Batched Type Coercion Warnings (HRD-06)

```r
# Collect character columns, emit single warning
char_cols <- character(0)
for (col in factor_cols) {
  if (!is.factor(data[[col]]) && !is.character(data[[col]])) {
    issues <- c(issues, sprintf(
      "Column `%s` must be a factor (found %s)", col, class(data[[col]])[1]
    ))
  } else if (is.character(data[[col]])) {
    char_cols <- c(char_cols, col)
  }
}

if (length(char_cols) > 0) {
  cli::cli_warn(
    c(
      "{length(char_cols)} column{?s} {?is/are} character instead of factor.",
      "i" = "Column{?s}: {.field {char_cols}}.",
      "i" = "{.fn rbmi::draws} will auto-coerce, but explicit conversion
             gives you control over level ordering.",
      "i" = "Example: {.code data${char_cols[1]} <- factor(data${char_cols[1]})}"
    ),
    class = c("rbmiUtils_warning_coercion", "rbmiUtils_warning")
  )
}
```

### Example 7: Converting Issues Vector to cli Bullets

```r
# At the end of validate_data(), replace:
# stop(paste(c("Data validation failed:", issues), collapse = "\n  - "), call. = FALSE)

# With:
if (length(issues) > 0) {
  bullets <- setNames(issues, rep("x", length(issues)))
  cli::cli_abort(
    c("Data validation failed:", bullets),
    class = c("rbmiUtils_error_validation", "rbmiUtils_error")
  )
}
```

### Example 8: Testing cli Errors by Class

```r
# New pattern (preferred):
test_that("validate_data errors on empty data frame", {
  dat <- make_test_data()[0, ]  # zero rows
  vars <- make_test_vars()
  expect_error(
    validate_data(dat, vars),
    class = "rbmiUtils_error_validation"
  )
})

# Hybrid pattern (class + keyword for specificity):
test_that("validate_data rejects malformed interaction term", {
  dat <- make_test_data()
  vars <- rbmi::set_vars(
    subjid = "USUBJID", visit = "AVISIT", group = "TRT",
    outcome = "CHG", covariates = c("A*")
  )
  expect_error(
    validate_data(dat, vars),
    "Malformed interaction",
    class = "rbmiUtils_error_validation"
  )
})
```

## State of the Art

| Old Approach (Current in data_helpers.R) | Current Best Practice (Used in newer files) | Impact |
|------------------------------------------|---------------------------------------------|--------|
| `stop(paste(...), call. = FALSE)` | `cli::cli_abort(c(...), class = ...)` | Formatted bullets, custom classes, catchable conditions |
| `assertthat::assert_that(cond, msg = ...)` | `if (!cond) cli::cli_abort(...)` | More control over message formatting, class hierarchy |
| `warning(sprintf(...), call. = FALSE)` | `cli::cli_warn(c(...), class = ...)` | Batching, bullets, consistent style |
| Individual `warning()` in loop | Single `cli::cli_warn()` after collection | Cleaner output, user sees one message instead of N |
| `expect_error(expr, "regex")` | `expect_error(expr, class = "pkg_error")` | Decoupled from message text, stable under rewording |

**Note:** assertthat is no longer actively developed. The `cli` + `rlang` pattern is the current tidyverse standard.

## Open Questions

1. **All-NA covariate "drop" semantics**
   - What we know: User decided "warn and automatically drop from the model." `validate_data()` does not modify data -- it only validates.
   - What is unclear: Does "drop from the model" mean (a) validate_data should not error on the all-NA covariate (skip it during validation), or (b) validate_data should actually remove the column from `vars$covariates` and return a modified vars object?
   - Recommendation: Interpret as (a) -- skip the all-NA covariate in validation, warn the user, and let them handle it. Changing the function's return type from `invisible(TRUE)` to a modified vars object would be a breaking API change. The warning message should tell the user to remove the column from covariates.
   - Confidence: MEDIUM -- the intent is clear but the mechanism requires a design choice.

2. **Single-subject / single-visit thresholds (Claude's Discretion)**
   - What we know: Edge cases with 1 subject or 1 visit should not crash. Currently, `validate_data()` does not check for minimum subject/visit counts.
   - What is unclear: Whether to warn, error, or simply allow these through.
   - Recommendation: Allow through without error -- these are technically valid inputs. The statistical model downstream may fail, but that is not validate_data's responsibility. Add test coverage to verify no crash. If desired, emit `cli::cli_inform()` as a hint: "Data contains only 1 subject" / "Data contains only 1 visit."
   - Confidence: HIGH -- defensive validation should ensure no crash, not enforce statistical adequacy.

3. **All-NA outcome severity (Claude's Discretion)**
   - What we know: An outcome column that is entirely NA means no observed data to model. This is likely a data error.
   - What is unclear: Whether this should be an error or a warning.
   - Recommendation: Make it an error (`cli::cli_abort()`). A completely missing outcome means the analysis cannot proceed -- there is nothing to impute or model. This is qualitatively different from missing some outcomes (which is normal in clinical trial data with dropouts).
   - Confidence: HIGH -- an all-NA outcome is unambiguously a problem.

4. **All-complete data informational message location**
   - What we know: When all data is complete, the user should see a confirmation message.
   - What is unclear: Whether this message belongs in `validate_data()` or `prepare_data_ice()`.
   - Recommendation: Put it in `prepare_data_ice()` when no ICE flags are found (nrow(ice_data) == 0). This is where the user is actively trying to prepare ICE data and discovering none is needed. Also add an informational message in `validate_data()` when no outcome values are missing.
   - Confidence: MEDIUM -- reasonable location but could go either way.

## Sources

### Primary (HIGH confidence)

- **Codebase inspection:** Direct reading of `R/data_helpers.R` (validate_data lines 64-235, prepare_data_ice lines 279-392), `tests/testthat/test-data_helpers.R` (437 lines, 37 test blocks), and 7 other R source files using cli patterns.
- **cli package documentation:** [cli_abort reference](https://cli.r-lib.org/reference/cli_abort.html) -- message formatting, named bullets, class parameter.
- **cli semantic markup:** [cli semantic CLI article](https://cli.r-lib.org/articles/semantic-cli.html) -- `{.arg}`, `{.cls}`, `{.field}`, `{.fn}`, `{.code}`, `{.val}` inline formatters.
- **cli bullets reference:** [cli_bullets reference](https://cli.r-lib.org/reference/cli_bullets.html) -- named vector mapping (x, i, !, v, >, *, space).
- **rbmi::set_vars() behavior:** Verified via direct R execution that `strategy` defaults to `"strategy"` (string literal), and `set_vars(strategy = NULL)` errors. NULL strategy only occurs with hand-built lists.
- **Interaction term parsing edge cases:** Verified via direct R execution that `"A*"`, `":B"`, `""`, `"*"` are silently accepted by current `strsplit()` logic.

### Secondary (MEDIUM confidence)

- **Phase 1 research (`01-RESEARCH.md`):** Established the custom error class hierarchy, cli migration pattern, and assertthat deprecation approach. This phase follows the same patterns.
- **Existing codebase conventions (`CONVENTIONS.md`, `TESTING.md`):** Test fixture patterns, assertion usage, error handling conventions documented from codebase analysis.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all libraries already in DESCRIPTION and verified in use
- Architecture patterns: HIGH -- patterns directly derived from existing codebase files written in earlier phases
- Pitfalls: HIGH -- all pitfalls identified from reading the actual current code and understanding the migration requirements
- Code examples: HIGH -- based on actual current code patterns and verified cli API

**Research date:** 2026-02-08
**Valid until:** 2026-03-08 (30 days -- all dependencies are stable, no version changes anticipated)
