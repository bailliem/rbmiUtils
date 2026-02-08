# Phase 2: Print/Summary & ARD Conversion - Research

**Researched:** 2026-02-08
**Domain:** S3 print/summary methods for rbmi objects + pharmaverse ARD (cards package) conversion
**Confidence:** HIGH

## Summary

This phase adds two capabilities: (1) enhanced print/summary methods for pool and analysis objects using cli formatting, and (2) a `pool_to_ard()` function that converts rbmi pool results into the pharmaverse ARD (Analysis Results Dataset) standard using the `cards` package.

The existing codebase already has `print.analysis()` and `summary.analysis()` registered in NAMESPACE (overriding rbmi's `print.analysis`). The pool class from rbmi already has `print.pool` (from rbmi) but does NOT have `summary.pool`. The requirements ask for enhanced versions of both pool print/summary and analysis print/summary.

For ARD conversion, the `cards` package (v0.7.1 installed) provides `ard_identity()` for pre-computed statistics and `as_card()` for direct construction. The recommended approach is direct construction with `as_card()` because it gives full control over grouping columns (visit, parameter_type, lsm_type/arm) and avoids per-row function call overhead. The ARD output must have `stat` as a list column, include `warning` and `error` list columns, and optionally include a `method` stat_name row.

**Primary recommendation:** Use `as_card()` for ARD construction (not `ard_identity()`), keep `cards` as Suggests (not Imports), and enhance existing print/summary methods with `cli` formatting. Override `print.pool` from rbmi in NAMESPACE (same pattern as existing `print.analysis` override).

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| cards | >= 0.4.0 | ARD data frame creation (`as_card()`, `check_ard_structure()`) | Pharmaverse standard for ARD; CDISC-endorsed |
| cli | >= 3.6.0 | Formatted console output for print/summary methods | Already in Imports from Phase 1; semantic formatting |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| dplyr | (existing) | Data manipulation for ARD construction | Already in Imports |
| tibble | (existing) | Tibble output from tidy_pool_obj() | Already in Suggests |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `as_card()` direct construction | `ard_identity()` per parameter | `ard_identity()` simpler for single stat lists but requires per-row calls + manual group column addition; `as_card()` is more efficient for batch conversion |
| `as_card()` direct construction | `tidy_as_ard()` | `tidy_as_ard()` is lifecycle "questioning" and designed for tidied model output, not pre-computed stats |
| Override `print.pool` | Standalone function `print_pool()` | S3 dispatch is idiomatic R; overriding is the same pattern used for `print.analysis` already |

**Installation:**
```r
# cards as Suggests (not Imports) - see dependency decision below
# Already available: cli, dplyr in Imports
```

## Architecture Patterns

### Recommended Project Structure
```
R/
  analyse_mi_data.R     # Existing: print.analysis, summary.analysis (ENHANCED)
  pool_methods.R        # NEW: print.pool, summary.pool (PRT-01, PRT-02)
  ard_conversion.R      # NEW: pool_to_ard() (ARD-01, ARD-02, ARD-03)
tests/testthat/
  test-pool_methods.R   # NEW: tests for print/summary pool
  test-ard_conversion.R # NEW: tests for pool_to_ard
  test-analyse_mi_data.R # Existing: add enhanced print/summary tests
```

### Pattern 1: S3 Method Override for Pool Objects
**What:** Register `print.pool` and `summary.pool` in rbmiUtils NAMESPACE, overriding rbmi's `print.pool`
**When to use:** Always - this is the standard S3 dispatch pattern
**Confidence:** HIGH - rbmiUtils already does this for `print.analysis`

```r
# In NAMESPACE (via roxygen2):
# S3method(print,pool)
# S3method(summary,pool)

#' @export
print.pool <- function(x, digits = 2, ...) {
  # Use cli for formatted output
  tidy_df <- tidy_pool_obj(x)

  cli::cli_h1("Pool Object")
  cli::cli_text("{.val {nrow(tidy_df)}} parameters across {.val {length(unique(tidy_df$visit))}} visit{?s}")
  # ... formatted table output
  invisible(x)
}

#' @export
summary.pool <- function(object, ...) {
  # More detailed output with visit-level breakdown
  # Returns summary list invisibly
}
```

### Pattern 2: ARD Construction via as_card()
**What:** Build long-format data frame with proper list columns, then call `cards::as_card()`
**When to use:** For `pool_to_ard()` function
**Confidence:** HIGH - verified with cards 0.7.1

```r
pool_to_ard <- function(pool_obj, ...) {
  # 1. Get tidy tibble
  tidy_df <- tidy_pool_obj(pool_obj)

  # 2. Build long-form rows with grouping columns
  stat_names <- c("estimate", "std.error", "conf.low", "conf.high", "p.value")
  stat_labels <- c("Estimate", "Std. Error", "95% CI Lower", "95% CI Upper", "p-value")

  rows <- lapply(seq_len(nrow(tidy_df)), function(i) {
    r <- tidy_df[i, ]
    data.frame(
      group1 = "visit",
      group1_level = I(as.list(rep(r$visit, 5))),
      group2 = "parameter_type",
      group2_level = I(as.list(rep(r$parameter_type, 5))),
      variable = r$parameter,
      context = "rbmi_pool",
      stat_name = stat_names,
      stat_label = stat_labels,
      stat = I(list(r$est, r$se, r$lci, r$uci, r$pval)),
      fmt_fun = I(as.list(rep(1L, 5))),
      warning = I(as.list(rep(list(NULL), 5))),
      error = I(as.list(rep(list(NULL), 5))),
      stringsAsFactors = FALSE
    )
  })

  ard_df <- do.call(rbind, rows)
  cards::tidy_ard_column_order(cards::as_card(ard_df))
}
```

### Pattern 3: Conditional Dependency Check
**What:** Check for cards availability at function call time (Suggests pattern)
**When to use:** When cards is in Suggests, not Imports
**Confidence:** HIGH - standard CRAN pattern

```r
pool_to_ard <- function(pool_obj, ...) {
  if (!requireNamespace("cards", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg cards} is required for ARD conversion.",
      class = "rbmiUtils_error_dependency"
    )
  }
  # ... proceed with conversion
}
```

### Anti-Patterns to Avoid
- **Using cat() directly in print methods:** Use cli formatting instead for consistent, colorized output. The existing print.analysis uses cat() - this should be modernized.
- **Returning the data frame from print:** print.pool must return `invisible(x)` (the original object), not the formatted output.
- **Hardcoding stat column names in ARD:** Use the standard names from cards conventions: `estimate`, `std.error`, `conf.low`, `conf.high`, `p.value` (matching broom-style names that cardx uses).
- **Omitting list columns in ARD:** The `stat`, `fmt_fun`, `warning`, and `error` columns MUST be list columns, not atomic vectors. This is a cards requirement.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| ARD data frame structure | Custom data frame with atomic columns | `cards::as_card()` with proper list columns | ARD requires list columns for stat/warning/error; cards validates structure |
| ARD column ordering | Manual column selection | `cards::tidy_ard_column_order()` | Standard ordering: group cols -> variable -> context -> stat_name -> stat -> fmt_fun -> warning/error |
| ARD structure validation | Custom validation logic | `cards::check_ard_structure()` | Checks required columns, ordering, method row presence |
| Console formatting | sprintf/paste for aligned output | `cli::cli_text()`, `cli::cli_h1()`, `cli::cli_bullets()` | Consistent formatting, automatic color support, width-aware |
| Significance flags | Custom p-value thresholds | Existing `format_pvalue()` function | Already handles thresholds, NA, HTML formatting |

**Key insight:** The ARD format has very specific structural requirements (list columns, specific column names, class 'card'). Building this manually without `as_card()` or `ard_identity()` will inevitably miss edge cases that downstream consumers (gtsummary) depend on.

## Common Pitfalls

### Pitfall 1: S3 Method Registration Conflict with rbmi
**What goes wrong:** Overriding `print.pool` from rbmi causes "method overwritten" warning on package load
**Why it happens:** Both rbmi and rbmiUtils register print.pool as S3 method
**How to avoid:** This is actually acceptable and expected -- rbmiUtils already does this for `print.analysis`. The warning `"Registered S3 method overwritten by 'rbmiUtils': method print.analysis from rbmi"` appears at load time. The same will happen for print.pool. This is standard behavior when extending a package.
**Warning signs:** If print.pool doesn't appear in NAMESPACE after roxygen2, the method won't be dispatched.

### Pitfall 2: ARD stat Column Must Be List Column
**What goes wrong:** Creating `stat` as an atomic numeric vector instead of a list column causes `check_ard_structure()` to fail and downstream gtsummary to break
**Why it happens:** Natural instinct is `stat = c(-2.5, 0.8)` but cards requires `stat = list(-2.5, 0.8)` or `stat = I(list(-2.5, 0.8))`
**How to avoid:** Always use `I(list(...))` or `I(as.list(...))` when building the data frame, or construct via `ard_identity()` which handles this automatically
**Warning signs:** `str(ard$stat)` shows `num` instead of `List`

### Pitfall 3: ARD group_level Must Be List Column Too
**What goes wrong:** Group level columns (group1_level, group2_level) must also be list columns
**Why it happens:** Same as stat column -- the ARD format stores heterogeneous types in list columns
**How to avoid:** Use `I(as.list(...))` for all `*_level` columns
**Warning signs:** `check_ard_structure()` warning about column types

### Pitfall 4: Missing method Row Warning
**What goes wrong:** `check_ard_structure()` emits note "Expecting a row with stat_name = 'method', but it is not present"
**Why it happens:** Many ARD objects include a "method" row describing the statistical method used
**How to avoid:** Include a method row with stat = "Rubin's pooling" (or appropriate method name). Alternatively, pass `method = FALSE` to `check_ard_structure()` in tests.
**Warning signs:** Note (not error) from check_ard_structure. Downstream tools (gtsummary) may or may not require this.

### Pitfall 5: Pool Object Without pars Structure
**What goes wrong:** `tidy_pool_obj()` depends on `as.data.frame.pool()` from rbmi which expects `$pars` list structure
**Why it happens:** Mock pool objects in tests must match rbmi's internal structure (decision 01-01-D3)
**How to avoid:** Use the `make_mock_pool()` helper from existing tests. Pool objects have: `$pars` (named list of parameter results), `$conf.level`, `$alternative`, `$N`, `$method`
**Warning signs:** "Error in as.data.frame.pool" when testing with simplified mock objects

### Pitfall 6: Overriding print.pool Changes Default Console Experience
**What goes wrong:** Users who rely on rbmi's default print.pool format see different output after loading rbmiUtils
**Why it happens:** S3 method override takes precedence
**How to avoid:** Ensure the enhanced print.pool includes all information from rbmi's version (N combined, method, confidence level, alternative) plus the new features (rounded estimates, formatted CIs, parameter labels). Consider adding a `style` argument or respecting `getOption("rbmiUtils.print_style")`.

## Code Examples

### Example 1: Enhanced print.pool with cli
```r
# Verified pattern using cli formatting
#' @export
print.pool <- function(x, digits = 2, ...) {
  tidy_df <- tidy_pool_obj(x)
  n_params <- nrow(tidy_df)
  visits <- unique(tidy_df$visit)
  n_visits <- length(visits[!is.na(visits)])

  cli::cli_h1("Pool Object")
  cli::cli_text("{n_params} parameter{?s} across {n_visits} visit{?s}")

  # Metadata
  if (!is.null(x$method)) cli::cli_text("{.field Method}: {x$method}")
  if (!is.null(x$N)) cli::cli_text("{.field N imputations}: {x$N}")
  if (!is.null(x$conf.level)) {
    cli::cli_text("{.field Confidence}: {x$conf.level * 100}%")
  }

  cli::cli_rule()

  # Format and display results table
  fmt_df <- format_results_table(tidy_df, est_digits = digits)
  # ... display formatted table

  invisible(x)
}
```

### Example 2: summary.pool with Visit Breakdown
```r
#' @export
summary.pool <- function(object, alpha = 0.05, ...) {
  tidy_df <- tidy_pool_obj(object)

  cli::cli_h1("Pool Object Summary")
  # ... metadata section

  # Visit-level breakdown with significance flags
  for (v in unique(tidy_df$visit)) {
    visit_df <- tidy_df[tidy_df$visit == v, ]
    cli::cli_h2(v)

    trt_rows <- visit_df[visit_df$parameter_type == "trt", ]
    for (i in seq_len(nrow(trt_rows))) {
      sig_flag <- if (!is.na(trt_rows$pval[i]) && trt_rows$pval[i] < alpha) "*" else ""
      cli::cli_text(
        "  {trt_rows$description[i]}: {round(trt_rows$est[i], 2)} ",
        "({round(trt_rows$lci[i], 2)}, {round(trt_rows$uci[i], 2)})",
        " p={format_pvalue(trt_rows$pval[i])} {sig_flag}"
      )
    }
  }

  # Return summary info invisibly
  invisible(list(
    n_parameters = nrow(tidy_df),
    visits = unique(tidy_df$visit),
    method = object$method,
    n_imputations = object$N
  ))
}
```

### Example 3: pool_to_ard() ARD Conversion
```r
# Verified working pattern with cards 0.7.1
pool_to_ard <- function(pool_obj, conf.level = 0.95) {
  if (!requireNamespace("cards", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg cards} is required for ARD conversion. Install with {.code install.packages('cards')}.",
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }

  if (!inherits(pool_obj, "pool")) {
    cli::cli_abort(
      "Input {.arg pool_obj} must be of class {.cls pool}, not {.cls {class(pool_obj)}}.",
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  tidy_df <- tidy_pool_obj(pool_obj)

  stat_names  <- c("estimate", "std.error", "conf.low", "conf.high", "p.value")
  stat_labels <- c("Estimate", "Std. Error",
                    paste0(conf.level * 100, "% CI Lower"),
                    paste0(conf.level * 100, "% CI Upper"),
                    "p-value")
  n_stats <- length(stat_names)

  rows <- lapply(seq_len(nrow(tidy_df)), function(i) {
    r <- tidy_df[i, ]

    base <- data.frame(
      group1       = rep("visit", n_stats),
      group1_level = I(as.list(rep(r$visit, n_stats))),
      group2       = rep("parameter_type", n_stats),
      group2_level = I(as.list(rep(r$parameter_type, n_stats))),
      group3       = rep("lsm_type", n_stats),
      group3_level = I(as.list(rep(
        if (is.na(r$lsm_type)) NA_character_ else r$lsm_type,
        n_stats
      ))),
      variable     = rep(r$parameter, n_stats),
      context      = rep("rbmi_pool", n_stats),
      stat_name    = stat_names,
      stat_label   = stat_labels,
      stat         = I(list(r$est, r$se, r$lci, r$uci, r$pval)),
      fmt_fun      = I(as.list(rep(1L, n_stats))),
      warning      = I(as.list(rep(list(NULL), n_stats))),
      error        = I(as.list(rep(list(NULL), n_stats))),
      stringsAsFactors = FALSE
    )
    base
  })

  ard_df <- do.call(rbind, rows)

  # Add method row (optional but recommended by check_ard_structure)
  method_label <- if (!is.null(pool_obj$method)) pool_obj$method else "unknown"
  # ... add method rows per unique variable grouping

  cards::tidy_ard_column_order(cards::as_card(ard_df))
}
```

### Example 4: Enhanced print.analysis with cli
```r
#' @export
print.analysis <- function(x, ...) {
  cli::cli_h1("Analysis Object")
  cli::cli_text("{length(x$results)} imputation{?s} analysed with {.fn {x$fun_name}}")

  cli::cli_rule()
  cli::cli_text("{.field Method}: {method_class}")
  cli::cli_text("{.field Pooling}: {class(x$results)[1]}")
  cli::cli_text("{.field Delta}: {if (!is.null(x$delta)) 'Applied' else 'None'}")

  # Parameter count and visit info
  if (length(x$results) > 0 && is.list(x$results[[1]])) {
    param_names <- names(x$results[[1]])
    n_params <- length(param_names)
    # Extract visits from parameter names
    visits <- unique(sub("^(trt|lsm)_(ref_|alt_)?", "", param_names))
    cli::cli_text("{.field Parameters}: {n_params}")
    cli::cli_text("{.field Visits}: {paste(visits, collapse = ', ')}")
  }

  cli::cli_text("")
  cli::cli_text("Next: {.code pool_obj <- rbmi::pool(analysis_obj)}")

  invisible(x)
}
```

### Example 5: summary.analysis with Parameter Preview
```r
#' @export
summary.analysis <- function(object, n_preview = 5, ...) {
  # ... detailed summary with parameter table preview
  if (length(object$results) > 0 && is.list(object$results[[1]])) {
    params <- object$results[[1]]
    preview <- head(names(params), n_preview)

    cli::cli_h2("Parameter Preview (from first imputation)")
    for (p in preview) {
      est <- params[[p]]$est
      se  <- params[[p]]$se
      cli::cli_text("  {p}: est={round(est, 3)}, se={round(se, 3)}")
    }
    if (length(params) > n_preview) {
      cli::cli_text("  ... and {length(params) - n_preview} more")
    }
  }

  invisible(summary_info)
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `tidy_as_ard()` for model output | `ard_identity()` for pre-computed stats | cards 0.4.0+ | `tidy_as_ard()` is lifecycle "questioning"; use `ard_identity()` or `as_card()` |
| `cat()` for print methods | `cli::cli_text()` with semantic markup | cli 3.0+ | Colorized, width-aware, pluralization-aware output |
| Custom ARD validation | `cards::check_ard_structure()` | cards 0.3.0+ | Standard validation with column order + method row checks |
| `stat` as atomic vector | `stat` as list column in ARD | cards design | List columns allow mixed types (numeric, character, NULL) per row |

**Deprecated/outdated:**
- `tidy_as_ard()`: Lifecycle "questioning" -- developers considering `ard_summary()` + `ard_formals()` as replacement. Do NOT use for new code.
- `cat()` in print methods: While functional, inconsistent with Phase 1's adoption of cli for error messages. Migrate to cli for consistency.

## Dependency Decision: cards as Suggests vs Imports

**Recommendation: cards as Suggests** (HIGH confidence)

Rationale:
1. `pool_to_ard()` is the ONLY function that uses cards
2. Most users may not need ARD conversion (many just want tidy tibble)
3. cards has its own dependency chain that would increase installation footprint
4. The `requireNamespace()` check pattern is standard CRAN practice
5. This matches the project's constraint from PROJECT.md: "cards, cardx, gtsummary, gt as Suggests vs Imports -- to be decided during planning"
6. cardx and gtsummary (Phase 3) will also be Suggests, keeping consistency

**Impact on code:** Only `pool_to_ard()` needs the `requireNamespace("cards")` guard. Print/summary methods use only cli (already Imports).

## Pool Object Structure Reference

Verified structure of rbmi `pool` class objects (HIGH confidence -- tested directly):

```
pool object (class: "pool"):
$pars       - Named list of parameter results
  $<param_name>  - Named list per parameter
    $est     - numeric: point estimate
    $se      - numeric: standard error
    $ci      - numeric[2]: c(lower, upper) confidence interval
    $pvalue  - numeric: p-value (can be NA)
$conf.level - numeric: confidence level (e.g., 0.95)
$alternative - character: "two.sided"
$N          - integer: number of results combined
$method     - character: pooling method ("rubin", "jackknife", "bootstrap", "bmlmi")
```

**as.data.frame.pool()** converts to tibble with columns: parameter, est, se, lci, uci, pval

**Existing S3 methods from rbmi:** `print.pool`, `as.data.frame.pool`
**Missing from rbmi:** `summary.pool` (does not exist -- rbmiUtils will create it)

## Analysis Object Structure Reference

Verified structure of rbmiUtils `analysis` class objects:

```
analysis object (class: c("analysis", "list")):
$results    - List of analysis results per imputation (class includes pooling method)
  [[i]]     - Named list of parameter estimates from imputation i
    $<param_name>
      $est   - numeric
      $se    - numeric
      $df    - numeric (often NA for Bayesian)
$delta      - data.frame or NULL
$fun        - function (the analysis function used)
$fun_name   - character (name of the analysis function)
$method     - method object (inherits from bayes/approxbayes/condmean/bmlmi)
```

**Existing S3 methods from rbmiUtils:** `print.analysis`, `summary.analysis` (both registered in NAMESPACE)
**Override note:** rbmiUtils overrides rbmi's `print.analysis` (warning at load time is expected)

## ARD Column Structure Reference

Standard ARD (cards) column structure (HIGH confidence -- verified with cards 0.7.1):

```
Required columns (in standard order):
  group1, group1_level    - First grouping variable name and value (list col for _level)
  group2, group2_level    - Second grouping variable name and value
  ...groupN, groupN_level - Additional grouping levels
  variable                - character: variable/parameter name
  variable_level          - list: variable level (often NA for continuous)
  context                 - character: analysis context identifier
  stat_name               - character: statistic name (e.g., "estimate", "p.value")
  stat_label              - character: human-readable label
  stat                    - list: statistic value (MUST be list column)
  fmt_fun                 - list: formatting function or integer digits
  warning                 - list: captured warnings (NULL if none)
  error                   - list: captured errors (NULL if none)
```

**Class:** `c("card", "tbl_df", "tbl", "data.frame")`

**check_ard_structure() validates:**
1. Presence of `warning` and `error` columns
2. Column ordering (when `column_order = TRUE`)
3. Presence of `stat_name = "method"` row (when `method = TRUE`)

## Open Questions

1. **Should print.pool preserve backward compatibility with rbmi's format?**
   - What we know: rbmi's print.pool shows an ASCII table with raw parameter names. The requirements ask for "rounded estimates, formatted CIs, and parameter labels." These are different outputs.
   - What's unclear: Whether users have scripts that parse rbmi's print output (unlikely, but possible).
   - Recommendation: Enhance the output (human-readable, formatted) but include all metadata rbmi's version shows (N, method, conf.level, alternative). This supersedes rbmi's version.

2. **How should lsm_type/arm be represented in ARD grouping?**
   - What we know: ARD uses group1/group2/group3 for grouping. Visit is clearly group1, parameter_type is group2.
   - What's unclear: Whether lsm_type should be group3 (adds NA for trt rows) or a separate variable_level.
   - Recommendation: Use group3 for lsm_type. NA values in group3_level for trt rows are acceptable and match the tidy_pool_obj output. This preserves the visit + parameter_type + arm metadata as required by ARD-03.

3. **Should pool_to_ard() accept a tidy_df directly (in addition to pool_obj)?**
   - What we know: Users might want to convert a pre-tidied data frame to ARD without re-tidying.
   - What's unclear: Whether this adds complexity for little benefit.
   - Recommendation: Accept pool_obj only for v1. The function internally calls tidy_pool_obj(). A future enhancement could accept tidy_df directly.

## Sources

### Primary (HIGH confidence)
- cards package (v0.7.1) -- tested locally: `ard_identity()`, `as_card()`, `check_ard_structure()`, `tidy_ard_column_order()`, column structure
- rbmi package -- tested locally: `print.pool` source, pool object structure, `as.data.frame.pool()`
- rbmiUtils codebase -- read directly: NAMESPACE, DESCRIPTION, analyse_mi_data.R, tidiers.R, test files
- [cards reference index](https://insightsengineering.github.io/cards/main/reference/index.html) -- full function listing
- [cards ard_identity docs](https://insightsengineering.github.io/cards/main/reference/ard_identity.html) -- function signature and examples
- [cards check_ard_structure docs](https://insightsengineering.github.io/cards/main/reference/check_ard_structure.html) -- validation parameters
- [cards as_card docs](https://insightsengineering.github.io/cards/main/reference/as_card.html) -- required columns

### Secondary (MEDIUM confidence)
- [cardx reference index](https://insightsengineering.github.io/cardx/main/reference/index.html) -- naming conventions, ARD conversion patterns from statistical models
- [R Consortium: Supercharging Statistical Analysis with ARDs](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/) -- ecosystem context
- [ARD PHUSE Workshop 2025](https://www.danieldsjoberg.com/ARD-PHUSE-workshop-2025/) -- ARD patterns and best practices

### Tertiary (LOW confidence)
- WebSearch results for pharmaverse ARD patterns -- general ecosystem context only

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - cards package tested locally, cli already in use
- Architecture: HIGH - S3 method patterns verified against existing codebase, ARD structure tested
- Pool/Analysis object structure: HIGH - examined directly via R sessions
- ARD conversion approach: HIGH - verified working prototype with cards 0.7.1
- Pitfalls: HIGH - discovered through direct testing (list columns, method row, S3 override)

**Research date:** 2026-02-08
**Valid until:** 2026-03-08 (cards package is stable; column structure unlikely to change)
