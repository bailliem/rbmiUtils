# Phase 8: MI Diagnostic Statistics - Research

**Researched:** 2026-02-10
**Domain:** Multiple Imputation diagnostic statistics (Rubin's rules, variance decomposition, ARD enrichment)
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Curated essentials: FMI, lambda, and Barnard-Rubin adjusted df (not full V_w/V_b/V_t breakdown)
- Include relative efficiency (RE = 1 / (1 + FMI/M)) alongside the essentials
- Barnard-Rubin small-sample adjustment for degrees of freedom
- Omit diagnostic rows entirely for non-Rubin pooling methods -- no NA rows, cleaner ARD
- Emit informative message (cli::cli_inform) when diagnostics are omitted due to non-Rubin pooling
- Auto-detect pooling method from pool object structure -- zero user burden, no method= parameter
- No convenience accessor function -- users filter the ARD data frame directly (e.g., dplyr::filter)
- No specific regulatory format alignment required (not CDISC-constrained)
- Degrees of freedom displayed as decimal (e.g., 42.7), not rounded integer

### Claude's Discretion
- stat_name naming convention (short lowercase vs descriptive vs statistical notation)
- Grouping approach (inline with parameters vs separate diagnostic section)
- Whether stat_label column is included
- Per-visit vs summary diagnostics for multi-visit parameters
- Opt-in flag vs automatic inclusion when analysis_obj provided
- FMI version (classical vs adjusted)
- Print layout for enriched ARD
- Decimal precision for FMI/lambda/RE display
- Whether enriched ARD must pass through efficacy_table() without modification

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

This phase enriches `pool_to_ard()` with MI diagnostic statistics by adding a new optional `analysis_obj` parameter. When provided, the function recomputes Rubin's rules variance components from the raw per-imputation estimates and standard errors stored in the analysis object, then appends diagnostic stat rows (FMI, lambda, RIV, Barnard-Rubin df, relative efficiency, and variance components V_w/V_b/V_t per MIDIAG-04) to the ARD for each parameter. When `analysis_obj` is not provided, the function returns the same base ARD as before (backward compatible).

The key technical insight is that rbmi's `pool()` function does NOT store intermediate Rubin's rules quantities (V_w, V_b, lambda, etc.) in the pool object -- it only stores final est/se/ci/pvalue per parameter. The diagnostic statistics must therefore be recomputed from the analysis object, which contains the per-imputation est/se/df vectors needed for Rubin's rules decomposition. This is why the `analysis_obj` parameter is required for MI diagnostics.

The `mice` package (version 3.19.0) provides an authoritative naming convention for MI diagnostics: `ubar` (within variance), `b` (between variance), `t` (total variance), `riv`, `lambda`, `fmi`, `df`, `dfcom`. This research recommends following mice's lowercase naming convention for stat_names in the ARD since it is the de facto standard in R for MI diagnostics.

**Primary recommendation:** Add optional `analysis_obj` parameter to `pool_to_ard()`. When provided and the pooling method is "rubin", compute diagnostics from raw analysis results and append as additional stat_name rows per parameter. When absent or non-Rubin, return the base ARD unchanged (with an informative cli message for non-Rubin).

## Standard Stack

### Core (already in project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| rbmi | >= 1.4 | Source of pool/analysis objects; contains `rubin_rules()` internally | Core dependency |
| cards | (Suggests) | ARD structure, `as_card()`, `check_ard_structure()` | pharmaverse ARD standard |
| cli | >= 3.6.0 | Informative messages for non-Rubin detection | Already used throughout codebase |
| dplyr | (Imports) | ARD data frame manipulation | Already used throughout codebase |

### Supporting (no new dependencies needed)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| stats | base R | `var()` for between-imputation variance | Used in Rubin's rules recomputation |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Recomputing Rubin's rules | Using rbmi:::rubin_rules() internal | Internal functions are not part of public API and could change; recomputation is simple (5 lines) and well-documented |
| Custom variance formulas | mice::pool.scalar() | Would add mice as dependency; the formulas are trivially reimplemented |

**No new dependencies required.** All diagnostic computations use base R arithmetic on vectors already available in the analysis object.

## Architecture Patterns

### Recommended Function Signature Change

The existing `pool_to_ard()` signature is:
```r
pool_to_ard <- function(pool_obj, conf.level = NULL)
```

The new signature should be:
```r
pool_to_ard <- function(pool_obj, analysis_obj = NULL, conf.level = NULL)
```

Adding `analysis_obj` as the second parameter (before `conf.level`) is the natural position since it is the most important new parameter. Making it `NULL` by default preserves backward compatibility.

### Pattern 1: Opt-In Diagnostics via Optional Parameter
**What:** When `analysis_obj` is provided and pooling is Rubin-based, append diagnostic rows. When absent, return base ARD.
**When to use:** Always -- this is the recommended approach.
**Why:** Success criterion 4 requires backward compatibility when `analysis_obj` is not provided. This is cleaner than a boolean flag because it naturally signals intent AND provides the data needed for computation.

```r
pool_to_ard <- function(pool_obj, analysis_obj = NULL, conf.level = NULL) {
  # ... existing base ARD construction ...

  # Append diagnostics if analysis_obj provided
  if (!is.null(analysis_obj)) {
    diag_rows <- compute_mi_diagnostics(pool_obj, analysis_obj)
    if (!is.null(diag_rows)) {
      ard_df <- rbind(ard_df, diag_rows)
    }
  }

  cards::tidy_ard_column_order(cards::as_card(ard_df))
}
```

### Pattern 2: Internal Diagnostic Computation Helper
**What:** Extract a `compute_mi_diagnostics()` internal function that encapsulates Rubin's rules recomputation.
**When to use:** Called by `pool_to_ard()` when analysis_obj is provided and method is Rubin.
**Why:** Separation of concerns -- ARD construction logic stays clean, diagnostic math is isolated and testable.

```r
#' @noRd
compute_mi_diagnostics <- function(pool_obj, analysis_obj) {
  # 1. Check if Rubin-based pooling
  is_rubin <- pool_obj$method == "rubin"
  if (!is_rubin) {
    cli::cli_inform(
      "MI diagnostic statistics are only available for Rubin's rules pooling (method={.val {pool_obj$method}}). Diagnostics omitted."
    )
    return(NULL)
  }

  # 2. Extract parameter names from analysis results
  param_names <- names(analysis_obj$results[[1]])
  M <- length(analysis_obj$results)

  # 3. For each parameter, compute diagnostics
  diag_list <- lapply(param_names, function(pname) {
    ests <- vapply(analysis_obj$results, function(x) x[[pname]]$est, numeric(1))
    ses  <- vapply(analysis_obj$results, function(x) x[[pname]]$se, numeric(1))
    dfs  <- vapply(analysis_obj$results, function(x) x[[pname]]$df, numeric(1))
    v_com <- unique(dfs)

    compute_rubin_diagnostics(ests, ses, v_com, M)
  })
  names(diag_list) <- param_names

  # 4. Build ARD rows from diagnostics
  build_diagnostic_ard_rows(diag_list, pool_obj, analysis_obj)
}
```

### Pattern 3: Rubin's Rules Diagnostic Formulas
**What:** Pure computational function implementing the Rubin's rules variance decomposition.
**Verified against:** rbmi:::rubin_rules source code AND mice::pool.scalar output.

```r
#' @noRd
compute_rubin_diagnostics <- function(ests, ses, v_com, M) {
  est_point <- mean(ests)
  var_w <- mean(ses^2)                    # Within-imputation variance (ubar in mice)
  var_b <- var(ests)                       # Between-imputation variance (b in mice)
  var_t <- var_w + var_b + var_b / M       # Total variance (t in mice)

  # Lambda: proportion of total variance due to missingness
  lambda <- (1 + 1/M) * var_b / var_t

  # RIV: relative increase in variance
  riv <- (1 + 1/M) * var_b / var_w

  # Barnard-Rubin adjusted degrees of freedom
  if (is.na(v_com) || (is.infinite(v_com) && var_b == 0)) {
    df_adj <- Inf
  } else {
    v_old <- (M - 1) / lambda^2
    if (!is.infinite(v_com)) {
      v_obs <- ((v_com + 1) / (v_com + 3)) * v_com * (1 - lambda)
    }
    if (lambda != 0) {
      df_adj <- if (is.infinite(v_com)) v_old else (v_old * v_obs) / (v_old + v_obs)
    } else {
      df_adj <- v_obs
    }
  }

  # FMI (adjusted, following mice convention)
  # fmi = (riv + 2/(df+3)) / (1 + riv)
  fmi <- (riv + 2 / (df_adj + 3)) / (1 + riv)

  # Relative efficiency
  re <- 1 / (1 + fmi / M)

  list(
    var_w = var_w,
    var_b = var_b,
    var_t = var_t,
    lambda = lambda,
    riv = riv,
    df_adj = df_adj,
    dfcom = v_com,
    fmi = fmi,
    re = re
  )
}
```

### Pattern 4: Diagnostic ARD Row Construction (Inline with Parameters)
**What:** Diagnostic stats appear as additional stat_name rows within the same variable grouping.
**Why inline:** The ARD standard expects one variable per parameter, with multiple stat_name rows. Adding diagnostics as additional stat_name rows (alongside "estimate", "std.error", etc.) is the natural ARD pattern. This also ensures diagnostics travel with their parameter when filtering.

```r
# For each parameter, the ARD will contain:
# stat_name = "estimate"     -> existing
# stat_name = "std.error"    -> existing
# stat_name = "conf.low"     -> existing
# stat_name = "conf.high"    -> existing
# stat_name = "p.value"      -> existing
# stat_name = "method"       -> existing
# stat_name = "fmi"          -> NEW diagnostic
# stat_name = "lambda"       -> NEW diagnostic
# stat_name = "riv"          -> NEW diagnostic
# stat_name = "df.adjusted"  -> NEW diagnostic
# stat_name = "df.complete"  -> NEW diagnostic
# stat_name = "var.within"   -> NEW diagnostic (MIDIAG-04)
# stat_name = "var.between"  -> NEW diagnostic (MIDIAG-04)
# stat_name = "var.total"    -> NEW diagnostic (MIDIAG-04)
# stat_name = "re"           -> NEW diagnostic
# stat_name = "m.imputations"-> NEW diagnostic (number of imputations, M)
```

### Anti-Patterns to Avoid
- **Modifying the pool object:** Do NOT try to store diagnostics in the pool object. The pool object comes from rbmi::pool() which we don't control.
- **Calling rbmi internal functions:** Do NOT use `rbmi:::rubin_rules()` or `rbmi:::rubin_df()`. These are internal and could change. The formulas are trivial to reimplement.
- **Adding a separate diagnostic function:** The decision is to enrich pool_to_ard(), not create a separate function. Users filter the ARD directly.
- **Creating NA rows for non-Rubin methods:** CONTEXT.md explicitly says "no NA rows, cleaner ARD" -- omit entirely with informative message.
- **Storing diagnostics in a separate group/variable:** Keep them inline with parameters for natural filtering.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Rubin's rules formulas | Complex statistical library | Simple arithmetic (5 formulas, base R) | The formulas are well-established: var_w = mean(ses^2), var_b = var(ests), etc. No library needed. |
| ARD structure validation | Custom validation | cards::check_ard_structure() | Already used; authoritative validation |
| ARD column ordering | Manual column reorder | cards::tidy_ard_column_order() | Already used; handles all ARD conventions |
| CLI messaging | Custom warning format | cli::cli_inform() | Already used; consistent with codebase pattern |

**Key insight:** The MI diagnostic formulas are simple arithmetic. The complexity is in the ARD integration (getting the data structure right), not the math. Do not over-engineer the computation side; focus on ARD structural correctness and backward compatibility.

## Common Pitfalls

### Pitfall 1: Analysis Object Parameter Names vs Pool Object Parameter Names
**What goes wrong:** The analysis object stores per-imputation results with parameter names like `"trt_Week 24"`. The pool object stores pooled results with the same names in `$pars`. The tidy_pool_obj() function parses these names to extract visit/parameter_type/lsm_type. Diagnostic rows must use the SAME parameter names as the base ARD rows, or filtering will break.
**Why it happens:** If you construct diagnostic rows from the analysis object parameter names but the pool object parameter names differ (e.g., due to sanitization), the ARD will have mismatched variable names.
**How to avoid:** Use the pool object's `names(pool_obj$pars)` as the authoritative parameter names, cross-referenced with analysis object. They should match because `rbmi::pool()` uses the same names. Verify in tests.
**Warning signs:** Diagnostic rows appear as separate "variables" in the ARD instead of additional stat_names within existing parameters.

### Pitfall 2: Cards List Column Requirements
**What goes wrong:** cards::check_ard_structure() requires stat, group_level, variable_level, warning, error, and fmt_fun to be LIST columns. If diagnostic rows are constructed with plain vectors instead of list columns, validation fails.
**Why it happens:** Using `data.frame()` without `I()` wrapping for list columns produces atomic vectors instead of list columns.
**How to avoid:** Use `I(as.list(...))` for all list columns when constructing diagnostic rows, matching the existing pattern in pool_to_ard().
**Warning signs:** `cards::check_ard_structure()` emits warnings about expected list columns.

### Pitfall 3: Backward Compatibility When analysis_obj Is NULL
**What goes wrong:** Adding a new parameter to pool_to_ard() could break existing code if the parameter ordering changes or if existing tests rely on positional arguments.
**Why it happens:** Adding analysis_obj as the second parameter shifts conf.level to position 3.
**How to avoid:** Verify that ALL existing calls to pool_to_ard() use named arguments for conf.level (review existing tests and examples). If any use positional, fix them first. Add analysis_obj with default NULL so omitting it produces identical behavior.
**Warning signs:** Existing tests fail after signature change.

### Pitfall 4: Assuming df Is Always Finite
**What goes wrong:** When complete-data df is Inf (infinite, as in large-sample approximation) or when lambda is exactly 0 (no between-imputation variance), the Barnard-Rubin formula has edge cases. Division by zero or Inf results can propagate.
**Why it happens:** LSM parameters with no real missingness impact can have var_b very close to 0, making lambda near 0.
**How to avoid:** Mirror rbmi's internal `rubin_df()` edge case handling exactly: check for NA v_com, infinite v_com, lambda == 0. Store Inf as Inf in the stat column (not NA).
**Warning signs:** NaN or NA values appearing in diagnostic stats for parameters with very low missingness.

### Pitfall 5: Confusing Classical FMI (lambda) with Adjusted FMI
**What goes wrong:** Using lambda as FMI when they are different quantities. In the MI literature, lambda = (1+1/M)*V_b/V_t is the "proportion of total variance due to missingness," while FMI = (RIV + 2/(df+3))/(1+RIV) is the "adjusted fraction of missing information" that accounts for finite M.
**Why it happens:** Many sources use "FMI" and "lambda" interchangeably. Rubin's original text uses gamma for what mice calls FMI.
**How to avoid:** Include BOTH lambda AND fmi as separate stat_names. Use mice's convention: `lambda` for the unadjusted proportion, `fmi` for the adjusted fraction.
**Warning signs:** Users comparing rbmiUtils FMI values with mice output and getting different numbers.

### Pitfall 6: Non-Rubin Method Detection
**What goes wrong:** The pool object's `$method` field is "rubin", "jackknife", "bootstrap (percentile)", or "bootstrap (normal)". Checking for exact string "rubin" works but could miss edge cases.
**Why it happens:** The method field is set by rbmi::pool() and is a free-text string, not a class-based dispatch.
**How to avoid:** Check `pool_obj$method == "rubin"` for the method field. Also cross-check with `"rubin" %in% class(analysis_obj$results)` if the analysis object is available. Both should agree.
**Warning signs:** Diagnostics being computed for non-Rubin methods, producing meaningless numbers.

### Pitfall 7: efficacy_table() Compatibility
**What goes wrong:** The existing `efficacy_table()` calls `tidy_pool_obj(pool_obj)` directly (not pool_to_ard). But if someone passes an enriched ARD to downstream code expecting only 6 stat_names per variable, things could break.
**Why it happens:** efficacy_table() does NOT use pool_to_ard() -- it calls tidy_pool_obj() directly. So efficacy_table() is naturally unaffected by pool_to_ard() changes.
**How to avoid:** Verify that efficacy_table() continues to work unchanged. The enriched ARD is a superset; downstream consumers that filter by stat_name will naturally skip diagnostic rows they don't recognize.
**Warning signs:** No actual risk here since efficacy_table() doesn't use pool_to_ard(), but document this for the planner.

## Code Examples

### Example 1: Recomputing Rubin's Rules Diagnostics (Verified Against rbmi Internals)

```r
# Source: Verified by running rbmi:::rubin_rules and rbmi:::rubin_df
# against manual computation on ADMI test data (see research session)

# Given: analysis_obj$results[[i]][[param]]$est, $se, $df
param <- "trt_Week 24"
ests <- vapply(analysis_obj$results, function(x) x[[param]]$est, numeric(1))
ses  <- vapply(analysis_obj$results, function(x) x[[param]]$se, numeric(1))
dfs  <- vapply(analysis_obj$results, function(x) x[[param]]$df, numeric(1))
v_com <- unique(dfs)  # Complete-data df (constant across imputations)
M <- length(ests)     # Number of imputations

# Rubin's rules
var_w <- mean(ses^2)                  # Within-imputation variance
var_b <- var(ests)                     # Between-imputation variance
var_t <- var_w + var_b + var_b / M     # Total variance

# Diagnostics
lambda <- (1 + 1/M) * var_b / var_t   # Proportion of variance due to missingness
riv <- (1 + 1/M) * var_b / var_w      # Relative increase in variance

# Barnard-Rubin adjusted df (mirrors rbmi:::rubin_df exactly)
v_old <- (M - 1) / lambda^2
v_obs <- ((v_com + 1) / (v_com + 3)) * v_com * (1 - lambda)
df_adj <- (v_old * v_obs) / (v_old + v_obs)

# Adjusted FMI (mice convention)
fmi <- (riv + 2 / (df_adj + 3)) / (1 + riv)

# Relative efficiency
re <- 1 / (1 + fmi / M)
```

### Example 2: Constructing Diagnostic ARD Rows (Following Existing Pattern)

```r
# Source: Pattern derived from existing pool_to_ard() in R/ard_conversion.R

# Diagnostic stat names and labels
diag_stat_names <- c("fmi", "lambda", "riv", "df.adjusted", "df.complete",
                     "var.within", "var.between", "var.total", "re", "m.imputations")
diag_stat_labels <- c(
  "Fraction of Missing Information",
  "Lambda (Proportion of Variance Due to Missingness)",
  "Relative Increase in Variance",
  "Barnard-Rubin Adjusted df",
  "Complete-Data df",
  "Within-Imputation Variance",
  "Between-Imputation Variance",
  "Total Variance",
  "Relative Efficiency",
  "Number of Imputations"
)

# For one parameter row:
n_diag <- length(diag_stat_names)
diag_values <- list(diag$fmi, diag$lambda, diag$riv, diag$df_adj, diag$dfcom,
                    diag$var_w, diag$var_b, diag$var_t, diag$re, M)

data.frame(
  group1       = rep("visit", n_diag),
  group1_level = I(as.list(rep(r$visit, n_diag))),
  group2       = rep("parameter_type", n_diag),
  group2_level = I(as.list(rep(r$parameter_type, n_diag))),
  group3       = rep("lsm_type", n_diag),
  group3_level = I(as.list(rep(lsm_val, n_diag))),
  variable       = rep(r$parameter, n_diag),
  variable_level = I(as.list(rep(NA, n_diag))),
  context        = rep("rbmi_pool", n_diag),
  stat_name  = diag_stat_names,
  stat_label = diag_stat_labels,
  stat       = I(diag_values),
  fmt_fun    = I(as.list(rep(1L, n_diag))),
  warning    = I(as.list(rep(list(NULL), n_diag))),
  error      = I(as.list(rep(list(NULL), n_diag))),
  stringsAsFactors = FALSE
)
```

### Example 3: Auto-Detecting Rubin Method

```r
# Source: Derived from rbmi::pool() source and as_analysis2() in codebase

# Pool object method field: "rubin", "jackknife", "bootstrap (percentile)", etc.
is_rubin_pool <- pool_obj$method == "rubin"

# Analysis object results class: c("rubin", "list"), c("jackknife", "list"), etc.
is_rubin_analysis <- "rubin" %in% class(analysis_obj$results)

# Use both checks for robustness
is_rubin <- is_rubin_pool && is_rubin_analysis
```

### Example 4: User-Facing API (Backward Compatible)

```r
# Existing usage (unchanged, backward compatible):
ard <- pool_to_ard(pool_obj)

# New enriched usage:
ard_enriched <- pool_to_ard(pool_obj, analysis_obj)

# Filtering diagnostics from enriched ARD:
diag_stats <- c("fmi", "lambda", "riv", "df.adjusted", "re")
ard_enriched |>
  dplyr::filter(stat_name %in% diag_stats)
```

## Discretion Recommendations

Based on research findings, here are recommendations for areas marked as "Claude's Discretion":

### stat_name Naming Convention: Use mice-aligned lowercase with dots
**Recommendation:** Use `fmi`, `lambda`, `riv`, `df.adjusted`, `df.complete`, `var.within`, `var.between`, `var.total`, `re`, `m.imputations`
**Rationale:** The `mice` package uses lowercase single-word names (`fmi`, `lambda`, `riv`). For multi-word stats, use dot-separated notation which aligns with existing stat_names in the ARD (`std.error`, `conf.low`, `conf.high`, `p.value`). This creates a consistent ARD where all stat_names follow the same pattern.

### Grouping Approach: Inline with parameters
**Recommendation:** Append diagnostic stat_name rows within each parameter's existing variable grouping.
**Rationale:** ARD convention is multiple stat_names per variable. Diagnostics are per-parameter quantities. Inline grouping means `dplyr::filter(ard, variable == "trt_Week 24")` returns both the estimate AND its diagnostics. A separate group would require a join to associate diagnostics with parameters.

### stat_label Inclusion: Yes, include descriptive labels
**Recommendation:** Include stat_label for all diagnostic stats.
**Rationale:** The existing ARD includes stat_label for all stats ("Estimate", "Std. Error", etc.). Diagnostic labels provide human-readable context ("Fraction of Missing Information", "Barnard-Rubin Adjusted df"). This is zero-cost and aids downstream formatting.

### Per-Visit vs Summary Diagnostics: Per-parameter (which is per-visit)
**Recommendation:** Compute diagnostics per parameter, which naturally means per-visit since each parameter corresponds to a visit.
**Rationale:** The analysis object stores per-imputation results per parameter. Each parameter (e.g., "trt_Week 24") has its own est/se/df vectors. FMI, lambda, etc. are computed per parameter. Since parameters map to visits, this naturally gives per-visit diagnostics.

### Opt-In vs Always-Include: Opt-in via analysis_obj parameter
**Recommendation:** Include diagnostics only when `analysis_obj` is provided. No boolean flag.
**Rationale:** The analysis_obj IS the data source for diagnostics. Its presence naturally signals intent. A boolean flag without the data source would be confusing. This design makes backward compatibility automatic (omit analysis_obj, get base ARD).

### FMI Version: Adjusted FMI (mice convention)
**Recommendation:** Use the adjusted FMI formula: `fmi = (riv + 2/(df+3)) / (1 + riv)`. Also include lambda as a separate stat for the unadjusted proportion.
**Rationale:** The mice package uses adjusted FMI, which accounts for finite M. Including both fmi and lambda as separate stats provides complete information. This follows the mice convention and the MI literature distinction between gamma (adjusted, what mice calls fmi) and lambda (unadjusted).

### Decimal Precision: 4 decimal places for proportions, full precision for df
**Recommendation:** Use `fmt_fun = 1L` (default formatting) in the ARD. Precision control is a display concern handled by downstream consumers. Store full numeric precision in the stat column.
**Rationale:** ARD is a data container, not a display format. Full precision preserves information. The fmt_fun value of 1L matches the existing convention. Users or downstream tools control display formatting.

### efficacy_table() Compatibility: No modification needed
**Recommendation:** efficacy_table() calls tidy_pool_obj() directly, NOT pool_to_ard(). It is naturally unaffected. The enriched ARD is a superset that works with any consumer that filters by stat_name.
**Rationale:** Verified by reading efficacy_table() source code. It takes pool_obj and calls tidy_pool_obj(), completely bypassing the ARD path.

## Tension: CONTEXT.md vs Requirements on V_w/V_b/V_t

The CONTEXT.md says "Curated essentials: FMI, lambda, and Barnard-Rubin adjusted df (not full V_w/V_b/V_t breakdown)." However, requirement **MIDIAG-04** explicitly states: "User can obtain within/between/total variance components per parameter in ARD output." Success criterion 2 also says: "User can extract within-imputation variance (V_w), between-imputation variance (V_b), total variance (V_t)."

**Recommendation:** Include V_w, V_b, V_t as stat_name rows to satisfy MIDIAG-04 and success criterion 2. The CONTEXT.md "curated essentials" language appears to describe the user-facing emphasis (reviewers care most about FMI/lambda/df/RE), not an exclusion of variance components. Since the requirements explicitly require them and they are trivially computed alongside the other diagnostics, include them.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Store FMI in pool object | Recompute from analysis object | rbmi design (always) | Must pass analysis_obj to get diagnostics |
| Classical FMI only (lambda) | Both lambda AND adjusted FMI | mice >= 3.0 (2018+) | Include both as separate stats |
| Old df formula (Rubin 1987) | Barnard-Rubin adjusted df | Barnard & Rubin 1999; rbmi uses it already | Already standard in rbmi |
| MI diagnostics as separate function | Inline in ARD output | pharmaverse ARD pattern (2024+) | Natural integration with cards ecosystem |

**Deprecated/outdated:**
- Rubin's 1987 original df formula (replaced by Barnard-Rubin 1999 adjustment) -- rbmi already uses the adjustment
- Treating lambda as FMI -- they are different quantities; modern practice distinguishes them

## Data Flow Summary

```
analysis_obj (per-imputation est/se/df)
       |
       v
pool_to_ard(pool_obj, analysis_obj)
       |
       +-- [1] Build base ARD (existing code, unchanged)
       |
       +-- [2] Check: is analysis_obj provided?
       |         NO  -> return base ARD (backward compatible)
       |         YES -> continue
       |
       +-- [3] Check: is pooling method Rubin?
       |         pool_obj$method == "rubin" && "rubin" %in% class(analysis_obj$results)
       |         NO  -> cli::cli_inform("diagnostics omitted"), return base ARD
       |         YES -> continue
       |
       +-- [4] For each parameter in names(analysis_obj$results[[1]]):
       |         - Extract est/se/df vectors across imputations
       |         - Compute: var_w, var_b, var_t, lambda, riv, df_adj, fmi, re
       |         - Build diagnostic ARD rows (matching existing group/variable structure)
       |
       +-- [5] rbind base ARD + diagnostic rows
       |
       +-- [6] cards::tidy_ard_column_order(cards::as_card(combined))
       |
       v
  Enriched ARD (passes check_ard_structure)
```

## Stat Names Complete Specification

| stat_name | stat_label | Value Source | Requirement |
|-----------|-----------|-------------|-------------|
| `fmi` | Fraction of Missing Information | (riv + 2/(df+3)) / (1 + riv) | MIDIAG-01 |
| `riv` | Relative Increase in Variance | (1 + 1/M) * V_b / V_w | MIDIAG-02 |
| `lambda` | Lambda | (1 + 1/M) * V_b / V_t | MIDIAG-03 |
| `var.within` | Within-Imputation Variance | mean(ses^2) | MIDIAG-04 |
| `var.between` | Between-Imputation Variance | var(ests) | MIDIAG-04 |
| `var.total` | Total Variance | V_w + V_b + V_b/M | MIDIAG-04 |
| `df.adjusted` | Barnard-Rubin Adjusted df | Barnard-Rubin formula | MIDIAG-05 |
| `re` | Relative Efficiency | 1 / (1 + fmi/M) | MIDIAG-06 |
| `method` | Method | (already exists) | MIDIAG-07 |
| `df.complete` | Complete-Data df | v_com from analysis results | Supporting |
| `m.imputations` | Number of Imputations | M = length(analysis_obj$results) | Supporting |

## Open Questions

1. **Parameter name alignment between pool and analysis objects**
   - What we know: Both use the same parameter names (verified with ADMI test data). Names come from the analysis function output (e.g., `rbmi::ancova`), and `rbmi::pool()` preserves them via `transpose_results()`.
   - What's unclear: Whether custom analysis functions (e.g., gcomp_responder_multi) could produce names that differ after pooling.
   - Recommendation: Add a validation check that `names(pool_obj$pars)` matches `names(analysis_obj$results[[1]])` and abort with a clear error if they don't.

2. **Edge case: All SEs are NA**
   - What we know: rbmi's `rubin_rules()` handles this by returning `var_t = NA` and `df = NA`. This happens when the analysis function returns NA SEs.
   - What's unclear: How often this occurs in practice.
   - Recommendation: Mirror rbmi's behavior: if all SEs are NA for a parameter, set all diagnostic values to NA for that parameter.

3. **BMLMI method classification**
   - What we know: BMLMI uses a modified pooling approach (not standard Rubin's rules). It has its own variance decomposition (SSW/SSB/MSW/MSB) that differs from Rubin's V_w/V_b.
   - What's unclear: Whether BMLMI diagnostics should be a future extension (DIAG-01 in requirements).
   - Recommendation: Classify BMLMI as non-Rubin for this phase. The requirements doc already defers "BMLMI-specific MI diagnostics" to future work (DIAG-01).

## Sources

### Primary (HIGH confidence)
- rbmi package source code (version >= 1.4): `rubin_rules()`, `rubin_df()`, `pool_internal.rubin()`, `pool()`, `get_pool_components()` -- examined directly via R
- rbmi pool object structure: verified with live ADMI data analysis
- Existing pool_to_ard() source: `/Users/bailliem/R-projects/rbmiUtils-gsd/R/ard_conversion.R`
- Existing analysis_obj structure: verified with `analyse_mi_data()` output
- cards::check_ard_structure() validation rules: examined function source

### Secondary (MEDIUM confidence)
- [mice package documentation](https://amices.org/mice/reference/pool.html) -- naming conventions for MI diagnostics (ubar, b, t, riv, lambda, fmi)
- [mice CRAN manual v3.19.0](https://cran.r-project.org/web/packages/mice/mice.pdf) -- verified December 2025 release
- [Measures of Missing Data Information (bookdown)](https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html) -- FMI formulas, lambda vs adjusted FMI
- [Rubin's Rules (bookdown)](https://bookdown.org/mwheymans/bookmi/rubins-rules.html) -- Barnard-Rubin df formula

### Tertiary (LOW confidence)
- None -- all findings verified against source code or authoritative documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies; all existing packages verified
- Architecture: HIGH -- pattern derived from existing codebase + verified analysis object structure
- Formulas: HIGH -- verified against rbmi source code AND cross-checked with mice output
- Pitfalls: HIGH -- derived from reading actual source code and running test computations
- Naming conventions: HIGH -- mice package naming verified via live R session

**Research date:** 2026-02-10
**Valid until:** 2026-04-10 (stable domain; formulas don't change)
