# Architecture Patterns: v3 ARD Enrichment & Polish

**Domain:** Integrating MI-specific diagnostics, enriched ARD metadata, and publication polish into an existing R package for clinical trial analysis
**Researched:** 2026-02-10
**Overall confidence:** HIGH (verified against rbmi source code, cards validation, existing codebase)

---

## Executive Summary

This document defines how three categories of new features -- MI-specific ARD metadata, imputation/draws diagnostic helpers, and publication styling refinements -- integrate into the existing 7-layer rbmiUtils architecture. The key architectural insight is that these features touch three distinct layers with clear boundaries: the ARD conversion layer (enriched metadata), a new diagnostics layer (describe helpers), and the existing reporting layer (styling refinements). None of these features require changes to core analysis or data preparation layers.

The most technically interesting challenge is computing MI diagnostics (FMI, lambda, relative increase in variance) that rbmi's `pool()` computes internally but discards. Two approaches exist: (a) re-derive from the analysis object's per-imputation estimates, or (b) re-derive from the pool object using the pool's `se` and the analysis object's per-imputation `est` values. Both are feasible; approach (a) is recommended because it uses the same mathematical path as rbmi itself.

---

## 1. Current Architecture (Baseline)

### Existing 7-Layer Structure

```
Layer 1: Data Preparation         validate_data(), prepare_data_ice(), summarise_missingness()
Layer 2: Analysis Execution       analyse_mi_data() -> analysis S3 class
Layer 3: Analysis Utilities       gcomp_responder(), gcomp_binary(), gcomp_responder_multi()
Layer 4: Results Processing       tidy_pool_obj(), format_results(), combine_results()
Layer 5: ARD Conversion           pool_to_ard() [requires cards]
Layer 6: Reporting                efficacy_table() [requires gt], plot_forest() [requires ggplot2+patchwork]
Layer 7: Storage Management       reduce_imputed_data(), expand_imputed_data()
```

Plus cross-cutting: S3 print/summary methods for pool and analysis classes, formatting utilities (format_pvalue, format_estimate, format_results_table).

### Key Architectural Properties

1. **Abstraction boundary at tidy tibble.** `tidy_pool_obj()` produces the stable contract between analysis and reporting. All downstream functions (ARD, tables, plots) consume this tibble.
2. **Soft dependencies.** cards, gt, ggplot2, patchwork are in Suggests with runtime guards.
3. **S3 method override pattern.** print.pool and summary.pool are registered by rbmiUtils, overriding rbmi's own methods. This was a deliberate decision (v1) that works because rbmiUtils loads after rbmi.
4. **Mock-based testing.** All reporting tests use hand-crafted mock pool objects rather than running real MCMC imputation.

---

## 2. New Feature Integration Map

### What Changes, What Does Not

| Feature | Layer(s) Affected | New Files | Modified Files | New Dependencies |
|---------|-------------------|-----------|----------------|------------------|
| MI-specific ARD metadata (FMI, lambda, r) | Layer 5 (ARD) | R/mi_diagnostics.R | R/ard_conversion.R | None new |
| describe_draws() | New cross-cutting | R/describe.R | None | None new |
| describe_imputation() | New cross-cutting | R/describe.R | None | None new |
| Efficacy table styling | Layer 6 (Reporting) | None | R/efficacy_table.R | None new |
| Forest plot styling | Layer 6 (Reporting) | None | R/plot_forest.R | None new |
| Documentation/examples | Cross-cutting | vignettes/, man/ | Multiple | None new |

**Critical constraint: No changes to Layers 1-4.** The analysis pipeline, data preparation, and results processing layers are stable and well-tested. All v3 features are additive on top of existing infrastructure.

---

## 3. MI-Specific ARD Metadata: Architecture

### The Problem

rbmi's `pool()` internally computes Rubin's rule components (within-variance, between-variance, lambda, degrees of freedom) but only returns the final pooled estimate, SE, CI, and p-value. The intermediate diagnostics are discarded:

```r
# Inside rbmi:::rubin_rules():
var_w <- mean(ses^2)          # within-imputation variance
var_b <- var(ests)             # between-imputation variance
var_t <- var_w + var_b + var_b/M  # total variance
# lambda and FMI are computed in rubin_df() but not returned
```

**Confidence: HIGH** -- verified by reading rbmi source code directly (`rbmi:::rubin_rules`, `rbmi:::rubin_df`, `rbmi:::pool_internal.rubin`).

### Approach A (Recommended): Re-derive from Analysis Object

The analysis object stores per-imputation parameter estimates and standard errors:

```
analysis_obj$results[[i]]$<parameter_name>$est  -- estimate from imputation i
analysis_obj$results[[i]]$<parameter_name>$se   -- standard error from imputation i
analysis_obj$results[[i]]$<parameter_name>$df   -- degrees of freedom from imputation i
```

We can re-derive all MI diagnostics by applying the same Rubin's rules formulas:

```r
compute_mi_diagnostics <- function(analysis_obj) {
  M <- length(analysis_obj$results)
  param_names <- names(analysis_obj$results[[1]])

  lapply(param_names, function(p) {
    ests <- vapply(analysis_obj$results, function(r) r[[p]]$est, numeric(1))
    ses  <- vapply(analysis_obj$results, function(r) r[[p]]$se, numeric(1))
    dfs  <- vapply(analysis_obj$results, function(r) r[[p]]$df, numeric(1))

    var_w <- mean(ses^2)        # within-imputation variance
    var_b <- var(ests)           # between-imputation variance
    var_t <- var_w + var_b + var_b / M  # total variance

    lambda <- (1 + 1/M) * var_b / var_t    # proportion of total variance from missingness
    fmi <- (var_b + var_b/M) / var_t         # fraction of missing information
    r_increase <- (1 + 1/M) * var_b / var_w  # relative increase in variance due to nonresponse

    list(
      parameter = p,
      m = M,
      var_within = var_w,
      var_between = var_b,
      var_total = var_t,
      lambda = lambda,
      fmi = fmi,
      rel_increase = r_increase,
      df_complete = unique(dfs)
    )
  })
}
```

**Advantages:**
- Uses the exact same mathematical formulation as rbmi
- No dependency on rbmi internals (only accesses documented analysis object structure)
- Works for all Rubin's-rule-based methods (bayes, approxbayes)

**Limitations:**
- Requires the analysis object, not just the pool object
- Does not apply to bootstrap/jackknife/condmean methods (those do not use Rubin's rules)

### Approach B (Alternative): Partial Re-derivation from Pool + Analysis

If only the pool object is available, we can extract the pooled SE (which equals sqrt(var_t)) and re-derive var_b by accessing per-imputation estimates from the analysis object. However, this still requires the analysis object for per-imputation values, making it equivalent to Approach A.

**Recommendation: Use Approach A exclusively.** Functions that need MI diagnostics should accept the analysis object, not the pool object.

### Integration with pool_to_ard()

The enriched `pool_to_ard()` should accept an optional `analysis_obj` parameter. When provided, MI diagnostics are added as additional stat_name rows in the ARD:

```r
pool_to_ard <- function(pool_obj, analysis_obj = NULL, conf.level = NULL) {
  # Existing logic: produce base ARD with estimate, std.error, conf.low, conf.high, p.value
  ard <- existing_pool_to_ard_logic(pool_obj, conf.level)

  # NEW: If analysis_obj provided and method is Rubin's, compute MI diagnostics
  if (!is.null(analysis_obj) && pool_obj$method == "rubin") {
    mi_diag <- compute_mi_diagnostics(analysis_obj)
    diag_rows <- build_ard_rows_from_diagnostics(mi_diag)
    ard <- rbind(ard, diag_rows)
    ard <- cards::tidy_ard_column_order(cards::as_card(ard))
  }

  ard
}
```

### New ARD stat_names for MI Diagnostics

| stat_name | stat_label | Value | When Present |
|-----------|-----------|-------|-------------|
| `"fmi"` | `"Fraction Missing Information"` | numeric [0,1] | Rubin's rule methods only |
| `"lambda"` | `"Proportion Variance from Missingness"` | numeric [0,1] | Rubin's rule methods only |
| `"rel_increase"` | `"Relative Increase in Variance"` | numeric >= 0 | Rubin's rule methods only |
| `"var_within"` | `"Within-Imputation Variance"` | numeric | Rubin's rule methods only |
| `"var_between"` | `"Between-Imputation Variance"` | numeric | Rubin's rule methods only |
| `"m_imputations"` | `"Number of Imputations"` | integer | Always |
| `"df_adjusted"` | `"Barnard-Rubin Degrees of Freedom"` | numeric | Rubin's rule methods only |
| `"method"` | `"Method"` | character | Always (already exists) |

**cards compatibility note:** The cards `check_ard_structure()` validates that a `stat_name == "method"` row exists and checks column types (list columns for stat, fmt_fn, warning, error, and level columns). It does NOT restrict stat_name values -- custom stat_names are fully supported. **Confidence: HIGH** -- verified by reading cards::check_ard_structure source code.

### New File: R/mi_diagnostics.R

This file should contain:

1. `compute_mi_diagnostics(analysis_obj)` -- internal function returning a list of diagnostic values per parameter
2. Helper to validate that the analysis object uses Rubin's rules (check `class(analysis_obj$results)[1] == "rubin"`)
3. Helper to build ARD rows from diagnostics in the same format as pool_to_ard()

---

## 4. describe_draws() and describe_imputation(): Architecture

### The Problem

rbmi's `print.draws()` and `print.imputation()` provide basic console output but are not designed for programmatic use or rich diagnostics. rbmiUtils users need:
- Structured summaries accessible as data (not just printed text)
- MCMC diagnostics for draws objects (convergence, ESS)
- Missingness summaries for imputation objects

### Approach: Standalone Functions (Not S3 Methods)

Per the architecture decision from v1 (Anti-Pattern 4 in the existing ARCHITECTURE.md), rbmiUtils must NOT register S3 print/summary methods for rbmi-owned classes (draws, imputation). Instead, provide standalone functions.

**Note:** rbmiUtils already overrides print.pool and summary.pool. This was a deliberate v1 decision that works because rbmiUtils loads after rbmi and replaces the methods. The same pattern should NOT be extended to draws/imputation because: (a) users may use draws/imputation objects without loading rbmiUtils, (b) the information density of describe_draws() would differ enough from print.draws() to justify a separate function name.

### describe_draws() Design

```r
describe_draws <- function(draws_obj) {
  # Input validation
  if (!inherits(draws_obj, "draws")) {
    cli::cli_abort("...")
  }

  # Extract structured information
  info <- list(
    n_samples = length(draws_obj$samples),
    n_failures = draws_obj$n_failures,
    formula = as.character(draws_obj$formula),
    method = class(draws_obj$method)[[2]],  # bayes, approxbayes, condmean
    method_details = extract_method_details(draws_obj$method),
    has_stan_fit = !is.null(draws_obj$fit)
  )

  # MCMC diagnostics (Bayesian methods only)
  if (!is.null(draws_obj$fit) && inherits(draws_obj$method, "bayes")) {
    info$mcmc <- list(
      # rhat, ESS from Stan fit -- use rstan::summary() if available
    )
  }

  class(info) <- c("rbmiUtils_draws_summary", "list")
  info
}

# Print method for the summary object (we OWN this class)
print.rbmiUtils_draws_summary <- function(x, ...) {
  cli::cli_h1("Draws Summary")
  cli::cli_text("{.field Samples}: {x$n_samples}")
  cli::cli_text("{.field Failures}: {x$n_failures}")
  cli::cli_text("{.field Method}: {x$method}")
  cli::cli_text("{.field Formula}: {x$formula[2]} ~ {x$formula[3]}")
  # ... more detail
  invisible(x)
}
```

### describe_imputation() Design

```r
describe_imputation <- function(impute_obj) {
  # Input validation
  if (!inherits(impute_obj, "imputation")) {
    cli::cli_abort("...")
  }

  # Extract structured information
  n_imp <- length(impute_obj$imputations)

  # Missingness by visit (same computation as rbmi's print.imputation)
  is_miss <- matrix(
    unlist(impute_obj$data$is_missing),
    ncol = length(impute_obj$data$visits),
    byrow = TRUE
  )
  miss_pct <- round((apply(is_miss, 2, sum) / nrow(is_miss)) * 100)

  info <- list(
    n_imputations = n_imp,
    visits = impute_obj$data$visits,
    missing_pct = setNames(miss_pct, impute_obj$data$visits),
    references = impute_obj$references,
    method = class(impute_obj$method)[[2]],
    n_subjects = nrow(is_miss),
    n_visits = ncol(is_miss),
    total_missing_pct = round(sum(is_miss) / length(is_miss) * 100, 1)
  )

  class(info) <- c("rbmiUtils_imputation_summary", "list")
  info
}
```

### Where These Live

Both functions belong in a single new file: **R/describe.R**. This is a cross-cutting diagnostics layer that sits alongside but does not interfere with the existing architecture layers. The file also contains the S3 print methods for the custom summary classes.

### Dependency Impact

- describe_draws() accesses draws_obj$fit (a Stan fit object) for MCMC diagnostics. rstan is already in Suggests for the package. Guard with `requireNamespace("rstan")`.
- describe_imputation() accesses impute_obj$data (an R6 longdata object). This is an rbmi internal but the access pattern is stable (rbmi's own print.imputation uses the same fields).
- Neither function requires any new package dependencies.

---

## 5. Publication Styling Refinements: Architecture

### efficacy_table() Modifications

The existing efficacy_table() in R/efficacy_table.R is self-contained and produces a gt table. Styling refinements are purely additive and do not change the data flow:

| Refinement | Type | Scope |
|-----------|------|-------|
| Font family control | New parameter `font_family` | efficacy_table() signature |
| Row padding/spacing | New parameter `row_padding` or gt theme | Internal gt::tab_options() calls |
| Alternating row shading | New parameter `zebra_stripe` | Internal gt::opt_row_striping() |
| Column alignment fine-tuning | Adjust existing cols_align() calls | Internal |
| Border styling | New theme function `theme_regulatory()` | New internal helper |
| Footnote formatting | Adjust existing tab_source_note() | Internal |

**Pattern: All styling through gt::tab_options().** No structural changes to the data pipeline. The tidy tibble contract is unchanged.

### Recommended: theme_regulatory() Internal Helper

Rather than adding many parameters to efficacy_table(), create an internal theming function:

```r
# Internal -- not exported
theme_regulatory <- function(gt_obj, font_family = "Times New Roman", font_size = 10) {
  gt_obj |>
    gt::tab_options(
      table.font.size = font_size,
      table.font.names = font_family,
      row_group.font.weight = "bold",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(2),
      table_body.border.bottom.style = "solid",
      # ... more regulatory styling
    )
}
```

Then `efficacy_table()` gains a single `theme` parameter:

```r
efficacy_table <- function(pool_obj, ..., theme = "regulatory") {
  # ... existing logic ...
  if (theme == "regulatory") tbl <- theme_regulatory(tbl)
  tbl
}
```

### plot_forest() Modifications

The existing plot_forest() in R/plot_forest.R uses ggplot2 + patchwork. Styling refinements:

| Refinement | Type | Scope |
|-----------|------|-------|
| Font family | New parameter or theme | theme_forest() internal helper |
| Panel spacing | Adjust patchwork widths parameter | Internal layout |
| Axis label formatting | Adjust labs() and scale_x_continuous() | Internal |
| Grid line weight | Adjust theme_forest() | Internal |
| Point/line aesthetics | Existing parameters (point_size, text_size) | Already parameterized |

**Pattern: Extend theme_forest().** The internal `theme_forest()` function already exists. Refinements extend it rather than adding parameters to the main function.

---

## 6. Component Diagram: v3 Architecture

```
Layer 1: Data Preparation         [UNCHANGED]
  validate_data(), prepare_data_ice(), summarise_missingness()
      |
      v
Layer 2: Analysis Execution       [UNCHANGED]
  analyse_mi_data() -> analysis S3 class
      |
      v
Layer 3: Analysis Utilities       [UNCHANGED]
  gcomp_responder(), gcomp_binary(), gcomp_responder_multi()
      |
      v
Layer 4: Results Processing       [UNCHANGED]
  tidy_pool_obj(), format_results(), combine_results()
      |
      v  (abstraction boundary)
      |
Layer 5: ARD Conversion           [ENRICHED]
  pool_to_ard()                   [MODIFIED - optional analysis_obj param]
  compute_mi_diagnostics()        [NEW - internal, R/mi_diagnostics.R]
  Depends: cards (Suggests)
      |
      v
Layer 6: Reporting                [REFINED]
  efficacy_table()                [MODIFIED - styling params]
  plot_forest()                   [MODIFIED - styling params]
  theme_regulatory()              [NEW - internal helper]
  Depends: gt, ggplot2, patchwork (Suggests)

Layer 7: Storage Management       [UNCHANGED]
  reduce_imputed_data(), expand_imputed_data()

Cross-cutting: Diagnostics        [NEW]
  describe_draws()                [NEW - R/describe.R]
  describe_imputation()           [NEW - R/describe.R]
  print.rbmiUtils_draws_summary   [NEW - S3 method for OWN class]
  print.rbmiUtils_imputation_summary [NEW - S3 method for OWN class]
```

### Data Flow for MI Diagnostics

```
rbmi::draws()
      |
      +---> describe_draws()  -----> rbmiUtils_draws_summary (print to console)
      |
      v
rbmi::impute()
      |
      +---> describe_imputation() -> rbmiUtils_imputation_summary (print to console)
      |
      v
analyse_mi_data()
      |
      +---> compute_mi_diagnostics() -> list of FMI/lambda/r per parameter
      |
      v
rbmi::pool()
      |
      v
pool_to_ard(pool_obj, analysis_obj)
      |
      +---> base ARD (existing stats) + MI diagnostic rows (if analysis_obj provided)
      |
      v
efficacy_table() / plot_forest()  (downstream, unchanged data contract)
```

### Data Flow for Styling Refinements

```
pool_obj
    |
    v
tidy_pool_obj()
    |
    +--> efficacy_table(pool_obj, theme = "regulatory")
    |        |
    |        +---> [existing gt pipeline] ---> theme_regulatory() ---> gt_tbl
    |
    +--> plot_forest(pool_obj)
             |
             +---> [existing ggplot2 pipeline] ---> theme_forest() ---> patchwork
```

No changes to data flow. Styling is applied as a final transformation on the output object.

---

## 7. Component Boundaries

| Component | Responsibility | Communicates With | Never Touches |
|-----------|---------------|-------------------|---------------|
| mi_diagnostics.R | Compute FMI/lambda/r from analysis object | analyse_mi_data output | Pool object internals, reporting |
| describe.R | Summarize draws/imputation objects | rbmi draws/impute objects | Analysis results, reporting |
| ard_conversion.R (modified) | Convert pool+diagnostics to ARD | tidy_pool_obj output, mi_diagnostics output | rbmi internals directly |
| efficacy_table.R (modified) | Apply publication styling | tidy_pool_obj output, gt | Analysis objects, rbmi |
| plot_forest.R (modified) | Apply publication styling | tidy_pool_obj output, ggplot2 | Analysis objects, rbmi |

---

## 8. New vs. Modified Files

### New Files

| File | Purpose | Exports | Tests |
|------|---------|---------|-------|
| R/mi_diagnostics.R | MI diagnostic computations (FMI, lambda, r) | `compute_mi_diagnostics()` (or keep internal) | test-mi_diagnostics.R |
| R/describe.R | describe_draws(), describe_imputation() + print methods | `describe_draws()`, `describe_imputation()` | test-describe.R |

### Modified Files

| File | What Changes | Backward Compatible? |
|------|-------------|---------------------|
| R/ard_conversion.R | pool_to_ard() gains optional `analysis_obj` parameter; additional stat_name rows when provided | Yes -- new param has default NULL |
| R/efficacy_table.R | New styling parameters (theme, font_family); internal theme_regulatory() helper | Yes -- new params have defaults |
| R/plot_forest.R | Refinements to theme_forest() and panel spacing | Yes -- visual-only changes |
| NAMESPACE | Export describe_draws, describe_imputation, print methods | Additive only |
| DESCRIPTION | Version bump; no new Imports (rstan already in Suggests) | Yes |
| tests/testthat/ | New test files + updates to existing ARD tests | Additive |
| vignettes/ | New/updated examples showing MI diagnostics | Additive |
| man/ | Roxygen-generated from new/modified functions | Auto-generated |

---

## 9. Suggested Build Order

### Phase ordering rationale

Build in dependency order. Each phase is independently testable.

### Phase 1: MI Diagnostics Core (R/mi_diagnostics.R)

**Build first because:** The ARD enrichment and describe helpers both depend on having the MI diagnostic computation available.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `compute_mi_diagnostics()` | New function computing FMI, lambda, r_increase, var_w, var_b | rbmi analysis object structure |
| Internal helper: validate Rubin method | Check `class(analysis_obj$results)[1] == "rubin"` | rbmi class structure |
| Tests | Mock analysis objects, verify formulas match rbmi | testthat, no optional deps |

**Test strategy:** Create mock analysis objects with known per-imputation estimates/SEs. Verify computed FMI/lambda/r match hand-calculated values. Compare against mice::pool() output for cross-validation.

### Phase 2: describe_draws() and describe_imputation() (R/describe.R)

**Build second because:** These are self-contained diagnostic functions with no downstream consumers in the package. They need MI diagnostics concepts but not the ARD integration.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `describe_draws()` | Extract draws metadata + optional MCMC diagnostics | rbmi draws structure, rstan (Suggests) |
| `describe_imputation()` | Extract imputation metadata + missingness summary | rbmi imputation structure |
| `print.rbmiUtils_draws_summary` | cli-formatted output | cli |
| `print.rbmiUtils_imputation_summary` | cli-formatted output | cli |
| NAMESPACE updates | Export new functions + S3 methods | roxygen2 |
| Tests | Mock draws/imputation objects | testthat |

**Test strategy:** Create mock draws objects (list with samples, n_failures, formula, method, fit=NULL) and mock imputation objects (list with imputations, data, references, method). Test print output with snapshot tests.

### Phase 3: Enriched ARD (modify R/ard_conversion.R)

**Build third because:** Depends on mi_diagnostics from Phase 1. Feeds into downstream reporting.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `pool_to_ard()` modification | Add optional `analysis_obj` parameter | Phase 1 (mi_diagnostics) |
| MI diagnostic ARD rows | Build additional stat_name rows for FMI, lambda, etc. | cards ARD format |
| Backward compatibility | Existing calls without analysis_obj unchanged | None |
| Tests | Verify enriched ARD passes check_ard_structure() | cards (Suggests) |

**Test strategy:** Extend existing test-ard_conversion.R. Add tests with and without analysis_obj. Verify new stat_names appear, values are correct, and check_ard_structure() still passes.

### Phase 4: Publication Styling (modify R/efficacy_table.R, R/plot_forest.R)

**Build fourth because:** Independent of Phases 1-3. Can be built in parallel.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `theme_regulatory()` | New internal helper for gt styling | gt (Suggests) |
| `efficacy_table()` | Add theme parameter, styling refinements | gt (Suggests) |
| `theme_forest()` | Extend with font and spacing options | ggplot2 (Suggests) |
| `plot_forest()` | Refinements to panel spacing and aesthetics | ggplot2, patchwork (Suggests) |
| Tests | Visual regression or snapshot tests | testthat |

**Test strategy:** Extend existing test-efficacy_table.R and test-plot_forest.R. Verify backward compatibility (default parameters produce same output). Add tests for new styling parameters.

### Phase 5: Documentation and Examples

**Build last because:** Depends on all functional code being complete.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| Vignette updates | Add MI diagnostics examples to pipeline vignette | All phases complete |
| README updates | Show FMI output, refined table/plot images | All phases complete |
| roxygen documentation | Document all new/modified functions | All phases complete |
| Pre-rendered images | Regenerate efficacy_table and plot_forest images | All phases complete |
| NEWS.md | Document v3 additions | All phases complete |
| pkgdown site rebuild | Incorporate new reference pages | All phases complete |

### Phase Dependency Graph

```
Phase 1 (MI Diagnostics Core)
    |
    +-------+
    |       |
    v       v
Phase 2   Phase 3
(describe) (ARD enrichment)
    |       |
    |       |
    v       v         Phase 4 (Styling)
    |       |              |
    +---+---+--------------+
        |
        v
    Phase 5 (Documentation)
```

Phases 2, 3, and 4 can run in parallel after Phase 1. Phase 4 has no dependency on Phases 1-3 and can technically start immediately, but is sequenced after Phase 1 for practical focus.

---

## 10. Patterns to Follow

### Pattern 1: Optional Parameter Enrichment

**What:** New parameters default to NULL; when provided, they enrich output without changing the base behavior.

**When:** pool_to_ard() gaining analysis_obj parameter.

```r
pool_to_ard <- function(pool_obj, analysis_obj = NULL, conf.level = NULL) {
  # Base ARD always produced (backward compatible)
  ard <- build_base_ard(pool_obj, conf.level)

  # Enrichment only when explicitly requested
  if (!is.null(analysis_obj)) {
    validate_analysis_obj(analysis_obj)
    diag_rows <- build_mi_diagnostic_rows(analysis_obj)
    ard <- rbind(ard, diag_rows)
  }

  cards::tidy_ard_column_order(cards::as_card(ard))
}
```

### Pattern 2: Own Your Classes

**What:** Create custom S3 classes for objects returned by your functions. Never override print methods for classes you do not own (except pool, which was a deliberate v1 decision).

**When:** describe_draws() and describe_imputation() return objects.

```r
describe_draws <- function(draws_obj) {
  # ... compute summary ...
  info <- list(n_samples = ..., n_failures = ..., ...)
  class(info) <- c("rbmiUtils_draws_summary", "list")
  info
}

# We OWN this class, so S3 method is safe
print.rbmiUtils_draws_summary <- function(x, ...) {
  cli::cli_h1("Draws Summary")
  # ...
}
```

### Pattern 3: Guard Against Non-Rubin Methods

**What:** MI diagnostics (FMI, lambda) only apply to Rubin's rule pooling. Always check the pooling method before computing.

**When:** Any function that computes MI diagnostics.

```r
compute_mi_diagnostics <- function(analysis_obj) {
  pooling_class <- class(analysis_obj$results)[[1]]

  if (pooling_class != "rubin") {
    cli::cli_inform(
      "MI diagnostics (FMI, lambda) are only available for Rubin's rule methods. Skipping."
    )
    return(NULL)
  }
  # ... compute diagnostics ...
}
```

### Pattern 4: Isolated Theming Functions

**What:** Keep visual styling in small, self-contained internal functions rather than embedding in main function bodies.

**When:** efficacy_table() and plot_forest() styling.

```r
# Good: isolated theme function
theme_regulatory <- function(gt_obj, font_family = "Times New Roman") {
  gt_obj |> gt::tab_options(...)
}

# Used in main function
efficacy_table <- function(pool_obj, ..., theme = "regulatory") {
  tbl <- build_table_logic(pool_obj, ...)
  if (theme == "regulatory") tbl <- theme_regulatory(tbl)
  tbl
}
```

---

## 11. Anti-Patterns to Avoid

### Anti-Pattern 1: Computing MI Diagnostics from Pool Object Alone

**What:** Attempting to back-derive FMI from the pool object's SE without per-imputation estimates.

**Why bad:** The pool object only stores sqrt(var_t). You cannot decompose var_t into var_w and var_b without the per-imputation values.

**Instead:** Always require the analysis object for MI diagnostics. Make analysis_obj a required parameter in compute_mi_diagnostics(), and optional in pool_to_ard().

### Anti-Pattern 2: Calling rbmi Internal Functions

**What:** Using `rbmi:::rubin_rules()` or `rbmi:::rubin_df()` directly.

**Why bad:** Triple-colon access to unexported functions is fragile. rbmi can change internal function signatures without notice.

**Instead:** Implement the Rubin's rules formulas directly in rbmiUtils. The formulas are simple, well-documented (Rubin 1987), and do not require access to rbmi internals.

### Anti-Pattern 3: Adding rstan to Imports

**What:** Making rstan a hard dependency for MCMC diagnostics in describe_draws().

**Why bad:** rstan is a heavy dependency (~80MB) with complex system requirements (C++ compiler). Many users will not have it installed.

**Instead:** Keep rstan in Suggests. Guard with `requireNamespace("rstan")`. When rstan is not available, describe_draws() should still return method/sample/failure information, just without ESS/Rhat diagnostics.

### Anti-Pattern 4: Modifying the Tidy Tibble Contract

**What:** Adding FMI/lambda columns to the output of tidy_pool_obj().

**Why bad:** tidy_pool_obj() is the stable abstraction boundary consumed by ALL downstream functions. Adding columns would change the contract and could break efficacy_table(), plot_forest(), and user code.

**Instead:** MI diagnostics are a separate concern. They flow through compute_mi_diagnostics() into pool_to_ard() as additional ARD rows. They do NOT contaminate the tidy tibble.

---

## 12. Scalability Considerations

| Concern | 10 parameters | 100 parameters | 1000 parameters |
|---------|--------------|----------------|-----------------|
| MI diagnostic computation | Instant (<1ms) | Instant (<10ms) | ~100ms (vectorizable) |
| ARD row count (with diagnostics) | ~120 rows | ~1200 rows | ~12000 rows |
| describe_draws() MCMC diagnostics | Instant (summary cached) | Same | Same |
| efficacy_table() rendering | <1s | 2-3s (gt rendering) | May need pagination |

The main scalability concern is gt table rendering for very large pool objects (many visits x parameters). This is a gt limitation, not an rbmiUtils concern. Consider documenting that efficacy_table() is designed for typical clinical trial output (5-20 rows).

---

## Sources

### HIGH confidence (verified against source code and official documentation)

- [rbmi pool.R source](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R) -- rubin_rules, rubin_df implementations
- [rbmi draws.R source](https://github.com/insightsengineering/rbmi/blob/main/R/draws.R) -- draws object structure
- [rbmi impute.R source](https://github.com/insightsengineering/rbmi/blob/main/R/impute.R) -- imputation object structure
- rbmi installed package: `getAnywhere("rubin_rules")`, `getAnywhere("print.draws")`, `getAnywhere("print.imputation")` -- verified in running R session
- cards installed package: `cards::check_ard_structure` source -- verified required columns and stat_name flexibility
- [cards check_ard_structure docs](https://search.r-project.org/CRAN/refmans/cards/html/check_ard_structure.html) -- validation rules
- Existing codebase: R/ard_conversion.R, R/efficacy_table.R, R/plot_forest.R, R/pool_methods.R -- current implementation

### MEDIUM confidence (verified cross-referencing multiple sources)

- [rbmi CRAN documentation](https://cran.r-project.org/web/packages/rbmi/rbmi.pdf) -- pool function reference
- [rbmi statistical specifications](https://cran.r-project.org/web/packages/rbmi/vignettes/stat_specs.html) -- Rubin's rules formulation
- [Rubin's Rules reference (bookdown)](https://bookdown.org/mwheymans/bookmi/rubins-rules.html) -- FMI formulas
- [mice pool documentation](https://amices.org/mice/reference/pool.html) -- cross-validation of FMI/lambda formulas

### Internal (existing codebase)

- `/Users/bailliem/R-projects/rbmiUtils-gsd/.planning/PROJECT.md` -- v3 milestone requirements
- `/Users/bailliem/R-projects/rbmiUtils-gsd/.planning/MILESTONES.md` -- v1/v2 history
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/ard_conversion.R` -- current pool_to_ard()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/efficacy_table.R` -- current efficacy_table()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/plot_forest.R` -- current plot_forest()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/pool_methods.R` -- current print.pool/summary.pool
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/tidiers.R` -- current tidy_pool_obj()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/tests/testthat/test-ard_conversion.R` -- ARD test patterns
- `/Users/bailliem/R-projects/rbmiUtils-gsd/tests/testthat/test-pool_methods.R` -- pool method test patterns

---

*Architecture research: 2026-02-10 (v3 ARD Enrichment & Polish)*
