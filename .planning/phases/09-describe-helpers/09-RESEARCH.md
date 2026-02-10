# Phase 9: Describe Helpers - Research

**Researched:** 2026-02-10
**Domain:** MI pipeline introspection -- structured summaries of rbmi draws and imputation objects
**Confidence:** HIGH

## Summary

Phase 9 adds two user-facing functions -- `describe_draws()` and `describe_imputation()` -- that let users inspect what happened during the MI pipeline without reading raw object internals. These functions return structured S3 objects with informative `print()` methods using `cli` formatting.

The key technical insight is that rbmi's `draws` and `imputation` objects already contain all the information needed. The `draws` object stores: method class, formula, number of samples, number of failures, and (for Bayesian methods only) a `stanfit` object from which MCMC convergence diagnostics (ESS, Rhat) can be extracted via `rstan::summary()`. The `imputation` object stores: an R6 `longdata` object with per-subject missingness flags and group assignments, method details, reference mappings, and the list of imputations. The `describe_*` functions simply extract, structure, and format this information.

The `describe_imputation()` function signature needs careful consideration. The requirements state "describe_imputation(imputed_data, original_data)" which implies it takes data.frames, but the richest information source is the `imputation` object itself (which contains missingness by subject, group assignments, visit structure, method, references, and number of imputations). Taking the `imputation` object directly would be simpler and more informative. Taking data.frames would require a `vars` argument to identify columns. The recommended approach is to accept the `imputation` object from `rbmi::impute()`, optionally with `original_data` for additional context, since this matches the describe_draws pattern of operating on pipeline objects.

**Primary recommendation:** Create two new exported functions (`describe_draws()`, `describe_imputation()`) that take rbmi pipeline objects, extract structured metadata, and return custom S3 classes with `print()` methods using `cli` formatting. For MCMC diagnostics, conditionally use `rstan::summary()` when `rstan` is available (it is already in Suggests).

## Standard Stack

### Core (already in project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| rbmi | >= 1.4 | Source of draws/imputation objects | Core dependency |
| cli | >= 3.6.0 | Formatted console output with `cli_h1`, `cli_text`, `cli_rule`, `cli_alert_*` | Already used for all formatted output in codebase |
| rlang | (Imports) | `%||%` operator, `check_installed()` for optional deps | Already used throughout |

### Supporting (no new dependencies needed)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| rstan | (Suggests) | Extract ESS, Rhat from stanfit objects via `rstan::summary()` | Only when draws used `method_bayes()` and `rstan` is installed |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| rstan::summary() | posterior::summarise_draws() | posterior is not in Suggests; rstan already is and works directly with stanfit |
| Custom S3 class | Plain list | S3 class enables print() dispatch; matches codebase pattern (analysis, pool classes) |
| cli for printing | base cat() | cli is already used everywhere; provides consistent formatting |

**No new dependencies required.** `rstan` is already in Suggests. All other libraries are already in Imports.

## Architecture Patterns

### Recommended File Structure
```
R/
  describe.R             # describe_draws(), describe_imputation(), their S3 print methods
tests/testthat/
  test-describe.R        # Tests for both describe functions
```

A single `R/describe.R` file is appropriate because the two functions are conceptually paired and share patterns. This mirrors how `R/pool_methods.R` groups the `print.pool` and `summary.pool` methods.

### Pattern 1: describe_draws() Structure
**What:** Extract draws metadata into a named list with custom S3 class.
**When to use:** User wants to understand what the imputation model did.

```r
describe_draws <- function(draws_obj) {
  # Validate input
  if (!inherits(draws_obj, "draws")) {
    cli::cli_abort(
      "{.arg draws_obj} must be a {.cls draws} object from {.fn rbmi::draws}.",
      class = c("rbmiUtils_error_type", "rbmiUtils_error")
    )
  }

  # Extract method info
  method_class <- class(draws_obj$method)
  method_name <- switch(method_class[2],
    "bayes" = "Bayesian (MCMC via Stan)",
    "approxbayes" = "Approximate Bayesian",
    "condmean" = "Conditional Mean",
    "unknown"
  )

  # Extract core info
  info <- list(
    method = method_name,
    method_class = method_class[2],
    n_samples = length(draws_obj$samples),
    n_failures = draws_obj$n_failures,
    formula = deparse(draws_obj$formula),
    covariance = draws_obj$method$covariance,
    same_cov = draws_obj$method$same_cov
  )

  # Extract MCMC diagnostics for Bayesian
  if (inherits(draws_obj$method, "bayes") && inherits(draws_obj$fit, "stanfit")) {
    if (requireNamespace("rstan", quietly = TRUE)) {
      stan_summary <- rstan::summary(draws_obj$fit)$summary
      info$mcmc <- list(
        rhat = stan_summary[, "Rhat"],
        ess = stan_summary[, "n_eff"],
        max_rhat = max(stan_summary[, "Rhat"], na.rm = TRUE),
        min_ess = min(stan_summary[, "n_eff"], na.rm = TRUE),
        n_params = nrow(stan_summary),
        converged = all(stan_summary[, "Rhat"] < 1.1, na.rm = TRUE)
      )
    }
  }

  # Method-specific details
  if (inherits(draws_obj$method, "bayes")) {
    info$bayes_control <- list(
      warmup = draws_obj$method$control$warmup,
      thin = draws_obj$method$control$thin,
      chains = draws_obj$method$control$chains,
      seed = draws_obj$method$control$seed
    )
  }
  if (inherits(draws_obj$method, "condmean")) {
    info$condmean_type <- draws_obj$method$type  # "jackknife" or "bootstrap"
  }

  structure(info, class = c("describe_draws", "list"))
}
```

### Pattern 2: describe_imputation() Structure
**What:** Extract imputation metadata including missingness by visit and treatment arm.
**When to use:** User wants to understand what happened during imputation.

```r
describe_imputation <- function(impute_obj, original_data = NULL) {
  # Validate input
  if (!inherits(impute_obj, "imputation")) {
    cli::cli_abort(
      "{.arg impute_obj} must be an {.cls imputation} object from {.fn rbmi::impute}.",
      class = c("rbmiUtils_error_type", "rbmiUtils_error")
    )
  }

  data <- impute_obj$data
  method_class <- class(impute_obj$method)

  # Number of imputations (M)
  n_imp <- length(impute_obj$imputations)

  # Method name
  method_name <- switch(method_class[2],
    "bayes" = "Bayesian (MCMC via Stan)",
    "approxbayes" = "Approximate Bayesian",
    "condmean" = "Conditional Mean",
    "unknown"
  )

  # Missingness from longdata R6 object
  visits <- data$visits
  ids <- data$ids
  groups <- data$group
  is_missing <- data$is_missing

  miss_mat <- matrix(unlist(is_missing), ncol = length(visits), byrow = TRUE)
  colnames(miss_mat) <- visits
  rownames(miss_mat) <- ids

  grp_vec <- vapply(ids, function(id) as.character(groups[[id]]), character(1))

  # Build missingness summary by visit x arm
  missingness <- expand.grid(visit = visits, group = unique(grp_vec),
                             stringsAsFactors = FALSE)
  missingness$n_total <- NA_integer_
  missingness$n_miss <- NA_integer_

  for (i in seq_len(nrow(missingness))) {
    subjs <- names(grp_vec[grp_vec == missingness$group[i]])
    missingness$n_total[i] <- length(subjs)
    missingness$n_miss[i] <- sum(miss_mat[subjs, missingness$visit[i]])
  }
  missingness$pct_miss <- round(100 * missingness$n_miss / missingness$n_total, 1)

  info <- list(
    method = method_name,
    method_class = method_class[2],
    n_imputations = n_imp,
    references = impute_obj$references,
    missingness = missingness,
    visits = visits,
    n_subjects = length(ids)
  )

  structure(info, class = c("describe_imputation", "list"))
}
```

### Pattern 3: S3 print() Methods with cli Formatting
**What:** Informative console output using `cli_h1`, `cli_text`, `cli_rule`, `cli_alert_*`.
**When to use:** Every describe function needs a print method.
**Pattern source:** Existing `print.analysis()` and `print.pool()` in the codebase.

```r
print.describe_draws <- function(x, ...) {
  cli::cli_h1("Draws Summary")
  cli::cli_text("{.field Method}: {x$method}")
  cli::cli_text("{.field Formula}: {x$formula}")
  cli::cli_text("{.field Samples}: {x$n_samples}")
  cli::cli_text("{.field Failures}: {x$n_failures}")
  cli::cli_text("{.field Covariance}: {x$covariance}")
  cli::cli_rule()

  if (!is.null(x$mcmc)) {
    cli::cli_h2("MCMC Convergence")
    if (x$mcmc$converged) {
      cli::cli_alert_success("All Rhat < 1.1 ({x$mcmc$n_params} parameters)")
    } else {
      n_bad <- sum(x$mcmc$rhat >= 1.1, na.rm = TRUE)
      cli::cli_alert_warning("{n_bad} parameter{?s} with Rhat >= 1.1")
    }
    cli::cli_text("{.field Max Rhat}: {round(x$mcmc$max_rhat, 3)}")
    cli::cli_text("{.field Min ESS}: {round(x$mcmc$min_ess, 1)}")
  }

  invisible(x)
}
```

### Anti-Patterns to Avoid
- **Overriding rbmi's print.draws or print.imputation:** Do NOT override these methods. The describe functions are separate, returning new S3 classes. Overriding would create confusion and potential conflicts.
- **Accepting data.frames instead of rbmi objects:** Taking `imputed_data` and `original_data` as data.frames would require a `vars` argument, would lose method/reference information, and would duplicate `summarise_missingness()` logic. Take the rbmi objects directly.
- **Requiring rstan for non-Bayesian draws:** MCMC diagnostics only apply to Bayesian draws. Check `inherits(draws_obj$method, "bayes")` AND `inherits(draws_obj$fit, "stanfit")` before attempting to extract.
- **Making rstan a hard dependency:** Use `requireNamespace("rstan", quietly = TRUE)` for graceful degradation. If rstan is not installed, skip MCMC diagnostics with an informative message.
- **Returning invisible(x) from print:** The print methods should return `invisible(x)` where `x` is the describe object, matching the base R convention and the pattern used in `print.analysis()`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| MCMC diagnostics (ESS, Rhat) | Custom ESS/Rhat computation | `rstan::summary(stanfit)$summary[, c("Rhat", "n_eff")]` | rstan is already in Suggests; computation is non-trivial for ESS |
| cli formatted output | Custom cat() formatting | cli::cli_h1, cli_text, cli_rule, cli_alert_* | Already used throughout codebase; handles width, color, formatting |
| Missingness by visit/arm | Custom aggregation from data.frames | Extract from `impute_obj$data$is_missing` + `$group` | The longdata R6 object already has pre-computed missingness flags |

**Key insight:** The rbmi objects already contain all the information needed. The describe functions are presentation layers, not computation layers. Do not re-derive information that is already stored in the objects.

## Common Pitfalls

### Pitfall 1: Confusing draws_obj$fit with draws_obj$samples
**What goes wrong:** Assuming all draws objects have a `$fit` (stanfit) field. Only `method_bayes()` draws have a non-NULL `$fit`. Condmean and approxbayes both have `$fit == NULL`.
**Why it happens:** The class hierarchy is `c("draws", "condmean", "list")` or `c("draws", "random", "list")`. The "random" class is shared by both bayes and approxbayes, but only bayes has a stanfit.
**How to avoid:** Always check `inherits(draws_obj$fit, "stanfit")` before attempting MCMC diagnostics, not just `inherits(draws_obj$method, "bayes")`.
**Warning signs:** Errors when calling `rstan::summary(NULL)`.

### Pitfall 2: condmean Sample Count Display
**What goes wrong:** The condmean method has 1 primary sample + N jackknife/bootstrap samples. `print.draws` displays this as "1 + N". If `describe_draws` just reports `length(draws_obj$samples)`, it shows "N+1" which is misleading without context.
**Why it happens:** For condmean, `draws_obj$samples[[1]]` is the primary (full-data) fit, and the rest are resampled. The user sees "1 + 10" in rbmi's output.
**How to avoid:** Check if the method is condmean and display sample count accordingly: `n_primary = 1, n_resampled = length(samples) - 1` for condmean vs `n_samples = length(samples)` for Bayesian.
**Warning signs:** User confused about why condmean shows "11 samples" when they specified `type = "jackknife"` with 10 subjects.

### Pitfall 3: longdata R6 Object Access Patterns
**What goes wrong:** The `$data` field in draws/impute objects is an R6 object, not a plain list. Field access works via `$` (e.g., `data$is_missing`) but the fields are not documented in the public API.
**Why it happens:** R6 objects use reference semantics. The `longdata` class is internal to rbmi and not part of the documented public API.
**How to avoid:** Access only the fields used by rbmi's own `print.draws()` and `print.imputation()` methods: `$visits`, `$ids`, `$is_missing`, `$group`, `$vars`, `$strategies`. These are stable because rbmi's own code depends on them.
**Warning signs:** R6 field access errors or unexpected NULL values.

### Pitfall 4: rstan Not Installed
**What goes wrong:** User has rbmi installed without rstan (possible for condmean-only workflows). Calling `rstan::summary()` without checking availability causes an error.
**Why it happens:** rstan is in Suggests, not Imports. It is only required for `method_bayes()`.
**How to avoid:** Use `requireNamespace("rstan", quietly = TRUE)` before any rstan call. If not available, skip MCMC diagnostics and optionally emit `cli::cli_inform()` suggesting installation.
**Warning signs:** `Error in loadNamespace("rstan")` errors.

### Pitfall 5: Naming Collision with Existing summarise_missingness()
**What goes wrong:** `describe_imputation()` computes missingness by visit and arm, which overlaps with the existing `summarise_missingness()` function. Users may be confused about when to use which.
**Why it happens:** Both functions show missingness patterns but from different perspectives: `summarise_missingness()` takes raw data + vars and provides detailed pattern classification (complete/monotone/intermittent), while `describe_imputation()` takes an imputation object and gives a high-level summary.
**How to avoid:** Document clearly that `describe_imputation()` is a quick summary of what happened during imputation, while `summarise_missingness()` is a detailed pre-imputation data exploration tool. Do not duplicate the pattern classification logic.
**Warning signs:** Users filing issues about "two functions that do the same thing."

### Pitfall 6: Registering S3 Methods Correctly
**What goes wrong:** If `print.describe_draws` and `print.describe_imputation` are not properly registered in NAMESPACE, they will not be dispatched.
**Why it happens:** roxygen2 requires `@export` or explicit `S3method()` registration.
**How to avoid:** Use `@export` on the print methods (roxygen2 detects S3 methods automatically from the naming convention `method.class`). Verify with `devtools::document()` that `S3method(print,describe_draws)` appears in NAMESPACE.
**Warning signs:** `print(describe_obj)` shows the raw list instead of formatted output.

## Code Examples

### Example 1: describe_draws() with Bayesian Method

```r
# User workflow
draws_obj <- rbmi::draws(data, data_ice, vars, method_bayes(...))
desc <- describe_draws(draws_obj)

# Console output via print(desc):
#
# -- Draws Summary -------------------------------------------------------
# Method: Bayesian (MCMC via Stan)
# Formula: CHG ~ 1 + TRT + AVISIT + BASE
# Samples: 100
# Failures: 0
# Covariance: us
# --------
# -- MCMC Convergence --
# v All Rhat < 1.1 (15 parameters)
# Max Rhat: 1.002
# Min ESS: 87.3

# Programmatic access
desc$mcmc$max_rhat    # 1.002
desc$mcmc$converged   # TRUE
desc$n_failures       # 0
```

### Example 2: describe_draws() with Conditional Mean Method

```r
draws_obj <- rbmi::draws(data, data_ice, vars, method_condmean(type = "jackknife"))
desc <- describe_draws(draws_obj)

# Console output:
#
# -- Draws Summary -------------------------------------------------------
# Method: Conditional Mean (jackknife)
# Formula: CHG ~ 1 + TRT + AVISIT + BASE
# Samples: 1 + 10 (1 primary + 10 jackknife)
# Failures: 0
# Covariance: us
```

### Example 3: describe_imputation()

```r
impute_obj <- rbmi::impute(draws_obj, references = c("Drug A" = "Placebo", "Placebo" = "Placebo"))
desc <- describe_imputation(impute_obj)

# Console output:
#
# -- Imputation Summary ---------------------------------------------------
# Method: Conditional Mean (jackknife)
# Number of Imputations (M): 11
# Subjects: 20
#
# -- References --
# Drug A  -> Placebo
# Placebo -> Placebo
#
# -- Missingness by Visit and Arm --
# Visit      Arm      N  Missing  %
# Week 4     Placebo  10   0       0.0%
# Week 4     Drug A   10   1      10.0%
# Week 8     Placebo  10   1      10.0%
# Week 8     Drug A   10   0       0.0%
# Week 12    Placebo  10   1      10.0%
# Week 12    Drug A   10   0       0.0%

# Programmatic access
desc$n_imputations    # 11
desc$missingness      # data.frame with visit/group/n_total/n_miss/pct_miss
desc$references       # named character vector
```

### Example 4: Testing Pattern (Following TDD Codebase Convention)

```r
test_that("describe_draws returns structured S3 object", {
  # Setup minimal draws object (mock or real)
  draws_obj <- make_mock_draws()  # helper

  desc <- describe_draws(draws_obj)

  expect_s3_class(desc, "describe_draws")
  expect_true(!is.null(desc$method))
  expect_true(!is.null(desc$n_samples))
  expect_true(!is.null(desc$n_failures))
  expect_true(!is.null(desc$formula))
})

test_that("describe_draws shows MCMC diagnostics for Bayesian", {
  skip_if_not_installed("rstan")
  draws_obj <- make_bayes_draws()  # helper with method_bayes

  desc <- describe_draws(draws_obj)

  expect_true(!is.null(desc$mcmc))
  expect_true(!is.null(desc$mcmc$rhat))
  expect_true(!is.null(desc$mcmc$ess))
  expect_true(is.logical(desc$mcmc$converged))
})

test_that("describe_draws omits MCMC for non-Bayesian", {
  draws_obj <- make_condmean_draws()

  desc <- describe_draws(draws_obj)

  expect_null(desc$mcmc)
})

test_that("describe_imputation returns structured S3 object", {
  impute_obj <- make_mock_impute()

  desc <- describe_imputation(impute_obj)

  expect_s3_class(desc, "describe_imputation")
  expect_true(is.data.frame(desc$missingness))
  expect_true(all(c("visit", "group", "n_total", "n_miss", "pct_miss")
                  %in% names(desc$missingness)))
})
```

## Design Decision: describe_imputation Signature

The requirements state:
> DESC-03: User can call describe_imputation(imputed_data, original_data) to get method, M, missingness by visit/arm

This implies data.frame inputs. However, there is a tension:

| Approach | Pros | Cons |
|----------|------|------|
| `describe_imputation(impute_obj)` | All info available (method, M, references, missingness by arm). No `vars` needed. Consistent with `describe_draws(draws_obj)` pattern. | Requires user to keep the impute_obj around, not just the imputed data.frame |
| `describe_imputation(imputed_data, original_data)` | Works with saved data.frames (common for large datasets). Matches requirement text. | Needs `vars` parameter to identify columns. Cannot determine method name or references from data.frames alone. Would duplicate `summarise_missingness()` logic. |
| `describe_imputation(impute_obj, original_data = NULL)` | Hybrid: uses impute_obj for method/references/missingness. `original_data` is optional for additional context. | Slightly different from requirement text |

**Recommendation:** Accept `impute_obj` as primary input. This provides method, M, references, and per-visit/arm missingness directly from the object. The `original_data` parameter is not needed because the impute object already contains the longdata with complete missingness information. If the requirement insists on data.frame inputs, a `vars` parameter would be mandatory and the function would lose access to method/reference information.

The planner should clarify this with the user or implement the `impute_obj` approach with clear documentation explaining why it takes the object rather than data.frames.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Raw object inspection (`str()`) | Structured describe functions | This phase | Users get formatted, informative summaries |
| rbmi's `print.draws()` via cat() | cli-formatted output with S3 objects | This phase | Consistent with rbmiUtils formatting; programmatic access to metadata |
| No MCMC diagnostics surfaced | ESS/Rhat from stanfit via rstan | This phase | Users can check convergence without knowing Stan internals |
| No per-arm missingness | Missingness by visit AND arm from impute_obj | This phase | Users see treatment-specific missing data patterns |

**Deprecated/outdated:**
- None -- this is new functionality

## Open Questions

1. **describe_imputation signature: impute_obj vs data.frames**
   - What we know: The impute_obj contains all needed information. Data.frames alone lack method/reference info and need a vars argument.
   - What's unclear: Whether the requirement text "describe_imputation(imputed_data, original_data)" is prescriptive about the signature or just describing the conceptual operation.
   - Recommendation: Implement with `impute_obj` parameter. If data.frame support is needed later, it can be added as an overloaded method or separate function.

2. **Should describe functions integrate with the ARD pipeline?**
   - What we know: Phase 8 added MI diagnostics to the ARD. describe functions are for human-readable summaries, not machine-readable data.
   - What's unclear: Whether describe output should also be convertible to ARD format.
   - Recommendation: No ARD integration. Describe functions serve a different purpose (human inspection during development) than ARD (machine-readable analysis results).

3. **Convergence threshold for Rhat**
   - What we know: The Stan community standard is Rhat < 1.01 (as of Stan 2.26+, 2021). Older standard was Rhat < 1.1. rbmi's own warning uses 1.1.
   - What's unclear: Which threshold to use for the `converged` boolean.
   - Recommendation: Use Rhat < 1.1 as the default threshold (matching rbmi's own convention). Include the max_rhat value so users can apply stricter thresholds if needed.

4. **Display of condmean sample count**
   - What we know: rbmi prints "1 + N" for condmean (1 primary + N resampled). This is unique to condmean.
   - What's unclear: Whether describe_draws should mirror this exact format.
   - Recommendation: Yes, mirror it for condmean. Display "1 + N (1 primary + N jackknife/bootstrap)" to match rbmi convention and provide clarity.

5. **Test strategy for Bayesian draws**
   - What we know: Bayesian draws require rstan and take time to run (Stan compilation + MCMC sampling). Tests should use `skip_if_not_installed("rstan")`.
   - What's unclear: Whether to use mock stanfit objects or actual minimal Stan runs in tests.
   - Recommendation: Use both. Mock objects for fast unit tests of the extraction logic. One integration test with a real (minimal) Bayesian draws for end-to-end verification, behind `skip_on_cran()` and `skip_if_not_installed("rstan")`.

## Sources

### Primary (HIGH confidence)
- rbmi package source code (>= 1.4): `print.draws()`, `print.imputation()`, `validate.draws()`, `longdata` R6 class -- examined via `getAnywhere()` and live R session
- rbmi `draws()` output structure: verified with condmean, approxbayes, and bayes methods on test data
- rbmi `impute()` output structure: verified with live data, including `$data$is_missing`, `$data$group`, `$references`
- rstan `summary()` output: verified on actual stanfit object from `method_bayes()` -- columns include `Rhat`, `n_eff`, `mean`, `se_mean`, `sd`, quantiles
- Existing codebase patterns: `print.analysis()`, `print.pool()`, `summary.analysis()`, `summary.pool()` in `R/analyse_mi_data.R` and `R/pool_methods.R`
- cli package usage: verified throughout codebase (`cli_h1`, `cli_h2`, `cli_text`, `cli_rule`, `cli_abort`, `cli_alert_*`)

### Secondary (MEDIUM confidence)
- Stan convergence diagnostics standards: Rhat < 1.01 is modern recommendation (Vehtari et al. 2021), but rbmi itself uses 1.1 threshold in `fit_mcmc()` warning messages
- mice package patterns: `mice::convergence()` function for MI convergence; naming conventions for MI diagnostics

### Tertiary (LOW confidence)
- None -- all findings verified against source code

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies; all existing packages verified via live R session
- Architecture: HIGH -- patterns derived from existing codebase (print.analysis, print.pool) and verified rbmi object structures
- Pitfalls: HIGH -- derived from live testing of all three method types (bayes, approxbayes, condmean) and edge case exploration
- MCMC diagnostics: HIGH -- verified rstan::summary() output on actual stanfit object; confirmed which methods produce stanfit (only bayes)
- Missingness extraction: HIGH -- verified longdata R6 fields match print.imputation() usage; tested group x visit aggregation

**Research date:** 2026-02-10
**Valid until:** 2026-04-10 (stable domain; rbmi object structure is stable across versions)
