# Technology Stack: v3 ARD Enrichment & Polish

**Project:** rbmiUtils v3 milestone
**Researched:** 2026-02-10
**Overall confidence:** HIGH (all library versions verified via CRAN; rbmi internals verified by source inspection; formulas verified against authoritative references)

**Scope:** This document covers ONLY what is needed for the v3 milestone features:
1. MI-specific ARD metadata (FMI, relative increase in variance, pooling method)
2. `describe_draws()` / `describe_imputation()` diagnostics helpers
3. Publication-quality table and plot polish

Existing stack (rbmi, cards, gt, ggplot2, patchwork, cli, lifecycle, dplyr, tidyr, purrr, rlang, assertthat, beeca) is validated and NOT re-researched.

---

## Key Finding: No New Dependencies Required

The v3 milestone can be delivered entirely with the current dependency set. No new packages need to be added to Imports or Suggests. Here is why, feature by feature.

---

## Feature 1: MI-Specific ARD Metadata

### What Needs to Happen

The current `pool_to_ard()` produces 6 stat rows per parameter: estimate, std.error, conf.low, conf.high, p.value, and method. The v3 goal is to add MI-specific pooling diagnostics:

| New stat_name | Definition | Source |
|---------------|-----------|--------|
| `fmi` | Fraction of Missing Information | Rubin (1987) formula 3.1.10 |
| `lambda` | Proportion of total variance due to missingness | Rubin (1987) formula 3.1.7 |
| `riv` | Relative increase in variance due to nonresponse | `(1 + 1/m) * var_b / var_w` |
| `var_within` | Pooled within-imputation variance (Ubar) | `mean(se_i^2)` |
| `var_between` | Between-imputation variance (B) | `var(est_i)` |
| `var_total` | Total variance under Rubin's rules (T) | `var_w + (1 + 1/m) * var_b` |
| `df_rubin` | Adjusted degrees of freedom | Barnard-Rubin (1999) |
| `m` | Number of imputations | `length(est_i)` |
| `pooling_method` | Pooling type label | "rubin", "bootstrap", "jackknife", "bmlmi" |

### Why No New Dependencies

**The computation requires only base R arithmetic.** The Rubin's rules formulas are:

```r
# Given per-imputation estimates and SEs from the analysis object:
M <- length(ests)
var_w <- mean(ses^2)                    # Within-imputation variance (Ubar)
var_b <- var(ests)                       # Between-imputation variance (B)
var_t <- var_w + (1 + 1/M) * var_b      # Total variance (T)
lambda <- (1 + 1/M) * var_b / var_t     # Proportion due to missingness
riv <- (1 + 1/M) * var_b / var_w        # Relative increase in variance
df_old <- (M - 1) / lambda^2            # Degrees of freedom (old formula)
# Barnard-Rubin (1999) small-sample correction:
df_obs <- ((v_com + 1)/(v_com + 3)) * v_com * (1 - lambda)
df <- (df_old * df_obs) / (df_old + df_obs)
fmi <- (riv + 2/(df + 3)) / (1 + riv)   # Fraction of missing information
```

These are the exact formulas used by mice::pool.scalar() (Rubin 1987, van Buuren 2018) and internally by rbmi's own `rubin_rules()` and `rubin_df()`.

**Verification of rbmi internals (HIGH confidence):** Direct source inspection of rbmi 1.6.0 confirms:
- `rbmi:::rubin_rules()` computes `var_w`, `var_b`, `var_t` and `df` but returns only `est_point`, `var_t`, `df`
- `rbmi:::rubin_df()` computes `lambda` internally but does not expose it
- `rbmi:::pool_internal.rubin()` calls `rubin_rules()` and then `parametric_ci()`, which returns only `est`, `ci`, `se`, `pvalue`
- The pool object at `pool_obj$pars[[i]]` contains only: `est`, `ci`, `se`, `pvalue`
- **Conclusion:** FMI/lambda/riv are NOT available from the pool object. They must be recomputed from per-imputation results.

### Data Access Path for Recomputation

The per-imputation estimates and SEs needed for Rubin's diagnostics are available from the **analysis object** (pre-pooling), not the pool object:

```r
# analysis_obj$results is a list of per-imputation named lists
# For rubin pooling, each element has: est, se, df
# Access pattern (verified from rbmi:::transpose_results):
analysis_obj$results[[imp_i]][[param_name]]$est  # estimate from imputation i
analysis_obj$results[[imp_i]][[param_name]]$se   # SE from imputation i
analysis_obj$results[[imp_i]][[param_name]]$df   # df from imputation i (complete-data df)
```

**Integration decision:** The enriched ARD function will need to accept EITHER:
- (a) The `analysis` object (to access per-imputation results for recomputation), OR
- (b) The `pool` object paired with the `analysis` object

Recommend option (a): a new function like `pool_to_ard(pool_obj, analysis_obj = NULL)` that computes MI diagnostics when the analysis object is provided, and falls back to basic ARD (current behavior) when it is not.

### cards ARD Integration

The MI diagnostics slot naturally into the existing ARD structure as additional `stat_name` rows. The `cards::as_card()` function (verified in cards 0.7.1) accepts any data frame with the standard ARD columns. No special cards API is needed -- just append rows with `stat_name = "fmi"`, `stat_name = "lambda"`, etc.

**Verified:** The existing `pool_to_ard()` builds the ARD manually with `data.frame()` + `cards::as_card()` + `cards::tidy_ard_column_order()`. Adding more stat rows is straightforward and requires no cards API changes.

### Technology Required

| What | Package | Version | Already in DESCRIPTION? |
|------|---------|---------|------------------------|
| Base arithmetic for Rubin's rules | base R | any | N/A (base) |
| ARD structure assembly | cards | >= 0.7.0 | YES (Suggests) |
| Per-imputation data access | rbmi | >= 1.4 | YES (Imports) |

**No new dependencies.**

---

## Feature 2: describe_draws() / describe_imputation() Helpers

### What Needs to Happen

New functions that summarize rbmi `draws` and `imputation` objects with cli-formatted output. These are diagnostic helpers for users to understand their imputation pipeline before pooling.

### rbmi Object Structures (verified by source inspection)

**draws object** (class `"draws"`):
- `x$formula` -- model formula
- `x$method` -- method object (class: bayes/approxbayes/condmean)
- `x$samples` -- list of posterior samples or bootstrap samples
- `x$n_failures` -- number of failed draws
- `x$data` -- input data metadata
- rbmi already has `print.draws()` which shows: number of samples, failures, formula, imputation type, method details

**imputation object** (class `"imputation"`):
- `x$references` -- named character vector (arm -> reference mapping)
- `x$data$visits` -- visit names
- `x$data$is_missing` -- list of logical vectors (per-subject missingness by visit)
- `x$imputations` -- list of imputed datasets
- `x$method` -- method object
- rbmi already has `print.imputation()` which shows: number of datasets, missingness percentages by visit, reference mapping

### What describe_draws() Should Add Beyond print.draws()

| Diagnostic | How to Compute | Package Needed |
|------------|---------------|----------------|
| Number of samples (total, successful, failed) | `length(x$samples)`, `x$n_failures` | base R |
| Failure rate (%) | `x$n_failures / (length(x$samples) + x$n_failures) * 100` | base R |
| Model formula (tidy display) | `deparse(x$formula)` | base R |
| Method class and parameters | `class(x$method)`, method slots | base R |
| Number of subjects | From data dimensions | base R |
| Number of visits | From data dimensions | base R |
| Covariate summary | From formula terms | base R |
| cli-formatted output | `cli::cli_h1()`, `cli::cli_text()`, `cli::cli_rule()` | cli (already Imports) |

### What describe_imputation() Should Add Beyond print.imputation()

| Diagnostic | How to Compute | Package Needed |
|------------|---------------|----------------|
| Number of imputed datasets | `length(x$imputations)` | base R |
| Missingness by visit (count and %) | `x$data$is_missing` matrix operations | base R |
| Overall missingness rate | Aggregate from per-visit | base R |
| Reference mapping display | `x$references` | base R |
| Strategy mapping (if available) | From imputation metadata | base R |
| Imputation range summary (min/max/mean of imputed values) | Loop over `x$imputations`, extract outcome | base R + dplyr |
| cli-formatted output with progress indicators | cli | cli (already Imports) |

### Technology Required

| What | Package | Version | Already in DESCRIPTION? |
|------|---------|---------|------------------------|
| Object introspection | base R | any | N/A |
| Formatted console output | cli | >= 3.6.0 | YES (Imports) |
| Data summarization | dplyr | any | YES (Imports) |

**No new dependencies.**

---

## Feature 3: Publication-Quality Table Polish

### What Needs to Happen

Refine `efficacy_table()` for stricter publication standards: better font control, row striping, border control, footnote placement, and decimal alignment.

### gt Features Available (verified in gt 1.3.0)

All styling features needed for publication polish are already in gt:

| Feature | gt Function | Status |
|---------|------------|--------|
| Row group styling (bold headers) | `tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups())` | Available |
| Alternating row colors | `tab_style(style = cell_fill(color = "#F9F9F9"), locations = cells_body(rows = even))` or `opt_row_striping()` | Available |
| Border control | `tab_style(style = cell_borders(...), ...)` | Available |
| Column width control | `cols_width(...)` | Available |
| Font family | `opt_table_font(font = "...")` or `tab_options(table.font.names = "...")` | Available |
| Decimal alignment | `cols_align(align = "right", ...)` combined with `fmt_number(decimals = ...)` | Available |
| Footnote placement | `tab_footnote(footnote = "...", locations = cells_column_labels(...))` | Available |
| Source notes | `tab_source_note(...)` | Already used |
| Table-wide options | `tab_options(...)` | Available |
| Pre-built styling themes | `opt_stylize(style = 1..6, color = "...")` | Available |
| Missing value display | `sub_missing(missing_text = "---")` | Available |

### Recommended Styling Pattern for Clinical Tables

```r
efficacy_table(pool_obj) |>
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups()
  ) |>
  gt::tab_options(
    table.font.size = gt::px(12),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    heading.border.bottom.style = "solid"
  )
```

### Technology Required

| What | Package | Version | Already in DESCRIPTION? |
|------|---------|---------|------------------------|
| Table styling | gt | >= 1.0.0 | YES (Suggests) |

**No new dependencies.**

---

## Feature 4: Publication-Quality Plot Polish

### What Needs to Happen

Refine `plot_forest()` for publication: consistent font sizing, better spacing, alignment control, and optional theming.

### ggplot2 Features Available (verified in ggplot2 4.0.1)

**Important note on ggplot2 4.0:** The major version bump (S3 to S7 migration) does NOT break standard plotting code. rbmiUtils uses only standard geom_*, scale_*, and theme() functions which are fully backward-compatible. No migration is needed.

| Feature | ggplot2 Approach | Status |
|---------|-----------------|--------|
| Font size control | `theme(text = element_text(size = ...))` | Available |
| Font family | `theme(text = element_text(family = "..."))` | Available |
| Panel spacing | `theme(plot.margin = margin(...))` | Already used |
| Axis label formatting | `scale_x_continuous(labels = ...)` | Available |
| Gridline control | `theme(panel.grid.major.x = element_line(...))` | Available |
| Legend positioning | `theme(legend.position = ...)` | Already used |
| Annotation text | `annotate("text", ...)` or `geom_text()` | Already used |
| Color palette | `scale_colour_manual(values = ...)` with Okabe-Ito | Already used |

### patchwork Features Available (verified in patchwork 1.3.2)

| Feature | patchwork Approach | Status |
|---------|-------------------|--------|
| Panel width ratios | `plot_layout(widths = c(...))` | Already used |
| Annotation | `plot_annotation(title = ..., caption = ...)` | Available |
| Theme propagation | `& theme(...)` | Already used |
| Tag panels (A, B, C) | `plot_annotation(tag_levels = "A")` | Available |
| Nested layouts | `(p1 | p2) / p3` or `wrap_plots()` | Available |
| Spacer | `plot_spacer()` | Available |

### Recommended Polish Additions

1. **Configurable theme function:** Expose `theme_forest()` as a public helper or add theme parameters to `plot_forest()` (font_size, font_family, background_color).
2. **Caption support:** Add `caption` parameter that passes to `patchwork::plot_annotation(caption = ...)`.
3. **Panel width customization:** Expose `widths` parameter (currently hardcoded as `c(3, 4, 1.5)` or `c(3, 5)`).

### Technology Required

| What | Package | Version | Already in DESCRIPTION? |
|------|---------|---------|------------------------|
| Plot engine | ggplot2 | >= 3.5.0 | YES (Suggests) |
| Panel composition | patchwork | >= 1.3.0 | YES (Suggests) |

**No new dependencies.**

---

## Complete Dependency Verdict

### Current DESCRIPTION (Unchanged)

```
Imports:
    assertthat,
    beeca,
    cli (>= 3.6.0),
    dplyr,
    lifecycle (>= 1.0.4),
    purrr,
    rbmi (>= 1.4),
    rlang,
    tidyr

Suggests:
    cards,
    ggplot2,
    gt,
    knitr,
    patchwork,
    readr,
    rmarkdown,
    rstan,
    spelling,
    testthat (>= 3.0.0),
    tibble
```

**Zero changes needed.** All v3 features are implementable with the existing dependency set.

### Why NOT to Add New Dependencies

| Candidate | Why NOT to Add |
|-----------|---------------|
| mice | Would add a heavy dependency (30+ transitive deps) just for pool.scalar(). The Rubin's rules formulas are 6 lines of arithmetic. |
| mitml | Same rationale as mice. The `testEstimates()` function is overkill when we only need Ubar, B, T, lambda, FMI. |
| cardx | Not needed. The MI diagnostics are custom statistics not covered by cardx's emmeans/survival ARD helpers. |
| gtsummary | Not needed for v3. The efficacy_table() already uses gt directly. Adding gtsummary would be for a future as_gtsummary() method (deferred to v4+). |
| scales | Would only be used for `scales::percent()` formatting. Can use `sprintf("%.1f%%", x * 100)` instead. |
| forestploter | Not needed. The existing ggplot2+patchwork forest plot is working well and just needs polish, not replacement. |
| gtExtras | Tempting for `gt_theme_*()` helpers, but adds a dependency for what amounts to ~10 lines of `tab_options()` calls. |

---

## Version Compatibility Matrix

All versions verified against CRAN as of 2026-02-10.

| Package | Installed | CRAN Latest | Min Required | Status |
|---------|-----------|-------------|-------------|--------|
| R | 4.4.x | 4.4.3 | >= 4.1 | OK |
| rbmi | 1.6.0 | 1.6.0 (2026-01-23) | >= 1.4 | OK |
| cards | 0.7.1 | 0.7.1 (2025-12-02) | >= 0.4.0 | OK |
| gt | 1.3.0 | 1.3.0 (2026-01-22) | >= 0.10.0 | OK |
| ggplot2 | 4.0.1 | 4.0.2 (2026-02-03) | >= 3.4.0 | OK |
| patchwork | 1.3.2 | 1.3.2 (2025-08-25) | >= 1.0.0 | OK |
| cli | 3.6.5 | 3.6.5 | >= 3.6.0 | OK |
| dplyr | 1.2.0 | 1.2.0 | >= 1.0.0 | OK |

### ggplot2 4.0 Compatibility Note

ggplot2 4.0 (released September 2025) migrated internals from S3 to S7. This does NOT affect rbmiUtils because the package uses only public API functions (geom_*, scale_*, theme*, aes, etc.), not custom Geom/Stat subclasses. No code changes needed.

---

## Architecture: How Features Integrate with Existing Code

### Feature 1: MI ARD Metadata -- Integration Points

```
analyse_mi_data() --> analysis_obj$results  [per-imputation est/se/df]
                          |
                          v
pool_to_ard(pool_obj, analysis_obj)  [enhanced signature]
    |-- tidy_pool_obj(pool_obj)      [existing: basic stats]
    |-- compute_rubin_diagnostics(analysis_obj)  [NEW: FMI/lambda/riv]
    |-- cards::as_card() + cards::tidy_ard_column_order()  [existing]
    v
ARD with 15 stat_names per parameter (was 6)
```

**Key design decision:** The `compute_rubin_diagnostics()` internal helper should handle all pooling methods:
- **Rubin (bayes/approxbayes):** Full diagnostics (FMI, lambda, riv, var_w, var_b, var_t, df)
- **Bootstrap:** Limited diagnostics (var_b only, no FMI -- bootstrap does not use within-variance)
- **Jackknife:** Limited diagnostics (similar to bootstrap)
- **BMLMI:** Rubin-like diagnostics (uses Rubin's rules within and between BMLMI levels)

The pooling method is determined by `class(analysis_obj$results)[[1]]` which is one of: "rubin", "bootstrap", "jackknife", "bmlmi".

### Feature 2: Describe Helpers -- Integration Points

```
rbmi::draws()   --> draws_obj   --> describe_draws(draws_obj)     [NEW]
rbmi::impute()  --> impute_obj  --> describe_imputation(impute_obj) [NEW]
```

These are standalone diagnostic functions. They do NOT modify existing functions. They consume rbmi objects and produce cli-formatted console output (same pattern as `print.pool()` and `summary.pool()` in pool_methods.R).

**Return value pattern:** Follow the `summary.pool()` pattern -- print to console, invisibly return a list of computed diagnostics for programmatic use.

### Feature 3 & 4: Table/Plot Polish -- Integration Points

```
efficacy_table(pool_obj, ...)  [enhanced with style parameters]
    |-- new params: theme, font_size, row_striping, ...
    |-- internally: gt::tab_style() + gt::tab_options()

plot_forest(pool_obj, ...)  [enhanced with polish parameters]
    |-- new params: font_family, widths, caption, ...
    |-- internally: ggplot2::theme() + patchwork::plot_annotation()
```

These are parameter additions to existing functions. No new functions needed for the table/plot polish.

---

## Rubin's Rules Reference Formulas

Included here for implementation reference. All formulas verified against:
- Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
- Barnard, J. and Rubin, D.B. (1999). Small-sample degrees of freedom with multiple imputation. Biometrika, 86(4), 948-955.
- van Buuren, S. (2018). Flexible Imputation of Missing Data, 2nd ed. Chapman & Hall/CRC.
- rbmi 1.6.0 source code (`rbmi:::rubin_rules`, `rbmi:::rubin_df`)
- mice package `pool.scalar()` documentation

### Definitions

| Symbol | Name | Formula |
|--------|------|---------|
| M | Number of imputations | `length(ests)` |
| Qbar | Pooled estimate | `mean(ests)` |
| Ubar | Within-imputation variance | `mean(ses^2)` |
| B | Between-imputation variance | `var(ests)` (note: uses `1/(M-1)` denominator) |
| T | Total variance | `Ubar + (1 + 1/M) * B` |
| r | Relative increase in variance | `(1 + 1/M) * B / Ubar` |
| lambda | Proportion of variance due to missingness | `(1 + 1/M) * B / T` |
| v_old | Old degrees of freedom | `(M - 1) / lambda^2` |
| v_com | Complete-data degrees of freedom | From analysis (typically N - p) |
| v_obs | Observed-data degrees of freedom | `((v_com + 1)/(v_com + 3)) * v_com * (1 - lambda)` |
| v | Adjusted degrees of freedom (Barnard-Rubin) | `(v_old * v_obs) / (v_old + v_obs)` |
| FMI | Fraction of missing information | `(r + 2/(v + 3)) / (1 + r)` |

### Edge Cases

| Condition | Behavior |
|-----------|----------|
| `v_com` is `Inf` (large sample) | Use `v_old` directly (skip Barnard-Rubin correction) |
| `lambda == 0` (no missing data impact) | `df = v_obs` |
| `B == 0` and `v_com` is `Inf` | `df = Inf` |
| All SEs are NA | Return NA for all diagnostics |
| Non-rubin pooling (bootstrap/jackknife) | FMI/lambda/riv not applicable; return NA with informative message |

---

## Confidence Assessment

| Recommendation | Confidence | Basis |
|---------------|------------|-------|
| No new dependencies needed | HIGH | All features verified against existing package APIs; formulas are base R arithmetic |
| Rubin's rules recomputation from analysis_obj | HIGH | rbmi source inspection confirms per-imputation est/se/df stored in analysis_obj$results |
| rbmi pool object lacks FMI/lambda | HIGH | Direct source inspection of rbmi 1.6.0: rubin_rules() discards intermediates |
| cards ARD extension for MI stats | HIGH | cards::as_card() accepts arbitrary stat_name values; verified in cards 0.7.1 |
| gt styling functions for table polish | HIGH | All functions verified present in gt 1.3.0 |
| ggplot2 4.0 compatibility | HIGH | Standard API only; S7 migration does not affect public geom/scale/theme functions |
| describe_draws() feasibility | HIGH | rbmi draws object structure verified by print.draws() source inspection |
| describe_imputation() feasibility | HIGH | rbmi imputation object structure verified by print.imputation() source inspection |

---

## Sources

### CRAN Version Verification (2026-02-10)
- [cards 0.7.1 -- CRAN](https://cran.r-project.org/web/packages/cards/index.html) (Published: 2025-12-02)
- [gt 1.3.0 -- CRAN](https://cran.r-project.org/web/packages/gt/index.html) (Published: 2026-01-22)
- [ggplot2 4.0.2 -- CRAN](https://cran.r-project.org/web/packages/ggplot2/index.html) (Published: 2026-02-03)
- [patchwork 1.3.2 -- CRAN](https://cran.r-project.org/web/packages/patchwork/index.html) (Published: 2025-08-25)
- [rbmi 1.6.0 -- CRAN](https://cran.r-project.org/web/packages/rbmi/index.html) (Published: 2026-01-23)

### rbmi Internals (source inspection)
- `rbmi:::rubin_rules` -- computes var_w, var_b, var_t, df; returns est_point, var_t, df
- `rbmi:::rubin_df` -- computes lambda internally; returns df only
- `rbmi:::pool_internal.rubin` -- calls rubin_rules then parametric_ci; returns est, ci, se, pvalue
- `rbmi:::transpose_results` -- transposes per-imputation results by parameter and component
- `rbmi:::get_pool_components("rubin")` -- returns c("est", "df", "se")
- `rbmi:::print.draws` -- shows formula, method, samples, failures
- `rbmi:::print.imputation` -- shows datasets count, missingness %, references

### Rubin's Rules References
- [Rubin's Rules -- Book of MI (Heymans)](https://bookdown.org/mwheymans/bookmi/rubins-rules.html)
- [Measures of Missing Data Information (Heymans)](https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html)
- [mice::pool documentation](https://amices.org/mice/reference/pool.html)
- [mice::pool.scalar documentation](https://rdrr.io/cran/mice/man/pool.scalar.html)
- [van Buuren: MI in a Nutshell](https://stefvanbuuren.name/fimd/sec-nutshell.html)

### ggplot2 4.0 Compatibility
- [ggplot2 4.0.0 announcement -- Tidyverse blog](https://tidyverse.org/blog/2025/09/ggplot2-4-0-0/)
- [Bioconductor ggplot2 4.0 migration guide](https://blog.bioconductor.org/posts/2025-07-07-ggplot2-update/)
