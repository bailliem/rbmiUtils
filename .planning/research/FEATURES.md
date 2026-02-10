# Feature Landscape

**Domain:** MI-specific ARD metadata, imputation diagnostics, and publication-quality output for clinical trial R packages
**Researched:** 2026-02-10
**Overall confidence:** HIGH (based on codebase inspection, rbmi internals inspection, mice package benchmarking, gt/ggplot2 official docs)

## Context: What rbmiUtils Has Today (v0.1.9, Post-v2 Milestone)

The package has a complete pipeline from rbmi pooled results to formatted regulatory output:

```
rbmi::pool() -> tidy_pool_obj() -> tibble(est, se, lci, uci, pval)
                                      |
                                      +-> pool_to_ard()       -> cards ARD (5 stats + method per param)
                                      +-> efficacy_table()    -> gt regulatory table (LS means, diffs, CIs, p-values)
                                      +-> plot_forest()       -> patchwork 3-panel forest plot (trt/lsm modes)
                                      +-> print.pool()        -> cli formatted console output
                                      +-> summary.pool()      -> visit-level breakdown with significance flags
                                      +-> format_results()    -> formatted strings
                                      +-> extract_trt_effects() / extract_lsm() -> filtered rows
```

### What's Missing for v3

1. **ARD lacks MI-specific metadata** -- pool_to_ard() emits 5 stats (estimate, std.error, conf.low, conf.high, p.value) plus a method row, but does NOT include FMI, relative increase in variance, within/between variance, or Rubin's degrees of freedom
2. **No diagnostics for draws/impute objects** -- users have no programmatic way to inspect imputation setup (number of draws, convergence, method, missingness)
3. **Table/plot styling is functional but not publication-polished** -- font sizing, spacing, alignment need refinement for journal/regulatory submission
4. **README and docs don't show realistic clinical workflow** -- examples use minimal mock data rather than realistic patterns

---

## Table Stakes

Features that are expected for a package claiming MI-aware ARD output and publication-quality reporting. Missing any of these would be a gap given the v3 milestone goals.

| # | Feature | Why Expected | Complexity | Dependencies | Notes |
|---|---------|-------------|------------|--------------|-------|
| TS-1 | **FMI (fraction of missing information) in ARD** | The single most important MI diagnostic statistic. mice::pool() reports it. Every MI textbook covers it. Regulatory reviewers expect it in CSR appendices. Users cannot compute it from pool objects because rbmi discards intermediate Rubin's rules quantities. | Medium | Existing analysis_obj results (per-imputation est, se, df) | rbmi's pool object does NOT store FMI -- it is computed inside `rubin_rules()` then discarded. Must re-derive from raw per-imputation estimates/SEs or extend pool_to_ard() to accept analysis objects alongside pool objects. See Architecture section. |
| TS-2 | **Relative increase in variance (RIV) in ARD** | The second standard MI diagnostic. Quantifies how much the variance increased due to missing data. mice reports it as `riv`. Always reported alongside FMI. | Low | Same computation path as FMI | RIV = (V_b + V_b/M) / V_w. Trivially computed once V_b and V_w are available. |
| TS-3 | **Lambda (proportion of total variance due to missingness) in ARD** | Third standard MI diagnostic. Lambda = (V_b + V_b/M) / V_t. Unadjusted version of FMI. mice reports it. | Low | Same computation path as FMI | Lambda is actually what rbmi computes internally in `rubin_df()` but discards. |
| TS-4 | **Pooling method identification in ARD** | Already partially present (method row in current pool_to_ard). But should be a proper stat_name row with standardized values matching Rubin/bootstrap/jackknife/bmlmi. | Low | Already in pool_obj$method | Current implementation stores method as a stat row. Verify standardization. |
| TS-5 | **describe_draws() -- imputation model summary** | Users need to report what imputation model was used. draws objects contain method type, number of samples, covariance structure, convergence info (for Bayes: ESS, HMC diagnostics). No existing function surfaces this. | Medium | rbmi draws object internals | draws object stores: method (class), samples (list), fit (stanfit or NULL), formula, n_failures. For Bayes: fit contains rstan convergence info. For condmean/approxbayes: fit is model summary. |
| TS-6 | **describe_imputation() -- imputation diagnostics** | Users need to report imputation statistics: method used, number of imputations, number of subjects imputed, visits with missing data, reference pattern. This is table-stakes for CSR appendix material. | Medium | rbmi impute object internals, original data missingness | The impute object has class "imputation" with imputed datasets. Need to extract: imputation method from draws$method, count of imputations from length of datasets, missingness from comparing imputed vs original data. |

---

## Differentiators

Features that elevate rbmiUtils beyond basic MI reporting. Not strictly expected but create significant value.

| # | Feature | Value Proposition | Complexity | Dependencies | Notes |
|---|---------|-------------------|------------|--------------|-------|
| D-1 | **Within/between/total variance components in ARD** | Full Rubin's rules decomposition (V_w, V_b, V_t) as additional ARD stat rows. Not just the derived ratios (FMI, lambda, RIV) but the raw components. Enables regulatory appendix tables showing full pooling diagnostics. mice::pool() reports `ubar` (V_w), `b` (V_b), `t` (V_t). | Low | Same computation as TS-1 | Once we compute V_w, V_b, V_t for FMI, adding them as stat rows is trivial. |
| D-2 | **Adjusted degrees of freedom in ARD** | Barnard-Rubin adjusted df for each parameter. Important for small-sample MI (few imputations or small studies). mice reports `df`. rbmi computes it internally but discards it. | Low | Same computation path | df = (v_old * v_obs) / (v_old + v_obs) where v_old = (M-1)/lambda^2. |
| D-3 | **Relative efficiency (RE) in ARD** | RE = 1 / (1 + FMI/M). Tells users how efficient their number of imputations was. A value of 0.99+ means enough imputations; lower values suggest more are needed. Unique diagnostic for assessing imputation adequacy. | Low | FMI computation from TS-1 | Trivial formula once FMI is computed. Great for automated quality checks. |
| D-4 | **Efficacy table font/spacing refinements** | Publication-quality gt tables need precise control over: font family (e.g., monospace for numbers), row group padding, header weight, column alignment, decimal alignment. Current efficacy_table() uses gt defaults. | Low-Medium | gt tab_options(), tab_style() | Use gt::tab_options() for global settings and gt::tab_style() + gt::cell_text() for targeted formatting. Key params: table.font.size, row_group.font.weight, column_labels.font.weight, data_row.padding. |
| D-5 | **Forest plot font/spacing/alignment refinements** | Publication-quality forest plots need: consistent font sizing across panels, proper left-alignment of visit labels, adequate spacing between rows, customizable point/line sizing, control over panel width ratios. Current plot_forest() works but has alignment imperfections. | Low-Medium | ggplot2 theme(), patchwork layout | Key refinements: text_size parameter already exists; add font_family, panel_widths, row_height parameters. Use ggplot2::element_text(family = ...) and patchwork::plot_layout(widths = ...). |
| D-6 | **README with realistic clinical trial workflow** | Current README shows the pipeline but uses minimal examples. A realistic workflow would show: full ADEFF -> draws -> impute -> analyse -> pool -> table/plot sequence with clinically meaningful data and interpretation notes. | Medium | None (documentation only) | High impact for adoption. Should include both continuous (ANCOVA) and binary (gcomp) endpoints. Show table and forest plot rendered output inline. |
| D-7 | **Documentation examples reflecting real usage** | Current function examples use `\donttest{}` wrappers with comments showing expected usage. Real examples with the package ADMI data would be more helpful, demonstrating actual output. | Medium | Existing ADMI/ADEFF sample data | Must balance example completeness with CRAN check time limits. Use `\donttest{}` for MCMC-based examples, run simpler examples. |

---

## Anti-Features

Features to explicitly NOT build in v3. These are common mistakes or scope-creep risks.

| # | Anti-Feature | Why Avoid | What to Do Instead |
|---|-------------|-----------|-------------------|
| AF-1 | **Re-implementing rbmi::pool() to store FMI** | Modifying rbmi internals is out of scope. rbmiUtils is a utility layer, not a fork. Patching pool() would create maintenance burden and diverge from upstream. | Compute MI diagnostics by re-deriving from per-imputation results stored in the analysis object. The formulas are straightforward (5 lines of R code). Accept analysis_obj as input alongside pool_obj. |
| AF-2 | **Generic imputation diagnostic plots** | Trace plots, density overlays, convergence diagnostics are specialized visualization tasks. The mice package has overplot(), stripplot(), densityplot(). Building equivalents for rbmi would be a major project. | describe_draws() and describe_imputation() should return structured data and print summaries. Users who need diagnostic PLOTS can extract the rstan fit object and use bayesplot/rstan directly. |
| AF-3 | **Automatic quality thresholds for MI diagnostics** | Automated "pass/fail" judgments on FMI values or ESS thresholds are context-dependent. A clinical statistician must interpret these in context. Automated warnings could create false confidence or alarm. | Report the numbers. Document interpretation guidelines in help pages. Let the user decide. |
| AF-4 | **gtsummary tbl_ard_*() integration** | While gtsummary can consume ARDs, its tbl_ard_summary()/tbl_ard_continuous() functions are designed for standard summary statistics, not the bespoke regulatory table layout (LS means by arm + treatment difference + p-value per visit row group). Forcing rbmi output into gtsummary templates would require awkward workarounds. | Keep using gt directly for efficacy_table(). The current approach (gt::gt() with row groups and styled columns) produces better regulatory table output than any gtsummary template. The ARD exists for interoperability, not for driving the table directly. |
| AF-5 | **Supporting non-Rubin pooling methods for MI diagnostics** | Bootstrap and jackknife pooling in rbmi do NOT use Rubin's rules. They have different variance decompositions. FMI/lambda/RIV are Rubin-specific concepts. Trying to compute "FMI equivalents" for bootstrap would be statistically misleading. | MI diagnostic stats (FMI, lambda, RIV) should only be computed when pool_obj$method == "rubin". For bootstrap/jackknife pooling, these fields should be NA with an informative note. Document this clearly. |
| AF-6 | **Sensitivity analysis comparison features** | Deferred to v4+. While valuable, side-by-side sensitivity tables and forest plot overlays are independent features that don't depend on the v3 MI metadata work. Including them would dilute focus. | Keep on v4+ roadmap per PROJECT.md. The v3 ARD enrichment provides the foundation that sensitivity comparison will later build on. |

---

## Feature Dependencies

```
                    Analysis Object (per-imputation results)
                              |
                              v
                    MI Diagnostics Computation
                    (V_w, V_b, V_t, lambda, riv, fmi, df, RE)
                              |
                    +---------+---------+
                    |                   |
                    v                   v
            TS-1-3: FMI/RIV/Lambda   D-1-3: V_w/V_b/df/RE
            in pool_to_ard()         in pool_to_ard()
                    |
                    v
            TS-4: Method standardization
                    (already present, verify)


            Draws Object              Impute Object + Original Data
                |                              |
                v                              v
            TS-5: describe_draws()      TS-6: describe_imputation()
            (method, n_samples,         (method, n_imputations,
             convergence, formula)       % missing, visits imputed)


            Existing efficacy_table()     Existing plot_forest()
                    |                              |
                    v                              v
            D-4: Font/spacing             D-5: Font/spacing/alignment
            refinements via               refinements via
            tab_options()/tab_style()     theme()/patchwork layout


            Existing README + docs
                    |
                    v
            D-6/D-7: Content overhaul
            (realistic examples, workflow)
```

### Critical Path

1. **MI diagnostics computation** is the foundation -- TS-1/2/3 and D-1/2/3 all depend on deriving Rubin's rules intermediate quantities from analysis object results
2. **describe_draws()** and **describe_imputation()** are independent of ARD enrichment -- can be built in parallel
3. **Table/plot refinements** (D-4, D-5) are independent of everything else -- pure styling changes
4. **Documentation overhaul** (D-6, D-7) should come last since it documents the new features

### Dependency on Existing Code

| New Feature | Consumes From Existing Code | Modifies Existing Code |
|-------------|----------------------------|----------------------|
| TS-1/2/3 (FMI/RIV/Lambda in ARD) | `pool_to_ard()` output structure, `tidy_pool_obj()` | **Yes**: pool_to_ard() needs new stat rows; may need new parameter to accept analysis_obj |
| D-1/2/3 (V_w/V_b/df/RE in ARD) | Same as TS-1/2/3 | Same as TS-1/2/3 |
| TS-5 (describe_draws) | rbmi draws object (accessed via analysis_obj$method or standalone) | **No**: new function, no existing code changes |
| TS-6 (describe_imputation) | rbmi impute object, original data, vars | **No**: new function, no existing code changes |
| D-4 (efficacy table styling) | `efficacy_table()` | **Yes**: add tab_options()/tab_style() calls inside efficacy_table(); possibly add style parameters |
| D-5 (forest plot styling) | `plot_forest()` | **Yes**: add theme() refinements, possibly new parameters for font_family, panel_widths |
| D-6/D-7 (README/docs) | All exported functions | **No**: documentation changes only |

---

## Deep Technical Analysis: MI Diagnostic Statistics

### The Core Problem

rbmi's `pool()` function computes MI diagnostics internally but discards them:

```r
# Inside rbmi:::rubin_rules()
var_w <- mean(ses^2)          # within-imputation variance (computed, discarded)
var_b <- var(ests)             # between-imputation variance (computed, discarded)
var_t <- var_w + var_b + var_b/M  # total variance (computed, discarded)

# Inside rbmi:::rubin_df()
lambda <- (1 + 1/M) * var_b / var_t  # proportion due to missingness (computed, discarded)

# Only est, se, ci, pvalue survive into pool_obj$pars[[param]]
```

### The Solution: Re-derive from Analysis Object

The analysis object stores per-imputation results:

```r
ana_obj$results[[i]]$`trt_Week 24`  # list(est = numeric, se = numeric, df = integer)
```

From these, ALL Rubin's rules quantities can be re-derived:

```r
ests <- sapply(ana_obj$results, function(r) r[[param]]$est)
ses  <- sapply(ana_obj$results, function(r) r[[param]]$se)
dfs  <- sapply(ana_obj$results, function(r) r[[param]]$df)
M    <- length(ests)

var_w  <- mean(ses^2)                    # within (mice: ubar)
var_b  <- var(ests)                      # between (mice: b)
var_t  <- var_w + var_b + var_b/M        # total (mice: t)
lambda <- (var_b + var_b/M) / var_t      # proportion of variance due to missingness
riv    <- (var_b + var_b/M) / var_w      # relative increase in variance
v_com  <- unique(dfs)
v_old  <- (M - 1) / lambda^2
v_obs  <- ((v_com + 1)/(v_com + 3)) * v_com * (1 - lambda)
df_adj <- (v_old * v_obs) / (v_old + v_obs)  # Barnard-Rubin adjusted df
fmi    <- (riv + 2/(df_adj + 3)) / (1 + riv)  # fraction of missing information
re     <- 1 / (1 + fmi/M)               # relative efficiency
```

### Verified Computation (HIGH Confidence)

Tested against rbmi internals with ADMI sample data:

```
MI Statistics for trt_Week 24 (M=5):
  var_w (within):  0.032542
  var_b (between): 0.000542
  var_t (total):   0.033192
  lambda:          0.0196
  riv:             0.0200
  fmi:             0.0238
  df_adj:          460.12
```

### Applicability by Pooling Method

| Pooling Method | pool_obj$method | MI Stats Available | Reason |
|---------------|----------------|-------------------|--------|
| Rubin's rules | "rubin" | FMI, lambda, RIV, V_w, V_b, V_t, df, RE | Full Rubin's rules decomposition from per-imputation est/se/df |
| Bootstrap (percentile) | "bootstrap (percentile)" | None | No within-imputation variance; CI from percentiles |
| Bootstrap (normal) | "bootstrap (normal)" | Partial (V_b only) | Between-imputation variance available, but no within |
| Jackknife | "jackknife" | None | Leave-one-out SE, not MI decomposition |
| BMLMI | "bmlmi" | Partial | Uses modified Rubin's rules with D parameter |

**Recommendation:** Implement MI diagnostics only for Rubin's rules pooling (method == "rubin"). Return NA for other methods with an informative attribute/note. This is honest -- FMI is defined only in the Rubin framework.

### Mice Package Benchmark

The mice::pool() output provides these columns per parameter (the gold standard for MI diagnostics):

| mice column | Meaning | Our equivalent | stat_name in ARD |
|------------|---------|----------------|-----------------|
| `estimate` | Pooled estimate | Already in ARD | estimate |
| `std.error` | Pooled SE | Already in ARD | std.error |
| `ubar` | Average within-imputation variance | var_w | var_within |
| `b` | Between-imputation variance | var_b | var_between |
| `t` | Total variance | var_t | var_total |
| `dfcom` | Complete-data degrees of freedom | v_com (from rbmi) | df_complete |
| `df` | Adjusted degrees of freedom | df_adj | df_adjusted |
| `riv` | Relative increase in variance | riv | riv |
| `lambda` | Proportion of variance due to missingness | lambda | lambda |
| `fmi` | Fraction of missing information | fmi | fmi |

**We should match mice's naming conventions** where possible to reduce confusion for users familiar with mice. The stat_name values in ARD should use recognizable names.

---

## Deep Technical Analysis: describe_draws()

### What the Draws Object Contains

```r
# draws object structure (from rbmi:::as_draws)
draws_obj <- list(
  data     = longdata,          # longDataConstructor R6 object
  method   = method,            # method object (class: "method" + "bayes"|"condmean"|etc.)
  samples  = sample_list,       # list of sample_single objects (beta, sigma per imputation)
  fit      = fit,               # For Bayes: stanfit object; for condmean: model fit summary
  n_failures = n_failures,      # integer: number of failed sample draws
  formula  = formula            # model formula used
)
# class: c("draws", "random"|"condmean", "list")
```

### What describe_draws() Should Report

| Information | Source | Available for |
|------------|--------|--------------|
| Imputation method | `class(draws_obj$method)[2]` | All |
| Number of draws/samples | `length(draws_obj$samples)` or `draws_obj$method$n_samples` | All |
| Covariance structure | `draws_obj$method$covariance` | Bayes/approxbayes |
| Same covariance across groups | `draws_obj$method$same_cov` | Bayes/approxbayes |
| Model formula | `draws_obj$formula` | All |
| Number of sample failures | `draws_obj$n_failures` | All |
| MCMC convergence (ESS) | `draws_obj$fit` via rstan | Bayes only |
| MCMC convergence (divergences) | `draws_obj$fit` via rstan | Bayes only |
| MCMC warmup/thin/chains | `draws_obj$method$control` | Bayes only |
| Bootstrap type | `draws_obj$method$type` | condmean only |

### Access Pattern Challenge

**Critical finding:** rbmiUtils does NOT have direct access to draws objects in its current workflow. The pipeline is:

```
rbmi::draws() -> rbmi::impute() -> get_imputed_data() -> analyse_mi_data() -> rbmi::pool()
```

The draws object is consumed by `rbmi::impute()` and not stored by rbmiUtils. The analysis object stores `method` but NOT the draws object or its fit.

**Options:**
1. **describe_draws() takes a draws object directly** -- user must retain it from the rbmi pipeline
2. **describe_draws() takes an analysis object** -- can only report method info, NOT convergence/fit

**Recommendation:** Option 1. Accept a draws object. This is the natural API: users produce draws objects, may want to inspect them before proceeding. The function signature would be `describe_draws(draws_obj)`. Document that users should inspect draws before calling impute/analyse.

---

## Deep Technical Analysis: describe_imputation()

### What Information Is Available

The impute object from rbmi has class "imputation" and contains the imputed datasets. However, rbmiUtils typically works with the stacked ADMI data (after `get_imputed_data()`), not the raw impute object.

| Information | Source | Notes |
|------------|--------|-------|
| Number of imputations | `length(unique(data$IMPID))` | From ADMI stacked data |
| Method used | From method object or analysis_obj$method | Available |
| Subjects with imputed values | Compare imputed vs original data | Requires original data reference |
| Visits with missing data | Compare imputed vs original data | Requires original data reference |
| Percentage missing overall | Original data missingness | Requires original data |
| Percentage missing by visit | Original data missingness by visit | Requires original data |
| Percentage missing by arm | Original data missingness by arm | Requires original data |
| ICE strategies applied | From data_ice if available | May not be retained |

### Design Options

1. **describe_imputation(imputed_data, original_data, vars)** -- Compare stacked imputed data against original to derive missingness summaries
2. **describe_imputation(impute_obj)** -- Take the rbmi impute object directly
3. **describe_imputation(analysis_obj, original_data, vars)** -- Extract from analysis object

**Recommendation:** Option 1. This matches the existing rbmiUtils pattern where functions work with ADMI-style stacked data frames (same signature pattern as `reduce_imputed_data()`). Users already have `imputed_data` and `original_data` in their workflow.

### Relationship to Existing Functions

rbmiUtils already has `summarise_missingness()` in the data preparation layer. describe_imputation() should complement it:

- `summarise_missingness(data, vars)` -- Pre-imputation: "what's missing in my data?"
- `describe_imputation(imputed_data, original_data, vars)` -- Post-imputation: "what was imputed and how?"

---

## Deep Technical Analysis: Publication-Quality Styling

### Efficacy Table (gt) Refinements Needed

Based on review of efficacy_table() and gt best practices:

| Current State | Refinement Needed | gt API | Complexity |
|--------------|-------------------|--------|------------|
| Default gt font | Specify font family for reproducibility | `tab_options(table.font.names = ...)` | Low |
| Default font size | Smaller size for dense regulatory tables (typically 9-10pt) | `tab_options(table.font.size = px(10))` | Low |
| Default row padding | Reduce for compact layout | `tab_options(data_row.padding = px(4))` | Low |
| Row group bold by default | Ensure row groups are visually distinct | `tab_options(row_group.font.weight = "bold")` | Low |
| Number alignment | Center-aligned currently; consider right-alignment for numbers | `cols_align(align = "right", columns = c("est", "se"))` | Low |
| Column width | Uncontrolled; CI column may be too wide | `cols_width()` | Low |
| Source notes styling | Smaller font for footnotes | `tab_options(source_notes.font.size = px(8))` | Low |
| Border styling | Subtle borders between row groups | `tab_style()` with `cell_borders()` | Low |

**Approach:** Add optional `style` parameter to efficacy_table() accepting a named list of overrides, with sensible publication defaults. Or provide a `theme_regulatory()` helper that returns a list of gt tab_options.

### Forest Plot (ggplot2/patchwork) Refinements Needed

| Current State | Refinement Needed | ggplot2 API | Complexity |
|--------------|-------------------|-------------|------------|
| text_size = 3 default | May be too small for some contexts; need scalable | Already parameterized | Low |
| No font_family control | Journals may require specific fonts | `theme(text = element_text(family = ...))` | Low |
| Fixed panel width ratio (3:4:1.5) | Allow customization | `patchwork::plot_layout(widths = ...)` | Low |
| Inconsistent left panel alignment | Visit labels and estimate text may not align perfectly | `hjust` consistency, `scale_x_continuous` limits | Medium |
| No horizontal gridlines in forest panel | Some journals prefer subtle horizontal guides | `panel.grid.major.y = element_line(...)` option | Low |
| plot_forest returns patchwork | Users can customize with `& theme()` but this is not documented prominently | Documentation improvement | Low |
| No option for landscape/portrait orientation | Regulatory documents may need specific aspect ratios | Not a parameter -- user controls with ggsave() | Documentation only |

**Approach:** Add `font_family` and `panel_widths` parameters to plot_forest(). Default panel_widths should remain c(3, 4, 1.5) or c(3, 5) for show_pvalues=FALSE. Alignment fixes are internal implementation changes, not new API.

---

## MVP Recommendation for v3

### Phase 1: MI Diagnostic Statistics (core value)

Priority order within this phase:

1. **Internal helper: compute_rubin_diagnostics(analysis_obj, param)** -- computes all MI stats for one parameter
2. **TS-1/2/3 + D-1/2/3: Enrich pool_to_ard()** -- add FMI, lambda, RIV, V_w, V_b, V_t, df_adjusted, RE as new stat rows
3. **TS-4: Verify method standardization** in existing ARD method row

This is the highest-value work because it fills a real gap: rbmi users currently cannot report standard MI diagnostics in a structured format.

### Phase 2: Describe Helpers (new functions)

4. **TS-5: describe_draws(draws_obj)** -- structured summary of imputation model
5. **TS-6: describe_imputation(imputed_data, original_data, vars)** -- structured summary of what was imputed

These are moderate complexity but high documentation value -- they make the package self-documenting.

### Phase 3: Visual Polish (refinements)

6. **D-4: Efficacy table styling** -- add publication defaults to efficacy_table()
7. **D-5: Forest plot alignment/fonts** -- add font_family, panel_widths; fix alignment

Lower complexity, moderate user-facing impact.

### Phase 4: Documentation (last)

8. **D-6: README overhaul** -- realistic workflow examples
9. **D-7: Function documentation examples** -- real usage patterns

Must come after features are finalized.

### Defer to v4+

- Sensitivity analysis comparison table (D-2 from v1 features list)
- Forest plot with sensitivity overlay
- Responder bar chart
- as_gt() / as_gtsummary() S3 methods
- BMLMI-specific MI diagnostics

---

## Complexity Estimates

| Feature | Lines of R Code (est.) | Test Lines (est.) | Documentation | Total Effort |
|---------|----------------------|-------------------|---------------|-------------|
| Internal MI diagnostics helper | 40-60 | 60-80 | Roxygen (internal) | Low-Medium |
| Enrich pool_to_ard() with MI stats | 60-100 | 80-120 | Update existing roxygen + vignette | Medium |
| describe_draws() | 80-120 | 60-80 | Full roxygen + vignette section | Medium |
| describe_imputation() | 80-120 | 60-100 | Full roxygen + vignette section | Medium |
| Efficacy table styling | 30-50 | 20-40 | Update existing roxygen | Low |
| Forest plot refinements | 30-50 | 20-40 | Update existing roxygen | Low |
| README overhaul | 0 R code | 0 tests | Full rewrite (~200 lines) | Medium |
| Documentation examples | 20-40 per function | 0 | Roxygen updates | Low-Medium |

**Total estimated new R code:** 340-540 lines
**Total estimated new test code:** 300-460 lines

---

## Sources

### HIGH Confidence (codebase inspection, official documentation, verified computation)

- rbmiUtils codebase: Direct inspection of all R/ source files, tests/, DESCRIPTION (2026-02-10)
- rbmi pool object structure: Inspected via `str(pool_obj)` -- contains pars (est, ci, se, pvalue per param), conf.level, alternative, N, method
- rbmi analysis object structure: Inspected via `str(ana_obj)` -- contains results (per-imputation est, se, df), delta, fun, fun_name, method
- rbmi `rubin_rules()` source: Verified var_w, var_b, var_t computation; confirmed intermediate values are discarded
- rbmi `rubin_df()` source: Verified lambda computation; confirmed Barnard-Rubin df formula
- rbmi `as_draws()` source: Verified draws object structure (data, method, samples, fit, n_failures, formula)
- rbmi diagnostic functions: `check_ESS()` and `check_hmc_diagn()` exist as internal (unexported) functions
- MI diagnostics computation: Verified against ADMI data, producing correct FMI/lambda/RIV values
- [mice::pool() documentation](https://amices.org/mice/reference/pool.html) -- columns: estimate, ubar, b, t, dfcom, df, riv, lambda, fmi
- [Rubin's rules formulas](https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html) -- lambda, RIV, FMI, RE formulas verified
- [gt tab_options() reference](https://gt.rstudio.com/reference/tab_options.html) -- confirmed API for font.size, row padding, font weight
- cards package v0.7.1: Confirmed check_ard_structure() validation, as_card() conversion

### MEDIUM Confidence (web search verified with official sources)

- [R Consortium: cards package and ARD standard](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/) -- ARD as emerging CDISC standard
- [Introduction to clinical tables with gt](https://www.r-bloggers.com/2024/02/introduction-to-clinical-tables-with-the-gt-package/) -- gt for regulatory submission tables
- [Forester package](https://github.com/rdboyes/forester) -- publication-ready forest plot reference implementation (font_family, justify, ggplot_width parameters)
- [forestploter package](https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html) -- forest plot as aligned table+plot composition
- [PharmaSUG 2023 QT-263: R Tables via GT for Regulatory Submissions](https://pharmasug.org/proceedings/2023/QT/PharmaSUG-2023-QT-263.pdf)

### LOW Confidence (needs validation during implementation)

- Whether rstan convergence diagnostics (ESS, divergences, BFMI) are reliably accessible from draws_obj$fit across all rbmi versions
- Whether describe_imputation() can reliably detect which values were imputed without the original pre-imputation dataset (may need original_data parameter)
- Optimal font sizes and padding values for regulatory table publication -- may need iterative visual testing

---

*Feature landscape research for v3 ARD Enrichment & Polish: 2026-02-10*
