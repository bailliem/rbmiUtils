# Architecture Patterns: Reporting Layer for rbmiUtils

**Domain:** Clinical trial reporting in R (ARD-based reporting layer for an MI analysis package)
**Researched:** 2026-02-07
**Overall confidence:** HIGH (verified against official package documentation and source code)

---

## Executive Summary

This document defines how a new reporting layer should be architectured on top of the existing rbmiUtils analysis pipeline. The core design problem is bridging rbmi's pooled results (a tibble with columns `parameter`, `est`, `se`, `lci`, `uci`, `pval`) into the pharmaverse's ARD ecosystem (cards/cardx/gtsummary) for publication-ready clinical trial tables and figures. The architecture must also address improved S3 print/summary methods for pipeline objects and a cleaner wrapper around `rbmi::analyse()`.

The recommended architecture adds three new layers on top of the existing six, with clear data flow boundaries: analysis results flow downward through the pipeline, and each layer transforms data without reaching back up.

---

## 1. The cards ARD Format: Structure and Columns

**Confidence: HIGH** (verified via official cards documentation and CDISC ARD onboarding materials)

An ARD (Analysis Results Dataset) is a data frame of class `"card"` with standardized columns. The canonical columns are:

| Column | Type | Purpose | Example |
|--------|------|---------|---------|
| `group1` | character | Stratification variable name | `"TRT"` |
| `group1_level` | list | Level of the grouping variable | `"Drug A"` |
| `variable` | character | Variable being analyzed | `"CHG"` |
| `variable_level` | list | Category level (categorical vars) | `"Week 24"` |
| `stat_name` | character | Statistic identifier (machine key) | `"estimate"`, `"conf.low"`, `"p.value"` |
| `stat_label` | character | Human-readable label | `"Estimate"`, `"95% CI Lower"`, `"P-value"` |
| `stat` | list | The actual value (stored as list column) | `list(-2.45)` |
| `context` | character | Analysis type/function origin | `"ancova_trt_effect"` |
| `fmt_fn` | list | Formatting function for display | `list(function(x) sprintf("%.2f", x))` |
| `warning` | list | Warnings captured during computation | `list(NULL)` |
| `error` | list | Errors captured during computation | `list(NULL)` |

**Key structural properties:**
- The `stat` column is a **list column**, not a simple numeric vector. Each cell wraps one value.
- `group1`/`group1_level` can extend to `group2`/`group2_level` etc. for multiple stratification levels.
- The `fmt_fn` column carries formatting instructions, allowing downstream table builders to format values without knowing the statistic type.
- The `"card"` class enables S3 dispatch in gtsummary and cards utilities.

**Sources:**
- [cards package overview](https://www.danieldsjoberg.com/ARD-onboarding/03-cards-examples.html)
- [cards CRAN reference](https://cran.r-project.org/web/packages/cards/refman/cards.html)
- [CDISC COSA ARD + gtsummary 2025 slides](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/slides/)

---

## 2. Converting rbmiUtils Results to ARD

**Confidence: HIGH** (verified via as_card documentation and tidy_pool_obj source code)

### The bridge: tidy tibble to ARD

The existing `tidy_pool_obj()` produces a tibble with columns:

```
parameter | description | visit | parameter_type | lsm_type | est | se | lci | uci | pval
```

The conversion to ARD requires **pivoting** each statistic (est, se, lci, uci, pval) from wide columns into long `stat_name`/`stat` rows, then adding ARD metadata columns.

### Recommended approach: `as_card()` + manual pivot

The cards package (v0.7.1) provides `as_card(x)` which converts a data frame with at minimum `stat_name`, `stat_label`, and `stat` columns to an ARD of class `"card"`. The conversion function should:

1. Take the tidy tibble from `tidy_pool_obj()`
2. Pivot the statistical columns to long format
3. Map rbmiUtils column names to ARD column names
4. Add metadata columns (context, fmt_fn, warning, error)
5. Apply `cards::as_card()` to get the class

### Concrete column mapping

| rbmiUtils tidy tibble | ARD `stat_name` | ARD `stat_label` |
|----------------------|----------------|------------------|
| `est` | `"estimate"` | `"Estimate"` |
| `se` | `"std.error"` | `"Std. Error"` |
| `lci` | `"conf.low"` | `"CI Lower Bound"` |
| `uci` | `"conf.high"` | `"CI Upper Bound"` |
| `pval` | `"p.value"` | `"P-value"` |

### Group columns

| rbmiUtils tidy tibble | ARD column |
|-----------------------|-----------|
| `parameter_type` = `"trt"` | `variable` = `"trt_effect"` |
| `parameter_type` = `"lsm"` | `variable` = `"lsm"` |
| Treatment group (from context) | `group1` = treatment variable name, `group1_level` = arm |
| `visit` | `variable_level` = visit name |

### Formatting functions

Attach formatting functions to each stat type:

| `stat_name` | Default `fmt_fn` |
|------------|------------------|
| `"estimate"` | `function(x) sprintf("%.2f", x)` |
| `"std.error"` | `function(x) sprintf("%.2f", x)` |
| `"conf.low"` | `function(x) sprintf("%.2f", x)` |
| `"conf.high"` | `function(x) sprintf("%.2f", x)` |
| `"p.value"` | `format_pvalue` (reuse existing) |

### Alternative: `tidy_as_ard()` for model-based results

For future extensions where raw model objects are available (not just pooled summaries), `cards::tidy_as_ard()` converts broom-style tidy data frames into ARD. This is more suitable for single-model results (e.g., g-computation models) than for pooled MI results.

**Sources:**
- [as_card documentation](https://insightsengineering.github.io/cards/v0.3.0/reference/as_card.html)
- [tidy_as_ard documentation](https://search.r-project.org/CRAN/refmans/cards/html/tidy_as_ard.html)

---

## 3. How gtsummary Consumes ARDs

**Confidence: HIGH** (verified via gtsummary ARD-first tables vignette)

### ARD-first table functions

gtsummary provides a family of `tbl_ard_*()` functions that accept pre-computed ARD objects:

| Function | Purpose | Suitable for rbmiUtils? |
|----------|---------|------------------------|
| `tbl_ard_summary()` | Descriptive statistics tables | Possible but designed for summary stats |
| `tbl_ard_wide_summary()` | Wide-format statistics in separate columns | **Best fit for efficacy tables** |
| `tbl_ard_continuous()` | Continuous variable summaries | Useful for baseline characteristics |
| `tbl_ard_hierarchical()` | Nested/hierarchical data | Not needed for efficacy |

### The recommended gtsummary entry point

For MI efficacy tables (treatment effects + LSMs by visit), `tbl_ard_wide_summary()` is the best fit because:
- It places each statistic in its own column (Estimate, 95% CI, P-value)
- This matches the standard regulatory table layout
- It accepts a pre-computed ARD without recalculating anything

### API pattern

```r
# Conceptual usage
ard_obj <- pool_to_ard(pool_obj)  # rbmiUtils function

gtsummary::tbl_ard_wide_summary(
  cards = ard_obj,
  statistic = list(
    all_continuous() ~ c("{estimate}", "{conf.low}, {conf.high}", "{p.value}")
  )
)
```

### Supplementary ARD objects

gtsummary works best when the ARD also contains:
- **Attributes ARD** (via `cards::ard_attributes()`) -- variable labels
- **Missing data ARD** (via `cards::ard_missing()`) -- not needed for pooled results
- **Total N ARD** (via `cards::ard_total_n()`) -- sample sizes

These can be combined with `cards::bind_ard()` before passing to gtsummary.

### Important constraint

gtsummary's `tbl_ard_*()` functions do **not** calculate new statistics. Every number that appears in the table must already exist in the ARD. This means the rbmiUtils ARD conversion must be complete before gtsummary sees it.

**Sources:**
- [ARD-first Tables vignette](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html)
- [tbl_ard_summary reference](https://www.danieldsjoberg.com/gtsummary//reference/tbl_ard_summary.html)
- [tbl_ard_wide_summary reference](https://www.danieldsjoberg.com/gtsummary//reference/tbl_ard_wide_summary.html)

---

## 4. Abstraction Boundary: Analysis Results vs. Reporting

### The boundary line

The clean boundary sits at the **tidy tibble** produced by `tidy_pool_obj()`. Everything above this line is analysis; everything below is reporting.

```
ANALYSIS DOMAIN                     REPORTING DOMAIN
==================                  ==================
rbmi::draws()                       pool_to_ard()
rbmi::impute()                        |
get_imputed_data()                  tbl_efficacy()
analyse_mi_data()                   tbl_responder()
rbmi::pool()                        plot_forest()
tidy_pool_obj()  ---- boundary ---  plot_responder_bar()
                                    print/summary methods
```

### Why this boundary

1. **The tidy tibble is a stable contract.** It has well-defined columns (est, se, lci, uci, pval, visit, parameter_type, lsm_type) that do not depend on rbmi internals or cards version.

2. **Analysis functions should not know about presentation.** `analyse_mi_data()` and `tidy_pool_obj()` should never import cards or gtsummary.

3. **Reporting functions should not know about imputation mechanics.** The ARD converter and table builders only need the tidy tibble, not the imputation method or number of draws.

4. **Testing is cleanest at this boundary.** The tidy tibble can be mocked as a simple data frame for all reporting layer tests, without needing rbmi at all.

### Dependency classification

| Layer | Dependency on cards/gtsummary/gt | Classification |
|-------|----------------------------------|----------------|
| Data Preparation | None | No change |
| Analysis Execution | None | No change |
| Analysis Utilities | None | No change |
| Results Processing (tidy/format) | None | No change |
| **ARD Conversion (new)** | cards | **Suggests** (soft dependency) |
| **Table Generation (new)** | gtsummary, gt | **Suggests** (soft dependency) |
| **Visualization (new)** | ggplot2 | **Suggests** (already there) |
| Storage Management | None | No change |

**Recommendation: All reporting packages go in `Suggests`, not `Imports`.** This keeps rbmiUtils installable without the full reporting stack. Functions that need cards/gtsummary should check availability at runtime with `rlang::check_installed()`.

---

## 5. S3 Print/Summary Methods: Structure and Strategy

**Confidence: HIGH** (verified from existing source code in analyse_mi_data.R)

### Current state

rbmiUtils already has `print.analysis()` and `summary.analysis()` methods. These are well-structured but minimal. The package does **not** provide print/summary methods for rbmi's `pool`, `draws`, or `imputation` classes (those are owned by rbmi).

### Strategy: Provide methods for rbmiUtils-created objects only

| Object | Class | Current methods | Proposed changes |
|--------|-------|----------------|------------------|
| `analyse_mi_data()` output | `"analysis"` | `print.analysis`, `summary.analysis` | Improve: show parameter names, first result preview |
| `tidy_pool_obj()` output | tibble (no custom class) | Default tibble print | Add: Custom class `"tidy_pool"` with print method showing visit/parameter summary |
| `pool_to_ard()` output | `"card"` (from cards) | cards print method | No custom method needed -- cards handles this |

### Print method principles

1. **Header with class and provenance:** Show what the object is and how it was created.
2. **Key dimensions:** Number of imputations, parameters, visits.
3. **Next steps:** Actionable guidance ("Use `rbmi::pool()` to obtain pooled estimates").
4. **No heavy computation:** Print methods should be instant, never re-run analysis.

### Do NOT create methods for rbmi classes

Providing `print.pool()` or `summary.draws()` from rbmiUtils would conflict with rbmi's own methods (S3 dispatch ambiguity). Instead:

- If improved display of pool objects is needed, provide a **separate function** like `describe_pool(pool_obj)` rather than an S3 method.
- For draws/impute objects, use wrapper functions: `describe_draws(draws_obj)`, `describe_imputation(impute_obj)`.

---

## 6. Wrapping rbmi::analyse() vs. Reimplementing

**Confidence: HIGH** (verified from both rbmi quickstart docs and analyse_mi_data.R source)

### What rbmi::analyse() does

```r
rbmi::analyse(imputations, fun = ancova, delta = NULL, ...)
```

- Takes an `imputations` object (output of `rbmi::impute()`)
- Extracts each imputed dataset internally
- Applies optional delta adjustments
- Runs the analysis function on each dataset
- Returns an `analysis` object (class `c("analysis", "list")`)

### What analyse_mi_data() currently reimplements

The current `analyse_mi_data()` in rbmiUtils:
1. Takes a **flat data frame with IMPID column** (not an imputations object)
2. Validates inputs (data, vars, method, fun, delta)
3. Splits data by IMPID using `group_split()`
4. Applies the analysis function to each split
5. Constructs an `analysis` object using `as_analysis2()` (a modified copy of rbmi's `as_analysis()`)

### Key difference: Input type

| Aspect | `rbmi::analyse()` | `analyse_mi_data()` |
|--------|-------------------|---------------------|
| Input | `imputations` object (from `rbmi::impute()`) | Flat data.frame with IMPID column |
| Data extraction | Internal (uses `extract_imputed_dfs()`) | Already done (user provides flat data) |
| Delta handling | Via imputations object metadata | Via separate delta parameter + left_join |
| Class construction | Internal `as_analysis()` | Copied `as_analysis2()` |

### Recommendation: Do NOT wrap rbmi::analyse() directly

Wrapping `rbmi::analyse()` would mean requiring the user to pass an `imputations` object, which defeats the purpose of `analyse_mi_data()`. The whole point of the rbmiUtils function is to work with flat data frames that have been extracted (potentially reduced/expanded, or created from non-rbmi imputations via `create_impid()`).

**Instead, refactor `analyse_mi_data()` to:**

1. **Remove the copied `as_analysis2()`** -- use `rbmi::as_analysis()` directly if it is exported, or keep a thin maintained copy with clear version pinning.
2. **Add delta uniqueness validation** (currently missing -- see CONCERNS.md).
3. **Improve IMPID validation** -- check for NA, non-unique within groups, etc.
4. **Keep the flat-data-frame API** -- this is the correct abstraction for the use case.

The refactoring is a **cleaning exercise**, not an API change. The function signature stays the same.

---

## 7. Recommended Architecture: 9-Layer Design

### Component Diagram

```
Layer 1: Data Preparation          [EXISTING - no changes]
  validate_data(), prepare_data_ice(), summarise_missingness()
      |
      v
Layer 2: Analysis Execution        [EXISTING - refactor internals only]
  analyse_mi_data() -> analysis S3 class
      |
      v
Layer 3: Analysis Utilities        [EXISTING - harden]
  gcomp_responder(), gcomp_binary(), gcomp_responder_multi()
      |
      v
Layer 4: Results Processing        [EXISTING - minor improvements]
  tidy_pool_obj(), format_results(), combine_results()
      |
      v  (this is the abstraction boundary)
      |
Layer 5: ARD Conversion            [NEW]
  pool_to_ard(), tidy_to_ard()
  Depends: cards (Suggests)
      |
      v
Layer 6: Table Generation          [NEW]
  tbl_efficacy(), tbl_responder()
  Depends: gtsummary, gt (Suggests)
      |
Layer 7: Visualization             [NEW]
  plot_forest(), plot_responder_bar()
  Depends: ggplot2 (Suggests, already present)

Layer 8: Storage Management        [EXISTING - harden]
  reduce_imputed_data(), expand_imputed_data()

Layer 9: Utility Functions         [EXISTING - no changes]
  get_imputed_data(), create_impid()
```

### Component Boundaries

| Component | Responsibility | Communicates With | Never Touches |
|-----------|---------------|-------------------|---------------|
| Data Preparation | Validate and prepare data for rbmi | User data, rbmi::set_vars() | Analysis results |
| Analysis Execution | Run analysis across imputations | Prepared data, user functions | Reporting, ARDs |
| Analysis Utilities | Specialized analysis functions | beeca, stats | Reporting, ARDs |
| Results Processing | Transform pool objects to tidy tibbles | rbmi::pool() output | cards, gtsummary |
| **ARD Conversion** | Transform tidy tibbles to ARD format | Results Processing output, cards | gtsummary, ggplot2 |
| **Table Generation** | Produce publication-ready tables from ARDs | ARD Conversion output, gtsummary | rbmi, analysis |
| **Visualization** | Produce publication-ready figures | Results Processing output, ggplot2 | cards, gtsummary |
| Storage Management | Compress/expand imputed datasets | Imputed data frames | Everything else |
| Utility Functions | Extract data, create IDs | rbmi imputation objects | Everything else |

### Data Flow

```
raw data
  |
  v
validate_data() -> prepare_data_ice() -> summarise_missingness()
  |
  v
[user calls rbmi::draws() -> rbmi::impute()]
  |
  v
get_imputed_data()  ------>  reduce_imputed_data()  (storage path)
  |                                    |
  v                                    v
analyse_mi_data()           expand_imputed_data()  (retrieval path)
  |
  v
rbmi::pool()
  |
  v
tidy_pool_obj()
  |
  +--- format_results()    (existing formatting path)
  |    combine_results()
  |
  +--- pool_to_ard()       (new ARD path, requires cards)
  |      |
  |      +--- tbl_efficacy()      (requires gtsummary + gt)
  |      +--- tbl_responder()     (requires gtsummary + gt)
  |
  +--- plot_forest()        (new visualization path, requires ggplot2)
  +--- plot_responder_bar() (new visualization path, requires ggplot2)
```

---

## 8. Suggested Build Order (Dependencies Between Components)

### Phase ordering rationale

Components must be built in dependency order. Each phase should be independently testable before the next begins.

### Phase 1: Foundation Hardening (no new dependencies)

**Build first because:** Everything else depends on these working correctly.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `analyse_mi_data()` | Refactor internals (remove copied `as_analysis2`, add delta validation) | rbmi only |
| `gcomp_*` functions | Add edge-case handling, input validation | beeca only |
| `reduce/expand_imputed_data()` | Fix type coercion, attribute preservation | dplyr only |
| `validate_data()`, `prepare_data_ice()` | Fill validation gaps | rbmi only |
| `tidy_pool_obj()` | Improve parameter parsing robustness | dplyr, tidyr only |
| `format_results()`, `format_results_table()` | Edge case handling | existing deps only |

**Test strategy:** Existing test infrastructure, no new test dependencies.

### Phase 2: Improved Print/Summary Methods (no new dependencies)

**Build second because:** Improves developer experience for all subsequent phases.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `print.analysis()` | Enhance with parameter preview | None new |
| `summary.analysis()` | Add result structure inspection | None new |
| Tidy pool class | Add `"tidy_pool"` class + print method to `tidy_pool_obj()` output | None new |
| Descriptive helpers | `describe_pool()`, `describe_draws()`, `describe_imputation()` | rbmi only |

**Test strategy:** Snapshot tests for print output; existing test infrastructure.

### Phase 3: ARD Conversion Layer (adds cards dependency)

**Build third because:** Tables and some figures depend on ARD objects.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `pool_to_ard()` | Convert pool object to ARD (pool -> tidy -> pivot -> as_card) | cards (Suggests) |
| `tidy_to_ard()` | Convert tidy tibble to ARD (for users who already have tidy results) | cards (Suggests) |
| `bind_ard` helpers | Combine efficacy ARD with attributes/N ARDs | cards (Suggests) |

**Test strategy:** Skip tests when cards not installed (`testthat::skip_if_not_installed("cards")`).

### Phase 4: Table Generation (adds gtsummary + gt dependencies)

**Build fourth because:** Depends on ARD conversion from Phase 3.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `tbl_efficacy()` | Efficacy summary table (treatment effects by visit) | gtsummary, gt (Suggests) |
| `tbl_responder()` | Responder analysis table (proportions by arm and visit) | gtsummary, gt (Suggests) |

**Test strategy:** Skip tests when gtsummary/gt not installed. Use snapshot tests for table output.

### Phase 5: Visualization (uses existing ggplot2 Suggests dependency)

**Build fifth because:** Can be built in parallel with Phase 4 since figures use the tidy tibble directly, not ARDs.

| Component | Action | Dependencies |
|-----------|--------|--------------|
| `plot_forest()` | Forest plot of treatment effects | ggplot2 (already Suggests) |
| `plot_responder_bar()` | Bar chart of responder proportions | ggplot2 (already Suggests) |

**Test strategy:** `vdiffr` for visual regression testing, or `testthat::skip_if_not_installed("ggplot2")`.

### Phase dependency graph

```
Phase 1 (Hardening)
    |
    v
Phase 2 (Print/Summary)
    |
    v
Phase 3 (ARD Conversion)  <-- adds cards dependency
    |           \
    v            \
Phase 4          Phase 5
(Tables)         (Figures)
adds gtsummary   uses ggplot2
```

Phases 4 and 5 can run in parallel after Phase 3. Phase 5 could even run in parallel with Phase 3 since figures use the tidy tibble directly, not ARDs.

---

## Patterns to Follow

### Pattern 1: Soft Dependency Check

**What:** Check for optional packages at function entry, not at package load time.

**When:** Any function in Layers 5-7 that depends on packages in `Suggests`.

```r
pool_to_ard <- function(pool_obj, ...) {
  rlang::check_installed("cards", reason = "to convert results to ARD format")
  # ... function body
}
```

### Pattern 2: Tidy-first, ARD-second

**What:** Always convert pool objects to tidy tibbles first, then to ARD.

**When:** Any reporting function that takes a pool object.

```r
# Good: Two-step conversion
tidy_df <- tidy_pool_obj(pool_obj)
ard <- tidy_to_ard(tidy_df)

# Bad: Direct pool-to-ARD without intermediate tidy step
ard <- magic_pool_to_ard(pool_obj)  # skips the stable contract
```

### Pattern 3: Mocked Input Tests for Reporting Layer

**What:** Test reporting functions with hand-crafted tidy tibbles, not real rbmi pipeline output.

**When:** All tests for Layers 5-7.

```r
test_that("pool_to_ard creates valid ARD", {
  skip_if_not_installed("cards")

  mock_tidy <- tibble::tibble(
    parameter = "trt_Week 24",
    description = "Treatment Comparison at Week 24",
    visit = "Week 24",
    parameter_type = "trt",
    lsm_type = NA_character_,
    est = -2.45, se = 0.89, lci = -4.20, uci = -0.70, pval = 0.006
  )

  ard <- tidy_to_ard(mock_tidy)
  expect_s3_class(ard, "card")
})
```

---

## Anti-Patterns to Avoid

### Anti-Pattern 1: Tight Coupling to ARD Column Names

**What:** Hardcoding ARD column names throughout the codebase.

**Why bad:** If cards changes column naming conventions (they have between v0.3 and v0.7), every function breaks.

**Instead:** Define column name mappings in a single internal constant, reference that constant everywhere. Isolate all cards-specific code to the ARD Conversion layer.

### Anti-Pattern 2: Importing cards/gtsummary in Imports

**What:** Adding cards/gtsummary to `Imports:` in DESCRIPTION.

**Why bad:** Forces every rbmiUtils user to install the full reporting stack even if they only need analysis utilities. Increases installation failures on restricted environments (e.g., validated pharma computing environments).

**Instead:** Keep all reporting packages in `Suggests:`. Use `rlang::check_installed()` at function entry.

### Anti-Pattern 3: Bypassing the Tidy Tibble Contract

**What:** Having table/figure functions reach back into pool or analysis objects directly.

**Why bad:** Creates hidden dependencies on rbmi internals. Breaks if rbmi changes pool object structure. Makes testing harder.

**Instead:** All reporting functions accept tidy tibbles (or ARDs derived from tidy tibbles). Never pass raw pool objects to gtsummary functions.

### Anti-Pattern 4: Registering S3 Methods for Classes You Do Not Own

**What:** Creating `print.pool()` or `summary.draws()` methods in rbmiUtils for rbmi classes.

**Why bad:** S3 method dispatch conflict. If both rbmi and rbmiUtils export `print.pool`, the user gets whichever package was loaded last. Creates maintenance coupling with rbmi releases.

**Instead:** Use standalone helper functions: `describe_pool()`, `describe_draws()`, `describe_imputation()`.

---

## Open Questions for Phase-Specific Research

| Question | When to resolve | Impact |
|----------|----------------|--------|
| Which `tbl_ard_*()` function best fits efficacy tables? | Phase 4 start | Table function design |
| Does `tbl_ard_wide_summary()` handle mixed continuous/treatment-effect rows? | Phase 4 start | May need custom gt assembly instead |
| Should forest plots use pooled CIs or recomputed CIs from ARD? | Phase 5 start | Determines figure function input |
| Is `cards::as_card()` stable across cards v0.5-v0.7? | Phase 3 start | May need version guard |
| Does `rbmi::as_analysis()` export publicly or only internally? | Phase 1 start | Determines refactoring approach for `as_analysis2` |

---

## Sources

### HIGH confidence (official documentation, verified)

- [cards package CRAN reference](https://cran.r-project.org/web/packages/cards/refman/cards.html) -- function listing and ARD structure
- [cards as_card() documentation](https://insightsengineering.github.io/cards/v0.3.0/reference/as_card.html) -- ARD conversion from data frames
- [cards tidy_as_ard() documentation](https://search.r-project.org/CRAN/refmans/cards/html/tidy_as_ard.html) -- model-based ARD construction
- [gtsummary ARD-first tables vignette](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html) -- how gtsummary consumes ARDs
- [gtsummary tbl_ard_summary reference](https://www.danieldsjoberg.com/gtsummary//reference/tbl_ard_summary.html) -- function signature and parameters
- [gtsummary tbl_ard_wide_summary reference](https://www.danieldsjoberg.com/gtsummary//reference/tbl_ard_wide_summary.html) -- wide-format table function
- [rbmi quickstart vignette](https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html) -- analyse() and pool() workflow
- [rbmi analyse() reference](https://search.r-project.org/CRAN/refmans/rbmi/html/analyse.html) -- function signature
- [rbmi pool() reference](https://search.r-project.org/CRAN/refmans/rbmi/html/pool.html) -- function signature and return structure
- [rbmi pool.R source](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R) -- pool object class structure

### MEDIUM confidence (workshop slides, verified with official docs)

- [CDISC COSA ARD + gtsummary 2025 slides](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/slides/) -- end-to-end ARD workflow
- [PHUSE ARD workshop 2025 - cardx slides](https://www.danieldsjoberg.com/ARD-PHUSE-workshop-2025/slides/02-intro-cardx/) -- custom ARD patterns with cardx
- [ARD onboarding tutorial](https://www.danieldsjoberg.com/ARD-onboarding/03-cards-examples.html) -- cards examples and column structure
- [R Consortium - ARDs and cards](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/) -- ecosystem overview

### Internal (existing codebase, verified by reading source)

- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/analyse_mi_data.R` -- analyse_mi_data() and as_analysis2()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/tidiers.R` -- tidy_pool_obj() column structure
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/formatting.R` -- format_results(), format_pvalue(), format_estimate()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/result_helpers.R` -- combine_results(), extract_trt_effects(), extract_lsm()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/R/analysis_utils.R` -- gcomp_responder(), gcomp_responder_multi()
- `/Users/bailliem/R-projects/rbmiUtils-gsd/.planning/codebase/ARCHITECTURE.md` -- existing architecture analysis
- `/Users/bailliem/R-projects/rbmiUtils-gsd/.planning/codebase/CONCERNS.md` -- known issues and tech debt
- `/Users/bailliem/R-projects/rbmiUtils-gsd/.planning/PROJECT.md` -- milestone requirements

---

*Architecture research: 2026-02-07*
