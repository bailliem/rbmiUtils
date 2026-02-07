# Domain Pitfalls

**Domain:** Clinical trial reporting layer for an R package (rbmiUtils) extending rbmi
**Researched:** 2026-02-07
**Overall confidence:** MEDIUM-HIGH (grounded in codebase analysis + verified ecosystem documentation)

---

## Critical Pitfalls

Mistakes that cause rewrites or major issues.

### Pitfall 1: ARD Column Contract Violations When Converting Custom Results

**What goes wrong:** The cards/ARD ecosystem requires a precise column contract: `variable`, `stat_name`, `stat_label`, `stat`, `group1`, `group1_level`, `context`, `fmt_fn`, `warning`, `error`. When converting rbmiUtils' tidy pool results (which use `est`, `se`, `lci`, `uci`, `pval`) into ARD format, developers commonly miss required metadata columns or use wrong column types. The `stat` column must be a list-column (each cell is a list element), not a simple numeric vector. The `fmt_fn` column must contain formatting functions (or NULL), not formatted strings. Missing the `context` column causes `tbl_ard_summary()` to fail silently or produce empty tables.

**Why it happens:** The tidy tibble from `tidy_pool_obj()` uses a flat, human-readable structure (one row per parameter, numeric columns for estimates). The ARD format is a long, machine-readable structure (one row per statistic, list-column for values). The mental model mismatch leads to surface-level column renaming instead of structural transformation.

**Consequences:** `gtsummary::tbl_ard_summary()` either errors with cryptic messages about missing columns, or silently drops statistics that lack proper `stat_name`/`stat_label` pairing. Tables render with missing cells and no warning. Debugging requires understanding the cards internals.

**Warning signs:**
- Unit tests that only check column names exist, not column types
- `tbl_ard_summary()` producing tables with fewer rows than expected
- Errors mentioning "Level column missing" or "`stat` must be a list"
- The ARD object prints without the `{cards} data frame` header

**Prevention:**
1. Use `cards::as_cards_fn()` to wrap the conversion function, which guarantees consistent output structure even on error (specifying `stat_names` upfront).
2. Validate the output with `cards::check_ard_structure()` (or equivalent) before passing to gtsummary.
3. Write explicit tests that the returned object inherits the `"card"` class and that `stat` is a list-column.
4. Map rbmiUtils statistics to ARD stat_names using a documented lookup: `est` -> `"estimate"`, `se` -> `"std.error"`, `lci` -> `"conf.low"`, `uci` -> `"conf.high"`, `pval` -> `"p.value"`.
5. Include `context = "rbmiUtils_pool"` to enable gtsummary to dispatch formatting correctly.

**Phase:** Reporting layer (ARD conversion). Address first before building gtsummary templates, since templates consume ARDs.

**Confidence:** MEDIUM -- cards API verified via official documentation; exact column requirements confirmed from [insightsengineering/cards](https://insightsengineering.github.io/cards/main/) and [gtsummary ARD-first article](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html). Specific error messages are inferred from package structure, not directly tested.

---

### Pitfall 2: Parameter Name Parsing Breaks on Underscored Visit or Treatment Names

**What goes wrong:** `tidy_pool_obj()` uses `tidyr::separate(parameter, into = c("parameter_type", "lsm_type", "visit"), sep = "_", fill = "right", extra = "merge")` to parse parameter names like `"trt_Week24"` or `"lsm_ref_Week 24"`. When visit names contain underscores (e.g., `"Week_24"`, `"Follow_Up_Visit"`) or treatment group names contain underscores passed through beeca (e.g., `"Drug_A_vs_Placebo"`), the separator becomes ambiguous. The `extra = "merge"` flag merges excess parts into the rightmost column, but the first split still consumes the wrong token.

**Why it happens:** The parameter naming convention uses `_` as both a structural delimiter (separating type/arm/visit) and a character that legitimately appears in visit and treatment names. This is a classic delimiter-collision problem.

**Consequences:** Silent data corruption: `parameter_type` gets assigned partial strings, `visit` contains truncated names, `lsm_type` is populated with fragments of the visit name. Downstream filtering in `extract_trt_effects()` and `extract_lsm()` returns empty results. `format_results()` produces tables with wrong visit labels. Because the errors are silent (no exception thrown), they may reach publication.

**Warning signs:**
- `tidy_pool_obj()` output has `NA` in the `visit` column for some rows
- `parameter_type` values other than `"trt"` or `"lsm"` appear
- `extract_trt_effects()` returns zero rows when data clearly has treatment comparisons
- Visit names in the output differ from visit names in the input data

**Prevention:**
1. Replace `tidyr::separate()` with structured regex parsing: `"^(trt|lsm)_(ref_|alt_)?(.+)$"` using `tidyr::extract()` or `stringr::str_match()`. This anchors the known prefixes and treats everything after as the visit.
2. Alternatively, store parameter metadata (type, arm, visit) as attributes during `analyse_mi_data()` so that parsing is unnecessary.
3. Add a validation step: after parsing, assert that `parameter_type %in% c("trt", "lsm")` and warn if any parsed values are unexpected.
4. Migrate away from the deprecated `tidyr::separate()` to `tidyr::separate_wider_regex()` which enables named capture groups.

**Phase:** Hardening phase (tidier/formatter functions). Must be fixed before ARD conversion relies on correctly parsed parameters.

**Confidence:** HIGH -- this fragility is directly observable in `R/tidiers.R` lines 98-106 and confirmed in CONCERNS.md. The `tidyr::separate()` supersession is documented at [tidyr reference](https://tidyr.tidyverse.org/reference/separate.html).

---

### Pitfall 3: rbmi Class Hierarchy Indexing with `class(method)[[2]]`

**What goes wrong:** `analyse_mi_data()` and `as_analysis2()` use `class(method)[[2]]` to determine the method type (bayes, approxbayes, condmean, bmlmi) and `switch()` to assign the pooling subclass. This assumes the class vector has a specific structure where the second element is the method type. If rbmi changes its class hierarchy -- adding, removing, or reordering class elements -- this code silently produces wrong pooling assignments or errors with "subscript out of bounds".

**Why it happens:** S3 class vectors in R are ordered character vectors, and their positional access is inherently fragile. rbmi moved from 1.4 to 1.5.x in 2025, adding covariance structure support and new method parameters. While the class vector structure has been stable so far, there is no API contract guaranteeing it will remain at position `[[2]]`.

**Consequences:** Wrong pooling method assignment means `rbmi::pool()` applies incorrect variance pooling rules. For clinical trial primary analyses, this directly corrupts reported treatment effects, confidence intervals, and p-values. The error is silent -- pool() may succeed but produce wrong numbers.

**Warning signs:**
- After upgrading rbmi, pooled estimates change unexpectedly
- `print.analysis()` shows unexpected `Method type:` or `Pooling method:` values
- `rbmi::pool()` errors with "method not found for class..."
- `class(method)` has a different length than expected

**Prevention:**
1. Replace `class(method)[[2]]` with `inherits()` checks: `if (inherits(method, "bayes")) ...` which is stable across class vector reordering.
2. Alternatively, wrap `rbmi::analyse()` directly instead of reimplementing the method-to-pooling mapping. This is already an active requirement in PROJECT.md.
3. Add an explicit assertion: `stopifnot(class(method)[[2]] %in% c("bayes", "approxbayes", "condmean", "bmlmi"))` to fail loud rather than silently.
4. Pin rbmi version more tightly: `rbmi (>= 1.4, < 2.0)` in DESCRIPTION to bound the compatibility window.
5. Add integration tests that verify the analysis-to-pool pipeline produces expected results for each method type.

**Phase:** Refactor of `analyse_mi_data()` to wrap `rbmi::analyse()`. This is the highest-priority structural fix because it eliminates the fragility entirely rather than patching it.

**Confidence:** HIGH -- the indexing pattern is directly visible in `R/analyse_mi_data.R` lines 148-155 and 271-281. rbmi version history confirms 1.4 -> 1.5 releases in 2025 from [CRAN rbmi news](https://cran.r-project.org/web/packages/rbmi/news/news.html).

---

### Pitfall 4: beeca Output Format Coupling in G-computation Functions

**What goes wrong:** `gcomp_responder()` and `gcomp_binary()` directly access columns named `STAT`, `STATVAL`, and `TRTVAL` from `beeca::get_marginal_effect()$marginal_results`. They also filter on specific `STAT` values like `"diff"`, `"diff_se"`, `"risk"`, `"risk_se"`. If beeca renames these columns or changes the stat naming convention, these functions break.

**Why it happens:** beeca has no version constraint in DESCRIPTION (`Imports: beeca` with no version). The output format is treated as a stable API, but beeca is a relatively young package (pharmaverse ecosystem, active development). Column-name coupling to upstream return values is a common integration fragility.

**Consequences:** `gcomp_responder()` returns empty data frames (when filter matches zero rows) or errors on missing column names. Since this function is used inside `analyse_mi_data()` across hundreds of imputations, a single beeca breaking change crashes the entire analysis pipeline with an unhelpful error deep in a lapply loop.

**Warning signs:**
- `beeca::get_marginal_effect()` output structure changes after package update
- `gcomp_responder()` returns lists with `NULL` or `numeric(0)` for `est` or `se`
- Errors mentioning "Column STAT not found" or "Column STATVAL not found"
- Binary outcome analysis produces different results after upgrading beeca

**Prevention:**
1. Add `beeca (>= X.Y)` version constraint to DESCRIPTION Imports.
2. Validate beeca output structure immediately after the `get_marginal_effect()` call: check that `STAT`, `STATVAL`, `TRTVAL` columns exist.
3. Write integration tests that call `beeca::get_marginal_effect()` on a fixed dataset and assert the output schema.
4. Abstract the column-name mapping into a single location (e.g., a named vector `BEECA_COLS <- c(stat = "STAT", value = "STATVAL", treatment = "TRTVAL")`) so that a beeca update requires changing only one place.
5. Consider wrapping the beeca call in a tryCatch that produces a clear error message mentioning the beeca version.

**Phase:** Hardening phase (gcomp functions). Should be addressed before ARD conversion attempts to consume gcomp results.

**Confidence:** HIGH -- directly observable in `R/analysis_utils.R` lines 83-100 and `R/utils.R` lines 163-171. beeca absence of version constraint confirmed in DESCRIPTION.

---

## Moderate Pitfalls

Mistakes that cause delays or technical debt.

### Pitfall 5: gtsummary Template Brittleness with Regulatory Table Structures

**What goes wrong:** Developers build gtsummary table templates assuming `tbl_ard_summary()` or `tbl_summary()` can handle the specific structure of clinical trial efficacy tables (treatment effect + LSM by visit, with specific header/footer formatting). Regulatory tables have non-standard layouts: multiple statistics per cell (estimate + CI on one line, p-value below), row grouping by visit with treatment effects and LSMs in the same table, footnotes with method descriptions. gtsummary's template system is designed for single-statistic-per-cell summaries.

**Why it happens:** gtsummary excels at demographic/baseline tables (Table 1) and simple comparison tables. Efficacy summary tables from clinical trials have more complex layouts. Developers start with `tbl_summary()` then discover they need heavy customization via `modify_*()` functions, custom column merging, or raw gt manipulation.

**Prevention:**
1. Prototype the exact target table layout (using a mock-up or existing SAS/R output) before writing code.
2. Use `tbl_ard_wide_summary()` for efficacy tables with multiple statistics rather than `tbl_ard_summary()`, since it supports wider layouts.
3. For complex regulatory layouts, consider building the gt table directly from the ARD rather than going through gtsummary's abstractions.
4. Build the ARD -> gtsummary -> gt pipeline incrementally: first get the data in, then customize formatting.
5. List gtsummary and gt as Suggests (not Imports) to avoid forcing the dependency on users who only need the analysis functions.

**Phase:** Reporting layer. Address after ARD conversion is working. Expect iteration.

**Confidence:** MEDIUM -- based on gtsummary documentation at [gtsummary tbl_ard functions](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html) and [CDISC COSA ARD talk](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/). The complexity of regulatory table layouts is domain knowledge from clinical trial reporting conventions.

---

### Pitfall 6: Forest Plot Scale, Ordering, and Reference Line Mistakes

**What goes wrong:** Clinical trial forest plots require specific conventions that ggplot2 does not enforce by default:

1. **Scale mismatch:** Using linear scale for odds/hazard ratios (should be log scale) or log scale for mean differences (should be linear). rbmiUtils will report mean differences from ANCOVA and risk differences from g-computation, which need linear scales, but if the package later supports odds ratios, the same function will produce misleading plots.
2. **Visit ordering reversed:** ggplot2 orders factor levels bottom-to-top on the y-axis. Clinical convention is chronological top-to-bottom. Developers forget `scale_y_discrete(limits = rev)` or fail to set factor levels explicitly.
3. **Reference line at wrong value:** Reference line should be at 0 for differences and 1 for ratios. Hardcoding either value makes the function wrong for the other scale.
4. **Label truncation:** Long visit names or subgroup labels overflow the plot area, especially when combining with CI annotations.
5. **Diamond vs. point size encoding:** Developers encode effect size in point size (common in meta-analysis forest plots for study weight), which is inappropriate for clinical trial visit-level effects where all visits have equal standing.

**Prevention:**
1. Parameterize `reference_line` with a sensible default (0 for difference scale) and document when to change it.
2. Accept a `scale` argument (`"identity"` vs `"log"`) that controls both the x-axis transform and the reference line default.
3. Always reverse y-axis factor levels for chronological top-to-bottom display.
4. Use `ggplot2::coord_cartesian(clip = "off")` with `ggplot2::geom_text(hjust = 0)` outside the plot area for labels, or use a two-panel layout (labels | plot) via patchwork or cowplot.
5. Use fixed point size for clinical trial forest plots; reserve size encoding for meta-analysis contexts.
6. Return the ggplot object (not a rendered plot) so users can customize further with standard ggplot2 operations.

**Phase:** Reporting layer (forest plot function). Can be developed in parallel with ARD conversion.

**Confidence:** MEDIUM -- clinical trial forest plot conventions are well-established domain knowledge. ggplot2 axis ordering behavior verified via [ggplot2 documentation](https://ggplot2.tidyverse.org/). Specific clinical trial conventions cross-referenced with [KHstats forest plot guide](https://www.khstats.com/blog/forest-plots/) and [forestplot vignette](https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html).

---

### Pitfall 7: print/summary S3 Method Dispatch Conflicts with rbmi's Own Methods

**What goes wrong:** rbmiUtils defines `print.analysis` and `summary.analysis` for its own `analysis` class. If rbmi also defines or later adds methods for the same class, R's S3 dispatch will use whichever package was loaded last (or whichever is found first on the search path). Additionally, if `rbmi::pool()` returns an object that inherits from multiple classes, adding `print.pool` in rbmiUtils could shadow rbmi's own `print.pool` method.

**Why it happens:** Both rbmiUtils and rbmi operate on the same S3 class names (`analysis`, `pool`). R's S3 dispatch is namespace-aware for registered methods but still susceptible to masking via `library()` load order. The `analysis` class in rbmiUtils is constructed by `as_analysis2()` to be compatible with rbmi's `pool()`, meaning the class vectors overlap.

**Consequences:** Users get different output depending on package load order. If rbmi adds a `print.analysis` method in a future version, there is a method conflict. Worse, if the method signatures differ (rbmi's method expects different internal structure), calling `print()` on an rbmiUtils analysis object could error or print misleading information.

**Warning signs:**
- `methods(print)` shows duplicate registrations for `print.analysis`
- R CMD check WARNING about "method masked" or "S3 method overwritten"
- `print()` output changes depending on whether rbmi or rbmiUtils is loaded first

**Prevention:**
1. For the `analysis` class: since rbmiUtils creates these objects via `as_analysis2()`, it legitimately owns the print/summary methods. Ensure they are registered in NAMESPACE with `S3method(print, analysis)` and `S3method(summary, analysis)`.
2. For the `pool` class: do NOT add print/summary methods directly. Instead, create a wrapper class (e.g., `"rbmiUtils_pool"`) that inherits from `"pool"` and add methods for the wrapper class.
3. Test for method conflicts: in test suite, check `methods(class = "analysis")` does not contain entries from rbmi.
4. Use the `...` argument in print/summary methods even if currently unused, to maintain forward compatibility.
5. Ensure `print()` returns `invisible(x)` and `summary()` returns a summary object (invisibly), following R conventions.

**Phase:** print/summary improvements. Address during that dedicated phase.

**Confidence:** HIGH -- the current `print.analysis` and `summary.analysis` implementations are visible in `R/analyse_mi_data.R` lines 330-442. S3 dispatch rules verified via [Advanced R S3 chapter](https://adv-r.hadley.nz/s3.html).

---

### Pitfall 8: Optional Dependencies (Suggests) Break CRAN Checks or Silently Degrade

**What goes wrong:** Adding cards, cardx, gtsummary, and gt as Suggests (rather than Imports) requires guarding every call to these packages with availability checks. Developers commonly miss guards in:
- Helper functions called by the main ARD conversion function
- Default argument values that reference Suggests packages
- Examples in roxygen documentation
- Vignette code that runs during `R CMD check`

CRAN policy requires that package functionality must not error if Suggests are not installed. Additionally, the cardx package has its own hidden Suggests dependencies (e.g., emmeans, broom.helpers) that `{renv}` will not record unless explicitly referenced.

**Why it happens:** During development, Suggests packages are always installed, so missing guards are invisible. The failure only manifests during `R CMD check --no-suggests` or when a user installs rbmiUtils without the optional packages.

**Consequences:** CRAN rejection if examples or tests call Suggests packages without guards. Users who install only rbmiUtils (without gtsummary) get cryptic "there is no package called 'gtsummary'" errors when they call any function, even unrelated ones, if the package-level imports are not properly conditional.

**Warning signs:**
- `R CMD check --no-suggests` produces errors
- `.Rd` examples fail when Suggests are not installed
- Tests fail on CI when Suggests installation is skipped
- `renv::snapshot()` misses cards/gtsummary/gt in lockfile

**Prevention:**
1. Use `rlang::check_installed("cards", reason = "to convert results to ARD format")` at the top of every function that uses a Suggests package. This gives clear error messages and offers installation in interactive sessions.
2. Wrap all examples using Suggests in `\donttest{}` or guard with `if (requireNamespace("gtsummary", quietly = TRUE))`.
3. Add `testthat::skip_if_not_installed("cards")` to every test that uses Suggests.
4. Run `R CMD check --no-suggests` in CI to catch unguarded calls early.
5. For cardx hidden dependencies: add explicit `library(emmeans)` or similar references in vignettes so renv captures them.
6. Consider making `cards` an Import (not Suggests) if the ARD conversion is core functionality, and keep `gtsummary`/`gt` as Suggests since they are only needed for table rendering.

**Phase:** Reporting layer and throughout. The dependency classification decision (Imports vs Suggests) should be made early in planning, as it affects every subsequent function.

**Confidence:** HIGH -- CRAN Suggests policy verified via [R Packages (2e) Dependencies chapter](https://r-pkgs.org/dependencies-in-practice.html). rlang guard pattern verified via [rlang is_installed reference](https://rlang.r-lib.org/reference/is_installed.html). cardx hidden dependency issue documented in [cardx README](https://github.com/insightsengineering/cardx).

---

### Pitfall 9: Formula Construction from Unsanitized User Input

**What goes wrong:** `as_simple_formula2()` in `R/analysis_utils.R` constructs formulas by pasting user-supplied covariate names: `paste0(outcome, " ~ 1 + ", paste0(covars, collapse = " + "))`. If covariate names contain special characters, spaces, or R formula operators (e.g., `"log(BASE)"`, `"I(AGE^2)"`, `"BASE * TRT"`), the resulting formula may be syntactically valid but semantically wrong, or may error with an unhelpful message from `stats::glm()`.

**Why it happens:** Clinical trial datasets sometimes use variable names with embedded transformations or interaction terms specified as single strings. The `extract_covariates2()` function splits on `:` and `*` but does not handle parenthesized expressions, `I()` wrappers, or spaces in variable names.

**Consequences:** Model fitting errors deep inside `stats::glm()` or `stats::lm()` with stack traces that point to the model formula rather than the user's input. For legitimate R formula operators embedded in covariate strings (like `"BASE:TRT"`), the intent may be an interaction, but the naive splitting in `extract_covariates2()` will try to validate `"BASE"` and `"TRT"` as separate columns.

**Prevention:**
1. Validate covariate names against `names(data)` before formula construction, and report unmatched names clearly.
2. For interaction terms: support them explicitly with documentation and testing, or reject them with a clear error message saying "use `*` syntax in the covariates argument".
3. For transformations: either support `I()` and function calls in covariates, or reject them early with "transform variables before passing to analysis functions".
4. Sanitize covariate names: reject names containing characters other than `[A-Za-z0-9_.]` unless they match known formula operators.
5. Environment scoping: `as_simple_formula2()` sets `environment(frm) <- globalenv()` which may fail if variables are defined in a non-global environment. Consider using `environment(frm) <- parent.frame()` or the calling environment.

**Phase:** Hardening phase (gcomp functions).

**Confidence:** HIGH -- directly visible in `R/analysis_utils.R` lines 205-214.

---

## Minor Pitfalls

Mistakes that cause annoyance but are fixable.

### Pitfall 10: Delta Data Without Subject-Visit Uniqueness Validation

**What goes wrong:** `analyse_mi_data()` validates that delta has the required columns but does not check for duplicate subject-visit combinations. A `dplyr::left_join()` with duplicated keys silently replicates rows, inflating the dataset. The analysis then runs on more rows than intended, producing biased estimates.

**Prevention:**
1. Add uniqueness check: `if (anyDuplicated(delta[, c(vars$subjid, vars$visit)])) stop("delta has duplicate subject-visit rows")`.
2. Place this check immediately after column validation, before the join.

**Phase:** Hardening phase (analyse_mi_data).

**Confidence:** HIGH -- directly observable in `R/analyse_mi_data.R` lines 195-214. Also flagged in CONCERNS.md.

---

### Pitfall 11: String-Based Key Matching with "|||" Separator Masking Type Mismatches

**What goes wrong:** `reduce_imputed_data()` and `expand_imputed_data()` convert subject and visit columns to character and concatenate with `"|||"` separator for matching. If a subject ID literally contains `"|||"`, key collisions occur. More commonly, the `as.character()` conversion masks type mismatches between `imputed_data` and `original_data` (e.g., factor vs. integer subject IDs that happen to have the same string representation but different semantics).

**Prevention:**
1. Use `interaction(col1, col2, drop = TRUE)` or `dplyr::left_join()` on the key columns directly instead of string concatenation.
2. Add explicit type equality checks before the join: warn if `class(imputed_data[[subjid]]) != class(original_data[[subjid]])`.
3. If keeping string concatenation, validate that no values in either key column contain the separator.

**Phase:** Hardening phase (storage functions).

**Confidence:** HIGH -- directly observable in `R/imputation_storage.R` lines 131-134 and lines 266-269.

---

### Pitfall 12: Responder Bar Chart Misrepresenting Proportions with Small Denominators

**What goes wrong:** When creating responder bar charts (proportion responding by arm and visit), developers use raw proportions without considering the denominator. Visits with high dropout have small denominators, making proportions volatile. Displaying proportions without confidence intervals or denominator annotation misleads readers into treating a 4/5 (80%) proportion the same as a 40/50 (80%) proportion.

**Prevention:**
1. Always annotate bars with n/N (count over denominator).
2. Include error bars showing binomial confidence intervals (Wilson or Clopper-Pearson).
3. If using ggplot2, use `geom_col()` (not `geom_bar(stat = "identity")`) to avoid confusion about count vs. proportion encoding.
4. Consider a table-under-the-figure showing n at risk per visit per arm.

**Phase:** Reporting layer (responder bar chart function).

**Confidence:** MEDIUM -- domain knowledge from clinical trial reporting conventions. Not verified against a specific source.

---

### Pitfall 13: tidyr::separate() Deprecation Creating Future Maintenance Burden

**What goes wrong:** `tidy_pool_obj()` uses `tidyr::separate()` which is superseded by `tidyr::separate_wider_delim()`, `tidyr::separate_wider_position()`, and `tidyr::separate_wider_regex()`. While superseded functions will continue to receive critical bug fixes, they will not receive new features and are excluded from tidyr documentation search. Over time, new tidyverse releases may produce deprecation warnings that surface to users.

**Prevention:**
1. Migrate to `tidyr::separate_wider_regex()` with named capture groups, which also solves the underscore-in-parameter-names problem (Pitfall 2).
2. If keeping `tidyr::separate()` for now, add a comment documenting the supersession and linking to the migration path.
3. Pin `tidyr (>= 1.3.0)` in DESCRIPTION to ensure the replacement functions are available.

**Phase:** Hardening phase (tidier/formatter functions). Combine with Pitfall 2 fix.

**Confidence:** HIGH -- `tidyr::separate()` supersession confirmed at [tidyr::separate reference](https://tidyr.tidyverse.org/reference/separate.html).

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation | Severity |
|---|---|---|---|
| ARD conversion (cards) | Column contract violations (Pitfall 1) | Use `cards::as_cards_fn()`, validate output class | Critical |
| ARD conversion (cards) | cardx hidden Suggests deps (Pitfall 8) | Explicit references for renv, `check_installed()` | Moderate |
| gtsummary table templates | Regulatory table layout mismatch (Pitfall 5) | Prototype layout first, consider direct gt | Moderate |
| Forest plot function | Scale/ordering/reference line errors (Pitfall 6) | Parameterize scale, always reverse y-axis | Moderate |
| Responder bar chart | Small denominator misrepresentation (Pitfall 12) | Annotate n/N, add CIs | Minor |
| print/summary methods | S3 dispatch conflicts with rbmi (Pitfall 7) | Use wrapper class for pool, register methods | Moderate |
| Refactor analyse_mi_data | class(method)[[2]] fragility (Pitfall 3) | Wrap rbmi::analyse() or use inherits() | Critical |
| Harden gcomp functions | beeca column coupling (Pitfall 4) | Version pin, validate output schema | Critical |
| Harden gcomp functions | Formula injection (Pitfall 9) | Sanitize covariate names, validate against data | Moderate |
| Harden storage functions | String key collision (Pitfall 11) | Use interaction() or direct join | Minor |
| Harden tidier/formatter | Parameter parsing fragility (Pitfall 2) | Regex parsing or metadata attributes | Critical |
| Harden tidier/formatter | tidyr::separate() supersession (Pitfall 13) | Migrate to separate_wider_regex() | Minor |
| Harden analyse_mi_data | Delta uniqueness gap (Pitfall 10) | Add uniqueness assertion | Minor |
| Dependency management | Suggests break CRAN check (Pitfall 8) | Guard all calls, run --no-suggests in CI | Moderate |

## Recommended Phase Ordering Based on Pitfalls

1. **Hardening first** -- Fix Pitfalls 2, 3, 4, 9, 10, 11 before building on top of fragile foundations. The parameter parsing (Pitfall 2) and class hierarchy (Pitfall 3) issues will propagate into ARD conversion if not fixed first.
2. **Refactor analyse_mi_data** -- Wrapping `rbmi::analyse()` eliminates Pitfall 3 entirely and reduces the internal helper surface area.
3. **ARD conversion** -- With hardened parsing and stable class hierarchy, convert tidy results to ARD format (Pitfall 1 guidance).
4. **Reporting layer** -- Build gtsummary templates (Pitfall 5), forest plots (Pitfall 6), and bar charts (Pitfall 12) on top of correct ARDs.
5. **print/summary improvements** -- Low risk of blocking other work (Pitfall 7).

## Sources

- [cards package documentation](https://insightsengineering.github.io/cards/main/) -- ARD structure, `as_cards_fn()`
- [cardx package](https://github.com/insightsengineering/cardx) -- hidden Suggests dependency warning
- [gtsummary ARD-first Tables](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html) -- `tbl_ard_summary()` requirements
- [CDISC COSA Spotlight ARD + gtsummary 2025](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/) -- clinical trial ARD workflow
- [R Packages (2e) Dependencies in Practice](https://r-pkgs.org/dependencies-in-practice.html) -- Suggests/Imports guidance, CRAN policy
- [rlang is_installed reference](https://rlang.r-lib.org/reference/is_installed.html) -- conditional dependency checks
- [Advanced R: S3](https://adv-r.hadley.nz/s3.html) -- S3 dispatch, method conflicts
- [tidyr::separate reference](https://tidyr.tidyverse.org/reference/separate.html) -- supersession status
- [rbmi CRAN news](https://cran.r-project.org/web/packages/rbmi/news/news.html) -- version history 1.4 -> 1.5
- [KHstats Forest Plots](https://www.khstats.com/blog/forest-plots/) -- clinical trial forest plot conventions
- [forestplot vignette](https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html) -- reference line and scale conventions
- [R Consortium: Supercharging Statistical Analysis with ARDs](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/) -- ARD ecosystem overview
- [PHUSE 2025 ARD Workshop](https://www.danieldsjoberg.com/ARD-PHUSE-workshop-2025/) -- ARD practical guidance
- Codebase analysis: `R/tidiers.R`, `R/analyse_mi_data.R`, `R/analysis_utils.R`, `R/imputation_storage.R`, `R/utils.R`, `R/data_helpers.R`, `R/formatting.R`
- Existing `.planning/codebase/CONCERNS.md` -- confirmed fragilities

---

*Pitfalls research: 2026-02-07*
