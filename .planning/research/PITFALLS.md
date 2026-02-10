# Domain Pitfalls: v3 ARD Enrichment & Polish

**Domain:** Adding MI-specific ARD metadata, imputation diagnostics, and publication polish to an existing clinical trial R package (rbmiUtils)
**Researched:** 2026-02-10
**Overall confidence:** MEDIUM-HIGH (grounded in codebase analysis, rbmi source code, cards validation source, and verified ecosystem documentation)

**Scope:** This document covers pitfalls specific to the v3 milestone features. For pitfalls related to the foundation layer (parameter parsing, class hierarchy, beeca coupling, etc.), see the v1/v2 PITFALLS.md. Many of those were resolved in v1 and v2 and are not repeated here.

---

## Critical Pitfalls

Mistakes that cause existing functionality to break or produce incorrect results.

### Pitfall 1: Adding FMI/Pooling Metadata Rows Breaks check_ard_structure() Validation

**What goes wrong:** The existing `pool_to_ard()` produces a valid ARD that passes `cards::check_ard_structure()`. When adding new MI-specific statistics (FMI, relative increase in variance, between-imputation variance, within-imputation variance, degrees of freedom), developers add new rows with custom `stat_name` values (e.g., `"fmi"`, `"riv"`, `"var_between"`, `"var_within"`, `"df"`). If these rows have incorrect column types -- particularly if `stat`, `fmt_fun`, `warning`, or `error` columns are not list-columns in the new rows -- the entire ARD fails validation. The `check_ard_structure()` function validates that `stat`, `fmt_fun`, `warning`, and `error` are all list-type columns. A single row with a non-list value in any of these columns corrupts the entire data frame.

**Why it happens:** The current `pool_to_ard()` constructs ARD rows using `I(list(...))` to create list-columns. When developers add new rows for FMI metadata, they may construct them separately and `rbind()` them to the existing ARD. If the new rows use plain numeric values instead of list-wrapped values for `stat`, or `NULL` instead of `list(NULL)` for warning/error, the `rbind()` silently coerces the list-columns to non-list types or produces inconsistent column types.

**Consequences:** The returned ARD object fails `cards::check_ard_structure()`, breaking the documented contract (the existing roxygen says "The resulting ARD passes cards::check_ard_structure() validation"). Downstream gtsummary consumers silently drop rows or error. Existing tests for `pool_to_ard()` that check `cards::check_ard_structure()` fail, alerting to the regression, but only if those tests run (they require `cards` to be installed).

**Warning signs:**
- `rbind()` producing data frames where `is.list(ard$stat)` returns `FALSE`
- `check_ard_structure()` errors mentioning column type mismatches
- New ARD has more rows than expected but gtsummary displays only the original statistics
- `str(ard$stat)` shows a mix of list and non-list entries

**Prevention:**
1. Build all new metadata rows using the same `data.frame(..., stat = I(list(...)))` pattern as the existing `pool_to_ard()` code (lines 100-117 in `R/ard_conversion.R`).
2. After adding new rows, verify the ARD structure before returning: `cards::check_ard_structure(ard_df)` should be called as the last step (it already is on line 122 via `cards::as_card()`).
3. Add a test specifically for the enriched ARD: after adding FMI rows, assert `is.list(ard$stat)`, `is.list(ard$warning)`, `is.list(ard$error)`, and `is.list(ard$fmt_fun)` remain TRUE.
4. Extra columns beyond the required set are permitted by `check_ard_structure()` -- it only validates that required columns exist and have correct types. So additional metadata can also be added as new columns rather than new rows if appropriate.
5. Use `cards::bind_ard()` instead of `rbind()` if available -- it is designed to safely combine ARD objects while preserving list-column types.

**Detection:** The existing test `test-ard_conversion.R` lines 25-49 already test `check_ard_structure()`. Add a parallel test for the enriched version.

**Phase:** ARD metadata enrichment. Must be the first feature implemented since it modifies the core ARD output.

**Confidence:** HIGH -- `check_ard_structure()` required columns and type enforcement verified from [cards source code](https://raw.githubusercontent.com/insightsengineering/cards/main/R/check_ard_structure.R): required columns are `"variable"`, `"stat_name"`, `"stat_label"`, `"stat"`, `"fmt_fun"`, `"warning"`, `"error"`, and it validates that `stat`, `fmt_fun`, `warning`, `error`, and level columns are list-type.

---

### Pitfall 2: FMI/Variance Components Not Available in rbmi Pool Object -- Must Be Recomputed

**What goes wrong:** Developers assume the rbmi `pool()` return object contains FMI (fraction of missing information), between-imputation variance (`var_b`), and within-imputation variance (`var_w`) as accessible fields. In reality, `rbmi::pool()` returns only: `pars` (with `est`, `ci`, `se`, `pvalue` per parameter), `conf.level`, `alternative`, `N`, and `method`. The variance components `var_w` and `var_b` are computed internally in `rbmi::rubin_rules()` but not returned. FMI is never computed at all.

The consequence: developers try `pool_obj$fmi` or `pool_obj$pars$trt_Week4$var_b` and get `NULL`, then either (a) skip the feature, (b) compute FMI incorrectly, or (c) try to access rbmi internals that are not part of the public API.

**Why it happens:** Other MI packages (notably `mice::pool()`) return FMI, lambda, RIV, and degrees of freedom as standard output fields. Developers familiar with `mice` expect `rbmi::pool()` to behave similarly. The rbmi documentation for `rubin_rules()` confirms it returns only `est_point`, `var_t`, and `df` -- not the decomposed variance components.

**Consequences:** If developers use `rbmi::rubin_rules()` directly to recompute the variance components, they need the per-imputation estimates and standard errors, which are available in the `analysis` object's `results` list but require knowledge of rbmi's internal structure. Getting this wrong means incorrect FMI values in the ARD, which would misrepresent the quality of the imputation to downstream consumers.

**Prevention:**
1. Accept that FMI must be computed from scratch using the `analysis` object (not the `pool` object). The formulas from Rubin (1987):
   - `var_w = mean(se_i^2)` (average within-imputation variance)
   - `var_b = var(est_i)` (between-imputation variance)
   - `var_t = var_w + (1 + 1/M) * var_b` (total variance)
   - `riv = (1 + 1/M) * var_b / var_w` (relative increase in variance)
   - `lambda = (1 + 1/M) * var_b / var_t` (proportion of total variance due to missingness)
   - `fmi = (riv + 2/(df+3)) / (riv + 1)` (fraction of missing information, adjusted)
2. Pass the `analysis` object alongside the `pool` object to the enriched `pool_to_ard()`, OR compute and store metadata before pooling.
3. Use `rbmi::rubin_rules()` (exported function) to validate the computation: it accepts vectors of estimates and SEs and returns `est_point`, `var_t`, and `df`. Compare your computed `var_t` against its output to verify correctness.
4. Document clearly that the FMI computation follows Barnard and Rubin (1999) with small-sample degrees of freedom adjustment, matching rbmi's internal approach.
5. Consider whether `pool_to_ard()` should accept an optional `analysis_obj` parameter for enrichment, or whether a separate `enrich_ard()` function should add metadata post-hoc.

**Warning signs:**
- `pool_obj$fmi` returns NULL
- FMI values outside [0, 1] (computation error)
- FMI values that don't correlate with missing data percentages (possible but check if extreme)
- `var_t` from manual computation doesn't match `pool_obj$pars[[x]]$se^2`

**Phase:** ARD metadata enrichment. Fundamental design decision -- must be resolved before implementation begins.

**Confidence:** HIGH -- rbmi `pool()` return structure verified from [rbmi pool.R source](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R). `rubin_rules()` return values verified from [rbmi documentation](https://www.rdocumentation.org/packages/rbmi/versions/1.5.2/topics/rubin_rules). FMI not computed anywhere in rbmi confirmed from source code review.

---

### Pitfall 3: describe_draws() Accesses rbmi Internal Object Structure That May Change

**What goes wrong:** Creating `describe_draws()` requires accessing the internal structure of rbmi's `draws` object. The draws object is a named list with: `data` (an R6 `longdata` object), `method` (method specification), `samples` (list of per-imputation parameter estimates containing `ids`, `beta`, `sigma`, `theta`, `failed`, `ids_samp`), `fit` (Stan fit object for Bayesian methods, NULL otherwise), `n_failures` (count of failed fits), and `formula` (the model formula). None of these fields have a stability guarantee -- they are implementation details of rbmi, not a documented public API.

**Why it happens:** rbmi exposes high-level functions (`draws()`, `impute()`, `analyse()`, `pool()`) but the intermediate objects are designed for internal consumption. The `draws` object is passed to `impute()`, not examined by users. Creating diagnostic helpers that peek inside these objects couples rbmiUtils to rbmi's internals.

**Consequences:** When rbmi releases a new version (e.g., 1.5 -> 2.0), the internal structure of the draws object may change. If `describe_draws()` accesses `draws_obj$samples[[1]]$sigma` and rbmi renames this field or restructures the samples list, `describe_draws()` silently returns wrong values or errors. Since rbmiUtils pins `rbmi (>= 1.4)` with no upper bound, this breakage happens automatically when users update rbmi.

**Prevention:**
1. Minimize the number of internal fields accessed. For `describe_draws()`, the safest information to extract is:
   - `draws_obj$method` (method class/type) -- relatively stable since it's user-provided
   - `draws_obj$n_failures` (number of failed fits) -- simple scalar
   - `draws_obj$formula` (model formula) -- standard R formula object
   - `length(draws_obj$samples)` (number of draws) -- list length
   - `draws_obj$fit` (Stan fit for Bayesian) -- presence/absence, not internals
2. For Bayesian methods: if `draws_obj$fit` is a stanfit object, use `rstan::summary()` or `rstan::get_sampler_params()` for convergence diagnostics (Rhat, n_eff). These are rstan public API, not rbmi internals.
3. Wrap every internal field access in a tryCatch or null-check: `tryCatch(draws_obj$samples[[1]]$sigma, error = function(e) NULL)`. Return "unavailable" for fields that don't exist rather than erroring.
4. Add defensive version checking: `if (packageVersion("rbmi") >= "2.0.0") cli::cli_warn("describe_draws() may need updating for rbmi >= 2.0")`.
5. Test with at least two rbmi versions in CI (the minimum pinned version and the latest CRAN version).
6. Use `inherits()` for method type detection (already established pattern from v1 hardening) rather than field-based detection.

**Warning signs:**
- Error messages like "$ operator is invalid for atomic vectors" or "subscript out of bounds"
- `describe_draws()` returning NULL or NA for fields that should have values
- Different behavior between rbmi 1.4 and 1.5

**Phase:** Diagnostic helpers. Implement after ARD metadata (which has a safer API surface).

**Confidence:** HIGH -- draws object structure verified from [rbmi draws.R source](https://github.com/insightsengineering/rbmi/blob/main/R/draws.R) and [as_draws documentation](https://www.rdocumentation.org/packages/rbmi/versions/1.4.1/topics/as_draws). The lack of stability guarantee is confirmed by the absence of any `@return` documentation specifying the fields as public API.

---

### Pitfall 4: describe_imputation() Coupling to rbmi extract_imputed_dfs() Internal Format

**What goes wrong:** Creating `describe_imputation()` to summarize imputation quality requires accessing individual imputed datasets from the `impute` object. The current `get_imputed_data()` function already does this via `rbmi::extract_imputed_dfs(impute_obj, idmap = TRUE)`, which returns a list of data frames. However, diagnostic summaries (number of imputations, imputed value ranges, convergence of imputed distributions) require iterating over potentially hundreds of complete datasets, which is memory-intensive. Developers may try to access `impute_obj` internals directly to avoid the full extraction cost.

**Why it happens:** The `impute` object stores imputations compactly. `rbmi::extract_imputed_dfs()` materializes all M complete datasets into memory. For diagnostic summaries, you only need summary statistics across imputations for originally-missing values, not the full datasets. But the rbmi API only provides the full extraction path.

**Consequences:** Memory exhaustion with large datasets and many imputations. For example, 200 subjects x 8 visits x 1000 imputations = 1.6M rows, each with all covariates. A `describe_imputation()` that extracts all datasets into memory just to compute summary statistics is wasteful and may fail on typical clinical trial analyst laptops.

**Prevention:**
1. For basic diagnostics (number of imputations, method used), access metadata from the `impute` object's top-level fields rather than extracting datasets: the number of imputations is `length(rbmi::extract_imputed_dfs(impute_obj))`, but a cheaper approach is inspecting the object structure directly.
2. For value-level diagnostics, extract a small sample of imputations (e.g., first 5) rather than all M: `imputed_dfs <- rbmi::extract_imputed_dfs(impute_obj)[1:5]`.
3. If `get_imputed_data()` has already been called and the full stacked data frame exists, compute diagnostics from that (it's already in memory). Design `describe_imputation()` to accept either an `impute` object or a stacked data frame.
4. Document memory requirements clearly: "For M=1000 imputations, describe_imputation() will extract N datasets into memory. Consider using describe_imputation(impute_obj, n_sample = 10) for large analyses."
5. Consider whether `describe_imputation()` should accept the pre-existing `ADMI` stacked data frame (which is already the standard workflow input to `analyse_mi_data()`) rather than the raw impute object.

**Phase:** Diagnostic helpers. Design the function signature carefully before implementation.

**Confidence:** MEDIUM -- the memory concern is inferred from the structure of `get_imputed_data()` (lines 52-82 in `R/utils.R`) and typical clinical trial dataset sizes. The `rbmi::extract_imputed_dfs()` API is stable and documented. Memory behavior under many imputations is not directly tested.

---

## Moderate Pitfalls

Mistakes that cause delays, visual regressions, or technical debt.

### Pitfall 5: Enriched pool_to_ard() API Signature Change Breaks Existing Callers

**What goes wrong:** To compute FMI and variance components, `pool_to_ard()` needs access to per-imputation results from the `analysis` object. The current signature is `pool_to_ard(pool_obj, conf.level = NULL)`. Adding a required `analysis_obj` parameter is a breaking API change. Making it optional means the enriched metadata is only available when both objects are passed, creating two code paths with different outputs. Users who upgrade rbmiUtils and call `pool_to_ard(pool_obj)` get different output (missing FMI rows) than users who call `pool_to_ard(pool_obj, analysis_obj)`.

**Why it happens:** The `pool` object is the natural input for ARD conversion (it has the pooled estimates). But FMI requires the pre-pooled per-imputation estimates. These live in the `analysis` object. The information was lost during `rbmi::pool()` -- pooling is a summarization that discards the per-imputation detail.

**Prevention:**
1. Make `analysis_obj` an optional parameter with a default of `NULL`. When NULL, return the existing 5-statistic ARD (backward compatible). When provided, add FMI/variance metadata rows.
2. Document the enrichment clearly: "Pass the analysis object to include MI diagnostic metadata (FMI, between/within variance). Without it, only pooled statistics are included."
3. Add a lifecycle badge: `@details ## MI Diagnostic Metadata \lifecycle{experimental}` to signal that the enrichment is new and the metadata row names may change.
4. Consider an alternative approach: a separate `ard_add_mi_metadata()` function that takes an existing ARD and an analysis object, adding rows. This keeps `pool_to_ard()` unchanged and follows the cards philosophy of composable transformations (similar to `cards::add_calculated_row()`).
5. Never add required parameters to an existing exported function -- always optional with backward-compatible defaults.

**Phase:** ARD metadata enrichment. Design decision needed before implementation.

**Confidence:** HIGH -- API design patterns verified from [R Packages (2e)](https://r-pkgs.org/). The `cards::add_calculated_row()` pattern confirmed from [cards documentation](https://cran.r-project.org/web/packages/cards/refman/cards.html).

---

### Pitfall 6: gt Table Styling Changes Break Existing Visual Output

**What goes wrong:** Modifying `efficacy_table()` styling (fonts, spacing, alignment, borders, header formatting) changes the visual output. Existing pre-rendered images in `man/figures/` (used in README, roxygen examples, and pkgdown site) no longer match the function output. Users who have screenshots or reports built with the old styling see visual regressions. More critically, snapshot tests (if any exist via `testthat::expect_snapshot()`) fail on the first run and need manual approval.

**Why it happens:** gt tables render to HTML, and their visual appearance depends on CSS, font availability, and rendering engine. The current `efficacy_table()` uses gt defaults with some customization (source notes, column alignment, number formatting). Adding styling refinements (row striping, header borders, font stacks, compact spacing) changes every rendered table.

**Consequences:**
- Pre-rendered images (`man/figures/efficacy_table-example.png`, `man/figures/README-efficacy-table-1.png`) become stale
- pkgdown site shows old images until rebuilt
- R CMD check may warn about changed man page content if images are regenerated but not committed
- Users' existing scripts produce different-looking tables after upgrading

**Warning signs:**
- `man/figures/` images look different from live `efficacy_table()` output
- pkgdown build produces different site than committed `docs/` directory
- User reports of "tables look different after update"

**Prevention:**
1. Add styling parameters to `efficacy_table()` rather than changing defaults: e.g., `compact = FALSE` (default preserves current look), `theme = "default"` (new themes available but not forced).
2. If changing defaults, increment the package version to signal the visual change and document it in NEWS.md.
3. Regenerate all pre-rendered images after styling changes: `man/figures/efficacy_table-example.png`, `man/figures/README-efficacy-table-1.png`, and any images in vignette outputs.
4. Use `gt::opt_table_font()` with system font stacks (e.g., `"rounded-sans"`) rather than specific font names to ensure cross-platform consistency.
5. Test rendered output on both macOS and Windows if possible -- font rendering differs significantly. gt tables rendered to PDF/LaTeX have different styling capabilities than HTML; be cautious about promising "publication quality" across all output formats.
6. Consider whether styling changes should be a new function (e.g., `efficacy_table_styled()`) or a `theme` parameter on the existing function.

**Phase:** Table/plot polish. Implement after ARD metadata (which is functionally more important than visual polish).

**Confidence:** MEDIUM -- gt cross-platform rendering differences verified from [gt font issues](https://github.com/rstudio/gt/issues/337) and [gt/LaTeX rendering discussion](https://forum.posit.co/t/rmarkdown-gt-package-rendering-different-in-html-and-pdf/140864). Pre-rendered image staleness confirmed from `man/figures/` directory listing.

---

### Pitfall 7: Forest Plot Patchwork Alignment Breaks with Styling Changes

**What goes wrong:** The current `plot_forest()` uses a three-panel patchwork composition with width ratios `c(3, 4, 1.5)`. Changing text size, font family, point size, or adding new annotation elements (e.g., FMI annotation, sample size labels) disrupts the panel alignment. Patchwork aligns panels based on plot area dimensions, but text elements in `theme_void()` panels (the left table and right p-value panels) don't participate in ggplot2's standard alignment system. When text becomes wider or taller, it overflows the allocated panel width.

**Why it happens:** The three-panel forest plot is a creative hack: the left and right panels are not true ggplot2 plots but text annotation panels using `geom_text()` on a `theme_void()` canvas. Patchwork aligns the plot areas, but `theme_void()` panels have no axes to align. When text content changes size, the layout breaks because patchwork cannot negotiate space between panels that have no alignment anchors. The `free()` function in patchwork can help but introduces its own complexity.

**Consequences:** Forest plot text overlaps with the forest panel, p-value labels are clipped, or excessive whitespace appears between panels. The issue is resolution/DPI-dependent: a plot that looks fine at screen resolution may have alignment issues when exported at 300 DPI for publication.

**Warning signs:**
- Text labels truncated or overlapping with the forest plot whiskers
- Plot looks correct in RStudio viewer but wrong in saved PNG/PDF
- Different appearance at different `width` and `height` values in `ggsave()`
- Alignment breaks when `text_size` or `point_size` parameters are changed from defaults

**Prevention:**
1. Any styling changes to `plot_forest()` should be tested at multiple output sizes: `ggsave("test.png", width = 8, height = 4)`, `ggsave("test.png", width = 12, height = 6)`, and `ggsave("test.png", width = 6, height = 3)`.
2. Expose `widths` as a parameter: `widths = c(3, 4, 1.5)` so users can adjust panel proportions without modifying the function.
3. For text panels, use `ggplot2::scale_x_continuous(expand = expansion(mult = c(0.3, 0.05)))` (already done in current code) -- verify this expansion is sufficient for longer text labels from v3 changes.
4. Prefer relative sizing (`text_size = 3` in points) over absolute sizing. Avoid pixel-based sizing.
5. If adding new annotation panels (e.g., FMI column), test with the two-panel layout (`show_pvalues = FALSE`) as well as three-panel.
6. Regenerate `man/figures/plot_forest-trt.png` and `man/figures/README-forest-plot-1.png` after any visual changes.

**Phase:** Plot polish. Can be implemented in parallel with table polish.

**Confidence:** MEDIUM -- patchwork alignment issues verified from [patchwork documentation](https://patchwork.data-imaginist.com/articles/guides/layout.html) and [patchwork issues](https://github.com/thomasp85/patchwork/issues/276). The specific three-panel layout fragility is directly observable in `R/plot_forest.R` lines 354-361.

---

### Pitfall 8: README/Documentation Changes Break pkgdown Build or GitHub Rendering

**What goes wrong:** Overhauling the README to add realistic clinical trial examples introduces several failure modes:

1. **Image path inconsistency:** pkgdown relocates `man/figures/` images to `reference/figures/` on the site. Image paths that work in GitHub README rendering (`man/figures/image.png`) may break on the pkgdown site, or vice versa.
2. **README.Rmd rendering:** If using `README.Rmd`, `usethis::use_readme_rmd()` sets `fig.path = "man/figures/README-"`. pkgdown expects a pre-rendered `README.md`. If `README.Rmd` is edited but not knitted, the site shows stale content.
3. **Example code in README runs during pkgdown build:** Code blocks in `README.Rmd` execute during knitting. If they depend on rbmi (which requires Stan for Bayesian methods), the README cannot be built on CI without a full Stan installation. This is why the current README uses pre-rendered images.
4. **HTML in README.md:** GitHub renders a subset of HTML. The current README uses `<figure>` tags. pkgdown renders the full README as HTML. Some HTML tags supported by pkgdown are not supported by GitHub's markdown renderer, and vice versa.
5. **Vignette code changes:** If documentation examples are moved from vignettes to README or vice versa, resource paths (images, data) may break because pkgdown processes `articles/` and `index.html` differently.

**Prevention:**
1. Keep using pre-rendered static images for README. Do NOT add code that calls `rbmi::draws()` in README.Rmd chunks -- it requires Stan/MCMC and is too slow for CI.
2. Test both rendering targets: `pkgdown::build_home()` for the site, and visually inspect the GitHub rendering of `README.md`.
3. Use relative paths from the package root for images: `man/figures/image.png`. This works for both GitHub and pkgdown.
4. After any README change, run `pkgdown::build_site()` locally and check `docs/index.html` for broken images or layout issues before committing.
5. Do NOT use markdown image resizing syntax (`![](image.png =WIDTHxHEIGHT)`) -- it breaks on pkgdown. Use HTML `<img>` tags with `width` attributes instead.
6. If switching from `README.Rmd` to `README.md`, remove the `README.Rmd` file entirely to avoid confusion about which file is authoritative.
7. Any new images must go in `man/figures/` -- pkgdown only copies images from `man/figures/` and `vignettes/` to the built site.

**Phase:** Documentation overhaul. Implement after functional changes (ARD metadata, diagnostic helpers) to avoid re-doing documentation.

**Confidence:** HIGH -- pkgdown image path handling verified from [pkgdown issues #1472](https://github.com/r-lib/pkgdown/issues/1472), [#1943](https://github.com/r-lib/pkgdown/issues/1943), [#2194](https://github.com/r-lib/pkgdown/issues/2194), and [#2389](https://github.com/r-lib/pkgdown/issues/2389). GitHub HTML rendering limitations confirmed from [R Packages (2e) chapter 18](https://r-pkgs.org/other-markdown.html). Current README structure confirmed from codebase (`README.md` lines 51-62 use `<figure>` tags, images in `man/figures/`).

---

### Pitfall 9: FMI Computation Requires Per-Parameter Estimates from analysis Object Internals

**What goes wrong:** To compute FMI for each pooled parameter, you need the per-imputation estimate (`est_i`) and standard error (`se_i`) for that parameter across all M imputations. These live in `analysis_obj$results` -- a list of M elements, each containing named parameter results with `$est` and `$se` fields. Extracting these requires iterating over `analysis_obj$results` and matching parameter names, which couples to rbmi's internal analysis result structure.

**Why it happens:** The `analysis` object is created by `as_analysis2()` (lines 295-345 in `R/analyse_mi_data.R`) which stores results as `list(results = rbmi::as_class(results, c(next_class, "list")))`. Each element of `results` is a named list where names are parameter names (e.g., `"trt_Week4"`, `"lsm_ref_Week4"`) and values are lists with `$est` and `$se`. This structure is how `rbmi::pool()` consumes the results -- so it's relatively stable, but it's still an internal contract.

**Consequences:** If the extraction logic has an off-by-one error or mismatches parameter names between the analysis and pool objects, FMI values are assigned to the wrong parameters. Since FMI values are plausible-looking numbers (between 0 and 1), the error is silent and may not be caught without explicit validation.

**Prevention:**
1. Extract per-imputation values by parameter name, not by position: `sapply(analysis_obj$results, function(r) r[["trt_Week4"]]$est)`.
2. Validate that parameter names in the analysis object match those in the pool object: `names(analysis_obj$results[[1]])` should match the parameter names in `tidy_pool_obj(pool_obj)$parameter`.
3. Compute FMI for one known parameter first and verify against manual calculation on the same data.
4. Add a sanity check: `var_t_manual = var_w + (1 + 1/M) * var_b` should equal `pool_obj$pars[[param]]$se^2` (within floating point tolerance). If they don't match, the extraction is wrong.
5. Consider wrapping the extraction in a standalone function (e.g., `extract_per_imputation_stats(analysis_obj, param_name)`) with clear error messages when the structure doesn't match expectations.

**Phase:** ARD metadata enrichment. Part of the FMI computation design.

**Confidence:** HIGH -- analysis object structure directly observable in `R/analyse_mi_data.R` and confirmed by `make_mock_pool()` helper in test file `test-ard_conversion.R` which shows the parameter naming convention.

---

### Pitfall 10: New pkgdown Reference Groups for Diagnostic Functions May Break Existing Site Layout

**What goes wrong:** Adding `describe_draws()` and `describe_imputation()` to the package requires adding them to `_pkgdown.yml` reference groups. If they are added to the wrong group or a new group is created without updating the group list, pkgdown's `build_reference()` either omits them (they appear in an "Other" section) or errors with "Topic not found" if the function name doesn't match the roxygen export.

**Why it happens:** The current `_pkgdown.yml` has 9 reference groups (Data Preparation, Analysis, Tidying, Reporting, Formatting, Storage, Utilities, Print & Summary Methods, Datasets). The new diagnostic functions don't fit neatly into any existing group. If a "Diagnostics" group is added, the YAML indentation must be exact -- a single space error causes silent drops.

**Prevention:**
1. Add a new "Diagnostics" reference group between "Tidying" and "Reporting" in `_pkgdown.yml`.
2. After modifying `_pkgdown.yml`, run `pkgdown::build_reference_index()` locally and verify all exported functions appear.
3. Check that the YAML is valid: `yaml::read_yaml("_pkgdown.yml")` should not error.
4. Match function names exactly: if the function is `describe_draws`, the contents entry must be `describe_draws` (not `describe_draws()` or `describe-draws`).

**Phase:** Documentation overhaul. Must coordinate with the implementation of diagnostic functions.

**Confidence:** HIGH -- `_pkgdown.yml` structure confirmed from codebase (lines 39-92). pkgdown reference grouping behavior documented at [pkgdown reference](https://pkgdown.r-lib.org/articles/pkgdown.html).

---

## Minor Pitfalls

Mistakes that cause annoyance or cosmetic issues but are fixable.

### Pitfall 11: gt Source Notes Ordering Affects Table Footer Readability

**What goes wrong:** The current `efficacy_table()` adds three source notes in sequence: pooling method, number of imputations, and confidence level (lines 224-226). Adding FMI summary or additional MI metadata to the footer requires careful ordering. `gt::tab_source_note()` adds notes in the order called, with the last call appearing at the bottom. If MI-specific metadata (e.g., "Average FMI: 0.23") is added first, it appears above the method note, which is confusing. If added last, it appears at the bottom, which may be overlooked.

**Prevention:**
1. Define a consistent note ordering: (1) Method, (2) N imputations, (3) Confidence level, (4) MI diagnostics (FMI, etc.).
2. Use `gt::tab_footnote()` for parameter-specific annotations (e.g., per-visit FMI) rather than `tab_source_note()` which applies to the whole table.
3. Consider using a single source note for MI metadata: "Pooling: Rubin's rules | M = 100 | 95% CI | Avg FMI = 0.23" to reduce footer clutter.

**Phase:** Table polish.

**Confidence:** HIGH -- gt source note ordering behavior confirmed from [gt tab_source_note documentation](https://gt.rstudio.com/reference/tab_options.html).

---

### Pitfall 12: Bootstrap and Conditional Mean Methods Return NA for Standard Errors

**What goes wrong:** When computing FMI for non-Bayesian pooling methods (bootstrap percentile, jackknife), the standard error may be `NA` in the pool object. The bootstrap percentile method computes CIs directly from quantiles without a standard error. If FMI computation assumes `se` is always available, it will produce `NaN` or error for bootstrap analyses.

**Prevention:**
1. Check `pool_obj$method` before attempting FMI computation. FMI from Rubin's rules is only meaningful for `method = "rubin"` (Bayesian and approximate Bayesian).
2. For bootstrap/jackknife methods, either skip FMI computation entirely or compute a bootstrap-based analog (variance of estimates across resamples).
3. Return `NA` for FMI when the pooling method doesn't support it, and document why.
4. Test FMI computation with all three method types: `method_bayes()`, `method_approxbayes()`, and `method_condmean()`.

**Phase:** ARD metadata enrichment.

**Confidence:** MEDIUM -- bootstrap SE behavior inferred from rbmi pool.R source code where method-specific handling produces NA for standard errors with percentile-based methods. Needs verification with actual rbmi output.

---

### Pitfall 13: Snapshot Tests for Visual Output Are Platform-Dependent

**What goes wrong:** Adding snapshot tests (`testthat::expect_snapshot()`) for `describe_draws()` or `describe_imputation()` console output captures platform-specific details like Unicode rendering, cli formatting codes, and terminal width. A snapshot created on macOS will fail on Windows CI if Unicode characters differ or terminal width assumptions change.

**Prevention:**
1. Use `expect_snapshot()` with `variant = .Platform$OS.type` if platform-specific output is expected.
2. Prefer structural assertions over snapshot tests: `expect_equal(result$n_imputations, 100)` rather than `expect_snapshot(print(result))`.
3. If using snapshot tests for cli-formatted output, set `cli::cli_inform()` width explicitly: `withr::local_options(cli.width = 80)`.
4. For gt table snapshots, use `gt::as_raw_html()` and compare string content rather than relying on visual rendering.

**Phase:** All phases with new console or visual output.

**Confidence:** MEDIUM -- testthat snapshot platform sensitivity documented in [testthat 3e snapshot testing](https://testthat.r-lib.org/articles/snapshotting.html).

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation | Severity |
|---|---|---|---|
| ARD metadata enrichment | List-column corruption from rbind (Pitfall 1) | Use same `I(list(...))` pattern, validate with `check_ard_structure()` | Critical |
| ARD metadata enrichment | FMI not in pool object (Pitfall 2) | Compute from analysis object, validate against `rubin_rules()` | Critical |
| ARD metadata enrichment | API signature change (Pitfall 5) | Optional `analysis_obj` param or separate `enrich_ard()` function | Moderate |
| ARD metadata enrichment | Per-imputation extraction fragility (Pitfall 9) | Extract by name not position, validate parameter name matching | Moderate |
| ARD metadata enrichment | Bootstrap/condmean SE is NA (Pitfall 12) | Check method before FMI, return NA for non-Rubin methods | Minor |
| Diagnostic helpers: describe_draws() | rbmi draws internal structure (Pitfall 3) | Defensive access with tryCatch, minimize fields accessed | Critical |
| Diagnostic helpers: describe_imputation() | Memory from full extraction (Pitfall 4) | Sample imputations, accept pre-extracted data | Moderate |
| Table polish | Visual regression from styling changes (Pitfall 6) | Add parameters not change defaults, regenerate images | Moderate |
| Table polish | Source note ordering (Pitfall 11) | Define consistent ordering convention | Minor |
| Plot polish | Patchwork alignment with new elements (Pitfall 7) | Test at multiple output sizes, expose `widths` param | Moderate |
| Documentation overhaul | Image paths, pkgdown build breakage (Pitfall 8) | Pre-rendered images, test both GitHub and pkgdown rendering | Moderate |
| Documentation overhaul | pkgdown reference group changes (Pitfall 10) | Validate YAML, run build_reference_index() locally | Minor |
| All phases | Snapshot test platform sensitivity (Pitfall 13) | Prefer structural assertions over snapshots | Minor |

---

## Recommended Phase Ordering Based on Pitfalls

1. **ARD metadata enrichment first** -- Pitfalls 1, 2, 5, 9, 12 cluster here. The fundamental design decision (how to get per-imputation stats for FMI, API signature for enriched `pool_to_ard()`) must be resolved before building on top. This is the highest-risk area because it modifies the core ARD output that passes `check_ard_structure()`.

2. **Diagnostic helpers second** -- Pitfalls 3 and 4 require careful API design. `describe_draws()` couples to rbmi internals (inherently fragile). `describe_imputation()` has memory considerations. Design function signatures before implementing.

3. **Table/plot polish third** -- Pitfalls 6 and 7 are moderate risk but visually impactful. Implement AFTER functional changes to avoid regenerating images twice. Add styling as new parameters, not changed defaults.

4. **Documentation overhaul last** -- Pitfalls 8, 10 are lower risk and should be done after all functional changes are complete. Regenerate all images, update `_pkgdown.yml`, rebuild site.

**Rationale:** Functional correctness (ARD, diagnostics) before visual polish (tables, plots) before documentation. Each layer depends on the previous one being stable.

---

## Sources

- [cards check_ard_structure source](https://raw.githubusercontent.com/insightsengineering/cards/main/R/check_ard_structure.R) -- required columns: `variable`, `stat_name`, `stat_label`, `stat`, `fmt_fun`, `warning`, `error`; list-column enforcement on `stat`, `fmt_fun`, `warning`, `error`
- [rbmi pool.R source](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R) -- pool object returns `pars`, `conf.level`, `alternative`, `N`, `method`; `var_w` and `var_b` computed internally but not returned; FMI never computed
- [rbmi rubin_rules documentation](https://www.rdocumentation.org/packages/rbmi/versions/1.5.2/topics/rubin_rules) -- returns `est_point`, `var_t`, `df` only
- [rbmi draws.R source](https://github.com/insightsengineering/rbmi/blob/main/R/draws.R) -- draws object structure: `data`, `method`, `samples`, `fit`, `n_failures`, `formula`
- [rbmi as_draws documentation](https://www.rdocumentation.org/packages/rbmi/versions/1.4.1/topics/as_draws) -- samples subcomponents: `ids`, `beta`, `sigma`, `theta`, `failed`, `ids_samp`
- [Rubin's Rules reference](https://bookdown.org/mwheymans/bookmi/rubins-rules.html) -- FMI, lambda, RIV formulas
- [Barnard-Rubin (1999) df adjustment](https://rdrr.io/cran/rbmi/man/rubin_df.html) -- small-sample degrees of freedom
- [cards package documentation](https://insightsengineering.github.io/cards/main/) -- `as_card()`, `add_calculated_row()`
- [gt font size issues](https://github.com/rstudio/gt/issues/337) -- `table.font.size` does not affect column labels
- [gt rendering HTML vs PDF](https://forum.posit.co/t/rmarkdown-gt-package-rendering-different-in-html-and-pdf/140864) -- cross-format differences
- [patchwork layout guide](https://patchwork.data-imaginist.com/articles/guides/layout.html) -- alignment, `free()`, panel proportions
- [patchwork title alignment issue](https://github.com/thomasp85/patchwork/issues/276) -- text angle alignment bugs
- [pkgdown image path issues](https://github.com/r-lib/pkgdown/issues/1472) -- `man/figures/` relocation to `reference/figures/`
- [pkgdown figure path relocation](https://github.com/r-lib/pkgdown/issues/1943) -- path rewriting behavior
- [pkgdown false missing image warnings](https://github.com/r-lib/pkgdown/issues/2194) -- build order: `build_home_index()` before `build_articles()`
- [pkgdown markdown image resizing](https://github.com/r-lib/pkgdown/issues/2389) -- resizing syntax breaks on pkgdown
- [R Packages (2e) Dependencies](https://r-pkgs.org/dependencies-in-practice.html) -- Suggests/Imports guidance
- [R Packages (2e) Other Markdown](https://r-pkgs.org/other-markdown.html) -- README/pkgdown interaction
- [testthat snapshot testing](https://testthat.r-lib.org/articles/snapshotting.html) -- platform-specific variants
- Codebase analysis: `R/ard_conversion.R`, `R/tidiers.R`, `R/pool_methods.R`, `R/efficacy_table.R`, `R/plot_forest.R`, `R/analyse_mi_data.R`, `R/utils.R`, `tests/testthat/test-ard_conversion.R`, `_pkgdown.yml`, `man/figures/`

---

*Pitfalls research for v3 ARD Enrichment & Polish: 2026-02-10*
