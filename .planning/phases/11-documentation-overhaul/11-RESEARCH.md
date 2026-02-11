# Phase 11: Documentation Overhaul - Research

**Researched:** 2026-02-11
**Domain:** R package documentation (roxygen2, pkgdown, vignettes, NEWS.md, pre-rendered images)
**Confidence:** HIGH

## Summary

Phase 11 is a documentation overhaul that brings all package documentation in line with the finalized v3 features delivered in phases 8, 9, and 10. The three new feature areas that need documenting are: (1) MI diagnostic statistics via `pool_to_ard(pool_obj, analysis_obj)` enrichment (phase 8), (2) `describe_draws()` and `describe_imputation()` helper functions (phase 9), and (3) publication styling parameters for `efficacy_table()` (`font_family`, `font_size`, `row_padding`) and `plot_forest()` (`font_family`, `panel_widths`) (phase 10).

The package already has a solid documentation foundation: four vignettes (pipeline, analyse2, data-preparation, efficient-storage), a README.Rmd that generates README.md with pre-rendered images, roxygen2 function documentation with `@examples` blocks, a `_pkgdown.yml` configuration, and a `data-raw/generate-doc-images.R` script for regenerating man/figures images. The existing NEWS.md covers v0.1.0 and v0.2.0 but has no v0.3.0 entries.

The primary work involves five streams: (a) updating README.Rmd to show a fuller ADEFF-through-pipeline workflow, (b) improving function `@examples` to use ADMI/ADEFF data with realistic patterns instead of commented-out pseudocode, (c) creating or updating vignettes for MI diagnostics and describe helpers, (d) regenerating images via `generate-doc-images.R` to reflect v3 styling, and (e) writing NEWS.md v0.3.0 entries.

**Primary recommendation:** Structure as 5 plan files matching DOCS-01 through DOCS-05, each independently executable, with the image regeneration (DOCS-04) last since it depends on any README/help-page styling changes from earlier plans.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| roxygen2 | 7.3.2 | Function documentation from inline R comments | Already in use; generates .Rd files from #' comments |
| pkgdown | 2.x | Package website generation | Already configured via `_pkgdown.yml` |
| knitr | current | Vignette rendering and README.Rmd knitting | Already VignetteBuilder in DESCRIPTION |
| rmarkdown | current | Vignette output format | Already in Suggests |
| gt | current | Table rendering for help page images | Already in Suggests; used by `efficacy_table()` |
| ggplot2 | current | Plot rendering for help page images | Already in Suggests; used by `plot_forest()` |
| patchwork | current | Forest plot panel composition | Already in Suggests; used by `plot_forest()` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| webshot2/chromote | current | gt table PNG export via `gt::gtsave()` | Image regeneration only (dev dependency) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Pre-rendered README images | Computed-on-render README.Rmd | Pre-rendered is better: avoids rbmi MCMC computation during knit, faster CI |
| New vignette for diagnostics | Extend existing pipeline vignette | Separate vignette is cleaner: pipeline.Rmd is already 390 lines, diagnostics deserve focused coverage |

## Architecture Patterns

### Current Documentation Structure
```
rbmiUtils/
  README.Rmd                      # Source for README.md, refs pre-rendered images
  README.md                       # Generated from README.Rmd
  NEWS.md                         # Changelog (currently v0.1.0 and v0.2.0)
  DESCRIPTION                     # Version 0.1.9 (needs bump to 0.3.0)
  _pkgdown.yml                    # pkgdown site configuration
  R/
    *.R                           # Source with roxygen2 @examples
  man/
    *.Rd                          # Generated from roxygen2
    figures/
      README-forest-plot-1.png    # Pre-rendered for README
      README-efficacy-table-1.png # Pre-rendered for README
      plot_forest-trt.png         # Pre-rendered for help page
      efficacy_table-example.png  # Pre-rendered for help page
      logo.png                    # Package logo
  vignettes/
    pipeline.Rmd                  # End-to-end walkthrough
    analyse2.Rmd                  # Storing and analyzing imputed data
    data-preparation.Rmd          # Pre-imputation validation
    efficient-storage.Rmd         # Reduce/expand imputed data
  data-raw/
    generate-doc-images.R         # Script to regenerate all man/figures PNGs
```

### Pattern 1: Pre-Rendered Images for README and Help Pages
**What:** Images referenced in README.md and roxygen `\figure{}` tags are pre-generated PNGs stored in `man/figures/`, produced by `data-raw/generate-doc-images.R`.
**When to use:** Any time the visual output of `efficacy_table()` or `plot_forest()` changes (e.g., after styling parameter additions in phase 10).
**Key details:**
- README.Rmd sets `fig.path = "man/figures/README-"` but uses static image references (not computed chunks) for the teaser images
- `generate-doc-images.R` runs the full analysis pipeline (ADMI -> analyse_mi_data -> pool -> outputs) and saves PNGs
- gt table PNG export requires `gt::gtsave()` which needs chromote/webshot2
- ggplot2 plot export uses `ggsave()`

### Pattern 2: Roxygen Examples with Package Datasets
**What:** Function `@examples` blocks use the bundled `ADMI` and `ADEFF` datasets to show real usage.
**When to use:** All exported functions should demonstrate usage with package data.
**Key considerations:**
- `ADMI` has 100,000 rows (100 imputations x 1,000 subject-visits) -- analysis pipeline examples need `\donttest{}` due to CRAN time limits
- `ADEFF` has 1,000 rows -- suitable for lightweight examples
- `describe_draws()` and `describe_imputation()` currently use `\dontrun{}` because they require rbmi objects that cannot be created without MCMC -- consider switching to `\donttest{}` with actual pipeline execution
- `pool_to_ard()` examples are commented-out pseudocode inside `\donttest{}` -- should be upgraded to executable code using ADMI

### Pattern 3: Vignette Cross-References
**What:** Vignettes link to each other via relative HTML paths and to function docs via `[function_name()]`.
**When to use:** Any new vignette should be cross-referenced from existing vignettes and from `_pkgdown.yml` navbar.

### Pattern 4: NEWS.md Section Structure
**What:** NEWS.md uses top-level headings per version with categorized entries.
**Current structure:**
```markdown
# rbmiUtils 0.2.0
## New Features
* Item
## Improvements
* Item

# rbmiUtils 0.1.0
## Breaking Changes
## New Features
## Dependencies
## Improvements
## Previous Releases
```

### Anti-Patterns to Avoid
- **Commented-out example code:** Examples like `# pool_obj <- rbmi::pool(analysis_obj)` provide no value -- they cannot be executed or tested. Use `\donttest{}` with actual executable code instead.
- **Synthetic toy data in examples when package data exists:** The `data-preparation.Rmd` vignette creates a 20-subject synthetic dataset. For DOCS-02, function examples should use `ADMI`/`ADEFF` which are realistic clinical trial data already bundled with the package.
- **Missing pkgdown reference entries:** `describe_draws()` and `describe_imputation()` are exported but not listed in `_pkgdown.yml` reference sections -- they will appear in an "Other" catch-all section instead of being properly categorized.
- **Stale images after styling changes:** The images in `man/figures/` were last regenerated 2026-02-10 (before phase 10 styling changes like `font_family` and `panel_widths`). If the default styling changed, images should be re-rendered.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| README image generation | Manual screenshot workflow | `data-raw/generate-doc-images.R` script | Already exists, ensures reproducibility |
| gt table to PNG | Custom HTML-to-image pipeline | `gt::gtsave(..., vwidth = 800)` | Handles Chrome rendering automatically |
| ggplot to PNG | Custom rendering | `ggplot2::ggsave()` | Standard, handles DPI and dimensions |
| ARD structure docs | Manual .Rd editing | roxygen2 `@details` and `@return` tags | Auto-generates .Rd from source comments |
| pkgdown site structure | Manual HTML | `_pkgdown.yml` reference sections | Declarative, maintained automatically |
| NEWS.md formatting | Custom changelog format | Standard R `# pkg version` heading format | Recognized by pkgdown, CRAN, usethis |

**Key insight:** All documentation tooling (roxygen2, pkgdown, knitr) is already configured. The task is updating content, not setting up infrastructure.

## Common Pitfalls

### Pitfall 1: README.Rmd / README.md Sync
**What goes wrong:** Editing README.md directly instead of README.Rmd, then losing changes on next knit.
**Why it happens:** README.md is the rendered output; edits belong in README.Rmd.
**How to avoid:** Always edit README.Rmd, then knit to regenerate README.md. The file includes a comment: "README.md is generated from README.Rmd. Please edit that file."
**Warning signs:** Git diff shows changes in README.md but not README.Rmd.

### Pitfall 2: Examples That Exceed CRAN Time Limits
**What goes wrong:** Executable examples that run rbmi MCMC pipelines take too long for R CMD check.
**Why it happens:** `draws()` with MCMC sampling is computationally expensive even with small n_samples.
**How to avoid:** Wrap MCMC-dependent examples in `\donttest{}`. Use `ADMI` dataset (pre-imputed) to skip the draws/impute steps. For `describe_draws()` and `describe_imputation()` which require rbmi objects, `\dontrun{}` is acceptable since these objects cannot be created cheaply.
**Warning signs:** R CMD check times out or takes > 10 minutes.

### Pitfall 3: Forgetting to Update _pkgdown.yml
**What goes wrong:** New exported functions appear in an "Other" section on the pkgdown site instead of being properly categorized.
**Why it happens:** `_pkgdown.yml` reference sections are manually maintained.
**How to avoid:** After adding new exports, update `_pkgdown.yml` reference sections. Currently missing: `describe_draws`, `describe_imputation`, `print.describe_draws`, `print.describe_imputation`.
**Warning signs:** Functions listed under "Other" on the pkgdown site.

### Pitfall 4: Stale NEWS.md Without Version Bump
**What goes wrong:** NEWS.md has v0.3.0 entries but DESCRIPTION still says Version: 0.1.9.
**Why it happens:** NEWS.md is updated but DESCRIPTION version is not bumped to match.
**How to avoid:** Coordinate NEWS.md entries with DESCRIPTION Version field. Both should reference v0.3.0.
**Warning signs:** Mismatch between DESCRIPTION version and NEWS.md top entry.

### Pitfall 5: Vignette Dependencies Not in Suggests
**What goes wrong:** New vignette code fails on CRAN check machines that lack optional packages.
**Why it happens:** Vignette uses a package not listed in Suggests.
**How to avoid:** Any package used in a vignette must be in Suggests (or Imports). Guard optional-package vignette chunks with `eval = requireNamespace("pkg", quietly = TRUE)`.
**Warning signs:** R CMD check fails on machines without optional packages.

### Pitfall 6: Image Regeneration Requires Chromium
**What goes wrong:** `gt::gtsave()` fails because no Chromium browser is available.
**Why it happens:** gt uses chromote for HTML-to-PNG rendering.
**How to avoid:** Document the chromote requirement in `generate-doc-images.R` (already done). Ensure the regeneration step mentions this prerequisite.
**Warning signs:** Error about missing Chrome/Chromium when running `generate-doc-images.R`.

## Code Examples

### README.Rmd Workflow Pattern (DOCS-01)
The README Quick Start currently shows a minimal ADMI-based example. The updated version should show the full ADEFF-through-pipeline workflow:

```r
# Current (minimal):
ana_obj  <- analyse_mi_data(data = ADMI, vars = vars, method = method)
pool_obj <- pool(ana_obj)
efficacy_table(pool_obj, arm_labels = c(ref = "Placebo", alt = "Drug A"))

# Target (ADEFF-through-pipeline):
data("ADEFF", package = "rbmiUtils")
# ... factor prep, vars definition, draws(), impute(), get_imputed_data() ...
# ... analyse_mi_data(), pool() ...
# ... efficacy_table() and plot_forest() with new styling params ...
```

### Executable Example Pattern for pool_to_ard (DOCS-02)
```r
#' @examples
#' \donttest{
#' library(rbmi)
#' data("ADMI", package = "rbmiUtils")
#' ADMI$TRT <- factor(ADMI$TRT, levels = c("Placebo", "Drug A"))
#' ADMI$USUBJID <- factor(ADMI$USUBJID)
#' ADMI$AVISIT <- factor(ADMI$AVISIT)
#'
#' vars <- set_vars(
#'   subjid = "USUBJID", visit = "AVISIT", group = "TRT",
#'   outcome = "CHG", covariates = c("BASE", "STRATA", "REGION")
#' )
#' method <- method_bayes(n_samples = 20,
#'   control = control_bayes(warmup = 20, thin = 1))
#'
#' ana_obj <- analyse_mi_data(ADMI, vars, method, fun = ancova)
#' pool_obj <- pool(ana_obj)
#'
#' # Base ARD
#' if (requireNamespace("cards", quietly = TRUE)) {
#'   ard <- pool_to_ard(pool_obj)
#'
#'   # Enriched ARD with MI diagnostics
#'   ard_diag <- pool_to_ard(pool_obj, analysis_obj = ana_obj)
#'   # Filter diagnostic stats
#'   ard_diag[ard_diag$stat_name == "fmi", ]
#' }
#' }
```

### describe_draws/describe_imputation Example Pattern (DOCS-02)
These functions require rbmi draws/impute objects which need MCMC -- `\dontrun{}` is appropriate, but examples should be realistic:
```r
#' @examples
#' \dontrun{
#' library(rbmi)
#' data("ADEFF", package = "rbmiUtils")
#' ADEFF$TRT <- factor(ADEFF$TRT01P, levels = c("Placebo", "Drug A"))
#' ADEFF$USUBJID <- factor(ADEFF$USUBJID)
#' ADEFF$AVISIT <- factor(ADEFF$AVISIT, levels = c("Week 24", "Week 48"))
#'
#' vars <- set_vars(
#'   subjid = "USUBJID", visit = "AVISIT", group = "TRT",
#'   outcome = "CHG", covariates = c("BASE", "STRATA", "REGION")
#' )
#' dat <- ADEFF[, c("USUBJID","STRATA","REGION","TRT","BASE","CHG","AVISIT")]
#' draws_obj <- draws(data = dat, vars = vars,
#'   method = method_bayes(n_samples = 100))
#'
#' # Inspect the draws object
#' desc <- describe_draws(draws_obj)
#' print(desc)
#' desc$method
#' desc$n_samples
#' }
```

### NEWS.md v0.3.0 Pattern (DOCS-05)
```markdown
# rbmiUtils 0.3.0

## New Features

* `describe_draws()` extracts structured metadata from rbmi draws objects,
  including method type, formula, sample count, and (for Bayesian methods)
  MCMC convergence diagnostics.
* `describe_imputation()` extracts imputation metadata including method,
  M, reference arm mappings, and a missingness breakdown by visit and arm.
* `pool_to_ard()` gains an `analysis_obj` parameter that enriches the ARD
  with MI diagnostic statistics (FMI, lambda, RIV, Barnard-Rubin adjusted
  df, relative efficiency) when the pooling method is Rubin's rules.

## Improvements

* `efficacy_table()` gains `font_family`, `font_size`, and `row_padding`
  parameters for publication-ready table styling.
* `plot_forest()` gains `font_family` and `panel_widths` parameters for
  customizable typography and panel layout.
* `plot_forest()` left panel now uses left-aligned text (hjust=0) for
  consistent positioning.
* Non-Rubin pooling methods now emit an informative `cli::cli_inform()`
  message when MI diagnostics are not applicable, rather than returning
  NA rows.
```

### Updated generate-doc-images.R Pattern (DOCS-04)
The existing script should be updated to use v3 styling parameters:
```r
# Updated forest plot with font_family demo
p_forest <- plot_forest(
  pool_obj,
  title = "Treatment Effect: Change from Baseline",
  arm_labels = c(ref = "Placebo", alt = "Drug A"),
  text_size = 4.5,
  point_size = 4,
  show_pvalues = FALSE
  # font_family and panel_widths use defaults for doc images
)

# Updated efficacy table with font_size and row_padding
tbl <- efficacy_table(
  pool_obj,
  title = "Table 14.2.1: ANCOVA of Change from Baseline",
  subtitle = "Mixed Model for Repeated Measures",
  arm_labels = c(ref = "Placebo", alt = "Drug A")
  # font_family, font_size, row_padding use defaults for doc images
)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `pool_to_ard(pool_obj)` only | `pool_to_ard(pool_obj, analysis_obj)` with MI diagnostics | Phase 8 (v0.3.0) | ARD now includes FMI/lambda/RIV/df stats |
| No pipeline introspection | `describe_draws()` and `describe_imputation()` | Phase 9 (v0.3.0) | Users can inspect rbmi objects before analysis |
| No typography control | `font_family`/`font_size`/`row_padding`/`panel_widths` | Phase 10 (v0.3.0) | Publication-ready output customization |
| Commented-out examples | Executable `\donttest{}` examples with ADMI/ADEFF | Phase 11 target | Examples testable and demonstrative |

**Deprecated/outdated:**
- None of the existing API has been deprecated. All v3 additions are backward-compatible (new optional parameters with NULL defaults).

## Inventory of Documentation Gaps

### Gap 1: README Missing Full ADEFF Pipeline
**Current:** README Quick Start shows 4 lines (analyse_mi_data -> pool -> outputs) using ADMI.
**Target:** Show realistic workflow starting from ADEFF data through factor prep, vars, draws, impute, get_imputed_data, analyse, pool, to efficacy_table and plot_forest.
**Requirement:** DOCS-01

### Gap 2: v3 Functions Not in README Key Features
**Current:** Key Features list does not mention `describe_draws()`, `describe_imputation()`, or enriched `pool_to_ard()`.
**Target:** Add these to the Key Features bullet list.
**Requirement:** DOCS-01

### Gap 3: Commented-Out Function Examples
**Current:** `pool_to_ard()`, `efficacy_table()`, and `plot_forest()` have commented-out pseudocode in `@examples`.
**Target:** Executable `\donttest{}` examples using ADMI/ADEFF.
**Requirement:** DOCS-02

### Gap 4: describe_* Examples Use Generic \dontrun
**Current:** `describe_draws()` and `describe_imputation()` examples use minimal `\dontrun{}` code.
**Target:** Realistic `\dontrun{}` examples using ADEFF data (these must remain `\dontrun` due to MCMC requirement).
**Requirement:** DOCS-02

### Gap 5: No Vignette Coverage for MI Diagnostics
**Current:** No vignette demonstrates `pool_to_ard(pool_obj, analysis_obj)` enrichment or describes the diagnostic stats.
**Target:** New vignette or extension of pipeline vignette showing MI diagnostics workflow.
**Requirement:** DOCS-03

### Gap 6: No Vignette Coverage for describe_* Helpers
**Current:** `describe_draws()` and `describe_imputation()` have no vignette coverage.
**Target:** Vignette section showing how to use these functions to inspect rbmi pipeline objects.
**Requirement:** DOCS-03

### Gap 7: _pkgdown.yml Missing New Functions
**Current:** `describe_draws`, `describe_imputation` not in `_pkgdown.yml` reference sections.
**Target:** Add to appropriate reference section (new "Introspection" or "Describe" category).
**Requirement:** DOCS-03 (or DOCS-01)

### Gap 8: Images May Not Reflect v3 Styling
**Current:** Images last regenerated 2026-02-10 (phase 10 styling changes were on 2026-02-11).
**Target:** Regenerate all images via updated `generate-doc-images.R`.
**Requirement:** DOCS-04

### Gap 9: No NEWS.md v0.3.0 Entries
**Current:** NEWS.md only has v0.1.0 and v0.2.0.
**Target:** v0.3.0 section documenting all phase 8/9/10 additions.
**Requirement:** DOCS-05

### Gap 10: DESCRIPTION Version Not Bumped
**Current:** Version: 0.1.9 in DESCRIPTION.
**Target:** Version: 0.3.0 to match NEWS.md entries.
**Requirement:** DOCS-05

### Gap 11: Learn More Section Missing Diagnostics Vignette
**Current:** README "Learn More" section links to pipeline and analyse2 vignettes only.
**Target:** Add link to new diagnostics/describe vignette.
**Requirement:** DOCS-01

## Recommended Plan Structure

| Plan | Requirement | Scope | Dependencies |
|------|-------------|-------|--------------|
| 11-01 | DOCS-01 | README.Rmd overhaul: ADEFF pipeline workflow, v3 features in Key Features, Learn More links, _pkgdown.yml updates | None |
| 11-02 | DOCS-02 | Function @examples: executable examples for pool_to_ard, efficacy_table, plot_forest; realistic ADEFF examples for describe_draws, describe_imputation | None |
| 11-03 | DOCS-03 | Vignettes: new MI diagnostics + describe helpers vignette | None |
| 11-04 | DOCS-04 | Image regeneration: update generate-doc-images.R, regenerate all man/figures PNGs | After 11-01 (README may change image requirements) |
| 11-05 | DOCS-05 | NEWS.md v0.3.0 entries + DESCRIPTION version bump | After 11-01 through 11-03 (needs complete feature inventory) |

## Open Questions

1. **DESCRIPTION version bump scope**
   - What we know: Current version is 0.1.9, NEWS.md targets v0.3.0
   - What's unclear: Whether to bump to 0.3.0 in this phase or defer to a release phase
   - Recommendation: Include version bump in DOCS-05 since NEWS.md entries reference v0.3.0 and they should be consistent

2. **README Quick Start length**
   - What we know: DOCS-01 requires showing ADEFF-through-pipeline workflow including draws/impute steps
   - What's unclear: How much pipeline detail to include in README vs. deferring to vignette
   - Recommendation: Show a streamlined but complete workflow (eval=FALSE code blocks) with link to pipeline vignette for the full walkthrough. The README should be a teaser, not the tutorial.

3. **New vignette vs. extending existing**
   - What we know: DOCS-03 requires covering MI diagnostics and describe helpers
   - What's unclear: Whether to create a new vignette or extend pipeline.Rmd
   - Recommendation: Create a new focused vignette (e.g., `diagnostics.Rmd`) since pipeline.Rmd is already 390 lines and the diagnostics topic is distinct enough to warrant its own article. However, this vignette will need rbmi MCMC objects, so chunks may need `eval=FALSE` or use pre-computed objects.

4. **Vignette eval strategy for MCMC-dependent code**
   - What we know: `describe_draws()` and `describe_imputation()` require rbmi draws/impute objects that need MCMC
   - What's unclear: Whether to use `eval=FALSE` with printed output, or pre-compute and load saved objects
   - Recommendation: Use the same pattern as pipeline.Rmd (run MCMC inline with small n_samples). The pipeline vignette already does this successfully. Alternatively, use `eval=FALSE` blocks with hand-formatted output comments for the describe_* sections only.

## Sources

### Primary (HIGH confidence)
- Codebase inspection: All R source files, vignettes, README.Rmd, _pkgdown.yml, generate-doc-images.R, NEWS.md, DESCRIPTION, NAMESPACE
- Phase 8/9/10 SUMMARY.md files documenting implemented features
- Phase 8 CONTEXT.md documenting locked decisions

### Secondary (MEDIUM confidence)
- [R Packages (2e) - Function documentation](https://r-pkgs.org/man.html) - roxygen2 best practices
- [R Packages (2e) - Vignettes](https://r-pkgs.org/vignettes.html) - vignette patterns
- [R-hub blog - Code examples](https://blog.r-hub.io/2020/01/27/examples/) - \donttest vs \dontrun guidance
- [R-hub blog - NEWS.md](https://blog.r-hub.io/2020/05/08/pkg-news/) - changelog formatting
- [rOpenSci - Precompute vignettes](https://ropensci.org/blog/2019/12/08/precompute-vignettes/) - pre-rendered vignette patterns

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All tools already in use in the package, verified by codebase inspection
- Architecture: HIGH - Documentation patterns well-established in the existing package
- Pitfalls: HIGH - Based on direct codebase analysis and known R package documentation patterns
- Content gaps: HIGH - Complete inventory derived from comparing current docs against phase 8/9/10 deliverables

**Research date:** 2026-02-11
**Valid until:** 2026-03-11 (stable domain -- documentation tooling and patterns change slowly)
