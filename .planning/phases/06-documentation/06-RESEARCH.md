# Phase 6: Documentation - Research

**Researched:** 2026-02-08
**Domain:** R package documentation (vignettes, README.Rmd, roxygen2, NEWS.md)
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Use package datasets (ADEFF/ADMI) -- no external dependencies, always reproducible
- Show full rbmi pipeline (draws/impute/analyse/pool) with explanation so reader understands the complete workflow
- Continuous analysis as primary walkthrough, binary/responder as a shorter appendix section showing how it differs
- Pipeline-focused framing (e.g., "From rbmi Analysis to Regulatory Tables") rather than rbmiUtils-centric title
- Keep both the existing analyse2 vignette and the new end-to-end vignette -- link between them
- README.Rmd renders to README.md (standard R package pattern, executable code examples)
- Rendered table and plot images appear in README as visual teasers
- Inline hyperlinks woven into prose (not callout boxes) in vignettes
- Add cross-references to existing vignettes (analyse2) where rbmi/beeca links add value
- Add @seealso sections in roxygen for functions that wrap or depend on rbmi/beeca
- Organize NEWS.md by version number: v1 milestone = 0.1.0, v2 milestone = 0.2.0
- Group entries with sub-bullets (New features, Bug fixes, Improvements, Breaking changes)
- NEWS.md already exists -- reorganize existing content into proper versioned structure

### Claude's Discretion
- Writing tone for vignette (tutorial vs reference -- pick what fits the audience)
- Whether to include data prep steps (validate_data, prepare_data_ice) in the pipeline narrative
- Customization depth in vignette (defaults only vs showing key customizations)
- README depth (minimal teaser vs quick-start snippet)
- Rendered output approach for function help pages (static images vs live rendering)
- Cross-reference link targets (pkgdown sites vs CRAN pages)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 6 is an R package documentation phase covering six requirements: an end-to-end vignette (DOC-01), enhanced README with visual teasers (DOC-02), rendered examples for `plot_forest()` and `efficacy_table()` (DOC-03, DOC-04), NEWS.md restructuring (DOC-05), and cross-references to rbmi/beeca docs in existing vignettes (DOC-06). The core technologies are standard R package documentation tooling: knitr/rmarkdown for vignettes, roxygen2 for function docs, and devtools::build_readme() for README rendering.

The main technical challenge is that the rbmi analysis pipeline (draws/impute/analyse/pool) is computationally expensive (MCMC sampling), making standard vignette building during `R CMD check` very slow. The existing analyse2 vignette handles this by running the code live (it uses the pre-computed ADMI dataset to skip the draws/impute step). The new end-to-end vignette needs to run the full pipeline from ADEFF data, which means either accepting the build time or using a precompute strategy. Since ADEFF has 1000 rows across 2 visits and the vignette uses `n_samples = 100`, this should complete in under 60 seconds on modern hardware -- acceptable for a vignette that builds only during `devtools::build_vignettes()` and `pkgdown::build_site()`, not during `R CMD check` (which uses the pre-built vignette from submission).

For rendered examples in help pages, the two viable approaches are: (1) static pre-rendered PNG images included via roxygen2 markdown `![](figure.png)` syntax in the `@details` or `@description` section, or (2) executable examples wrapped in `\donttest{}` that render live when the user runs `example(plot_forest)`. The static image approach is more reliable for display but adds maintenance burden. The recommendation is to use static images in `man/figures/` referenced via markdown in the roxygen documentation, because `efficacy_table()` produces a gt table (not a plot) which cannot be rendered in a standard R help page, and `plot_forest()` output requires ggplot2+patchwork which are Suggests dependencies.

**Primary recommendation:** Use a tutorial-style vignette tone with the precomputed ADEFF pipeline, include data prep steps briefly, show key customizations, and use static PNG images in man/figures for function help pages.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| knitr | >= 1.40 | Vignette engine | Already in Suggests + VignetteBuilder |
| rmarkdown | >= 2.20 | Vignette rendering | Already in Suggests, standard html_vignette output |
| roxygen2 | >= 7.3.0 | Function documentation with markdown | Already configured in DESCRIPTION (RoxygenNote: 7.3.2) |
| devtools | (dev tool) | `build_readme()`, `build_vignettes()` | Standard R package development workflow |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| gt | (Suggests) | Render efficacy_table() output | Already in Suggests, needed for table screenshots |
| ggplot2 | (Suggests) | Render plot_forest() output | Already in Suggests, needed for plot screenshots |
| patchwork | (Suggests) | Multi-panel forest plot layout | Already in Suggests, required by plot_forest() |
| webshot2 | (dev tool) | Screenshot gt tables to PNG | Only needed once to generate static images |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Static images in man/figures | Live `\donttest{}` examples | Static images work everywhere; live examples require all Suggests deps and long runtime |
| html_vignette output | html_document output | html_vignette is smaller (CRAN-friendly), already used in existing vignettes |
| Precompute vignette (.Rmd.orig) | Live-build vignette | Live build is simpler; precompute only needed if build time exceeds ~2 min |

## Architecture Patterns

### Documentation File Structure
```
rbmiUtils/
+-- README.Rmd              # Source (renders to README.md)
+-- README.md               # Generated (committed to git)
+-- NEWS.md                 # Version history
+-- man/
|   +-- figures/
|   |   +-- rbmiUtils.png             # Existing hex logo
|   |   +-- README-forest-plot-1.png  # Generated by README.Rmd
|   |   +-- README-efficacy-table-1.png # Generated by README.Rmd
|   |   +-- plot_forest-trt.png       # Static image for help page
|   |   +-- efficacy_table-example.png # Static image for help page
|   +-- plot_forest.Rd
|   +-- efficacy_table.Rd
+-- vignettes/
|   +-- .gitignore
|   +-- analyse2.Rmd                  # Existing (enhanced with cross-refs)
|   +-- data-preparation.Rmd          # Existing (enhanced with cross-refs)
|   +-- efficient-storage.Rmd         # Existing (enhanced with cross-refs)
|   +-- pipeline.Rmd                  # NEW end-to-end vignette
```

### Pattern 1: End-to-End Vignette Structure
**What:** A tutorial-style vignette walking through the complete rbmi + rbmiUtils pipeline
**When to use:** This is the primary getting-started guide for new users

The vignette should follow this narrative structure:
```r
# Typical vignette structure for pipeline.Rmd:

# 1. Introduction (what the reader will learn)
# 2. Setup & Data Loading (library calls, data("ADEFF"))
# 3. Data Preparation (validate_data, prepare_data_ice, summarise_missingness)
# 4. rbmi Pipeline (set_vars, method_bayes, draws, impute)
# 5. Analysis (analyse_mi_data with ancova)
# 6. Results (pool, tidy_pool_obj)
# 7. Reporting Output
#    a. Efficacy Table (efficacy_table())
#    b. Forest Plot (plot_forest())
# 8. Appendix: Binary/Responder Analysis
#    (brief section showing how the pipeline differs for binary endpoints)
```

### Pattern 2: README.Rmd with Visual Teasers
**What:** README.Rmd that renders executable code producing images saved to man/figures/
**When to use:** Standard R package README pattern

```r
# In README.Rmd YAML and setup:
# ---
# output: github_document
# ---
# ```{r, include = FALSE}
# knitr::opts_chunk$set(
#   collapse = TRUE,
#   comment = "#>",
#   fig.path = "man/figures/README-",
#   out.width = "100%"
# )
# ```

# For the forest plot teaser:
# ```{r forest-plot, fig.width=10, fig.height=5}
# p <- plot_forest(pool_obj, title = "Treatment Effect by Visit")
# p
# ```

# For the efficacy table teaser (save as image):
# ```{r efficacy-table, eval=FALSE}
# tbl <- efficacy_table(pool_obj,
#   title = "Table 14.2.1",
#   arm_labels = c(ref = "Placebo", alt = "Drug A")
# )
# ```
# Then include a pre-rendered image:
# ![Efficacy Table](man/figures/README-efficacy-table-1.png)
```

Note: gt tables cannot be rendered inline in github_document output. The README must either: (a) save a screenshot of the gt table using `gt::gtsave()` and include it as a static image, or (b) show the tidy data frame output instead. The recommended approach is to generate the gt table image via a separate script and include it as a static image.

### Pattern 3: Static Images in Roxygen Help Pages
**What:** Pre-rendered PNG images included in function documentation via markdown
**When to use:** For `plot_forest()` and `efficacy_table()` help pages

```r
# In R/plot_forest.R roxygen:
#' @details
#' ...existing details...
#'
#' **Example output (treatment difference mode):**
#'
#' \if{html}{\figure{plot_forest-trt.png}{options: width=80\%}}
#' \if{latex}{\figure{plot_forest-trt.pdf}{options: width=5in}}
```

Alternative markdown approach (simpler, works in HTML help):
```r
#' **Example output:**
#'
#' ![Treatment difference forest plot](plot_forest-trt.png)
```

The image files must be placed in `man/figures/` directory.

### Pattern 4: NEWS.md Tidyverse Convention
**What:** Version-organized changelog following tidyverse style guide
**When to use:** Standard NEWS.md structure

```markdown
# rbmiUtils 0.2.0

## New Features

* `efficacy_table()` creates regulatory-style gt tables from pool objects.
* `plot_forest()` creates three-panel forest plots from pool objects.
* `pool_to_ard()` converts pool objects to pharmaverse ARD format.

## Improvements

* `tidy_pool_obj()` now uses regex-based parameter parsing (#HRD-01).

## Breaking Changes

* `tidy_pool_obj()` output columns changed from truncated to full visit names.

# rbmiUtils 0.1.0

## New Features

* Initial release with core analysis and data preparation utilities.
* ...
```

### Pattern 5: Cross-References in Vignette Prose
**What:** Inline hyperlinks to rbmi and beeca documentation woven into narrative
**When to use:** In all vignettes where rbmi/beeca concepts are mentioned

```r
# Inline link pattern in vignette Rmd:
# "The [rbmi](https://cran.r-project.org/package=rbmi) package implements
# reference-based multiple imputation. The core pipeline consists of
# [`draws()`](https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html),
# `impute()`, `analyse()`, and `pool()` -- see the
# [rbmi quickstart vignette](https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html)
# for details."

# For beeca references:
# "Binary responder analysis uses
# [`beeca::get_marginal_effect()`](https://openpharma.github.io/beeca/reference/get_marginal_effect.html)
# for covariate-adjusted marginal effects."
```

### Anti-Patterns to Avoid
- **Running draws() in README.Rmd:** The MCMC sampling is too slow for README rendering. Use the pre-computed ADMI dataset in README examples instead.
- **Using `\dontrun{}` for working examples:** CRAN prefers `\donttest{}` for examples that work but are slow. `\dontrun{}` is for examples that literally cannot run without special setup.
- **Hardcoding pkgdown URLs that may change:** Use CRAN URLs as the stable link target (e.g., `https://cran.r-project.org/package=rbmi`) rather than GitHub Pages URLs which could move.
- **Evaluating gt table rendering in github_document:** gt tables produce HTML; `github_document` output cannot display them inline. Use static images instead.
- **Forgetting to commit man/figures/ images:** Generated plot PNGs must be committed to git or they will not appear on GitHub/CRAN.

## Discretion Recommendations

Based on research, here are recommendations for the areas marked as Claude's Discretion:

### Writing Tone: Tutorial
**Recommendation:** Use tutorial tone, not reference tone. The audience is clinical trial statisticians who know rbmi but are new to rbmiUtils. A tutorial walks them through step-by-step, explaining "why" at each step. The existing analyse2 vignette already serves as a focused reference; the new vignette should be the approachable "getting started" guide.

### Data Prep Steps: Include Briefly
**Recommendation:** Include `validate_data()` and `summarise_missingness()` as a brief "Data Preparation" section (3-5 paragraphs). Skip `prepare_data_ice()` in the main narrative since the ADEFF dataset's default imputation does not use ICE flags in the simplest pipeline. Mention it in a "See Also" note linking to the data-preparation vignette.

### Customization Depth: Defaults Plus Key Customizations
**Recommendation:** Show defaults first (the simple path), then show 1-2 key customizations for each output function. For `efficacy_table()`: show `arm_labels` and `title`/`subtitle`. For `plot_forest()`: show `display = "lsm"` as an alternative to the default `"trt"` mode.

### README Depth: Quick-Start Snippet with Visual Teasers
**Recommendation:** Keep the README focused on installation + a compact pipeline snippet (using pre-computed ADMI to avoid draws() runtime) + rendered forest plot and efficacy table images as visual teasers + link to the full end-to-end vignette. The current README is already quite long; the enhancement should replace the verbose utility listing with visual output and a "Learn more" link.

### Rendered Output for Help Pages: Static Images
**Recommendation:** Use pre-rendered static PNG images stored in `man/figures/`. Reasons:
1. `efficacy_table()` returns a gt object which renders as HTML -- it cannot display in base R `?` help without an image
2. `plot_forest()` requires ggplot2 + patchwork (Suggests, not Imports) so `\donttest{}` examples would fail if those packages are not installed
3. Static images provide consistent visual documentation regardless of user environment
4. Images display well in both pkgdown site and base R HTML help

Generate images via a maintenance script (e.g., `data-raw/generate-doc-images.R`) that creates the pool object from ADMI, renders the outputs, and saves to `man/figures/`.

### Cross-Reference Link Targets: CRAN for Stability, pkgdown for Rich Content
**Recommendation:** Use CRAN URLs for package-level references (stable, permanent). Use pkgdown site URLs for specific function references and vignette deep-links (richer content). Both rbmi and beeca have pkgdown sites:
- rbmi: `https://cran.r-project.org/package=rbmi` (CRAN) / vignettes on CRAN
- beeca: `https://openpharma.github.io/beeca/` (pkgdown)
- rbmiUtils: `https://openpharma.github.io/rbmiUtils/` (pkgdown)

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| gt table screenshots | Manual HTML-to-PNG pipeline | `gt::gtsave("file.png")` | Handles CSS, fonts, cropping automatically; requires webshot2 |
| README rendering | `knitr::knit()` manually | `devtools::build_readme()` | Ensures rendering against current package source |
| Vignette cross-links | Raw HTML links | pkgdown auto-linking with backtick syntax | pkgdown resolves `function()` references to correct URLs automatically |
| NEWS.md parsing | Custom changelog | tidyverse NEWS.md conventions | pkgdown's `build_news()` parses standard format into HTML changelog |
| Roxygen figure inclusion | Custom Rd tags | `\if{html}{\figure{file.png}{options: width=80\%}}` | Standard Rd figure tag, works in R help viewer and pkgdown |

**Key insight:** R package documentation tooling (roxygen2, pkgdown, devtools) has mature, well-tested patterns for every documentation task in this phase. The risk is not in finding the right tool but in misusing the existing tools (e.g., wrong image paths, forgetting to rebuild, inconsistent cross-reference patterns).

## Common Pitfalls

### Pitfall 1: README.Rmd and README.md Out of Sync
**What goes wrong:** README.Rmd is edited but README.md is not re-rendered; GitHub shows stale content.
**Why it happens:** Developer edits README.Rmd but forgets to run `devtools::build_readme()`.
**How to avoid:** The `usethis::use_readme_rmd()` pre-commit hook prevents commits when README.Rmd is newer than README.md. Verify this hook is installed. Always run `devtools::build_readme()` after editing README.Rmd.
**Warning signs:** `git status` shows changes to README.Rmd but not README.md.

### Pitfall 2: man/figures Images Not Committed
**What goes wrong:** Images render locally but are missing on GitHub/CRAN because they were not committed to git.
**Why it happens:** `.gitignore` may exclude PNG files, or developer forgets to `git add` new image files in `man/figures/`.
**How to avoid:** After generating images, explicitly `git add man/figures/*.png`. Verify the `.gitignore` does not exclude `man/figures/`.
**Warning signs:** Broken image links on GitHub README preview.

### Pitfall 3: gt Table in github_document Output
**What goes wrong:** `efficacy_table()` output renders as raw HTML tags in the GitHub README instead of a formatted table.
**Why it happens:** github_document format does not support HTML widget output from gt. The gt package produces HTML which GitHub Markdown cannot render.
**How to avoid:** Save gt table as PNG image using `gt::gtsave()`, then include the image in README with `![](man/figures/...)`. Do not try to render gt inline in README.Rmd chunks.
**Warning signs:** Raw `<table>` HTML appearing in README.md output.

### Pitfall 4: Vignette Build Time Too Long for CRAN
**What goes wrong:** `R CMD check` takes too long because vignettes run expensive MCMC computations.
**Why it happens:** CRAN rebuilds vignettes during checks. If the vignette calls `rbmi::draws()` with many samples, it can exceed time limits.
**How to avoid:** CRAN distributes the pre-built vignettes from submission, but still evaluates vignette code during incoming checks. Use small sample sizes in vignette examples (n_samples = 100, warmup = 200) and monitor total build time. If build time exceeds 2 minutes, switch to the precompute pattern (.Rmd.orig + knitr::knit()).
**Warning signs:** `devtools::build_vignettes()` takes more than 2 minutes.

### Pitfall 5: Broken Cross-Reference Links
**What goes wrong:** Links to rbmi/beeca functions or vignettes return 404.
**Why it happens:** URLs change when packages move organizations or update pkgdown sites.
**How to avoid:** Prefer CRAN URLs for stability (`https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html`). These are maintained by CRAN infrastructure and are unlikely to break. Verify all external links before release.
**Warning signs:** Link rot detected during manual review or automated link checking.

### Pitfall 6: Forgetting VignetteIndexEntry in New Vignette
**What goes wrong:** New vignette does not appear in the package's vignette list.
**Why it happens:** Missing `%\VignetteIndexEntry{}` in the vignette YAML header.
**How to avoid:** Use `usethis::use_vignette("pipeline")` to create the skeleton, which includes all required metadata. Verify with `devtools::build_vignettes()`.
**Warning signs:** `browseVignettes("rbmiUtils")` does not show the new vignette.

### Pitfall 7: NEWS.md Format Not Parsed by pkgdown
**What goes wrong:** pkgdown changelog page does not render correctly -- missing version groupings or broken formatting.
**Why it happens:** NEWS.md uses inconsistent heading levels or non-standard formatting.
**How to avoid:** Follow tidyverse NEWS.md conventions exactly: `# package version` for level 1 headings, `## Section` for groupings. Test with `pkgdown::build_news()`.
**Warning signs:** pkgdown build warnings about NEWS.md parsing.

## Code Examples

### Example 1: Vignette YAML Header for New Pipeline Vignette
```yaml
---
title: "From rbmi Analysis to Regulatory Tables"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
    number_sections: true
vignette: >
  %\VignetteIndexEntry{From rbmi Analysis to Regulatory Tables}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```

### Example 2: README.Rmd Forest Plot Teaser
```r
# After pool_obj is created from ADMI data:
# ```{r forest-plot, fig.width=10, fig.height=4.5}
# library(ggplot2)
# library(patchwork)
# p <- plot_forest(pool_obj,
#   title = "Treatment Effect: Change from Baseline",
#   arm_labels = c(ref = "Placebo", alt = "Drug A")
# )
# p
# ```
```
Source: README.Rmd fig.path convention from R Packages (2e), Chapter 18

### Example 3: gt Table Screenshot for README
```r
# In a separate script (data-raw/generate-doc-images.R):
# library(rbmiUtils)
# library(gt)
#
# # Create pool_obj from ADMI (same as README example)
# # ... pipeline code ...
#
# tbl <- efficacy_table(pool_obj,
#   title = "Table 14.2.1: ANCOVA of Change from Baseline",
#   subtitle = "Mixed Model for Repeated Measures",
#   arm_labels = c(ref = "Placebo", alt = "Drug A")
# )
# gt::gtsave(tbl, "man/figures/README-efficacy-table-1.png", vwidth = 800)
```
Source: gt::gtsave() documentation

### Example 4: Static Image in Roxygen Documentation
```r
#' @details
#' The plot displays treatment effects with confidence intervals across visits.
#'
#' **Example output (treatment difference mode):**
#'
#' \if{html}{\figure{plot_forest-trt.png}{options: width=80\%}}
#' \if{latex}{\figure{plot_forest-trt.pdf}{options: width=5in}}
```
Source: roxygen2 rd-formatting vignette (https://roxygen2.r-lib.org/articles/rd-formatting.html)

### Example 5: Cross-Reference Pattern in Vignette Prose
```markdown
The analysis pipeline begins with the
[rbmi](https://cran.r-project.org/package=rbmi) package, which implements
reference-based multiple imputation for continuous longitudinal endpoints.
For statistical background, see the
[rbmi statistical specifications vignette](https://cran.r-project.org/web/packages/rbmi/vignettes/stat_specs.html).

For binary responder endpoints, `rbmiUtils` integrates with
[beeca](https://openpharma.github.io/beeca/) for covariate-adjusted
marginal effect estimation using the method of Ge et al. (2011).
```

### Example 6: @seealso Addition for Existing Functions
```r
#' @seealso
#' * [analyse_mi_data()] to analyse imputed datasets
#' * [rbmi::draws()] for the imputation model fitting step
#' * [rbmi::pool()] for pooling analysis results using Rubin's rules
#' * The [rbmi quickstart vignette](https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html)
#'   for the complete rbmi pipeline
```

### Example 7: NEWS.md Restructured Format
```markdown
# rbmiUtils 0.2.0

## Breaking Changes

* `tidy_pool_obj()` now uses regex-based parameter parsing instead of
  splitting on `_`. Output columns (`parameter_type`, `lsm_type`, `visit`)
  now contain the full visit name rather than truncated fragments.

## New Features

* `efficacy_table()` creates regulatory-style gt tables from pool objects.
* `plot_forest()` creates publication-quality three-panel forest plots.
* `pool_to_ard()` converts pool objects to pharmaverse ARD format.
* `print()` and `summary()` S3 methods for pool and analysis objects.
* `create_impid()` converts lists of imputed datasets into stacked
  data.frames.
* `combine_results()` combines tidy results from multiple analyses.
* `format_results()` for publication-ready formatting.
* `extract_trt_effects()` and `extract_lsm()` convenience filters.

## Improvements

* Added `cli` and `lifecycle` dependencies for improved messaging.
* Enhanced input validation in `analyse_mi_data()` and
  `prepare_data_ice()`.
* Expanded test coverage and integration tests.

# rbmiUtils 0.1.0

## New Features

* Initial release with core utilities.
* `analyse_mi_data()` applies analysis functions to imputed datasets.
* `tidy_pool_obj()` tidies pooled results for reporting.
* `get_imputed_data()` extracts imputed data from rbmi objects.
* `gcomp_responder()` and `gcomp_responder_multi()` for binary endpoints.
* `gcomp_binary()` for g-computation on binary outcomes.
* `validate_data()` for pre-flight validation before `rbmi::draws()`.
* `prepare_data_ice()` builds `data_ice` from flag columns.
* `summarise_missingness()` tabulates missing data patterns.
* `reduce_imputed_data()` and `expand_imputed_data()` for storage.
* `format_pvalue()`, `format_estimate()`, `format_results_table()`.
* `ADEFF` and `ADMI` example datasets.
* Documentation via pkgdown.

# rbmiUtils 0.1.8

# rbmiUtils 0.1.7

* Moved `tidyr` from Suggests to Imports.
* Standardized examples to use native pipe operator `|>`.

# rbmiUtils 0.1.6

* Added additional tests for all utility functions.

# rbmiUtils 0.1.4

* First release. Preparation for CRAN submission.
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `\dontrun{}` for slow examples | `\donttest{}` for slow examples | roxygen2 7.0+ / CRAN policy | CRAN prefers `\donttest{}`; `\dontrun{}` reserved for truly non-runnable code |
| Manual Rd `\figure{}` tags | roxygen2 markdown `\if{html}{\figure{...}}` | roxygen2 7.1.0 (2020) | Markdown-based figure inclusion in roxygen comments |
| Inline HTML tables in README | Static image screenshots of gt tables | gt 0.8+ with gtsave | github_document cannot render HTML widgets; images work everywhere |
| `knitr::knit()` for README | `devtools::build_readme()` | devtools 2.4+ | Renders against current package source, ensures consistency |

**Deprecated/outdated:**
- `\dontrun{}` for slow-but-working examples: CRAN now flags this; use `\donttest{}` instead
- `rmarkdown::render()` for README: `devtools::build_readme()` is the current standard

## Open Questions

1. **Vignette build time with ADEFF pipeline**
   - What we know: ADEFF has 1000 rows (200 subjects x 5 visits), and the vignette would use `n_samples = 100, warmup = 200, thin = 2`. The existing analyse2 vignette runs in acceptable time using the pre-computed ADMI dataset.
   - What's unclear: Exact build time for `draws()` + `impute()` on ADEFF with these parameters. May vary by machine.
   - Recommendation: Build the vignette with live code first. If build time exceeds 2 minutes, switch to precompute pattern (rename to `.Rmd.orig`, add `precompile.R` script). Monitor with `system.time(devtools::build_vignettes())`.

2. **Image size impact on package**
   - What we know: CRAN has a 5MB tarball size guideline. Each PNG image adds ~50-200KB.
   - What's unclear: Current package tarball size and remaining budget.
   - Recommendation: Use reasonable image dimensions (800x400 for forest plots, 600x400 for tables) and PNG compression. Monitor with `devtools::build()` and check tarball size. If too large, convert the vignette to a pkgdown article instead.

3. **rbmi pkgdown site URL stability**
   - What we know: rbmi was historically at `insightsengineering.github.io/rbmi` but the package is under `openpharma` organization. The CRAN page is stable.
   - What's unclear: Whether the rbmi pkgdown site URL is still active or has moved.
   - Recommendation: Use CRAN URLs for rbmi vignette links (`https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html`) which are guaranteed stable. For beeca, the pkgdown site at `https://openpharma.github.io/beeca/` appears current.

## Sources

### Primary (HIGH confidence)
- R Packages (2e), Chapter 17: Vignettes (https://r-pkgs.org/vignettes.html) - vignette structure, articles vs vignettes, CRAN considerations
- R Packages (2e), Chapter 18: Other Markdown Files (https://r-pkgs.org/other-markdown.html) - README.Rmd patterns, fig.path, build_readme
- roxygen2 rd-formatting vignette (https://roxygen2.r-lib.org/articles/rd-formatting.html) - figure inclusion in man pages
- Tidyverse Style Guide, Chapter 10: NEWS (https://style.tidyverse.org/news.html) - NEWS.md formatting conventions
- gt::gtsave() documentation (https://gt.rstudio.com/reference/gtsave.html) - saving gt tables as PNG
- pkgdown auto-linking (https://pkgdown.r-lib.org/articles/linking.html) - cross-package reference syntax

### Secondary (MEDIUM confidence)
- rOpenSci precompute vignettes blog post (https://ropensci.org/blog/2019/12/08/precompute-vignettes/) - .Rmd.orig pattern
- R-hub blog on examples (https://blog.r-hub.io/2020/01/27/examples/) - donttest vs dontrun guidance
- R-hub blog on NEWS (https://blog.r-hub.io/2020/05/08/pkg-news/) - NEWS.md conventions

### Tertiary (LOW confidence)
- None -- all findings verified with official documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all tools are standard R package infrastructure already in use
- Architecture: HIGH - patterns verified against R Packages (2e) and existing package structure
- Pitfalls: HIGH - well-documented issues in R packaging community
- Discretion recommendations: MEDIUM - based on judgment informed by best practices, but subjective

**Research date:** 2026-02-08
**Valid until:** 2026-04-08 (stable domain; R packaging patterns change slowly)
