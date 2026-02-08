# Phase 7: Site Polish - Research

**Researched:** 2026-02-08
**Domain:** pkgdown site configuration (navbar, reference grouping, open graph, footer, logo)
**Confidence:** HIGH

## Summary

Phase 7 covers five requirements for making the rbmiUtils pkgdown site professional: hex logo display (SITE-01), custom navbar (SITE-02), grouped function reference (SITE-03), social/open graph cards (SITE-04), and a custom footer with openpharma/pharmaverse links (SITE-05). All five requirements are implemented through a single `_pkgdown.yml` configuration file, with minor supporting tasks for the logo file and favicons.

The project already has a hex logo (`man/figures/rbmiUtils.png`) but it is not named following the pkgdown convention (`logo.png` or `logo.svg`). The existing `_pkgdown.yml` is minimal (only URL and bootstrap 5). The package exports 22 functions across 7 clear functional layers plus 2 datasets, the vignettes are already built, and NEWS.md exists. The pkgdown version in use is 2.1.1 with Bootstrap 5.

All five requirements are achievable purely through `_pkgdown.yml` YAML configuration plus one file rename (logo). No custom R code, template packages, or external dependencies are needed. The entire phase is configuration work, not development work.

**Primary recommendation:** Rename the existing hex logo to `logo.png`, then build out `_pkgdown.yml` with navbar, reference groups, opengraph, and footer sections. Run `pkgdown::build_favicons()` after the logo rename to generate the favicon set.

## Standard Stack

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| pkgdown | 2.1.1 | Site generation, all configuration | Already in use, all SITE-0x requirements are pkgdown config |
| `_pkgdown.yml` | N/A | Single configuration file for all 5 requirements | Standard pkgdown configuration mechanism |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| `usethis::use_logo()` | any | Properly scale and place logo at `man/figures/logo.png` | If logo needs resizing (current logo may already be correct size) |
| `pkgdown::build_favicons()` | 2.1.1 | Generate favicon set from logo for browser tab icons | After logo is in place as `logo.png` |
| `pkgdown::build_reference_index()` | 2.1.1 | Quickly iterate on reference grouping | During development of reference groups |
| `pkgdown::build_site()` | 2.1.1 | Full site rebuild to verify all changes | Final verification |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Renaming rbmiUtils.png to logo.png | Keeping current name | pkgdown auto-detects `logo.png`/`logo.svg` for navbar and favicons; using non-standard name breaks auto-detection |
| usethis::use_logo() | Manual file copy/rename | usethis handles resizing to hex spec (181x209 at 2x retina) and generates README markdown; manual approach works if dimensions are already correct |
| Custom template package | Inline _pkgdown.yml config | Template packages are for organizations with many packages needing consistent branding; overkill for a single package |

## Architecture Patterns

### _pkgdown.yml Complete Structure
The final `_pkgdown.yml` should have this overall structure:

```yaml
url: https://openpharma.github.io/rbmiUtils/

template:
  bootstrap: 5
  opengraph:
    image:
      src: man/figures/logo.png
      alt: "rbmiUtils hex logo - Utilities and Extensions for rbmi"
    twitter:
      card: summary

home:
  title: "rbmiUtils: Utility Functions to Support and Extend rbmi"
  description: >
    Bridges rbmi analysis results into publication-ready regulatory tables
    and forest plots for clinical trial workflows.

navbar:
  structure:
    left: [intro, reference, articles, news]
    right: [search, github, lightswitch]
  components:
    intro:
      text: Get Started
      href: articles/pipeline.html
    articles:
      text: Articles
      menu:
        - text: "End-to-End Pipeline"
        - text: "From rbmi Analysis to Regulatory Tables"
          href: articles/pipeline.html
        - text: "-------"
        - text: "Workflow Guides"
        - text: "Storing and Analyzing Imputed Data"
          href: articles/analyse2.html
        - text: "Data Preparation and Validation"
          href: articles/data-preparation.html
        - text: "Efficient Storage of Imputed Data"
          href: articles/efficient-storage.html

reference:
  - title: "Data Preparation"
    desc: "Validate and prepare data before imputation"
    contents:
      - validate_data
      - prepare_data_ice
      - summarise_missingness
  - title: "Analysis"
    desc: "Apply analysis functions across imputed datasets"
    contents:
      - analyse_mi_data
      - gcomp_responder
      - gcomp_responder_multi
      - gcomp_binary
  - title: "Tidying"
    desc: "Tidy and extract results from pooled objects"
    contents:
      - tidy_pool_obj
      - extract_trt_effects
      - extract_lsm
  - title: "Reporting"
    desc: "Publication-ready tables and plots"
    contents:
      - efficacy_table
      - plot_forest
  - title: "Formatting"
    desc: "Format values for publication"
    contents:
      - format_pvalue
      - format_estimate
      - format_results_table
      - format_results
  - title: "Storage"
    desc: "Efficient storage of imputed datasets"
    contents:
      - reduce_imputed_data
      - expand_imputed_data
  - title: "Utilities"
    desc: "Helper functions for imputed data workflows"
    contents:
      - get_imputed_data
      - create_impid
      - combine_results
      - pool_to_ard
  - title: "Print & Summary Methods"
    desc: "S3 methods for rbmiUtils objects"
    contents:
      - starts_with("print.")
      - starts_with("summary.")
  - title: "Datasets"
    desc: "Example datasets for demonstrations"
    contents:
      - ADEFF
      - ADMI

footer:
  structure:
    left: developed_by
    right: [built_with, community]
  components:
    community: >
      Part of the [openpharma](https://openpharma.github.io/) ecosystem
      | [pharmaverse](https://pharmaverse.org/)
```

### Pattern 1: Logo as `logo.png` in man/figures
**What:** pkgdown auto-detects `man/figures/logo.png` (or `logo.svg`) and uses it for: (1) navbar display, (2) favicon generation, (3) default Open Graph image, (4) README display.
**When to use:** Always -- this is the standard convention.
**Current state:** The project has `man/figures/rbmiUtils.png` which must be renamed (or copied) to `man/figures/logo.png`.
**Impact on README:** The README.Rmd and README.md both reference `man/figures/rbmiUtils.png`. After renaming, update the reference to `man/figures/logo.png`.

### Pattern 2: Reference Grouping with Explicit Function Lists
**What:** List each function explicitly in the `reference:` section rather than using `starts_with()` selectors.
**When to use:** When the package has a modest number of exports (22 functions + 2 datasets) and the grouping does not follow naming conventions.
**Why:** rbmiUtils function names do not follow a prefix convention (e.g., they are not named `dp_validate_data`, `dp_prepare_data_ice`), so `starts_with()` selectors would not work for most groups. Explicit listing is clearer and easier to maintain.
**Exception:** `starts_with("print.")` and `starts_with("summary.")` work well for S3 methods.

### Pattern 3: Footer with Custom Components
**What:** Define custom `components` in the footer alongside the built-in `developed_by` and `built_with` components. Each component is markdown text that gets converted to HTML.
**When to use:** When you need links beyond the default author/pkgdown attribution.

### Anti-Patterns to Avoid
- **Overriding navbar defaults unnecessarily:** Only override the specific components you need to change. If you override the entire navbar, you lose future pkgdown default improvements.
- **Using `has_concept()` without roxygen `@family` tags:** The `has_concept()` selector only works if functions have `@family` or `@concept` tags in their roxygen documentation. This package does not use them, so use explicit function lists instead.
- **Putting logo in a non-standard location:** Always use `man/figures/logo.png` -- pkgdown, usethis, and the R package ecosystem all expect this location.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Favicon generation | Manual favicon creation with image editing tools | `pkgdown::build_favicons()` | Generates complete favicon set (all sizes/formats) via realfavicongenerator.net API |
| Logo resizing | Manual image resizing | `usethis::use_logo("path/to/logo.png")` | Handles hex spec dimensions (181x209 retina), copies to `man/figures/logo.png`, generates README markdown |
| Open Graph meta tags | Manual `<meta>` tag HTML in templates | `template: opengraph:` in `_pkgdown.yml` | pkgdown generates correct og:image, og:title, og:description meta tags automatically |
| Navbar HTML | Custom HTML template for navigation | `navbar:` section in `_pkgdown.yml` | pkgdown handles responsive design, Bootstrap 5 integration, and accessibility |
| Social media card preview image | Custom image for social cards | Default logo-based card from pkgdown | pkgdown uses the logo automatically if no custom image is specified |

**Key insight:** Every requirement in this phase is a pkgdown configuration concern. There is zero need for custom HTML templates, custom CSS, or custom R code. The entire phase is `_pkgdown.yml` editing plus a logo file rename.

## Common Pitfalls

### Pitfall 1: Logo File Not Named `logo.png`
**What goes wrong:** pkgdown cannot auto-detect the logo for navbar, favicons, or Open Graph cards. The site builds without a logo in the navbar, favicons are missing, and social cards show no image.
**Why it happens:** The package was created with a custom logo filename (`rbmiUtils.png`) before adopting pkgdown conventions.
**How to avoid:** Rename (or copy) `man/figures/rbmiUtils.png` to `man/figures/logo.png`. Update any references in README.Rmd.
**Warning signs:** No logo in navbar after `build_site()`, no favicon in browser tab.

### Pitfall 2: Missing Functions in Reference Groups
**What goes wrong:** pkgdown throws an error if any exported function is not included in at least one reference group, or warns if a function appears in multiple groups.
**Why it happens:** When explicitly listing functions, it is easy to miss one, especially after adding new exports.
**How to avoid:** After writing the reference section, run `pkgdown::build_reference_index()` which will report any unmatched topics. Every topic in the NAMESPACE must appear in exactly one group.
**Warning signs:** Build warnings like "Topics missing from index" or "Topics found in multiple groups".

### Pitfall 3: Open Graph Image URL Must Be Absolute
**What goes wrong:** Social media platforms cannot resolve relative image paths. The og:image meta tag must contain a full URL for cards to display properly.
**Why it happens:** Using a relative path like `man/figures/logo.png` without the `url:` field being set in `_pkgdown.yml`.
**How to avoid:** Ensure `url: https://openpharma.github.io/rbmiUtils/` is set in `_pkgdown.yml` (it already is). pkgdown will construct the absolute URL from this base URL plus the image path. Alternatively, use a fully qualified URL in the `opengraph.image.src` field.
**Warning signs:** Social media debugger tools (Twitter Card Validator, Facebook Sharing Debugger, opengraph.xyz) show no preview image.

### Pitfall 4: Footer Markdown Not Rendering as Expected
**What goes wrong:** Footer components are joined with spaces and then converted from markdown to HTML. If you use block-level markdown (headers, lists), the rendering may break.
**Why it happens:** The footer rendering expects inline markdown (links, bold, italic) not block-level elements.
**How to avoid:** Use only inline markdown in footer components: links `[text](url)`, bold `**text**`, italic `*text*`. Use ` | ` as a visual separator between items.
**Warning signs:** Footer displays raw markdown text instead of rendered HTML.

### Pitfall 5: Navbar Articles Menu Order
**What goes wrong:** Articles appear in an unexpected order in the dropdown menu.
**Why it happens:** When you override the `articles` component with a custom menu, the order is exactly as specified in the YAML. If you don't override, pkgdown uses alphabetical order of vignette filenames.
**How to avoid:** Explicitly define the `articles` component with the desired order using `menu:` items.
**Warning signs:** Pipeline vignette not appearing first in the articles dropdown.

### Pitfall 6: Forgetting to Rebuild After Config Changes
**What goes wrong:** Changes to `_pkgdown.yml` are not reflected in the built site.
**Why it happens:** The `docs/` directory contains previously built output. Config changes require a site rebuild.
**How to avoid:** After modifying `_pkgdown.yml`, run `pkgdown::build_site()` (or targeted builds like `pkgdown::build_reference_index()`, `pkgdown::build_home()`).
**Warning signs:** Old site content still showing after YAML changes.

## Code Examples

### Example 1: Minimal Complete `_pkgdown.yml`
```yaml
# Source: pkgdown.r-lib.org/articles/customise.html + metadata.html
url: https://openpharma.github.io/rbmiUtils/

template:
  bootstrap: 5
  opengraph:
    image:
      src: man/figures/logo.png
      alt: "rbmiUtils hex logo"
    twitter:
      card: summary
```

### Example 2: Navbar with Custom Articles Menu
```yaml
# Source: pkgdown.r-lib.org/articles/customise.html
navbar:
  structure:
    left: [intro, reference, articles, news]
    right: [search, github, lightswitch]
  components:
    intro:
      text: Get Started
      href: articles/pipeline.html
    articles:
      text: Articles
      menu:
        - text: "End-to-End Pipeline"
        - text: "From rbmi Analysis to Regulatory Tables"
          href: articles/pipeline.html
        - text: "-------"
        - text: "Workflow Guides"
        - text: "Storing and Analyzing Imputed Data"
          href: articles/analyse2.html
        - text: "Data Preparation and Validation"
          href: articles/data-preparation.html
        - text: "Efficient Storage of Imputed Data"
          href: articles/efficient-storage.html
```

### Example 3: Reference Groups
```yaml
# Source: pkgdown.r-lib.org/articles/pkgdown.html, r-pkgs.org/website.html
reference:
  - title: "Data Preparation"
    desc: "Validate and prepare data before imputation"
    contents:
      - validate_data
      - prepare_data_ice
      - summarise_missingness
  - title: "Analysis"
    desc: "Apply analysis functions across imputed datasets"
    contents:
      - analyse_mi_data
      - gcomp_responder
      - gcomp_responder_multi
      - gcomp_binary
```

### Example 4: Custom Footer with Community Links
```yaml
# Source: pkgdown.r-lib.org/articles/customise.html
footer:
  structure:
    left: developed_by
    right: [built_with, community]
  components:
    community: >
      Part of the [openpharma](https://openpharma.github.io/) ecosystem
      | [pharmaverse](https://pharmaverse.org/)
```

### Example 5: Logo Rename and Favicon Generation (R commands)
```r
# Option A: If logo dimensions are already correct
file.copy("man/figures/rbmiUtils.png", "man/figures/logo.png")
pkgdown::build_favicons(overwrite = TRUE)

# Option B: Use usethis to handle resizing
usethis::use_logo("man/figures/rbmiUtils.png")
# This creates man/figures/logo.png at correct dimensions
# and prints README markdown snippet
pkgdown::build_favicons(overwrite = TRUE)
```

### Example 6: Open Graph Verification
```r
# After building site, verify Open Graph tags
# Check generated HTML for meta tags:
readLines("docs/index.html") |>
  grep("og:|twitter:", value = TRUE)

# Or use external tools:
# - https://www.opengraph.xyz/
# - https://cards-dev.twitter.com/validator
# - https://developers.facebook.com/tools/debug/
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Bootstrap 3 templates | Bootstrap 5 (default since pkgdown 2.0) | pkgdown 2.0.0 (Dec 2021) | Modern responsive design, dark mode support |
| `template: default` | `template: bootstrap: 5` | pkgdown 2.0.0 | Must explicitly opt in (or is default for new sites) |
| No footer customization | `footer: structure: + components:` | pkgdown 2.0.0 | Full control over footer content |
| No lightswitch | `lightswitch` navbar component | pkgdown 2.0.0+ | Dark/light/auto mode toggle |
| Manual og:image meta | `template: opengraph:` YAML | pkgdown 2.0.0 | Automatic Open Graph metadata generation |

**Deprecated/outdated:**
- Bootstrap 3 templates: pkgdown 2.0+ defaults to Bootstrap 5. The project already uses BS5.
- `_pkgdown.yml` `home: strip_header: true`: No longer needed in BS5 templates.

## Existing Project State

Critical facts about the current project state that the planner must account for:

| Item | Current State | Action Needed |
|------|--------------|---------------|
| Logo file | `man/figures/rbmiUtils.png` (hex sticker, exists) | Rename/copy to `man/figures/logo.png` |
| `_pkgdown.yml` | Minimal: URL + bootstrap 5 only | Expand with all 5 SITE-0x configurations |
| README.Rmd | References `man/figures/rbmiUtils.png` | Update to reference `man/figures/logo.png` |
| README.md | References `man/figures/rbmiUtils.png` | Regenerate from README.Rmd after update |
| Vignettes | 4 vignettes exist: pipeline, analyse2, data-preparation, efficient-storage | Map into navbar articles menu |
| NEWS.md | Exists with versioned entries | Will auto-detect for navbar news link |
| Exported functions | 22 functions + 2 datasets | All must be assigned to reference groups |
| S3 methods | `print.pool`, `summary.pool`, `print.analysis`, `summary.analysis` | Include in reference groups (use `starts_with()`) |
| Favicons | None generated yet (no `pkgdown/favicon/` directory) | Generate with `pkgdown::build_favicons()` |
| pkgdown version | 2.1.1 | Current, all features available |
| Site URL | `https://openpharma.github.io/rbmiUtils/` | Already set in `_pkgdown.yml` |
| GitHub URL | `https://github.com/openpharma/rbmiUtils` | In DESCRIPTION, auto-detected by pkgdown |
| `.Rbuildignore` | Already includes `^_pkgdown.yml$` | No changes needed |

### Complete Exported Function Inventory

For reference grouping, here is the full list of 22 exported functions organized by proposed layer:

**Data Preparation (3):**
- `validate_data` - Pre-flight validation before imputation
- `prepare_data_ice` - Prepare intercurrent event data
- `summarise_missingness` - Summarize missing data patterns

**Analysis (4):**
- `analyse_mi_data` - Apply analysis across imputed datasets
- `gcomp_responder` - G-computation for binary outcomes (single visit)
- `gcomp_responder_multi` - G-computation for binary outcomes (multiple visits)
- `gcomp_binary` - G-computation binary analysis wrapper

**Tidying (3):**
- `tidy_pool_obj` - Tidy pooled results into tibble
- `extract_trt_effects` - Extract treatment effects from results
- `extract_lsm` - Extract least squares means

**Reporting (2):**
- `efficacy_table` - Publication-ready gt efficacy table
- `plot_forest` - Three-panel forest plot

**Formatting (4):**
- `format_pvalue` - Format p-values for publication
- `format_estimate` - Format estimates with CIs
- `format_results_table` - Format full results table
- `format_results` - Format results with customization

**Storage (2):**
- `reduce_imputed_data` - Compress imputed data for storage
- `expand_imputed_data` - Restore full imputed data

**Utilities (4):**
- `get_imputed_data` - Extract imputed datasets from rbmi objects
- `create_impid` - Add IMPID column to imputed list
- `combine_results` - Combine multiple analysis results
- `pool_to_ard` - Convert pool object to pharmaverse ARD format

**S3 Methods (4, may or may not need explicit listing):**
- `print.pool` / `summary.pool`
- `print.analysis` / `summary.analysis`

**Datasets (2):**
- `ADEFF` - Example efficacy trial dataset
- `ADMI` - Example multiple imputation dataset

## Open Questions

1. **Should `rbmiUtils.png` be deleted after creating `logo.png`?**
   - What we know: pkgdown only auto-detects `logo.png`/`logo.svg`. Having both files wastes space but causes no functional issues.
   - What's unclear: Whether any downstream references (GitHub, CRAN) point to the old filename.
   - Recommendation: Keep both for now to avoid breaking any external links. Add a note for future cleanup.

2. **Should S3 print/summary methods appear in reference groups?**
   - What we know: These are exported and will cause pkgdown warnings if not listed. They are documented in man pages.
   - What's unclear: Whether users would look for them in the reference index.
   - Recommendation: Include them in a "Print & Summary Methods" group or fold them into their respective functional groups (pool methods with Tidying, analysis methods with Analysis).

3. **Twitter/X card configuration**
   - What we know: pkgdown supports `twitter: creator:` and `twitter: site:` fields for Twitter card metadata.
   - What's unclear: Whether the package maintainer has a Twitter/X handle to use, and whether the openpharma organization has one.
   - Recommendation: Use `card: summary` without creator/site handles. Can be added later if desired.

4. **Favicon generation requires network access**
   - What we know: `pkgdown::build_favicons()` calls the realfavicongenerator.net API. This requires network access during the build step.
   - What's unclear: Whether the build environment has unrestricted network access.
   - Recommendation: Run `build_favicons()` locally and commit the generated `pkgdown/favicon/` directory. The favicons only need regenerating when the logo changes.

## Sources

### Primary (HIGH confidence)
- [pkgdown customization docs](https://pkgdown.r-lib.org/articles/customise.html) - navbar structure, footer components, template includes
- [pkgdown metadata docs](https://pkgdown.r-lib.org/articles/metadata.html) - Open Graph configuration, social media cards
- [pkgdown introduction](https://pkgdown.r-lib.org/articles/pkgdown.html) - reference section grouping
- [pkgdown build_home](https://pkgdown.r-lib.org/reference/build_home.html) - logo auto-detection from man/figures/logo.png
- [pkgdown build_favicons](https://pkgdown.r-lib.org/reference/build_favicons.html) - favicon generation from logo
- [R Packages book, Chapter 19](https://r-pkgs.org/website.html) - reference grouping patterns
- [usethis::use_logo](https://usethis.r-lib.org/reference/use_logo.html) - logo sizing and placement

### Secondary (MEDIUM confidence)
- [admiral _pkgdown.yml](https://github.com/pharmaverse/admiral/blob/main/_pkgdown.yml) - pharmaverse reference grouping patterns, footer with organization links
- [pkgdown GitHub issue #1518](https://github.com/r-lib/pkgdown/issues/1518) - custom footer behavior

### Tertiary (LOW confidence)
- None -- all findings verified with official documentation.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all configuration via official pkgdown YAML, verified with official docs
- Architecture: HIGH - `_pkgdown.yml` structure verified with pkgdown.r-lib.org docs and real-world examples
- Pitfalls: HIGH - common issues documented in pkgdown docs and GitHub issues
- Function grouping: HIGH - based on direct inspection of NAMESPACE and R source files

**Research date:** 2026-02-08
**Valid until:** 2026-04-08 (pkgdown configuration is stable; 60-day validity)
